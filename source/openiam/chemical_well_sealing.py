# @author: Jaisree Iyer
# iyer5@llnl.gov

# Import libraries
import os
import sys
import logging
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import SystemModel, ComponentModel
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

from openiam.cfi.commons import process_parameters


class ChemicalWellSealing(ComponentModel):
    """
    The Chemical Well Sealing component is based on the model described in
    :cite:`IyerEtAl2020`. It predicts whether a fracture
    at the cement caprock interface, upon exposure to |CO2|, would self-seal
    or not due to calcite precipitation. The model couples two-phase flow
    of supercritical |CO2| and brine through fractures, advective and
    diffusive transport along the fracture, diffusive transport within
    the cement, and chemical reactions between cement and carbonated brine.
    If the fracture is predicted to self-seal, the time required for sealing
    is also computed.
    The original model is described in :cite:`WalshEtAl2013`,
    :cite:`WalshEtAl2014a`, and :cite:`IyerEtAl2017` and was
    calibrated using experimental data presented in :cite:`WalshEtAl2013`,
    :cite:`WalshEtAl2014a`, and :cite:`WalshEtAl2014b`.

    Component model input definitions:

    * **fractureAperture** [|m|] (1.0e-5 to 2.0e-3) - aperture of the fractured
      leakage path (default: 2.0e-5). Any input aperture that is lower than 10 micron
      (lower bound) is set to 10 micron.

    * **fractureLength** [|m|] (10.0 to 400.0) - length of the fractured leakage path
      (default: 20)

    * **maxOverpressure** [|Pa|] (1.0e+6 to 1.5e+7) - maximum overpressure
      the base of the fracture is expected to experience (default: 5.0e+6).

    The output from the Chemical Well Sealing component informs about
    the sealing ability of the fractured leakage pathway. In case the fracture
    seals, the component would also report sealing time.

    * **seal_flag** [-] - flag informing whether a fracture would seal (1)
      or not (0) due to calcite precipitation

    * **seal_time** [|s|] - predicted time for sealing a fracture by calcite
      recipitation. If fracture doesn't seal this variable is set to 0.0.

    """
    def __init__(self, name, parent):
        """
        Constructor method of ChemicalWellSealing class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :returns: ChemicalWellSealing class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 0.0}

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'ChemicalWellSealing'

        # Define output dictionary labels
        self.output_labels = ['seal_flag', 'seal_time']

        # Set default parameters of the component model
        self.add_default_par('fractureAperture', value=2.0e-5)
        self.add_default_par('fractureLength', value=20)
        self.add_default_par('maxOverpressure', value=5.0e+6)

        # Define dictionary of parameters boundaries
        self.pars_bounds = dict()
        self.pars_bounds['fractureAperture'] = [1.0e-5, 2.0e-3]
        self.pars_bounds['fractureLength'] = [10.0, 400.0]
        self.pars_bounds['maxOverpressure'] = [1.0e+6, 1.5e+7]

        # Indicate how often the component should be run
        self.default_run_frequency = 1
        self.run_frequency = 1 # run only once

    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        debug_msg = 'Input parameters of component {} are {}'.format(self.name, p)
        logging.debug(debug_msg)

        for key, val in p.items():
            if key in self.pars_bounds:
                if val > self.pars_bounds[key][1]:
                    error_msg = ''.join([
                        'Parameter {} of ChemicalWellSealing component is out of ',
                        'boundaries. The parameter values do not satisfy ',
                        'the assumptions of the component.']).format(key)
                    logging.error(error_msg)
                    raise ValueError(error_msg)

                if val < self.pars_bounds[key][0]:
                    if key in ['fractureAperture']:
                        warn_msg = ''.join([
                            'Parameter {} of ChemicalWellSealing component ',
                            'will be modified to satisfy the component ',
                            'requirements.']).format(key)
                        logging.warning(warn_msg)
                    else:
                        warn_msg = ''.join([
                            'Parameter {} of ChemicalWellSealing component ',
                            'is out of boundaries. Reported results are ',
                            'obtained by extrapolation and may be ',
                            'subject to errors.']).format(key)
                        logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter of ',
                    'ChemicalWellSealing component.']).format(key)
                logging.warning(warn_msg)


    def simulation_model(self, p, time_point=0.0):
        """
        Predicts whether a fractured pathway along the cement-caprock interface
        in a well would self-seal due to calcite precipitation. If the fracture
        is predicted to self-seal, the sealing time will also be reported.

        :param p: input parameters of ChemicalWellSealing component
        :type p: dict

        :returns: dictionary of predictions of the chemical well sealing component
             keys: ['seal_flag', 'seal_time']
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Convert the units of the input parameters
        aperture = actual_p["fractureAperture"]*1.0e6
        length = actual_p["fractureLength"]
        overpressure = actual_p["maxOverpressure"]/1.0e6

        # Run the component only once to identify which leakage pathways
        # will seal and which ones will not.
        if time_point == 0.0:
            # run the ROM
            output = run_coupled_model(aperture, length, overpressure)
            out = {}
            out['seal_flag'] = output['seal_flag']
            out['seal_time'] = output['seal_time']*24*60*60

            # Exit method
            return out

    def connect_with_system(self, component_data, name2obj_dict, *args, **kwargs):
        """
        Code to add chemical well sealing to system model for control file interface.

        :param component_data: Dictionary of component data
        :type component_data: dict

        :param name2obj_dict: Dictionary to map component names to objects
        :type name2obj: dict

        :returns: None
        """
        # Process parameters data if any
        process_parameters(self, component_data, name2obj_dict)

        if 'Outputs' in component_data:
            comp_outputs = component_data['Outputs']
            for obs_nm in comp_outputs:
                self.add_obs(obs_nm, index=[0])


def get_critical_brine_residence_time(logAperture):
    """
    Function to calculate the critical brine residence time for a given fracture
    aperture above which a fracture at the cement-caprock interface would
    self-seal. The critical residence time has been numerically computed
    for several aperture values (Table S3 in the supplementary information of the
    paper https://doi.org/10.1021/acs.est.9b05039). For apertures not listed
    in the table, the critical brine residence time is computed by taking
    the natural logarithm of the values in the table and performing linear
    interpolation. For smaller apertures, the aperture value is set to 10 microns.

    Input:
        logAperture: natural log of the aperture in microns

    Output:
        criticalResidenceTime: critical brine residence time in seconds.
    """
    # Import interpolation function
    from scipy import interpolate

    # Natural log of Table S3 in the supplementary information of the paper
    # https://doi.org/10.1021/acs.est.9b05039
    fractureAperture = [np.log(10), np.log(50), np.log(100), np.log(250),
                        np.log(500), np.log(1000), np.log(2000)]
    residenceTime = [np.log(0.04), np.log(0.8), np.log(4), np.log(16),
                     np.log(80), np.log(200), np.log(1000)]

    # Linear interpolation
    f = interpolate.interp1d(fractureAperture, residenceTime)

    # Critical brine residence time in seconds
    if logAperture < np.log(10):
        logAperture = np.log(10)
    criticalResidenceTime = np.exp(f(logAperture))*60
    return criticalResidenceTime


def run_coupled_model(aperture, length, overpressure):
    """
    Function to predict whether a fracture would self-seal or not. If the fracture
    is predicted to self-seal, the time required for sealing is also computed.
    The calculations are based on the study published at:
        https://doi.org/10.1021/acs.est.9b05039

    Inputs:

    :aperture: fracture aperture in microns (between 10 and 2000 microns)
    :length: length of the fracture in metres (between 10 and 400 m)
    :overpressure: maximum overpressure (in MPa) in the reservoir
        (bottom of the leakage path) due to CO2 injection (between 1 and 15 MPa)

    Output:

    :outputDict:
        seal_flag - Flag informing whether a the fracture would seal (1) or not (0)
        seal_time - Time in days for a fracture to seal. If fracture doesn't seal
        this variable is set to 0
    """
    # Initialize output
    outputDict = {"seal_flag": 0, "seal_time": 0.0}

    # Calculations to determine the brine residence time
    viscosity = 5.154*0.0001 # viscosity of 1m NaCl solution at 60 C
    aperture = max(aperture, 10.0)
    maximumBrineVelocity = (aperture*10**(-6))**2*overpressure*10**6/12/viscosity/length
    minimumBrineResidenceTime = length/maximumBrineVelocity

    # Function to calculate the critical brine residence time for sealing
    critcalBrineResidenceTime = get_critical_brine_residence_time(np.log(aperture))

    # Comparison between brine residence time and critical brine residence time.
    if minimumBrineResidenceTime > critcalBrineResidenceTime:
        outputDict["seal_flag"] = 1
        outputDict["seal_time"] = 9.54*10**(-3)*aperture**1.79

    return outputDict


if __name__ == "__main__":
    # Create system model
    time_array = 365.25*np.arange(0.0, 2.0)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add chemical well sealing ROM and define parameters
    cws = sm.add_component_model_object(ChemicalWellSealing(name='cws', parent=sm))
    cws.add_par('fractureAperture', value=0.000043)
    cws.add_par('fractureLength', value=160)
    cws.add_par('maxOverpressure', value=14000000)

    # Add observations (output) from the chemical well sealing ROM
    cws.add_obs('seal_flag', index=[0])
    cws.add_obs('seal_time', index=[0])

    # Run the system model
    sm.forward()

    # Print the observations
    print('seal_flag:', sm.obs['cws.seal_flag_0'].sim)
    print('seal_time:', sm.obs['cws.seal_time_0'].sim)

    # Expected output
    # seal_flag: 1
    # seal_time: 691784.2222

    #
    # Test Cases
    # fractureAperture (m)  fractureLength (m)  maxOverpressure (Pa)  seal_flag  seal_time (s)  Error/Warning  Message
    # 0.00066               25                  6800000               0         0                No
    # 0.0014                70                  1600000               0         0                No
    # 0.00021               200                 13000000              0         0                No
    # 0.000043              160                 14000000              1         691784.2222      No
    # 0.000018              97                  3600000               1         145545.875       No
    # 0.000015              396                 1800000               1         105018.4178      No
    # 0.0000046             77                  11000000              1         50823.21299      Warning        Parameter fractureAperture of ChemicalSealingWellbore component will be modified to satisfy the component requirements.
    # 0.000042              5                   5200000               0         0                Warning        Parameter fractureLength of ChemicalSealingWellbore component is out of boundaries. Reported results are obtained by extrapolation and may be subject to errors.
    # 0.001979              14                  500000                0         0                Warning        Parameter maxOverpressure of ChemicalSealingWellbore component is out of boundaries. Reported results are obtained by extrapolation and may be subject to errors.
    # 0.003                 223                 2100000                                          Error          Parameter fractureAperture of ChemicalSealingWellbore component is out of boundaries. The parameter values do not satisfy the assumptions of the component.
    # 0.000159              600                 6410000                                          Error          Parameter fractureLength of ChemicalSealingWellbore component is out of boundaries. The parameter values do not satisfy the assumptions of the component.
    # 0.000246              378                 18100000                                         Error          Parameter maxOverpressure of ChemicalSealingWellbore component is out of boundaries. The parameter values do not satisfy the assumptions of the component.
    #
