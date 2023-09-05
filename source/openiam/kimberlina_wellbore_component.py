# -*- coding: utf-8 -*-
# Created on Nov 13, 2018
# @author: Kayyum Mansoor
# mansoor1@llnl.gov
import os
import sys
import logging
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import SystemModel, ComponentModel
except ImportError as err:
    print('Unable to load IAM class module: '+str(err))


class KimberlinaWellbore(ComponentModel):
    """
    The module contains the solution class for the Kimberlina wellbore model
    based on deep Kimberlina reservoir simulations

    Component model input definitions:

    * **time** * [years] (0 to 200) - simulation time (default: 10)

    * **pressure** * [|Pa|] (17130510 to 37483630) - bottom hole pressure
      (default: 28228195.8)

    * **CO2saturation** * [-] (0.0 to 0.634324) - bottom hole |CO2| saturation
      (default: 0.45741)

    * **logK_well_permeability** * [|log10| |m^2|] (-11.8468 to -10.001) - well permeability
      (default: -11.5926)

    * **previous_co2_mass** * [|kg|] (1.58e-3 to 1.62e+9) - cumulative brine mass
      from previous timestep (default: 3.49e+5)

    * **previous_brine_mass** * [|kg|] (2.37e+1 to 8.70e+6) - cumulative |CO2| mass
      from previous timestep (default: 2.75e+3)

    Observations of the Kimberlina Wellbore ROM component are:

    * **CO2_mass** [|kg|] - cumulative |CO2| mass

    * **brine_mass** [|kg|] - cumulative brine mass

    * **CO2_rate** [|kg/s|] - |CO2| rate

    * **brine_rate** [|kg/s|] - brine rate.
    """
    def __init__(self, name, parent):
        """
        Constructor method of KimberlinaWellbore class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :returns: KimberlinaWellbore class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25, 'time_step': 365.25}

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'KimberlinaWellbore'

        # Default parameters
        self.add_default_par('time', value=10)
        self.add_default_par('dt', value=365.25)
        self.add_default_par('pressure', value=28228195.8)
        self.add_default_par('CO2saturation', value=0.45741)
        self.add_default_par('logK_well_permeability', value=-11.5926)
        self.add_default_par('previous_co2_mass', value=348953.5334519825)
        self.add_default_par('previous_brine_mass', value=2754.645894381841)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['logK_well_permeability'] = [-11.8468, -10.001]
        self.pars_bounds['time'] = [0.0, 200.0]
        self.pars_bounds['dt'] = [1, 200.0]

        # Define dictionary of temporal data limits
        self.temp_data_bounds = dict()
        self.temp_data_bounds['pressure'] = [17130510, 37483630]
        self.temp_data_bounds['CO2saturation'] = [0, 0.634324]
        self.temp_data_bounds['previous_co2_mass'] = [10**(-2.80192602858),
                                                      10**(9.20977325986)]
        self.temp_data_bounds['previous_brine_mass'] = [10**(1.37524597484),
                                                        10**(6.93933611696)]
        self.temp_data_bounds['time'] = [1, 200]

        # Initialize previous state data here
        self._previous_co2_mass = 0.01
        self._previous_brine_mass = 0.01

        # Define output dictionary labels
        self.output_labels = ['CO2_mass', 'brine_mass', 'CO2_rate', 'brine_rate']

        # Check if another Kimberlina wellbore component is already part of the system model
        orig_kwc_cmpnt = self._parent.ml_models_components.get(
            'KimberlinaWellbore', None)

        try:
            import components.wellbore.kimberlina_wellbore.kimberlina_wellbore_rom as kwbrom
        except ImportError:
            print('\nERROR: Unable to load ROM for Kimberlina Wellbore component\n')
            sys.exit()

        if orig_kwc_cmpnt is None:
            # Register the component as the one using ml models
            # The component will be registered only if no other Kimberlina wellbore
            # components were added to the same system model
            self._parent.register_ml_model_component(self, 'KimberlinaWellbore')

            # Initiate solution object
            self.sol = kwbrom.Solution()
        else:
            # Get solution class object from existing component of the same type
            self.sol = orig_kwc_cmpnt.sol

        debug_msg = 'KimberlinaWellbore created with name {}'.format(name)
        logging.debug(debug_msg)

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
                if ((val < self.pars_bounds[key][0]) or
                        (val > self.pars_bounds[key][1])):
                    warn_msg = 'Parameter {} is out of bounds.'.format(key)
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {key} not recognized as a Kimberlina ',
                    'wellbore component input parameter.']).format(key=key)
                logging.warning(warn_msg)

    def check_temporal_inputs(self, time, temp_inputs):
        """
        Check whether temporal data fall within specified boundaries.

        :param temp_inputs: temporal input data of component model
        :type temp_inputs: dict
        """
        for key, val in temp_inputs.items():
            if ((val < self.temp_data_bounds[key][0]) or
                    (val > self.temp_data_bounds[key][1])):
                warn_msg = ''.join([
                    '{0} {1} is outside the model range {2} to {3} ',
                    'at time t = {4}.']).format(key, val, self.temp_data_bounds[key][0],
                                                self.temp_data_bounds[key][1], time)
                logging.warning(warn_msg)

    def simulation_model(self, p, pressure=0.0, CO2saturation=0.0,
                         previous_co2_mass=0.0, previous_brine_mass=0.0,
                         time_point=365.25, time_step=365.25):
        """
        Return cumulative CO2 and brine mass from wellbores based on several metrics.

        :param p: input parameters of Kimberlina wellbore model
        :type p: dict

        :param pressure: Bottom hole pressure at reservoir/aquifer interface, Pa
        :type pressure: float

        :param CO2saturation: Bottom hole CO2 saturation at reservoir/aquifer interface, (-)
        :type CO2saturation: float

        :param time_point: time point in days at which the leakage rates are
            to be calculated; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param time_step: difference between the current and previous
            time points in days; by default, its value is 365.25 (1 year in days)
        :type time_step: float

        :returns: vol - dictionary of observations (leakage rates) of Kimberlina
            wellbore model; keys:['CO2_mass', 'brine_mass', 'CO2_rate', 'brine_rate']
        """
        # Check assign defailt values
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        actual_p.update(p)

        # Check whether co2, brine, rates and mass inputs satisfy the model requirements
        if (time_point >= 365.25) and ((
                previous_co2_mass >= self.temp_data_bounds['previous_co2_mass'][0]) or (
                    previous_brine_mass >= self.temp_data_bounds['previous_brine_mass'][0])):
            self.check_temporal_inputs(
                time_point, dict(list(zip(
                    ['CO2saturation', 'pressure',
                     'previous_co2_mass', 'previous_brine_mass', 'time'],
                    [CO2saturation, pressure,
                     previous_co2_mass, previous_brine_mass, time_point/365.25]))))

        p['CO2saturation'] = CO2saturation
        p['pressure'] = pressure

        # Use previous state variables
        if self._previous_co2_mass == 0.01:
            self._previous_co2_mass = previous_co2_mass
        if self._previous_brine_mass == 0.01:
            self._previous_brine_mass = previous_brine_mass

        p['previous_co2_mass'] = self._previous_co2_mass
        p['previous_brine_mass'] = self._previous_brine_mass

        # assign time
        actual_p['time'] = time_point/365.25
        actual_p['dt'] = time_step

        actual_p.update(p)

        labels = ['CO2_mass', 'brine_mass', 'CO2_rate', 'brine_rate']

        inputArray1 = [
            actual_p['pressure'], actual_p['CO2saturation'],
            actual_p['logK_well_permeability'], actual_p['previous_co2_mass'],
            actual_p['time'], actual_p['dt']]
        inputArray2 = [
            actual_p['pressure'], actual_p['CO2saturation'],
            actual_p['logK_well_permeability'], actual_p['previous_brine_mass'],
            actual_p['time'], actual_p['dt']]

        if (time_point > 0) and ((
                previous_co2_mass >= self.temp_data_bounds['previous_co2_mass'][0]) or (
                    previous_brine_mass >= self.temp_data_bounds['previous_brine_mass'][0])):
            self.sol.find(inputArray1, inputArray2)
            out = dict(list(zip(labels, self.sol.Outputs)))
        else:
            out = dict(list(zip(labels, np.zeros(len(labels)))))

        self._previous_co2_mass = out['CO2_mass']
        self._previous_brine_mass = out['brine_mass']

        return out

    # Attributes for system connections
    adapters = ['RateToMass']
    needsXY = False

    def connect_with_system(self):
        """
        Placeholder for the method needed for the control file interface
        """
        pass


if __name__ == "__main__":
    # Create system model
    time_array = 365.25*np.arange(0.0, 5.0)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    sm = SystemModel(model_kwargs=sm_model_kwargs)
    # Add Kimberlina wellbore model object and define parameters
    kwb = sm.add_component_model_object(KimberlinaWellbore(name='kwb', parent=sm))

    kwb.add_par('logK_well_permeability', value=-10.58833)
    #kwb.add_par('dt', value=365.25)

    kwb.model_kwargs['pressure'] = 28228195.8 # Pa
    kwb.model_kwargs['CO2saturation'] = 0.45741
    kwb.model_kwargs['previous_co2_mass'] = 348953.5334519825 # kg
    kwb.model_kwargs['previous_brine_mass'] = 2754.645894381841 # kg

    kwb.add_obs('CO2_mass')
    kwb.add_obs('brine_mass')
    kwb.add_obs('CO2_rate')
    kwb.add_obs('brine_rate')

    # Run the system model
    sm.forward()
    # Print the observations
    print('CO2 Mass:', sm.collect_observations_as_time_series(kwb, 'CO2_mass'))
    print('Brine Mass:', sm.collect_observations_as_time_series(kwb, 'brine_mass'))
    print('CO2 rate:', sm.collect_observations_as_time_series(kwb, 'CO2_rate'))
    print('Brine rate:', sm.collect_observations_as_time_series(kwb, 'brine_rate'))

    # The answer should be
    #  CO2 Mass: [0.00000000e+00 6.13886995e-01 1.36926020e+01 1.52343591e+02
    #   1.02803378e+03 4.78264325e+03 1.68179063e+04 4.77191207e+04
    #   1.14111262e+05 2.37647556e+05 4.42056955e+05 7.45573527e+05
    #   1.15958582e+06 1.68613190e+06 2.31780518e+06 3.03941245e+06
    #   3.83048357e+06 4.66806044e+06 5.52910487e+06 6.39242818e+06
    #   7.23971960e+06 8.05635167e+06 8.83131941e+06 9.55705401e+06
    #   1.02289227e+07 1.08449158e+07 1.14048977e+07 1.19104093e+07
    #   1.23639121e+07 1.27686270e+07 1.31280690e+07 1.34461047e+07
    #   1.37265581e+07 1.39731379e+07 1.41893445e+07 1.43786186e+07
    #   1.45439235e+07 1.46880592e+07 1.48135358e+07 1.49227061e+07
    #   1.50175221e+07 1.50998910e+07 1.51713025e+07 1.52332045e+07
    #   1.52868892e+07 1.53333296e+07 1.53735231e+07 1.54083343e+07
    #   1.54383879e+07 1.54643696e+07]
    #  Brine Mass: [   0.           39.96582597   52.93537418   68.05336453   85.31527973
    #    104.68967813  126.12447189  149.5523741   174.89591591  202.07045001
    #    230.98771363  261.55750683  293.69108466  327.30022321  362.30036031
    #    398.6086219   436.14889072  474.84685099  514.63351516  555.44373892
    #    597.21652591  639.89678384  683.43223338  727.775093    772.88098484
    #    818.70915403  865.22193899  912.38796202  960.17498587 1008.55613206
    #   1057.50637992 1107.00410634 1157.02804323 1207.56375748 1258.59574
    #   1310.10751666 1362.08938724 1414.53309455 1467.43195012 1520.77764422
    #   1574.56611666 1628.7948963  1683.46118782 1738.56497193 1794.10871213
    #   1850.09561893 1906.52471627 1963.4030871  2020.7342048  2078.52478327]
    #  CO2 rate: [0.00000000e+00 1.91360241e-08 4.14439469e-07 4.39358472e-06
    #   2.77489475e-05 1.18976395e-04 3.81374474e-04 9.79200394e-04
    #   2.10384001e-03 3.91462893e-03 6.47734298e-03 9.61785978e-03
    #   1.31192581e-02 1.66852382e-02 2.00165186e-02 2.28663543e-02
    #   2.50675313e-02 2.65412094e-02 2.72848516e-02 2.73570647e-02
    #   2.68490449e-02 2.58775087e-02 2.45572460e-02 2.29971417e-02
    #   2.12902352e-02 1.95196431e-02 1.77447533e-02 1.60186970e-02
    #   1.43706379e-02 1.28246408e-02 1.13900292e-02 1.00779438e-02
    #   8.88703230e-03 7.81364307e-03 6.85117273e-03 5.99773406e-03
    #   5.23819730e-03 4.56738313e-03 3.97611430e-03 3.45939875e-03
    #   3.00453520e-03 2.61011261e-03 2.26289628e-03 1.96155550e-03
    #   1.70116371e-03 1.47160773e-03 1.27365509e-03 1.10310036e-03
    #   9.52341350e-04 8.23310366e-04]
    #  Brine rate: [0.00000000e+00 1.26612372e-06 4.10980183e-07 4.79060206e-07
    #   5.46997085e-07 6.13937638e-07 6.79227627e-07 7.42385422e-07
    #   8.03088378e-07 8.61109023e-07 9.16332789e-07 9.68698291e-07
    #   1.01825164e-06 1.06500933e-06 1.10908742e-06 1.15053938e-06
    #   1.18957933e-06 1.22626436e-06 1.26076331e-06 1.29319795e-06
    #   1.32369974e-06 1.35245576e-06 1.37955515e-06 1.40514043e-06
    #   1.42931946e-06 1.45220705e-06 1.47390121e-06 1.49460108e-06
    #   1.51427941e-06 1.53310601e-06 1.55113975e-06 1.56848830e-06
    #   1.58516290e-06 1.60138015e-06 1.61710594e-06 1.63230970e-06
    #   1.64720608e-06 1.66184080e-06 1.67626358e-06 1.69042304e-06
    #   1.70445384e-06 1.71840633e-06 1.73227025e-06 1.74613355e-06
    #   1.76007492e-06 1.77411802e-06 1.78813019e-06 1.80236681e-06
    #   1.81671349e-06 1.83127293e-06]
