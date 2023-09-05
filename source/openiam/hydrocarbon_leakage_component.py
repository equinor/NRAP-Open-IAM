# @author: Nate Mitchell
# Created on February 10th, 2023
# Last Modified: March 2nd, 2023
# Nathaniel.Mitchell@netl.doe.gov
# Backend ROMs and ML models developed by Mohamed Mehana, Los Alamos National Laboratory
# mzm@lanl.gov

import os
import sys
import logging
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import ComponentModel
except ImportError as err:
    print('Unable to import ComponentModel class:', err)

import components.wellbore.hydrocarbon_leakage.hydrocarbon_leakage_ROM as hydrcarbrom
from openiam.cfi.commons import process_parameters
import openiam.cfi.strata as iam_strata

class HydrocarbonLeakage(ComponentModel):
    """
    The HydrocarbonLeakage component model is a reduced order model predicting
    liquid and gas leakage to the shallow aquifer between 100 to 410 years
    post-injection of |CO2| to depleted hydrocarbon field. The model output
    begins 100 years and extends to 410 years after injection stops. The component
    is based on a machine learning regression model fitted to the results
    of compositional multiphase flow transport simulations using a neural network.
    Total 192,000 data from ~1,000 numerical simulations were used to develop
    the model. The model predicts |CO2| and methane leakage in liquid and gas
    phases to a shallow aquifer. The model also predicts the total liquid (oil)
    and gas leakage to the shallow aquifer, where these total masses include
    all hydrocarbons (light, intermediate, and heavy) as well as |CO2|.
    The depth to the bottom of the shallow aquifer is assumed to be 60 ft
    (18.288 m) below the surface. The top of the shallow aquifer extends to the
    surface. Input parameters were sampled using Latin Hypercube Sampling
    across wide ranges.

    Values of input parameters FCO2, FC1, FC4, and FC7Plus must sum to one. To
    allow some leniency (e.g., for issues related to rounding errors), the sum
    of these values must be within 0.0001 (0.01%) of one (0.9999 to 1.0001).
    This option was created for cases when the sum of the values is different
    from 1 by a small value (e.g., 1.0e-6). If the sum of the provided values
    is not sufficiently close to one, then a warning message is printed.

    Since the temporal bounds for the HydrocarbonLeakage component are years
    100 to 410 after injection stops, the component produces zero results
    for any times outside of this range. Any results outside the
    applicable time range should not be considered valid, however.

    The description of the component's parameters is presented below:

    * **reservoirDepth** [|m|] (914.4 to 2743.2) - depth to the top of the
      reservoir (default: 2000); *linked to Stratigraphy*

    * **NTG** [-] (0.4 to 1.0) - net-to-gross ratio representing the fraction of
      reservoir contributing to the flow (default: 0.6)

    * **logResPerm** [log10 m^2] (-14.0057 to -13.0057) - logarithm of
      reservoir permeability (default: -13.5)

    * **reservoirPressureMult** [-] (1.0 to 1.2) - factor used to represent a
      state of the reservoir pressurization post-injection (relative to the
      reservoir pressure calculated lithostatically) (default: 1.1)

    * **logWellPerm** [log10 m^2] (-17.0057 to -12.0057) - logarithm of
      wellbore permeability (default: -13.0)

    * **avgWaterSaturation** [-] (0.471 to 0.792) - average water saturation
      in the reservoir (default: 0.5)

    * **FCO2** [-] (0.432 to 0.693) - mole fraction of |CO2| in the reservoir
      post-injection (default: 0.55)

    * **FC1** [-] (0.010 to 0.113) - mole fraction of methane in the reservoir
      post-injection (default: 0.05)

    * **FC4** [-] (0.010 to 0.111) - mole fraction of intermediate hydrocarbons
      in the reservoir post-injection  (default: 0.05)

    * **FC7Plus** [-] (0.123 to 0.500) - mole fraction of heavy hydrocarbons in
      the reservoir post-injection  (default: 0.35)

    Component model outputs:

    * **mass_oil_aquifer** [|kg|] - cumulative mass of oil leaked to the aquifer.
      This output includes all hydrocarbons (light, intermediate, and heavy
      hydrocarbons) and |CO2| in the liquid phase.

    * **mass_gas_aquifer** [|kg|] - cumulative mass of gas leaked to the aquifer.
      This output includes all hydrocarbons (light, intermediate, and heavy
      hydrocarbons) and |CO2| in the gas phase.

    * **mass_methane_gas_aquifer** [|kg|] - cumulative mass of methane gas leaked
      to the aquifer

    * **mass_methane_oil_aquifer** [|kg|] - cumulative mass of the methane in oil
      phase leaked to the aquifer

    * **mass_CO2_aquifer** [|kg|] - cumulative mass of liquid |CO2| leakage to
      aquifer

    * **mass_CO2_gas_aquifer** [|kg|] - cumulative mass of |CO2| gas leaked to the
      aquifer

    """
    def __init__(self, name, parent, **kwargs):
        """
        Constructor method of HydrocarbonLeakage class

        :param name: name of component model
        :type name: [str]

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel [object]

        """
        super().__init__(
            name, parent, model=self.simulation_model)

        self.time_array = parent.time_array

        # Add type attribute
        self.class_type = 'HydrocarbonLeakage'

        # Define default parameters
        self.add_default_par('reservoirDepth', value=2000)
        self.add_default_par('NTG', value=0.6)
        self.add_default_par('logResPerm', value=-13.5)
        self.add_default_par('reservoirPressureMult', value=1.1)
        self.add_default_par('logWellPerm', value=-13.0)
        self.add_default_par('avgWaterSaturation', value=0.5)
        self.add_default_par('FCO2', value=0.55)
        self.add_default_par('FC1', value=0.05)
        self.add_default_par('FC4', value=0.05)
        self.add_default_par('FC7Plus', value=0.35)

        # Define dictionary of parameter boundaries
        self.pars_bounds['reservoirDepth'] = [915.6, 2742.2]
        self.pars_bounds['NTG'] = [0.4, 1.0]
        self.pars_bounds['logResPerm'] = [-14.00, -13.00]
        self.pars_bounds['reservoirPressureMult'] = [1.0, 1.2]
        self.pars_bounds['logWellPerm'] = [-15.61, -12.01]
        self.pars_bounds['avgWaterSaturation'] = [0.471, 0.792]
        self.pars_bounds['FCO2'] = [0.432, 0.693]
        self.pars_bounds['FC1'] = [0.010, 0.113]
        self.pars_bounds['FC4'] = [0.010, 0.111]
        self.pars_bounds['FC7Plus'] = [0.123, 0.500]

        # The model only needs to be run once, as output for all times can be
        # obtained at once.
        self.run_frequency = 1
        self.default_run_frequency = 1

        # Initiate solution object
        self.sol = hydrcarbrom.Solution()

        debug_msg = 'HydrocarbonLeakage component created with name {}'.format(name)
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
                if ((val < self.pars_bounds[key][0])or
                        (val > self.pars_bounds[key][1])):
                    warn_msg = ''.join([
                        'Parameter {} of HydrocarbonLeakage component {} with ',
                        'value {} is out of boundaries.']).format(key, self.name, val)
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of HydrocarbonLeakage component {}.'
                    ]).format(key, self.name)
                logging.warning(warn_msg)

        # Make sure that FCO2, FC1, FC4, and FC7Plus sum to 1.
        FCO2_val = p['FCO2']
        FC1_val = p['FC1']
        FC4_val = p['FC4']
        FC7Plus_val = p['FC7Plus']

        if not np.abs((FCO2_val + FC1_val + FC4_val + FC7Plus_val) - 1) <= 0.0001:
            warn_msg = ''.join([
                'The sum of the HydrocarbonLeakage component parameters FCO2, ',
                'FC1, FC4, and FC7Plus has to be sufficiently close to 1 ',
                '(from 0.9999 to 1.0001). The values provided to the ',
                'HydrocarbonLeakage component {} do not sum to 1. ',
                'Check your input.']).format(self.name)
            logging.warning(warn_msg)

    def connect_with_system(self, component_data, name2obj_dict, *args, **kwargs):
        """
        Code to add hydrocarbon leakage component to system model for control file interface.

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
            new_comp_outputs = []
            for obs_nm in comp_outputs:
                new_comp_outputs.append(obs_nm)
                self.add_obs(obs_nm)

            component_data['Outputs'] = new_comp_outputs

        # Determine number of shale layers in the stratigraphy
        strata = name2obj_dict['strata']
        if 'numberOfShaleLayers' in strata.deterministic_pars:
            num_shale_layers = strata.deterministic_pars['numberOfShaleLayers'].value
        elif 'numberOfShaleLayers' in strata.default_pars:
            num_shale_layers = strata.default_pars['numberOfShaleLayers'].value
        else:
            num_shale_layers = 3

        # Consider an option when reservoirDepth is already added
        # as referring to spatially varying stratigraphy
        if 'reservoirDepth' not in self.pars and \
                'reservoirDepth' not in self.deterministic_pars:
            res_depth_expr = ' + '.join(['{}.shale{}Thickness'.format(
                strata.name, ind) for ind in range(1, num_shale_layers+1)])+' + '+\
                    ' + '.join(['{}.aquifer{}Thickness'.format(
                        strata.name, ind) for ind in range(1, num_shale_layers)])

            # Depth to the top of reservoir (usually)
            self.add_composite_par('reservoirDepth', res_depth_expr)

    def simulation_model(self, p, **kwargs):
        """
        Return leaked masses (kg) for gas and oil phases of CO2 and CH4, as
        well as total gas and total oil.

        :param p: input parameters of the HydrocarbonLeakage ROM
        :type p: dict

        :param time_point: time point in days at which the leakage rates are
            to be calculated; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :returns: out - dictionary of observations of leaked masses of CH4 and
            CO2 gas and oil (liquid) as well as total gas and total liquid.
            keys: ['mass_CO2_aquifer', 'mass_CO2_gas_aquifer',
                   'mass_methane_oil_aquifer', 'mass_methane_gas_aquifer',
                   'mass_oil_aquifer', 'mass_gas_aquifer']
        """
        # Assign default values
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        actual_p.update(p)

        output_names = hydrcarbrom.OUTPUT_NAMES
        output_check = []

        # If one of the outputs is not being used, output_check is used to skip
        # the ROM for that output.
        for label in output_names:
            if label + '_0' in self.obs:
                output_check.append(True)
            else:
                output_check.append(False)

        input_array = hydrcarbrom.set_up_input(self.time_array, actual_p)

        self.sol.find(input_array, output_check)

        output_names = hydrcarbrom.OUTPUT_NAMES

        out = {}

        for count, label in enumerate(output_names):
            for time_ind in range(len(self.time_array)):
                out['{}_{}'.format(label, time_ind)] = self.sol.outputs[
                    time_ind, count]

        # Return dictionary of outputs
        return out


if __name__ == "__main__":
    from openiam import SystemModel
    import matplotlib.pyplot as plt

    __spec__ = None

    # Setup logging level
    logging.basicConfig(level=logging.WARNING)

    # Define keyword arguments of the system model
    num_years = 410
    t0 = 0.0  # initial time point
    time_array = 365.25*np.arange(t0, t0+num_years+10, 10)

    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add hydrocarbon leakage component
    hcl_comp = sm.add_component_model_object(HydrocarbonLeakage(
        name='hcl', parent=sm))

    # Add parameters of hydrocarbon component model
    hcl_comp.add_par('reservoirDepth', value=1046.88)
    hcl_comp.add_par('NTG', value=0.444)
    hcl_comp.add_par('logResPerm', value=-13.23)
    hcl_comp.add_par('reservoirPressureMult', value=1.13)
    hcl_comp.add_par('logWellPerm', value=-12.44)
    hcl_comp.add_par('avgWaterSaturation', value=0.719)
    hcl_comp.add_par('FCO2', value=0.647)
    hcl_comp.add_par('FC1', value=0.062)
    hcl_comp.add_par('FC4', value=0.082)
    hcl_comp.add_par('FC7Plus', value=0.209)

    # Add observations (output) of the hydrocarbon leakage component
    hcl_observations = ['mass_CO2_aquifer', 'mass_CO2_gas_aquifer',
                        'mass_methane_oil_aquifer', 'mass_methane_gas_aquifer',
                        'mass_oil_aquifer', 'mass_gas_aquifer']
    for obs_nm in hcl_observations:
        hcl_comp.add_obs(obs_nm)

    # Run forward simulation
    sm.forward()

    # Collect outputs into a dictionary
    outputs = {obs_nm: sm.collect_observations_as_time_series(hcl_comp, obs_nm)
               for obs_nm in hcl_observations}

    # Setup plot parameters
    label_size = 13
    font_size = 16
    ticks_size = 12
    line_width = 3

    # Plot results
    # First group of plots
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for ind, obs_nm in enumerate(['mass_CO2_aquifer', 'mass_CO2_gas_aquifer']):
        ax[ind].plot(time_array/365.25, outputs[obs_nm]/1.0e+6, '-',
                     color='steelblue', linewidth=line_width)
        ax[ind].set_xlabel('Time, years', fontsize=label_size)
    ax[0].set_ylabel(r'Mass of liquid CO$_2$ in aquifer, [Kt]', fontsize=label_size)
    ax[1].set_ylabel(r'Mass of gas CO$_2$ in aquifer, [Kt]', fontsize=label_size)
    fig.suptitle('Results of simulation', fontsize=font_size)

    plt.show()

    # Second group of plots
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for ind, obs_nm in enumerate(['mass_methane_oil_aquifer', 'mass_methane_gas_aquifer']):
        ax[ind].plot(time_array/365.25, outputs[obs_nm]/1.0e+6, '-',
                     color='steelblue', linewidth=line_width)
        ax[ind].set_xlabel('Time, years', fontsize=label_size)
    ax[0].set_ylabel('Mass of methane in oil phase in aquifer, [Kt]', fontsize=label_size)
    ax[1].set_ylabel('Mass of methane in aquifer, [Kt]', fontsize=label_size)
    fig.suptitle('Results of simulation', fontsize=font_size)

    plt.show()

    # Third group of plots
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for ind, obs_nm in enumerate(['mass_oil_aquifer', 'mass_gas_aquifer']):
        ax[ind].plot(time_array/365.25, outputs[obs_nm]/1.0e+6, '-',
                     color='steelblue', linewidth=line_width)
        ax[ind].set_xlabel('Time, years', fontsize=label_size)
    ax[0].set_ylabel('Mass of oil leaked in aquifer, [Kt]', fontsize=label_size)
    ax[1].set_ylabel('Mass of gas leaked in aquifer, [Kt]', fontsize=label_size)
    fig.suptitle('Results of simulation', fontsize=font_size)

    plt.show()
