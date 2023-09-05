# -*- coding: utf-8 -*-
# Created on Jul 16, 2018
# @author: Kayyum Mansoor
# mansoor1@llnl.gov
import os
import sys
import logging

import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from openiam import SystemModel, ComponentModel, IAM_DIR
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

from openiam.cfi.commons import process_parameters, process_dynamic_inputs

try:
    import components.aquifer.deep_alluvium.deep_alluvium_aquifer_ROM as daarom
except ImportError as err:
    print('\nERROR: Unable to load ROM for Deep Alluvium Aquifer component\n')
    sys.exit()


class DeepAlluviumAquifer(ComponentModel):
    """
    The Deep Alluvium Aquifer component model is a reduced order model which can be
    used to predict the changes in diluted groundwater chemistry if |CO2| and brine
    were to leak into a deep alluvium aquifer similar to the one located below the
    Kimberlina site, in the Southern San Joaquin Valley, California. The protocol
    allows uncertainty and variability in aquifer heterogeneity, fluid transport,
    and potential |CO2| and brine leakage rates from abandoned or damaged oil and gas
    wells to be collectively evaluated to assess potential changes in groundwater
    pH, total dissolved solids (TDS), and changes in the aquifer pressure resulting
    from leakage.

    Although the Deep Alluvium Aquifer model was developed using site-specific data
    from the LLNL's Kimberlina Model (version 1.2), the model accepts aquifer
    characteristics as variable inputs and, therefore, may have broader
    applicability. Careful consideration should be given to the hydrogeochemical
    character of the aquifer before using this model at a new site.

    Model was created using the ``py-earth`` Python package :cite:`pyearth2013`.
    Simulation data used to build this model was created by Mansoor et al. :cite:`DAA2018`.
    In the NRAP-Open-IAM control file, the type name for the Deep Alluvium Aquifer
    component is ``DeepAlluviumAquifer``.

    Component model input definitions:

    * **logK_sand1** [|log10| |m^2|] (-12.92 to -10.92) - permeability of layer 1
      at depth between 10 and 546 |m| (default: -11.92)

    * **logK_sand2** [|log10| |m^2|] (-12.72 to -10.72) - permeability of layer 2
      at depth between 546 and 1225 |m| (default: -11.72)

    * **logK_sand3** [|log10| |m^2|] (-12.7 to -10.7) - permeability of layer 3
      at depth between 1225 and 1411 |m| (default: -11.70)

    * **logK_caprock** [|log10| |m^2|] (-16.699 to -14.699) - permeability of caprock
      at depth between 0 and 10 |m| (default: -15.70)

    * **correlationLengthX** [|m|] (200 to 2000) - correlation length in x-direction
      (default: 1098.99)

    * **correlationLengthZ** [|m|] (10 to 150) - correlation length in z-direction
      (default: 79.81)

    * **sandFraction** [-] (0.7 to 0.9) - sand volume fraction (default: 0.8)

    * **groundwater_gradient** [-] (0.001000 to 0.001667) - regional groundwater
      gradient (dh/dx=change in hydraulic head/distance) (default: 0.001333)

    * **leak_depth** [|m|] (424.36 to 1341.48) - depth of leakage interval
      (default: 885.51).

    Component model dynamic inputs:

    * **brine_rate** [|kg/s|] (0 to 0.017) - brine rate (default: 0.0003)

    * **brine_mass** [|kg|] (238.14419 to 8689604.29) - cumulative brine mass
      (default: 84722.74=10**4.928)

    * **co2_rate** [|kg/s|] (0 to 0.385) - |CO2| rate (default: 0.045)

    * **co2_mass** [|kg|] (1.002 to 1.621e+9) - cumulative |CO2| mass
      (default: 1.636e+7=10**7.214).

    Observations from the Deep Alluvium Aquifer component are:

    * **TDS_volume** [|m^3|] - volume of plume above baseline TDS change in |millig/L|
      (change in TDS > 100 |millig/L|)

    * **TDS_dx** [|m|] - length of plume above baseline TDS change in |millig/L|
      (change in TDS > 100 |millig/L|)

    * **TDS_dy** [|m|] - width of plume above baseline TDS change in |millig/L|
      (change in TDS > 100 |millig/L|)

    * **TDS_dz** [|m|] - height of plume above baseline TDS change in |millig/L|
      (change in TDS > 100 |millig/L|)

    * **Pressure_volume** [|m^3|] - volume of plume above baseline pressure change
      in |Pa| (change in pressure > 500 |Pa|)

    * **Pressure_dx** [|m|] - length of plume above baseline pressure change in |Pa|
      (change in pressure > 500 |Pa|)

    * **Pressure_dy** [|m|] - width of plume above baseline pressure change in |Pa|
      (change in pressure > 500 |Pa|)

    * **Pressure_dz** [|m|] - height of plume above baseline pressure change in |Pa|
      (change in pressure > 500 |Pa|)

    * **pH_volume** [|m^3|] - volume of plume below pH threshold (pH < 6.75)

    * **pH_dx** [|m|] - length of plume below pH threshold (pH < 6.75)

    * **pH_dy** [|m|] - width of plume below pH threshold (pH < 6.75)

    * **pH_dz** [|m|] - height of plume below pH threshold (pH < 6.75).

    """
    def __init__(self, name, parent):
        """
        Constructor method of DeepAlluviumAquifer class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :returns: DeepAlluviumAquifer class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25}   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'DeepAlluviumAquifer'

        # Initiate solution object
        self.sol = daarom.Solution()

        self.add_default_par('logK_sand1', value=-11.92098495)
        self.add_default_par('logK_sand2', value=-11.7198002)
        self.add_default_par('logK_sand3', value=-11.70137252)
        self.add_default_par('logK_caprock', value=-15.69758676)
        self.add_default_par('correlationLengthX', value=1098.994284)
        self.add_default_par('correlationLengthZ', value=79.8062668)
        self.add_default_par('sandFraction', value=0.800121364)
        self.add_default_par('groundwater_gradient', value=0.001333374)
        self.add_default_par('leak_depth', value=885.5060281)

        # Define dictionary of boundaries
        self.pars_bounds = dict()

        self.pars_bounds['logK_sand1'] = [-12.92, -10.92]
        self.pars_bounds['logK_sand2'] = [-12.72, -10.72]
        self.pars_bounds['logK_sand3'] = [-12.7, -10.7]
        self.pars_bounds['logK_caprock'] = [-16.699, -14.699]
        self.pars_bounds['correlationLengthX'] = [200, 2000]
        self.pars_bounds['correlationLengthZ'] = [10, 150]
        self.pars_bounds['sandFraction'] = [0.7, 0.9]
        self.pars_bounds['groundwater_gradient'] = [0.001000, 0.001667]
        self.pars_bounds['leak_depth'] = [424.36, 1341.48]

        # Define dictionary of temporal data limits
        self.temp_data_bounds = dict()
        self.temp_data_bounds['brine_rate'] = ['brine rate', 0.00000, 0.01700]
        self.temp_data_bounds['brine_mass'] = ['brine mass', 10**2.37684, 10**6.93934]
        self.temp_data_bounds['co2_rate'] = ['CO2 rate', 0.00000, 0.38498]
        self.temp_data_bounds['co2_mass'] = ['CO2 mass', 10**0.00107, 10**9.20977]
        self.temp_data_bounds['time'] = ['simulation time', 0, 200]

        # Define observations names
        self.out_labels = [
            'TDS_volume', 'TDS_dx', 'TDS_dy', 'TDS_dz',
            'Pressure_volume', 'Pressure_dx', 'Pressure_dy', 'Pressure_dz',
            'pH_volume', 'pH_dx', 'pH_dy', 'pH_dz']

        debug_msg = 'DeepAlluviumAquifer component created with name {}'.format(name)
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
                if ((val < self.pars_bounds[key][0]) or (
                        val > self.pars_bounds[key][1])):
                    warn_msg = ''.join([
                        'Parameter {} of DeepAlluviumAquifer component {} ',
                        'is out of bounds.']).format(key, self.name)
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of a DeepAlluviumAquifer component {}.']).format(key, self.name)
                logging.warning(warn_msg)

    def check_temporal_inputs(self, time, temp_inputs):
        """
        Check whether temporal data fall within specified boundaries.

        :param temp_inputs: temporal input data of component model
        :type temp_inputs: dict
        """
        for key, val in temp_inputs.items():
            if (val < self.temp_data_bounds[key][1]) or (
                    val > self.temp_data_bounds[key][2]):
                warn_msg = ''.join([
                    'Temporal input {} (value: {}) of DeepAlluviumAquifer component {} ',
                    'is outside the model range [{}, {}] at time t = {}']).format(
                        self.temp_data_bounds[key][0].lower(), val, self.name,
                        self.temp_data_bounds[key][1],
                        self.temp_data_bounds[key][2], time)
                logging.warning(warn_msg)

    # deep alluvium aquifer model function
    def simulation_model(self, p, co2_rate=0.045, brine_rate=0.0003,
                         co2_mass=10**7.214, brine_mass=10**4.928,
                         time_point=365.25):
        """
        Return volume of impacted aquifer based on several metrics.

        :param p: input parameters of deep alluvium aquifer model
        :type p: dict

        :param co2_rate: CO2 leakage rate, kg/s
        :type co2_rate: float

        :param brine_mass: Brine leakage rate, kg/s
        :type brine_mass: float

        :param co2_mass: Cumulative CO2 mass leaked, kg
        :type co2_mass: float

        :param brine_mass: Cumulative brine mass leaked, kg
        :type brine_mass: float

        :param time_point: time point in days at which the leakage rates are
            to be calculated; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :returns: vol - dictionary of observations of deep alluvium aquifer
            impacted volumes model;
            keys: ['TDS', 'pH', 'As', 'Pb', 'Cd', 'Ba', 'Benzene', 'Naphthalene', 'Phenol']
        """
        # Check assign default values
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        actual_p.update(p)

        # Check whether co2, brine, rates and mass inputs satisfy the model requirements
        # Conditional statement used to supress excessive error messsages
        if co2_mass <= 0:
            co2_mass = 1e-2
        if brine_mass <= 0:
            brine_mass = 1e-2

        if ((time_point/365.25 >= 1) and (
                np.log10(co2_mass) >= self.temp_data_bounds['co2_mass'][1])) or (
                    (time_point/365.25 >= 1) and (
                        np.log10(brine_mass) >= self.temp_data_bounds['brine_mass'][1])):
            self.check_temporal_inputs(
                time_point, dict(list(zip(
                    ['co2_rate', 'brine_rate', 'co2_mass', 'brine_mass', 'time'],
                    [co2_rate, brine_rate, co2_mass, brine_mass, time_point/365.25]))))

        # Apply nominal value for brine/co2
        p['co2_rate'] = co2_rate
        p['brine_rate'] = brine_rate
        p['co2_mass'] = co2_mass
        p['brine_mass'] = brine_mass

        # Update default values of parameters with the provided ones
        actual_p.update(p)

        inputArray = np.array([
            time_point/365.25, actual_p['brine_rate'], np.log10(actual_p['brine_mass']),
            actual_p['co2_rate'], np.log10(actual_p['co2_mass']),
            actual_p['logK_sand1'], actual_p['logK_sand2'],
            actual_p['logK_sand3'], actual_p['logK_caprock'],
            actual_p['correlationLengthX'], actual_p['correlationLengthZ'],
            actual_p['sandFraction'], actual_p['groundwater_gradient'],
            actual_p['leak_depth']])

        if ((time_point > 0) and (
                np.log10(co2_mass) >= self.temp_data_bounds['co2_mass'][1])) or (
                    (time_point > 0) and (
                        np.log10(brine_mass) >= self.temp_data_bounds['brine_mass'][1])):

            self.sol.find(inputArray)
            out = dict(list(zip(self.out_labels, self.sol.Outputs)))
        else:
            out = dict(list(zip(self.out_labels, np.zeros(len(self.out_labels)))))
        return out

    # Attributes for system connections
    adapters = ['RateToMass']
    needsXY = False

    def connect_with_system(self, component_data, name2obj_dict,
                            locations, system_adapters, **kwargs):
        """
        Code to add deep alluvium aquifer to system model for control file interface.

        :param component_data: Dictionary of component data including 'Connection'
        :type component_data: dict

        :param name2obj_dict: Dictionary to map component names to objects
        :type name2obj: dict

        :param locations: Dictionary with keys being the names of wellbore components
            added to the system and values being a dictionary with 3 keys:
            'number' (number of locations), 'coordx' (x-coordinates of locations)
            and 'coordy' (y-coordinates of locations).
        :type locations: dict

        :param system_adapters: Dictionary with keys being the names of adapters
            added to the system and values being a dictionary with 3 keys:
            'Type' (type of adapter), 'Connection' (name of wellbore component
            to which the given adapter is connected) and 'AquiferName' (name of
            aquifer for which the leakage rates of brine and CO2 from
            the connected wellbore component are transformed into the masses)
        :type system_adapters: dict

        :returns: None
        """
        # Process parameters data if any
        process_parameters(self, component_data, name2obj_dict)

        # Process dynamic inputs if any
        process_dynamic_inputs(self, component_data)

        # Make model connections
        if 'Connection' in component_data:
            connection = None
            try:
                connection = name2obj_dict[component_data['Connection']]
            except KeyError:
                pass

            if 'AquiferName' in component_data:
                aq_name = component_data['AquiferName']
            else:
                aq_name = 'aquifer1'

            for ad_nm in system_adapters:
                if (system_adapters[ad_nm]['Connection'] == component_data['Connection']) and (
                        system_adapters[ad_nm]['AquiferName'] == aq_name):
                    aname = ad_nm

            adapter = name2obj_dict[aname]

            connection.add_obs_to_be_linked(
                'CO2_{aquifer_name}'.format(aquifer_name=aq_name))
            adapter.add_obs_to_be_linked(
                'mass_CO2_{aquifer_name}'.format(aquifer_name=aq_name))
            connection.add_obs_to_be_linked(
                'brine_{aquifer_name}'.format(aquifer_name=aq_name))
            adapter.add_obs_to_be_linked(
                'mass_brine_{aquifer_name}'.format(aquifer_name=aq_name))

            self.add_kwarg_linked_to_obs(
                'co2_rate', connection.linkobs['CO2_{aquifer_name}'.format(
                    aquifer_name=aq_name)])
            self.add_kwarg_linked_to_obs(
                'co2_mass', adapter.linkobs['mass_CO2_{aquifer_name}'.format(
                    aquifer_name=aq_name)])
            self.add_kwarg_linked_to_obs(
                'brine_rate', connection.linkobs['brine_{aquifer_name}'.format(
                    aquifer_name=aq_name)])
            self.add_kwarg_linked_to_obs(
                'brine_mass', adapter.linkobs['mass_brine_{aquifer_name}'.format(
                    aquifer_name=aq_name)])


if __name__ == "__main__":
    # Create system model
    time_array = 365.25*np.arange(0.0, 2.0)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add deep alluvium aquifer model object and define parameters
    daa = sm.add_component_model_object(DeepAlluviumAquifer(name='daa', parent=sm))

    daa.add_par('logK_sand1', value=-11.92098495)
    daa.add_par('logK_sand2', value=-11.7198002)
    daa.add_par('logK_sand3', value=-11.70137252)
    daa.add_par('logK_caprock', value=-15.69758676)
    daa.add_par('correlationLengthX', value=1098.994284)
    daa.add_par('correlationLengthZ', value=79.8062668)
    daa.add_par('sandFraction', value=0.800121364)
    daa.add_par('groundwater_gradient', value=0.001333374)
    daa.add_par('leak_depth', value=885.5060281)

    daa.model_kwargs['brine_rate'] = 3.21903E-05     # kg/s
    daa.model_kwargs['brine_mass'] = 10**4.71081307  # kg
    daa.model_kwargs['co2_rate'] = 0.060985038       # kg/s
    daa.model_kwargs['co2_mass'] = 10**6.737803184   # kg

    daa.add_obs('TDS_volume')
    daa.add_obs('TDS_dx')
    daa.add_obs('TDS_dy')
    daa.add_obs('TDS_dz')

    daa.add_obs('Pressure_volume')
    daa.add_obs('Pressure_dx')
    daa.add_obs('Pressure_dy')
    daa.add_obs('Pressure_dz')

    daa.add_obs('pH_volume')
    daa.add_obs('pH_dx')
    daa.add_obs('pH_dy')
    daa.add_obs('pH_dz')

    # Run the system model
    sm.forward()

    # Print the observations
    # Since the observations at the particular time points are different variables,
    # method collect_observations_as_time_series creates lists of
    # values of observations belonging to a given component (e.g. ca) and having the same
    # common name (e.g. 'pH', 'TDS', etc) but differing in indices.
    # More details are given in the docstring and documentation to the method
    # collect_observations_as_time_series of SystemModel class.
    print('{:15}'.format('TDS_volume'),
          sm.collect_observations_as_time_series(daa, 'TDS_volume'))
    print('{:15}'.format('TDS_dx'),
          sm.collect_observations_as_time_series(daa, 'TDS_dx'))
    print('{:15}'.format('TDS_dy'),
          sm.collect_observations_as_time_series(daa, 'TDS_dy'))
    print('{:15}'.format('TDS_dz'),
          sm.collect_observations_as_time_series(daa, 'TDS_dz'))

    print('{:15}'.format('Pressure_volume'),
          sm.collect_observations_as_time_series(daa, 'Pressure_volume'))
    print('{:15}'.format('Pressure_dx'),
          sm.collect_observations_as_time_series(daa, 'Pressure_dx'))
    print('{:15}'.format('Pressure_dy'),
          sm.collect_observations_as_time_series(daa, 'Pressure_dy'))
    print('{:15}'.format('Pressure_dz'),
          sm.collect_observations_as_time_series(daa, 'Pressure_dz'))

    print('{:15}'.format('pH_volume'),
          sm.collect_observations_as_time_series(daa, 'pH_volume'))
    print('{:15}'.format('pH_dx'),
          sm.collect_observations_as_time_series(daa, 'pH_dx'))
    print('{:15}'.format('pH_dy'),
          sm.collect_observations_as_time_series(daa, 'pH_dy'))
    print('{:15}'.format('pH_dz'),
          sm.collect_observations_as_time_series(daa, 'pH_dz'))

#    # Expected output
#    TDS_volume      [      0.         2436214.13519926]
#    TDS_dx          [  0.         121.07792897]
#    TDS_dy          [  0.         139.62493095]
#    TDS_dz          [  0.         257.07734713]
#    Pressure_volume [      0.         2945213.33234456]
#    Pressure_dx     [  0.         176.71076959]
#    Pressure_dy     [  0.         152.39886469]
#    Pressure_dz     [  0.        260.9642596]
#    pH_volume       [      0.         2656675.37737679]
#    pH_dx           [ 0.         95.07087112]
#    pH_dy           [ 0.         88.08772333]
#    pH_dz           [  0.         249.01291277]
