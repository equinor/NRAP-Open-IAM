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
    print('Unable to load IAM class module: {}'.format(err))

from openiam.cfi.commons import process_parameters, process_dynamic_inputs


class DeepAlluviumAquiferML(ComponentModel):
    """
    The Deep Alluvium Aquifer ML component model is a reduced order model which can be
    used to predict the changes in diluted groundwater chemistry if |CO2| and brine
    were to leak into a deep alluvium aquifer similar to the one located below the
    Kimberlina site, in the Southern San Joaquin Valley, California. The protocol
    allows uncertainty and variability in aquifer heterogeneity, fluid transport,
    and potential |CO2| and brine leakage rates from abandoned or damaged oil and gas
    wells to be collectively evaluated to assess potential changes in groundwater
    pH, total dissolved solids (TDS), and changes in the aquifer pressure resulting
    from leakage.

    In the NRAP-Open-IAM control file, the type name for the Deep Alluvium
    Aquifer ML component is ``DeepAlluviumAquiferML``. The description
    of the possible component's parameters are provided below:

    * **logK_sand1** * [|log10| |m^2|] (-12.92 to -10.92) - permeability
      of layer 1 at the depth between 10 and 546 m (default: -11.91)

    * **logK_sand2** * [|log10| |m^2|] (-12.72 to -10.72) - permeability
      of layer 2 at the depth between 546 and 1225 m (default: -11.71)

    * **logK_sand3** * [|log10| |m^2|] (-12.7 to -10.7) - permeability
      of layer 3 at the depth between 1225 and 1411 m (default: -11.69)

    * **logK_caprock** * [|log10| |m^2|] (-16.7 to -14.7) - caprock permeability
      at the depth between 0 and 5 m (default: -15.7)

    * **correlationLengthX** * [|m|] (200 to 2000) - correlation length in x-direction
      (default: 1098.235)

    * **correlationLengthZ** * [|m|] (10 to 150) - correlation length in z-direction
      (default: 79.827)

    * **sandFraction** * [-] (0.7 to 0.9) - sand volume fraction (default: 0.8)

    * **groundwater_gradient** * [-] (0.001 to 0.0017) - regional groundwater gradient
      (dh/dx=change in hydraulic head over distance) (default: 0.0013)

    * **leak_depth** * [-] (424.4 to 1341.5) - depth of leakage interval (default: 883.3).

    Component model dynamic inputs:

    * **simtime** * [years] (0 to 200) - simulation time (default: 20)

    * **brine_rate** * [|kg/s|] (3.4951e-7 to 0.017) - brine rate (default: 0.0003)

    * **brine_mass** * [|kg|] (23.727 to 8.689e+6) - cumulative brine mass
      (default: 8.472e+4)

    * **co2_rate** * [|kg/s|] (5.0e-11 to 0.385) - |CO2| rate (default: 0.045)

    * **co2_mass** * [|kg|] (1.577e-3 to 1.621e+9) - cumulative |CO2| mass
      (default: 1.636e+7)

    Observations from the Deep Alluvium Aquifer component are:

    * **TDS_volume** [|m^3|] - volume of plume above baseline TDS change in mg/L
      (change in TDS > 100 mg/L)

    * **TDS_dx** [|m|] - length of plume above baseline TDS change in mg/L
      (change in TDS > 100 mg/L)

    * **TDS_dy** [|m|] - width of plume above baseline TDS change in mg/L
      (change in TDS > 100 mg/L)

    * **TDS_dz** [|m|] - height of plume above baseline TDS change in mg/L
      (change in TDS > 100 mg/L)

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
        Constructor method of DeepAlluviumAquiferML class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :returns: DeepAlluviumAquiferML class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25}  # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'DeepAlluviumAquiferML'

        # Default parameters here
        self.add_default_par('logK_sand1', value=-11.3845)
        self.add_default_par('logK_sand2', value=-11.9252)
        self.add_default_par('logK_sand3', value=-10.8862)
        self.add_default_par('logK_caprock', value=-14.731)
        self.add_default_par('correlationLengthX', value=520.721)
        self.add_default_par('correlationLengthZ', value=112.442)
        self.add_default_par('sandFraction', value=0.743)
        self.add_default_par('groundwater_gradient', value=0.00105408)
        self.add_default_par('leak_depth', value=715.99)

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
        self.temp_data_bounds['brine_rate'] = [3.4951E-07, 0.017002313]
        self.temp_data_bounds['brine_mass'] = [23.72739171, 8696409.865] # 10**1.37525, 10**6.93934
        self.temp_data_bounds['co2_rate'] = [5.000E-11, 0.384984283]
        self.temp_data_bounds['co2_mass'] = [0.001577866, 1620951423] # 10**-2.80193, 10**9.20977
        self.temp_data_bounds['time'] = [1, 200]

        # Define output dictionary labels
        self.output_labels = [
            'TDS_volume', 'TDS_dx', 'TDS_dy', 'TDS_dz',
            'Pressure_volume', 'Pressure_dx', 'Pressure_dy', 'Pressure_dz',
            'pH_volume', 'pH_dx', 'pH_dy', 'pH_dz']

        # Check if another deep alluvium aquifer ml component is already part of the system model
        orig_daac_cmpnt = self._parent.ml_models_components.get(
            'DeepAlluviumAquiferML', None)

        try:
            import components.aquifer.deep_alluvium_ml.deep_alluvium_aquifer_ml_rom as daamlrom
        except ImportError:
            print('\nERROR: Unable to load ROM for Deep Alluvium Aquifer (ML) component\n')
            sys.exit()

        if orig_daac_cmpnt is None:
            # Register the component as the one using ml models
            # The component will be registered only if no other deep alluvium aquifer ml
            # components were added to the same system model
            self._parent.register_ml_model_component(self, 'DeepAlluviumAquiferML')

            # Initiate solution object
            self.sol = daamlrom.Solution()
        else:
            # Get solution class object from existing component of the same type
            self.sol = orig_daac_cmpnt.sol

        debug_msg = 'DeepAlluviumAquiferML created with name {name}'.format(name=name)
        logging.debug(debug_msg)

    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        debug_msg = 'Input parameters of {name} component are {p}'.format(name=self.name, p=p)
        logging.debug(debug_msg)

        for key, val in p.items():
            if key in self.pars_bounds:
                if ((val < self.pars_bounds[key][0])or
                        (val > self.pars_bounds[key][1])):
                    warn_msg = 'Parameter {} is out of bounds.'.format(key)
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} not recognized as ',
                    'a Deep Alluvium Aquifer input parameter.']).format(key)
                logging.warning(warn_msg)

    def check_temporal_inputs(self, time, temp_inputs):
        """
        Check whether temporal data fall within specified boundaries.

        :param temp_inputs: temporal input data of component model
        :type temp_inputs: dict
        """
        debug_msg = 'Temporal Inputs at time {0} are {1}'.format(time, temp_inputs)
        logging.debug(debug_msg)
        for key, val in temp_inputs.items():
            if ((val < self.temp_data_bounds[key][0]) or
                    (val > self.temp_data_bounds[key][1])):
                warn_msg = ''.join([
                    '{0} {1} is outside the model range {2} to {3} ',
                    'at time t = {4}.']).format(key, val, self.temp_data_bounds[key][0],
                                                self.temp_data_bounds[key][1], time)
                logging.warning(warn_msg)

    def simulation_model(self, p, co2_rate=0.045, brine_rate=0.0003,
                         co2_mass=10**7.214, brine_mass=10**4.928,
                         time_point=365.25):
        """
        Return volume of impacted aquifer based on several metrics.

        :param p: input parameters of deep alluvium aquifer model
        :type p: dict

        :param co2_rate: CO2 leakage rate, kg/s
        :type co2_rate: float

        :param brine_rate: Brine leakage rate, kg/s
        :type brine_rate: float

        :param co2_mass: Cumulative CO2 mass leaked, kg
        :type co2_mass: float

        :param brine_mass: Cumulative brine mass leaked, kg
        :type brine_mass: float

        :param time_point: time point in days at which the leakage rates are
            to be calculated; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param kwargs: additional keyword arguments for the model
        :type kwargs: dict (if provided)

        :returns: vol - dictionary of observations of deep alluvium aquifer
            impacted volumes; keys: ['TDS', 'Pressure', 'pH']
        """
        # Check assign defailt values
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        actual_p.update(p)

        # Check whether CO2, brine, rates and mass inputs satisfy the model requirements
        # Conditional statement used to supress excessive error messsages
        if co2_mass <= 0:
            co2_mass = self.temp_data_bounds['co2_mass'][0]
        if brine_mass <= 0:
            brine_mass = self.temp_data_bounds['brine_mass'][0]

        if (time_point >= 365.25) and ((
                co2_mass >= self.temp_data_bounds['co2_mass'][0]) or (
                    brine_mass >= self.temp_data_bounds['brine_mass'][0])):
            self.check_temporal_inputs(
                time_point, dict(list(zip(
                    ['co2_rate', 'brine_rate', 'co2_mass', 'brine_mass', 'time'],
                    [co2_rate, brine_rate, co2_mass, brine_mass, time_point/365.25]))))

        # Apply nominal value for brine/CO2
        if (co2_rate is not None) and (brine_rate is not None):
            p['co2_rate'] = co2_rate
            p['brine_rate'] = brine_rate
            p['co2_mass'] = co2_mass
            p['brine_mass'] = brine_mass
        else:
            error_msg = 'CO2 and/or brine rate arguments of the model method are not provided.'
            logging.error(error_msg)

        # Update default values of parameters with the provided ones
        actual_p.update(p)
        inputArray = np.array([
            time_point/365.25, actual_p['brine_rate'], actual_p['brine_mass'],
            actual_p['co2_rate'], actual_p['co2_mass'],
            actual_p['logK_sand1'], actual_p['logK_sand2'], actual_p['logK_sand3'],
            actual_p['logK_caprock'], actual_p['correlationLengthX'],
            actual_p['correlationLengthZ'], actual_p['sandFraction'],
            actual_p['groundwater_gradient'], actual_p['leak_depth']])

        labels = ['TDS_volume', 'TDS_dx', 'TDS_dy', 'TDS_dz',
                  'Pressure_volume', 'Pressure_dx', 'Pressure_dy', 'Pressure_dz',
                  'pH_volume', 'pH_dx', 'pH_dy', 'pH_dz']

        if time_point > 0:
            self.sol.find(inputArray)
            out = dict(list(zip(labels, self.sol.Outputs)))
        else:
            out = dict(list(zip(labels, np.zeros(len(labels)))))

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
    time_array = 365.25*np.arange(0.0, 4.0)

    sm_model_kwargs = {'time_array': time_array} # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add deep alluvium aquifer model object and define parameters
    daaml = sm.add_component_model_object(DeepAlluviumAquiferML(name='daaml', parent=sm))
    daaml.add_par('logK_sand1', value=-11.3845)
    daaml.add_par('logK_sand2', value=-11.9252)
    daaml.add_par('logK_sand3', value=-10.8862)
    daaml.add_par('logK_caprock', value=-14.731)
    daaml.add_par('correlationLengthX', value=520.721)
    daaml.add_par('correlationLengthZ', value=112.442)
    daaml.add_par('sandFraction', value=0.743)
    daaml.add_par('groundwater_gradient', value=0.00105408)
    daaml.add_par('leak_depth', value=715.99)

    daaml.model_kwargs['brine_rate'] = 6.24E-05    # kg/s
    daaml.model_kwargs['brine_mass'] = 2646900     # 10**6.422737534  kg
    daaml.model_kwargs['co2_rate'] = 0.150771      # kg/s
    daaml.model_kwargs['co2_mass'] = 171604000.1   # 10**8.234527407 kg

    # Add observations (output) of the deep alluvium aquifer model
    daaml.add_obs('TDS_volume')
    daaml.add_obs('TDS_dx')
    daaml.add_obs('TDS_dy')
    daaml.add_obs('TDS_dz')

    daaml.add_obs('Pressure_volume')
    daaml.add_obs('Pressure_dx')
    daaml.add_obs('Pressure_dy')
    daaml.add_obs('Pressure_dz')

    daaml.add_obs('pH_volume')
    daaml.add_obs('pH_dx')
    daaml.add_obs('pH_dy')
    daaml.add_obs('pH_dz')

    # Run the system model
    sm.forward()

    # Print the observations
    print('{:15}'.format('TDS_volume'),
          sm.collect_observations_as_time_series(daaml, 'TDS_volume'))
    print('{:15}'.format('TDS_dx'),
          sm.collect_observations_as_time_series(daaml, 'TDS_dx'))
    print('{:15}'.format('TDS_dy'),
          sm.collect_observations_as_time_series(daaml, 'TDS_dy'))
    print('{:15}'.format('TDS_dz'),
          sm.collect_observations_as_time_series(daaml, 'TDS_dz'))

    print('{:15}'.format('Pressure_volume'),
          sm.collect_observations_as_time_series(daaml, 'Pressure_volume'))
    print('{:15}'.format('Pressure_dx'),
          sm.collect_observations_as_time_series(daaml, 'Pressure_dx'))
    print('{:15}'.format('Pressure_dy'),
          sm.collect_observations_as_time_series(daaml, 'Pressure_dy'))
    print('{:15}'.format('Pressure_dz'),
          sm.collect_observations_as_time_series(daaml, 'Pressure_dz'))

    print('{:15}'.format('pH_volume'),
          sm.collect_observations_as_time_series(daaml, 'pH_volume'))
    print('{:15}'.format('pH_dx'),
          sm.collect_observations_as_time_series(daaml, 'pH_dx'))
    print('{:15}'.format('pH_dy'),
          sm.collect_observations_as_time_series(daaml, 'pH_dy'))
    print('{:15}'.format('pH_dz'),
          sm.collect_observations_as_time_series(daaml, 'pH_dz'))

# Expected output
# TDS_volume      [0.  1602563.55016136   5438865.18640179   8059917.15962449]
# TDS_dx          [0.  65.69592945        122.73450292       213.45443131]
# TDS_dy          [0.  182.430035         225.89221453       252.38214608]
# TDS_dz          [0.  46.41604274        104.95408265       161.20099329]
# Pressure_volume [0.  14855869.2546651   11719586.45504402  10645598.30293803]
# Pressure_dx     [0.  278.70072389       335.24981865       379.97843937]
# Pressure_dy     [0.  269.11794614       302.68000236       314.19433156]
# Pressure_dz     [0.  348.06364686       370.36887966       389.19799003]
# pH_volume       [0.  26076001.86850136  54044350.52776746  67349593.38404109]
# pH_dx           [0.  519.86851253       656.91480486       664.48539853]
# pH_dy           [0.  362.15100201       523.36799963       582.61274082]
# pH_dz           [0.  842.09423191       1037.53341233      1082.90360577]
