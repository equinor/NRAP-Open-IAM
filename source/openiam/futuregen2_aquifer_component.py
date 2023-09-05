# -*- coding: utf-8 -*-
# Created on Oct 16, 2018
# @author: Diana Bacon
# diana.bacon@pnnl.gov

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
    import components.aquifer.FutureGen2.aquifer.futuregen2_aquifer_ROM as fgarom
except ImportError:
    print('\nERROR: Unable to load ROM for FutureGen2 Aquifer component\n')
    sys.exit()


class FutureGen2Aquifer(ComponentModel):
    """
    The FutureGen 2.0 Aquifer component model is a reduced order model
    that can be used to predict the impact that carbon dioxide (|CO2|) and
    brine leakage from the |CO2| storage reservoir at the FutureGen 2.0 site
    might have on overlying aquifers or monitoring units. The model predicts
    the size of “impact plumes” according to four metrics:
    pH, total dissolved solids (TDS), pressure and dissolved |CO2|.

    The FutureGen 2.0 Aquifer model is a regression model fitted to the results of
    STOMP-CO2E-R multiphase flow and reactive transport simulations of |CO2| and
    brine leakage using the ``py-earth`` Python package (:cite:`pyearth2013`).
    The ``py-earth`` package is a Python implementation of the Multivariate Adaptive
    Regression Splines algorithm (:cite:`Friedman1991`), in the style of ``scikit-learn``
    (:cite:`scikit`), a library of machine-learning methods.

    The aquifer simulations used to train the FutureGen 2.0 Aquifer component model were
    based on modeling done for monitoring program design at the FutureGen 2.0 site
    (:cite:`Vermeul2016`), as well as porosity and permeability values from the ELAN logs
    and core samples taken from the characterization well. Isothermal simulations
    were performed for training the aquifer component model.

    In the NRAP-Open-IAM control file, the type name for the FutureGen 2.0 Aquifer component is
    ``FutureGen2Aquifer``. The description of the possible component's parameters are
    provided below:

    * **aqu_thick** [|m|] (30 to 90) - thickness of unit (default: 33.2);
      *linked to Stratigraphy*

    * **depth** [|m|] (100 to 700) - depth to bottom of unit (default: 590.1);
      *linked to Stratigraphy*

    * **por** [-] (0.02 to 0.2) - porosity of unit (default: 0.118)

    * **log_permh** [|log10| |m^2|] (-14 to -11) - horizontal permeability
      (default: -13.39)

    * **log_aniso** [|log10|] (0 to 3) - anisotropy ratio (default: 0.3)

    * **rel_vol_frac_calcite** [-] (0 to 1) - relative volume fraction of calcite
      in solid phase (default: 0.01).

    Component model dynamic inputs:

    * **brine_rate** [|kg/s|] (0 to 31.622) - brine rate

    * **brine_mass** [|kg|] (0 to 6.985e+10) - cumulative brine mass

    * **co2_rate** [|kg/s|] (0 to 31.622) - |CO2| rate

    * **co2_mass** [|kg|] (0 to 6.985e+10) - cumulative |CO2| mass.

    Observations from the FutureGen 2.0 Aquifer component are:

    * **Pressure_volume** [|m^3|] - volume of plume where relative change
      in pressure > 0.065%

    * **Pressure_dx** [|m|] - length of plume where relative change
      in pressure > 0.065%

    * **Pressure_dy** [|m|] - width of plume where relative change
      in pressure > 0.065%

    * **Pressure_dz** [|m|] - height of plume where relative change
      in pressure > 0.065%

    * **pH_volume** [|m^3|] - volume of plume where absolute change in pH > 0.2

    * **pH_dx** [|m|] - length of plume where absolute change in pH > 0.2

    * **pH_dy** [|m|] - width of plume where absolute change in pH  > 0.2

    * **pH_dz** [|m|] - height of plume where absolute change in pH  > 0.2

    * **TDS_volume** [|m^3|] - volume of plume where relative change in TDS > 10%

    * **TDS_dx** [|m|] - length of plume where relative change in TDS > 10%

    * **TDS_dy** [|m|] - width of plume where relative change in TDS > 10%

    * **TDS_dz** [|m|] - height of plume where relative change in TDS > 10%

    * **Dissolved_CO2_volume** [|m^3|] - volume of plume where relative change
      in dissolved |CO2| concentration > 20%

    * **Dissolved_CO2_dx** [|m|] - length of plume where relative change
      in dissolved |CO2| concentration > 20%

    * **Dissolved_CO2_dy** [|m|] - width of plume where relative change
      in dissolved |CO2| concentration > 20%

    * **Dissolved_CO2_dz** [|m|] - height of plume where relative change
      in dissolved |CO2| concentration > 20%

    """
    def __init__(self, name, parent):
        """
        Constructor method of FutureGen2Aquifer class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :returns: Futuregen2Aquifer class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25}

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'FutureGen2Aquifer'

        # Define default parameters
        self.add_default_par('aqu_thick', value=33.2)
        self.add_default_par('depth', value=590.1)
        self.add_default_par('por', value=0.118)
        self.add_default_par('log_permh', value=-13.39)
        self.add_default_par('log_aniso', value=0.3)
        self.add_default_par('rel_vol_frac_calcite', value=0.01)

        # Define dictionary of parameter boundaries
        self.pars_bounds = dict()
        self.pars_bounds['aqu_thick'] = [30., 90.]
        self.pars_bounds['depth'] = [100., 700.]
        self.pars_bounds['por'] = [0.02, 0.2]
        self.pars_bounds['log_permh'] = [-14., -11.]
        self.pars_bounds['log_aniso'] = [0., 3.]
        self.pars_bounds['rel_vol_frac_calcite'] = [0., 1.]

        # Define dictionary of temporal data limits
        self.temp_data_bounds = dict()
        self.temp_data_bounds['co2_rate'] = ['CO2 rate', 0., 10**1.5]
        self.temp_data_bounds['brine_rate'] = ['Brine rate', 0., 10**1.5]
        self.temp_data_bounds['co2_mass'] = ['CO2 mass', 0., (10**1.5)*70*365.25*86400]
        self.temp_data_bounds['brine_mass'] = ['Brine mass', 0., (10**1.5)*70*365.25*86400]

        # Initiate solution object
        self.sol = fgarom.Solution()

        debug_msg = 'FutureGen2Aquifer component created with name {}'.format(name)
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
                        'Parameter {} of FutureGen2Aquifer component {} ',
                        'with value {} is out of boundaries.']).format(
                            key, self.name, val)
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of FutureGen2Aquifer component {}.']).format(key, self.name)
                logging.warning(warn_msg)

    def check_temporal_inputs(self, time, temp_inputs):
        """
        Check whether temporal data fall within specified boundaries.

        :param temp_inputs: temporal input data of component model
        :type temp_inputs: dict
        """
        for key, val in temp_inputs.items():
            if ((val < self.temp_data_bounds[key][1]) or
                    (val > self.temp_data_bounds[key][2])):
                warn_msg = ''.join([
                    'Input data {} for component {} (value {}) is outside the model range ',
                    '{} to {} at time t = {}.']).format(
                        self.temp_data_bounds[key][0], self.name, val,
                        self.temp_data_bounds[key][1],
                        self.temp_data_bounds[key][2], time)
                logging.warning(warn_msg)

    def simulation_model(self, p, co2_rate=None, brine_rate=None,
                         co2_mass=None, brine_mass=None, time_point=365.25):
        """
        Return volume of impacted aquifer based on several metrics.

        :param p: input parameters of FutureGen2 aquifer model
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

        :returns: out - dictionary of observations of aquifer impacted volumes
            model;
            keys: ['TDS_volume', 'TDS_dx', 'TDS_dy', 'TDS_dz',
                   'pH_volume', 'pH_dx', 'pH_dy', 'pH_dz',
                   'Pressure_volume', 'Pressure_dx', 'Pressure_dy', 'Pressure_dz',
                   'Dissolved_CO2_volume', 'Dissolved_CO2_dx',
                   'Dissolved_CO2_dy', 'Dissolved_CO2_dz']
        """
        # Check assign default values
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        actual_p.update(p)

        # Check whether CO2 and brine rate inputs satisfy the model requirements
        self.check_temporal_inputs(time_point, dict(list(zip(
            ['co2_rate', 'brine_rate', 'co2_mass', 'brine_mass'],
            [co2_rate, brine_rate, co2_mass, brine_mass]))))

        if (co2_rate is not None) and (brine_rate is not None) and (
                co2_mass is not None) and (brine_mass is not None):
            co2_rate = max(0, co2_rate)
            co2_mass = max(0, co2_mass)
            brine_rate = max(0, brine_rate)
            brine_mass = max(0, brine_mass)
        else:
            error_msg = ''.join(['Dynamic input data is not provided for ',
                                 'FutureGen2Aquifer component {}']).format(self.name)
            logging.error(error_msg)
            raise ValueError(error_msg)

        log_co2_rate = np.log1p(co2_rate)
        log_brine_rate = np.log1p(brine_rate)
        log_co2_mass = np.log1p(co2_mass)
        log_brine_mass = np.log1p(brine_mass)

        inputArray = np.array([actual_p['aqu_thick'], actual_p['depth'], actual_p['por'],
                               actual_p['log_permh'], actual_p['log_aniso'],
                               actual_p['rel_vol_frac_calcite'],
                               log_co2_rate, log_brine_rate, time_point/365.25,
                               log_co2_mass, log_brine_mass])

        labels = ['TDS_volume', 'TDS_dx', 'TDS_dy', 'TDS_dz',
                  'pH_volume', 'pH_dx', 'pH_dy', 'pH_dz',
                  'Pressure_volume', 'Pressure_dx', 'Pressure_dy', 'Pressure_dz',
                  'Dissolved_CO2_volume', 'Dissolved_CO2_dx',
                  'Dissolved_CO2_dy', 'Dissolved_CO2_dz']

        if time_point/365.25 > 0:
            self.sol.find(inputArray)
            out = dict(zip(labels, self.sol.Outputs))
        else:
            out = dict(zip(labels, np.zeros(len(labels))))
        return out

    # Attributes for system connections
    system_params = ['{aquifer_name}Thickness']
    system_collected_inputs = {'co2_rate': 'CO2_{aquifer_name}',
                               'brine_rate': 'brine_{aquifer_name}',
                               'co2_mass': 'mass_CO2_{aquifer_name}',
                               'brine_mass': 'mass_brine_{aquifer_name}'}
    adapters = ['RateToMass']
    needsXY = False

    def connect_with_system(self, component_data, name2obj_dict,
                            locations, system_adapters, **kwargs):
        """
        Code to add FutureGen2 aquifer to system model for control file interface.

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
            aquifer to which the leakage rates of brine and CO2 from
            the connected wellbore component are transformed into the masses)
        :type system_adapters: dict

        :returns: None
        """
        # Process parameters data if any
        process_parameters(self, component_data, name2obj_dict)

        # Process dynamic inputs if any
        process_dynamic_inputs(self, component_data)

        if 'AquiferName' not in component_data:
            err_msg = ''.join(["Required argument 'AquiferName' is missing ",
                               "in the setup of FutureGen 2 Aquifer component. ",
                               "It should be setup to 'aquifer#' where # ",
                               "is an index of aquifer of interest."])
            logging.error(err_msg)
            raise KeyError(err_msg)
        else:
            aq_name = component_data['AquiferName']


        # Make model connections
        if 'Connection' in component_data:
            connection = None
            try:
                connection = name2obj_dict[component_data['Connection']]
            except KeyError:
                pass

            for ad_nm in system_adapters:
                if (system_adapters[ad_nm]['Connection'] == component_data['Connection']) and (
                        system_adapters[ad_nm]['AquiferName'] == aq_name):
                    aname = ad_nm

            adapter = name2obj_dict[aname]

            connection.add_obs_to_be_linked('CO2_{aquifer_name}'.format(
                aquifer_name=aq_name))
            adapter.add_obs_to_be_linked('mass_CO2_{aquifer_name}'.format(
                aquifer_name=aq_name))
            connection.add_obs_to_be_linked('brine_{aquifer_name}'.format(
                aquifer_name=aq_name))
            adapter.add_obs_to_be_linked('mass_brine_{aquifer_name}'.format(
                aquifer_name=aq_name))
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
        # End Connection if statement

        # Take care of parameter depth (to the bottom of the aquifer)
        strata = name2obj_dict['strata']
        if 'depth' not in self.pars and \
                'depth' not in self.deterministic_pars:
            self.add_par_linked_to_par(
                'depth', strata.composite_pars['{}Depth'.format(aq_name)])

        # Take care of parameter thick (possibly) defined by stratigraphy component
        if 'aqu_thick' not in self.pars and \
                'aqu_thick' not in self.deterministic_pars:
            sparam = '{aq}Thickness'.format(aq=aq_name)
            connect = None
            if sparam in strata.pars:
                connect = strata.pars
            elif sparam in strata.deterministic_pars:
                connect = strata.deterministic_pars
            elif sparam in strata.default_pars:
                connect = strata.default_pars
            if connect is None:
                sparam = 'aquiferThickness'
                if sparam in strata.pars:
                    connect = strata.pars
                elif sparam in strata.deterministic_pars:
                    connect = strata.deterministic_pars
                elif sparam in strata.default_pars:
                    connect = strata.default_pars
                else:
                    err_msg = ''.join([
                        'Unable to find "{}" or "{}" parameters. Please ',
                        'check setup of the stratigraphy.']).format(
                            sparam, '{aq}Thickness'.format(aq=aq_name))
                    logging.error(err_msg)
                    raise KeyError(err_msg)

            self.add_par_linked_to_par('aqu_thick', connect[sparam])


if __name__ == "__main__":
    # Create system model
    time_array = 365.25*np.arange(0, 2)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add FutureGen aquifer model object and define parameters
    fga = sm.add_component_model_object(FutureGen2Aquifer(name='fga', parent=sm))

    # St Peter Sandstone
    fga.add_par('aqu_thick', value=61.6)
    fga.add_par('depth', value=591.9)
    fga.add_par('por', value=0.18)
    fga.add_par('log_permh', value=-11.92)
    fga.add_par('log_aniso', value=0.30)
    fga.add_par('rel_vol_frac_calcite', value=0.01)

    fga.model_kwargs['brine_rate'] = 1.0e-3                    # kg/s
    fga.model_kwargs['co2_rate'] = 1.0e-2                      # kg/s
    fga.model_kwargs['brine_mass'] = 1.0e-3 * 1 *86400*365.25  # kg
    fga.model_kwargs['co2_mass'] = 1.0e-2 * 1 *86400*365.25    # kg

    # Add observations (output) from the aquifer model
    fga.add_obs('TDS_volume')
    fga.add_obs('TDS_dx')
    fga.add_obs('TDS_dy')
    fga.add_obs('TDS_dz')

    fga.add_obs('Pressure_volume')
    fga.add_obs('Pressure_dx')
    fga.add_obs('Pressure_dy')
    fga.add_obs('Pressure_dz')

    fga.add_obs('pH_volume')
    fga.add_obs('pH_dx')
    fga.add_obs('pH_dy')
    fga.add_obs('pH_dz')

    fga.add_obs('Dissolved_CO2_volume')
    fga.add_obs('Dissolved_CO2_dx')
    fga.add_obs('Dissolved_CO2_dy')
    fga.add_obs('Dissolved_CO2_dz')

    # Run the system model
    sm.forward()

    # Print the observations
    np.set_printoptions(precision=2)
    print('TDS_volume', sm.collect_observations_as_time_series(fga, 'TDS_volume'))
    print('TDS_dx', sm.collect_observations_as_time_series(fga, 'TDS_dx'))
    print('TDS_dy', sm.collect_observations_as_time_series(fga, 'TDS_dy'))
    print('TDS_dz', sm.collect_observations_as_time_series(fga, 'TDS_dz'))

    print('Pressure_volume',
          sm.collect_observations_as_time_series(fga, 'Pressure_volume'))
    print('Pressure_dx',
          sm.collect_observations_as_time_series(fga, 'Pressure_dx'))
    print('Pressure_dy',
          sm.collect_observations_as_time_series(fga, 'Pressure_dy'))
    print('Pressure_dz',
          sm.collect_observations_as_time_series(fga, 'Pressure_dz'))

    print('pH_volume', sm.collect_observations_as_time_series(fga, 'pH_volume'))
    print('pH_dx', sm.collect_observations_as_time_series(fga, 'pH_dx'))
    print('pH_dy', sm.collect_observations_as_time_series(fga, 'pH_dy'))
    print('pH_dz', sm.collect_observations_as_time_series(fga, 'pH_dz'))

    print('Dissolved_CO2_volume',
          sm.collect_observations_as_time_series(fga, 'Dissolved_CO2_volume'))
    print('Dissolved_CO2_dx',
          sm.collect_observations_as_time_series(fga, 'Dissolved_CO2_dx'))
    print('Dissolved_CO2_dy',
          sm.collect_observations_as_time_series(fga, 'Dissolved_CO2_dy'))
    print('Dissolved_CO2_dz',
          sm.collect_observations_as_time_series(fga, 'Dissolved_CO2_dz'))

    # Expected output
    # TDS_volume [  0.         221.65165841]
    # TDS_dx [0.         7.01228794]
    # TDS_dy [0.         7.01228794]
    # TDS_dz [ 0.         10.70068883]
    # Pressure_volume [ 0.         24.96894867]
    # Pressure_dx [0.         2.36196175]
    # Pressure_dy [0.         2.36196175]
    # Pressure_dz [0.         1.76299383]
    # pH_volume [     0.         283401.12242706]
    # pH_dx [  0.         100.29078768]
    # pH_dy [  0.         100.29078768]
    # pH_dz [ 0.        23.9857807]
    # Dissolved_CO2_volume [     0.         118108.91206191]
    # Dissolved_CO2_dx [ 0.         59.56185527]
    # Dissolved_CO2_dy [ 0.         59.56185527]
    # Dissolved_CO2_dz [ 0.         34.14447187]
