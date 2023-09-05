# -*- coding: utf-8 -*-
# Created on Sep 21, 2021
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

import components.aquifer.generic.generic_aquifer_ROM as genrom
from openiam.cfi.commons import process_parameters, process_dynamic_inputs

GA_SCALAR_OBSERVATIONS = [
    'Dissolved_salt_volume', 'Dissolved_salt_dr', 'Dissolved_salt_dz',
    'Dissolved_CO2_volume', 'Dissolved_CO2_dr', 'Dissolved_CO2_dz']

GA_GRID_OBSERVATIONS = ['Dissolved_CO2_mass_fraction', 'Dissolved_salt_mass_fraction']
GA_GRID_COORDINATES = ['r_coordinate', 'z_coordinate']


class GenericAquifer(ComponentModel):
    """
    The Generic Aquifer component model is a surrogate model
    that can be used to predict the leakage of carbon dioxide (|CO2|) and
    brine from a |CO2| storage reservoir. The model predicts
    the mass fraction of |CO2| and salt on a 100x10 radial grid surrounding the
    leaky well and outputs these as gridded observations.  The model also predicts
    the volume and dimensions of aquifer where pore water concentrations exceed
    specified threshold values of dissolved |CO2| and salt.

    The Generic Aquifer model is a machine learning regression model fitted to the
    results of STOMP-CO2E-R multiphase flow and reactive transport simulations of
    |CO2| and brine leakage using Tensorflow 2.4.  50,000 nonisothermal multiphase flow
    simulations were used to train the Generic Aquifer component model.  Input
    parameters were varied using Latin Hypercube Sampling across wide ranges.

    In the NRAP-Open-IAM control file, the type name for the Generic Aquifer component
    is ``GenericAquifer``. The description of the component's parameters are
    provided below:

    * **aqu_thick** [|m|] (25 to 250) - thickness of unit (default: 33.2);
      *linked to Stratigraphy*

    * **top_depth** [|m|] (100 to 4100) - depth to the top of the aquifer (default: 590.1);
      *linked to Stratigraphy*

    * **por** [-] (0.02 to 0.25) - porosity of unit (default: 0.118)

    * **log_permh** [|log10| |m^2|] (-14 to -10) - horizontal permeability
      (default: -13.39)

    * **log_aniso** [-] (0 to 3) - anisotropy ratio Kh/Kv (default: 0.3)

    * **aquifer_salinity** [-] (0.0 to 0.015) - salt mass fraction
      in aquifer water (default: 0.005)

    * **reservoir_salinity** [-] (0.015 to 0.05) - salt mass fraction
      in leak water (default: 0.03)

    * **dissolved_salt_threshold** [-] (0.0 to 1.0) - threshold for salt mass
      fraction (default: 0.02)

    * **dissolved_co2_threshold** [-] (0.0 to 1.0) - threshold for CO2 mass
      fraction (default: 0.01)

    Component model dynamic inputs:

    * **brine_mass** [|kg|] (0 to 6.985e+10) - cumulative brine mass

    * **co2_mass** [|kg|] (0 to 6.985e+10) - cumulative |CO2| mass.

    Observations from the Generic Aquifer component are:

    * **Dissolved_salt_volume** [|m^3|] - volume of plume where relative change
      in salt mass fraction > dissolved_salt_threshold

    * **Dissolved_salt_dr** [|m|] - radius of plume where relative change
      in salt mass fraction > dissolved_salt_threshold

    * **Dissolved_salt_dz** [|m|] - height of plume where relative change
      in salt mass fraction > dissolved_salt_threshold

    * **Dissolved_CO2_volume** [|m^3|] - volume of plume where dissolved
      |CO2| mass fraction > dissolved_co2_threshold

    * **Dissolved_CO2_dr** [|m|] - radius of plume where dissolved
      |CO2| mass fraction > dissolved_co2_threshold

    * **Dissolved_CO2_dz** [|m|] - height of plume where dissolved
      |CO2| mass fraction > dissolved_co2_threshold

    Gridded observations from the Generic Aquifer component are:

    * **Dissolved_CO2_mass_fraction** [-] - mass fraction of |CO2| in aquifer
      pore water on a 100x10 radial grid surrounding the leaky well

    * **Dissolved_salt_mass_fraction** [-] - mass fraction of salt in aquifer
      pore water on a 100x10 radial grid surrounding the leaky well

    * **r_coordinate** [m] - radial coordinates of the points in the
      100x10 radial grid surrounding the leaky well. The 100 radii are within
      range from 1.62 m to about 77.5 km.

    * **z_coordinate** [m] - depth coordinates of the points in the 100x10
      radial grid surrounding the leaky well. The 10 depths used are within
      the aquifer modeled by the GenericAquifer. The minimum depth is 5%
      of the aquifer's thickness above the base of the aquifer, while
      the maximum depth is 95% of the aquifer's thickness above the base
      of the aquifer. The increment used between depth values is 10%
      of the aquifer's thickness.

    """
    def __init__(self, name, parent):
        """
        Constructor method of GenericAquifer class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :returns: GenericAquifer class object
        """
        if not self.model_data_check():
            model_files_folder = os.sep.join([IAM_DIR, 'source', 'components',
                                              'aquifer', 'generic'])
            error_msg = ''.join(['Generic Aquifer component {} cannot be created ',
                                 'as required model files are missing ',
                                 'in the folder\n {}.']).format(
                                     name, model_files_folder)
            # TODO Add instructions on where to download required model files once
            # we put them on EDX
            logging.error(error_msg)
            raise FileNotFoundError(error_msg)

        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25}

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'GenericAquifer'

        # Define default parameters
        self.add_default_par('aqu_thick', value=33.2)
        self.add_default_par('top_depth', value=590.1)
        self.add_default_par('por', value=0.118)
        self.add_default_par('log_permh', value=-13.39)
        self.add_default_par('log_aniso', value=0.3)
        self.add_default_par('aquifer_salinity', value=0.005)
        self.add_default_par('reservoir_salinity', value=0.03)
        self.add_default_par('dissolved_salt_threshold', value=0.02)
        self.add_default_par('dissolved_co2_threshold', value=0.01)

        # Define dictionary of parameter boundaries
        self.pars_bounds = dict()
        self.pars_bounds['aqu_thick'] = [25., 250.]
        self.pars_bounds['top_depth'] = [100., 4100.]
        self.pars_bounds['por'] = [0.02, 0.25]
        self.pars_bounds['log_permh'] = [-14., -10.]
        self.pars_bounds['log_aniso'] = [0., 3.]
        self.pars_bounds['aquifer_salinity'] = [0., 0.015]
        self.pars_bounds['reservoir_salinity'] = [0.015, 0.05]
        self.pars_bounds['dissolved_salt_threshold'] = [0.0, 1.0]
        self.pars_bounds['dissolved_co2_threshold'] = [0.0, 1.0]

        # Define dictionary of temporal data limits
        self.temp_data_bounds = dict()
        self.temp_data_bounds['co2_mass'] = ['CO2 mass', 0., (10**1.5)*70*365.25*86400]
        self.temp_data_bounds['brine_mass'] = ['Brine mass', 0., (10**1.5)*70*365.25*86400]

        # Define gridded observations names
        self.grid_obs_keys = GA_GRID_OBSERVATIONS

        # Initiate solution object
        self.sol = genrom.Solution()

        debug_msg = 'GenericAquifer component created with name {}'.format(name)
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
                        'Parameter {} of GenericAquifer component {} with value {}',
                        'is out of boundaries.']).format(key, self.name, val)
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of GenericAquifer component {}.']).format(key, self.name)
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
                    'Input data {} (value {}) is outside the model range ',
                    '{} to {} at time t = {}.']).format(
                        self.temp_data_bounds[key][0], val,
                        self.temp_data_bounds[key][1],
                        self.temp_data_bounds[key][2], time)
                logging.warning(warn_msg)

    # Attributes for system connections
    system_params = ['{aquifer_name}Thickness']
    system_collected_inputs = {'co2_mass': 'mass_CO2_{aquifer_name}',
                               'brine_mass': 'mass_brine_{aquifer_name}'}
    adapters = ['RateToMass']
    needsXY = False

    def connect_with_system(self, component_data, name2obj_dict,
                            locations, system_adapters, **kwargs):
        """
        Code to add Generic Aquifer aquifer to system model for control file interface.

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
                               "in the setup of Generic Aquifer component. ",
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

            adapter.add_obs_to_be_linked('mass_CO2_{aquifer_name}'.format(
                aquifer_name=aq_name))
            adapter.add_obs_to_be_linked('mass_brine_{aquifer_name}'.format(
                aquifer_name=aq_name))
            self.add_kwarg_linked_to_obs(
                'co2_mass', adapter.linkobs['mass_CO2_{aquifer_name}'.format(
                    aquifer_name=aq_name)])
            self.add_kwarg_linked_to_obs(
                'brine_mass', adapter.linkobs['mass_brine_{aquifer_name}'.format(
                    aquifer_name=aq_name)])
        # End Connection if statement

        # Take care of parameter top_depth (to the top of the aquifer)
        strata = name2obj_dict['strata']
        if 'top_depth' not in self.pars and \
                'top_depth' not in self.deterministic_pars:
            aq_index = int(aq_name[7:])  # index of aquifer of interest
            # depth to bottom of shale layer above the aquifer of interest
            above_shale_name = 'shale{}Depth'.format(aq_index+1)
            self.add_par_linked_to_par(
                'top_depth', strata.composite_pars[above_shale_name])

        # Take care of parameter aqu_thick (possibly) defined by stratigraphy component
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
            if not connect:
                sparam = 'aquiferThickness'
                if sparam in strata.pars:
                    connect = strata.pars
                elif sparam in strata.deterministic_pars:
                    connect = strata.deterministic_pars
                elif sparam in strata.default_pars:
                    connect = strata.default_pars
                else:
                    print('Unable to find parameter ' + sparam)

            self.add_par_linked_to_par('aqu_thick', connect[sparam])

        if 'Outputs' in component_data:
            comp_outputs = component_data['Outputs']
            new_comp_outputs = []
            for obs_nm in comp_outputs:
                if obs_nm in GA_SCALAR_OBSERVATIONS:
                    new_comp_outputs.append(obs_nm)
                    self.add_obs(obs_nm)
                elif obs_nm in GA_GRID_OBSERVATIONS:
                    for ind1 in range(100):
                        for ind2 in range(10):
                            augm_obs_nm = obs_nm+'_coord_{}_{}'.format(ind1+1, ind2+1)
                            self.add_local_obs(augm_obs_nm, obs_nm, 'matrix', (ind1, ind2))
                            new_comp_outputs.append(augm_obs_nm)
                    # Add gridded observations
                    self.add_grid_obs(
                        obs_nm, constr_type='matrix', output_dir=kwargs['output_dir'])
                elif obs_nm in GA_GRID_COORDINATES:
                    # Add the radii and depths for gridded observations
                    self.add_grid_obs(
                        obs_nm, constr_type='matrix', output_dir=kwargs['output_dir'])
                else:
                    warn_msg = ''.join([
                        '{} is not recognised as observation name ',
                        'of Generic Aquifer component {}.']).format(obs_nm, self.name)
                    logging.warning(warn_msg)

            component_data['Outputs'] = new_comp_outputs


    def simulation_model(self, p, co2_mass=None, brine_mass=None,
                         time_point=365.25, **kwargs):
        """
        Return volume of impacted aquifer based on several metrics.

        :param p: input parameters of Generic aquifer model
        :type p: dict

        :param co2_mass: Cumulative CO2 mass leaked, kg
        :type co2_mass: float

        :param brine_mass: Cumulative brine mass leaked, kg
        :type brine_mass: float

        :param time_point: time point in days at which the leakage rates are
            to be calculated; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :returns: out - dictionary of observations of aquifer impacted volumes
            and gridded observations of dissolved salt mass fraction and
            dissolved co2 mass fraction;
            keys: ['Dissolved_salt_volume', 'Dissolved_salt_dr',
                   'Dissolved_salt_dz', 'Dissolved_CO2_volume',
                   'Dissolved_CO2_dr', 'Dissolved_CO2_dz',
                   'Dissolved_salt_mass_fraction',
                   'Dissolved_CO2_mass_fraction']
        """
        # Assign default values
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        actual_p.update(p)

        # Check whether CO2 and brine rate inputs satisfy the model requirements
        self.check_temporal_inputs(time_point, dict(zip(
            ['co2_mass', 'brine_mass'], [co2_mass, brine_mass])))
        if (co2_mass is not None) and (brine_mass is not None):
            co2_mass = max(0, co2_mass)
            brine_mass = max(0, brine_mass)
        else:
            error_msg = ''.join(['Dynamic input data is not provided for ',
                                 'GenericAquifer component {}']).format(self.name)
            logging.error(error_msg)
            raise ValueError(error_msg)

        inputArray = np.array([actual_p['aqu_thick'], actual_p['top_depth'], actual_p['por'],
                               actual_p['log_permh'], actual_p['log_aniso'],
                               co2_mass, brine_mass,
                               actual_p['aquifer_salinity'], actual_p['reservoir_salinity'],
                               time_point/365.25, actual_p['dissolved_salt_threshold'],
                               actual_p['dissolved_co2_threshold']])

        self.sol.find(inputArray)

        scalar_outputs = GA_SCALAR_OBSERVATIONS

        gridded_outputs = GA_GRID_COORDINATES + GA_GRID_OBSERVATIONS

        out = {}

        for count, label in enumerate(scalar_outputs):
            out[label] = self.sol.Outputs[count]

        for count, label in enumerate(gridded_outputs):
            out[label] = self.sol.GriddedOutputs[count]

        return out


    @staticmethod
    def model_data_check():
        """
        Check whether model data files for Generic Aquiferr component were downloaded
        and placed to the right place.
        """
        model_data_dir = os.sep.join([
            IAM_DIR, 'source', 'components', 'aquifer', 'generic'])

        subfolders = ['{}-{}m'.format(100*ind, 100*(ind+4)-1) \
                      for ind in range(1, 38, 4)]

        data_files = ['co2_generator.h5', 'co2_input_transformer.pkl',
                      'co2_output_transformer.pkl', 'salt_generator.h5',
                      'salt_input_transformer.pkl', 'salt_output_transformer.pkl']

        check_flag = 0   # no issues: data is present

        for subfold in subfolders:
            for file in data_files:
                file_path = os.sep.join([model_data_dir, subfold, file])
                if not os.path.isfile(file_path):
                    check_flag = check_flag + 1

        return (check_flag == 0)


if __name__ == "__main__":
    # Create system model
    time = 1
    time_array = 365.25*np.arange(0, time+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add FutureGen aquifer model object and define parameters
    ga = sm.add_component_model_object(GenericAquifer(name='ga', parent=sm))

    # parameters from run 3263 of 6250 (2100-2499m)
    par = [1.0260E+02, -2.3066E+03, 7.9001E-02, -1.2434E+01, 2.2254E+00,
           5.7080E-03, 4.0938E-01, -5.0149E+00, 4.0617E-02]

    ga.add_par('aqu_thick', value=par[0])
    ga.add_par('top_depth', value=-par[1])
    ga.add_par('por', value=par[2])
    ga.add_par('log_permh', value=par[3])
    ga.add_par('log_aniso', value=par[4])
    ga.add_par('aquifer_salinity', value=par[5])
    ga.add_par('reservoir_salinity', value=par[8])
    ga.add_par('dissolved_salt_threshold', value=0.02)
    ga.add_par('dissolved_co2_threshold', value=0.01)

    log_co2_rate = par[6] # kg/s
    log_brine_rate = par[7] # kg/s
    ga.model_kwargs['co2_mass'] = 10**log_co2_rate * time * 365.25 * 86400   # kg
    ga.model_kwargs['brine_mass'] = 10**log_brine_rate * time * 365.25 * 86400   # kg

    # Add observations (output) from the aquifer model
    ga.add_obs('Dissolved_salt_volume')
    ga.add_obs('Dissolved_salt_dr')
    ga.add_obs('Dissolved_salt_dz')
    ga.add_obs('Dissolved_CO2_volume')
    ga.add_obs('Dissolved_CO2_dr')
    ga.add_obs('Dissolved_CO2_dz')

    # Run the system model
    sm.forward()

    # Print the observations
    np.set_printoptions(precision=2)
    print('Dissolved_salt_volume',
          sm.collect_observations_as_time_series(ga, 'Dissolved_salt_volume'))
    print('Dissolved_salt_dr',
          sm.collect_observations_as_time_series(ga, 'Dissolved_salt_dr'))
    print('Dissolved_salt_dz',
          sm.collect_observations_as_time_series(ga, 'Dissolved_salt_dz'))
    print('Dissolved_CO2_volume',
          sm.collect_observations_as_time_series(ga, 'Dissolved_CO2_volume'))
    print('Dissolved_CO2_dr',
          sm.collect_observations_as_time_series(ga, 'Dissolved_CO2_dr'))
    print('Dissolved_CO2_dz',
          sm.collect_observations_as_time_series(ga, 'Dissolved_CO2_dz'))

    # Expected output
    # Dissolved_salt_volume [   0.   2706.93]
    # Dissolved_salt_dr [0.   0.26]
    # Dissolved_salt_dz [0.   8.21]
    # Dissolved_CO2_volume [      0.   5373325.51]
    # Dissolved_CO2_dr [ 0.   12.31]
    # Dissolved_CO2_dz [  0.   180.58]
