# -*- coding: utf-8 -*-
# Created on July 2, 2019
# @author: Kayyum Mansoor
# mansoor1@llnl.gov
import logging
import os
import sys
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import SystemModel, ComponentModel
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

from openiam.cfi.commons import process_parameters, process_dynamic_inputs



class AlluviumAquiferLF(ComponentModel):
    """
    The Alluvium Aquifer LF (low flux) component model is a reduced order model
    which can be used to predict the changes in diluted groundwater chemistry if |CO2|
    and brine were to leak into an overlying alluvium aquifer similar to the High
    Plains aquifer, Haskel County, Kansas, USA. The protocol allows uncertainty and
    variability in aquifer heterogeneity, fluid transport and geochemical reactions
    to be collectively evaluated to assess potential changes in groundwater pH,
    total dissolved solids (TDS), As, Ba, Cd, Pb, benzene, naphthalene, and phenol
    concentrations

    Although the alluvium aquifer ROM was developed using site-specific data
    from the High Plains aquifer, the model accepts aquifer characteristics as
    variable inputs and, therefore, may have more broad applicability. Careful
    consideration should be given to the hydrogeochemical character of the
    aquifer before using this model at a new site.

    Component model input parameters:

    * **sandFraction** [-] (0.6 to 0.9) - sand fraction (default: 0.7087)

    * **correlationLengthX** [|m|] (200 to 2500) - correlation length in x-direction
      (default: 473.9740)

    * **correlationLengthZ** [|m|] (0.5 to 25) - correlation length in  z-direction
      (default: 22.5966)

    * **permeabilitySand** [|log10| |m^2|] (-14 to -10) - permeability of sand
      (default: -13.3634)

    * **permeabilityClay** [|log10| |m^2|] (-18 to -15) - permeability of clay
      (default: -15.5075)

    * **NaMolality** [|log10| molality] (-3 to 1) - sodium molality
      (default: 0.9279)

    * **PbMolality** [|log10| molality] (-9.1836 to -5.6836) - lead molality
      (default: -6.2021)

    * **benzeneMolality** [|log10| molality] (-10 to -6) - benzene molality
      (default:  -7.7138)

    * **tMitigation** [years] (1 to 100) - mitigation time (default: 7.7387)

    Component model dynamic inputs:

    * **brine_rate** [|kg/s|] (0 to 0.00075) - brine rate (default: 0.000145878)

    * **brine_mass** [|kg|] (1.584e+3 to 1.698e+6) - cumulative brine mass
      (default: 2.301e+4)

    * **co2_rate** [|kg/s|] (0 to 0.005001) - |CO2| rate (default: 0.6883e-5)

    * **co2_mass** [|kg|] (3.388 to 5.623e+6) - cumulative |CO2| mass
      (default: 543.01)

    * **simtime** [years] (1 to 200) - simulation time (default: 5)

    Observations of the Alluvium Aquifer (low flux) component are:

    * **TDS_volume** [|m^3|] - volume of aquifer above TDS threshold in |millig/L|
      (TDS > 1300 |millig/L|)

    * **pH_volume** [|m^3|] - volume of aquifer below pH threshold (pH < 7.0)

    * **As_volume** [|m^3|] - volume of aquifer above arsenic threshold in |microg/L|
      (arsenic > 9.3 |microg/L|)

    * **Ba_volume** [|m^3|] - volume of aquifer above barium threshold in |microg/L|
      (barium > 140 |microg/L|)

    * **Cd_volume** [|m^3|] - volume of aquifer above cadmium threshold in |microg/L|
      (cadmium > 0.25 |microg/L|)

    * **Pb_volume** [|m^3|] - volume of aquifer above lead threshold in |microg/L|
      (lead > 0.63 |microg/L|)

    * **Benzene_volume** [|m^3|] - volume of aquifer above benzene threshold
      (benzene > 0.03 |microg/L|)

    * **Naphthalene_volume** [|m^3|] - volume of aquifer above naphthalene threshold
      (napthalene > 0.2 |microg/L|)

    * **Phenol_volume** [|m^3|] - volume of aquifer above phenol threshold
      (phenol > 0.003 |microg/L|).
    """
    def __init__(self, name, parent):
        """
        Constructor method of AlluviumAquiferLF class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :returns: AlluviumAquiferLF class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25}  # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'AlluviumAquiferLF'

        # Check if another alluvium aquifer ml component is already part of the system model
        orig_aac_cmpnt = self._parent.ml_models_components.get(
            'AlluviumAquiferLF', None)

        try:
            import components.aquifer.alluvium_lf.alluvium_aquifer_lf_rom as aalfrom
        except ImportError:
            print('\nERROR: Unable to load ROM for Alluvium Aquifer (low flux) component\n')
            sys.exit()

        if orig_aac_cmpnt is None:
            # Register the component as the one using ml models
            # The component will be registered only if no other alluvium aquifer lf
            # components were added to the same system model
            self._parent.register_ml_model_component(self, 'AlluviumAquiferLF')

            # Initiate solution object
            self.sol = aalfrom.Solution()
        else:
            # Get solution class object from existing component of the same type
            self.sol = orig_aac_cmpnt.sol

        self.add_default_par('sandFraction', value=0.7087)
        self.add_default_par('correlationLengthX', value=473.9740)
        self.add_default_par('correlationLengthZ', value=22.5966)
        self.add_default_par('logK_sand', value=-13.3634)
        self.add_default_par('logK_clay', value=-15.5075)
        self.add_default_par('NaMolality', value=0.9279)
        self.add_default_par('PbMolality', value=-6.2021)
        self.add_default_par('benzeneMolality', value=-7.7138)
        self.add_default_par('tMitigation', value=7.7387)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['sandFraction'] = [0.60, 0.90]
        self.pars_bounds['correlationLengthX'] = [200.0, 2500.0]
        self.pars_bounds['correlationLengthZ'] = [0.5, 25.0]
        self.pars_bounds['logK_sand'] = [-14.0, -10.0]
        self.pars_bounds['logK_clay'] = [-18.0, -15.0]
        self.pars_bounds['NaMolality'] = [-3.0, 1.0]
        self.pars_bounds['PbMolality'] = [-9.1836, -5.6836]
        self.pars_bounds['benzeneMolality'] = [-10, -6]
        self.pars_bounds['tMitigation'] = [1.0, 100.0]

        # Define dictionary of temporal data limits
        self.temp_data_bounds = dict()
        self.temp_data_bounds['co2_rate'] = ['CO2 rate', 0.0, 0.005001]
        self.temp_data_bounds['brine_rate'] = ['brine rate', 0.0, 0.00075]
        self.temp_data_bounds['co2_mass'] = ['CO2 mass', 10**0.53, 10**6.75]
        self.temp_data_bounds['brine_mass'] = ['brine mass', 10**3.20, 10**6.23]
        self.temp_data_bounds['simtime'] = ['simulation time', 1.0, 200.0]

        # Define observations names
        self.out_labels = ['TDS', 'pH', 'As', 'Pb', 'Cd', 'Ba',
                           'Benzene', 'Naphthalene', 'Phenol']
        self.out_labels = [val+'_volume' for val in self.out_labels]

        debug_msg = 'AlluviumAquiferLF component created with name {name}'.format(name=name)
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
                    warn_msg = ''.join([
                        'Parameter {} of AlluviumAquiferLF component {} ',
                        'is out of bounds.']).format(key, self.name)
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of a AlluviumAquiferLF component {}.']).format(key, self.name)
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
                    'Temporal input {} of AlluviumAquiferLF component {} ',
                    'is outside the model range [{}, {}] at time t = {}']).format(
                        self.temp_data_bounds[key][0].lower(), self.name,
                        self.temp_data_bounds[key][1],
                        self.temp_data_bounds[key][2], time)
                logging.warning(warn_msg)

    def simulation_model(self, p, co2_rate=0.6883e-5, brine_rate=0.000145878,
                         co2_mass=10**2.73481, brine_mass=10**4.36206,
                         time_point=365.25):
        """
        Return of impacted aquifer based on several metrics.

        :param p: input parameters of alluvium aquifer model
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

        :returns: vol - dictionary of observations of alluvium aquifer impacted
            volumes; keys:
            ['TDS_volume', 'pH_volume', 'As_volume', 'Pb_volume', 'Cd_volume',
            'Ba_volume', 'Benzene_volume', 'Naphthalene_volume', 'Phenol_volume']
        """
        # Update assigned default values
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        actual_p.update(p)

        # Check whether CO2, brine, rates and mass inputs satisfy the model requirements
        # Conditional statement used to supressess excessive error messsages
        if (co2_mass is not None) and (brine_mass is not None):
            if co2_mass <= 0:
                warn_msg = ''.join([
                    'Temporal input CO2 mass of AlluviumAquiferLF component {} ',
                    'is below the model range [{}, {}] at time t = {}.',
                    'The current value will be replaced by ',
                    'the lower boundary value.']).format(
                        self.name, self.temp_data_bounds['co2_mass'][1],
                        self.temp_data_bounds['co2_mass'][2], time_point)
                logging.warning(warn_msg)
                co2_mass = self.temp_data_bounds['co2_mass'][1]
            if brine_mass <= 0:
                warn_msg = ''.join([
                    'Temporal input brine mass of AlluviumAquiferLF component {} ',
                    'is below the model range [{}, {}] at time t = {}.',
                    'The current value will be replaced by ',
                    'the lower boundary value.']).format(
                        self.name, self.temp_data_bounds['brine_mass'][1],
                        self.temp_data_bounds['brine_mass'][2], time_point)
                logging.warning(warn_msg)
                brine_mass = self.temp_data_bounds['brine_mass'][1]
        else:
            error_msg = 'CO2 and/or brine mass arguments of the model method are not provided.'
            logging.error(error_msg)

        if (time_point >= 365.25) and ((
                np.log10(co2_mass) >= self.temp_data_bounds['co2_mass'][1]) or (
                    np.log10(brine_mass) >= self.temp_data_bounds['brine_mass'][1])):
            self.check_temporal_inputs(
                time_point, dict(list(zip(
                    ['co2_rate', 'brine_rate', 'co2_mass', 'brine_mass', 'simtime'],
                    [co2_rate, brine_rate, np.log10(co2_mass),
                     np.log10(brine_mass), time_point/365.25]))))

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
            actual_p['co2_rate'], actual_p['co2_mass'], actual_p['sandFraction'],
            actual_p['correlationLengthX'], actual_p['correlationLengthZ'],
            actual_p['logK_sand'], actual_p['logK_clay'],
            actual_p['NaMolality'], actual_p['PbMolality'],
            actual_p['benzeneMolality'], actual_p['tMitigation']])

        if time_point/365.25 > 0:
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
        Code to add alluvium aquifer to system model for control file interface.

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
    time_array = 365.25*np.arange(0.0, 6.0)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add alluvium aquifer model object and define parameters
    aalf = sm.add_component_model_object(AlluviumAquiferLF(name='aalf', parent=sm))
    aalf.add_par('sandFraction', value=0.7087)
    aalf.add_par('correlationLengthX', value=473.9740)
    aalf.add_par('correlationLengthZ', value=22.5966)
    aalf.add_par('logK_sand', value=-13.3634)
    aalf.add_par('logK_clay', value=-15.5075)
    aalf.add_par('NaMolality', value=0.9279)
    aalf.add_par('PbMolality', value=-6.2021)
    aalf.add_par('benzeneMolality', value=-7.7138)
    aalf.add_par('tMitigation', value=7.7387)

    aalf.model_kwargs['co2_rate'] = 0.000006883      # kg/s
    aalf.model_kwargs['co2_mass'] = 10**2.73481      # kg
    aalf.model_kwargs['brine_rate'] = 0.000145878    # kg/s
    aalf.model_kwargs['brine_mass'] = 10**4.36206    # kg

    # Add observations (output) of the alluvium aquifer model
    aalf.add_obs('TDS_volume')
    aalf.add_obs('pH_volume')
    aalf.add_obs('As_volume')
    aalf.add_obs('Pb_volume')
    aalf.add_obs('Cd_volume')
    aalf.add_obs('Ba_volume')
    aalf.add_obs('Benzene_volume')
    aalf.add_obs('Naphthalene_volume')
    aalf.add_obs('Phenol_volume')

    # Run the system model
    sm.forward()

    # Print the observations
    print('{:20}'.format('TDS volume'),
          sm.collect_observations_as_time_series(aalf, 'TDS_volume'))
    print('{:20}'.format('pH volume'),
          sm.collect_observations_as_time_series(aalf, 'pH_volume'))
    print('{:20}'.format('As volume'),
          sm.collect_observations_as_time_series(aalf, 'As_volume'))
    print('{:20}'.format('Pb volume'),
          sm.collect_observations_as_time_series(aalf, 'Pb_volume'))
    print('{:20}'.format('Cd volume'),
          sm.collect_observations_as_time_series(aalf, 'Cd_volume'))
    print('{:20}'.format('Ba volume'),
          sm.collect_observations_as_time_series(aalf, 'Ba_volume'))
    print('{:20}'.format('Benzene volume'),
          sm.collect_observations_as_time_series(aalf, 'Benzene_volume'))
    print('{:20}'.format('Naphthalene volume'),
          sm.collect_observations_as_time_series(aalf, 'Naphthalene_volume'))
    print('{:20}'.format('Phenol volume'),
          sm.collect_observations_as_time_series(aalf, 'Phenol_volume'))

    # Expected output
    # TDS volume [0., 243.60472506, 315.42085149, 681.51000667, 1464.27143244, 2817.58361957]
    # pH volume [0., 135.55234468, 547.96051906, 605.7159738, 606.48105472, 678.17416717]
    # As volume [0., 10., 17.32817682, 25.83795442, 53.71518848, 97.3487519]
    # Pb volume [0., 179.47389199, 1554.30216956, 5972.94243509, 13346.66092414, 24999.17432054]
    # Cd volume [0., 105.02705854, 2362.73319012, 9383.10921071, 11061.49479283, 8258.10081888]
    # Ba volume [0., 3973.21653329, 8082.43147709, 10288.42580438, 12179.40618369, 13643.12399265]
    # Benzene volume [0., 62.15296578, 353.59654512, 714.91883593, 1223.97332442, 2131.52468324]
    # Naphthalene volume [0., 38.76136875, 133.63125123, 538.19287788, 1639.28247224, 3103.13930929]
    # Phenol volume [0., 4436.94112399, 4309.9468996, 4803.11442366, 7642.13761712, 12014.17764713]
