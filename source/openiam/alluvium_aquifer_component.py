# -*- coding: utf-8 -*-
# Created on June 7, 2018
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

try:
    import components.aquifer.alluvium.alluvium_aquifer_ROM as aarom
except ImportError:
    print('\nERROR: Unable to load ROM for Alluvium Aquifer component\n')
    sys.exit()

class AlluviumAquifer(ComponentModel):
    """
    The Alluvium Aquifer component modelis a reduced order model which can be used
    to predict the changes in diluted groundwater chemistry if |CO2| and brine were to
    leak into an overlying alluvium aquifer similar to the High Plains aquifer,
    Haskel County, Kansas, USA. The protocol allows uncertainty and variability in
    aquifer heterogeneity, fluid transport and geochemical reactions to be
    collectively evaluated to assess potential changes in groundwater pH, total
    dissolved solids (TDS), As, Ba, Cd, Pb, benzene, naphthalene, and phenol
    concentrations by developing a scaling function that can be applied to correct
    the output from the hydrology ROM for geochemical reactions.

    Although the Alluvium Aquifer ROM was developed using site-specific data
    from the High Plains aquifer, the model accepts aquifer characteristics as
    variable inputs and, therefore, may have more broad applicability. Careful
    consideration should be given to the hydrogeochemical character of the
    aquifer before using this model at a new site.

    The Alluvium Aquifer will reject a number of simulation results, typically
    around 20 to 30 percent of results will be rejected.  When results are rejected
    the model will return NaN (Not a Number) values.  Rejected results need to
    be filtered out for analysis. For more information on the Alluvium Aquifer
    model development see Carroll et al., :cite:`Carroll2016`.

    Component model input definitions:

    * **sandFraction** [-] (0.35 to 0.65) - sand volume fraction (default: 0.422)

    * **correlationLengthX** [|m|] (200 to 2500) - correlation length in x-direction
      (default: 361.350)

    * **correlationLengthZ** [|m|] (0.5 to 25) - correlation length in z-direction
      (default: 17.382)

    * **permeabilityClay** [|log10| (|m^2|)] (-18 to -15) - permeability of clay
      (default: -16.340)

    * **NaMolality** [|log10| molality] (-3 to 1) - sodium molality (default: 0.121)

    * **PbMolality** [|log10| molality] (-8.5 to -5) - lead molality
      (default: -5.910)

    * **benzeneMolality** [|log10| molality] (0 to 1) - benzene molality
      (default: 0.109)

    * **tMitigation** [years] (1 to 200) - mitigation time (default: 87.914)

    * **CEC** [meq/100g] (0.1 to 40) - cation exchange capacity (default: 32.073)

    * **Asbrine** [|log10| molality] (-9 to -5) - arsenic concentration in the
      leaking brine (default: -5.397)

    * **Babrine** [|log10| molality] (-5.1 to -2.3) - barium concentration in
      the leaking brine (default: -3.397)

    * **Cdbrine** [|log10| molality] (-9 to -6) - cadmium concentration in the
      leaking brine (default: -8.574)

    * **Pbbrine** [|log10| molality] (-8.5 to -5) - lead concentration in the
      leaking brine (default: -7.719)

    * **Benzene_brine** [|log10| molality] (-8.8927 to -4.8927) - benzene
      concentration in the leaking brine (default: -8.610)

    * **Benzene_kd** [|L/kg|] (-4.5 to 0.69) - benzene distribution coefficient
      (default: -3.571)

    * **Benzene_decay** [|1/s|] (-6.1 to 0) - benzene degradation constant
      (default: -2.732)

    * **PAH_brine** [|mol/L|] (-10 to -4.1) - naphthalene concentration in the
      leaking brine (default: -7.118)

    * **PAH_kd** [|L/kg|] (-3.1 to 1.98) - naphthalene distribution coefficient
      (default: -0.985)

    * **PAH_decay** [|1/s|] (-6.45 to 0) - naphthalene degradation constant
      (default: -3.371)

    * **phenol_brine** [|mol/L|] (-10 to -3.7) - phenol concentration in the leaking
      brine (default: -6.666)

    * **phenol_kd** [|L/kg|] (-6 to 0.15) - phenol distribution coefficient
      (default: -1.342)

    * **phenol_decay** [|1/s|] (-5.62999999999999 to 0) - phenol degradation
      constant (default: -3.546)

    * **porositySand** [-] (0.25 to 0.5) - porosity (default: 0.375)

    * **densitySand** [|kg/m^3|] (1500 to 2500) - density (default: 2165.953)

    * **VG_mSand** [-] (0.52 to 0.79) - Van Genuchten parameter m (default: 0.627)

    * **VG_alphaSand** [|1/m|] (-4.69 to -3.81) - Van Genuchten parameter alpha
      (default: -4.341)

    * **permeabilitySand** [|log10| |m^2|] (-14 to -10) - permeability of sand
      (default: -12.430)

    * **Clbrine** [|mol/L|] (-2 to 0.73) - chlorine concentration in the leaking
      brine (default: -0.339)

    * **calcitevol** [-] (0.035 to 0.2) - volume fraction of calcite
      (default: 0.165)

    * **V_goethite** [-] (0 to 0.2) - volume fraction of goethite (default: 0.004)

    * **V_illite** [-] (0 to 0.3) - volume fraction of illite (default: 0.006)

    * **V_kaolinite** [-] (0 to 0.2) - volume fraction of kaolinite (default: 0.004)

    * **V_smectite** [-] (0 to 0.5) - volume fraction of smectite (default: 0.010)

    Component model dynamic inputs:

    * **co2_rate** [|kg/s|] (0 to 0.5) - |CO2| rate (default: 0.25005)

    * **co2_mass** [|kg|] (169.824 to 1.142e+9) - cumulative |CO2| mass
      (default: 4.405e+5)

    * **brine_rate** [|kg/s|] (0 to 0.075) - brine rate (default: 0.0375)

    * **brine_mass** [|kg|] (3.162e+4 to 1.330e+8) - cumulative brine mass
      (default: 2.051e+6)

    * **simtime** [years] (1 to 200) - simulation time (default: 1)

    Possible observations from the Alluvium Aquifer component are:

    * **TDS_volume** [|m^3|] - volume of aquifer above TDS threshold in mg/L
      (TDS > 1300 mg/L)

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
      (naphthalene > 0.2 |microg/L|)

    * **Phenol_volume** [|m^3|] - volume of aquifer above phenol threshold
      (phenol  > 0.003 |microg/L|).
    """
    def __init__(self, name, parent):
        """
        Constructor method of AlluviumAquifer class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :returns: AlluviumAquifer class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25}   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'AlluviumAquifer'

        # Initiate solution object
        self.sol = aarom.Solution()

        self.add_default_par('sandFraction', value=0.422)
        self.add_default_par('correlationLengthX', value=361.350)
        self.add_default_par('correlationLengthZ', value=17.382)
        self.add_default_par('permeabilityClay', value=-16.340)
        self.add_default_par('NaMolality', value=0.121)
        self.add_default_par('PbMolality', value=-5.910)
        self.add_default_par('benzeneMolality', value=0.109)
        self.add_default_par('tMitigation', value=87.914)
        self.add_default_par('CEC', value=32.073)
        self.add_default_par('Asbrine', value=-5.397)
        self.add_default_par('Babrine', value=-3.397)
        self.add_default_par('Cdbrine', value=-8.574)
        self.add_default_par('Pbbrine', value=-7.719)
        self.add_default_par('Benzene_brine', value=-8.610)
        self.add_default_par('Benzene_kd', value=-3.571)
        self.add_default_par('Benzene_decay', value=-2.732)
        self.add_default_par('PAH_brine', value=-7.118)
        self.add_default_par('PAH_kd', value=-0.985)
        self.add_default_par('PAH_decay', value=-3.371)
        self.add_default_par('phenol_brine', value=-6.666)
        self.add_default_par('phenol_kd', value=-1.342)
        self.add_default_par('phenol_decay', value=-3.546)
        self.add_default_par('porositySand', value=0.468)
        self.add_default_par('densitySand', value=2165.953)
        self.add_default_par('VG_mSand', value=0.627)
        self.add_default_par('VG_alphaSand', value=-4.341)
        self.add_default_par('permeabilitySand', value=-12.430)
        self.add_default_par('Clbrine', value=-0.339)
        self.add_default_par('calcitevol', value=0.165)
        self.add_default_par('V_goethite', value=0.004)
        self.add_default_par('V_illite', value=0.006)
        self.add_default_par('V_kaolinite', value=0.004)
        self.add_default_par('V_smectite', value=0.010)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['sandFraction'] = [0.35, 0.65]
        self.pars_bounds['correlationLengthX'] = [200, 2500]
        self.pars_bounds['correlationLengthZ'] = [0.5, 25.0]
        self.pars_bounds['permeabilityClay'] = [-18.0, -15.0]
        self.pars_bounds['NaMolality'] = [-3.0, 1.0]
        self.pars_bounds['PbMolality'] = [-8.5, -5.0]
        self.pars_bounds['benzeneMolality'] = [0.0, 1.0]
        self.pars_bounds['tMitigation'] = [1.0, 200.0]
        self.pars_bounds['CEC'] = [0.1, 40.0]
        self.pars_bounds['Asbrine'] = [-9.0, -5.0]
        self.pars_bounds['Babrine'] = [-5.1, -2.3]
        self.pars_bounds['Cdbrine'] = [-9.0, -6.0]
        self.pars_bounds['Pbbrine'] = [-8.5, -5.0]
        self.pars_bounds['Benzene_brine'] = [-8.8927, -4.8927]
        self.pars_bounds['Benzene_kd'] = [-4.5, 0.69]
        self.pars_bounds['Benzene_decay'] = [-6.1, 0.0]
        self.pars_bounds['PAH_brine'] = [-10.0, -4.1]
        self.pars_bounds['PAH_kd'] = [-3.1, 1.98]
        self.pars_bounds['PAH_decay'] = [-6.45, 0]
        self.pars_bounds['phenol_brine'] = [-10.0, -3.7]
        self.pars_bounds['phenol_kd'] = [-6.0, 0.15]
        self.pars_bounds['phenol_decay'] = [-5.62999999999999, 0.0]
        self.pars_bounds['porositySand'] = [0.25, 0.5]
        self.pars_bounds['densitySand'] = [1500, 2500]
        self.pars_bounds['VG_mSand'] = [0.52, 0.79]
        self.pars_bounds['VG_alphaSand'] = [-4.69, -3.81]
        self.pars_bounds['permeabilitySand'] = [-14.0, -10.0]
        self.pars_bounds['Clbrine'] = [-2.0, 0.73]
        self.pars_bounds['calcitevol'] = [0.035, 0.2]
        self.pars_bounds['V_goethite'] = [0, 0.2]
        self.pars_bounds['V_illite'] = [0, 0.3]
        self.pars_bounds['V_kaolinite'] = [0, 0.2]
        self.pars_bounds['V_smectite'] = [0, 0.5]

        # Define dictionary of temporal data limits
        self.temp_data_bounds = dict()
        self.temp_data_bounds['co2_rate'] = ['CO2 rate', 0., 0.5]
        self.temp_data_bounds['brine_rate'] = ['brine rate', 0., 0.075]
        self.temp_data_bounds['co2_mass'] = ['CO2 mass', 2.23, 9.058]    # Given in |log10| kg, so 169.8 kg to 1.14e+9 kg
        self.temp_data_bounds['brine_mass'] = ['brine mass', 4.5, 8.124] # Given in |log10| kg, so 3.16e+4 kg to 1.33e+8 kg
        self.temp_data_bounds['simtime'] = ['simulation time', 1, 200]

        # Define observations names
        self.out_labels = ['TDS', 'pH', 'As', 'Pb', 'Cd', 'Ba',
                           'Benzene', 'Naphthalene', 'Phenol']
        self.out_labels = [val+'_volume' for val in self.out_labels]

        debug_msg = 'AlluviumAquifer component created with name {}'.format(name)
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
                if (val < self.pars_bounds[key][0]) or (
                        val > self.pars_bounds[key][1]):
                    warn_msg = ''.join([
                        'Parameter {} of AlluviumAquifer component {} ',
                        'is out of boundaries.']).format(key, self.name)
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of a AlluviumAquifer component {}.']).format(key, self.name)
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
                    'Temporal input {} of AlluviumAquifer component {} ',
                    'is outside the model range [{}, {}] at time t = {}'])
                if key in ['co2_mass', 'brine_mass']:
                    warn_msg = warn_msg.format(
                        self.temp_data_bounds[key][0].lower(), self.name,
                        10**self.temp_data_bounds[key][1],
                        10**self.temp_data_bounds[key][2], time)
                else:
                    warn_msg = warn_msg.format(
                        self.temp_data_bounds[key][0].lower(), self.name,
                        self.temp_data_bounds[key][1],
                        self.temp_data_bounds[key][2], time)
                logging.warning(warn_msg)

    def normalize_temp_inpar(self, temp_inputs):
        """
        Normalize (0-1) temporal input parameters

        :param temp_inputs: temporal input data of component model
        :type temp_inputs: dict

        """
        k, v = list(temp_inputs.keys())[0], list(temp_inputs.values())[0]
        nv = (v - self.temp_data_bounds[k][1])/(
            self.temp_data_bounds[k][2] - self.temp_data_bounds[k][1])
        return  nv

    def normalize_inpar(self, inputs):
        """
        Normalize (0-1) input parameters

        :param inputs: input data of component model
        :type inputs: dict
        """
        k, v = list(inputs.keys())[0], list(inputs.values())[0]
        nv = (v - self.pars_bounds[k][0])/(
            self.pars_bounds[k][1] - self.pars_bounds[k][0])
        return  nv

    def simulation_model(self, p, co2_rate=0.25, brine_rate=0.0375,
                         co2_mass=10**5.644, brine_mass=10**6.312,
                         time_point=365.25):

        """
        Return volume of impacted aquifer based on several metrics.

        :param p: input parameters of alluvium aquifer model
        :type p: dict

        :param co2_rate: CO2 leakage rate, kg/s
        :type co2_rate: [float]

        :param brine_mass: Brine leakage rate, kg/s
        :type brine_mass: [float]

        :param co2_mass: Cumulative CO2 mass leaked, kg
        :type co2_mass: [float]

        :param brine_mass: Cumulative brine mass leaked, kg
        :type brine_mass: [float]

        :param time_point: time point in days at which the leakage rates are
            to be calculated; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param time_step: difference between the current and previous
            time points in days; by default, its value is 365.25 (1 year in days)
        :type time_step: float

        :param kwargs: additional keyword arguments for the model
        :type kwargs: dict (if provided)

        :returns: vol - dictionary of observations of carbonate aquifer impacted volumes
            model; keys:
            ['TDS_volume', 'pH_volume', 'As_volume', 'Pb_volume', 'Cd_volume',
            'Ba_volume', 'Benzene_volume', 'Naphthalene_volume', 'Phenol_volume']

        """
        # Check assign defailt values
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        actual_p.update(p)

        # Normalize input paramters
        for key, val in actual_p.items():
            p[key] = self.normalize_inpar(dict([(key, val)]))

        # Check whether co2, brine, rates and mass inputs satisfy the model requirements
        # Conditional statement used to supressess excessive error messsages
        if co2_mass <= 0:
            co2_mass = 1e-12
        if brine_mass <= 0:
            brine_mass = 1e-12

        if ((time_point/365.25 >= 1) and (
                np.log10(co2_mass) >= self.temp_data_bounds['co2_mass'][1])) or (
                    (time_point/365.25 >= 1)  and (
                        np.log10(brine_mass) >= self.temp_data_bounds['brine_mass'][1])):
            self.check_temporal_inputs(time_point, dict(list(zip(
                ['co2_rate', 'brine_rate', 'co2_mass', 'brine_mass', 'simtime'],
                [co2_rate, brine_rate, np.log10(co2_mass),
                 np.log10(brine_mass), time_point/365.25]))))

        #apply nominal value for brine/bo2
        p['norm_co2_rate'] = self.normalize_temp_inpar(
            dict([("co2_rate", co2_rate)]))
        p['norm_brine_rate'] = self.normalize_temp_inpar(
            dict([("brine_rate", brine_rate)]))
        p['norm_co2_mass'] = self.normalize_temp_inpar(
            dict([("co2_mass", np.log10(co2_mass))]))
        p['norm_brine_mass'] = self.normalize_temp_inpar(
            dict([("brine_mass", np.log10(brine_mass))]))

        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # normalize simtime
        actual_p['simtime'] = self.normalize_temp_inpar(
            dict([("simtime", time_point/365.25)]))

        oparam = np.array([
            actual_p['norm_co2_rate'], actual_p['norm_co2_mass'], actual_p['sandFraction'],
            actual_p['correlationLengthX'], actual_p['correlationLengthZ'],
            actual_p['permeabilityClay'], actual_p['NaMolality'], actual_p['PbMolality'],
            actual_p['benzeneMolality'], actual_p['tMitigation'], actual_p['norm_brine_rate'],
            actual_p['norm_brine_mass'], actual_p['simtime'], actual_p['CEC'], actual_p['Asbrine'],
            actual_p['Babrine'], actual_p['Cdbrine'], actual_p['Pbbrine'],
            actual_p['Benzene_brine'], actual_p['Benzene_kd'], actual_p['Benzene_decay'],
            actual_p['PAH_brine'], actual_p['PAH_kd'], actual_p['PAH_decay'],
            actual_p['phenol_brine'], actual_p['phenol_kd'], actual_p['phenol_decay'],
            actual_p['porositySand'], actual_p['densitySand'], actual_p['VG_mSand'],
            actual_p['VG_alphaSand'], actual_p['permeabilitySand'], actual_p['Clbrine'],
            actual_p['calcitevol'], actual_p['V_goethite'], actual_p['V_illite'],
            actual_p['V_kaolinite'], actual_p['V_smectite']])

        if ((time_point/365.25 >= 1) and (
                np.log10(co2_mass) >= self.temp_data_bounds['co2_mass'][1])) or (
                    (time_point/365.25 >= 1)  and (
                        np.log10(brine_mass) >= self.temp_data_bounds['brine_mass'][1])):

            self.sol.find(oparam)
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
    time_array = 365.25*np.arange(0.0, 2.0)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add alluvium aquifer model object and define parameters
    aa = sm.add_component_model_object(AlluviumAquifer(name='aa', parent=sm))

    aa.add_par('sandFraction', value=0.422)
    aa.add_par('correlationLengthX', value=361.350)
    aa.add_par('correlationLengthZ', value=17.382)
    aa.add_par('permeabilityClay', value=-16.340)
    aa.add_par('NaMolality', value=0.121)
    aa.add_par('PbMolality', value=-5.910)
    aa.add_par('benzeneMolality', value=0.109)
    aa.add_par('tMitigation', value=87.914)
    aa.add_par('CEC', value=32.073)
    aa.add_par('Asbrine', value=-5.397)
    aa.add_par('Babrine', value=-3.397)
    aa.add_par('Cdbrine', value=-8.574)
    aa.add_par('Pbbrine', value=-7.719)
    aa.add_par('Benzene_brine', value=-8.610)
    aa.add_par('Benzene_kd', value=-3.571)
    aa.add_par('Benzene_decay', value=-2.732)
    aa.add_par('PAH_brine', value=-7.118)
    aa.add_par('PAH_kd', value=-0.985)
    aa.add_par('PAH_decay', value=-3.371)
    aa.add_par('phenol_brine', value=-6.666)
    aa.add_par('phenol_kd', value=-1.342)
    aa.add_par('phenol_decay', value=-3.546)
    aa.add_par('porositySand', value=0.468)
    aa.add_par('densitySand', value=2165.953)
    aa.add_par('VG_mSand', value=0.627)
    aa.add_par('VG_alphaSand', value=-4.341)
    aa.add_par('permeabilitySand', value=-12.430)
    aa.add_par('Clbrine', value=-0.339)
    aa.add_par('calcitevol', value=0.165)
    aa.add_par('V_goethite', value=0.004)
    aa.add_par('V_illite', value=0.006)
    aa.add_par('V_kaolinite', value=0.004)
    aa.add_par('V_smectite', value=0.010)

    aa.model_kwargs['co2_rate'] = 1.90e-02    # kg/s
    aa.model_kwargs['co2_mass'] = 6.00e+05    # kg
    aa.model_kwargs['brine_rate'] = 4.62e-03  # kg/s
    aa.model_kwargs['brine_mass'] = 1.46e+05  # kg

    # Add observations (output) from the alluvium aquifer model
    # When add_obs method is called, system model automatically
    # (if no other options of the method is used) adds a list of observations
    # with names aa.obsnm_0, aa.obsnm_1,.... The final index is determined by the
    # number of time points in system model time_array attribute.
    # For more information, see the docstring of add_obs of ComponentModel class.

    aa.add_obs('TDS_volume')
    aa.add_obs('pH_volume')
    aa.add_obs('As_volume')
    aa.add_obs('Pb_volume')
    aa.add_obs('Cd_volume')
    aa.add_obs('Ba_volume')

    # Run the system model
    sm.forward()

    # Print the observations
    # Since the observations at the particular time points are different variables,
    # method collect_observations_as_time_series creates lists of
    # values of observations belonging to a given component (e.g. ca) and having the same
    # common name (e.g. 'pH', 'TDS', etc) but differing in indices.
    # More details are given in the docstring and documentation to the method
    # collect_observations_as_time_series of SystemModel class.
    print('TDS:', sm.collect_observations_as_time_series(aa, 'TDS_volume'))
    print('pH:', sm.collect_observations_as_time_series(aa, 'pH_volume'))
    print('As:', sm.collect_observations_as_time_series(aa, 'As_volume'))
    print('Pb:', sm.collect_observations_as_time_series(aa, 'Pb_volume'))
    print('Cd:', sm.collect_observations_as_time_series(aa, 'Cd_volume'))
    print('Ba:', sm.collect_observations_as_time_series(aa, 'Ba_volume'))

    # Expected output
    # ('TDS', array([       0.        ,  6878026.69860782]))
    # ('pH', array([  0.00000000e+00,   5.74256378e+09]))
    # ('As', array([  0.00000000e+00,   3.64713147e+08]))
    # ('Pb', array([  0.00000000e+00,   6.26222833e+09]))
    # ('Cd', array([  0.00000000e+00,   2.14112532e+09]))
    # ('Ba', array([  0.00000000e+00,   6.04804066e+09]))
