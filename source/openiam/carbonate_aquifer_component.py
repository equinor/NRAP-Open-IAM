# -*- coding: utf-8 -*-
import os
import sys
import logging
import ctypes
from sys import platform

import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import IAM_DIR, SystemModel, ComponentModel
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

from openiam.cfi.commons import process_parameters, process_dynamic_inputs

try:
    import components.aquifer as aqmodel
except ImportError:
    print('\nERROR: Unable to load ROM for Carbonate Aquifer component\n')
    sys.exit()


EDX_WINDOWS_CAAQ_LIB_HTTP = 'https://edx.netl.doe.gov/resource/fc330991-7a17-46b6-8995-17feefdf414c'
GITLAB_WINDOWS_CAAQ_LIB_HTTP = 'https://gitlab.com/NRAP/OpenIAM/-/blob/master/source/components/aquifer/carbonate/carbonate.dll'


class CarbonateAquifer(ComponentModel):
    """
    The Carbonate Aquifer component model is a reduced-order model that can be
    used to predict the impact that carbon dioxide (|CO2|) and brine leaks from a
    |CO2| storage reservoir might have on overlying aquifers. The model predicts
    the size of “impact plumes” according to nine water quality metrics, see
    :cite:`RN1087`, :cite:`RN631`, :cite:`RN1605`.

    Although the Carbonate Aquifer model was developed using site-specific data
    from the Edwards aquifer, the model accepts aquifer characteristics as
    variable inputs and, therefore, may have more broad applicability. Careful
    consideration should be given to the hydrogeochemical character of the
    aquifer before using this model at a new site. Guidelines and examples are
    presented in :cite:`RN1635`.

    The size of “impact plumes” are calculated using two alternative
    definitions of “impact” which should be selected by user: 1) changes that cause an
    exceedance of a drinking water standard or maximum contaminant level (MCL);
    and 2) changes that are above and beyond “natural background variability” in
    the aquifer, :cite:`RN1608`.

    Component model input definitions:

    * **ithresh** [-] (1 or 2) - threshold, either 1: MCL or 2: No-impact
      (default: 2)

    * **rmin** [|m|] (0 to 100) - maximum distance between leaks for them to be
      considered one leak (default: 15)

    * **perm_var** [|log10| |m^4|] (0.017 to 1.89) - logarithm of permeability
      variance (default: 0.9535)

    * **corr_len** [|m|] (1 to 3.95) - correlation length (default: 2.475)

    * **aniso** [-] (1.1 to 49.1) - anisotropy factor: ratio of horizontal to
      vertical permeability (default: 25.1)

    * **mean_perm** [|log10| |m^2|] (-13.8 to -10.3) - logarithm of mean
      permeability (default: -12.05)

    * **hyd_grad** [-] (2.88e-4 to 1.89e-2) - horizontal hydraulic gradient
      (default: 9.59e-03)

    * **calcite_ssa** [|m^2/g|] (0 to 1.0e-2) - calcite surface area
      (default: 5.5e-03)

    * **organic_carbon** [-] (0 to 1.0e-2) - organic carbon volume fraction
      (default: 5.5e-03)

    * **benzene_kd** [|log10| K_oc] (1.49 to 1.73) - benzene distribution
      coefficient (default: 1.61)

    * **benzene_decay** [|log10| day] (0.15 to 2.84) - benzene decay constant
      (default: 0.595)

    * **nap_kd** [|log10| K_oc] (2.78 to 3.18) - naphthalene distribution
      coefficient (default: 2.98)

    * **nap_decay** [|log10| day] (-0.85 to 2.04) - naphthalene decay constant
      (default: 0.595)

    * **phenol_kd** [|log10| K_oc] (1.21 to 1.48) - phenol distribution coefficient
      (default: 1.35)

    * **phenol_decay** [|log10| day] (-1.22 to 2.06) - phenol decay constant
      (default: 0.42)

    * **cl** [|log10| molality] (0.1 to 6.025) - brine salinity (default: 0.776)

    * **logf** [-] (0 or 1) - type of transform of output plume volume;
      0: linear, 1: log (default: 0)

    * **aqu_thick** [|m|] (100 to 500) - aquifer thickness (default: 300);
      *linked to Stratigraphy*

    Component model dynamic inputs:

    * **brine_rate** [|kg/s|] (0 to 0.075) - brine rate

    * **brine_mass** [|kg|] (0 to 2.0e+8) - cumulative brine mass

    * **co2_rate** [|kg/s|] (0 to 0.5) - |CO2| rate

    * **co2_mass** [|kg|] (0 to 2.0e+9) - cumulative |CO2| mass.

    Possible observations from the Carbonate Aquifer component are:

    * **pH_volume** [|m^3|] - volume of aquifer below pH threshold

    * **Flux** [|kg/s|] - |CO2| leakage rate to atmosphere

    * **dx** [|m|] - length of impacted aquifer volume in x-direction

    * **dy** [|m|] - width of impacted aquifer volume in y-direction

    * **TDS_volume** [|m^3|] - volume of aquifer above TDS threshold in |millig/L|

    * **As_volume** [|m^3|] - volume of aquifer above arsenic threshold in |microg/L|

    * **Pb_volume** [|m^3|] - volume of aquifer above lead threshold in |microg/L|

    * **Cd_volume** [|m^3|] - volume of aquifer above cadmium threshold in |microg/L|

    * **Ba_volume** [|m^3|] - volume of aquifer above barium threshold in |microg/L|

    * **Benzene_volume** [|m^3|] - volume of aquifer above benzene threshold

    * **Naphthalene_volume** [|m^3|] - volume of aquifer above naphthalene threshold

    * **Phenol_volume** [|m^3|] - volume of aquifer above phenol threshold.

    """

    # Attributes for system connections
    system_params = ['{aquifer_name}Thickness']
    system_collected_inputs = {'co2_rate': 'CO2_{aquifer_name}',
                               'brine_rate': 'brine_{aquifer_name}',
                               'co2_mass': 'mass_CO2_{aquifer_name}',
                               'brine_mass': 'mass_brine_{aquifer_name}'}
    adapters = ['RateToMass']
    needsXY = True

    def __init__(self, name, parent, header_file_dir=None):
        """
        Constructor method of CarbonateAquifer class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param header_file_dir: name of directory with Fortran DLL
        :type header_file_dir: str

        :returns: CarbonateAquifer class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25}   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'CarbonateAquifer'

        # Determine whether library is in the default location
        if header_file_dir is None:
            self.header_file_dir = aqmodel.__path__[0]
        else:
            self.header_file_dir = header_file_dir

        # Specify name of and path to dynamic library
        if platform in ("linux", "linux2"):
            # linux
            library_name = "carbonate.so"
        elif platform == "darwin":
            # OS X
            library_name = "carbonate.dylib"
        elif platform == "win32":
            # Windows
            library_name = "carbonate.dll"

        self.library = os.sep.join([self.header_file_dir, 'carbonate', library_name])
        library_folder = os.sep.join([self.header_file_dir, 'carbonate'])

        if not self.model_data_check():
            if platform == "win32":
                error_msg = ''.join([
                    'Carbonate Aquifer component "{}" cannot be created as ',
                    'the required library file "{}" is missing in the folder\n {}.\n\n',
                    'For Windows OS the required library file "{}" can be downloaded ',
                    'from a corresponding folder on EDX here:\n{}\n or directly ',
                    'from GitLab here:\n{}.\nThe downloaded file should be placed into ',
                    'the folder\n{}.']).format(name, library_name, library_folder,
                                               library_name, EDX_WINDOWS_CAAQ_LIB_HTTP,
                                               GITLAB_WINDOWS_CAAQ_LIB_HTTP,
                                               library_folder)
            else:
                error_msg = ''.join([
                    'Carbonate Aquifer component "{}" cannot be created as ',
                    'the required library file "{}" is missing in the folder\n {}.\n\n',
                    'For non-Windows OS the required library file "{}" should ',
                    'be compiled using the files provided in the folder\n{}\nand ',
                    'placed into the same folder.']).format(name, library_name,
                                                          library_folder,
                                                          library_name, library_folder)
            logging.error(error_msg)
            sys.exit()

        # Set default parameters of the component model
        self.add_default_par('ithresh', value=2)
        self.add_default_par('rmin', value=15.0)
        self.add_default_par('perm_var', value=0.9535)
        self.add_default_par('corr_len', value=2.475)
        self.add_default_par('aniso', value=25.1)
        self.add_default_par('mean_perm', value=-12.05)
        self.add_default_par('aqu_thick', value=300.)
        self.add_default_par('hyd_grad', value=9.59e-3)
        self.add_default_par('calcite_ssa', value=5.5e-03)
        self.add_default_par('organic_carbon', value=5.5e-03)
        self.add_default_par('benzene_kd', value=1.61)
        self.add_default_par('benzene_decay', value=1.5)
        self.add_default_par('nap_kd', value=2.98)
        self.add_default_par('nap_decay', value=0.595)
        self.add_default_par('phenol_kd', value=1.35)
        self.add_default_par('phenol_decay', value=0.42)
        self.add_default_par('cl', value=0.776)
        self.add_default_par('logf', value=0)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['ithresh'] = [1, 2]
        self.pars_bounds['rmin'] = [0, 100]
        self.pars_bounds['perm_var'] = [0.017, 1.89]
        self.pars_bounds['corr_len'] = [1.0, 3.95]
        self.pars_bounds['aniso'] = [1.1, 49.1]
        self.pars_bounds['mean_perm'] = [-13.8, -10.3]
        self.pars_bounds['aqu_thick'] = [100., 500.]
        self.pars_bounds['hyd_grad'] = [2.88e-4, 1.89e-2]
        self.pars_bounds['calcite_ssa'] = [0, 1.e-2]
        self.pars_bounds['organic_carbon'] = [0, 1.e-2]
        self.pars_bounds['benzene_kd'] = [1.49, 1.73]
        self.pars_bounds['benzene_decay'] = [0.15, 2.84]
        self.pars_bounds['nap_kd'] = [2.78, 3.18]
        self.pars_bounds['nap_decay'] = [-0.85, 2.04]
        self.pars_bounds['phenol_kd'] = [1.21, 1.48]
        self.pars_bounds['phenol_decay'] = [-1.22, 2.06]
        self.pars_bounds['cl'] = [0.1, 6.025]
        self.pars_bounds['logf'] = [0, 1]

        # Define dictionary of temporal data limits
        self.temp_data_bounds = dict()
        self.temp_data_bounds['co2_rate'] = ['CO2 rate', 0., 0.5]
        self.temp_data_bounds['brine_rate'] = ['brine rate', 0., 0.075]
        self.temp_data_bounds['co2_mass'] = ['CO2 mass', 0., 2.0e+9]
        self.temp_data_bounds['brine_mass'] = ['brine mass', 0., 2.0e+8]

        self.output_labels = [
            'pH_volume', 'Flux', 'dx', 'dy', 'TDS_volume', 'As_volume',
            'Pb_volume', 'Cd_volume', 'Ba_volume', 'Benzene_volume',
            'Naphthalene_volume', 'Phenol_volume']

        debug_msg = 'CarbonateAquifer component created with name {}'.format(name)
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
                        'Parameter {} (value: {}) of CarbonateAquifer component {} ',
                        'is out of bounds.']).format(key, val, self.name)
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of a CarbonateAquifer component {}.']).format(key, self.name)
                logging.warning(warn_msg)

    def check_temporal_inputs(self, time, temp_inputs):
        """
        Check whether temporal data fall within specified boundaries.

        :param temp_inputs: temporal input data of component model
        :type temp_inputs: dict
        """
        for key, vals in temp_inputs.items():
            for val in vals:
                if ((val < self.temp_data_bounds[key][1]) or
                        (val > self.temp_data_bounds[key][2])):
                    warn_msg = ''.join([
                        'Temporal input {} (value: {}) of CarbonateAquifer component {} ',
                        'is outside the model range [{}, {}] at time t = {}']).format(
                            self.temp_data_bounds[key][0].lower(),
                            val, self.name,
                            self.temp_data_bounds[key][1],
                            self.temp_data_bounds[key][2], time)
                    logging.warning(warn_msg)

    # Carbonate aquifer model function
    def simulation_model(self, p, x=None, y=None, co2_rate=None, brine_rate=None,
                         co2_mass=None, brine_mass=None, time_point=365.25):
        """
        Return volume of impacted aquifer based on several metrics.

        :param p: input parameters of carbonate aquifer model
        :type p: dict

        :param x: horizontal coordinate of leaking well (m)
        :type x: [float]

        :param y: horizontal coordinate of leaking well (m)
        :type y: [float]

        :param co2_rate: |CO2| leakage rate, kg/s
        :type co2_rate: [float]

        :param brine_rate: brine leakage rate, kg/s
        :type brine_rate: [float]

        :param co2_mass: cumulative |CO2| mass leaked, kg
        :type co2_mass: [float]

        :param brine_mass: cumulative brine mass leaked, kg
        :type brine_mass: [float]

        :param time_point: time point in days at which the leakage rates are
            to be calculated; by default its value is 365.25 (1 year in days)
        :type time_point: float

        :returns: vol - dictionary of observations of carbonate aquifer impacted volumes
            model; keys:
            ['pH_volume', 'Flux', 'dx', 'dy', 'TDS_volume', 'As_volume',
            'Pb_volume', 'Cd_volume', 'Ba_volume', 'Benzene_volume',
            'Napthalene_volume', 'Phenol_volume']

        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}

        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Determine number of leaks
        # x,y,co2_rate,brine_rate,co2_mass,brine_mass must all be same length
        if x is not None:
            nleak = len(x)
        else:
            error_msg = ''.join([
                'Locations of the leaking wells are not provided for ',
                'Carbonate Aquifer component {}']).format(self.name)
            logging.error(error_msg)

        # Check whether pressure and saturation inputs satisfy the model requirements
        if (co2_mass is not None) and (brine_mass is not None):
            self.check_temporal_inputs(time_point, dict(list(zip(
                ['co2_rate', 'brine_rate'], [co2_rate, brine_rate]))))
        else:
            error_msg = ''.join([
                'CO2 and/or brine mass arguments of the Carbonate Aquifer ',
                'component model are not provided.'])
            logging.error(error_msg)

        # convert rates and masses to correct units
        co2_rate_conv = [val * 1000. for val in co2_rate] # kg/s -> g/s
        brine_rate_conv = [val * 1000. for val in brine_rate]
        co2_mass_conv = [val * 1e-6 for val in co2_mass] # kg -> kTon
        brine_mass_conv = [val * 1e-6 for val in brine_mass]

        # get time in days from keyword arguments, convert to years
        time = time_point/365.25

        # 1:MCL, 2:No-Impact
        ithresh = actual_p['ithresh']

        # This is not defined in AIM
        rmin = actual_p['rmin']

        # aquifer parameters
        aquifer = np.array([actual_p['perm_var'], actual_p['corr_len'],
                            actual_p['aniso'], actual_p['mean_perm'],
                            actual_p['aqu_thick']/100., actual_p['hyd_grad']])
        # ROM expects aquifer thickness in hundreds of meters

        # geochemical parameters
        chem = np.array([actual_p['calcite_ssa'], actual_p['organic_carbon'],
                         actual_p['benzene_kd'], actual_p['nap_kd'],
                         actual_p['phenol_kd'], actual_p['benzene_decay'],
                         actual_p['nap_decay'], actual_p['phenol_decay']])

        # log brine salinity
        cl = [actual_p['cl'] for i in range(nleak)]

        # log transform of output plume volume 0=linear, 1=log
        logf = actual_p['logf']*np.ones((12,), dtype=int)

        # Load the dynamic library
        external_lib = ctypes.cdll.LoadLibrary(os.path.join(
            os.getcwd(), self.library))

        # Specify Fortran function name
        functionName = "gw_plumes"

        # Get needed function as attribute of the library
        function = getattr(external_lib, functionName)

        # Define C classes
        INT = ctypes.c_int
        DOUBLE = ctypes.c_double
        EightDoubles = DOUBLE*8
        SixDoubles = DOUBLE*6
        TwelveDoubles = DOUBLE*12
        LEAKS = DOUBLE*nleak
        TwelveInts = INT*12

        # Set argument types for value and array pointers
        function.argtypes = [ctypes.POINTER(INT),
                             ctypes.POINTER(DOUBLE),
                             ctypes.POINTER(SixDoubles),
                             ctypes.POINTER(DOUBLE),
                             ctypes.POINTER(INT),
                             ctypes.POINTER(LEAKS),
                             ctypes.POINTER(LEAKS),
                             ctypes.POINTER(LEAKS),
                             ctypes.POINTER(LEAKS),
                             ctypes.POINTER(LEAKS),
                             ctypes.POINTER(LEAKS),
                             ctypes.POINTER(LEAKS),
                             ctypes.POINTER(TwelveDoubles),
                             ctypes.POINTER(EightDoubles),
                             ctypes.POINTER(TwelveInts)]
        function.restype = None

        # Define values of the input parameters that will be passed to Fortran function
        par_ithresh = INT(ithresh)
        par_rmin = DOUBLE(rmin)
        par_aquifer = SixDoubles(*aquifer)
        par_time = DOUBLE(time)
        par_nleak = INT(nleak)
        par_x = LEAKS(*x)
        par_y = LEAKS(*y)
        par_cl = LEAKS(*cl)
        par_co2rate = LEAKS(*co2_rate_conv)
        par_co2mass = LEAKS(*co2_mass_conv)
        par_brinerate = LEAKS(*brine_rate_conv)
        par_brinemass = LEAKS(*brine_mass_conv)
        par_chem = EightDoubles(*chem)
        par_logf = TwelveInts(*logf)

        out_vol = TwelveDoubles()    # initialize output array

        # Change to the directory where the dynamic library is located
        currentWorkDir = os.getcwd()
        targetWorkDir = os.path.join(os.getcwd(), self.header_file_dir)
        os.chdir(targetWorkDir)

        # Call the Fortran function
        function(par_ithresh, par_rmin, par_aquifer, par_time, par_nleak,
                 par_x, par_y, par_cl, par_co2rate, par_co2mass,
                 par_brinerate, par_brinemass, out_vol, par_chem, par_logf)

        # Return to original directory
        os.chdir(currentWorkDir)

        # Return a dictionary of the output volumes
        vol = dict(list(zip(self.output_labels, out_vol)))
        return vol

    def connect_with_system(self, component_data, name2obj_dict,
                            locations, system_adapters, **kwargs):
        """
        Code to add carbonate aquifer to system model for control file interface.

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
        system_collected_inputs = {'co2_rate': 'CO2_{aquifer_name}',
                                   'brine_rate': 'brine_{aquifer_name}',
                                   'co2_mass': 'mass_CO2_{aquifer_name}',
                                   'brine_mass': 'mass_brine_{aquifer_name}'}
        # Leaking wells locations
        self.model_kwargs['x'] = component_data['locX']
        self.model_kwargs['y'] = component_data['locY']

        # Process parameters data if any
        process_parameters(self, component_data, name2obj_dict)

        # Process dynamic inputs if any
        process_dynamic_inputs(
            self, component_data, array_like=True, check_second_dim=True,
            size=len(self.model_kwargs['x']),
            quantity_to_compare='Number of leaking wells')

        # Make model connections
        if 'Connection' in component_data:
            connections = np.array([component_data['Connection']]).flatten()
            if len(connections) == 1 and not locations[connections[0]]['well']:
                # Connect to seal or fault component
                # Get component providing outputs
                connector = name2obj_dict[connections[0]]
                # Add observations of the connector and link to input of the aquifer
                for key, sci_item in system_collected_inputs.items():
                    obs_nm = sci_item.format(aquifer_name='aquifer')
                    connector.add_obs_to_be_linked(obs_nm, obs_type='grid')
                    self.add_kwarg_linked_to_obs(
                        key, connector.linkobs[obs_nm], obs_type='grid')
            else:
                # collectors is a dictionary containing inputs collected into arrays
                collectors = self._parent.collectors
                for collect in collectors:  # e.g. collect can be 'co2_rate'
                    cdict = collectors[collect][self.name]
                    argument = cdict['argument'].format(
                        aquifer_name=component_data['AquiferName'])
                    # connections is a list of wellbore components connected to the aquifer
                    connections = cdict['Connection']

                    for connect in connections:
                        for ind in range(locations[connect]['number']):
                            cname = connect + '_{0:03}'.format(ind)
                            for ad_nm in system_adapters:
                                if (system_adapters[ad_nm]['Connection'] == cname) and (
                                        system_adapters[ad_nm]['AquiferName'] == \
                                            component_data['AquiferName']):
                                    aname = ad_nm

                            adapter = name2obj_dict[aname]
                            connector = name2obj_dict[cname]
                            if argument in adapter.linkobs:
                                cdict['data'].append(adapter.linkobs[argument])
                            else:
                                cdict['data'].append(connector.linkobs[argument])

                for sinput in self.system_collected_inputs:
                    self.add_kwarg_linked_to_collection(
                        sinput, collectors[sinput][self.name]['data'])
        # End Connection if statement

        # Take care of parameter aqu_thick (possibly) defined by stratigraphy component
        sparam = '{aq}Thickness'.format(aq=component_data['AquiferName'])
        strata = name2obj_dict['strata']
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


    @staticmethod
    def model_data_check():
        """
        Check whether required library file for Carbonate Aquifer component was
        downloaded/compiled and placed to the right folder.
        """

        check_flag = 0

        # Specify name of and path to dynamic library
        if platform in ("linux", "linux2"):
            # linux
            library_name = "carbonate.so"
        elif platform == "darwin":
            # OS X
            library_name = "carbonate.dylib"
        elif platform == "win32":
            # Windows
            library_name = "carbonate.dll"

        library = os.sep.join([
            IAM_DIR, 'source', 'components', 'aquifer', 'carbonate', library_name])

        if not os.path.isfile(library):
            check_flag = check_flag + 1

        return (check_flag == 0)


if __name__ == "__main__":
    # Create system model
    time_array = 365.25*np.arange(0.0, 2.0)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add carbonate aquifer model object and define parameters
    ca = sm.add_component_model_object(CarbonateAquifer(name='ca', parent=sm))
    ca.add_par('perm_var', value=1.24e+00)
    ca.add_par('corr_len', value=2.01E+00)
    ca.add_par('aniso', value=4.36e+01)
    ca.add_par('mean_perm', value=-1.14E+01)
    ca.add_par('aqu_thick', value=4.71E+02)
    ca.add_par('hyd_grad', value=1.63e-02)
    ca.add_par('calcite_ssa', value=4.13E-03)
    ca.add_par('cl', value=4.59E+00)
    ca.add_par('organic_carbon', value=5.15e-03)
    ca.add_par('benzene_kd', value=1.53e+00)
    ca.add_par('nap_kd', value=2.78e+00)
    ca.add_par('phenol_kd', value=1.29e+00)
    ca.add_par('benzene_decay', value=1.58e+00)
    ca.add_par('nap_decay', value=3.65e-01)
    ca.add_par('phenol_decay', value=2.61e-01)
    ca.model_kwargs['x'] = [0.]
    ca.model_kwargs['y'] = [0.]
    ca.model_kwargs['co2_rate'] = [1.90e-02]     # kg/s
    ca.model_kwargs['co2_mass'] = [6.00e+05]     # kg
    ca.model_kwargs['brine_rate'] = [4.62e-03]   # kg/s
    ca.model_kwargs['brine_mass'] = [1.46e+05]   # kg

    # Add observations (output) from the carbonate aquifer model
    # When add_obs method is called, system model automatically
    # (if no other options of the method is used) adds a list of observations
    # with names ca.obsnm_0, ca.obsnm_1,.... The final index is determined by the
    # number of time points in system model time_array attribute.
    # For more information, see the docstring of add_obs of ComponentModel class.
    ca.add_obs('pH_volume')
    ca.add_obs('Flux')
    ca.add_obs('dx')
    ca.add_obs('dy')
    ca.add_obs('TDS_volume')
    ca.add_obs('As_volume')
    ca.add_obs('Pb_volume')
    ca.add_obs('Cd_volume')
    ca.add_obs('Ba_volume')
    ca.add_obs('Benzene_volume')
    ca.add_obs('Naphthalene_volume')
    ca.add_obs('Phenol_volume')

    # Run the system model
    sm.forward()

    # Print the observations
    # Since the observations at the particular time points are different variables,
    # method collect_observations_as_time_series creates lists of
    # values of observations belonging to a given component (e.g. ca) and having the same
    # common name (e.g. 'pH', 'Flux', etc) but differing in indices.
    # More details are given in the docstring and documentation to the method
    # collect_observations_as_time_series of SystemModel class.
    print('pH:', sm.collect_observations_as_time_series(ca, 'pH_volume'))
    print('Flux:', sm.collect_observations_as_time_series(ca, 'Flux'))
    print('dx:', sm.collect_observations_as_time_series(ca, 'dx'))
    print('dy:', sm.collect_observations_as_time_series(ca, 'dy'))
    print('TDS:', sm.collect_observations_as_time_series(ca, 'TDS_volume'))
    print('As:', sm.collect_observations_as_time_series(ca, 'As_volume'))
    print('Pb:', sm.collect_observations_as_time_series(ca, 'Pb_volume'))
    print('Cd:', sm.collect_observations_as_time_series(ca, 'Cd_volume'))
    print('Ba:', sm.collect_observations_as_time_series(ca, 'Ba_volume'))
    print('Benzene:',
          sm.collect_observations_as_time_series(ca, 'Benzene_volume'))
    print('Naphthalene:',
          sm.collect_observations_as_time_series(ca, 'Naphthalene_volume'))
    print('Phenol:',
          sm.collect_observations_as_time_series(ca, 'Phenol_volume'))

    # Expected output
    #    pH:          [0. 0.]
    #    Flux:        [0. 0.]
    #    dx:          [97.81492584 97.81492584]
    #    dy:          [61.16387968 61.16387968]
    #    TDS:         [18940957.05166382 18925088.43222356]
    #    As:          [318247.57621476 318247.57621476]
    #    Pb:          [0. 0.]
    #    Cd:          [2301435.20549617 2301637.01148604]
    #    Ba:          [0. 0.]
    #    Benzene:     [36434231.34681176 36434231.34681176]
    #    Naphthalene: [391423.9405283 391423.9405283]
    #    Phenol:      [30248015.99834153 30248015.99834153]
