# -*- coding: utf-8 -*-
import ctypes
import logging
import os
import sys
from sys import platform

import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from openiam import SystemModel, ComponentModel, IAM_DIR
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

from openiam.cfi.commons import process_parameters, process_dynamic_inputs

try:
    import components.atmosphere as atmModel
except ImportError:
    print('\nERROR: Unable to load ROM for Atmospheric ROM component\n')
    sys.exit()


EDX_WINDOWS_ATMROM_LIB_HTTP = 'https://edx.netl.doe.gov/resource/fc330991-7a17-46b6-8995-17feefdf414c'
GITLAB_WINDOWS_ATMROM_LIB_HTTP = 'https://gitlab.com/NRAP/OpenIAM/-/blob/master/source/components/atmosphere/atmdisrom.dll'


class AtmosphericROM(ComponentModel):
    """
    The Atmospheric model is meant to be used for performing scoping studies for |CO2|
    dispersion after leakage out of the ground. The employed method is an
    extension of the nomograph approach of Britter and McQuaid (1988)
    :cite:`BritterMcQuaid1988` developed for estimating dense gas plume length from
    a single or multiple leakage sources. The method is very fast
    and, therefore, amenable to general system-level geologic carbon sequestration
    (GCS) risk assessment. The method is conservative: it assumes the wind
    could be from any direction and handles multiple sources by a simple
    superposition approach :cite:`ZHANG2016323`. A user's manual for the standalone
    model is available at :cite:`Zhang2016`.

    The model is intended to be used for large |CO2| leakage rates
    (e.g., leakage from an open wellbore).  It may not be suitable for
    very small leakage rate, as, in general, small release rates
    (e.g., less than 1.0e-5 kg/s) do not form a dense gas release due to ambient
    mixing. The inputs to the model are leakage rate(s) from leaky well(s),
    location(s) of leaky well(s), ambient conditions (wind speed), and receptor
    locations (home or business locations where people are present). The outputs
    from the model are flags at receptors indicating whether the |CO2| concentration
    at the location exceeds a pre-defined critical value, and the critical
    downwind distance from the sources.

    Within the control file interface, receptor locations can be specified with
    the ``receptors`` keyword argument assigned a full path (including a name) to
    a csv file containing x- and y-coordinates of the receptors.  Alternatively,
    the ``x_receptor`` and ``y_receptor`` keywords can be assigned a list of x- and
    y-coordinates of the receptors, respectively. In the NRAP-Open-IAM control file,
    the type name for the Atmospheric model component is ``AtmosphericROM``.

    Component model input definitions:

    * **T_amb** [|C|] (5 to 40) - ambient temperature (default: 15)

    * **P_amb** [atmosphere] (0.7 to 1.08) - ambient pressure (default: 1)

    * **V_wind** [|m/s|]  (1.e-10 to 20) - wind velocity (default: 5)

    * **C0_critical** [-] (0.002 to 0.1) - critical concentration (default: 0.01)

    * **T_source** [|C|] (5 to 50) - released |CO2| temperature (default: 15)

    * **x_receptor** [|m|] - x-coordinate of receptor

    * **y_receptor** [|m|] - y-coordinate of receptor

    Possible observations from the Atmospheric Model component are:

    * **outflag_r###** [-] - count of critical distances receptor is within
      from original leak points; here, ### is a receptor number starting at 000

    * **num_sources** [-] - number of sources. The possible maximum is a number
      of leakage points; could be less as leakage sources can potentially coalesce.

    * **x_new_s###** [|m|] - x-coordinate of leakage source; here ### is a source number
      starting at 000

    * **y_new_s###** [|m|] - y-coordinate of leakage source

    * **critical_distance_s###** [|m|] - critical downwind distance from each source.

    """
    def __init__(self, name, parent, header_file_dir=None):
        """
        Constructor method of AtmosphericROM class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param header_file_dir: name of directory with Fortran DLL
        :type header_file_dir: str

        :returns: AtmosphericROM class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25, 'time_step': 365.25} # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'AtmosphericROM'

        # Set up library
        if platform in ("linux", "linux2"):
            # linux
            library_name = "atmdisrom.so"
        elif platform == "darwin":
            # OS X
            library_name = "atmdisrom.dylib"
        elif platform == "win32":
            # Windows...
            library_name = "atmdisrom.dll"

        if header_file_dir is None:
            header_file_dir = atmModel.__path__[0]

        self.library = os.sep.join([header_file_dir, library_name])
        library_folder = header_file_dir

        if not self.model_data_check():
            if platform == "win32":
                error_msg = ''.join([
                    'Atmospheric ROM component "{}" cannot be created as the ',
                    'required library file "{}" is missing in the folder\n {}.\n',
                    'For Windows OS the required library file "{}" can be downloaded ',
                    'from a corresponding folder on EDX here:\n{}\n or directly ',
                    'from GitLab here:\n{}.\nThe downloaded file should be ',
                    'placed into the folder\n{}.']).format(
                        name, library_name, library_folder, library_name,
                        EDX_WINDOWS_ATMROM_LIB_HTTP, GITLAB_WINDOWS_ATMROM_LIB_HTTP,
                        library_folder)
            else:
                error_msg = ''.join([
                    'Atmospheric ROM component "{}" cannot be created as the ',
                    'required library file "{}" is missing in the folder\n {}.\n',
                    'For non-Windows OS the required library file "{}" should ',
                    'be compiled using the files provided in the folder\n{}\nand ',
                    'placed into the same folder.']).format(
                        name, library_name, library_folder, library_name, library_folder)
            logging.error(error_msg)
            sys.exit()

        # Placeholder for keyword arguments of the 'model' method:
        # to let the system model know that this component needs the specified keyword arguments
        self.model_kwargs['time_point'] = 365.25  # default value of 365.25 days
        self.model_kwargs['time_step'] = 365.25   # default value of 365.25 days

        # Set default values for user defined parameters of the component model
        self.add_default_par('T_amb', value=15.0)
        self.add_default_par('P_amb', value=1.0)
        self.add_default_par('V_wind', value=5.0)
        self.add_default_par('C0_critical', value=0.01)
        self.add_default_par('T_source', value=15.0)
        self.model_kwargs['x_receptor'] = [0.0]
        self.model_kwargs['y_receptor'] = [0.0]

        # Define dictionary of model parameters boundaries
        self.pars_bounds = dict()
        self.pars_bounds['T_amb'] = [5.0, 40.0]
        self.pars_bounds['P_amb'] = [0.7, 1.08]
        self.pars_bounds['V_wind'] = [1.0e-10, 20.0]
        self.pars_bounds['C0_critical'] = [0.002, 0.1]
        self.pars_bounds['T_source'] = [5.0, 50.0]

        debug_msg = 'AtmosphericROM component created with name {}'.format(name)
        logging.debug(debug_msg)

    # should not this be something written once but all component models can use?
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
                        'Parameter {} of AtmosphericROM component {} ',
                        'is out of bounds.']).format(key, self.name)
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of AtmosphericROM component {}.']).format(key, self.name)
                logging.warning(warn_msg)

    def simulation_model(self, p, x_coor=None, y_coor=None, co2_leakrate=None,
                         time_point=365.25, time_step=365.25,
                         x_receptor=None, y_receptor=None):

        """
        Return flag at receptors if |CO2| concentration is above critical values.

        :param p: dictionary of input parameters
        :type p: dict

        :param time_point: current time point in days (at which the output is
            to be provided); by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param time_step: difference between current and previous time points
            in days; by default, its value is 365.25 (1 year in days)
        :type time_step: float

        :param x_coor: horizontal x-coordinate of leaking well (m)
        :type x_coor: [float]

        :param y_coor: horizontal y-coordinate of leaking well (m)
        :type y_coor: [float]

        :param x_receptor: x-coordinate of receptor (m)
        :type x_receptor: [float]

        :param y_receptor: y-coordinate of receptor (m)
        :type y_receptor: [float]

        :param co2_leakrate: CO2 leakage rate, kg/s
        :type co2_leakrate: [float]

        :returns: dictionary of locations and flags at receptors; keys:['out_flag']
                  new source locations if CO2 leakage sources are combined;
                  keys:['x_new', 'y_new']; and critical downwind
                  distance: ['critical_distance']
        """
        msg = '{sr} model call input {p}, leak {leak}'.format(
            sr=self.name, p=p, leak=co2_leakrate)
        logging.debug(msg)

        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Set parameters
        T_amb = actual_p['T_amb']
        P_amb = actual_p['P_amb']
        V_wind = actual_p['V_wind']
        C0_critical = actual_p['C0_critical']
        T_source = actual_p['T_source']

        if (x_receptor is not None) and (y_receptor is not None):
            N_receptor = len(x_receptor)
        else:
            err_msg = 'Receptor locations are not provided.'
            logging.error(err_msg)

        input_array = np.array([T_amb, P_amb, V_wind,
                                C0_critical, T_source])
        Ninput = np.size(input_array)

        # determine number of leaks
        # x,y,co2_rate,brine_rate,co2_mass,brine_mass must all be same length
        if (x_coor is not None) and (y_coor is not None):
            nleak = len(x_coor)
        else:
            err_msg = 'Locations of CO2 source (leaking wells) are not provided.'
            logging.error(err_msg)

        # get time in days from keyword arguments, convert to years
        Nstep = int(time_point/time_step)

        # Specify function name
        functionName = "atmdisrom"

        # Define c classes
        INT = ctypes.c_int
        DOUBLE = ctypes.c_double
        LEAKS = DOUBLE*nleak
        RECEPTORS = DOUBLE*N_receptor
        NDOUBLE = DOUBLE*Ninput
        NINT = INT*N_receptor

        # Load DLL
        external_lib = ctypes.cdll.LoadLibrary(os.path.join(os.getcwd(), self.library))

        # Get needed function as attribute of the library
        function = getattr(external_lib, functionName)

        # Set argument types for value and pointers
        function.argtypes = [ctypes.POINTER(INT),  # nleak
                             ctypes.POINTER(INT),  # n_receptor
                             ctypes.POINTER(INT),  # n_step
                             ctypes.POINTER(NDOUBLE),  # input_array
                             ctypes.POINTER(LEAKS),  # x_coor
                             ctypes.POINTER(LEAKS),  # y_coor
                             ctypes.POINTER(RECEPTORS),  # x_receptor
                             ctypes.POINTER(RECEPTORS),	 # y_receptor
                             ctypes.POINTER(LEAKS),  # co2_leakrate
                             ctypes.POINTER(NINT),  # output_array
                             ctypes.POINTER(INT),  # num_source
                             ctypes.POINTER(LEAKS),  # combined source location, x
                             ctypes.POINTER(LEAKS),  # combined source location, y
                             ctypes.POINTER(LEAKS),  # critical downwind distance
                             ctypes.POINTER(INT)]  # message if using B&M method is valid
        function.restype = None

        par_nleak = INT(nleak)
        par_nreceptor = INT(N_receptor)
        par_nstep = INT(Nstep)
        par_inputarray = NDOUBLE(*input_array)
        par_x = LEAKS(*x_coor)
        par_y = LEAKS(*y_coor)
        par_co2rate = LEAKS(*co2_leakrate)
        par_xre = RECEPTORS(*x_receptor)
        par_yre = RECEPTORS(*y_receptor)

        out_flag = NINT()    # initialize output array
        num_source = INT()
        x_new = LEAKS()
        y_new = LEAKS()
        critical_distance = LEAKS()
        log_message = INT()

        function(par_nleak, par_nreceptor, par_nstep, par_inputarray,
                 par_x, par_y, par_xre, par_yre, par_co2rate, out_flag,
                 num_source, x_new, y_new, critical_distance, log_message)

        out = dict()

        if log_message == 1:
            logging.warning(''.join(['Dense gas criteria is not satisfied, ',
                                     'results from B&M may not be valid.']))
        if log_message == 2:
            logging.warning(''.join(['Alpha value in the B&M is out of valid ',
                                     'range, results may not be correct.']))
        if log_message == 3:
            logging.warning(''.join(['Dense gas criteria is not satisfied, ',
                                     'alpha value in the B&M is out of valid ',
                                     'range, results may not be correct.']))

        for i, outflag_i in enumerate(out_flag):
            # Extract value from the float type of output
            out['outflag_r{0:03}'.format(i)] = outflag_i

        out['num_sources'] = num_source.value

        for i, x_new_i in enumerate(x_new):
            out['x_new_s{0:03}'.format(i)] = x_new_i

        for i, y_new_i in enumerate(y_new):
            out['y_new_s{0:03}'.format(i)] = y_new_i
        for i, critical_distance_i in enumerate(critical_distance):
            out['critical_distance_s{0:03}'.format(i)] = critical_distance_i

        # Return output dictionary
        return out

    system_collected_inputs = {'co2_leakrate': 'CO2_atm'}

    def connect_with_system(self, component_data, name2obj_dict,
                            locations, adapters, **kwargs):
        """
        Add AtmosphericROM to system model for control file interface.

        :param component_data: Dictionary of component data including 'Connection'
        :type component_data: dict

        :param name2obj_dict: Dictionary to map component names to objects
        :type name2obj: dict

        :param locations: Dictionary with keys being the names of wellbore components
            added to the system and values being a dictionary with 3 keys:
            'number' (number of locations), 'coordx' (x-coordinates of locations)
            and 'coordy' (y-coordinates of locations).
        :type locations: dict

        :param adapters: Dictionary with keys being the names of adapters
            added to the system and values being a dictionary with 3 keys:
            'Type' (type of adapter), 'Connection' (name of wellbore component
            to which the given adapter is connected) and 'AquiferName' (name of
            aquifer for which the leakage rates of brine and CO2 from
            the connected wellbore component are transformed into the masses)
        :type adapters: dict

        :returns: None
        """
        # Leaking wells locations
        self.model_kwargs['x_coor'] = component_data['locX']
        self.model_kwargs['y_coor'] = component_data['locY']

        # Process parameters data if any
        process_parameters(self, component_data, name2obj_dict)

        # Process dynamic inputs if any
        process_dynamic_inputs(
            self, component_data, array_like=True, check_second_dim=True,
            size=len(self.model_kwargs['x_coor']),
            quantity_to_compare='Number of leaking wells')

        if ('receptors' in component_data) and (component_data['receptors']):
            xy_data = np.genfromtxt(
                os.sep.join([IAM_DIR, component_data['receptors']]),
                delimiter=',')
            self.model_kwargs['x_receptor'] = xy_data[:, 0]
            self.model_kwargs['y_receptor'] = xy_data[:, 1]
        elif ('x_receptor' in component_data) and (component_data['x_receptor']):
            self.model_kwargs['x_receptor'] = component_data['x_receptor']
            self.model_kwargs['y_receptor'] = component_data['y_receptor']

        # Make model connections
        if 'Connection' in component_data:
            collectors = self._parent.collectors
            for collect in collectors:
                cdict = collectors[collect][self.name]
                argument = cdict['argument']
                connections = cdict['Connection']
                for connect in connections:
                    for ind in range(locations[connect]['number']):
                        cname = connect + '_{0:03}'.format(ind)
                        connector = name2obj_dict[cname]
                        connector.add_obs_to_be_linked(argument)
                        cdict['data'].append(connector.linkobs[argument])

            for sinput in self.system_collected_inputs:
                self.add_kwarg_linked_to_collection(
                    sinput, collectors[sinput][self.name]['data'])

        # End Connection if statement
        if 'Outputs' in component_data:
            if 'outflag' in component_data['Outputs']:
                component_data['Outputs'].remove('outflag')
                for ir in range(len(self.model_kwargs['x_receptor'])):
                    component_data['Outputs'].append('outflag_r{0:03}'.format(ir))
            source_outputs = ['x_new', 'y_new', 'critical_distance']
            for source in source_outputs:
                if source in component_data['Outputs']:
                    component_data['Outputs'].remove(source)
                    for ip in range(len(component_data['locX'])):
                        component_data['Outputs'].append(source + '_s{0:03}'.format(ip))


    @staticmethod
    def model_data_check():
        """
        Check whether required library file for Atmospheric ROM component was
        downloaded/compiled and placed to the right folder.
        """

        check_flag = 0

        # Specify name of and path to dynamic library
        if platform in ("linux", "linux2"):
            # linux
            library_name = "atmdisrom.so"
        elif platform == "darwin":
            # OS X
            library_name = "atmdisrom.dylib"
        elif platform == "win32":
            # Windows...
            library_name = "atmdisrom.dll"

        library = os.sep.join([
            IAM_DIR, 'source', 'components', 'atmosphere', library_name])

        if not os.path.isfile(library):
            check_flag = check_flag + 1

        return (check_flag == 0)


if __name__ == "__main__":
    # Create system model
    time_array = 365.25*np.arange(0.0, 2.0)
    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # input parameters by user
    satm = sm.add_component_model_object(AtmosphericROM(name='satm', parent=sm))
    satm.add_par('T_amb', value=20.0)
    satm.add_par('P_amb', value=1.01E+00)
    satm.add_par('V_wind', value=5.0)
    satm.add_par('C0_critical', value=0.01)

    # input parameters passed by other models
    satm.model_kwargs['x_coor'] = [0.1, 100.2]
    satm.model_kwargs['y_coor'] = [0.2, 100.3]
    co2_leakrate = [[1.90e-01, 1.90e-01],[0.1, 0.1]]  # kg/s
    satm.add_dynamic_kwarg('co2_leakrate', co2_leakrate)

    satm.model_kwargs['x_receptor'] = [100., 10.]
    satm.model_kwargs['y_receptor'] = [110., 10.]

    n_receptors = len(satm.model_kwargs['x_receptor'])
    n_sources = len(co2_leakrate[0])
    # Add observations (output) from the atmospheric ROM model
    for ind in range(n_receptors):
        satm.add_obs('outflag_r{0:03}'.format(ind))

    for ind in range(n_sources):
        satm.add_obs('x_new_s{0:03}'.format(ind))
        satm.add_obs('y_new_s{0:03}'.format(ind))
        satm.add_obs('critical_distance_s{0:03}'.format(ind))

    satm.add_obs('num_sources')

    sm.forward()

    num_sources = sm.collect_observations_as_time_series(satm, 'num_sources')

    # Print the observations
    for ind in range(n_receptors):
        print('Receptor_flag_r{0:03}'.format(ind),
              satm.obs['outflag_r{0:03}_1'.format(ind)].sim)

    print('Number of sources', num_sources[0])
    for ind in range(num_sources[0]):
        print('New_Source_x{0:03}'.format(ind),
              satm.obs['x_new_s{0:03}_1'.format(ind)].sim)
        print('New_Source_y{0:03}'.format(ind),
              satm.obs['y_new_s{0:03}_1'.format(ind)].sim)
        print('New_Source_r{0:03}'.format(ind),
              satm.obs['critical_distance_s{0:03}_1'.format(ind)].sim)
