# -*- coding: utf-8 -*-
from math import sqrt
import sys
import os
import logging
import csv
import numpy as np
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import SystemModel, ComponentModel, IAM_DIR
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

from openiam.cfi.commons import process_parameters

try:
    import components.reservoir.theis.theis_reservoir_ROM as tresrom
except ImportError:
    print('\nERROR: Unable to load ROM for Theis Reservoir component\n')
    sys.exit()

# Units for injection rates and time in the control file interface and GUI.
INJECTION_DATA_UNITS = {'injRates': 'm^3/s^{-1}', 'injTimes': 'year'}


class TheisReservoir(ComponentModel):
    """
    The Theis Reservoir component model is an analytical model for pressure in the
    reservoir.

    In the NRAP-Open-IAM control file, the type name for the Theis Reservoir component is
    ``TheisReservoir``. The description of the component's parameters are
    provided below:

    * **initialPressure** [|Pa|] (8.0e+4 to 1.0e+7) - initial pressure at the top
      of the reservoir (default: 1.0e+6);

    * **reservoirThickness** [|m|] (1 to 1600) - thickness of reservoir (default: 50);
      *linked to Stratigraphy*

    * **logResPerm** [|log10| |m^2|] (-14 to -9) - base 10 logarithm of reservoir
      permeability (default: -12)

    * **reservoirPorosity** [-] (0.01 to 1) - porosity of reservoir (default: 0.3)

    * **compressibility** [|Pa^-1|] (1.0e-11 to 1.0e-6) - compressibility
      of porous media (default: 1.0e-10)

    * **CO2Density** [|kg/m^3|] (100 to 1500) - density of |CO2| phase
      (default: 479)

    * **brineDensity** [|kg/m^3|] (900 to 1500) - density of brine phase
      (default: 1000)

    * **brineViscosity** [|Pa*s|] (1.0e-4 to 5.0e-3) - viscosity of brine phase
      (default: 2.535e-3).

    Possible observation from the Theis Reservoir component is:

    * **pressure** [|Pa|] - pressure at top of the reservoir at the user
      defined location(s)

    * **CO2saturation** [-] - |CO2| saturation at the top of the reservoir at the
      user defined location(s)

    For compatibility with wellbore components in NRAP-Open-IAM
    observation CO2saturation of Theis Reservoir component has fixed values
    of 0 at each time point for any simulation.
    """
    def __init__(self, name, parent, injX=0., injY=0., locX=100., locY=100.,
                 injTimes=None, injRates=None):
        """
        Constructor method of TheisReservoir class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param injX: x-coordinate of the injector/producer location
        :type injX: float

        :param injY: y-coordinate of the injector/producer location
        :type injY: float

        :param locX: x-coordinate of the location for pressure to be calculated
        :type locX: float

        :param locY: y-coordinate of the location for pressure to be calculated
        :type locY: float

        :param injTimes: injection time(s) in days
        :type injTimes: list of floats

        :param injRates: injection rate(s) (+ injection, - withdrawal)
        :type injRates: list of floats

        :returns: TheisReservoir class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {}

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'TheisReservoir'

        # Set default parameters of the component model
        self.add_default_par('reservoirThickness', value=50.0)
        self.add_default_par('logResPerm', value=-12.0)
        self.add_default_par('reservoirPorosity', value=0.3)
        self.add_default_par('compressibility', value=1.0e-10)
        self.add_default_par('brineDensity', value=1000.0)
        self.add_default_par('CO2Density', value=479.0)
        self.add_default_par('brineViscosity', value=2.535e-3)
        self.add_default_par('initialPressure', value=1.0e+6)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['reservoirThickness'] = [1.0, 1600.0]
        self.pars_bounds['logResPerm'] = [-14.0, -9.0]
        self.pars_bounds['reservoirPorosity'] = [0.01, 1.0]
        self.pars_bounds['compressibility'] = [1.0e-11, 1.0e-6]
        self.pars_bounds['brineDensity'] = [900.0, 1500.0]
        self.pars_bounds['CO2Density'] = [100.0, 1500.0]
        self.pars_bounds['brineViscosity'] = [1.0e-4, 5.0e-3]
        self.pars_bounds['initialPressure'] = [8.0e+4, 1.0e+7]

        # Locations and injection setup
        self.process_injection_and_location_args(injX, injY, locX, locY,
                                                 injTimes, injRates)

        # Indicate to the system model that this component will be run only once:
        # at the very first time point
        self.run_frequency = 1

        debug_msg = 'TheisReservoir component created with name {}'.format(self.name)
        logging.debug(debug_msg)


    def process_injection_and_location_args(self, injX, injY, locX, locY,
                                            injTimes, injRates, setup=1):
        """
        Setup locations of injection and observation wells, as well as injection
        rates and injection times.
        injTimes are provided in days similar to system model time_array.
        """
        # Setup locations of injection wells
        if injX is not None:
            self.injX = np.asarray(injX).flatten()
        if injY is not None:
            self.injY = np.asarray(injY).flatten()
        if len(self.injX) != len(self.injY):
            err_msg = 'injTimes and injRates must be of the same length'
            logging.error(err_msg)
            raise ValueError(err_msg)

        # Setup locations of interest at which pressure is to be calculated
        if locX is not None:
            self.locX = locX
        if locY is not None:
            self.locY = locY

        self.num_inj_wells = len(self.injX)  # default is 1

        # Number of injection times coincide with the number of columns in array/matrix injRates
        if injTimes is not None:
            self.injTimes = np.asarray(injTimes)
            if self.injTimes.ndim == 1:
                self.injTimes = np.tile(self.injTimes, (self.num_inj_wells, 1))
        elif setup == 1:  # setting for the first time
            self.injTimes = np.tile(np.asarray(self._parent.time_array),
                                    (self.num_inj_wells, 1))

        # Specify default injection rate(s)(+ injection, - withdrawal) and time(s)
        # dimension of self.injRates is num_of_locs x num_inj_times
        if injRates is not None:
            self.injRates = np.asarray(injRates)
            # the same injection rates for all injection well locations or
            # for a single injection well
            if self.injRates.ndim == 1:
                self.injRates = np.tile(self.injRates, (self.num_inj_wells, 1))
            if self.injRates.shape != self.injTimes.shape:
                err_msg = 'injTimes and injRates must be of the same length'
                logging.error(err_msg)
                raise ValueError(err_msg)
        elif setup == 1:  # setting for the first time
            self.injRates = 0.1*np.ones(self.injTimes.shape)


    def check_injection_data(self, comp_name, component_data):
        """
        Code to read and check input for injection times and injection rates
        within the control file interface.

        :param component_data: Dictionary of component data
        :type component_data: dict

        :returns: None
        """
        injRatesInput = None
        if 'injRates' in component_data:
            injRatesInput = self.get_injection_data(
                comp_name, component_data, 'injRates')

        if injRatesInput is None:
            warning_msg = ''.join([
                'Argument injRates was not properly (or at all) setup for ',
                'the Theis Reservoir component {}. The default ',
                'injection rate of 0.1 {} will be used.']).format(comp_name, 'm^3/s^{-1}')
            logging.warning(warning_msg)

        injTimesInput = None
        if 'injTimes' in component_data:
            injTimesInput = self.get_injection_data(
                comp_name, component_data, 'injTimes')

        if injTimesInput is None:
            warning_msg = ''.join([
                'Argument injTimes was not properly (or at all) setup for ',
                'the Theis Reservoir component {}. The injection times will ',
                'be set to the time points defined in ModelParams ',
                'section.']).format(comp_name)
            logging.warning(warning_msg)

        return injTimesInput, injRatesInput


    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        debug_msg = 'Input parameters of component {} are {}'.format(self.name, p)
        logging.debug(debug_msg)

        for key, val in p.items():
            warn_msg = ''.join([
                'Parameter {} of TheisReservoir component {} ',
                'is out of boundaries.']).format(key, self.name)
            if key in self.pars_bounds:
                if (val < self.pars_bounds[key][0]) or (
                        val > self.pars_bounds[key][1]):
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of TheisReservoir component {}.']).format(key, self.name)
                logging.warning(warn_msg)


    def connect_with_system(self, component_data, name2obj_dict, *args, **kwargs):
        """
        Code to add Theis Reservoir to system model for control file interface.

        :param component_data: Dictionary of component data including injTimes and injRates
        :type component_data: dict

        :param name2obj_dict: Dictionary to map component names to objects
        :type name2obj: dict

        :returns: None
        """
        # Process parameters data if any
        process_parameters(self, component_data, name2obj_dict)

        # Get original component name
        orig_comp_name = self.name.split('_')[0]

        # Get locations at which pressure is to be calculated
        try:
            locX = component_data['coordx']
            locY = component_data['coordy']
        except KeyError:
            warn_msg = ''.join([
                'No (or incomplete) data about locations of interest ',
                'is provided for component {}. Default location (100, 100) ',
                'will be used']).format(orig_comp_name)
            logging.warning(warn_msg)
            locX = 100.0
            locY = 100.0

        # Get injection well locations
        try:
            injX = component_data['injX']
            injY = component_data['injY']
        except KeyError:
            warn_msg = ''.join([
                'No (or incomplete) data about injection well locations ',
                'is provided for component {}. Default location (0, 0) ',
                'will be used']).format(orig_comp_name)
            logging.warning(warn_msg)
            injX = 0.0
            injY = 0.0

        # Get injection times and rates specified in the component_data
        injTimes, injRates = self.check_injection_data(orig_comp_name, component_data)

        if isinstance(injTimes, list):
            injTimes = np.array(injTimes)

        # Convert the injTimes from years to days
        injTimes *= 365.25

        # Setup attributes of the component corresponding to the injection and locations
        # setup=1 since if one of injTimes, injRates is None the other still has
        # to be setup and in correct way according to the (possibly) new number
        # of injection wells
        self.process_injection_and_location_args(
            injX, injY, locX, locY, injTimes, injRates, setup=1)

        # Take care of parameter reservoirThickness
        strata = name2obj_dict['strata']
        par_name = 'reservoirThickness'
        if (par_name not in self.pars) and (par_name not in self.deterministic_pars):
            connect = None
            if par_name in strata.pars:
                connect = strata.pars
            elif par_name in strata.deterministic_pars:
                connect = strata.deterministic_pars
            elif par_name in strata.default_pars:
                connect = strata.default_pars
            else:
                err_msg = ''.join(['Unable to find "reservoirThickness" ',
                                   'parameter. Please check setup of the stratigraphy.'])
                logging.error(err_msg)
                raise KeyError(err_msg)

            self.add_par_linked_to_par(par_name, connect[par_name])

    @staticmethod
    def get_injection_data(comp_name, component_data, input_type):
        """
        Used for injRates and injTimes in the control file interface. This
        method checks if the input is a list or a string.

        If string inputs are provided for injRates and/or injTimes, this method
        checks if the string is a path to a .csv file containing injRates or injTimes.
        If so, the corresponding input is returned.
        """
        input_from_yaml = None

        # Check in what form the injection rates or injection times are provided
        if isinstance(component_data[input_type], list):
            input_from_yaml = component_data[input_type]

        elif isinstance(component_data[input_type], str):  # possibly a path to file
            fileDirName = os.path.join(IAM_DIR, component_data[input_type])

            if os.path.isfile(fileDirName):
                data = np.genfromtxt(
                    fileDirName, delimiter=",", dtype='f8', comments='#')

                if len(data.shape) == 1:  # data for a single injection well provided
                    # Check if the first row has a string label
                    if np.isnan(data[0]):
                        with open(fileDirName, newline='') as csvfile:
                            reader = csv.reader(csvfile)
                            label = next(reader)
                        csvfile.close()

                        if label[0] == input_type:  # input_type is 'injTimes' or 'injRates'
                            input_from_yaml = data[1:None].tolist()
                        else:
                            warning_msg = ''.join([
                                'The file provided for {inp_type} of the ',
                                'Theis Reservoir component {nm} has a non-numeric ',
                                'value in the first row. The value is not ',
                                'the label {inp_type}, however, which causes a lack ',
                                'of clarity. The {inp_type} will be set to the ',
                                'remaining entries (2nd row and on), which ',
                                'may result in errors. Check your input.']).format(
                                    inp_type=input_type, nm=comp_name)
                            logging.debug(warning_msg)
                            input_from_yaml = data[1:None].tolist()
                    else:
                        input_from_yaml = data.tolist()

                else:  # data for multiple injection wells is provided

                    # If the first row doesn't start with '#' and the first row
                    # has labels, that row will be read as NaN values.
                    # For cases when two columns are given for injRates and injTimes
                    if np.isnan(data[0, 0]) and data.shape[1] == 2:
                        with open(fileDirName, newline='') as csvfile:
                            reader = csv.reader(csvfile)
                            labels = next(reader)
                        csvfile.close()

                        selected_col = None
                        for column_ind in range(0, len(labels)):
                            if labels[column_ind] == input_type:
                                selected_col = column_ind
                                break

                        if selected_col is not None:
                            input_from_yaml = data[1:None, selected_col].tolist()

                        else:
                            warning_msg = ''.join([
                                'The file provided for {} of the ',
                                'TheisReservoir component {} has more than ',
                                'one column and the columns have string ',
                                'labels {}. When searching ',
                                'for the {} input, the label {} was not ',
                                'found among present in the header ',
                                'of the file {}. The {} will be set to the ',
                                'remaining entries (2nd row and on) of the ',
                                'first column, which may result in errors.',
                                'Check your input.']).format(
                                    input_type, comp_name, labels, input_type,
                                    input_type, fileDirName, input_type)

                            logging.debug(warning_msg)
                            input_from_yaml = data[1:None, 0].tolist()

                    elif not np.isnan(data[0, 0]):
                        # Get injection well locations
                        injX = component_data['InjectionWell']['coordx']
                        num_inj_wells = len(np.asarray(injX))
                        if data.shape[0] == num_inj_wells:
                            input_from_yaml = data
                        elif input_type == 'injRates':
                            err_msg = ''.join([
                                'Number of rows ({}) in the file {} provided ',
                                'for the setup of injRates does not equal to ',
                                'the number of injection wells ({}). Check your ',
                                'input.']).format(
                                    data.shape[0], fileDirName, num_inj_wells)
                            logging.error(err_msg)
                            raise TypeError(err_msg)

        return input_from_yaml


    def simulation_model(self, p, injX=None, injY=None, locX=None, locY=None,
                         injTimes=None, injRates=None):
        """
        Return pressure and CO2 saturation at the bottom of leaking well.
        Note that the names of parameters contained in the input dictionary
        coincide with those defined in this modules docstring.

        :param p: input parameters of simple reservoir model
        :type p: dict

        :param time_step: difference between the current and previous
            time points in days; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param injX: x-coordinate(s) of the injector/producer locations
        :type injX: list of floats

        :param injY: y-coordinate(s) of the injector/producer locations
        :type injY: list of floats

        :param locX: x-coordinate of the location for pressure and CO2 saturation
            to be calculated
        :type locX: float

        :param locY: y-coordinate of the location for pressure and CO2 saturation
            to be calculated
        :type locY: float

        :param injTimes: injection time(s) time points (in days) for which the pressure
            are to be calculated; by default, its value is 365.25
            (1 year in days)
        :type injTimes: list of floats

        :param injRates: injection rate(s) (+ injection, - withdrawal)
        :type injRates: list of floats

        :returns: out - dictionary of observations of Theis reservoir
            component model;
            keys: ['pressure']

        """
        debug_msg = ''.join([
            'TheisReservoir component {} model call input: ',
            'parameters {}, location (x, y) ({}, {})']).format(
                self.name, p, self.locX, self.locY)
        logging.debug(debug_msg)

        # Locations and injection setup
        # For setup=2 parameters injTimes, injRates are not resetup if the value
        # is None
        self.process_injection_and_location_args(injX, injY, locX, locY,
                                                 injTimes, injRates, setup=2)

        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}

        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Create instance of Parameters class
        inputParameters = tresrom.Parameters()

        # Set up parameters
        inputParameters.reservoirThickness = actual_p['reservoirThickness']
        inputParameters.reservoirPermeability = 10**actual_p['logResPerm']
        inputParameters.reservoirPorosity = actual_p['reservoirPorosity']
        inputParameters.initialPressure = actual_p['initialPressure']
        inputParameters.brineDensity = actual_p['brineDensity']
        inputParameters.CO2Density = actual_p['CO2Density']
        inputParameters.brineViscosity = actual_p['brineViscosity']
        inputParameters.compressibility = actual_p['compressibility']

        # Get injection rates at all times and corresponding times
        injTimes = self.injTimes
        injRates = self.injRates

        if self._parent.time_array is not None:
            time_array = self._parent.time_array

            injRates = np.zeros((self.num_inj_wells, len(time_array)))
            injTimes = np.tile(time_array, (self.num_inj_wells, 1))

            for ind in range(self.num_inj_wells):
                if (len(time_array) != len(self.injTimes[ind, :])) or (
                        np.not_equal(time_array, self.injTimes[ind, :]).any()):
                    injRates[ind, :] = np.interp(
                        time_array, self.injTimes[ind, :], self.injRates[ind, :])
                else:
                    injRates[ind, :] = self.injRates[ind, :]

        # Create solution object with defined input parameters
        sol = tresrom.Solution(inputParameters)
        dp = np.zeros(injTimes.shape)

        for ind in range(self.num_inj_wells):
            # Setup distance between injection and leaking well
            sol.update_distance(sqrt((self.locX-self.injX[ind])**2+(self.locY-self.injY[ind])**2))

            # Setup new injection rates
            sol.update_inj_rates(injRates[ind, :])

            # Setup new injection times
            sol.update_inj_times(injTimes[ind, :])

            # Find solution corresponding to the inputParameters
            sol.find()

            # Get pressure
            dp[ind, :] = sol.dp

        if self.num_inj_wells > 1:
            pressure = np.sum(dp, axis=0) + inputParameters.initialPressure
        else:
            pressure = dp[0, :] + inputParameters.initialPressure

        # Initialize dictionary of pressures
        out = dict()

        for i, t in enumerate(time_array):
            out['time_{}'.format(i)] = t

        for i, val in enumerate(pressure):
            out['pressure_{}'.format(i)] = val
            out['CO2saturation_{}'.format(i)] = 0.0

        debug_msg = 'TheisReservoir component {} model call output: {}'.format(
            self.name, out)
        logging.debug(debug_msg)

        # Return dictionary of outputs
        return out

    def reset(self):
        pass

    # Attributes for system connections
    needs_locXY = True
    needs_injXY = True

    system_params = ['reservoirThickness']


if __name__ == "__main__":
    __spec__ = None

    logging.basicConfig(level=logging.WARNING)
    # Define keyword arguments of the system model
    time_array = np.array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.])
    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    rates = [5.0e-3, 6.0e-3, 7.0e-3, 8.0e-3, 0, 0, 0, 0, 0, 0, 0, 0,]
    tres = sm.add_component_model_object(
        TheisReservoir(name='tres', parent=sm,
                       injX=100, injY=100, locX=100, locY=45,
                       injTimes=time_array, injRates=rates))

    tres.add_obs('pressure')
    tres.add_obs('CO2saturation')

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    # Set input parameters
    tres.add_par('initialPressure', value=1.0e6)      # Pa
    tres.add_par('reservoirThickness', value=30)      # m
    tres.add_par('logResPerm', value=-10.69897)       # m^2
    tres.add_par('reservoirPorosity', value=.2)
    tres.add_par('brineDensity', value=1000)          # kg/m^3
    tres.add_par('brineViscosity', value=2.535e-3)    # Pa*s
    tres.add_par('CO2Density', value=800)             # kg/m^3
    tres.add_par('compressibility', value=2.46e-9)    # 1/Pa

    # Run system model
    sm.forward()

    pressures = sm.collect_observations_as_time_series(tres, 'pressure')
    print(repr(pressures))
    saturation = sm.collect_observations_as_time_series(tres, 'CO2saturation')
    print(repr(saturation))

    def twolines(x, y1, y2, label1, label2):
        _, ax = plt.subplots()

        # Plot linear sequence, and set tick labels to the same color
        ax.plot(x, y1, color='red')
        ax.tick_params(axis='y', labelcolor='red')
        ax.set_ylabel(label1, color='red')

        # Generate a new Axes instance, on the twin-X axes (same position)
        ax2 = ax.twinx()
        ax2.plot(x, y2, color='blue', linestyle='dashed')
        ax2.tick_params(axis='y', labelcolor='blue')
        ax2.set_ylabel(label2, color='blue')

        plt.tight_layout()
        plt.show()

    twolines(time_array, rates, pressures, 'Injection Rate', 'Pressure')
