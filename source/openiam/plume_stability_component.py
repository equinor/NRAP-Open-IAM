# -*- coding: utf-8 -*-
import logging
from re import sub, split
import sys
import os

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import SystemModel, ComponentModel, IAM_DIR
    from openiam.mesh2D import read_Mesh2D_data
    from openiam.reservoir_data_interpolator import (check_file_format,
                                                     read_time_points,
                                                     read_data_headers)
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))


class PlumeStability(ComponentModel):
    """
    The Plume Stability component model produces quantitative metrics of the area,
    change in area over time, mobility and spreading :cite:`harp2019development`.
    Plume mobility is the effective centroid velocity including the speed and
    direction of movement. Plume spreading is the effective longitudinal dispersion
    of the plume along its direction of maximum elongation. This direction is
    returned by the model as well. The mobility and spreading metrics are
    comprehensive in that they can effectively handle and account for complex
    continuous and discontinuous plumes and intra‐plume migration. The metrics are
    calculated using 2D-scalar attribute field values as inputs. The model can read
    in field values formatted in the NRAP-Open-IAM dataset format. In order to
    process 3D scalar field data for the model, it is recommended to collapse the 3D
    data to 2D using a maximization approach as described in
    :cite:`harp2019development`.

    In the NRAP-Open-IAM control file, the type name for the Plume Stability component
    is ``PlumeStability``. The data files providing input for the Plume Stability component
    need to satisfy the same requirements imposed on the data files used as input
    for the Lookup Table Reservoir component. In particular, for control file
    setup of the Plume Stability component the following three keywords have
    the same meaning:

    * ``FileDirectory`` is a directory where files with the simulation data for the
      component are located, and which contains files described below;

    * ``TimeFile`` keyword specifies a name of the .csv file that stores the
      time points (in years) at which the results in the data files are provided;

    * ``ParametersFilename`` keyword contains a name of the .csv file containing
      the names and values of the parameters used to create the given set
      of data files; in addition, it lists the names of the .csv files in the folder
      ``FileDirectory`` containing simulation data for each of the data
      file in the set.

    The additional keywords of the component's control file interface are:

    * ``Variables`` is a list of observations names provided in the data files and
      for which some (or all) metrics will be calculated;

    * ``Thresholds`` is a dictionary of pairs (observation name, value) providing
      threshold value above which the change in the observation value should be
      taken into account for the calculation of the plume stability metrics.

    The only component model input parameter is ``index`` which indicates
    the index of the data file from the list in the last column of ``ParametersFilename``
    to be used to produce plume stability metrics. The minimum and maximum value
    of the parameter is defined by the indices of data files provided in the list.

    Possible observations from the Plume Stability component are

    * **{obs}_areas** - area of the plume above the predefined threshold

    * **{obs}_areas_dt** - change in the area of the plume above the predefined threshold

    * **{obs}_mobility** - velocity of centroid of the plume above the predefined threshold

    * **{obs}_mobility_angles** - angles/direction at which the centroid of the plume
      above the predefined threshold is changing

    * **{obs}_spreading** - longitudinal dispersion of the plume above
      the predefined threshold along its direction of maximum elongation

    * **{obs}_spreading_angles** - angles/direction at which the dispersion
      of the plume occurs.

    Above, {obs} determines the name of observation for which the plume stability
    metrics are to be calculated and for which the data is provided in the data files
    used as input for the component, e.g. **pressure**, **CO2saturation**, etc.
    """
    def __init__(self, name, parent, file_directory=None,
                 variable_names=None, thresholds=None,
                 parameter_filename='parameters_and_filenames.csv',
                 time_file='time_points.csv'):
        """
        Constructor method of PlumeStability class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param file_directory: Location of reservoir simulation data in NRAP-Open-IAM format
        :type file_directory: string

        :param variable_names: Names of lookup table variables to perform
            plume stability on (e.g., pressure)
        :type variable_names: list of strings

        :param thresholds: Dictionary of cutoff value for the variable value to define
            if the plume exists at a location; keys are names of the variables
        :type thresholds: {}

        :param parameter_filename: name of *.csv file containing information
            about sets of model parameters and files corresponding
            to the reservoir simulation data
        :type parameter_filename: str

        :param time_file: name of *.csv file to read information about time
            points (in years) at which observations (e.g., pressure and |CO2|
            saturation) in the simulation data were calculated
        :type filename: str

        :returns: PlumeStability class object
        """
        super().__init__(name, parent, model=self.simulation_model)

        # Add type attribute
        self.class_type = 'PlumeStability'

        # Setup instance attributes
        self.filenames = {}
        self.parameter_values = []
        self.file_directory = None
        self.parameter_filename = None

        self.time_file = None
        self.time_points = None

        # List of variables available in the provided lookup tables
        self.variable_list = []
        # List of variables for which the observations will be provided
        self.variable_names = []
        # Dictionary of thresholds for variables in the self.variable_names list
        self.thresholds = {}

        # Check whether file_directory is provided
        if file_directory is not None:
            # Save data folder
            self.file_directory = file_directory
            # Save parameter name
            self.parameter_filename = parameter_filename
            # Read file with info about data files and corresponding parameters
            self.read_parameter_filename()

            # Check and save name of the file with the time points to use in the model method
            self.hdf5_time_file, self.time_file = check_file_format(
                time_file, data_type='time')

            # Read time points
            _, self.time_points = read_time_points(self.hdf5_time_file,
                                                   self.file_directory,
                                                   self.time_file)

            # Read names of variables provided in the lookup tables
            self.read_possible_variables()

            # Set variable_names attribute
            self.set_variable_names(variable_names)

            # Set thresholds attribute
            self.set_thresholds(thresholds)

            # Define default value of index parameter
            # Default value indicates the first dataset in the provided list of filenames
            self.add_default_par('index', value=self.parameter_values[0][0])
        else :
            # Define default value of index parameter
            # Default value indicates the first dataset in the provided list of filenames
            self.add_default_par('index', value=1)

        # Indicate to the system model that this component will be run only once:
        # at the very first time point
        self.run_frequency = 1
        self.default_run_frequency = 1

        msg = 'PlumeStability object created with name {name}'.format(name=name)
        logging.debug(msg)

    def read_parameter_filename(self):
        """ Read file with names of data files and corresponding parameter values."""

        # By default parameter_filename is 'parameters_and_filenames.csv'
        # For the control file interface the name of the file can be different
        with open(os.path.join(self.file_directory, self.parameter_filename)) as fh:
            # Get names of the parameters from the first row (the first name one should be index)
            names = fh.readline().strip().split(',')
            if names[0] != 'index':
                err_msg = "".join([
                "The first column of the file {} should contain the indices ",
                "of the lookup tables specified in the last column (filename) and ",
                "be named 'index'. Please update the file by adding ",
                "the corresponding column and try again. Alternatively, ",
                "one can download the updated data set from its original ",
                "location."]).format(self.parameter_filename)
                logging.error(err_msg)
                raise NameError(err_msg)

            for l in fh:
                vs = l.strip().split(',')
                # Assumes parameter values are specified in the first columns except the last
                self.parameter_values.append([float(v) for v in vs[0:-1]])
                # Assumes filename will be the last column
                self.filenames[int(vs[0])] = vs[-1]

    def set_variable_names(self, variable_names):
        """ Set variable_names attribute content. """
        for var_name in variable_names:
            if var_name in self.variable_list:
                self.variable_names.append(var_name)
            else:
                err_msg = '{} is not a valid variable\n'.format(var_name)
                err_msg += 'Valid variables include: {}'.format(self.variable_list)
                logging.error(err_msg)

    def set_thresholds(self, thresholds):
        """ Set thresholds attribute. """
        for var_name in thresholds:
            self.thresholds[var_name] = thresholds[var_name]

        # Check whether threshold is setup for every variable for which
        # plume metrics have to be calculated
        for var_name in self.variable_names:
            if var_name not in self.thresholds:
                err_msg = 'Threshold is not setup for variable {}'.format(var_name)
                logging.error(err_msg)

    def read_possible_variables(self):
        """ Read names of variables provided in the lookup tables."""

        # Assume that all data files have identical variables
        # Assume that 'x', 'y', 'z', and 'area' are the only possible
        # non-variable header entries
        # Get file with the first index
        ind1 = list(self.filenames.keys())[0]

        # Check whether file is of appropriate type
        self.hdf5_data_format, _ = check_file_format(self.filenames[ind1],
                                                     data_type='data')

        # Read headers in the data file
        self.data_headers = read_data_headers(
            self.hdf5_data_format, self.file_directory, self.filenames[ind1])

        if not "area" in self.data_headers:
            wng_msg = ''.join([
                'Areas are not provided in the file directory {}, ',
                'uniform areas will be approximated based on ',
                'the assumption that the mesh is uniformly spaced.']).format(
                    self.file_directory)
            logging.warning(wng_msg)
        # Remove # symbol from first name if it exists
        self.data_headers[0] = sub("^#", '', self.data_headers[0]).strip()
        # Remove non-time varying variables assuming that
        # they will not end in '_[0-9]+' (underscore followed by an integer)
        tsnames = []
        for nm in self.data_headers:
            split_res = split("_[0-9]+$", nm)
            if len(split_res) == 2:
                tsnames.append(split_res[0])
        # Make sure that the number of variable names is divisible by the number of time_points
        if (len(tsnames))/(len(self.time_points))%1 == 0.:
            self.variable_list = np.unique(tsnames)
        else:
            err_msg = ''.join([
                'The number of variable names is not consistent with ',
                'the number of time points in {}.']).format(
                    os.path.join.join(self.file_directory, self.filenames[ind1]))
            logging.error(err_msg)

    def add_variable(self, variable, threshold):
        """ Add variable to the variables list if it is not there already."""

        if variable in self.variable_list:
            if variable not in self.variable_names:
                self.variable_names.append(variable)
                self.thresholds[variable] = threshold
            else:
                msg = ''.join([
                    '{} variable is already in the list of variables. ',
                    'Threshold is updated to be {}.']).format(variable, threshold)
                logging.debug(msg)
                # Update the threshold
                self.thresholds[variable] = threshold
        else:
            err_msg = '{} is not a valid variable. '.format(variable)
            err_msg += 'Valid variables include: {}.'.format(self.variable_list)
            logging.error(err_msg)

    @property
    def realization_weight(self):
        return self.pars['index'].discrete_vals[1]

    @realization_weight.setter
    def realization_weight(self, values):
        values = np.array(values)
        if values.shape == self.pars['index'].discrete_vals[1].shape:
            # Ensure probabilities sum to one
            values /= np.sum(values)
            self.pars['index'].discrete_vals[1] = values
            warn_msg = ''.join([
                "Realization weights have been ",
                "normalized to sum to one. Proportionally, they are ",
                "consistent with values entered.\n"])
            logging.warning(warn_msg)

    def connect_with_system(self, component_data, *args, **kwargs):
        """
        Code to add PlumeStability component to system model for control file interface.

        :param component_data: Dictionary of component data
        :type component_data: dict

        :returns: None
        """
        if 'FileDirectory' in component_data:
            self.file_directory = os.path.join(IAM_DIR,
                                               component_data['FileDirectory'])
        else:
            err_msg = ''.join(['FileDirectory keyword must be provided for ',
                               'PlumeStability Component.'])
            logging.error(err_msg)
            raise IOError(err_msg)

        if 'ParameterFilename' in component_data:
            self.parameter_filename = component_data['ParameterFilename']
        else:
            self.parameter_filename = 'parameters_and_filenames.csv'

        # Read file with info about data files and corresponding parameters
        self.read_parameter_filename()

        if 'TimeFile' in component_data:
            self.time_file = component_data['TimeFile']
            self.hdf5_time_file, self.time_file = check_file_format(
                self.time_file, data_type='time')
        else:
            self.time_file = 'time_points.csv'
            self.hdf5_time_file = False

        # Read in time points
        _, self.time_points = read_time_points(self.hdf5_time_file,
                                               self.file_directory,
                                               self.time_file)

        # Read names of variables provided in the lookup tables
        self.read_possible_variables()

        if 'Variables' in component_data:
            variable_names = component_data['Variables']
            # Set variable_names attribute
            self.set_variable_names(variable_names)
        else:
            err_msg = ''.join(['Variables keyword must be provided for ',
                               'PlumeStability Component.'])
            logging.error(err_msg)
            raise IOError(err_msg)

        if 'Thresholds' in component_data:
            thresholds = component_data['Thresholds']
            # Set thresholds attribute
            self.set_thresholds(thresholds)
        else:
            err_msg = ''.join(['Thresholds keyword must be provided for ',
                               'PlumeStability Component.'])
            logging.error(err_msg)
            raise IOError(err_msg)

        if 'Parameters' in component_data:
            par_name = 'index'
            if par_name in component_data['Parameters']:

                if not isinstance(component_data['Parameters'][par_name], dict):
                    component_data['Parameters'][par_name] = {
                        'value': component_data['Parameters'][par_name],
                        'vary': False}
                    self.add_par(par_name, **component_data['Parameters'][par_name])

                elif 'value' in component_data['Parameters'][par_name]:
                    value = component_data['Parameters'][par_name]['value']
                    self.add_par(par_name, value=value, vary=False)
                    debug_msg = ''.join([
                        'Parameter {} of PlumeStability component ',
                        'created with value {}.']).format(par_name, value)
                    logging.debug(debug_msg)

                else:
                    if 'values' in component_data['Parameters'][par_name]:
                        component_data['Parameters'][par_name]['Values'] = (
                            component_data['Parameters'][par_name]['values'])
                    values = component_data['Parameters'][par_name]['Values']
                    # Number of provided values
                    nv = len(values)

                    if 'Weights' not in component_data['Parameters'][par_name]:
                        if 'weights' in component_data['Parameters'][par_name]:
                            weights = component_data['Parameters'][par_name]['weights']
                        else:
                            weights = [1.0/nv for v in values]
                    else:
                        weights = component_data['Parameters'][par_name]['Weights']

                    self.add_par(par_name, value=values[0],
                                 vary=True, discrete_vals=(values, weights))
                    debug_msg = ''.join([
                        'Parameter {name} of PlumeStability component created with ',
                        'values {values} and weights {weights}.']).format(
                            name=par_name, values=values, weights=weights)
                    logging.debug(debug_msg)

    def simulation_model(self, p):
        """
        Calculate plume stability metrics based on moment analysis method
        described in Harp et al. (GG:S&T, 2019).

        :param p: input parameters
        :type p: dict

        :returns: out - dictionary of plume stability metrics and times

        """
        msg = '{sr} model call input {p}'.format(sr=self.name, p=p)
        logging.debug(msg)

        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}

        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # index is limited to be within the set of indices of linked data sets
        index = int(actual_p['index'])

        # verify file type for data file to be used and
        # change if necessary
        self.hdf5_data_format, _ = check_file_format(self.filenames[index],
                                                     data_type='data')

        # Create Mesh2D object
        M = read_Mesh2D_data(
            os.path.join(self.file_directory, self.filenames[index]),
            self.variable_list, self.time_points, self.hdf5_data_format)

        # Check that specified variables are valid
        for var_name in self.variable_names:
            if not var_name in self.variable_list:
                err_msg = '{} is not a valid variable. '.format(var_name)
                err_msg += 'Valid variables include: {}.'.format(self.variable_list)
                logging.error(err_msg)
            if var_name not in self.thresholds:
                err_msg = 'Threshold is not setup for variable {}.'.format(var_name)
                logging.error(err_msg)

        # By default, interpolation is not required
        to_interpolate = False
        # By default, time points are defined by data set time poins
        time_array = self.time_points

        # Check that parent system model has time array defined
        if self._parent.time_array is not None:
            # Time points in time_array of the system model are defined in days
            # so we need to convert them to years to compare
            t_array = self._parent.time_array/365.25
            # Check whether data set points coincide with system model time points
            if (len(t_array) != len(self.time_points)) or (
                    np.not_equal(t_array, self.time_points).any()):
                to_interpolate = True
                time_array = t_array
                warn_msg = ''.join(['Time points of the dataset do not coincide ',
                                    'with the points of the system ',
                                    'time array. Temporal interpolation ',
                                    'will be performed for the system ',
                                    'time array.'])
                logging.warning(warn_msg)
        else:
            err_msg = 'Argument time_array has to be defined for system model.'
            logging.error(err_msg)

        out = {}
        # Cycle over all variables for which the metrics were requested
        for var_name in self.variable_names:
            # Get the specified variable
            V = M.variables[var_name]
            # Get the corresponding threshold
            threshold = self.thresholds[var_name]

            # Collect some outputs
            times = V.times

            # Calculate area and its derivative
            areas = V.plume_areas(threshold)
            areas_dt = V.plume_areas_dt(threshold)
            # Calculate mobility and angles
            mobility, mobility_angles = V.mobility(threshold)
            mobility_angles[np.where(np.isnan(mobility_angles))[0]] = -9999
            # Calculate spreading and angles
            spreading, spreading_angles = V.spreading(threshold)
            spreading_angles[np.where(np.isnan(spreading_angles))[0]] = -9999

            if to_interpolate:
                areas = np.interp(time_array, times, areas)
                areas_dt = np.interp(time_array, times, areas_dt)
                mobility = np.interp(time_array, times, mobility)
                mobility_angles = np.interp(time_array, times, mobility_angles)
                spreading = np.interp(time_array, times, spreading)
                spreading_angles = np.interp(time_array, times, spreading_angles)

            for i, t in enumerate(areas):
                out['{}_areas_{}'.format(var_name, i)] = t
            for i, t in enumerate(areas_dt):
                out['{}_areas_dt_{}'.format(var_name, i)] = t
            for i, t in enumerate(mobility):
                out['{}_mobility_{}'.format(var_name, i)] = t
            for i, t in enumerate(mobility_angles):
                out['{}_mobility_angles_{}'.format(var_name, i)] = t
            for i, t in enumerate(spreading):
                out['{}_spreading_{}'.format(var_name, i)] = t
            for i, t in enumerate(spreading_angles):
                out['{}_spreading_angles_{}'.format(var_name, i)] = t

        for i, t in enumerate(time_array):
            out['times_{}'.format(i)] = t

        msg = '{sr} model call output {out}'.format(sr=self.name, out=out)
        logging.debug(msg)

        # Return dictionary of outputs
        return out

    def reset(self):
        pass

if __name__ == "__main__":
    # For multiprocessing in Spyder
    __spec__ = None
    logging.basicConfig(level=logging.WARNING)
    file_directory = os.sep.join(['..', 'components', 'reservoir',
                                  'lookuptables', 'Kimb_54_sims'])

    if not os.path.exists(os.sep.join([file_directory, 'Reservoir_data_sim01.csv'])):
        msg = ''.join(['\nKimberlina data set can be downloaded ',
                       'from one of the following places:\n',
                       '1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam \n',
                       '2. https://gitlab.com/NRAP/Kimberlina_data \n'])
        logging.error(msg)

    option = 2
    if option == 1:
        num_years = 200
        # Time array different from time points provided in the data set
        time_array = 365.25*np.arange(0.0, num_years+1)
        data_dim = "2D"
        parameter_filename = 'parameters_and_filenames.csv'
        nvals = 10
    elif option == 2:
        # Time array is the same as in the data set but converted to days
        time_array = np.genfromtxt(
            os.sep.join([file_directory, 'time_points.csv']), delimiter=',')*365.25
        data_dim = "2D"
        parameter_filename = 'parameters_and_filenames.csv'
        nvals = 10
    elif option == 3:
        # Time array is the same as in the data set but converted to days
        time_array = np.genfromtxt(
            os.sep.join([file_directory, 'time_points.csv']), delimiter=',')*365.25
        data_dim = "3D"
        parameter_filename = 'parameters_and_filenames_h5.csv'
        nvals = 2

    # Time array argument has to be defined for use of plume stability component.
    # It has to be defined as in the data set or whatever (in the latter case
    # interpolation will be used).
    sm_model_kwargs = {'time_array': time_array}
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    sps = sm.add_component_model_object(
        PlumeStability(name='sps', parent=sm, file_directory=file_directory,
                       variable_names=['pressure'],  # can be a list of names
                       thresholds={'pressure': 1.0e6},  # must be a dictionary
                       # containing thresholds for all variables in variable_names
                       parameter_filename=parameter_filename,
                       time_file='time_points.csv'))
    # Variable(s) can be modified
    # Print available variables
    print('Variables available in data set', sps.variable_list)
    # Print names of variables for which metrics has to be calculated
    print('Metrics to be calculated for', sps.variable_names)

    # Reset variable and threshold
    sps.add_variable('CO2saturation', 0.01)
    sps.variable_names.remove('pressure')
    print('Variables are updated')

    # Add discrete parameter datafile_index. For 2D processing, each file with data has the same probability
    # of being used for sampling the observation.
    # When parameter is added as discrete, it has an associated with it value which is used
    # for forward simulations. For range of values from 1 to 54 the value of 28
    # is used as the value located at the next position from the center of the list
    # or right at the center depending whether the number of values is even or odd.
    # For example, for list of values [3, 7, 1, 2, 9] value of 1 would be used as value
    # for forward simulation.
    # For 3D example, either index 55 (CSV file) or
    # 56 (HDF5 file) will be used with equal probability.
    if data_dim == '2D':
        sps.add_par('index', discrete_vals=(
        list(range(1, len(sps.filenames)+1)),
        [1./len(sps.filenames)]*len(sps.filenames)))

    elif data_dim == '3D':
        sps.add_par('index', discrete_vals=([1, 2], [0.5, 0.5]))


    # Observations have to be added explicitly
    sps.add_obs('times')
    for var_name in sps.variable_names:
        sps.add_obs('{}_areas'.format(var_name))
        sps.add_obs('{}_areas_dt'.format(var_name))
        sps.add_obs('{}_mobility'.format(var_name))
        sps.add_obs('{}_spreading'.format(var_name))

    # Print names of variables for which metrics has to be calculated
    print('Metrics to be calculated for', sps.variable_names)
    # One can also do the same as above by
    # sps.variable_names = ['CO2saturation']
    # sps.thresholds['CO2saturation'] = 0.01

    # 'index' parameter is within the set of indices for lookup tables.
    # To change the dataset probability, one can do, for example
    # Collect current weights
    wts = sps.realization_weight
    # Modify weight of the 1st dataset, for example
    wts[0] *= 2
    # Assign modified weights back to 'index' parameter
    # Note that this will produce a warning that the weights have been renormalized to sum to one.
    sps.realization_weight = wts

    out = sm.forward()
    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    for var_name in sps.variable_names:
        print('{} area'.format(var_name),
              sm.collect_observations_as_time_series(
                  sps, '{}_areas'.format(var_name)),
              sep='\n')
        print('------------------------------------------------------------------')
        print('{} area deriv'.format(var_name),
              sm.collect_observations_as_time_series(
                  sps, '{}_areas_dt'.format(var_name)),
              sep='\n')
        print('------------------------------------------------------------------')
        print('{} mobility'.format(var_name),
              sm.collect_observations_as_time_series(
                  sps, '{}_mobility'.format(var_name)),
              sep='\n')

    # One can run system model with a different value of index using
    # out = sm.forward(pardict={'sps.index': 1})

    # Create sampleset with 10 evenly spaced realization ids
    s = sm.parstudy(nvals=nvals)
    # One can create a sampleset of all realizations in Kimb_54_sims
    #    s = sm.create_sampleset([[v] for v in range(1, 55)])
    # Run sampleset
    s.run(cpus=nvals//2+1, verbose=False)
    # Postprocess output into easy to use dictionary
    out = s.collect_observations_as_time_series()

    # Plot results
    font_size = 14
    f1, ax = plt.subplots(4, sharex=True, figsize=(10, 10))
    for i in range(len(s.indices)):
        # Convert output from m^2 to km^2 or from m to km
        ax[0].plot(out['sps.times'][i], out['sps.CO2saturation_areas'][i]/1.0e+6, '-', linewidth=3)
        ax[1].plot(out['sps.times'][i], out['sps.CO2saturation_areas_dt'][i]/1.0e+6, '-', linewidth=3)
        ax[2].plot(out['sps.times'][i], out['sps.CO2saturation_mobility'][i]/1.0e+3, '-', linewidth=3)
        ax[3].plot(out['sps.times'][i], out['sps.CO2saturation_spreading'][i]/1.0e+6, '-', linewidth=3)

    f1.subplots_adjust(left=0.2, bottom=0.1, right=0.95, top=0.95, wspace=0.1)
    ax[0].set_ylabel(r'Plume area,{}[km$^2$]'.format('\n'), fontsize=font_size)
    ax[1].set_ylabel(r'Change in plume area,{}[km$^2$/year]'.format('\n'), fontsize=font_size)
    ax[2].set_ylabel(r'Mobility,{}[m/year]'.format('\n'), fontsize=font_size)
    ax[3].set_ylabel(r'Spreading,{}[km$^2$/year]'.format('\n'), fontsize=font_size)
    ax[3].set_xlabel('Time [years]', fontsize=font_size)

    for ind in range(4):
        ax[ind].tick_params(axis='both', which='major', labelsize=font_size-2)
        ax[ind].set_xlim(0, 150)

    # Align y-labels
    labelx = -0.1
    for ind in range(4):
        ax[ind].yaxis.set_label_coords(labelx, 0.5)
    # The mobility and spreading angles could be plotted as well,
    # but it may not be that useful to plot an entire ensemble of them.

    if data_dim == "2D":
        f1.savefig(file_directory + '\\plume_stability_graph_2Da.jpg', format='jpg')
    elif data_dim == "3D":
        f1.savefig(file_directory + '\\plume_stability_graph_3Da.jpg', format='jpg')
    f1.show()

    # Plot one of the realizations
    font_size = 14
    f2, ax = plt.subplots(2, 2, figsize=(10, 10))
    for i in [0]:
        # Convert output from m^2 to km^2 or from m to km
        ax[0, 0].plot(out['sps.times'][i], out['sps.CO2saturation_areas'][i]/1.0e+6, '-k', linewidth=3)
        ax[0, 1].plot(out['sps.times'][i], out['sps.CO2saturation_areas_dt'][i]/1.0e+6, '-k', linewidth=3)
        ax[1, 0].plot(out['sps.times'][i], out['sps.CO2saturation_mobility'][i]/1.0e+3, '-k', linewidth=3)
        ax[1, 1].plot(out['sps.times'][i], out['sps.CO2saturation_spreading'][i]/1.0e+6, '-k', linewidth=3)

    f2.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95, wspace=0.3)
    ax[0, 0].set_ylabel(r'Plume area,{}[km$^2$]'.format('\n'), fontsize=font_size)
    ax[0, 1].set_ylabel(r'Change in plume area,{}[km$^2$/year]'.format('\n'), fontsize=font_size)
    ax[1, 0].set_ylabel(r'Mobility,{}[m/year]'.format('\n'), fontsize=font_size)
    ax[1, 1].set_ylabel(r'Spreading,{}[km$^2$/year]'.format('\n'), fontsize=font_size)

    for ind1 in range(2):
        for ind2 in range(2):
            ax[ind1, ind2].set_xlabel('Time [years]', fontsize=font_size)
            ax[ind1, ind2].tick_params(axis='both', which='major', labelsize=font_size-2)
            ax[ind1, ind2].set_xlim(0, 150)

    # Align y-labels
    labelx = -0.12
    for ind1 in range(2):
        for ind2 in range(2):
            ax[ind1, ind2].yaxis.set_label_coords(labelx, 0.5)

    if data_dim == "2D":
        f2.savefig(file_directory + '\\plume_stability_graph_2Db.jpg', format='jpg')
    elif data_dim == "3D":
        f2.savefig(file_directory + '\\plume_stability_graph_3Db.jpg', format='jpg')
    f2.show()

    # Remove all handlers from the logger for proper work in the consecutive runs
    while logging.getLogger('').handlers:
        logging.getLogger('').handlers.pop()
