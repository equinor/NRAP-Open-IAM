# -*- coding: utf-8 -*-
import os
import sys
import logging
import numpy as np
from numpy.random import default_rng

source_folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(source_folder)
sys.path.append(os.path.join(source_folder, 'components', 'seal'))


try:
    from openiam import SystemModel, SamplerModel
except ImportError as err:
    print('Unable to load IAM class module: '+ str(err))

try:
    import components.seal.seal_intro as intro
    import components.seal.seal_perm as perm
    import components.seal.frac_random as fran
    import components.seal.seal_units as sunit
    from components.seal.seal_setup import SEAL_SETUP_DICT
except ImportError:
    print('\nERROR: Unable to load Permeability Sampler for Seal Horizon\n')
    sys.exit()


SHPS_PARAMETERS = ['avePermeability', 'stdDevPermeability',
                   'minPermeability', 'maxPermeability', 'heterFactor']
SCPS_PARAMETERS = ['perm_mean', 'perm_std',
                   'perm_min', 'perm_max', 'perm_heter_factor']  # seal_controls pars


def create_perm_array(num_cells, mu_d, perm_scale,
                      min_val, max_val, seed=None):
    """
    Create permeability array.
    """
    if seed is None:
        rng = default_rng()
    else:
        rng = default_rng(seed)


    perm_array = np.zeros(num_cells)
    for ind in range(num_cells):
        # Change seed for each new cell: seedx = seed+ind
        perm_array[ind] = fran.evaluate_lognorm(
            mu_d, perm_scale, min_val, max_val, rng)

    return perm_array

def evaluate_areal_heter(perm_array, perm_heter_factor, percent=perm.PERCENT,
                         seed=None):
    """
    Update portion of a provided permeability array by a factor.
    """
    rng = default_rng(seed=seed)

    num_elems = len(perm_array)
    heter_limit = int(percent*num_elems/100.0) # how many elements need to be adjusted
    index_set = np.arange(num_elems)
    rng.shuffle(index_set) # shuffle indices

    for ind in range(heter_limit):
        heter_index = index_set[ind]
        perm_array[heter_index] = perm_array[heter_index]*perm_heter_factor

    return perm_array


class SHPermeabilitySampler(SamplerModel):
    """
    Class for model providing permeability parameter for Seal Horizon component.

    The class can be used as provider of permeability parameter values
    for any other component.

    Distribution parameters for permeability of the seal horizon layer

    * **avePermeability** [|m^2|] (1.0e-22 to 1.0e-16) - mean total vertical
      permeability of a lognormal distribution; equivalent value
      for fractured rock (default: 2.5e-16)

    * **stdDevPermeability** [|m^2|] (0 to 1.0e-17) - standard deviation
      of the total vertical permeability distribution (default: 5.0e-17)

    * **minPermeability** [|m^2|] (1.0e-24 to 1.0e-17) - minimum total vertical
      permeability; this value truncates (censors) the vertical random
      distribution and limits lower values (default: 1.0e-18)

    * **maxPermeability** [|m^2|] (1.0e-21 to 1.0e-12) - maximum total vertical
      permeability; this value truncates (censors) the random distribution and
      limits higher values (default: 1.0e-15)

    * **heterFactor** [-] (1.0e-2 to 100) - increase factor of the permeability
      of cells selected for heterogeneity, if the heterogeneity approach is used
      (default: 0.5).

    The possible output from the SH Permeability Sampler component is

    * **permeability** [|m^2|] - effective permeability of a layer or cells
    over a given area.
    """
    def __init__(self, name, parent, grid_shape,
                 constr_type='array', coordinates=None,
                 reproducible=True, heterogeneity_approach=False):
        """
        Constructor method of SHPermeabilitySampler class

        :param name: name of component model
        :type name: [str]

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel [object]

        :param grid_shape: shape of grid for which the output is produced;
            for constr_type='matrix' shape of output is grid_shape, i.e.,
            (grid_shape[0], grid_shape[1]);
            for constr_type='array' shape of output is (grid_shape[0]*grid_shape[1],)
        :type shape: tuple of one or two integers

        :param constr_type: type of the gridded observation. Possible values:
            'array', 'matrix'. The parameter indicates what type of gridded
            observation, a given component will provide. 'array' option means that
            the returned observation obs will be an array-like object and
            every element can be accessed as obs[ind]; 'matrix' means that the returned
            observation is a matrix-like object and every element can be accessed
            by providing 2 (or 3) indices of elements: either like
            obs[ind1][ind2] or obs[ind1,ind2] for 2d case, or
            obs[ind1][ind2][ind3] or obs[ind1,ind2,ind3] for 3d case.
        :type constr_type: str

        :param coordinates: x, y (and z) coordinates of the grid associated
            with the observation; by default, it is None. In general,
            coordinates should be a dictionary with keys 'x', 'y' (and 'z').
            coordinates[key] should be an array-like object
            for any key in coordinates.
        :type coordinates: dict of array-like objects

        :returns: SHPermeabilitySampler class object
        """
        super().__init__(name=name, parent=parent, reproducible=reproducible)

        # Add type attribute
        self.class_type = 'SHPermeabilitySampler'

        # Define gridded observations names
        self.grid_obs_keys = ['permeability']

        self.controls = {}
        # Add default parameters specific to this sampler
        # Note that seed parameter is already parameter of the sampler
        mD2m = sunit.microd_to_metersq()
        for ind, par_name in enumerate(SHPS_PARAMETERS[0:4]):
            self.add_default_par(
                par_name, value=mD2m*SEAL_SETUP_DICT['Permeability'][par_name]['valu'])
            self.controls[SCPS_PARAMETERS[ind]] = mD2m*SEAL_SETUP_DICT['Permeability'][
                par_name]['valu']
        self.controls['heterogeneityApproach'] = heterogeneity_approach
        self.add_default_par('heterFactor', value=SEAL_SETUP_DICT[
            'Permeability']['heterFactor'])
        self.controls['perm_heter_factor'] = SEAL_SETUP_DICT[
            'Permeability']['heterFactor']

        # Define parameters boundaries dictionary self.pars_bounds
        self.define_pars_bounds()

        # Logging output of default parameters and their values
        debug_msg = ''.join([par_nm+': '+str(self.default_pars[par_nm].value) \
            for par_nm in self.default_pars.keys()])
        logging.debug(debug_msg)

        # Define additional controls
        if len(grid_shape) == 1:
            self.controls['grid_approach'] = False
            self.controls['grid_rows'] = grid_shape[0]
            self.controls['grid_cols'] = 1
            self.controls['num_cells'] = grid_shape[0]
            if constr_type != 'array':
                warn_msg = ''.join([
                    'Argument constr_type has a wrong value of {}. ',
                    'Default value of "array" will be used ',
                    'instead.']).format(constr_type)
                logging.warning(warn_msg)
            self.controls['constr_type'] = 'array'
        elif len(grid_shape) == 2: # length is 2
            self.controls['grid_approach'] = True
            self.controls['grid_rows'] = grid_shape[0]
            self.controls['grid_cols'] = grid_shape[1]
            self.controls['num_cells'] = grid_shape[0]*grid_shape[1]
            self.controls['constr_type'] = constr_type
        else:
            err_msg = ''.join([
                'Argument grid_shape is of wrong length. It should have ',
                'length 1 or 2 but has length {}.']).format(len(grid_shape))
            logging.error(err_msg)
            raise ValueError(err_msg)

        # Add gridded observation to be linked
        # It makes sense to add it here since the model has only one observation
        # and it's assumed that if the object of SHPermeabilitySampler class is
        # created then it's observation will be used by some other component.
        # Otherwise, beyond test purposes there are no other reasons to create
        # the object
        self.add_obs_to_be_linked(
            'permeability', obs_type='grid', constr_type=self.controls['constr_type'],
            coordinates=coordinates, index=[0])

    def define_pars_bounds(self):
        """
        Define dictionary of parameter boundaries.

        The method uses boundaries imposed in the standalone Seal Flux code.
        """
        #
        param_bounds = intro.define_input_limits()

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        scalar = sunit.microd_to_metersq()
        self.pars_bounds['avePermeability'] = [
            scalar*param_bounds['perm_mean'][0], scalar*param_bounds['perm_mean'][1]]
        self.pars_bounds['stdDevPermeability'] = [
            scalar*param_bounds['perm_std'][0], scalar*param_bounds['perm_std'][1]]
        self.pars_bounds['minPermeability'] = [
            scalar*param_bounds['perm_min'][0], scalar*param_bounds['perm_min'][1]]
        self.pars_bounds['maxPermeability'] = [
            scalar*param_bounds['perm_max'][0], scalar*param_bounds['perm_max'][1]]
        self.pars_bounds['heterFactor'] = param_bounds['perm_heter_factor']

        # np.inf is used to ephasize that no upper limit is imposed
        self.pars_bounds['seed'] = [0, np.inf]

    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        debug_msg = 'Input parameters of component {} are {}'.format(self.name, p)
        logging.debug(debug_msg)

        # First check seed parameter if it's provided
        if 'seed' in p:
            # check_seed_parameter is a default method of SamplerModel Class
            # checking whether seed is a positive integer
            self.check_seed_parameter(p['seed'])

        # Check the rest of the sampler parameters
        for key in SHPS_PARAMETERS:
            if key in p:
                val = p[key]
                if ((val < self.pars_bounds[key][0]) or
                        (val > self.pars_bounds[key][1])):
                    warn_msg = ''.join([
                        'Parameter {} of value {} of SH Permeability Sampler {} ',
                        'is out of boundaries [{}, {}].']).format(
                            key, val, self.name, self.pars_bounds[key][0],
                            self.pars_bounds[key][1])
                    logging.warning(warn_msg)

    def sample(self, p):
        """
        :param p: input parameters of SHPermeabilitySampler model
        :type p: dict

        """
        # Obtain the default values of the parameters from dictionary of
        # default parameters.
        actual_p = {k: v.value for k, v in self.default_pars.items()}

        if not self.controls['heterogeneityApproach'] and 'heterFactor' in p:
            warn_msg = ''.join([
                'Parameter heterFactor of component {} will not be used ',
                'since argument heterogeneity_approach provided to the ',
                'constructor method has value of False.']).format(self.name)
            logging.warning(warn_msg)
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Update controls with parameter values
        m2mD = sunit.metersq_to_microd()
        for ind, par_name in enumerate(SCPS_PARAMETERS[0:4]):
            self.controls[par_name] = m2mD*actual_p[SHPS_PARAMETERS[ind]]

        # Get new parameters
        self.controls['perm_location'], self.controls['perm_scale'] = \
            intro.convert_lognorm_terms(self.controls['perm_mean'],
                                        self.controls['perm_std'])
        # Define output
        mD2m = sunit.microd_to_metersq()
        if self.reproducible:
            permeability = mD2m*create_perm_array(
                self.controls['num_cells'], self.controls['perm_location'],
                self.controls['perm_scale'], self.controls['perm_min'],
                self.controls['perm_max'], seed=actual_p['seed'])
        else:
            permeability = mD2m*create_perm_array(
                self.controls['num_cells'], self.controls['perm_location'],
                self.controls['perm_scale'], self.controls['perm_min'],
                self.controls['perm_max'])

        if self.controls['heterogeneityApproach']:
            if self.reproducible:
                permeability = evaluate_areal_heter(
                    permeability, actual_p['heterFactor'], seed=actual_p['seed']+1234)
            else:
                permeability = evaluate_areal_heter(permeability, actual_p['heterFactor'])

        # Check whether matrix-type output is needed
        if self.controls['constr_type'] == 'matrix':
            permeability = permeability.reshape(self.controls['grid_rows'],
                                                self.controls['grid_cols'])

        out = {'permeability': permeability}

        # Return output
        return out


def test_scenario1():
    """ Run script example illustrating work of SH Permeability Sampler model.

    """
    # Setup logging and constants.
    logging.basicConfig(level=logging.DEBUG)

    # Define keyword arguments of the system model.
    model_kwargs = {'time_point': 0.0} # time is given in days

    output_dir = '../../../output/samplers/'
    # Create system model.
    sm = SystemModel(model_kwargs=model_kwargs)
    perm_sampler = sm.add_component_model_object(
        SHPermeabilitySampler(name='ps', parent=sm, grid_shape=(10, 10)))
    perm_sampler.add_par('seed', value=234, vary=False)
    perm_sampler.add_par('avePermeability', value=1.0e-18, vary=False)
    perm_sampler.add_par('stdDevPermeability', value=1.0e-19, vary=False)
    perm_sampler.add_par('minPermeability', value=1.0e-19, vary=False)
    perm_sampler.add_par('maxPermeability', value=1.0e-17, vary=False)
    perm_sampler.add_grid_obs('permeability', constr_type='array', index=[0],
                              output_dir=output_dir)

    sm.forward()

    # Collect saved observation
    data = sm.collect_gridded_observations_as_time_series(
        perm_sampler, 'permeability', output_dir,
        indices=[0], rlzn_number=0)

    print(data[0].shape)
    print(data[0])


def test_scenario2():
    """ Run script example illustrating work of SH Permeability Sampler model.

    """
    # Setup logging and constants.
    logging.basicConfig(level=logging.DEBUG)

    # Define keyword arguments of the system model.
    model_kwargs = {'time_point': 0.0} # time is given in days

    output_dir = '../../../output/samplers/'
    # Create system model.
    sm = SystemModel(model_kwargs=model_kwargs)
    perm_sampler = sm.add_component_model_object(
        SHPermeabilitySampler(name='ps', parent=sm, grid_shape=(10, 10),
                              constr_type='matrix'))
    perm_sampler.add_par('seed', value=234, vary=False)
    perm_sampler.add_par('avePermeability', value=1.0e-18, vary=False)
    perm_sampler.add_par('stdDevPermeability', value=1.0e-19, vary=False)
    perm_sampler.add_par('minPermeability', value=1.0e-19, vary=False)
    perm_sampler.add_par('maxPermeability', value=1.0e-17, vary=False)
    perm_sampler.add_grid_obs('permeability', constr_type='matrix', index=[0],
                              output_dir=output_dir)

    sm.forward()

    # Collect saved observation
    data = sm.collect_gridded_observations_as_time_series(
        perm_sampler, 'permeability', output_dir,
        indices=[0], rlzn_number=0)

    print(data[0].shape)
    print(data[0])

    from matplotlib import cm
    import matplotlib.pyplot as plt
    fig1 = plt.figure(figsize=(8, 8))
    ax1 = fig1.add_subplot(111)
    im1 = ax1.imshow(np.flipud(data[0]))
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    plt.grid(False)
    fig1.colorbar(im1, ax=ax1)

    from mpl_toolkits.mplot3d import axes3d

    x = np.arange(1, 11)
    y = x
    xx, yy = np.meshgrid(x, y)
    xx = xx.reshape(100,)
    yy = yy.reshape(100,)
    zz = np.zeros(100,)

    dx = np.ones(100,)
    dy = np.ones(100,)
    dz = data[0].reshape(100,)

    cmap = cm.get_cmap('Spectral')
    max_height = np.max(dz)
    min_height = np.min(dz)
    colors = [cmap((k-min_height)/max_height) for k in dz]

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.bar3d(xx, yy, zz, dx, dy, dz, shade=True, color=colors)
    ax2.set_xlabel('x [m]')
    ax2.set_ylabel('y [m]')
    ax2.set_zlabel(r'k [$m^2$]')


def test_scenario3():
    """ Run script example illustrating work of SH Permeability Sampler model.

    """
    # Setup logging and constants.
    logging.basicConfig(level=logging.DEBUG)

    # Define keyword arguments of the system model.
    model_kwargs = {'time_point': 0.0} # time is given in days

    output_dir = '../../../output/samplers/'
    # Create system model.
    sm = SystemModel(model_kwargs=model_kwargs)
    perm_sampler = sm.add_component_model_object(
        SHPermeabilitySampler(name='ps', parent=sm, grid_shape=(100,)))
    perm_sampler.add_par('seed', value=234, vary=False)
    perm_sampler.add_par('avePermeability', value=1.0e-18, vary=False)
    perm_sampler.add_par('stdDevPermeability', value=1.0e-19, vary=False)
    perm_sampler.add_par('minPermeability', value=1.0e-19, vary=False)
    perm_sampler.add_par('maxPermeability', value=1.0e-17, vary=False)
    perm_sampler.add_grid_obs('permeability', constr_type='array', index=[0],
                              output_dir=output_dir)

    sm.forward()

    # Collect saved observation
    data = sm.collect_gridded_observations_as_time_series(
        perm_sampler, 'permeability', output_dir,
        indices=[0], rlzn_number=0)

    print(data[0].shape)
    print(data[0])

    import matplotlib.pyplot as plt
    x = np.arange(1, 101)
    plt.figure(figsize=(5, 5))
    plt.plot(x, data[0], '-o')
    plt.title('Layer permeability')
    plt.xlabel('x [m]')
    plt.ylabel('k [m]')


def test_scenario4():
    """ Run script example illustrating work of SH Permeability Sampler model.

    """
    # Setup logging and constants.
    logging.basicConfig(level=logging.DEBUG)

    # Define keyword arguments of the system model.
    model_kwargs = {'time_point': 0.0} # time is given in days

    output_dir = '../../../output/samplers/'
    # Create system model.
    sm = SystemModel(model_kwargs=model_kwargs)

    # Define samplers parameters
    mean_val = 1.0e-18
    std_dev = 1.0e-19
    min_val = 1.0e-19
    max_val = 1.0e-17
    seed = 2345
    perm_sampler1 = sm.add_component_model_object(
        SHPermeabilitySampler(name='ps1', parent=sm, grid_shape=(10, 10),
                              constr_type='matrix'))
    perm_sampler1.add_par('seed', value=seed, vary=False)
    perm_sampler1.add_par('avePermeability', value=mean_val, vary=False)
    perm_sampler1.add_par('stdDevPermeability', value=std_dev, vary=False)
    perm_sampler1.add_par('minPermeability', value=min_val, vary=False)
    perm_sampler1.add_par('maxPermeability', value=max_val, vary=False)
    perm_sampler1.add_grid_obs('permeability', constr_type='matrix', index=[0],
                               output_dir=output_dir)

    # Adding second sampler with the same parameters but heterogeneity approach set to True
    perm_sampler2 = sm.add_component_model_object(
        SHPermeabilitySampler(name='ps2', parent=sm, grid_shape=(10, 10),
                              constr_type='matrix', heterogeneity_approach=True))
    perm_sampler2.add_par('seed', value=seed, vary=False)
    perm_sampler2.add_par('avePermeability', value=mean_val, vary=False)
    perm_sampler2.add_par('stdDevPermeability', value=std_dev, vary=False)
    perm_sampler2.add_par('minPermeability', value=min_val, vary=False)
    perm_sampler2.add_par('maxPermeability', value=max_val, vary=False)
    perm_sampler2.add_par('heterFactor', value=0.5, vary=False)
    perm_sampler2.add_grid_obs('permeability', constr_type='matrix', index=[0],
                               output_dir=output_dir)

    sm.forward()

    # Collect saved observation
    data1 = sm.collect_gridded_observations_as_time_series(
        perm_sampler1, 'permeability', output_dir,
        indices=[0], rlzn_number=0)
    data2 = sm.collect_gridded_observations_as_time_series(
        perm_sampler2, 'permeability', output_dir,
        indices=[0], rlzn_number=0)

    print(data1[0].shape, data2[0].shape)
    data_min = min(np.min(data1[0]), np.min(data2[0]))
    data_max = max(np.max(data1[0]), np.max(data2[0]))

    import matplotlib.pyplot as plt
    fig1 = plt.figure(figsize=(17, 8))
    ax1 = fig1.add_subplot(121)
    im1 = ax1.pcolor(data1[0], vmin=data_min, vmax=data_max)
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    ax1.set_title('Original 100 cells')
    fig1.colorbar(im1, ax=ax1)
    ax2 = fig1.add_subplot(122)
    im2 = ax2.pcolor(data2[0], vmin=data_min, vmax=data_max)
    ax2.set_xlabel('x [m]')
    ax2.set_ylabel('y [m]')
    ax2.set_title('8 cells have updated permeability')
    plt.grid(False)
    fig1.colorbar(im2, ax=ax2)

    from mpl_toolkits.mplot3d import axes3d
    x = np.arange(1, 11)
    y = x
    xx, yy = np.meshgrid(x, y)
    xx = xx.reshape(100,)
    yy = yy.reshape(100,)
    zz = np.zeros(100,)

    dx = np.ones(100,)
    dy = np.ones(100,)
    dz1 = data1[0].reshape(100,)
    dz2 = data2[0].reshape(100,)

    fig2 = plt.figure(figsize=(16, 8))
    ax11 = fig2.add_subplot(121, projection='3d')
    ax11.bar3d(xx, yy, zz, dx, dy, dz1)
    ax11.set_xlabel('x [m]')
    ax11.set_ylabel('y [m]')
    ax11.set_zlabel(r'k [$m^2$]')
    ax11.set_title('Without heterogeneity')
    ax12 = fig2.add_subplot(122, projection='3d')
    ax12.bar3d(xx, yy, zz, dx, dy, dz2)
    ax12.set_xlabel('x [m]')
    ax12.set_ylabel('y [m]')
    ax12.set_zlabel(r'k [$m^2$]')
    ax12.set_title('With heterogeneity')


if __name__ == "__main__":

    test = 4
    if test == 1:
        test_scenario1()
    elif test == 2:
        test_scenario2()
    elif test == 3:
        test_scenario3()
    else:
        test_scenario4()
