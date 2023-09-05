# -*- coding: utf-8 -*-
import os
import sys
import logging
import numpy as np
from numpy.random import default_rng
np.set_printoptions(threshold=np.inf)
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

try:
    from openiam import SystemModel, SamplerModel
except ImportError as err:
    print('Unable to load IAM class module: '+ str(err))

try:
    import components.seal.seal_refresh as sref
    import components.seal.seal_intro as intro
    from components.seal.seal_setup import SEAL_SETUP_DICT
except ImportError:
    print('\nERROR: Unable to load Thickness Sampler for Seal Horizon\n')
    sys.exit()

SHTS_PARAMETERS = ['aveThickness', 'stdDevThickness',
                   'minThickness', 'maxThickness']
SCTS_PARAMETERS = ['thickness_ave', 'thickness_std',
                   'thickness_min', 'thickness_max']  # seal_controls pars


class SHThicknessSampler(SamplerModel):
    """
    Class for model providing thickness parameter for Seal Horizon component.

    The class can be used as provider of thickness parameter values
    for any other component.

    Distribution parameters for thickness of the seal horizon layer

    * **aveThickness** [|m|] (10 to 1000) - mean of the truncated normal
      distribution for thickness (default: 100)

    * **stdDevThickness** [|m|] (0 to 500) - standard deviation
      of the thickness distribution (default: 20)

    * **minThickness** [|m|] (5 to 100) - minimum thickness; this value
      truncates the distribution and limits lower values (default: 75)

    * **maxThickness** [|m|] (10 to 1000) - maximum thickness; this value
      truncates the distribution and limits higher values (default: 125).

    The possible output from the SH Thickness Sampler component is

    * **thickness** [|m|] - thickness of a layer or cells over a given area.
    """

    def __init__(self, name, parent, grid_shape, constr_type='array', coordinates=None,
                 reproducible=True):
        """
        Constructor method of SHThicknessSampler class

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

        :returns: SHThicknessSampler class object
        """
        super().__init__(name=name, parent=parent, reproducible=reproducible)

        # Add type attribute
        self.class_type = 'SHThicknessSampler'

        # Define gridded observations names
        self.grid_obs_keys = ['thickness']

        self.controls = {}
        # Add default parameters specific to this sampler
        # Note that seed parameter is already parameter of the sampler
        for ind, par_name in enumerate(SHTS_PARAMETERS):
            self.add_default_par(
                par_name, value=SEAL_SETUP_DICT['Thickness'][par_name]['valu'])
            self.controls[SCTS_PARAMETERS[ind]] = SEAL_SETUP_DICT['Thickness'][
                par_name]['valu']

        # Define parameters boundaries dictionary self.pars_bounds
        self.define_pars_bounds()

        # Logging output of default parameters and their values
        debug_msg = ''.join([par_nm+': '+str(par_obj.value) \
            for par_nm, par_obj in self.default_pars.items()])
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
        # and it's assumed that if the object of SHThicknessSampler class is
        # created then it's observation will be used by some other component.
        # Otherwise, beyond test purposes there are no other reasons to create
        # the object
        self.add_obs_to_be_linked(
            'thickness', obs_type='grid', constr_type=self.controls['constr_type'],
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
        self.pars_bounds['aveThickness'] = param_bounds['thickness_ave']
        self.pars_bounds['stdDevThickness'] = param_bounds['thickness_std']
        self.pars_bounds['minThickness'] = param_bounds['thickness_min']
        self.pars_bounds['maxThickness'] = param_bounds['thickness_max']

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
        for key in SHTS_PARAMETERS:
            if key in p:
                val = p[key]
                if ((val < self.pars_bounds[key][0]) or
                        (val > self.pars_bounds[key][1])):
                    warn_msg = ''.join([
                        'Parameter {} of value {} of SH Thickness Sampler {} ',
                        'is out of boundaries [{}, {}].']).format(
                            key, val, self.name, self.pars_bounds[key][0],
                            self.pars_bounds[key][1])
                    logging.warning(warn_msg)

    def sample(self, p):
        """
        :param p: input parameters of SHThicknessSampler model
        :type p: dict

        """
        # Obtain the default values of the parameters from dictionary of
        # default parameters.
        actual_p = {k: v.value for k, v in self.default_pars.items()}

        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Update controls with parameter values
        for ind, par_name in enumerate(SCTS_PARAMETERS):
            self.controls[par_name] = actual_p[SHTS_PARAMETERS[ind]]

        # Define output
        if self.reproducible:
            thickness = sref.thickness_variability(
                self.controls, actual_p['seed'])
        else:
            rng = default_rng()
            seed = rng.integers(1, 123456789)
            thickness = sref.thickness_variability(self.controls, seed)
        # Check whether matrix-type output is needed
        if self.controls['constr_type'] == 'matrix':
            thickness = thickness.reshape(self.controls['grid_rows'],
                                          self.controls['grid_cols'])

        out = {'thickness': thickness}

        # Return output
        return out


def test_scenario1():
    """ Run script example illustrating work of SH Thickness Sampler model.

    """
    # Setup logging and constants.
    logging.basicConfig(level=logging.DEBUG)

    # Define keyword arguments of the system model.
    model_kwargs = {'time_point': 0.0} # time is given in days

    output_dir = '../../../output/samplers/'
    # Create system model.
    sm = SystemModel(model_kwargs=model_kwargs)
    t_sampler = sm.add_component_model_object(
        SHThicknessSampler(name='ts', parent=sm, grid_shape=(3, 3)))
    t_sampler.add_par('seed', value=234, vary=False)
    t_sampler.add_par('aveThickness', value=500, vary=False)
    t_sampler.add_par('stdDevThickness', value=50, vary=False)
    t_sampler.add_par('minThickness', value=450, vary=False)
    t_sampler.add_par('maxThickness', value=575, vary=False)
    t_sampler.add_grid_obs('thickness', constr_type='array', index=[0],
                           output_dir=output_dir)

    sm.forward()

    # Collect saved observation
    data = sm.collect_gridded_observations_as_time_series(
        t_sampler, 'thickness', output_dir,
        indices=[0], rlzn_number=0)

    print(data[0].shape)
    print(data[0])


def test_scenario2():
    """ Run script example illustrating work of SH Thickness Sampler model.

    """
    # Setup logging and constants.
    logging.basicConfig(level=logging.DEBUG)

    # Define keyword arguments of the system model.
    model_kwargs = {'time_point': 0.0} # time is given in days

    output_dir = '../../../output/samplers/'
    # Create system model.
    sm = SystemModel(model_kwargs=model_kwargs)
    t_sampler = sm.add_component_model_object(
        SHThicknessSampler(name='ts', parent=sm, grid_shape=(10, 10),
                           constr_type='matrix'))
    t_sampler.add_par('seed', value=234, vary=False)
    t_sampler.add_par('aveThickness', value=500, vary=False)
    t_sampler.add_par('stdDevThickness', value=50, vary=False)
    t_sampler.add_par('minThickness', value=450, vary=False)
    t_sampler.add_par('maxThickness', value=575, vary=False)
    t_sampler.add_grid_obs('thickness', constr_type='matrix', index=[0],
                           output_dir=output_dir)

    sm.forward()

    # Collect saved observation
    data = sm.collect_gridded_observations_as_time_series(
        t_sampler, 'thickness', output_dir,
        indices=[0], rlzn_number=0)

    print(data[0].shape)
    print(data[0])

    x = np.arange(1, 11)
    y = x
    xx, yy = np.meshgrid(x, y)
    xx = xx.reshape(100,)
    yy = yy.reshape(100,)
    plt.figure(figsize=(5, 5))
    plt.scatter(xx, yy, c=data[0].reshape(100, ), marker='s', s=200)
    plt.title('Layer thickness')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')

    plt.figure(figsize=(8, 8))
    plt.imshow(np.flipud(data[0]))
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.colorbar()


def test_scenario3():
    """ Run script example illustrating work of SH Thickness Sampler model.

    """
    # Setup logging and constants.
    logging.basicConfig(level=logging.DEBUG)

    # Define keyword arguments of the system model.
    model_kwargs = {'time_point': 0.0} # time is given in days

    output_dir = '../../../output/samplers/'
    # Create system model.
    sm = SystemModel(model_kwargs=model_kwargs)
    t_sampler = sm.add_component_model_object(
        SHThicknessSampler(name='ts', parent=sm, grid_shape=(100,)))
    t_sampler.add_par('seed', value=234, vary=False)
    t_sampler.add_par('aveThickness', value=500, vary=False)
    t_sampler.add_par('stdDevThickness', value=12, vary=False)
    t_sampler.add_par('minThickness', value=450, vary=False)
    t_sampler.add_par('maxThickness', value=575, vary=False)
    t_sampler.add_grid_obs('thickness', constr_type='array', index=[0],
                           output_dir=output_dir)

    sm.forward()

    # Collect saved observation
    data = sm.collect_gridded_observations_as_time_series(
        t_sampler, 'thickness', output_dir,
        indices=[0], rlzn_number=0)

    print(data[0].shape)
    print(data[0])

    x = np.arange(1, 101)
    plt.figure(figsize=(5, 5))
    plt.plot(x, data[0], '-o')
    plt.title('Layer thickness')
    plt.xlabel('x [m]')
    plt.ylabel('h [m]')


if __name__ == "__main__":

    test = 1
    if test == 1:
        test_scenario1()
    elif test == 2:
        test_scenario2()
    else:
        test_scenario3()
