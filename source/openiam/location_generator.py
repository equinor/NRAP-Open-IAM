"""
Created on Fri Sep 20 12:54:26 2019

"""
from math import ceil
import sys
import os
import logging
import numpy as np
from numpy.random import RandomState
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import SystemModel, SamplerModel
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))


def generate_locations(x_min, x_max, y_min, y_max, num_locations, seed=None,
                       z_min=None, z_max=None):
    """ Generate x- and y-coordinates uniformly within the defined domain.

    :param x_min: lower boundary of the domain along x-axis
    :type x_min: float

    :param x_max: upper boundary of the domain along x-axis
    :type x_max: float

    :param y_min: lower boundary of the domain along y-axis
    :type y_min: float

    :param y_max: upper boundary of the domain along y-axis
    :type y_max: float

    :param z_min: lower boundary of the domain along z-axis. None if not used.
    :type z_min: float

    :param z_max: upper boundary of the domain along z-axis. None if not used.
    :type z_max: float

    :param num_locations: number of locations to be generated
    :type num_locations: int

    :param seed: value of the seed argument to create pseudo-random number
        generator. By default, the value is None: this allows to produce the random
        sequence every time the function is called.
    :type seed: float

    :returns: 2 numpy.arrays of shape (num_locations,) with
        x- and y-coordinates of the generated locations
    """
    if seed is None:
        prng = RandomState()  # pseudo-random number generator
    else:
        prng = RandomState(ceil(seed))  # using ceil since seed has to be an integer

    # Generate uniform random variables and transform according to the defined setup
    x_coords = x_min + (x_max - x_min) * prng.rand(num_locations)
    y_coords = y_min + (y_max - y_min) * prng.rand(num_locations)

    if z_min is not None and z_max is not None:
        z_coords = z_min + (z_max - z_min) * prng.rand(num_locations)
    else:
        z_coords = None

    return x_coords, y_coords, z_coords


class LocationGenerator(SamplerModel):
    """ NRAP-Open-IAM LocationGenerator component class. """
    def __init__(self, name, parent, x_min=0, x_max=1000, y_min=0, y_max=1000,
                 num_locations=1, reproducible=True, z_min=None, z_max=None):
        """
        Constructor method of LocationGenerator class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param x_min: minimal x-coordinate (in meters) of the x-interval along
            which the (well) locations should be placed
        :type x_min: float

        :param x_max: maximal x-coordinate (in meters) of the x-interval along
            which the (well) locations should be placed
        :type x_max: float

        :param y_min: minimal y-coordinate (in meters) of the y-interval along
            which the (well) locations should be placed
        :type y_min: float

        :param y_max: maximal x-coordinate (in meters) of the y-interval along
            which the (well) locations should be placed
        :type y_max: float

        :param num_locations: number of (well) locations to be generated
        :type num_locations: int

        :param reproducible: flag variable indicating whether the produced output
            should be reproducible, i.e. if seed parameter of the model method
            should be used to generate the sequence of random locations.
        :type reproducible: boolean

        :returns: LocationGenerator class object
        """
        super().__init__(name, parent, reproducible=reproducible,
                         default_seed_value=3)

        # Add type attribute
        self.class_type = 'LocationGenerator'

        # Save setup attributes of the component
        self.check_ranges(num_locations, x_min, x_max, y_min, y_max,
                          z_min=z_min, z_max=z_max)

        msg = 'LocationGenerator created with name {}'.format(self.name)
        logging.debug(msg)

    def invalid_ranges_msg(self, var, min_val, max_val):
        msg = ''.join([
            'Provided minimum value ({}) for {}-coordinates is greater than ',
            'the provided maximum value ({}) for LocationGenerator {}. ',
            'Please check your setup.']).format(min_val, var, max_val, self.name)
        logging.error(msg)
        raise ValueError(msg)

    def check_ranges(self, num_locations, x_min, x_max, y_min, y_max, z_min, z_max):

        if x_min <= x_max:
            self.x_min = x_min
            self.x_max = x_max
        else:
            self.invalid_ranges_msg('x', x_min, x_max)
        if y_min <= y_max:
            self.y_min = y_min
            self.y_max = y_max
        else:
            self.invalid_ranges_msg('y', y_min, y_max)

        self.num_locations = num_locations

        # The z coordinates produced will be None if z_min and/or z_max are None.
        if z_min is not None and z_max is not None:
            if z_min <= z_max:
                self.z_min = z_min
                self.z_max = z_max
            else:
                self.invalid_ranges_msg('z', z_min, z_max)
        else:
            self.z_min = None
            self.z_max = None

    def sample(self, p):
        """ Return x-, y-, and z-coordinates of well locations."""
        # Obtain the default values of the parameters from dictionary
        # of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # The generated 'random' (well) locations will be the same
        # for the runs with the same seed and on the same machine
        if self.reproducible:    # if reproducible results are needed
            x_coords, y_coords, z_coords = generate_locations(
                self.x_min, self.x_max, self.y_min, self.y_max,
                self.num_locations, seed=actual_p['seed'],
                z_min=self.z_min, z_max=self.z_max)
        else:
            x_coords, y_coords, z_coords = generate_locations(
                self.x_min, self.x_max, self.y_min, self.y_max,
                self.num_locations, z_min=self.z_min, z_max=self.z_max)

        if z_coords is None:
            out = {'locX': x_coords, 'locY': y_coords}
        else:
            out = {'locX': x_coords, 'locY': y_coords, 'locZ': z_coords}

        return out


def test_case1():
    # For multiprocessing in Spyder
    __spec__ = None

    from openiam import SimpleReservoir

    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Define number of wells
    num_wells = 5
    # Define boundaries of random wells domains
    x_min = 100.0
    x_max = 300.0
    y_min = 150.0
    y_max = 300.0
    # Add generator component
    gen = sm.add_component_model_object(LocationGenerator(
        name='gen', parent=sm, x_min=x_min, x_max=x_max,
        y_min=y_min, y_max=y_max, num_locations=num_wells,
        reproducible=True))

    # Add parameters of generator component
    gen.add_par('seed', value=11, min=3, max=12345)
    gen.add_obs_to_be_linked('locX', obs_type='grid')
    gen.add_obs_to_be_linked('locY', obs_type='grid')
    gen.add_local_obs('locX1', grid_obs_name='locX',
                      constr_type='array', loc_ind=0, index=[0])
    gen.add_local_obs('locY1', grid_obs_name='locY',
                      constr_type='array', loc_ind=0, index=[0])
    # Locations can be obtained as observations of the generator component
    # Since the observations considered to be of gridded type the outputs
    # can only be obtained by adding them using add_obs method and
    # then reading the files with data after simulation is complete

    # Add simple reservoir component
    # We don't specify the location for the reservoir component, although
    # there is a default reservoir setup with locX=100, locY=100
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))
    # Add parameters of reservoir component model
    # All parameters of the reservoir component are deterministic
    # so all uncertainty in the simulation comes from the uncertainty
    # of the well location
    sres.add_par('numberOfShaleLayers', value=3, vary=False)
    sres.add_par('injRate', value=0.5, vary=False)
    sres.add_par('shale1Thickness', value=40.0, vary=False)
    sres.add_par('shale2Thickness', value=50.0, vary=False)

    # Simple reservoir component has keyword arguments of the model method:
    # locX and locY which we would link to the output of the generator component
    # Generator component outputs 5 random locations according to the setup above.
    # We can link reservoir component locX and locY to any of these produced locations
    # using arguments constr_type='array' and loc_ind=[2]; 2 is chosen arbitrarily
    sres.add_kwarg_linked_to_obs('locX', gen.linkobs['locX'],
                                 obs_type='grid', constr_type='array', loc_ind=[2])
    sres.add_kwarg_linked_to_obs('locY', gen.linkobs['locY'],
                                 obs_type='grid', constr_type='array', loc_ind=[2])

    # Add observations of reservoir component model
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')

    # Run forward simulation
    sm.forward()

    print('__________________Results of forward simulation___________________')
    # Print pressure and saturation which should be different every time the code is run
    # if attribute reproducible is False
    print('Pressure',
          sm.collect_observations_as_time_series(sres, 'pressure'), sep='\n')
    print('CO2 saturation',
          sm.collect_observations_as_time_series(sres, 'CO2saturation'), sep='\n')

    num_samples = 100
    ncpus = 5
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=458)

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)

    print('__________________Results of stochastic simulations_______________')

    # Extract results from stochastic simulations
    pressure = np.ones((num_samples, len(time_array)))
    CO2_saturation = np.ones((num_samples, len(time_array)))

    for ind in range(len(time_array)):
        pressure[:, ind] = s.recarray['sres.pressure_{}'.format(ind)]
        CO2_saturation[:, ind] = s.recarray['sres.CO2saturation_{}'.format(ind)]

    print('Average pressure', np.mean(pressure, axis=0), sep='\n')
    print('Standard deviation of pressure', np.std(pressure, axis=0), sep='\n')
    print('Average CO2 saturation', np.mean(CO2_saturation, axis=0), sep='\n')
    print('Standard deviation of CO2 saturation', np.std(CO2_saturation, axis=0), sep='\n')

    # Plot results
    fontsize = 16
    fig = plt.figure(figsize=(13, 5))
    ax = fig.add_subplot(121)
    for j in range(num_samples):
        plt.plot(time_array/365.25, pressure[j]/1.0e+6, color="maroon", linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=fontsize-2)
    plt.ylabel('Pressure, P (MPa)', fontsize=fontsize-2)
    plt.title('Pressure: Randomly allocated leaking well', fontsize=fontsize)
    plt.tight_layout()
    plt.tick_params(labelsize=fontsize-4)
    plt.xlim([0, 50])
    plt.ylim([20, 45])
    ax.get_yaxis().set_label_coords(-0.10, 0.5)

    ax = fig.add_subplot(122)
    for j in range(num_samples):
        plt.plot(time_array/365.25, CO2_saturation[j], color="green", linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=fontsize-2)
    plt.ylabel('Saturation, S (-)', fontsize=fontsize-2)
    plt.title(r'CO$_2$ saturation: Randomly allocated leaking well', fontsize=fontsize)
    plt.tight_layout()
    plt.tick_params(labelsize=fontsize-4)
    plt.xlim([0, 50])
    plt.ylim([0.0, 1.0])
    ax.get_yaxis().set_label_coords(-0.10, 0.5)

if __name__ == "__main__":

    test_case1()
