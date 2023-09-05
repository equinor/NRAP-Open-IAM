# -*- coding: utf-8 -*-
"""
Example illustrates system model with two linked simple components.
The first component returns gridded observation utilized by the second component.
Example also illustrates forward and lhs instance methods of SystemModel class.

Example of run:
$ python iam_sys_cmpnt1_cmpnt2_gridded_obs.py
"""

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel

def comp_model1(p, x=None, y=None, time_point=0.0):
    # Check whether keyword arguments x and y were provided
    if (x is not None) and (y is not None):
        xx = x
        yy = y
    else:
        N = 5
        xx = np.linspace(0., 10., N+1)
        yy = np.linspace(0., 10., N+1)

    # Check which parameters are provided in the input dictionary
    if 'var1' in p:
        var1 = p['var1']
    else:
        var1 = 2.
    if 'var2' in p:
        var2 = p['var2']
    else:
        var2 = 1.0

    # Setup output dictionary
    out = {}
    out['plane'] = np.zeros((len(x)*len(y),))
    ind = 0
    for xv in xx:
        for yv in yy:
            out['plane'][ind] = np.log10(time_point+1)*(xv+yv)/(10*var1) + \
                var2*np.log10(time_point+1)
            ind = ind + 1
    out['time_squared'] = var2*np.sin(time_point)**2/1.0e+3 + var1*time_point/10.
    return out

def comp_model2(p, time_point=0.0, kw1=1.0, kw2=None, kw3=1.0, kw4=None):

    # Check keyword arguments
    if kw2 is None:
        kw2 = [1.0]

    if kw4 is None:
        kw4 = [1.3, 1.5]

    # Check which parameters are provided in the input dictionary
    if 'var1' in p:
        var1 = p['var1']
    else:
        var1 = -1.0  # default values if the desired parameter value was not passed

    if 'var2' in p:
        var2 = p['var2']
    else:
        var2 = 2.0

    if 'var3' in p:
        var3 = p['var3']
    else:
        var3 = -3.0

    # Define output of the function (output is also a dictionary)
    output = dict()
    output['output1'] = np.abs((kw1*var1**2 + kw3*var2*var3)*np.sin(time_point))
    output['output2'] = np.abs(
        kw2[0]*var2**2 + kw4[0]/100.+ kw3*var1*var3*np.cos(time_point))

    # Component model should return a dictionary of outputs
    return output


if __name__ == "__main__":
    # For multiprocessing in Spyder
    __spec__ = None
    # Define keyword arguments of the system model
    num_years = 10
    time_array = 365.25*np.arange(0.0, num_years+1, step=1./4.)
    num_tp = len(time_array)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add component 1
    N = 100
    x = np.linspace(0., 10., N+1)

    cm1 = sm.add_component_model(name='cm1', model=comp_model1,
                                 model_kwargs={'x': x, 'y': x, 'time_point': 365.25})
    cm1.add_par('var1', min=2., max=5., value=3.)
    cm1.add_par('var2', min=1., max=6.5, value=4.)

    # Define indices of time points at which the gridded observation should be saved
    grid_obs_ind = list(range(0, 9))

    # Add gridded observations
    grid_save_type = 'csv'
    cm1.add_grid_obs('plane', constr_type='array',
                     output_dir='test_output',
                     index=grid_obs_ind, save_type=grid_save_type)

    # Add local (scalar) observations
    cm1.add_local_obs('plane_loc1', 'plane', constr_type='array', loc_ind=65)
    cm1.add_local_obs('plane_loc2', 'plane', constr_type='array', loc_ind=74)
    cm1.add_obs('time_squared')

    # In order to be used as arguments of the subsequent components,
    # observations should be added using add_obs_to_be_linked method
    # of the component providing the observations
    cm1.add_obs_to_be_linked('time_squared')
    cm1.add_obs_to_be_linked('plane', obs_type='grid')

    # Add component 2
    cm2 = sm.add_component_model(name='cm2', model=comp_model2,
                                 model_kwargs={'time_point': 0.0})
    cm2.add_par('var1', value=1., vary=False)
    cm2.add_par('var2', value=1., vary=False)
    cm2.add_par('var3', value=1., vary=False)

    # Link keyword arguments of component 2 to the respective observations of component 1
    cm2.add_kwarg_linked_to_obs('kw1', cm1.linkobs['time_squared'])
    cm2.add_kwarg_linked_to_obs('kw2', cm1.linkobs['plane'], obs_type='grid')
    cm2.add_kwarg_linked_to_obs('kw3', cm1.linkobs['plane'], obs_type='grid',
                                constr_type='array', loc_ind=[45])
    cm2.add_kwarg_linked_to_obs('kw4', cm1.linkobs['plane'], obs_type='grid',
                                constr_type='array', loc_ind=[23, 33, 43, 53])
    cm2.add_obs('output1')
    cm2.add_obs('output2')

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    sm.forward()

    # Check saved files. The following piece of code provides example of
    # how to read the content of file with gridded observation
    print('Selected content of saved files')
    data = {}

    for ind, val in enumerate(grid_obs_ind):
        # Load gridded observations based on save data type
        filename = '_'.join(['cm1', 'plane', 'sim_0', 'time_'+str(val)])+'.'+grid_save_type
        file_path = os.sep.join(['test_output', filename])
        d = sm.get_gridded_observation_file(file_path)
        data[ind] = d[64]
        print('Data from file ' + filename, data[ind])

    # Extract data using collect_gridded_observations_as_time_series method
    data_to_comp = sm.collect_gridded_observations_as_time_series(
        cm1, 'plane', os.sep.join([os.getcwd(), 'test_output']),
        indices=grid_obs_ind, save_type=grid_save_type)

    # Compare data
    for ind in range(len(grid_obs_ind)):
        print('Data extracted from file:', data[ind])
        print('Data extracted with the method:', data_to_comp[ind, 64])

    plane_loc1 = sm.collect_observations_as_time_series(cm1, 'plane_loc1')
    plane_loc2 = sm.collect_observations_as_time_series(cm1, 'plane_loc2')
    time_squared = sm.collect_observations_as_time_series(cm1, 'time_squared')
    print('------------------------------------------------------------------')
    print('Observations of the first component')
    print('------------------------------------------------------------------')
    print(plane_loc1, plane_loc2, time_squared, sep='\n')

    output1 = sm.collect_observations_as_time_series(cm2, 'output1')
    output2 = sm.collect_observations_as_time_series(cm2, 'output2')
    print('------------------------------------------------------------------')
    print('Observations of the second component')
    print('------------------------------------------------------------------')
    print(output1, output2, sep='\n')

    print('------------------------------------------------------------------')
    print('                          UQ illustration ')
    print('------------------------------------------------------------------')

    import random
    num_samples = 15
    ncpus = 5
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=random.randint(500, 1100))

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)

    fig1, ax1 = plt.subplots(nrows=2, ncols=2, figsize=(12, 7), sharex=True)
    for j in range(num_samples):
        ax1[0, 0].semilogy(time_array/365.25, results[j, 0:(1*num_tp)],
                           color='#551166', linewidth=1)
    ax1[0, 0].set_ylabel('Obs: plane_loc1', fontsize=14)

    for j in range(num_samples):
        ax1[0, 1].semilogy(time_array/365.25, results[j, (1*num_tp):(2*num_tp)],
                           color='#AA0022', linewidth=1)
    ax1[0, 1].set_ylabel('Obs: plane_loc2', fontsize=14)

    for j in range(num_samples):
        ax1[1, 0].semilogy(time_array/365.25, results[j, (2*num_tp):(3*num_tp)],
                           color='#000066', linewidth=1)
    ax1[1, 0].set_xlabel('Time, t [years]', fontsize=14)
    ax1[1, 0].set_ylabel('Obs: time_squared', fontsize=14)

    for j in range(num_samples):
        ax1[1, 1].semilogy(time_array/365.25, results[j, (3*num_tp):(4*num_tp)],
                           color='#00ABBB', linewidth=1)
    ax1[1, 1].set_xlabel('Time, t [years]', fontsize=14)
    ax1[1, 1].set_ylabel('Obs: output1', fontsize=14)

    # Extract data for gridded observation plane and stochastic simulations
    data = {}
    for real_num in range(1, 6):
        data[real_num] = sm.collect_gridded_observations_as_time_series(
            cm1, 'plane', os.sep.join([os.getcwd(), 'test_output']),
            indices=grid_obs_ind, rlzn_number=real_num, save_type=grid_save_type)

    from mpl_toolkits.mplot3d import Axes3D
    fig2 = plt.figure(figsize=(12, 7))
    ax2 = fig2.add_subplot(111, projection='3d')
    colors = ['blue', 'red']
    ds = data[1].shape[1]
    sub_indices = np.arange(0, ds, 500)
    for real_num in range(1, 2):
        for t_ind in grid_obs_ind[1:]:
            ax2.plot(time_array[t_ind]*np.ones(len(sub_indices))/365.25,
                     np.arange(1, ds+1)[sub_indices],
                     data[real_num][t_ind][sub_indices],
                     '-', color=colors[real_num-1], linewidth=10)
    ax2.set_xlabel('Time, t [years]', fontsize=14)
    ax2.set_ylabel('Location index, [-]', fontsize=14)
    ax2.set_zlabel('Observation plane, [-]', fontsize=14)
