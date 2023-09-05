"""
Example illustrates sampling of the lookup data file index rather
than the corresponding set of parameters for lookup table reservoir component.
Example also shows two ways of obtaining data from the lookup table reservoir component
at several points of interest.

This example requires the additional Kimberlina data set.
Kimberlina data set (pressure-in-pa-kimb-54-sims.zip or
Pressure_In_Pa_Kimb_54_sims.zip) can be downloaded from one of the following places:
1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam
2. https://gitlab.com/NRAP/Kimberlina_data

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/Kimb_54_sims

Example of run:
$ python iam_sys_lutreservoir_discrete_pars.py
"""

import sys
import os
import logging
import time
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, LookupTableReservoir


if __name__ == '__main__':
    start_time = time.time()
    # For multiprocessing in Spyder
    __spec__ = None
    logging.basicConfig(level=logging.WARNING)

    file_directory = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                  'lookuptables', 'Kimb_54_sims'])

    if not os.path.exists(os.sep.join([file_directory, 'Reservoir_data_sim01.csv'])):
        msg = ''.join(['\nKimberlina data set can be downloaded ',
                       'from one of the following places:\n',
                       '1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam \n',
                       '2. https://gitlab.com/NRAP/Kimberlina_data \n',
                       'Check this example description for more information.'])
        logging.error(msg)

    locX = [37478.0, 40360.3, 36757.1]
    locY = [48333.0, 52035.9, 46481.5]
    num_locs = len(locX)
    num_years = 50

    # Define keyword arguments of the system model
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Read file with signatures of interpolators and names of files with the corresponding data
    signature_data = np.genfromtxt(
        os.path.join(file_directory, 'parameters_and_filenames.csv'),
        delimiter=",", dtype='str')

    # The first row (except the last element) of the file contains names of the parameters.
    # The last element of the first row is supposed to be named 'filename'
    # The last column is for the filenames containing data for a particular realization
    par_names = signature_data[0, 1:-1]
    num_pars = len(par_names)

    # Add reservoir component
    ltres = sm.add_component_model_object(LookupTableReservoir(
        name='ltres', parent=sm, locX=locX, locY=locY))

    # It makes sense to use build_on_the_fly=False if during the simulation
    # all linked lookup tables will be used. If only subset of the lookup tables
    # will be used then user might consider using option build_on_the_fly=True
    # Option build_on_the_fly=True can be used selectively for particular interpolators
    # when interpolators are created before creating reservoir component
    # when it's known in advance which subset of lookup tables will be used.
    # That is, build_on_the_fly=True option should be used for interpolators which
    # won't be used, and build_on_the_fly=False option should be used
    # for the interpolators which will be used.
    ltres.build_and_link_interpolators(
        file_directory=file_directory, intr_family='reservoir',
        default_values={'salinity': 0.1, 'temperature': 50.0},
        recompute_triangulation=False, build_on_the_fly=False)

    # Add parameter of reservoir component model
    ltres.add_par('index', value=5, discrete_vals=(list(range(3, 12)), 9*[1/9]))
#    ltres.add_par('index', value=5, discrete_vals=(list(range(1, 55)), 54*[1/54]))

    # Add observations of reservoir component model
    ltres.add_grid_obs('pressure', constr_type='array', output_dir='test_output')
    ltres.add_grid_obs('CO2saturation', constr_type='array', output_dir='test_output')

    # We want to compare results returned by an individual reservoir component
    # and three grouped reservoir components: they should be the same
    # Initialize list of reservoir components
    ltress = []
    for ind in range(num_locs):
        # Add reservoir component
        # Interpolator family 'reservoir' was created above when we called method
        # build_and_link_interpolators on the ltres object
        ltress.append(sm.add_component_model_object(
            LookupTableReservoir(
                name='ltress'+str(ind+1), parent=sm,
                intr_family='reservoir', locX=locX[ind], locY=locY[ind])))

        # Add parameter of reservoir component model from a group of components
        ltress[-1].add_par_linked_to_par('index', ltres.pars['index'])

    # For convenience of post-processing we're adding first the pressure observations
    # then saturation observations
    for ind in range(len(locX)):
        ltress[ind].add_obs('pressure')
    for ind in range(len(locX)):
        ltress[ind].add_obs('CO2saturation')

    # Run system model using current values of its parameters
    sm.forward()

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    # Collect pressure and saturation observations from saved files
    saved_pressure = np.zeros((len(time_array), 3))
    saved_saturation = np.zeros((len(time_array), 3))
    grid_obs_ind = range(len(time_array))
    for ind in grid_obs_ind:
        filename = 'ltres_pressure_sim_0_time_{}.npz'.format(ind)
        d = np.load(os.sep.join(['test_output', filename]))
        saved_pressure[ind, :] = d['data']
        d.close()

        filename = 'ltres_CO2saturation_sim_0_time_{}.npz'.format(ind)
        d = np.load(os.sep.join(['test_output', filename]))
        saved_saturation[ind, :] = d['data']
        d.close()

    # Collect pressure and saturation observations from a group of components
    # Plot results
    fig, ax = plt.subplots(1, 3, figsize=(16, 4))
    for ind in range(len(locX)):
        pressure = sm.collect_observations_as_time_series(
            ltress[ind], 'pressure')
        saturation = sm.collect_observations_as_time_series(
            ltress[ind], 'CO2saturation')
        print('Pressure:', pressure, sep='\n')
        print('CO2saturation:', saturation, sep='\n')

        ax[0].plot(time_array/365.25, pressure/1.0e+6,
                   color="maroon", linewidth=1)
        ax[0].plot(time_array/365.25, saved_pressure/1.0e+6,
                   color="maroon", marker='x', linewidth=1)
        ax[0].set_xlabel('Time, t (years)', fontsize=14)
        ax[0].set_ylabel('Pressure, P (MPa)', fontsize=14)
        ax[0].set_title('Pressure: leaking well', fontsize=18)
        plt.tight_layout()
        plt.tick_params(labelsize=12)
        ax[0].set_xlim([0, 50])
        ax[0].get_yaxis().set_label_coords(-0.12, 0.5)

        ax[1].plot(time_array/365.25, saturation,
                   color="green", linewidth=1)
        ax[1].plot(time_array/365.25, saved_saturation,
                   color="green", marker='x', linewidth=1)
        ax[1].set_xlabel('Time, t (years)', fontsize=14)
        ax[1].set_ylabel('Saturation, S (-)', fontsize=14)
        ax[1].set_title(r'CO$_2$ saturation: leaking well', fontsize=18)
        plt.tight_layout()
        plt.tick_params(labelsize=12)
        ax[1].set_xlim([0, 50])
        ax[1].get_yaxis().set_label_coords(-0.12, 0.5)

        ax[2].plot(time_array/365.25, abs(saturation-saved_saturation[:, ind]),
                   color="green", linewidth=1)
        ax[2].plot(time_array/365.25, abs(pressure - saved_pressure[:, ind]),
                   color="maroon", linewidth=1)
        ax[2].set_xlabel('Time, t (years)', fontsize=14)
        ax[2].set_ylabel('Absolute difference', fontsize=14)
        ax[2].set_title('Difference between observations', fontsize=18)
        plt.tight_layout()
        plt.tick_params(labelsize=12)
        ax[2].set_xlim([0, 50])
        ax[2].get_yaxis().set_label_coords(-0.12, 0.5)

    print('------------------------------------------------------------------')
    print('                    LHS method illustration ')
    print('------------------------------------------------------------------')

    # Run stochastic simulation
    num_samples = 10
    ncpus = 2
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=450)

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)

    end_time = time.time()
    print('Passed time', end_time-start_time)

    # Compare results of a particular realization
    sim_num = 4
    # Collect pressure and saturation observations from saved files
    saved_pressure = np.zeros((len(time_array), 3))
    saved_saturation = np.zeros((len(time_array), 3))
    grid_obs_ind = range(len(time_array))
    for ind in grid_obs_ind:
        filename = 'ltres_pressure_sim_{}_time_{}.npz'.format(sim_num, ind)
        d = np.load(os.sep.join(['test_output', filename]))
        saved_pressure[ind, :] = d['data']
        d.close()

        filename = 'ltres_CO2saturation_sim_{}_time_{}.npz'.format(sim_num, ind)
        d = np.load(os.sep.join(['test_output', filename]))
        saved_saturation[ind, :] = d['data']
        d.close()

    # Collect pressure and saturation observations from a group of components
    # Plot results
    fig, ax = plt.subplots(1, 3, figsize=(16, 4))
    for ind in range(len(locX)):
        pressure = results[sim_num-1, ind*len(time_array):(ind+1)*len(time_array)]
        saturation = results[sim_num-1, (ind+3)*len(time_array):(ind+4)*len(time_array)]
        print('Pressure:', pressure, sep='\n')
        print('CO2saturation:', saturation, sep='\n')

        ax[0].plot(time_array/365.25, pressure/1.0e+6,
                   color="maroon", linewidth=1)
        ax[0].plot(time_array/365.25, saved_pressure/1.0e+6,
                   color="maroon", marker='x', linewidth=1)
        ax[0].set_xlabel('Time, t (years)', fontsize=14)
        ax[0].set_ylabel('Pressure, P (MPa)', fontsize=14)
        ax[0].set_title('Pressure: leaking well', fontsize=18)
        plt.tight_layout()
        plt.tick_params(labelsize=12)
        ax[0].set_xlim([0, 50])
        ax[0].get_yaxis().set_label_coords(-0.12, 0.5)

        ax[1].plot(time_array/365.25, saturation, color="green", linewidth=1)
        ax[1].plot(time_array/365.25, saved_saturation, color="green", marker='x', linewidth=1)
        ax[1].set_xlabel('Time, t (years)', fontsize=14)
        ax[1].set_ylabel('Saturation, S (-)', fontsize=14)
        ax[1].set_title(r'CO$_2$ saturation: leaking well', fontsize=18)
        plt.tight_layout()
        plt.tick_params(labelsize=12)
        ax[1].set_xlim([0, 50])
        ax[1].get_yaxis().set_label_coords(-0.12, 0.5)

        ax[2].plot(time_array/365.25, abs(saturation-saved_saturation[:, ind]),
                   color="green", linewidth=1)
        ax[2].plot(time_array/365.25, abs(pressure - saved_pressure[:, ind]),
                   color="maroon", linewidth=1)
        ax[2].set_xlabel('Time, t (years)', fontsize=14)
        ax[2].set_ylabel('Absolute difference', fontsize=14)
        ax[2].set_title('Difference between observations', fontsize=18)
        plt.tight_layout()
        plt.tick_params(labelsize=12)
        ax[2].set_xlim([0, 50])
        ax[2].get_yaxis().set_label_coords(-0.12, 0.5)
