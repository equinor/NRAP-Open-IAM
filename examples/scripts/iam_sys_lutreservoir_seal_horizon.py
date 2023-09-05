"""
Example illustrates simple linking of lookup table reservoir and seal horizon component.

This example requires the additional Kimberlina data set.
Kimberlina data set (pressure-in-pa-kimb-54-sims.zip or
Pressure_In_Pa_Kimb_54_sims.zip) can be downloaded from one of the following places:
1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam
2. https://gitlab.com/NRAP/Kimberlina_data

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/Kimb_54_sims

Example of run:
$ python iam_sys_lutreservoir_seal_horizon.py
"""
import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, LookupTableReservoir, SealHorizon

if __name__ == "__main__":

    # Setup location of lookup table data set
    file_directory = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                  'lookuptables', 'Kimb_54_sims'])

    # Setup save type for gridded observations
    grid_save_type = 'npz'

    if not os.path.exists(os.sep.join([file_directory, 'Reservoir_data_sim01.csv'])):
        msg = ''.join(['\nKimberlina data set can be downloaded ',
                       'from one of the following places:\n',
                       '1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam \n',
                       '2. https://gitlab.com/NRAP/Kimberlina_data  \n',
                       'Check this example description for more information.'])
        logging.error(msg)

    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25 * np.arange(0.0, num_years+1)
    seal_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=seal_model_kwargs)

    # Read file with signatures of interpolators and names of files
    # with the corresponding data.
    signature_data = np.genfromtxt(
        os.path.join(file_directory, 'parameters_and_filenames.csv'),
        delimiter=",", dtype='str')

    # Setup locations at which pressure and saturation will be calculated
    loc_X = [35315.8, 36036.4, 36757.1, 37477.7, 38198.4, 38919, 39639.7,
             40360.3, 41081, 41801.6, 42522.3, 43242.9, 43963.6, 44684.2, 45404.8]

    loc_Y = 15*[51110.2]

    # Add reservoir component
    ltres = sm.add_component_model_object(
        LookupTableReservoir(name='ltres', parent=sm, locX=loc_X, locY=loc_Y))

    ltres.build_and_link_interpolators(file_directory=file_directory,
                                       intr_family='reservoir',
                                       default_values={'salinity': 0.1,
                                                       'temperature': 50.0},
                                       recompute_triangulation=False,
                                       build_on_the_fly=False)
    # Add parameter of the reservoir component indicating what lookup table
    # will be used
    ltres.add_par('index', value=5, vary=False)

    # Define output folder to keep data files with gridded observations
    output_dir = 'sim_data'

    # Add gridded observations of the reservoir component
    ltres.add_grid_obs('pressure', constr_type='array',
                       output_dir=output_dir, save_type=grid_save_type)
    ltres.add_grid_obs('CO2saturation', constr_type='array',
                       output_dir=output_dir, save_type=grid_save_type)

    # Add observations of reservoir component model to be used as input
    # of the SealHorizon component
    ltres.add_obs_to_be_linked('pressure', obs_type='grid',
                               constr_type='array')
    ltres.add_obs_to_be_linked('CO2saturation', obs_type='grid',
                               constr_type='array')

    # Setup coordinates of the cell centers
    cell_coord_x = loc_X
    cell_coord_y = loc_Y
    num_cells = len(cell_coord_x)

    shc = sm.add_component_model_object(
        SealHorizon(name='shc', parent=sm,
                    locX=cell_coord_x, locY=cell_coord_y, area=6.3e+4))

    # Add time varying keyword arguments of the component
    shc.add_kwarg_linked_to_obs('pressure', ltres.linkobs['pressure'],
                                obs_type='grid', constr_type='array')
    shc.add_kwarg_linked_to_obs('CO2saturation',
                                ltres.linkobs['CO2saturation'],
                                obs_type='grid', constr_type='array')

    # Add gridded observations of the SealHorizon component
    shc.add_grid_obs('CO2_aquifer', constr_type='array',
                     output_dir=output_dir, save_type=grid_save_type)
    shc.add_grid_obs('brine_aquifer', constr_type='array',
                     output_dir=output_dir, save_type=grid_save_type)

    # Add scalar observations of the SealHorizon component
    shc.add_obs('CO2_aquifer_total')
    shc.add_obs('brine_aquifer_total')

    print('Starting simulation...')
    # Run forward simulation
    sm.forward()
    print('Simulation is finished.')

    # Collect scalar observations
    total_CO2_rate = sm.collect_observations_as_time_series(
        shc, 'CO2_aquifer_total')
    print('Cumulative CO2 rate at all 15 locations:', total_CO2_rate, sep='\n')

    total_brine_rate = sm.collect_observations_as_time_series(
        shc, 'brine_aquifer_total')
    print('Cumulative brine rate at all 15 locations:', total_brine_rate, sep='\n')

    # Plot scalar observations
    label_size = 13
    title_size = 16
    f1, ax = plt.subplots(1, 2, figsize=(12, 5))
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.9, wspace=0.25)
    ax[0].semilogy(time_array/365.25, total_CO2_rate)
    ax[0].set_xlabel('Time, [years]', fontsize=label_size)
    ax[0].set_ylabel('Leakage rates, [kg/s]', fontsize=label_size)
    ax[0].set_title(r'Total leakage rates of CO$_2$', fontsize=title_size)

    ax[1].semilogy(time_array/365.25, total_brine_rate)
    ax[1].set_xlabel('Time, [years]', fontsize=label_size)
    ax[1].set_ylabel('Leakage rates, [kg/s]', fontsize=label_size)
    ax[1].set_title('Total leakage rates of brine', fontsize=title_size)

    # Collect gridded observations
    grid_outputs = {}
    # From reservoir component
    for nm in ['pressure', 'CO2saturation']:
        grid_outputs[nm] = []
        for ind in range(num_years+1):
            filename = 'ltres_{}_sim_0_time_{}.{}'.format(nm, ind, grid_save_type)
            file_loc = os.sep.join([output_dir, filename])
            d = sm.get_gridded_observation_file(file_loc)
            grid_outputs[nm].append(d)
        grid_outputs[nm] = np.array(grid_outputs[nm])

    # From SealHorizon component
    for nm in ['CO2_aquifer', 'brine_aquifer']:
        grid_outputs[nm] = []
        for ind in range(num_years+1):
            filename = 'shc_{}_sim_0_time_{}.npz'.format(nm, ind)
            file_loc = os.sep.join([output_dir, filename])
            d = sm.get_gridded_observation_file(file_loc)
            grid_outputs[nm].append(d)
        grid_outputs[nm] = np.array(grid_outputs[nm])

    # Plot gridded observations
    line_width = 1
    f2, ax = plt.subplots(2, 2, figsize=(12, 8))
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95, wspace=0.25)
    for ind in range(num_cells):
        # Plot pressure
        ax[0, 0].plot(time_array/365.25, grid_outputs['pressure'],
                      '-r', linewidth=line_width)
        ax[0, 0].set_xlabel('Time, [years]', fontsize=label_size)
        ax[0, 0].set_ylabel('Pressure, [Pa]', fontsize=label_size)
        ax[0, 0].get_yaxis().set_label_coords(-0.16, 0.5)

        # Plot saturation
        ax[0, 1].plot(time_array/365.25, grid_outputs['CO2saturation'],
                      '-b', linewidth=line_width)
        ax[0, 1].set_xlabel('Time, [years]', fontsize=label_size)
        ax[0, 1].set_ylabel(r'CO$_2$ saturation, [-]', fontsize=label_size)
        ax[0, 1].get_yaxis().set_label_coords(-0.16, 0.5)

        # Plot CO2 leakage rates
        ax[1, 0].plot(time_array/365.25, grid_outputs['CO2_aquifer'],
                      '-g', linewidth=line_width)
        ax[1, 0].set_xlabel('Time, [years]', fontsize=label_size)
        ax[1, 0].set_ylabel(r'CO$_2$ leakage rates, [kg/s]', fontsize=label_size)
        ax[1, 0].get_yaxis().set_label_coords(-0.16, 0.5)

        # Plot brine leakage rates
        ax[1, 1].plot(time_array/365.25, grid_outputs['brine_aquifer'],
                      '-k', linewidth=line_width)
        ax[1, 1].set_xlabel('Time, [years]', fontsize=label_size)
        ax[1, 1].set_ylabel('Brine leakage rates, [kg/s]', fontsize=label_size)
        ax[1, 1].get_yaxis().set_label_coords(-0.16, 0.5)
