import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, ReservoirDataInterpolator,
                     LookupTableReservoir, SealHorizon,
                     SHThicknessSampler, SHPermeabilitySampler)

if __name__ == "__main__":

    # Setup location of lookup table data set
    file_directory = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                  'lookuptables', 'Kimb_54_sims'])

    if not os.path.exists(os.sep.join([file_directory, 'Reservoir_data_sim01.csv'])):
        msg = ''.join(['\nKimberlina data set can be downloaded ',
                       'from one of the following places:\n',
                       '1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam \n',
                       '2. https://gitlab.com/NRAP/Kimberlina_data  \n',
                       'Check this example description for more information.'])
        logging.error(msg)

    # Define output directory
    output_directory = os.sep.join(['..', '..', 'output', 'scripts'])

    # Define file type for saving gridded observations
    grid_save_type = 'npz'

    # Create output directory if needed
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Define domain of the problem: we focus on area with positive CO2 saturation
    xa = 38500.0
    xb = 41500.0
    yc = 47500.0
    yd = 50500.0
    nx = 21
    ny = 21

    loc_x = np.linspace(xa, xb, num=nx, endpoint=True)
    loc_y = np.linspace(yc, yd, num=ny, endpoint=True)

    loc_xx, loc_yy = np.meshgrid(loc_x, loc_y)
    loc_xx = loc_xx.T.reshape((nx*ny, ))
    loc_yy = loc_yy.T.reshape((nx*ny, ))

    cell_area = (loc_x[1]-loc_x[0])*(loc_y[1]-loc_y[0])
    index_value = 3

    # Define keyword arguments of the system model.
    num_years = 10  # reduce simulation time
    time_array = 365.25 * np.arange(0.0, num_years+1)
    seal_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model.
    sm = SystemModel(model_kwargs=seal_model_kwargs)

    # Add reservoir component.
    interpr1 = sm.add_interpolator(ReservoirDataInterpolator(
        name='int1', parent=sm,
        header_file_dir=file_directory,
        time_file='time_points.csv',
        data_file='Reservoir_data_sim0{}.csv'.format(index_value), # for index value < 10
        interp_2d=True,
        index=index_value,
        signature={}), intr_family='reservoir1')

    ltres = sm.add_component_model_object(
        LookupTableReservoir(name='ltres', parent=sm,
                             intr_family='reservoir1',
                             parameter_filename='parameters_and_filenames.csv',
                             locX=loc_xx, locY=loc_yy))
    # Add parameters and observations
    ltres.add_par('index', value=index_value, vary=False)
    ltres_obs = ['pressure', 'CO2saturation']
    for obs_nm in ltres_obs:
        ltres.add_grid_obs(obs_nm, constr_type='array',
                           output_dir=output_directory, save_type=grid_save_type)
        # Add observations of reservoir component model to be linked
        ltres.add_obs_to_be_linked(obs_nm, obs_type='grid', constr_type='array')

    # Calculate pressure and saturation at initial time point
    res_data = ltres.model({'index': index_value}, time_point=0)
    cell_init_pressure = res_data['pressure']
    brine_density = 1080.0
    # Approximate depth of cell bases
    cell_depth = (cell_init_pressure-101325.0)/(9.8*brine_density)
    ave_depth = np.mean(cell_depth)
    ave_pressure = np.mean(cell_init_pressure)

    # Add samplers components
    # Setup thickness sampler
    t_sampler = sm.add_component_model_object(
        SHThicknessSampler(name='ts', parent=sm, grid_shape=(nx*ny,)))
    t_sampler.add_par('seed', value=111, vary=False)  # arbitrary value
    # The chosen samples are just for illustration purposes
    t_sampler.add_par('aveThickness', value=450, vary=False)
    t_sampler.add_par('stdDevThickness', value=12, vary=False)
    t_sampler.add_par('minThickness', value=430, vary=False)
    max_thickness = 460.0
    t_sampler.add_par('maxThickness', value=max_thickness, vary=False)
    t_sampler.add_grid_obs('thickness', constr_type='array', index=[0],
                           output_dir=output_directory, save_type=grid_save_type)

    # Setup permeability sampler
    perm_sampler = sm.add_component_model_object(
        SHPermeabilitySampler(name='ps', parent=sm, grid_shape=(nx*ny,)))
    perm_sampler.add_par('seed', value=123, vary=False)
    perm_sampler.add_par('avePermeability', value=1.0e-17, vary=False)
    perm_sampler.add_par('stdDevPermeability', value=1.0e-18, vary=False)
    perm_sampler.add_par('minPermeability', value=0.7e-17, vary=False)
    perm_sampler.add_par('maxPermeability', value=1.0e-15, vary=False)
    perm_sampler.add_grid_obs('permeability', constr_type='array', index=[0],
                              output_dir=output_directory, save_type=grid_save_type)

    # Setup Seal Horizon component
    shc = sm.add_component_model_object(
        SealHorizon(name='shc', parent=sm,
                    locX=loc_xx, locY=loc_yy, area=cell_area))
    # Setup some parameters
    shc.add_par('brineDensity', value=brine_density, vary=False)
    shc.add_par('CO2Density', value=590.0, vary=False)

    shc.add_par('aveBaseDepth', value=ave_depth, vary=False)
    shc.add_par('aveBasePressure', value=ave_pressure, vary=False)
    shc.add_par('staticDepth', value=ave_depth-max_thickness, vary=False)
    stat_pressure = 101325.0+9.8*brine_density*(ave_depth-460.0)
    shc.add_par('staticPressure', value=stat_pressure, vary=False)

    shc.add_par('influenceModel', value=0, vary=False)
    shc.model_kwargs['relativeModel'] = 'BC'
    shc.add_par('lambda', value=2.0, vary=False)

    # Link reservoir component outputs
    shc.add_kwarg_linked_to_obs(
        'pressure', ltres.linkobs['pressure'],
        obs_type='grid', constr_type='array')
    shc.add_kwarg_linked_to_obs(
        'CO2saturation', ltres.linkobs['CO2saturation'],
        obs_type='grid', constr_type='array')

    # Link thickness sampler output
    shc.add_kwarg_linked_to_obs(
        'thickness', t_sampler.linkobs['thickness'],
        obs_type='grid', constr_type='array')

    # Link permeability sampler output
    shc.add_kwarg_linked_to_obs(
        'permeability', perm_sampler.linkobs['permeability'],
        obs_type='grid', constr_type='array')

    # Define SH observations
    shc_scalar_obs = ['CO2_aquifer_total', 'brine_aquifer_total',
                      'mass_CO2_aquifer_total', 'mass_brine_aquifer_total']
    shc_grid_obs = ['CO2_aquifer', 'brine_aquifer',
                    'mass_CO2_aquifer', 'mass_brine_aquifer']

    # Add scalar observations
    for obs_nm in shc_scalar_obs:
        shc.add_obs(obs_nm)

    # Add gridded observations
    for obs_nm in shc_grid_obs:
        shc.add_grid_obs(obs_nm, constr_type='array',
                         output_dir=output_directory, save_type=grid_save_type)

    print('Starting simulation...')
    sm.forward()
    print('Simulation is finished.')

    # Collect gridded observations
    outputs = {}
    for obs_nm in ltres_obs:
        outputs[obs_nm] = sm.collect_gridded_observations_as_time_series(
            ltres, obs_nm, output_directory, rlzn_number=0, save_type=grid_save_type)

    for obs_nm in shc_grid_obs:
        outputs[obs_nm] = sm.collect_gridded_observations_as_time_series(
            shc, obs_nm, output_directory, rlzn_number=0, save_type=grid_save_type)

    outputs['thickness'] = sm.collect_gridded_observations_as_time_series(
        t_sampler, 'thickness', output_directory, indices=[0], rlzn_number=0,
        save_type=grid_save_type)

    outputs['permeability'] = sm.collect_gridded_observations_as_time_series(
        perm_sampler, 'permeability', output_directory, indices=[0], rlzn_number=0,
        save_type=grid_save_type)

    # Collect scalar observations
    for obs_nm in shc_scalar_obs:
        outputs[obs_nm] = sm.collect_observations_as_time_series(shc, obs_nm)

    # Plot pressure profiles
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 7))
    for ind in range(nx*ny): # go over all cells
        ax.plot(time_array/365.25, outputs['pressure'][:, ind]/1.0e+6)
        ax.set_xlabel('Time, [years]')
        ax.set_ylabel('Pressure, [MPa]')
    fig.savefig(os.sep.join([
        output_directory, 'samplers_script_pressure.png']), dpi=200)

    # Plot CO2 saturation profiles
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 7))
    for ind in range(nx*ny): # go over all cells
        ax.plot(time_array/365.25, outputs['CO2saturation'][:, ind])
        ax.set_xlabel('Time, [years]')
        ax.set_ylabel(r'CO$_2$ saturation, [-]')
    fig.savefig(os.sep.join([
        output_directory, 'samplers_script_CO2saturation.png']), dpi=200)

    # Plot CO2 leakage rates profiles
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 7))
    for ind in range(nx*ny): # go over all cells
        ax.plot(time_array/365.25, outputs['CO2_aquifer'][:, ind])
        ax.set_xlabel('Time, [years]')
        ax.set_ylabel(r'CO$_2$ leakage rates, [kg/s]')
    fig.savefig(os.sep.join([
        output_directory, 'samplers_script_CO2_leakage_rates.png']), dpi=200)

    # Plot brine leakage rates profiles
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 7))
    for ind in range(nx*ny): # go over all cells
        ax.plot(time_array/365.25, outputs['brine_aquifer'][:, ind])
        ax.set_xlabel('Time, [years]')
        ax.set_ylabel('Brine leakage rates, [kg/s]')
    fig.savefig(os.sep.join([
        output_directory, 'samplers_script_brine_leakage_rates.png']), dpi=200)

    # Plot CO2 leaked masses profiles
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 7))
    for ind in range(nx*ny): # go over all cells
        ax.plot(time_array/365.25, outputs['mass_CO2_aquifer'][:, ind])
        ax.set_xlabel('Time, [years]')
        ax.set_ylabel(r'CO$_2$ leaked masses, [kg]')
    fig.savefig(os.sep.join([
        output_directory, 'samplers_script_CO2_leaked_mass.png']), dpi=200)

    # Plot brine leaked masses profiles
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 7))
    for ind in range(nx*ny): # go over all cells
        ax.plot(time_array/365.25, outputs['mass_brine_aquifer'][:, ind])
        ax.set_xlabel('Time, [years]')
        ax.set_ylabel(r'Brine leaked masses, [kg]')
    fig.savefig(os.sep.join([
        output_directory, 'samplers_script_brine_leaked_mass.png']), dpi=200)

    # Plot total leakage rates
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 7))
    ax[0].plot(time_array/365.25, outputs['CO2_aquifer_total'])
    ax[1].plot(time_array/365.25, outputs['brine_aquifer_total'])
    for j in range(2):
        ax[j].set_xlabel('Time, [years]')
    ax[0].set_ylabel(r'Total CO$_2$ leakage rates, [kg/s]')
    ax[1].set_ylabel(r'Total brine leakage rates, [kg/s]')
    fig.savefig(os.sep.join([
        output_directory, 'samplers_script_total_leakage_rates.png']), dpi=200)

    # Plot total leaked masses
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 7))
    ax[0].plot(time_array/365.25, outputs['mass_CO2_aquifer_total'])
    ax[1].plot(time_array/365.25, outputs['mass_brine_aquifer_total'])
    for j in range(2):
        ax[j].set_xlabel('Time, [years]')
    ax[0].set_ylabel(r'Total CO$_2$ leaked mass, [kg]')
    ax[1].set_ylabel(r'Total brine leaked mass, [kg]')
    fig.savefig(os.sep.join([
        output_directory, 'samplers_script_total_leaked_masses.png']), dpi=200)
    plt.close('all')

    # Plot generated thickness
    x = np.arange(1, nx+1)
    y = np.arange(1, ny+1)
    xx, yy = np.meshgrid(x, y)
    xx = xx.T.reshape(nx*ny,)
    yy = yy.T.reshape(nx*ny,)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 7))
    im = ax.scatter(xx, yy, c=outputs['thickness'][0], marker='s', s=200)
    ax.set_title('Seal layer thickness')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    fig.colorbar(im, ax=ax)
    fig.savefig(os.sep.join([
        output_directory, 'samplers_script_thickness.png']), dpi=200)

    # Plot generated permeability
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 7))
    im = ax.scatter(xx, yy, c=outputs['permeability'][0], marker='s', s=200)
    ax.set_title('Seal layer permeability')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    fig.colorbar(im, ax=ax)
    fig.savefig(os.sep.join([
        output_directory, 'samplers_script_permeability.png']), dpi=200)
    plt.close('all')
