"""
Example illustrates linking of lookup table reservoir, Kimberlina wellbore and
deep alluvium aquifer (based on ML) components.

This example requires the additional Kimberlina data set.
Kimberlina data set (kimb-closed-200-sims.zip or Kimb_closed_200_sims.zip)
can be downloaded from one of the following places:
1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam
2. https://gitlab.com/NRAP/Kimberlina_data

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/Kimb_closed_200_sims

Example of run:
$ python iam_sys_lutreservoir_Kimbwell_deepalluvium_aquifer.py
"""
import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, ReservoirDataInterpolator, LookupTableReservoir,
                     KimberlinaWellbore, DeepAlluviumAquiferML)


if __name__ == '__main__':
    __spec__ = None  # for multiprocessing in Spyder (Windows)
    logging.basicConfig(level=logging.ERROR)

    # Define keyword arguments of the system model
    num_years = 10
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    data_folder = os.path.join('..', '..', 'source', 'components', 'reservoir',
                               'lookuptables', 'Kimb_closed_200_sims')
    if not os.path.exists(os.sep.join([data_folder, 'Reservoir_data_sim0001.csv'])):
        msg = ''.join(['\nKimberlina data set can be downloaded ',
                       'from one of the following places:\n',
                       '1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam \n',
                       '2. https://gitlab.com/NRAP/Kimberlina_data \n',
                       'Check this example description for more information.'])
        logging.error(msg)

    # Read file with signatures of interpolators and names of files with the corresponding data
    signature_data = np.genfromtxt(
        os.path.join(data_folder, 'parameters_and_filenames.csv'),
        delimiter=",", dtype='str')

    # The first row (except the last element) of the file contains names of the parameters
    # First column (column=0) in the signature_data is for index so the parameter
    # names from the second column (column=1)
    par_names = signature_data[0, 1:-1]
    num_pars = len(par_names)

    # Create and add interpolators to the system model
    # corresponding to the first 4 lookup tables
    num_luts = 4
    for ind in range(num_luts):
        # Create signature for the interpolator to be created
        signature = {par_names[j]: float(signature_data[ind+1, j+1]) for j in range(num_pars)}

        intpr = sm.add_interpolator(ReservoirDataInterpolator(
            name='intpr{}'.format(ind+1), parent=sm, header_file_dir=data_folder,
            time_file='time_points.csv',
            data_file='Reservoir_data_sim{ind1:04}.csv'.format(ind1=ind+1),
            index=int(signature_data[ind+1, 0]),
            signature=signature), intr_family='reservoir')

        debug_msg = 'Signature of the created interpolator is {}'.format(signature)
        logging.debug(debug_msg)

    logging.debug('All interpolators are created')

    # Coordinates of the injection well
    inj_well_x = 300025.19
    inj_well_y = 3934175.82

    # Leaking wells at the distance of 0.2 km (200 m)
    num_wells = 5
    theta = np.linspace(0, 2*np.pi, num=num_wells, endpoint=False)
    dist_radius = 200
    leak_well_x = inj_well_x + dist_radius*np.cos(theta)
    leak_well_y = inj_well_y + dist_radius*np.sin(theta)

    # Plot grid points, injection and leaking wells
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 6))
    ax.plot(intpr.points[:, 0], intpr.points[:, 1], 'ok',
            label='grid points', markersize=4)
    ax.plot(inj_well_x, inj_well_y, '*r', label='injection well')
    ax.plot(leak_well_x, leak_well_y, 'sb', label='leaking well')
    for ind in range(num_wells):
        ax.annotate(str(ind+1), (leak_well_x[ind], leak_well_y[ind]),
                    color='red', fontsize=14)
    ax.set_xlabel('x [km]')
    ax.set_ylabel('y [km]')
    ax.set_aspect('equal')
    ax.set_xlim(298500, 301500)
    ax.set_ylim(3932500, 3935500)
    ax.set_xticklabels(np.linspace(298.5, 301.5, num=7, endpoint=True))
    ax.set_yticklabels(np.linspace(3932.5, 3935.5, num=7, endpoint=True))
    ax.legend()

    # Initialize list of reservoir components
    ltress = []
    # List of observations of reservoir component to be added to the system model
    res_obs = ['pressure', 'CO2saturation']

    # Loop over all wells
    for ind in range(num_wells):
        # Add reservoir component
        ltress.append(sm.add_component_model_object(LookupTableReservoir(
            name='ltres{}'.format(ind+1), parent=sm, intr_family='reservoir',
            locX=leak_well_x[ind], locY=leak_well_y[ind])))

        # Add index parameter of reservoir component model
        if ind == 0:
            ltress[-1].add_par('index', value=1,
                               discrete_vals=(list(range(1, num_luts+1)),
                                              num_luts*[1/num_luts]))
        else:
            # All reservoir interpolators need to have the same signature index
            # during simulation
            ltress[-1].add_par_linked_to_par(
                'index', ltress[0].pars['index'])

        for obs_nm in res_obs:
            # Add observations of reservoir component model
            ltress[-1].add_obs(obs_nm)
            # Add observations to be used as input for Kimberlina wellbore component
            ltress[-1].add_obs_to_be_linked(obs_nm)

    # Initialize list of wellbore components
    kws = []
    # List of observations of wellbore component to be added to the system model
    well_obs = ['CO2_mass', 'brine_mass', 'CO2_rate', 'brine_rate']

    # Loop over all wells
    for ind in range(num_wells):
        # Add Kimberlina wellbore component
        kws.append(sm.add_component_model_object(KimberlinaWellbore(
            name='kw{}'.format(ind+1), parent=sm)))

        # Assume that some CO2 and brine are already present in the system
        # to avoid warning messages about zero masses
        kws[-1].model_kwargs['previous_co2_mass'] = 0.01
        kws[-1].model_kwargs['previous_brine_mass'] = 50.0

        # Add parameters of Kimberlina wellbore component
        kws[-1].add_par('logK_well_permeability',
                        min=-11.5, max=-10.1, value=-11.1 + ind/10.0)

        # Add keyword arguments linked to the output provided by reservoir model
        for obs_nm in res_obs:
            kws[-1].add_kwarg_linked_to_obs(obs_nm, ltress[ind].linkobs[obs_nm])

        for obs_nm in well_obs:
            # Add observations of Kimberlina wellbore component
            kws[-1].add_obs(obs_nm)
            # Add observations to be used as input for aquifer components
            kws[-1].add_obs_to_be_linked(obs_nm)

    # Initialize list of aquifer components
    daas = []
    # List of observations of aquifer component to be added to the system model
    aq_obs = ['TDS_volume', 'TDS_dx', 'TDS_dy', 'TDS_dz',
              'Pressure_volume', 'Pressure_dx', 'Pressure_dy', 'Pressure_dz',
              'pH_volume', 'pH_dx', 'pH_dy', 'pH_dz']

    # Loop over all wells
    for ind in range(num_wells):
        # Add deep alluvium aquifer component
        daas.append(sm.add_component_model_object(DeepAlluviumAquiferML(
            name='daa{}'.format(ind+1), parent=sm)))

        # Add parameters of the aquifer component
        daas[-1].add_par('logK_sand1', value=-11.3, vary=False)
        daas[-1].add_par('logK_sand2', value=-11.3, vary=False)
        daas[-1].add_par('logK_sand3', value=-11.3, vary=False)
        daas[-1].add_par('logK_caprock', value=-16.5, vary=False)
        daas[-1].add_par('correlationLengthX', value=520.721, vary=False)
        daas[-1].add_par('correlationLengthZ', value=112.442, vary=False)
        daas[-1].add_par('sandFraction', value=0.743, vary=False)
        daas[-1].add_par('groundwater_gradient', value=0.00105408, vary=False)
        daas[-1].add_par('leak_depth', value=715.99, vary=False)

        # Add keyword arguments linked to the output provided by wellbore component
        for obs_nm in well_obs:
            daas[-1].add_kwarg_linked_to_obs(obs_nm.lower(), kws[ind].linkobs[obs_nm])

        # Add observations (output) of the aquifer component
        for obs_nm in aq_obs:
            daas[-1].add_obs(obs_nm)

    # Run system model using current values of its parameters
    print('--------------------------------------------------------------')
    print('Started forward simulation...')
    print('--------------------------------------------------------------')
    sm.forward()  # system model is run deterministically

    # Setup plots parameters
    line_width = 1
    to_print_results = True
    # ------------------------------------------------------------------------
    print('--------------------------------------------------------------')
    print('-------------- Reservoir component observations --------------')
    # Create figure
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    collected_obs = {}
    for ind in range(num_wells):
        # Collect observations
        for obs_nm in res_obs:
            collected_obs[obs_nm] = sm.collect_observations_as_time_series(
                ltress[ind], obs_nm)
        if to_print_results:
            print('--------------------------------------------------------------')
            print('Pressure at well {}'.format(ind+1),
                  collected_obs['pressure'], sep='\n')
            print('CO2 saturation at well {}'.format(ind+1),
                  collected_obs['CO2saturation'], sep='\n')

        # Plot data
        ax[0].plot(time_array/365.25, collected_obs['pressure']/1.0e+6, '-',
                   linewidth=line_width, label='well {}'.format(ind+1))
        ax[1].plot(time_array/365.25, collected_obs['CO2saturation'], '-',
                   linewidth=line_width, label='well {}'.format(ind+1))

    res_y_labels = ['Pressure, P [MPa]', r'CO$_2$ saturation, S [-]']
    # Add plot labels and legend
    for ind in range(2):
        ax[ind].set_xlabel('Time, t [years]')
        ax[ind].set_ylabel(res_y_labels[ind])
        ax[ind].legend()
    fig.subplots_adjust(left=0.1, top=0.9, right=0.95, bottom=0.15, wspace=0.25)
    fig.suptitle('Reservoir component observations: single realization')

    # ------------------------------------------------------------------------
    print('--------------------------------------------------------------')
    print('-------------- Wellbore component observations ---------------')
    # Define y-axis labels
    well_y_labels = [r'CO$_2$ mass, m [kg]', 'Brine mass, m [kg]',
                     r'CO$_2$ rate, q [kg/s]', 'Brine rate, q [kg/s]']
    # Create figure
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    collected_obs = {}
    for ind in range(num_wells):
        # Collect observations
        for obs_nm in well_obs:
            collected_obs[obs_nm] = sm.collect_observations_as_time_series(
                kws[ind], obs_nm)
        if to_print_results:
            print('--------------------------------------------------------------')
            print('CO2 mass through well {}'.format(ind+1),
                  collected_obs['CO2_mass'], sep='\n')
            print('Brine mass through well {}'.format(ind+1),
                  collected_obs['brine_mass'], sep='\n')
            print('CO2 rate through well {}'.format(ind+1),
                  collected_obs['CO2_rate'], sep='\n')
            print('Brine rate through well {}'.format(ind+1),
                  collected_obs['brine_rate'], sep='\n')

        # Plot data
        for nm_ind, obs_nm in enumerate(well_obs):
            ax[nm_ind//2, nm_ind%2].plot(time_array/365.25, collected_obs[obs_nm],
                                         '-', linewidth=line_width,
                                         label='well {}'.format(ind+1))
            # Add y-axis label and adjust its position
            ax[nm_ind//2, nm_ind%2].set_ylabel(well_y_labels[nm_ind])
            ax[nm_ind//2, nm_ind%2].get_yaxis().set_label_coords(-0.2, 0.5)

    # Add x-axis label
    for ind in range(2):
        ax[1, ind].set_xlabel('Time, t [years]')

    # Add legend
    ax[1, 1].legend()
    # Adjust space between subplots
    fig.subplots_adjust(left=0.1, top=0.95, right=0.95, bottom=0.08, wspace=0.25)
    fig.suptitle('Wellbore component observations: single realization')

    # ------------------------------------------------------------------------
    print('--------------------------------------------------------------')
    print('--------------- Aquifer component observations ---------------')
    # Define y-axis labels
    aq_y_labels = ['TDS volume, V [m$^3$]', 'TDS dx, [m]', 'TDS dy, [m]', 'TDS dz, [m]',
                   'Pressure volume, V [m$^3$]', 'Pressure dx, [m]',
                   'Pressure dy, [m]', 'Pressure dz, [m]',
                   'pH volume, V [m$^3$]', 'pH dx, [m]', 'pH dy, [m]', 'pH dz, [m]']
    # Create figure
    fig, ax = plt.subplots(3, 4, figsize=(18, 9), sharex=True)
    collected_obs = {}
    for ind in range(num_wells):
        # Collect observations
        for obs_nm in aq_obs:
            collected_obs[obs_nm] = sm.collect_observations_as_time_series(
                daas[ind], obs_nm)
        if to_print_results:
            print('--------------------------------------------------------------')
            print('TDS volume affected because of well {}'.format(ind+1),
                  collected_obs['TDS_volume'], sep='\n')
            print('TDS dx affected because of well {}'.format(ind+1),
                  collected_obs['TDS_dx'], sep='\n')
            print('TDS dy affected because of well {}'.format(ind+1),
                  collected_obs['TDS_dy'], sep='\n')
            print('TDS dz affected because of well {}'.format(ind+1),
                  collected_obs['TDS_dz'], sep='\n')
            print('--------------------------------------------------------------')
            print('Pressure volume affected because of well {}'.format(ind+1),
                  collected_obs['Pressure_volume'], sep='\n')
            print('Pressure dx affected because of well {}'.format(ind+1),
                  collected_obs['Pressure_dx'], sep='\n')
            print('Pressure dy affected because of well {}'.format(ind+1),
                  collected_obs['Pressure_dy'], sep='\n')
            print('Pressure dz affected because of well {}'.format(ind+1),
                  collected_obs['Pressure_dz'], sep='\n')
            print('--------------------------------------------------------------')
            print('pH volume affected because of well {}'.format(ind+1),
                  collected_obs['pH_volume'], sep='\n')
            print('pH dx affected because of well {}'.format(ind+1),
                  collected_obs['pH_dx'], sep='\n')
            print('pH dy affected because of well {}'.format(ind+1),
                  collected_obs['pH_dy'], sep='\n')
            print('pH dz affected because of well {}'.format(ind+1),
                  collected_obs['pH_dz'], sep='\n')

        # Plot data
        for nm_ind, obs_nm in enumerate(aq_obs):
            ax[nm_ind//4, nm_ind%4].plot(time_array/365.25, collected_obs[obs_nm],
                                         '-', linewidth=line_width,
                                         label='well {}'.format(ind+1))
            ax[nm_ind//4, nm_ind%4].set_ylabel(aq_y_labels[nm_ind])
            if nm_ind%4 == 0:
                ax[nm_ind//4, nm_ind%4].get_yaxis().set_label_coords(-0.2, 0.5)
            else:
                ax[nm_ind//4, nm_ind%4].get_yaxis().set_label_coords(-0.14, 0.5)

    # Add x-axis label
    for ind in range(4):
        ax[2, ind].set_xlabel('Time, t [years]')

    # Add legend
    ax[2, 3].legend()
    # Adjust space between subplots
    fig.subplots_adjust(left=0.1, top=0.95, right=0.95, bottom=0.08, wspace=0.25)
    fig.suptitle('Aquifer component observations: single realization')

    # ------------------------------------------------------------------------
    # Setup LHS based study
    num_samples = 24
    ncpus = 4
    s = sm.lhs(siz=num_samples, seed=765432)

    # Run model using values in samples for parameter values
    print('--------------------------------------------------------------')
    print('Started LHS based simulations...')
    print('--------------------------------------------------------------')
    s.run(cpus=ncpus, verbose=False)

    # Extract results from stochastic simulations
    out = s.collect_observations_as_time_series()

    for ind in range(num_wells):
        collected_obs = {}
        fig, ax = plt.subplots(3, 4, figsize=(18, 9), sharex=True)
        for nm_ind, nm in enumerate(aq_obs):
            collected_obs[nm] = out['daa{}.{}'.format(ind+1, nm)]
            for j in range(num_samples):
                ax[nm_ind//4, nm_ind%4].plot(time_array/365.25, collected_obs[nm][j],
                                             '-', color='grey', linewidth=line_width)
            ax[nm_ind//4, nm_ind%4].set_ylabel(aq_y_labels[nm_ind])
            if nm_ind%4 == 0:
                ax[nm_ind//4, nm_ind%4].get_yaxis().set_label_coords(-0.2, 0.5)
            else:
                ax[nm_ind//4, nm_ind%4].get_yaxis().set_label_coords(-0.14, 0.5)

        # Add x-axis label
        for j in range(4):
            ax[2, j].set_xlabel('Time, t [years]')

        # Adjust space between subplots
        fig.subplots_adjust(left=0.1, top=0.95, right=0.95, bottom=0.08, wspace=0.25)
        fig.suptitle('Aquifer impact results associated with well {}'.format(ind+1))

    plt.show()

    # Remove all handlers from the logger for proper work in the consecutive runs
    while logging.getLogger('').handlers:
        logging.getLogger('').handlers.pop()
