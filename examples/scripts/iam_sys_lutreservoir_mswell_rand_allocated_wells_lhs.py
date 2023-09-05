"""
Example illustrates system model with multiple lookup table reservoir components
each linked to a separate multisegmented wellbore component. Additionally, example
illustrates application of LocationGenerator component for a random placement
of the predetermined number of wells within a defined spatial domain.
The original domain is divided into 4 subdomains: one well is placed randomly
into each subdomain.

This example requires the additional Kimberlina data set.
Kimberlina data set (pressure-in-pa-kimb-54-sims.zip or
Pressure_In_Pa_Kimb_54_sims.zip) can be downloaded from one of the following places:
1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam
2. https://gitlab.com/NRAP/Kimberlina_data

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/Kimb_54_sims
"""
import sys
import os
import logging
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, ReservoirDataInterpolator,
                     LookupTableReservoir, MultisegmentedWellbore,
                     LocationGenerator)


if __name__ == '__main__':
    logging.basicConfig(level=logging.ERROR)
    __spec__ = None

    file_directory = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                  'lookuptables', 'Kimb_54_sims'])

    if not os.path.exists(os.sep.join([file_directory, 'Reservoir_data_sim01.csv'])):
        msg = ''.join(['\nKimberlina data set can be downloaded ',
                       'from one of the following places:\n',
                       '1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam \n',
                       '2. https://gitlab.com/NRAP/Kimberlina_data  \n',
                       'Check this example description for more information.'])
        logging.error(msg)

    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Read file with signatures of interpolators and names of files with the corresponding data
    signature_data = np.genfromtxt(
        os.path.join(file_directory, 'parameters_and_filenames.csv'),
        delimiter=",", dtype='str')

    # The first row (except the last element) of the file contains names of the parameters
    par_names = signature_data[0, 1:-1]
    num_pars = len(par_names)

    # Create only two interpolators
    num_interpolators = 2

    # Create and add interpolators to the system model
    for ind in range(num_interpolators):
        signature = {par_names[j]: float(signature_data[ind+1, j+1]) for j in range(num_pars)}

        intpr = sm.add_interpolator(
            ReservoirDataInterpolator(
                name='int'+str(ind+1), parent=sm,
                header_file_dir=file_directory,
                time_file='time_points.csv',
                data_file='Reservoir_data_sim{ind1:02}.csv'.format(ind1=ind+1),
                index=int(signature_data[ind+1, 0]),
                signature=signature),
            intr_family='reservoir')

        debug_msg = 'Signature of the created interpolator is {}'.format(signature)
        logging.debug(debug_msg)

    logging.debug('All needed interpolators are created')

    # Define number of wells
    num_wells = 4
    # Define boundaries of random wells domain
    # We split domain into number of subdomains (equal to num_wells)
    # of the same area and randomly place one well within each subdomain
    x_min = 35800.0
    x_max = 40700.0
    y_min = 47800.0
    y_max = 49800.0
    x_coords = np.linspace(x_min, x_max, num=3, endpoint=True)
    y_coords = np.linspace(y_min, y_max, num=3, endpoint=True)

    # Add generator components
    gens = []
    # Set generator index
    w_ind = 1
    for ind1 in range(2):
        for ind2 in range(2):
            gens.append(sm.add_component_model_object(LocationGenerator(
                name='gen{}'.format(w_ind), parent=sm,
                x_min=x_coords[ind1], x_max=x_coords[ind1+1],
                y_min=y_coords[ind2], y_max=y_coords[ind2+1],
                num_locations=1, reproducible=True)))

            # Add parameters of generator component
            # We use generator index to set the initial seed value
            gens[-1].add_par('seed', value=2*w_ind+10, min=3*w_ind, max=900*w_ind)
            gens[-1].add_obs_to_be_linked('locX', obs_type='grid')
            gens[-1].add_obs_to_be_linked('locY', obs_type='grid')

            # Update index for next generator
            w_ind = w_ind+1

    # Initialize list of reservoir components
    ltress = []

    # Add LUT components: one for each well
    for ind in range(num_wells):
        ltress.append(sm.add_component_model_object(
            LookupTableReservoir(name='ltres{}'.format(ind+1), parent=sm,
                                 intr_family='reservoir')))

        # Link locations of the LUT components to the output of generator
        ltress[ind].add_kwarg_linked_to_obs(
            'locX', gens[ind].linkobs['locX'], obs_type='grid', constr_type='array', loc_ind=[0])
        ltress[ind].add_kwarg_linked_to_obs(
            'locY', gens[ind].linkobs['locY'], obs_type='grid', constr_type='array', loc_ind=[0])

        # Add observations of reservoir component model
        ltress[ind].add_obs('pressure')
        ltress[ind].add_obs('CO2saturation')
        # Add observations of reservoir component model to be use for the input
        # of multisegmented wellbore components
        ltress[ind].add_obs_to_be_linked('pressure')
        ltress[ind].add_obs_to_be_linked('CO2saturation')

        # Add signature index parameter of LUT reservoir components
        if ind == 0:
            ltress[0].add_par('index', value=2, vary=False)
        else:
            # All reservoir components will correspond to the same lookup table
            ltress[ind].add_par_linked_to_par(
                'index', ltress[0].deterministic_pars['index'])

    # Initialize list of wellbore components
    mss = []
    for i in range(num_wells):
        mss.append(sm.add_component_model_object(
            MultisegmentedWellbore(name='ms{}'.format(i+1), parent=sm)))

        # Add parameters of multisegmented wellbore components
        mss[-1].add_par('wellRadius', value=0.015, vary=False)
        mss[-1].add_par('numberOfShaleLayers', value=4, vary=False)
        mss[-1].add_par('shale1Thickness', value=220., vary=False)
        mss[-1].add_par('shale2Thickness', value=500., vary=False)
        mss[-1].add_par('shale3Thickness', value=350., vary=False)
        mss[-1].add_par('shale4Thickness', value=340., vary=False)
        mss[-1].add_par('aquifer1Thickness', value=150., vary=False)
        mss[-1].add_par('aquifer2Thickness', value=720., vary=False)
        mss[-1].add_par('aquifer3Thickness', value=400., vary=False)
        mss[-1].add_par('reservoirThickness', value=400., vary=False)
        mss[-1].add_par('logWellPerm', value=-12.5, vary=False)

        # Add observations of multisegmented wellbore components
        mss[-1].add_obs('CO2_aquifer1')
        mss[-1].add_obs('brine_aquifer1')

        # Add keyword arguments linked to the output provided by reservoir model
        mss[-1].add_kwarg_linked_to_obs('pressure', ltress[i].linkobs['pressure'])
        mss[-1].add_kwarg_linked_to_obs('CO2saturation', ltress[i].linkobs['CO2saturation'])

    sm.forward()

    print('__________________Results of forward simulation___________________')
    print('____________________ Pressure and saturation _____________________')
    print('')

    # Print pressure and saturation
    linespec = ['-k', '-r', '-b', '-g', '-y', '-m']

    f1, ax = plt.subplots(1, 2, figsize=(13, 5))
    for i, ltres in enumerate(ltress):
        pressure = sm.collect_observations_as_time_series(ltres, 'pressure')
        saturation = sm.collect_observations_as_time_series(ltres, 'CO2saturation')

        print('Pressure at location {}:'.format(i+1), pressure, sep='\n')
        print('CO2 saturation at location {}:'.format(i+1), saturation, sep='\n')
        print('__________________________________________________________________')
        ax[0].plot(time_array/365.25, pressure/1.0e+6,
                   linespec[i], label='wellbore {}'.format(i+1))
        ax[0].set_xlabel('Time, [years]')
        ax[0].set_ylabel('Pressure, [MPa]')
        ax[0].set_title(r'Pressure at the wellbore')
        ax[1].plot(time_array/365.25, saturation,
                   linespec[i], label='wellbore {}'.format(i+1))
        ax[1].set_xlabel('Time, [years]')
        ax[1].set_ylabel(r'CO$_2$ saturation, [-]')
        ax[1].set_title(r'CO$_2$ at the wellbore')

    # Setup legend on sublots
    ax[0].legend()
    ax[1].legend()
    plt.tight_layout()

    print('____________ CO2 and brine leakage rates to aquifer 1 ____________')
    print('')

    f2, ax = plt.subplots(1, 2, figsize=(13, 5))
    for i, ms in enumerate(mss):
        CO2_aquifer1 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer1')
        brine_aquifer1 = sm.collect_observations_as_time_series(ms, 'brine_aquifer1')

        print('CO2 leakage rates to aquifer 1 at location {}:'.format(i+1),
              CO2_aquifer1, sep='\n')
        print('Brine leakage rates to aquifer 1 at location {}:'.format(i),
              brine_aquifer1, sep='\n')
        print('__________________________________________________________________')
        ax[0].plot(time_array/365.25, CO2_aquifer1,
                   linespec[i], label='wellbore '+str(i+1))
        ax[0].set_xlabel('Time, [years]')
        ax[0].set_ylabel('Leakage rates, [kg/s]')
        ax[0].set_title(r'CO$_2$ leakage to aquifer 1')
        ax[1].plot(time_array/365.25, brine_aquifer1,
                   linespec[i], label='wellbore '+str(i+1))
        ax[1].set_xlabel('Time, [years]')
        ax[1].set_ylabel('Leakage rates, [kg/s]')
        ax[1].set_title('Brine leakage to aquifer 1')
        for j in range(2):
            ax[j].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    # Setup legend on sublots
    ax[0].legend()
    ax[1].legend()
    plt.tight_layout()
    plt.show()

    num_samples = 100
    ncpus = 5
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=291)

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)

    print('_________________ Results of stochastic simulations ______________')

    # Extract results from stochastic simulations
    pressure = {}
    CO2_saturation = {}

    for ind in range(num_wells):
        pressure[ind] = np.ones((num_samples, len(time_array)))
        CO2_saturation[ind] = np.ones((num_samples, len(time_array)))

        for t_ind in range(len(time_array)):
            pressure[ind][:, t_ind] = s.recarray[
                'ltres{}.pressure_{}'.format(ind+1, t_ind)]
            CO2_saturation[ind][:, t_ind] = s.recarray[
                'ltres{}.CO2saturation_{}'.format(ind+1, t_ind)]

    # Plot results
    fontsize = 16
    fig = plt.figure(figsize=(13, 5))
    ax = fig.add_subplot(121)
    well_colors = ['blue', 'maroon', 'orange', 'green']
    for ind in range(num_wells):
        plt.plot(time_array/365.25, pressure[ind][0]/1.0e+6,
                 color=well_colors[ind], linewidth=1, label='subdomain {}'.format(ind+1))
        for j in range(1, num_samples):
            plt.plot(time_array/365.25, pressure[ind][j]/1.0e+6,
                     color=well_colors[ind], linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=fontsize-2)
    plt.ylabel('Pressure, P (MPa)', fontsize=fontsize-2)
    plt.title('Pressure: Randomly allocated wells', fontsize=fontsize)
    plt.legend()
    plt.tight_layout()
    plt.tick_params(labelsize=fontsize-4)
    plt.xlim([0, 50])
    plt.ylim([22, 36])
    ax.get_yaxis().set_label_coords(-0.10, 0.5)

    ax = fig.add_subplot(122)
    for ind in range(num_wells):
        plt.plot(time_array/365.25, CO2_saturation[ind][0],
                 color=well_colors[ind], linewidth=1, label='subdomain {}'.format(ind+1))
        for j in range(1, num_samples):
            plt.plot(time_array/365.25, CO2_saturation[ind][j],
                     color=well_colors[ind], linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=fontsize-2)
    plt.ylabel('Saturation, S (-)', fontsize=fontsize-2)
    plt.title(r'CO$_2$ saturation: Randomly allocated wells', fontsize=fontsize)
    plt.legend()
    plt.tight_layout()
    plt.tick_params(labelsize=fontsize-4)
    plt.xlim([0, 50])
    plt.ylim([0.0, 0.5])
    ax.get_yaxis().set_label_coords(-0.10, 0.5)
