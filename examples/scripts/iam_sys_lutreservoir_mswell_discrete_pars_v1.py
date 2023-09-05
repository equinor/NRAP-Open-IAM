"""
Example illustrates simple linking of lookup table reservoir and multisegmented
wellbore component. Example also shows use of discrete parameters for
lookup table reservoir component.

This example requires the additional Kimberlina data set.
Kimberlina data set (pressure-in-pa-kimb-54-sims.zip or
Pressure_In_Pa_Kimb_54_sims.zip) can be downloaded from one of the following places:
1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam
2. https://gitlab.com/NRAP/Kimberlina_data

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/Kimb_54_sims

Example of run:
$ python iam_sys_lutreservoir_mswell_discrete_pars_v1.py
"""

import sys
import os
import logging
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, ReservoirDataInterpolator,
                     LookupTableReservoir, MultisegmentedWellbore)


if __name__ == '__main__':
    # For multiprocessing in Spyder
    __spec__ = None
    logging.basicConfig(level=logging.WARNING)

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

    num_interpolators = signature_data.shape[0]-1  # -1 since the first line is a header

    # Create and add interpolators to the system model
    for ind in range(num_interpolators):

        signature = {par_names[j]: float(signature_data[ind+1, j+1]) for j in range(num_pars)}

        intpr = sm.add_interpolator(ReservoirDataInterpolator(
            name='int'+str(ind+1), parent=sm,
            header_file_dir=file_directory,
            time_file='time_points.csv',
            data_file='Reservoir_data_sim{ind1:02}.csv'.format(ind1=ind+1),
            index=int(signature_data[ind+1, 0]),
            signature=signature), intr_family='reservoir')

        debug_msg = 'Signature of the created interpolator is {}'.format(signature)
        logging.debug(debug_msg)

    logging.debug('All interpolators are created')

    # Add reservoir component
    ltres1 = sm.add_component_model_object(LookupTableReservoir(
        name='ltres1', parent=sm,
        intr_family='reservoir', locX=37478.0, locY=48333.0))

    par_values = [[-13.3, -12.3], [0.215, 0.338], [-18.7, -16.7]]

    # Add parameters of reservoir component model; all are assumed to have discrete distribution
    for j in range(num_pars):
        ltres1.add_par(par_names[j], value=par_values[j][0],
                       discrete_vals=(par_values[j], [0.5, 0.5]))

    # Add observations of reservoir component model
    ltres1.add_obs('pressure')
    ltres1.add_obs('CO2saturation')
    ltres1.add_obs_to_be_linked('pressure')
    ltres1.add_obs_to_be_linked('CO2saturation')

    # Add second reservoir component
    ltres2 = sm.add_component_model_object(LookupTableReservoir(
        name='ltres2', parent=sm,
        intr_family='reservoir', locX=37578.0, locY=49333.0))
    for j in range(num_pars):
        ltres2.add_par_linked_to_par(par_names[j], ltres1.pars[par_names[j]])
    ltres2.add_obs('pressure')
    ltres2.add_obs('CO2saturation')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))

    # Add parameters of multisegmented wellbore component
    ms.add_par('wellRadius', min=0.05, max=0.09, value=0.06)
    ms.add_par('numberOfShaleLayers', value=4, vary=False)
    ms.add_par('shale1Thickness', value=250., vary=False)
    ms.add_par('shale2Thickness', value=550., vary=False)
    ms.add_par('shale3Thickness', value=400., vary=False)
    ms.add_par('shale4Thickness', value=400., vary=False)
    ms.add_par('aquifer1Thickness', value=150., vary=False)
    ms.add_par('aquifer2Thickness', value=720., vary=False)
    ms.add_par('aquifer3Thickness', value=400., vary=False)
    ms.add_par('reservoirThickness', value=400., vary=False)
    ms.add_par('logWellPerm', min=-14.0, max=-11.0, value=-13.0)

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', ltres1.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', ltres1.linkobs['CO2saturation'])

    # Add observations of multisegmented wellbore component model
    ms.add_obs('CO2_aquifer1')
    ms.add_obs('CO2_aquifer2')
    ms.add_obs('CO2_aquifer3')
    ms.add_obs('brine_aquifer1')
    ms.add_obs('brine_aquifer2')
    ms.add_obs('brine_aquifer3')

    import random
    num_samples = 100
    ncpus = 10

    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=random.randint(500, 1100))

    debug_msg = 'Created sample set with samples {}'.format(s.samples.values)
    logging.debug(debug_msg)

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)

    fig = plt.figure(figsize=(15, 12))
    ax = fig.add_subplot(221)
    ind_list = list(range(len(time_array)))
    pressures = np.array([s.recarray['ltres1.pressure_'+str(indd)] for indd in ind_list])
    plt.plot(time_array/365.25, pressures, linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel('Pressure (Pa)', fontsize=14)
    plt.title(r'Bottomhole Pressure', fontsize=18)
    plt.tick_params(labelsize=12)
    plt.xlim([0, 50])

    ax = fig.add_subplot(222)
    sats = np.array([s.recarray['ltres1.CO2saturation_'+str(indd)] for indd in ind_list])
    plt.plot(time_array/365.25, sats, linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel('CO$_2$ Saturation', fontsize=14)
    plt.title(r'CO$_2$ Saturation', fontsize=18)
    plt.tick_params(labelsize=12)
    plt.xlim([0, 50])

    ax = fig.add_subplot(223)
    lrq = np.array([s.recarray['ms.CO2_aquifer1_'+str(indd)] for indd in ind_list])
    plt.semilogy(time_array/365.25, lrq, linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel('Leakage rates, q (kg/s)', fontsize=14)
    plt.title(r'Leakage of CO$_2$: aquifer 1', fontsize=18)
    plt.tick_params(labelsize=12)
    plt.xlim([0, 50])

    ax = fig.add_subplot(224)
    lbrq = np.array([s.recarray['ms.brine_aquifer1_'+str(indd)] for indd in ind_list])
    plt.plot(time_array/365.25, lbrq, linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel('Leakage rates, q (kg/s)', fontsize=14)
    plt.title(r'Leakage of brine: aquifer 1', fontsize=18)
    plt.tick_params(labelsize=12)
    plt.xlim([0, 50])

    fig = plt.figure(figsize=(15, 5))
    ax = fig.add_subplot(121)
    ind_list = list(range(len(time_array)))
    pressures = np.array([s.recarray['ltres2.pressure_'+str(indd)] for indd in ind_list])
    plt.plot(time_array/365.25, pressures, linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel('Pressure (Pa)', fontsize=14)
    plt.title(r'Observation Well Pressure', fontsize=18)
    plt.tick_params(labelsize=12)
    plt.xlim([0, 50])

    ax = fig.add_subplot(122)
    ind_list = list(range(len(time_array)))
    sats = np.array([s.recarray['ltres2.CO2saturation_'+str(indd)] for indd in ind_list])
    plt.plot(time_array/365.25, sats, linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel('CO$_2$ Saturation', fontsize=14)
    plt.title(r'Observation Well CO$_2$ Saturation', fontsize=18)
    plt.tick_params(labelsize=12)
    plt.xlim([0, 50])
