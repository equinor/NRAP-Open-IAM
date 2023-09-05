"""
Example illustrates simple linking of lookup table reservoir and multisegmented
wellbore component.

This example requires the additional Kimberlina data set.
Kimberlina data set (pressure-in-pa-kimb-54-sims.zip or
Pressure_In_Pa_Kimb_54_sims.zip) can be downloaded from one of the following places:
1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam
2. https://gitlab.com/NRAP/Kimberlina_data

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/Kimb_54_sims

Example of run:
$ python iam_sys_lutreservoir_mswell.py
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
    ltres = sm.add_component_model_object(LookupTableReservoir(
        name='ltres', parent=sm,
        intr_family='reservoir', locX=37478.0, locY=48333.0))

    # Add parameters of reservoir component model
    for j in range(num_pars):
        # Add arbitrary line (here, 2) of values from signature_file
        ltres.add_par(par_names[j], value=float(signature_data[2, j+1]), vary=False)

    # Add observations of reservoir component model
    ltres.add_obs('pressure')
    ltres.add_obs('CO2saturation')

    # Add observations to be used as input for multisegmented wellbore component
    ltres.add_obs_to_be_linked('pressure')
    ltres.add_obs_to_be_linked('CO2saturation')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))

    # Add parameters of multisegmented wellbore component
    ms.add_par('wellRadius', min=0.01, max=0.02, value=0.015)
    ms.add_par('numberOfShaleLayers', value=4, vary=False)
    ms.add_par('shale1Thickness', value=200., vary=False)
    ms.add_par('shale2Thickness', value=550., vary=False)
    ms.add_par('shale3Thickness', value=400., vary=False)
    ms.add_par('shale4Thickness', value=400., vary=False)
    ms.add_par('aquifer1Thickness', value=150., vary=False)
    ms.add_par('aquifer2Thickness', value=720., vary=False)
    ms.add_par('aquifer3Thickness', value=400., vary=False)
    ms.add_par('reservoirThickness', value=400., vary=False)
    ms.add_par('logWellPerm', min=-14.0, max=-11.0, value=-13.0)

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', ltres.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', ltres.linkobs['CO2saturation'])

    # Add observations of multisegmented wellbore component model
    ms.add_obs('CO2_aquifer1')
    ms.add_obs('CO2_aquifer2')
    ms.add_obs('CO2_aquifer3')
    ms.add_obs('brine_aquifer1')
    ms.add_obs('brine_aquifer2')
    ms.add_obs('brine_aquifer3')

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    # Collect pressure and saturation observations
    pressure = sm.collect_observations_as_time_series(ltres, 'pressure')
    saturation = sm.collect_observations_as_time_series(ltres, 'CO2saturation')
    CO2_aquifer1 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer1')

    print('Pressure', pressure, sep='\n')
    print('CO2saturation', saturation, sep='\n')
    print('CO2_aquifer1', CO2_aquifer1, sep='\n')
    print('CO2_aquifer2', sm.collect_observations_as_time_series(ms, 'CO2_aquifer2'), sep='\n')
    print('CO2_aquifer3', sm.collect_observations_as_time_series(ms, 'CO2_aquifer3'), sep='\n')

    # Plot results
    fig = plt.figure(figsize=(18, 4))
    ax = fig.add_subplot(131)
    plt.plot(time_array/365.25, pressure/1.0e+6, color="maroon", linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14, fontweight='bold')
    plt.ylabel('Pressure, P (MPa)', fontsize=14, fontweight='bold')
    plt.title('Pressure: leaking well', fontsize=18, fontweight='bold')
    plt.tight_layout()
    plt.tick_params(labelsize=12)
    plt.xlim([0, 50])
    plt.ylim([26, 36])
    ax.get_yaxis().set_label_coords(-0.12, 0.5)

    ax = fig.add_subplot(132)
    plt.plot(time_array/365.25, saturation, color="green", linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14, fontweight='bold')
    plt.ylabel('Saturation, S (-)', fontsize=14, fontweight='bold')
    plt.title(r'CO$_2$ saturation: leaking well', fontsize=18, fontweight='bold')
    plt.tight_layout()
    plt.tick_params(labelsize=12)
    plt.xlim([0, 50])
    plt.ylim([0.0, 0.5])
    ax.get_yaxis().set_label_coords(-0.12, 0.5)

    ax = fig.add_subplot(133)
    plt.plot(time_array/365.25, CO2_aquifer1, color="green", linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14, fontweight='bold')
    plt.ylabel('Rate, (kg/s)', fontsize=14, fontweight='bold')
    plt.title(r'CO$_2$ leakage rates to aquifer 1', fontsize=18, fontweight='bold')
    plt.tight_layout()
    plt.tick_params(labelsize=12)
    plt.xlim([0, 50])
    plt.ylim([0.0, 1.4e-5])
    ax.get_yaxis().set_label_coords(-0.12, 0.5)
