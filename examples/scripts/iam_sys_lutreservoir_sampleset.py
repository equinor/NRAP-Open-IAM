"""
Example illustrates setup of lookup table reservoir component and custom sample
set for discrete parameters.

This example requires the additional Kimberlina data set.
Kimberlina data set (pressure-in-pa-kimb-54-sims.zip or
Pressure_In_Pa_Kimb_54_sims.zip) can be downloaded from one of the following places:
1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam
2. https://gitlab.com/NRAP/Kimberlina_data

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/Kimb_54_sims

Example of run:
$ python iam_sys_lutreservoir_sampleset.py
"""

import sys
import os
import logging
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, ReservoirDataInterpolator, LookupTableReservoir)


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

    # Add index parameter of reservoir component model
    num_samples = 5
    ltres.add_par('index', value=1,
                   discrete_vals=[list(range(1, num_samples+1)), 5*[1.0]])

    # Add observations of reservoir component model
    ltres.add_obs('pressure')
    ltres.add_obs('CO2saturation')

    samples = np.arange(1, num_samples+1).reshape(num_samples, 1)
    s = sm.create_sampleset(samples)

    results = s.run(cpus=5, verbose=False)

    # Extract results from stochastic simulations
    out = s.collect_observations_as_time_series()

    # Setup plot parameters
    label_size = 13
    font_size = 16
    ticks_size = 12
    line_width = 2

    # Plot results
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for ind, obs_nm in enumerate(['CO2saturation', 'pressure']):
        for sind in range(num_samples):
            ax[ind].plot(time_array/365.25, out['ltres.'+obs_nm][sind], '-',
                         linewidth=line_width, label='scenario {}'.format(sind+1))
        ax[ind].set_xlabel('Time, [years]', fontsize=label_size)
    ax[0].set_ylabel(r'CO$_2$ saturation, [-]',
                     fontsize=label_size)
    ax[1].set_ylabel(r'Pressure, [Pa]',
                     fontsize=label_size)
    ax[1].legend()
    fig.suptitle('Results of simulation', fontsize=font_size)
