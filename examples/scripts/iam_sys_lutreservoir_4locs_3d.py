"""
Example illustrates system model with multiple lookup table reservoir components
placed at four different locations.
Example illustrates use of 3d reservoir data.

Example of run:
$ python iam_sys_lutreservoir_4locs_3d.py
"""

import sys
import os
import logging
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, ReservoirDataInterpolator,
                     LookupTableReservoir)


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)

    file_directory = os.sep.join(['..', '..', 'source', 'components',
                                  'reservoir', 'lookuptables', 'Test_3d'])

    # Define keyword arguments of the system model
    num_years = 20
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Read file with signatures of interpolators and names of files with the corresponding data
    signature_data = np.genfromtxt(
        os.path.join(file_directory, 'parameters_and_filenames.csv'),
        delimiter=",", dtype='str')

    # Determine number of interpolators (should be 2)
    num_interpolators = signature_data.shape[0]-1  # -1 since the first line is a header

    # Create and add interpolators to the system model
    # Note an option interp_2d setup to False to indicate that the data is
    # to be treated as 3d
    for ind in range(num_interpolators):
        signature = {}

        intpr = sm.add_interpolator(ReservoirDataInterpolator(
            name='int'+str(ind+1), parent=sm,
            header_file_dir=file_directory,
            time_file='time_points.csv',
            data_file='Reservoir_data_sim{ind1:02}.csv'.format(ind1=ind+1),
            index=int(signature_data[ind+1, 0]),
            signature=signature,
            interp_2d=False),   # new parameter interp_2d
                                    intr_family='reservoir')

        msg = 'Signature of the created interpolator is {}'.format(signature)
        logging.debug(msg)

    logging.debug('All interpolators are created')

    # Setup locations
    well_xyzs = np.array([
        [200.0, 200.0, 15.0],
        [100.0, 300.0, 11.0],
        # location on the grid boundary
        [169.636959999974, 194.519359999802, 18.2546239999997],
        # location, the internal grid point
        [169.636959999974, 113.259679999668, 13.210184]])

    # Determine number of different locations
    num_locs = well_xyzs.shape[0]

    # Initialize list of reservoir components
    ltress = []

    # Add reservoir component for each location
    for i in range(num_locs):
        # Add reservoir component with scalar observations
        ltress.append(sm.add_component_model_object(LookupTableReservoir(
            name='ltres'+str(i+1), parent=sm,
            intr_family='reservoir',
            locX=well_xyzs[i, 0], locY=well_xyzs[i, 1], locZ=well_xyzs[i, 2],
            interp_2d=False)))

        # Add observations of reservoir component model
        ltress[-1].add_obs('pressure')
        ltress[-1].add_obs('CO2saturation')
        ltress[-1].add_obs_to_be_linked('pressure')
        ltress[-1].add_obs_to_be_linked('CO2saturation')

        # Add parameters of reservoir component models
        ltress[-1].add_par('index', value=1, vary=False)

    sm.forward()

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    print('====================== PRESSURE AND SATURATION ===================')
    print('')

    # Print and plot pressure and saturation at all locations
    linespec = ['-k', '-r', '-b', '-g', '-y', '-m']
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    fig.subplots_adjust(left=0.08, right=0.98, wspace=0.22)

    for i, ltres in enumerate(ltress):
        pressure = sm.collect_observations_as_time_series(ltres, 'pressure')
        saturation = sm.collect_observations_as_time_series(ltres, 'CO2saturation')

        print('Pressure at location {}:'.format(i+1), pressure, sep='\n')
        print('CO2 saturation at location {}:'.format(i+1), saturation, sep='\n')
        print('-----------------------------------')
        ax[0].plot(time_array/365.25, pressure/1.0e+6, linespec[i],
                   label='wellbore '+str(i+1))
        ax[0].set_xlabel('Time, [years]')
        ax[0].set_ylabel('Pressure, [MPa]')
        ax[1].plot(time_array/365.25, saturation, linespec[i],
                   label='wellbore '+str(i+1))
        ax[1].set_xlabel('Time, [years]')
        ax[1].set_ylabel(r'CO$_2$ saturation, [-]')

    # Setup legend on sublots
    ax[0].legend()
    ax[1].legend()
