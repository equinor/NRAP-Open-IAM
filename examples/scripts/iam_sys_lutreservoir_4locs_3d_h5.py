"""
Example illustrates system model with multiple lookup table reservoir components
placed at five different locations.
Example illustrates use of 3d reservoir data provided through h5 lookup table
for Illinois Basin Decatur Project.

Example of run:
$ python iam_sys_lutreservoir_4locs_3d_h5.py
"""

import sys
import os
import logging
# import h5py
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, ReservoirDataInterpolator,
                     LookupTableReservoir)


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)

    file_directory = os.sep.join(['..', '..', 'source', 'components',
                                  'reservoir', 'lookuptables', 'IBDP'])

    # with h5py.File(os.path.join(file_directory, 'Reservoir_sim_BCIBDP.h5'), 'r') as hf:
    #     a = list(hf.keys())
    #     x = hf['x'][()]
    #     y = hf['y'][()]
    #     z = hf['z'][()]
    #     saturation = hf['CO2saturation_60'][()]

    # xmin = 95152.64220897049     # np.min(x)
    # xmax = 113873.17170033549    # np.max(x)
    # ymin = 347431.0728592322     # np.min(y)
    # ymax = 365799.01135136606    # np.max(y)
    # zmin = -2412.664119727896    # np.min(z)
    # zmax = -1294.2823751943602   # np.max(z)

    # indices = np.where(saturation>0.1)[0]
    # selected_indices = [684363, 684364, 684365, 684366, 684468, 684469, 684470,
    #                     684471, 684472, 684473, 684474]
    # for ind in range(len(selected_indices)):
    #     print(x[ind], y[ind], z[ind])

    # Printed coordinates
    # 100000.21436359595 347431.0728592322 -2248.862731934832
    # 100662.55668743292 347672.1455554775 -2255.7614736333844
    # 101271.19572039068 347893.67205823516 -2262.6160133791923
    # 101879.83460452373 348115.1985609927 -2269.9009640640243
    # 102442.78980150161 348320.09768752794 -2276.960551758384
    # 102921.25169780018 348494.2435648072 -2283.175167772104
    # 103327.9656367273 348642.27516561415 -2288.58547939596
    # 103673.67254434386 348768.10220489383 -2293.259538280488
    # 103967.537703319 348875.0602186705 -2297.2766698235523
    # 104217.29452088742 348965.96414217906 -2300.716385448552
    # 104429.53042769543 349043.21189423034 -2303.654066602896

    # Define keyword arguments of the system model
    num_years = 5
    time_array = 365.25*np.linspace(0.0, num_years, num=6*num_years)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Read file with signatures of interpolators and names of files with the corresponding data
    signature_data = np.genfromtxt(
        os.path.join(file_directory, 'parameters_and_filenames.csv'),
        delimiter=",", dtype='str')

    # Create and add interpolators to the system model
    # Note an option interp_2d set to False to indicate that the data is 3d
    signature = {}

    intpr = sm.add_interpolator(
        ReservoirDataInterpolator(
            name='int1', parent=sm,
            header_file_dir=file_directory,
            time_file='time_points.csv',
            data_file='Reservoir_sim_BCIBDP.h5',
            index=int(signature_data[1, 0]),
            signature=signature,
            interp_2d=False),
        intr_family='reservoir')


    msg = 'Signature of the created interpolator is {}'.format(signature)
    logging.debug(msg)

    logging.debug('All interpolators are created')

    # Setup locations
    well_xyzs = np.array([
        [100001.0, 347440.0, -2248.0],
        [103967.0, 348875.0, -2297.0],
        [101880.0, 358115.0, -2269.0],
        [100921.0, 358494.0, -2283.0],
        [100000.21436359595, 347431.0728592322, -2248.862731934832]])

    # xmin = 95152.64220897049     # np.min(x)
    # xmax = 113873.17170033549    # np.max(x)
    # ymin = 347431.0728592322     # np.min(y)
    # ymax = 365799.01135136606    # np.max(y)
    # zmin = -2412.664119727896    # np.min(z)
    # zmax = -1294.2823751943602   # np.max(z)

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
