"""
Example illustrates system model with multiple lookup table reservoir components
each linked to a separate multisegmented wellbore component. 4 of the 5 wellbore
locations are random; remaining wellbore is placed at the predetermined location.
Example also illustrates an option of using a gridded observation of
lookup table reservoir component.

This example also illustrates how to save all the outputs produced by the simulation.
Changing variable 'save_output' at the beginning of simulation (line 43)
from True to False allows to cancel saving of the outputs.

This example requires the additional Kimberlina data set.
Kimberlina data set (pressure-in-pa-kimb-54-sims.zip or
Pressure_In_Pa_Kimb_54_sims.zip) can be downloaded from one of the following places:
1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam
2. https://gitlab.com/NRAP/Kimberlina_data

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/Kimb_54_sims

Example of run:
$ python iam_sys_lutreservoir_mswell_5locs_gridded_obs.py
"""

import sys
import os
import logging
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, ReservoirDataInterpolator,
                     LookupTableReservoir, MultisegmentedWellbore)
from matk import pyDOE


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)

    # Change the variable value to False if saving the outputs is not needed.
    # By default, the results will be saved in the folder 'output/csv_files' within root
    # folder of NRAP-Open-IAM
    save_output = True

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
            header_file_dir=os.path.join('..', '..', 'source', 'components',
                                         'reservoir', 'lookuptables', 'Kimb_54_sims'),
            time_file='time_points.csv',
            data_file='Reservoir_data_sim{ind1:02}.csv'.format(ind1=ind+1),
            index=int(signature_data[ind+1, 0]),
            signature=signature), intr_family='reservoir')
        debug_msg = 'Signature of the created interpolator is {}'.format(signature)
        logging.debug(debug_msg)

    logging.debug('All interpolators are created')

    # Create randomly located leaky well locations within box defined by
    # xmin,xmax,ymin,ymax
    num_rand_locs = 4
    xymins = np.array([37300., 48200.])
    xymaxs = np.array([37600., 48500.])
    well_xys = xymins + pyDOE.lhs(2, samples=num_rand_locs)*(xymaxs-xymins)

    # Initialize list of reservoir components
    ltress = []

    # Extra reservoir (named ltres1) will correspond to the LUT component with gridded observation
    for i in range(num_rand_locs+1):
        if i > 0:
            # Add reservoir component with scalar observations
            ltress.append(sm.add_component_model_object(LookupTableReservoir(
                name='ltres'+str(i+1), parent=sm,
                intr_family='reservoir', locX=well_xys[i-1, 0], locY=well_xys[i-1, 1])))

            # Add observations of reservoir component model
            ltress[-1].add_obs('pressure')
            ltress[-1].add_obs('CO2saturation')
            ltress[-1].add_obs_to_be_linked('pressure')
            ltress[-1].add_obs_to_be_linked('CO2saturation')
        else:
            # Add reservoir component with gridded observations
            ltress.append(sm.add_component_model_object(
                LookupTableReservoir(name='ltres'+str(i+1), parent=sm,
                                     intr_family='reservoir')))   # no location is specified

            # Optional keyword argument constr_type is used to illustrate
            # the type of the gridded observation to be linked. Since the source
            # of gridded data is a lookup table (already known data)
            # we won't use add_grid_obs method, i.e., we will not be saving
            # the gridded observation in the file
            ltress[-1].add_obs_to_be_linked(
                'pressure', obs_type='grid', constr_type='array')
            ltress[-1].add_obs_to_be_linked(
                'CO2saturation', obs_type='grid', constr_type='array')
            # We'll use add_local_obs to track saturation and pressure
            # at the predefined wellbore location
            loc_ind = 2535
            ltress[-1].add_local_obs('loc_pressure', grid_obs_name='pressure',
                                     constr_type='array', loc_ind=loc_ind)
            ltress[-1].add_local_obs('loc_CO2saturation', grid_obs_name='CO2saturation',
                                     constr_type='array', loc_ind=loc_ind)

        # Add parameters of reservoir component models
        for j in range(num_pars):
            # 2 is an index of an arbitrary row in the signature file
            ltress[-1].add_par(par_names[j], value=float(signature_data[2, j+1]), vary=False)

    # Initialize list of wellbore components
    mss = []
    for i in range(num_rand_locs+1):
        mss.append(sm.add_component_model_object(MultisegmentedWellbore(
            name='ms'+str(i+1), parent=sm)))

        # Add parameters of multisegmented wellbore components
        mss[-1].add_par('wellRadius', min=0.01, max=0.02, value=0.015)
        mss[-1].add_par('numberOfShaleLayers', value=4, vary=False)
        mss[-1].add_par('shale1Thickness', value=250., vary=False)
        mss[-1].add_par('shale2Thickness', value=550., vary=False)
        mss[-1].add_par('shale3Thickness', value=400., vary=False)
        mss[-1].add_par('shale4Thickness', value=400., vary=False)
        mss[-1].add_par('aquifer1Thickness', value=150., vary=False)
        mss[-1].add_par('aquifer2Thickness', value=720., vary=False)
        mss[-1].add_par('aquifer3Thickness', value=400., vary=False)
        mss[-1].add_par('reservoirThickness', value=400., vary=False)
        mss[-1].add_par('logWellPerm', min=-14.0, max=-11.0, value=-13.0+np.random.rand())

        # Add observations of multisegmented wellbore component model
        mss[-1].add_obs('CO2_aquifer1')
        mss[-1].add_obs('CO2_aquifer2')
        mss[-1].add_obs('CO2_aquifer3')
        mss[-1].add_obs('brine_aquifer1')
        mss[-1].add_obs('brine_aquifer2')
        mss[-1].add_obs('brine_aquifer3')

        if i > 0:
            # Add keyword arguments linked to the output provided by reservoir model
            mss[-1].add_kwarg_linked_to_obs('pressure', ltress[i].linkobs['pressure'])
            mss[-1].add_kwarg_linked_to_obs('CO2saturation', ltress[i].linkobs['CO2saturation'])
        else:
            # Add keyword arguments linked to the output provided by reservoir model
            # constr_type and loc_ind are used to specify location at which
            # the scalar observation is extracted. These arguments are
            # component specific: they depend on the kind
            # of gridded observation returned by the supplying component.
            # Even for a single location loc_ind should be a list.
            mss[-1].add_kwarg_linked_to_obs('pressure', ltress[i].linkobs['pressure'],
                                            obs_type='grid', constr_type='array',
                                            loc_ind=[loc_ind])
            mss[-1].add_kwarg_linked_to_obs('CO2saturation', ltress[i].linkobs['CO2saturation'],
                                            obs_type='grid', constr_type='array',
                                            loc_ind=[loc_ind])

    sm.forward(save_output=save_output)

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    print('====================== PRESSURE AND SATURATION ===================')
    print('')

    # Print pressure and saturation
    linespec = ['-k', '-r', '-b', '-g', '-y', '-m']
    f1, ax = plt.subplots(1, 2, figsize=(12, 4))

    for i, ltres in enumerate(ltress):
        if i > 0:
            pressure = sm.collect_observations_as_time_series(ltres, 'pressure')
            saturation = sm.collect_observations_as_time_series(ltres, 'CO2saturation')
        else:
            pressure = sm.collect_observations_as_time_series(ltres, 'loc_pressure')
            saturation = sm.collect_observations_as_time_series(ltres, 'loc_CO2saturation')

        print('Pressure at location {}:'.format(i+1), pressure, sep='\n')
        print('CO2 saturation at location {}:'.format(i+1), saturation, sep='\n')
        print('-----------------------------------')
        ax[0].plot(time_array/365.25, pressure/1.0e+6, linespec[i], label='wellbore '+str(i+1))
        ax[0].set_xlabel('Time, [years]')
        ax[0].set_ylabel('Pressure, [MPa]')
        ax[1].plot(time_array/365.25, saturation, linespec[i], label='wellbore '+str(i+1))
        ax[1].set_xlabel('Time, [years]')
        ax[1].set_ylabel(r'CO$_2$ saturation, [-]')

    # Setup legend on sublots
    ax[0].legend()
    ax[1].legend()

    print('============= CO2 and BRINE LEAKAGE RATES TO AQUIFER 1 ===========')
    print('')
    f2, ax = plt.subplots(1, 2, figsize=(12, 4))

    for i, ms in enumerate(mss):
        CO2_aquifer1 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer1')
        brine_aquifer1 = sm.collect_observations_as_time_series(ms, 'brine_aquifer1')

        print('CO2 leakage rates to aquifer 1 at location {}:'.format(i+1),
              CO2_aquifer1, sep='\n')
        print('Brine leakage rates to aquifer 1 at location {}:'.format(i),
              brine_aquifer1, sep='\n')
        print('-----------------------------------')
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
    plt.show()
