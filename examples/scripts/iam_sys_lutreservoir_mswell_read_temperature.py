"""
Example illustrates application of LookuptableReservoir to read pressure,
CO2saturation, and temperature data from a .csv file.

This example uses a modified version of the file Reservoir_data_sim01.csv from the
Kimberlina data set (https://gitlab.com/NRAP/Kimberlina_data). The file
should be located in the following location:
    source/components/reservoir/lookuptables/Test_metrics
The file from the Kimberlina data set was modified to include columns
for the fields temperature_1 through temperature_10 (corresponding to the times
in columns 1 through 10 of the time_points.csv file in the same folder).
These temperature columns were made to all have a constant value of 50 (degrees Celcius).
This example also uses a modified version of the parameters_and_filenames.csv file.
This edited version (parameters_and_filenames_read_temperature.csv) should
also be located in the Test_metrics folder and only has one entry corresponding
to the file Reservoir_data_sim01_read_temperature.csv. The parameters
are unchanged from the Kimberlina data set, however (logResPerm: -13.3,
reservoirPorosity: 0.215, logShalePerm: -18.7).
If the above files are not present, this example will produce an error.
Although the temperature data are constant, this approach would also work
on metrics that vary over space and time.

Example of run:
$ python iam_sys_lutreservoir_mswell_read_temperature.py

Last Modified: August 5th, 2022
Author: Nate Mitchell
nathaniel.mitchell@netl.doe.gov
"""

import sys
import os
import logging
import numpy as np
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, ReservoirDataInterpolator,
                     LookupTableReservoir, MultisegmentedWellbore)


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)

    file_directory = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                  'lookuptables', 'Test_metrics'])

    if not os.path.exists(os.sep.join([
            file_directory, 'Reservoir_data_sim01_read_temperature.csv'])):
        msg = 'Required data files are missing. Please check example description.'
        logging.error(msg)

    # Define keyword arguments of the system model
    num_years = 30
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Read file with signatures of interpolators and names of files with the corresponding data
    signature_data = np.genfromtxt(
        os.path.join(file_directory, 'parameters_and_filenames_read_temperature.csv'),
        delimiter=",", dtype='str')

    # The first row (except the last element) of the file contains names of the parameters
    par_names = signature_data[0, 1:-1]
    num_pars = len(par_names)

    num_interpolators = signature_data.shape[0]-1  # -1 since the first line is a header

    # Create and add interpolator to the system model
    # This example only uses the one .csv file that has been modified to include temperature
    # data (Reservoir_data_sim01_read_temperature.csv), so ind is only evaluated as 0.
    for ind in range(0,1):

        signature = {par_names[j]: float(signature_data[ind+1, j+1]) for j in range(num_pars)}

        intpr = sm.add_interpolator(ReservoirDataInterpolator(
            name='int'+str(ind+1), parent=sm,
            header_file_dir=file_directory,
            time_file='time_points.csv',
            data_file='Reservoir_data_sim{ind1:02}_read_temperature.csv'.format(ind1=ind+1),
            index=int(signature_data[ind+1, 0]),
            signature=signature), intr_family='reservoir')

        debug_msg = 'Signature of the created interpolator is {}'.format(signature)
        logging.debug(debug_msg)

    logging.debug('All interpolators are created')

    # Add reservoir component
    ltres = sm.add_component_model_object(LookupTableReservoir(
        name='ltres', parent=sm, intr_family='reservoir',
        locX=37478.0, locY=48333.0,
        parameter_filename='parameters_and_filenames_read_temperature.csv'))

    # Add parameters of reservoir component model
    for j in range(num_pars):
        # Add arbitrary line (here, 1) of values from signature_file
        ltres.add_par(par_names[j], value=float(signature_data[1, j+1]), vary=False)

    # Add observations of reservoir component model
    ltres.add_obs('pressure')
    ltres.add_obs('CO2saturation')
    ltres.add_obs('temperature')

    # Add observations to be used as input for multisegmented wellbore component
    ltres.add_obs_to_be_linked('pressure')
    ltres.add_obs_to_be_linked('CO2saturation')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))

    # Add parameters of multisegmented wellbore component
    ms.add_par('wellRadius', value=0.015, vary=False)
    ms.add_par('numberOfShaleLayers', value=4, vary=False)
    ms.add_par('shale1Thickness', value=200., vary=False)
    ms.add_par('shale2Thickness', value=550., vary=False)
    ms.add_par('shale3Thickness', value=400., vary=False)
    ms.add_par('shale4Thickness', value=400., vary=False)
    ms.add_par('aquifer1Thickness', value=150., vary=False)
    ms.add_par('aquifer2Thickness', value=720., vary=False)
    ms.add_par('aquifer3Thickness', value=400., vary=False)
    ms.add_par('reservoirThickness', value=400., vary=False)
    ms.add_par('logWellPerm', value=-13.0, vary=False)

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
    temp = sm.collect_observations_as_time_series(ltres, 'temperature')
    CO2_aquifer1 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer1')
    CO2_aquifer2 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer2')
    CO2_aquifer3 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer3')

    print('Year', time_array / 365.25, sep='\n')
    print('Pressure', pressure, sep='\n')
    print('CO2saturation', saturation, sep='\n')
    print('CO2_aquifer1', CO2_aquifer1, sep='\n')
    print('CO2_aquifer2', CO2_aquifer2, sep='\n')
    print('CO2_aquifer3', CO2_aquifer3, sep='\n')
    print('Temperature', temp, sep='\n')
