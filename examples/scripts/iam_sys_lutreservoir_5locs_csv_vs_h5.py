"""
Example illustrates system model with multiple lookup table reservoir components
placed at five different locations.

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

from openiam import SystemModel, ReservoirDataInterpolator, LookupTableReservoir
from matk import pyDOE

def create_system_model(well_xys, time_array, option='csv'):
    """
    Create system model with provided inputs using csv or h5 lookup tables
    for Lookup Table Reservoir component.

    Parameters
    ----------
    well_xys : TYPE
        DESCRIPTION.
    time_array : TYPE
        DESCRIPTION.
    option : TYPE, optional
        DESCRIPTION. The default is 'csv'.

    Returns
    -------
    sm : Instance of SystemModel class
        DESCRIPTION.
    ltress : list of instances of LookupTableReservoir class
        DESCRIPTION.

    """
    sm_model_kwargs = {'time_array': time_array}

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Read file with signatures of interpolators and names of files with the corresponding data
    # Signature defines set of parameters used when creating a particular lookup table
    if option == 'csv':
        signature_data = np.genfromtxt(
            os.path.join(
                '..', '..', 'source', 'components', 'reservoir',
                'lookuptables', 'Kimb_54_sims', 'parameters_and_filenames.csv'),
            delimiter=",", dtype='str')
    else:
        signature_data = np.genfromtxt(
            os.path.join(
                '..', '..', 'source', 'components', 'reservoir',
                'lookuptables', 'Kimb_54_sims', 'parameters_and_filenames_h5.csv'),
            delimiter=",", dtype='str')

    # The first row (except the last element) of the file contains names of the parameters
    # First column (column=0) in the signature_data is for index so the parameter
    # names start from the second column (column=1)
    par_names = signature_data[0, 1:-1]
    num_pars = len(par_names)

    num_interpolators = signature_data.shape[0]-1  # -1 since the first line is a header

    if option == 'csv':
        file_ext = 'csv'
    else:
        file_ext = 'h5'

    # Create and add interpolators to the system model
    for ind in range(num_interpolators):
        signature = {par_names[j]: float(signature_data[ind+1, j+1]) for j in range(num_pars)}

        if ind == 0:
            intrp1 = sm.add_interpolator(
                ReservoirDataInterpolator(
                    name='int'+str(ind+1), parent=sm,
                    header_file_dir=os.path.join(
                        '..', '..', 'source', 'components', 'reservoir',
                        'lookuptables', 'Kimb_54_sims'),
                    time_file='time_points.csv',
                    data_file='Reservoir_data_sim{ind1:02}.{f_ext}'.format(
                        ind1=ind+1, f_ext=file_ext),
                    index=int(signature_data[ind+1, 0]),
                    signature=signature), intr_family='reservoir')
        else:
            _ = sm.add_interpolator(
                ReservoirDataInterpolator(
                    name='int'+str(ind+1), parent=sm,
                    header_file_dir=os.path.join(
                        '..', '..', 'source', 'components', 'reservoir',
                        'lookuptables', 'Kimb_54_sims'),
                    time_file='time_points.csv',
                    data_file='Reservoir_data_sim{ind1:02}.{f_ext}'.format(
                        ind1=ind+1, f_ext=file_ext),
                    index=int(signature_data[ind+1, 0]),
                    signature=signature,
                    triangulation=intrp1.triangulation), intr_family='reservoir')

        msg = 'Signature of the created interpolator is {}'.format(signature)
        logging.debug(msg)

    logging.debug('All interpolators are created')

    # Initialize list of reservoir components
    ltress = []
    for i, crds in enumerate(well_xys):
        # Add reservoir component
        ltress.append(sm.add_component_model_object(
            LookupTableReservoir(
                name='ltres'+str(i+1), parent=sm,
                intr_family='reservoir', locX=crds[0], locY=crds[1])))

        # Add parameters of reservoir component model
        for j in range(num_pars):
            ltress[-1].add_par(par_names[j],
                               value=float(signature_data[8+i, j+1]),
                               vary=False)  # 8+i is some arbitrary row index

        # Add observations of reservoir component model
        ltress[-1].add_obs('pressure')
        ltress[-1].add_obs('CO2saturation')

    return sm, ltress

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
    time_array = 365.25*np.arange(0.0, num_years+1)   # time is defined in days

    # Create 5 randomly located leaky well locations within box defined by
    # xmin,xmax,ymin,ymax
    num_locations = 5
    xymins = np.array([37300., 48200.])
    xymaxs = np.array([37600., 48500.])
    well_xys = xymins + pyDOE.lhs(2, samples=num_locations)*(xymaxs-xymins)

    # Setup system model with lookup table reservoir component based on csv lookup tables
    csv_sm, csv_ltress = create_system_model(well_xys, time_array, option='csv')
    print('csv based system model is created')

    # Time the execution
    csv_sm.forward()
    print('csv based simulation is finished')

    # Setup system model with lookup table reservoir component based on h5 lookup tables
    h5_sm, h5_ltress = create_system_model(well_xys, time_array, option='h5')
    print('h5 based system model is created')

    # Time the execution
    h5_sm.forward()
    print('h5 based simulation is finished')

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    print('====================== PRESSURE AND SATURATION ===================')
    print('')
    # Print pressure and saturation
    linespec = ['r', 'b', 'g', 'k', 'm']
    f1, ax = plt.subplots(1, 2, figsize=(12, 4))

    for ind in range(num_locations):
        pressure1 = csv_sm.collect_observations_as_time_series(csv_ltress[ind], 'pressure')
        saturation1 = csv_sm.collect_observations_as_time_series(csv_ltress[ind], 'CO2saturation')
        pressure2 = h5_sm.collect_observations_as_time_series(h5_ltress[ind], 'pressure')
        saturation2 = h5_sm.collect_observations_as_time_series(h5_ltress[ind], 'CO2saturation')
        print('Difference in pressure:', pressure1-pressure2, sep='\n')
        print('--------------------------------------------------------------')
        print('Difference in saturation:', saturation1-saturation2, sep='\n')
        print('--------------------------------------------------------------')
        ax[0].plot(time_array/365.25, pressure1/1.0e+6, '-'+linespec[ind],
                    label='location {}'.format(ind+1))
        ax[0].plot(time_array/365.25, pressure2/1.0e+6, 'x'+linespec[ind])
        ax[0].set_xlabel('Time, [years]')
        ax[0].set_ylabel('Pressure, [MPa]')
        ax[1].plot(time_array/365.25, saturation1, '-'+linespec[ind],
                    label='location {}'.format(ind+1))
        ax[1].plot(time_array/365.25, saturation2, 'x'+linespec[ind])
        ax[1].set_xlabel('Time, [years]')
        ax[1].set_ylabel(r'CO$_2$ saturation, [-]')

    # Setup legend on sublots
    ax[0].legend()
    ax[1].legend()
    plt.show()
