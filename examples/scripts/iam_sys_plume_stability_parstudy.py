'''
Example illustrates use of plume stability component.

This example requires the additional Kimberlina data set.
Kimberlina data set (pressure-in-pa-kimb-54-sims.zip or
Pressure_In_Pa_Kimb_54_sims.zip) can be downloaded from one of the following places:
1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam
2. https://gitlab.com/NRAP/Kimberlina_data

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/Kimb_54_sims
'''
import sys
import os
import logging
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import SystemModel, PlumeStability

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

    option = 1
    if option == 1:
        num_years = 200
        time_array = 365.25*np.arange(0.0, num_years+1)
        # Time array different from time points provided in the data set
        sm_model_kwargs = {'time_array': time_array}   # time is given in days
    else:
        # Time array is the same as in the data set but converted to days
        time_array = np.genfromtxt(
            os.sep.join([file_directory, 'time_points.csv']), delimiter=',')*365.25

    # Time array argument has to be defined for use of plume stability component.
    # It has to be defined as in the data set or whatever (in the latter case
    # interpolation will be used).
    sm_model_kwargs = {'time_array': time_array}
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    sps = sm.add_component_model_object(
        PlumeStability(name='sps', parent=sm, file_directory=file_directory,
                       # variable_names must be a list (even it contains only one name)
                       variable_names=['pressure', 'CO2saturation'],
                       # thresholds must be a dictionary containing thresholds
                       # for all variables in variable_names
                       thresholds={'pressure': 1.0e6, 'CO2saturation': 0.01},
                       parameter_filename='parameters_and_filenames.csv',
                       time_file='time_points.csv'))

    # Variable(s) can be modified
    # Print available variables
    print('Variables available in data set:', sps.variable_list)
    # Print names of variables for which metrics will be calculated
    print('Metrics to be calculated for', sps.variable_names)

    # It is possible to modify variables and thresholds after the plume component
    # is added. For example, to add new variable or update threshold one can use
    #    sps.add_variable('CO2saturation', 0.01). The method will check whether the
    # variable is already in the list and then adds it if it's not there.
    # If the variable was already added, only threshold will be updated.
    # The threshold can be modified by referring directly as
    #    sps.thresholds['CO2saturation'] = 0.01
    # To remove variable from the list the following method can be used
    #    sps.variable_names.remove('pressure')

    # Add discrete parameter sample number. Each sample has the same probability
    # When parameter is added as discrete, it has an associated with it value which is used
    # for forward simulations. For range of values from 1 to 54 the value of 28
    # is used as the value located at the next position from the center of the list
    # or right at the center depending whether the number of values is even or odd.
    # For example, for list of values [3, 7, 1, 2, 9] value of 1 would be used as value
    # for forward simulation.
    sps.add_par('index', discrete_vals=(
        list(range(1, len(sps.filenames)+1)),
        [1./len(sps.filenames)]*len(sps.filenames)))
    print('Value of {} will be used for index parameter in forward simulation.'.format(
        int(sps.pars['index'].value)))

    # Observations have to be added explicitly
    sps.add_obs('times')
    for var_name in sps.variable_names:
        sps.add_obs('{}_areas'.format(var_name))
        sps.add_obs('{}_areas_dt'.format(var_name))
        sps.add_obs('{}_mobility'.format(var_name))
        sps.add_obs('{}_spreading'.format(var_name))

    # 'index' parameter is bounded from above
    # by the number of datasets in the 'file_directory'.
    # To change the dataset probability, one can do, for example
    # Collect current weights
    #    wts = sps.realization_weight
    # Modify weight of the 1st dataset
    #    wts[0] *= 2
    # Modify weight of the 10th dataset
    #    wts[9] *= 3
    # Assign modified weights back to 'index' parameter
    # Note that this will produce a warning that the weights have been renormalized to sum to one.
    #    sps.realization_weight = wts

    out = sm.forward()
    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    outputs = {}

    for var_name in sps.variable_names:
        # Get data
        outputs[var_name] = {}
        outputs[var_name]['area'] = sm.collect_observations_as_time_series(
            sps, '{}_areas'.format(var_name))
        outputs[var_name]['area_dt'] = sm.collect_observations_as_time_series(
            sps, '{}_areas_dt'.format(var_name))
        outputs[var_name]['mobility'] = sm.collect_observations_as_time_series(
            sps, '{}_mobility'.format(var_name))
        outputs[var_name]['spreading'] = sm.collect_observations_as_time_series(
            sps, '{}_spreading'.format(var_name))

        # Print collected data
        print('{} area'.format(var_name), outputs[var_name]['area'], sep='\n')
        print('------------------------------------------------------------------')
        print('{} area deriv'.format(var_name), outputs[var_name]['area_dt'], sep='\n')
        print('------------------------------------------------------------------')
        print('{} mobility'.format(var_name), outputs[var_name]['mobility'], sep='\n')
        print('------------------------------------------------------------------')
        print('{} spreading'.format(var_name), outputs[var_name]['spreading'], sep='\n')

    # Plot results
    font_size = 16
    labelx = -0.1
    title_part = {'pressure': 'Pressure', 'CO2saturation': r'CO$_2$ saturation'}
    times = time_array/365.25
    for var_name in sps.variable_names:
        f, ax = plt.subplots(4, sharex=True, figsize=(10, 10))
        f.subplots_adjust(left=0.15, bottom=0.07, right=0.95, top=0.93, wspace=0.1)
        ax[0].plot(times, outputs[var_name]['area']/1000**2)
        ax[1].plot(times, outputs[var_name]['area_dt']/1000**2)
        ax[2].plot(times, outputs[var_name]['mobility']/1000)
        ax[3].plot(times, outputs[var_name]['spreading']/1000**2)
        ax[0].set_ylabel(r'Plume area,{}[km$^2$]'.format('\n'), fontsize=font_size)
        ax[1].set_ylabel(r'Change in plume area,{}[km$^2$/year]'.format('\n'), fontsize=font_size)
        ax[2].set_ylabel(r'Mobility,{}[km/year]'.format('\n'), fontsize=font_size)
        ax[3].set_ylabel(r'Spreading,{}[km$^2$/year]'.format('\n'), fontsize=font_size)
        ax[3].set_xlabel('Time [years]', fontsize=font_size)

        for ind in range(4):
            ax[ind].tick_params(axis='both', which='major', labelsize=font_size-2)

        # Align y-labels
        for ind in range(4):
            ax[ind].yaxis.set_label_coords(labelx, 0.5)

        plt.suptitle('{} metrics'.format(title_part[var_name]), fontsize=font_size+3)
        plt.show()

    # One can run system model with a different value of parameter index using
    #    out = sm.forward(pardict={'sps.index': 1})

    # Create sampleset with 20 evenly spaced realization ids
    s = sm.parstudy(nvals=20)
    # One can create a sampleset of all realizations in Kimb_54_sims
    #    s = sm.create_sampleset([[v] for v in range(1, 55)])

    # Run sampleset
    s.run(cpus=4, verbose=False)
    # Postprocess output into easy to use dictionary
    out = s.collect_observations_as_time_series()

    # Plot results
    for var_name in sps.variable_names:
        f, ax = plt.subplots(4, sharex=True, figsize=(10, 10))
        for i in range(len(s.indices)):
            ax[0].plot(times, out['sps.{}_areas'.format(var_name)][i]/1.0e+6)
            ax[1].plot(times, out['sps.{}_areas_dt'.format(var_name)][i]/1.0e+6)
            ax[2].plot(times, out['sps.{}_mobility'.format(var_name)][i]/1.0e+3)
            ax[3].plot(times, out['sps.{}_spreading'.format(var_name)][i]/1.0e+6)

        f.subplots_adjust(left=0.15, bottom=0.07, right=0.95, top=0.93, wspace=0.1)
        ax[0].set_ylabel(r'Plume area,{}[km$^2$]'.format('\n'), fontsize=font_size)
        ax[1].set_ylabel(r'Change in plume area,{}[km$^2$/year]'.format('\n'), fontsize=font_size)
        ax[2].set_ylabel(r'Mobility,{}[km/year]'.format('\n'), fontsize=font_size)
        ax[3].set_ylabel(r'Spreading,{}[km$^2$/year]'.format('\n'), fontsize=font_size)
        ax[3].set_xlabel('Time [years]', fontsize=font_size)

        for ind in range(4):
            ax[ind].tick_params(axis='both', which='major', labelsize=font_size-2)

        # Align y-labels
        for ind in range(4):
            ax[ind].yaxis.set_label_coords(labelx, 0.5)

        plt.suptitle('{} metrics'.format(title_part[var_name]), fontsize=font_size+3)
        plt.show()

    # Remove all handlers from the logger for proper work in the consecutive runs
    while logging.getLogger('').handlers:
        logging.getLogger('').handlers.pop()
