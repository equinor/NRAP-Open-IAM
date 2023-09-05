"""
The module contains the solution class for hydrocarbon leakage component.

Last Modified: March 2nd, 2023
"""

import sys
import os
import logging
import numpy as np
import tensorflow as tf
from sklearn.preprocessing import MinMaxScaler

# This suppresses the excessive messages produced by tensorflow
tf.get_logger().setLevel(logging.ERROR)

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from openiam import IAM_DIR

OUTPUT_NAMES = ['mass_CO2_aquifer', 'mass_CO2_gas_aquifer',
                'mass_methane_oil_aquifer', 'mass_methane_gas_aquifer',
                'mass_oil_aquifer', 'mass_gas_aquifer']

# The ROM filenames correspond to the order of outputs in OUTPUT_NAMES
OUTPUT_ROM_NAMES = ['CO2_o.h5', 'CO2_g.h5', 'CH4_o.h5',
                    'CH4_g.h5', 'oil.h5', 'gas.h5']

# These are the maximum output values in the training dataset. They are used to
# scale the output from the TensorFlow model. The values correspond with the
# names in OUTPUT_NAMES.
y_train_max_values = [
    878092.4019262992, 1721264.4714578395,
    112336.38179058068, 445995.58449341706,
    2811463.6629231963, 2219098.3053520545]

input_parameter_names = [
    'postInjectionTime', 'reservoirDepth', 'NTG', 'logResPerm',
    'reservoirPressureMult', 'logWellPerm', 'avgWaterSaturation',
    'FCO2', 'FC1', 'FC4', 'FC7Plus']

# The minimum and maximum training values used for the models. The 11 columns
# correspond to the entries in input_parameter_names. Here, postInjectionTime
# is in years, reservoirDepth is in feet, both the reservoir and well
# permeabilities are in mD (millidarcies).
x_training_min = [100, 3003.92, 0, 10.0855, 1.0001, 0.24313, 0, 0, 0, 0, 0]
x_training_max = [410, 8996.64, 1, 99.911, 1.19996, 999.68, 1, 1, 1, 1, 1]

MIN_TIME = 100
MAX_TIME = 410

# Assuming gravitational acceleration is 9.81 m/s^2
LBS_TO_KG = 0.45359237

# This result is given for model times outside the range of MIN_TIME to MAX_TIME.
# Such results should not be considered valid, however.
VALUE_FOR_TIMES_OUTSIDE_RANGE = 0

set_neg_results_to_zero = True

class Solution():
    def __init__(self):
        """ Create an instance of the Solution class. """
        self.outputs = np.zeros(6) # placeholder

    def find(self, input_array, output_check):
        """ Find output of the ROM for the provided input parameters."""
        self.outputs = np.zeros((input_array.shape[0], len(OUTPUT_NAMES)))

        no_results = False

        input_array_copy = np.copy(input_array)
        input_array_clip = np.where(np.logical_and(
            input_array_copy >= MIN_TIME, input_array_copy <= MAX_TIME))

        if input_array_clip[0].shape[0] > 0:
            input_values_scaled = scale_input(input_array_copy)

            model_directory = os.path.join(IAM_DIR, 'source', 'components',
                                           'wellbore', 'hydrocarbon_leakage')
            os.chdir(model_directory)

            rom_output = np.zeros((input_array_copy.shape[0],
                                   input_array_copy.shape[1]))

            # Load models
            roms = {'mass_CO2_aquifer': tf.keras.models.load_model('CO2_o.h5'),
                    'mass_CO2_gas_aquifer': tf.keras.models.load_model('CO2_g.h5'),
                    'mass_methane_oil_aquifer': tf.keras.models.load_model('CH4_o.h5'),
                    'mass_methane_gas_aquifer': tf.keras.models.load_model('CH4_g.h5'),
                    'mass_oil_aquifer': tf.keras.models.load_model('oil.h5'),
                    'mass_gas_aquifer': tf.keras.models.load_model('gas.h5')
                    }

            # Scaled predictions that range from 0 to 1 converted to tensor
            input_values_tensor = tf.convert_to_tensor(input_values_scaled,
                                                       dtype=tf.float64)

            for ind, obs_name in enumerate(OUTPUT_NAMES):
                if output_check[ind]:

                    predicts = roms[obs_name](
                        input_values_tensor, training=False).numpy().reshape(-1)

                    # The ROMs can produce values that are negative but very
                    # close to zero. In such cases, the ROM developer recommends
                    # setting the values to zero.
                    if set_neg_results_to_zero:
                        predicts[predicts < 0] = 0

                    # Scale the output. All results are produced in 1.0e+3 lbs, so
                    # they need to be converted to kg.
                    rom_output[:, ind] = predicts * (
                        y_train_max_values[ind] * 1.0e3 * LBS_TO_KG)

        else:
            warning_msg = ''.join([
                'The time range did not include any years between {} and {}. ',
                'HydrocarbonLeakage components only produce output for years ',
                'in this range. The HydrocarbonLeakage output will be NaN ',
                'values.']).format(MIN_TIME, MAX_TIME)
            logging.debug(warning_msg)
            no_results = True

        if no_results:
            self.outputs[:, :] = np.ones((
                input_array.shape[0], input_array.shape[1])) * np.NaN
        else:
            for ind in range(len(y_train_max_values)):
                if not output_check[ind]:
                    self.outputs[:, ind] = (np.ones(
                        (input_array.shape[0], 1))
                        * VALUE_FOR_TIMES_OUTSIDE_RANGE).reshape(-1)
                else:
                    for time_ind in range(input_array.shape[0]):
                        if not (MIN_TIME <= input_array[time_ind, 0] <= MAX_TIME):
                            self.outputs[time_ind, ind] = \
                                VALUE_FOR_TIMES_OUTSIDE_RANGE
                        else:
                            self.outputs[time_ind, ind] = rom_output[time_ind, ind]

        return self.outputs

def set_up_input(time_array, input_dict):
    """
    This function sets up the input from input_array. The reservoirDepth is
    converted from m to ft, and both logResPerm and logWellPerm are converted
    from log10 m^2 to mD (milliDarcies).
    """
    days_per_year = 365.25

    ft_per_m = (1 / 0.3048)

    mD_per_m2 = (1 / 9.869223e-13) * 1000

    input_values = np.zeros((len(time_array), len(input_parameter_names)))

    for time_ind, time_val in enumerate(time_array):
        input_row_temp = [
            time_val / days_per_year,
            input_dict['reservoirDepth'] * ft_per_m,
            input_dict['NTG'], ((10 ** (input_dict['logResPerm'])) * mD_per_m2),
            input_dict['reservoirPressureMult'],
            ((10 ** (input_dict['logWellPerm'])) * mD_per_m2),
            input_dict['avgWaterSaturation'], input_dict['FCO2'],
            input_dict['FC1'], input_dict['FC4'], input_dict['FC7Plus']]

        input_values[time_ind, :] = np.array(input_row_temp)

    return input_values

def scale_input(input_array):
    """
    Function that scales the input provided by the range used in training the
    ROM (x_training_min to x_training_max). The 11 columns in the x_training_min
    and x_training_max lists correspond to the order of parameters in
    input_parameter_names.
    """
    # The first column has the time_array. The times outside of the allowed range
    # receive NaN results, but they also have to be removed so they don't change
    # the scaling.
    for time_ind in range(input_array.shape[0]):
        if input_array[time_ind, 0] < MIN_TIME:
            input_array[time_ind, 0] = MIN_TIME
        elif input_array[time_ind, 0] > MAX_TIME:
            input_array[time_ind, 0] = MAX_TIME

    input_values_for_scaling = np.zeros((input_array.shape[0] + 2,
                                         input_array.shape[1]))
    input_values_for_scaling[2:None, :] = input_array

    # Insert rows with the min and max training values so that it scales properly
    input_values_for_scaling[0, :] = x_training_min
    input_values_for_scaling[1, :] = x_training_max

    input_values_for_scaling2 = np.copy(input_values_for_scaling)
    input_values_for_scaling[:, (0, 1, 3, 4, 5)] = MinMaxScaler().fit_transform(
        input_values_for_scaling2[:, (0, 1, 3, 4, 5)])

    # Get the scaled input while excluding the inserted rows
    input_values_scaled = np.zeros((input_array.shape[0], input_array.shape[1]))

    input_values_scaled = input_values_for_scaling[2:None, :]

    return input_values_scaled
