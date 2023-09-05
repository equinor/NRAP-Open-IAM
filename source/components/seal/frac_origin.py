#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions for fracture analysis start, input from YAML file and checking.

Author: Ernest N. Lindner
Date: 08/19/2022

Module Name
    frac_origin

Contents (18)
    create_perm_storage(frac_controls)
    create_entry_storage(frac_controls)
    transfer_dictionary_values(seal_controls, frac_controls)
    define_param_limits()
    check_frac_params(frac_controls, params_bounds)
    check_logic(parameter_name, ave_valu, min_valu, max_valu)
    cross_check_input(frac_controls)
    check_pressure(frac_controls)
    check_blanks(seal_controls)
    setup_aperture_lognormal(frac_controls)
    ----
    setup_matrix_lognormal(frac_controls)
    setup_threshold(frac_controls)
    setup_kappa(frac_controls)
    frac_launch(alone, yaml_file, seal_controls)
    assign_values(perm_sum, threshold_sum, grid)
    compute_frac_time(flag_time, frac_controls)
    obtain_fracture_data(file_name):
    evaluate_fracs(grid, frac_controls, sim_step, alive)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import math                         # For pow
import logging                      # For reporting errors
import time                         # For calculating runtime
from datetime import timedelta      # For calculating runtime
import numpy as np                  # Array operations

import frac_random as fran          # Random number generators
import frac_stochastic as fsch      # Computations for random fractures
import frac_user as fuse            # Matrix and fracture computations
import frac_view as view            # Save fracture summary data
import seal_config as scfg          # IO directory names
import seal_file as sfile           # Error in file operations
import seal_units as sun            # For unit conversion
import frac_decipher as fde         # YAML parameters

# Constants
DEG_TO_KAPPA = 0.00875              # Conversion for 2-sigma deg. to Kappa
HOBO = True                         # Debug flag - print fracture results

logging.basicConfig(format='\n    > %(levelname)s: %(message)s',
                    level=logging.WARNING)


def create_perm_storage(frac_controls):
    """Create NumPy arrays for permeability loops.

    Parameters
    ----------
    frac_controls = (dict) fracture parameters

    Returns
    -------
    random_perm = (array) permeabilities for stochastic fractures
    user_threshold = (array) threshold pressure for user fractures
    matrix_perm = (array) permeabilities for matrix
    """
    # Create 1D storage arrays for a realization.
    num_cells = frac_controls['num_cells']

    # Setup/zero-out permeability arrays.
    random_perm = np.zeros(num_cells)
    user_perm = np.zeros(num_cells)
    matrix_perm = np.zeros(num_cells)

    return random_perm, user_perm, matrix_perm


def create_entry_storage(frac_controls):
    """Create NumPy arrays for threshold pressure loops.

    Parameters
    ----------
    frac_controls = (dict) fracture parameters

    Returns
    -------
    random_threshold = (array) entry pressures for stochastic fractures
    user_threshold = (array) entry pressures for user fractures
    matrix_threshold = (array) entry pressures for matrix
    """
    # Create 1D storage arrays for a realization.
    num_cells = frac_controls['num_cells']

    # Setup/zero-out permeability arrays.
    random_threshold = np.zeros(num_cells)
    user_threshold = np.zeros(num_cells)
    matrix_threshold = np.zeros(num_cells)

    return random_threshold, user_threshold, matrix_threshold


def transfer_dictionary_values(seal_controls, frac_controls):
    """Convert seal_controls to frac_controls as needed.

    Parameters
    ----------
    seal_controls = (dict) seal control parameters
    frac_controls = (dict) frac_control parameters

    Returns
    -------
    1. frac_controls = (dict) updated frac_control parameters
    2. These parameters are not checked - so routine is after param check.

    """
    # Get run params needed for fractures.
    frac_controls['title'] = seal_controls['title']
    frac_controls['simulations'] = seal_controls['realizations']
    frac_controls['num_cells'] = seal_controls['num_cells']
    frac_controls['grid_approach'] = seal_controls['grid_approach']
    if frac_controls['grid_approach']:
        frac_controls['cell_height'] = seal_controls['cell_height']
        frac_controls['cell_width'] = seal_controls['cell_width']
        frac_controls['static_depth'] = seal_controls['static_depth']
        frac_controls['entry_pressure'] = seal_controls['entry_pressure']
    frac_controls['ave_base_depth'] = seal_controls['ave_base_depth']
    frac_controls['ave_base_pressure'] = seal_controls['ave_base_pressure']
    frac_controls['version'] = seal_controls['version']

    return frac_controls


def define_param_limits():
    """Define numeric limits/ranges of input parameters.

    Parameters
    ----------
    N/A

    Returns
    -------
    params_bounds = (dict) dictionary with min/max indices
    """
    # Define dictionary of input limits.
    params_bounds = {}

    # Run Params:
    params_bounds['simulations'] = [1, 200]

    # Connectivity:
    params_bounds['connect_factor'] = [0.0, 1.01]

    # Random Fracs - Density (per m2):
    params_bounds['density_ave'] = [1.0E-10, 1.0E+02]
    params_bounds['density_min'] = [1.0E-12, 5.0E+01]
    params_bounds['density_max'] = [1.0E-10, 5.0E+02]

    # Random Fracs - Orientation (deg):
    params_bounds['orient_mu'] = [-180.0, 180.0]
    params_bounds['orient_sigma'] = [0.0, 180.0]

    # Random Fracs - Length (m):
    params_bounds['length_ave'] = [0.1, 1.0E+02]
    params_bounds['length_dev'] = [0.001, 5.0E+01]
    params_bounds['length_eta'] = [-5.00, 5.00]
    params_bounds['length_min'] = [0.1, 50.0]
    params_bounds['length_max'] = [1.0, 5.0E+03]

    # Random Fracs - Vary Aperture (mm):
    params_bounds['aperture_ave'] = [1.0E-04, 5.0E+01]
    params_bounds['aperture_dev'] = [0.0, 1.0E+01]
    params_bounds['aperture_min'] = [0.0, 1.0E+01]
    params_bounds['aperture_max'] = [1.0E-03, 1.0E+02]

    # Random Fracs - Correlate Aperture:
    params_bounds['aperture_alpha'] = [0.4, 1.0]
    params_bounds['aperture_beta'] = [0.00001, 1.0]

    # Random Threshold (Pa):
    params_bounds['entry_pressure'] = [1.0, 5.0E+06]

    # Pressure-Aperture Correction:
    params_bounds['residual_aperture'] = [0.0001, 0.1]
    params_bounds['wide_aperture'] = [0.001, 10.0]
    params_bounds['stress_limit'] = [1.0E+06, 5.0E+07]
    params_bounds['theta_aperture'] = [0.1, 10.0]

    # Reference aperture:
    params_bounds['ref_user_aperture'] = [1.0E-04, 5.0E+01]
    params_bounds['ref_user_pressure'] = [1.0, 5.0E+06]

    # Matrix Permeability:
    params_bounds['rock_perm_ave'] = [1.0E-08, 1.0E+04]
    params_bounds['rock_perm_dev'] = [0.0, 1.0E+02]
    params_bounds['rock_perm_min'] = [1.0E-10, 1.0E+02]
    params_bounds['rock_perm_max'] = [1.0E-06, 1.0E+03]

    # Matrix Reference:
    params_bounds['ref_matrix_perm'] = [1.0E-10, 1.0E+04]
    params_bounds['ref_matrix_threshold'] = [1.0, 5.0E+07]

    return params_bounds


def check_frac_params(frac_controls, params_bounds):
    """Check whether input parameters fall within specified limits.

    Parameters
    ----------
    frac_controls = (dict) input parameters of seal component model
    params_bounds = (dict) limits on seal parameters (w. 2 values)

    Returns
    -------
    statx = (int) error code, default = 0, if warning = -1

    Notes
    -----
    1. "logging" module still uses old string format - added statements
            to construct string messages prior to logging.
    2. Some parameters are checked using a defined list.
    """
    # Define error control flag.
    statx = 0

    # Combine control strings into single list.
    command_parameters = ['user_input_file']

    # Combine True/False parameters into single list.
    truth_list = [True, False]
    true_false_parameters = ['pressure_approach',
                             'random_approach', 'user_approach',
                             'plot_fractures', 'correlate_approach']

    # Combine other parameters into single list.
    method_parameters = ['length_approach']
    method_list = ['LOGNORM', 'POWER']

    # Check all items if values are OK or within defined range.
    for key, val in frac_controls.items():
        # Check items in functions list have some input.
        if key in command_parameters:
            if val == '':
                txt_wrong = f'Parameter <{key}> is Not Defined!'
                logging.error(txt_wrong)
                statx += -1

        # Check items in method list.
        elif [ele for ele in method_parameters if ele in key]:
            if val not in method_list:
                txt_wrong = (f'Parameter <{key}> is Not Defined Correctly '
                             + 'and Has a Value of "{val}"!')
                logging.error(txt_wrong)
                statx += -1

        # Check True/False parameters.
        elif [ele for ele in true_false_parameters if ele in key]:
            if val not in truth_list:
                txt_wrong = (f'Parameter <{key}> is Not True/False '
                             + 'and Has a Value of "{val}"!')
                logging.error(txt_wrong)
                statx += -1

        # Check numeric values to be within limits.
        elif key in params_bounds:
            if ((val < params_bounds[key][0])
                    or (val > params_bounds[key][1])):
                txt_outside = (f'\n   -->> Parameter <{key}> is outside of '
                               + 'the recommended limits with a value of '
                               + f'{val}!')
                logging.error(txt_outside)
                statx += -1

        else:
            # Missing parameters in frac_controls.
            txt_not = (f'\n   -->> Parameter <{key}> not recognized as a '
                       + 'Seal frac_controls input parameter.')
            logging.error(txt_not)
            statx += -1

    return statx


def check_logic(parameter_name, ave_valu, min_valu, max_valu):
    """Check if stochastic input is set up correctly.

    Parameters
    ----------
    parameter_name = (str) name of value being checked
    ave_valu = (float) average value
    min_valu = (float) minimum value
    max_valu = (float) maximum value

    Returns
    -------
    error_flag = error code, default = 0, if warning = -1

    Notes
    -----
    1. "logging" module still uses string format;
        create string messages prior to calling <logging>.
    """
    # Set error flag as OK.
    error_flag = 0

    # Check relationships of data.
    if min_valu > max_valu:
        txt_cross = f'\n   -->> Minimum {min_valu} exceeds maximum of '
        txt_cross += 'f{max_valu} for + parameter_name'
        logging.error(txt_cross)
        error_flag = -1

    if ave_valu > max_valu or ave_valu < min_valu:
        txt_outside = f'\n   -->> Mean {parameter_name} is outside bounds! '
        logging.error(txt_outside)
        error_flag = -1

    return error_flag


def cross_check_input(frac_controls):
    """Check if input for logic flaws.

    Parameters
    ----------
    frac_controls = (dict) input parameters of component

    Returns
    -------
    statx = (int) error code, default = 0, if warning = -1

    Notes
    -----
    1. "logging" module still uses string format -
        create string messages prior to calling <logging>.
    """
    # Run checks on statistic inputs - ensure mean is inside bounds!
    # -------------------------------------------------------------------------
    # Run checks if random fractures are generated
    statx = 0
    if frac_controls['random_approach']:
        # Check Density.
        statx = check_logic("density", frac_controls['density_ave'],
                            frac_controls['density_min'],
                            frac_controls['density_max'])

        # Check aperture.
        statx += check_logic("aperture", frac_controls['aperture_ave'],
                             frac_controls['aperture_min'],
                             frac_controls['aperture_max'])

        # Check length - only if using a lognormal approach.
        if "LOGNORM" in frac_controls['length_approach']:
            statx += check_logic("length", frac_controls['length_ave'],
                                 frac_controls['length_min'],
                                 frac_controls['length_max'])
        else:
            min_valu = frac_controls['length_min']
            max_valu = frac_controls['length_max']
            if min_valu > max_valu:
                msg = (f'\n   -->> Minimum {min_valu} exceeds maximum '
                       f'{max_valu}! for length!')
                logging.error(msg)
                statx += -1

    return statx


def check_pressure(frac_controls):
    """Check if input for pressure is OK.

    Parameters
    ----------
    frac_controls = (dict) input parameters of component

    Returns
    -------
    statx = (int) error code, default = 0, if warning = -1

    Notes
    -----
    Calculation is in Pa.

    """
    # Get minimum lithostatic pressure (in Pa) & input pressure.
    depth = frac_controls['ave_base_depth']
    min_pressure = sun.brine_pressure(depth)
    input_pressure = frac_controls['ave_base_pressure']

    # Check pressure and report if in error.
    if input_pressure < min_pressure:
        # Error - message.
        msg = (f'\n   -->> Ave Base Pressure of {input_pressure} is less than '
               + f' the Minimum from Gradient of {min_pressure}! ')
        logging.error(msg)
        statx = -1
    else:
        statx = 0

    return statx


def check_blanks(frac_controls):
    """Check whether if any input parameter is blank (= "None"").

    Parameters
    ----------
    frac_controls = (dict) input parameters of component model

    Returns
    -------
    None

    Notes
    -----
    1."logging" module still uses string format - added statements
        to create string messages prior to logging.
    """
    # Define error flag.
    err_flag = 0

    # Check all values for "None".
    for key, val in frac_controls.items():
        if val is None:
            txt_blank = f'\n   -->> Parameter <{key}> is Not Defined <None>!'
            logging.error(txt_blank)
            err_flag = -1

    # Report error if found. Missing data is Fatal!
    if err_flag == -1:
        sfile.opx_problem("Input Parameter(s) are Missing!")

    # return None


def setup_aperture_lognormal(frac_controls):
    """Define log-normal aperture input variables for distribution.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters

    Returns
    -------
    frac_controls = (dict) with new frac_controls parameters
    """
    # Convert aperture parameters for log-normal distribution.
    frac_controls['aperture_mu'], frac_controls['aperture_scale'] = \
        fran.convert_lognorm_terms(frac_controls['aperture_ave'],
                                   frac_controls['aperture_dev'])

    # Convert length parameters for log-normal distribution.
    if "LOGNORM" in frac_controls['length_approach']:
        frac_controls['length_mu'], frac_controls['length_scale'] = \
            fran.convert_lognorm_terms(frac_controls['length_ave'],
                                       frac_controls['length_dev'])

    return frac_controls


def setup_matrix_lognormal(frac_controls):
    """Define log-normal matrix input variables for distribution.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters

    Returns
    -------
    frac_controls = (dict) with new frac_controls parameters
    """
    # Convert matrix parameters for log-normal distribution.
    frac_controls['rock_perm_mu'], frac_controls['rock_perm_scale'] = \
        fran.convert_lognorm_terms(frac_controls['rock_perm_ave'],
                                   frac_controls['rock_perm_dev'])

    return frac_controls


def setup_threshold(frac_controls):
    """Define reference threshold pressure factor for fractures.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters

    Returns
    -------
    frac_controls = (dict) parameters - adjusted
    """
    # Define variables depending on approach.
    if frac_controls['random_approach']:
        actual_b = frac_controls['aperture_ave']
        pressure = frac_controls['entry_pressure']
        frac_controls['threshold_factor'] = pressure * actual_b

    elif frac_controls['user_approach']:
        actual_b = frac_controls['ref_user_aperture']
        pressure = frac_controls['ref_user_pressure']
        frac_controls['threshold_factor'] = pressure * actual_b

    return frac_controls


def setup_kappa(frac_controls):
    """Define kappa for von Mises orientation distribution.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters

    Returns
    -------
    frac_controls = (dict) parameters - adjusted

    Notes
    -----
    1. sigma is the spread in degrees - must convert to kappa;
        kappa is "concentration" -> larger k is a smaller spread.
    2. Rough empirical fit for kappa = (sigma*DEG_TO_KAPPA)^-2;
    3. sigma should be > 0.1 for numerical stability of kappa.
    """
    sigma = frac_controls['orient_sigma']
    if sigma < 0.1:
        sigma = 0.1 * DEG_TO_KAPPA
    else:
        sigma *= DEG_TO_KAPPA
    frac_controls['orient_disperse'] = math.pow(sigma, -2.0)

    return frac_controls


def frac_launch(alone, yaml_file, seal_controls):
    """Open input control file and get/check data.

    Parameters
    ----------
    alone = (bool) status of code operation; True = if stand-alone
    yaml_file = (str) name of input file
    seal_controls = (dict) seal control parameters
    frac_controls = (dict) fracture control parameters

    Returns
    -------
    frac_controls = (dict) fracture control parameters - new
    """
    frac_controls = {}
    if seal_controls['fracture_approach']:
        # Input control parameters from text YAML file and seal_controls.
        sfile.echo_status(alone, "READING FRACTURE CONTROL FILE.")
        frac_controls = fde.acquire_frac_parameters(yaml_file, frac_controls)

        # Check for blanks in data.
        check_blanks(frac_controls)

        # Identify input limits for data.
        param_bounds = define_param_limits()

        # Check control values.
        status = check_frac_params(frac_controls, param_bounds)
        status += cross_check_input(frac_controls)

        # Add prior seal parameters to frac list and final check.
        frac_controls = transfer_dictionary_values(seal_controls,
                                                   frac_controls)
        status += check_pressure(frac_controls)

        # Consider all warnings fatal! -> exit program if error found.
        if status < 0:
            sfile.opx_problem("Input Parameter(s) Outside of Range or Wrong!")

        # Add matrix parameters for threshold & lognormal
        frac_controls = setup_threshold(frac_controls)
        frac_controls = setup_matrix_lognormal(frac_controls)

        # Define other parameters only if random approach.
        if frac_controls['random_approach']:

            # Setup log-normal parameters.
            frac_controls = setup_aperture_lognormal(frac_controls)

            # Convert sigma to kappa for von Mises distribution.
            frac_controls = setup_kappa(frac_controls)

            # Define orientation of second fracture set.
            frac_controls['orient_set_2'] = frac_controls['orient_mu'] + 90.0
            if frac_controls['orient_set_2'] > 360.0:
                frac_controls['orient_set_2'] -= 360.0

    else:
        # No fracture generation or input - set dictionary to error code.
        frac_controls['contents'] = -999.999

    return frac_controls


def assign_values(perm_sum, threshold_sum, grid):
    """Check whether if any input parameter is blank (= "None"").

    Parameters
    ----------
    perm_sum = (array) permeability results
    threshold_sum = (array) entry pressure results
    grid = (list of class) collection of cells

    Returns
    -------
    grid = (list of class) updated values
    """
    # Adjust grid for results.
    for indx, cell in enumerate(grid):
        cell.permeability = perm_sum[indx]
        cell.entry = threshold_sum[indx]

    return grid


def compute_frac_time(flag_time, frac_controls):
    """Obtain seal parameters from a formatted YAML file.

    Parameters
    ----------
    flag_time = (time) start time of computations
    frac_controls = (dict) dictionary of fracture controls (empty)

    Returns
    -------
    frac_controls = (dict) dictionary of fracture parameters - with new values
    """
    checker_flag = time.monotonic()
    elapsed_time = timedelta(seconds=(checker_flag - flag_time))
    frac_controls['elapsed'] = elapsed_time

    return frac_controls


def obtain_fracture_data(file_name):
    """Input fracture data from CSV file as list.

    Parameters
    ----------
    file_name = (str) name of user input file (w.o. destination directory)

    Returns
    -------
    user_list = (list) list of fracture data

    Notes
    -----
    1. Used only for fracture data input.
    """
    # Construct path name to file.
    subdirectory_path, destination = sfile.get_path_name(scfg.INPUT_DIRECTORY,
                                                         file_name,
                                                         scfg.EXTENSION_CSV)
    # Get data from csv file.
    user_list = sfile.acquire_csv_list(destination, subdirectory_path,
                                       file_name)

    return user_list


def evaluate_fracs(grid, frac_controls, sim_numbr, alive):
    """Evaluate the effect of fracturing effect on permeability & pressure.

    Parameters
    ----------
    grid = (list) (class) a collection of cells
    frac_controls = (dict) fracture control parameters
    sim_numbr = (int) current simulation number
    alive = (bool) stand-alone operations flag

    Returns
    -------
    grid = (list) updated cell values
    """
    # -------------------------------------------------------------------------
    # A. ARRAY SETUP
    # -------------------------------------------------------------------------
    # Start clock for frac run time.
    flag_time = time.monotonic()

    # Setup/zero-out NumPy permeability arrays.
    random_perm, user_perm, matrix_perm = \
        create_perm_storage(frac_controls)

    # Setup/zero-out NumPy threshold pressure arrays.
    random_threshold, user_threshold, matrix_threshold = \
        create_entry_storage(frac_controls)

    # Setup NumPy arrays for combined output.
    combo_list = [[], [], []]

    # -------------------------------------------------------------------------
    # B. COMPUTATION LOOP
    # -------------------------------------------------------------------------
    # >> Case 1. **** Evaluate random fractures. ****
    if frac_controls['random_approach']:
        # Generate random fractures and calculate permeabilities for grid.
        random_perm, random_threshold, combo_list[0] = \
            fsch.compute_random_permeability(frac_controls, grid, alive)

    # >> Case 2. **** Evaluate user fractures. ****
    if frac_controls['user_approach']:
        # Input user fractures.
        user_data = obtain_fracture_data(frac_controls['user_input_file'])

        # Process user fractures to obtain permeability.
        user_perm, user_threshold, combo_list[1] = \
            fuse.process_user_lines(frac_controls, grid, user_data)

    # >> Case 3. **** Evaluate matrix characteristics ****
    # Use probability distribution to get matrix values.
    matrix_perm, matrix_threshold, combo_list[2] = \
        fuse.define_matrix_perm(frac_controls)

    # >> 4. **** Compute final equivalent properties. ****
    # Sum/process case arrays to get equivalent values for grid.
    perm_sum = random_perm + user_perm + matrix_perm
    entry_sum = fuse.sum_threshold_arrays(frac_controls, random_threshold,
                                          user_threshold, matrix_threshold)

    # Compute elapsed time of run and save.
    frac_controls = compute_frac_time(flag_time, frac_controls)

    # Debugging option printout.
    if HOBO:
        view.write_data_results(frac_controls, combo_list, sim_numbr)

    # Set grid values and return.
    grid = assign_values(perm_sum, entry_sum, grid)

    return grid


#
# -----------------------------------------------------------------------------
# - End of module
