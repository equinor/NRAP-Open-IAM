#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions for setup and data checking.

Author: Ernest N. Lindner
Date: 08/18/2022

Module Name
    seal_intro

Contents (14)
    check_missing_values(seal_controls)
    set_output_config(seal_controls)
    define_let_limits(param_bounds)
    define_input_limits()
    check_input_parameters(seal_controls, param_bounds)
    check_input_functions(key, val)
    cross_check_input(seal_controls)
    check_grid_area(seal_controls)
    create_cells(seal_controls)
    create_layout(seal_controls, grid)
    ----
    input_coordinates(grid, param_bounds)
    convert_lognorm_terms(mean_val, std_dev)
    preliminaries(alone, yaml_file)
    finish_definitions(seal_controls, param_bounds)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import sys                          # For buffer flush
import math                         # For Tanh + other functions
import logging                      # For reporting errors

import seal_file as sfile           # For file operations
import seal_model as mod            # For Class definitions
import seal_upload as sup           # Input file data
import seal_config as scfg           # IO directory and file names

# Other constants
ECHO = False                        # Control to print for debugging

logging.basicConfig(format='\n    > %(levelname)s: %(message)s',
                    level=logging.WARNING)


def check_missing_values(seal_controls):
    """Check if any input parameter is blank (= None) or Insufficient.

    Parameters
    ----------
    seal_controls = (dict) seal control parameters

    Returns
    -------
    stats = (int) status

    Notes
    -----
    1. "logging" module still uses old string format - added statements
        to create string messages prior to logging.
    2. Error messages -> explicitly format string before <logging>
        as function still uses old string format!
    """
    # Define error flag.
    stats = 0

    # Check all values for None.
    for key, val in seal_controls.items():
        if val is None:
            txt_blank = f'Parameter <{key}> is "None"!'
            logging.error(txt_blank)
            stats = -1

    # Check all function values.
    for key, val in seal_controls.items():
        stats += check_input_functions(key, val)

    # Check array values for correct range.
    if seal_controls['time_points'] < 2:
        msg = 'Defined Number Of Time Points Is Incorrect! (Must Be >1)'
        logging.error(msg)
        stats += -1

    if seal_controls['num_cells'] < 1:
        msg = 'Defined Number Of Cells is Incorrect (Must Be >0)!'
        logging.error(msg)
        stats += -1

    # Report error if found. Missing data is Fatal!
    if stats < 0:
        sys.stdout.flush()

    return stats


def set_output_config(seal_controls):
    """Reset output directory for code config file.

    Parameters
    ----------
    seal_controls = (dict) seal control parameters

    Returns
    -------
    None
    """
    test = seal_controls['output_directory']

    if test and test.strip():
        scfg.OUTPUT_DIRECTORY = seal_controls['output_directory']

    # return None


def define_let_limits(param_bounds):
    """Define numeric limits/ranges of LET model parameters.

    Parameters
    ----------
    param_bounds = (dict) parameter limits with min/max indices

    Returns
    -------
    param_bounds = (dict) parameter limits with min/max indices - updated
    """
    # Two-phase model parameters for L-E-T model.
    param_bounds['l_wetting'] = [0.5, 5.0]
    param_bounds['e_wetting'] = [0.1, 30.0]
    param_bounds['t_wetting'] = [0.0, 3.0]
    param_bounds['l_nonwet'] = [0.5, 5.0]
    param_bounds['e_nonwet'] = [0.1, 30.0]
    param_bounds['t_nonwet'] = [0.0, 3.0]

    # Other parameters for two-phase L-E-T model.
    param_bounds['l_capillary'] = [0.01, 5.0]
    param_bounds['e_capillary'] = [0.01, 30.0]
    param_bounds['t_capillary'] = [0.01, 3.0]
    param_bounds['max_capillary'] = [1.0E+02, 2.0E+08]

    return param_bounds


def define_input_limits():
    """Define numeric limits/ranges of seal input parameters.

    Parameters
    ----------
    N/A

    Returns
    -------
    param_bounds = (dict) parameter limits with min/max indices
    """
    # Define dictionary of input bounds.
    param_bounds = {}

    # Define controls (sec).
    param_bounds['start_time'] = [0.0, 5.0E+3]
    param_bounds['end_time'] = [1.0E+1, 1.0E+5]
    param_bounds['time_points'] = [1, 100]
    param_bounds['realizations'] = [1, 10000]
    param_bounds['total_inject'] = [1.0E+0, 5.0E+6]

    # Define permeability limits (in microdarcys!).
    param_bounds['perm_mean'] = [1.0E-04, 1.0E+02]
    param_bounds['perm_std'] = [0.0, 1.0E+01]
    param_bounds['perm_min'] = [1.0E-06, 1.0E+01]
    param_bounds['perm_max'] = [1.0E-03, 1.0E+06]  # max = 100 darcys
    param_bounds['perm_heter_factor'] = [1.0E-02, 1.0E+02]

    # Define grid limits.
    param_bounds['grid_rows'] = [1, 100]
    param_bounds['grid_cols'] = [1, 100]
    param_bounds['cell_height'] = [1.0, 1.0E+04]
    param_bounds['cell_width'] = [1.0, 1.0E+04]

    # Define geometry limits (m/m2/Pa).
    param_bounds['num_cells'] = [1, 10000]
    param_bounds['area'] = [1.0, 2.6E+05]
    param_bounds['static_depth'] = [0.8E+02, 9.5E+03]
    param_bounds['static_pressure'] = [1.0E+06, 6.0E+07]

    # Define fluid parameter limits.
    param_bounds['temperature'] = [3.1E+01, 1.8E+02]
    param_bounds['salinity'] = [0.0, 8.0E+04]
    param_bounds['co2_density'] = [9.3E+01, 1.05E+03]
    param_bounds['co2_viscosity'] = [1.8E-05, 1.4E-04]
    param_bounds['brine_density'] = [8.8E+02, 1.08E+03]
    param_bounds['brine_viscosity'] = [1.5E-04, 1.6E-03]
    param_bounds['co2_solubility'] = [0.0, 2.0]
    param_bounds['ave_base_depth'] = [0.8E+02, 9.5E+03]
    param_bounds['ave_base_pressure'] = [1.0E+06, 6.0E+07]

    # Define thickness limits.
    param_bounds['thickness_ave'] = [10.0, 1000.0]
    param_bounds['thickness_std'] = [0.00, 500.0]
    param_bounds['thickness_min'] = [5.0, 600.0]
    param_bounds['thickness_max'] = [10.0, 1000.0]

    # Define relative permeability limits.
    param_bounds['resid_brine'] = [0.01, 0.35]
    param_bounds['resid_co2'] = [0.00, 0.35]
    param_bounds['perm_ratio'] = [0.00, 1.5]
    param_bounds['entry_pressure'] = [1.0E+02, 2.0E+06]

    # Define two-phase model parameters for BC model.
    param_bounds['zeta'] = [0.00, 5.0]

    # Define LET Parameters.
    param_bounds = define_let_limits(param_bounds)

    # Define Time-model parameters.
    param_bounds['influence'] = [0.00, 1.00]
    param_bounds['total_effect'] = [0.01, 200.0]
    param_bounds['rate_effect'] = [0.01, 0.65]
    param_bounds['reactivity'] = [0.0, 10.0]
    param_bounds['clay_content'] = [0.0, 100.0]
    param_bounds['carbonate_content'] = [0.0, 100.0]

    # Define repository input limits.
    param_bounds['reservoir_pressure'] = [1.0E+05, 1.0E+08]
    param_bounds['reservoir_saturation'] = [0.0, 1.0]
    param_bounds['coordinates'] = [-1.0E+05, 1.0E+05]

    # Define axis control of plot.
    param_bounds['max_draw_time'] = [1.0, 1.0E+04]

    return param_bounds


def check_input_parameters(seal_controls, param_bounds):
    """Check whether input parameters fall within specified limits.

    Parameters
    ----------
    seal_controls = (dict) seal control parameters
    param_bounds = (dict) bounds on seal parameters (w. 2 values)

    Returns
    -------
    statx = (int) error code, default = 0, if warning/error = -1

    Notes
    -----
    1. "logging" module still uses string format - added statements
        to create string messages prior to logging.
    2. Some parameters are checked using a defined list.
    """
    # Define parameter for checking.
    statx = 0

    # Combine True/False parameters into single list.
    true_false_parameters = ['approach', 'initialize', 'plot',
                             'time_input', 'choice', 'use']
    truth_list = [True, False]

    # Combine function types into a single list.
    function_parameters = ['clay_type', 'model', 'relative_model', 'title',
                           'version', 'output_directory']

    # Check if values are within defined range.
    for key, val in seal_controls.items():

        # Check if code_start value is positive.
        if 'code_start' in key:
            if val <= 0:
                txt_wrong = (f'Parameter <{key}> is Not Defined Correctly '
                             f'and Has a Value of "{val}"!')
                logging.error(txt_wrong)
                statx += -1

        # Check items in functions list have some input.
        elif key in function_parameters:
            if val == '':
                txt_wrong = f'Parameter <{key}> is Not Defined!'
                logging.error(txt_wrong)
                statx += -1

        # Check True/False controls.
        elif [ele for ele in true_false_parameters if ele in key]:
            if val not in truth_list:
                txt_wrong = (f'Parameter <{key}> is Not Defined Correctly '
                             f'and Has a Value of "{val}"!')
                logging.error(txt_wrong)
                statx += -1

        # Check numeric values to be within bounds.
        elif key in param_bounds:
            if ((val < param_bounds[key][0])
                    or (val > param_bounds[key][1])):
                txt_outside = (f'Parameter <{key}> is Outside of Recommended '
                               f'Bounds with a Value of {val}!')
                logging.error(txt_outside)
                statx += -1

        else:
            # Missing parameter in seal_controls.
            txt_not_seal = (f'Parameter <{key}> Not Recognized as a '
                            'Seal Input Parameter.')
            logging.error(txt_not_seal)
            statx += -1

    return statx


def check_input_functions(key, val):
    """Check whether input parameters have specified function values.

    Parameters
    ----------
    key = (str) input parameter title from dictionary
    val = (float) item in parameter list
    error_msg = (str) error message to output if error

    Returns
    -------
    function_statx = error code, if warning it is set to (-1)
    """
    # Define parameter lists for checking.
    clay_list = ["smectite", "illite", "chlorite"]
    model_list = [0, 1, 2]
    relative_list = ["BC", "LET"]
    function_statx = 0

    # Clay types.
    if key == 'clay_type':
        if val not in clay_list:
            msg = (f'Parameter <{key}> is Not Defined Correctly '
                   f'and Has a Value of "{val}"!')
            logging.error(msg)
            function_statx += -1

    # Model types.
    elif key == 'model':
        if val not in model_list:
            msg = (f'Parameter <{key}> is Not Defined Correctly '
                   f'and Has a Value of "{val}"!')
            logging.error(msg)
            function_statx += -1

    # Relative permeability types.
    else:
        if key == 'relative_model':
            if val not in relative_list:
                msg = (f'Parameter <{key}> is Not Defined Correctly '
                       f'and Has a Value of "{val}"!')
                logging.error(msg)
                function_statx += -1

    sys.stdout.flush()
    return function_statx


def cross_check_input(seal_controls):
    """Check if input is within bounds for permeability and thickness.

    Parameters
    ----------
    seal_seal_controls = (dict) seal control parameters

    Returns
    -------
    statx = (int) error code, default = 0, if warning/error = -1

    Notes
    -----
    1. "logging" module still uses string format -
        so code creates string messages prior to calling <logging>.
    2. Flush added as current environment not pushing all comments.
    """
    # Set default return.
    statx = 0

    # Run check only if no file input for permeability.
    if not seal_controls['perm_input_approach']:
        # Run checks on order of min and max permeability values.
        if seal_controls['perm_min'] > seal_controls['perm_max']:
            causative = "permeability"
            txt_cross = f'Minimum {causative} Exceeds Maximum {causative}!'
            logging.error(txt_cross)
            statx += -1
            sys.stdout.flush()

        # Check mean to be inside bounds.
        if (seal_controls['perm_mean'] > seal_controls['perm_max'] or
                seal_controls['perm_mean'] < seal_controls['perm_min']):
            causative = "permeability"
            txt_outside = f'Mean {causative} is Outside Min./Max. Bounds!'
            logging.error(txt_outside)
            statx += -1
            sys.stdout.flush()

    # Run check only if no file input for thickness.
    if not seal_controls['thickness_approach']:
        # Run checks on thickness values.
        if seal_controls['thickness_min'] > seal_controls['thickness_max']:
            causative = "thickness"
            txt_cross = f'Minimum {causative} Exceeds Maximum {causative}!'
            logging.error(txt_cross)
            statx += -1
            sys.stdout.flush()

        # Check average thickness against bounds.
        condition_a = (seal_controls['thickness_ave']
                       < seal_controls['thickness_min'])
        condition_b = (seal_controls['thickness_ave']
                       > seal_controls['thickness_max'])
        if condition_a or condition_b:
            causative = "thickness"
            txt_outside = f'Mean {causative} is Outside Min./Max. Bounds!'
            logging.error(txt_outside)
            statx += -1
            sys.stdout.flush()

    return statx


def check_grid_area(seal_controls):
    """Check if grid area-related input is consistent.

    Parameters
    ----------
    seal_controls = (dict) seal control parameters

    Returns
    -------
    statx = (int) error code, default = 0, if warning or error = -1

    """
    # Set default return (statx = -1 will terminate program).
    statx = 0

    # Run check only if no file input for permeability.
    if seal_controls['grid_approach']:

        # Provide simple definitions.
        high = seal_controls['cell_height']
        wide = seal_controls['cell_width']
        box = high * wide
        input_area = seal_controls['area']

        # Check area - log warning.
        if (box < 0.99 * input_area) or (box > 1.01 * input_area):
            run_txt = ('Yaml Input Error: '
                       f'Product of Cell Dimensions {high} and {wide} '
                       f'Is Not Equal to Defined Area of {input_area}!')
            logging.error(run_txt)
            statx = -1
            sys.stdout.flush()

    return statx


def create_cells(seal_controls):
    """Create grid (a list of cells) with default coordinates.

    Parameters
    ----------
    seal_controls = (dict) seal control parameters

    Returns
    -------
    grid = (list of class) - empty class cells
    """
    # Define variables for clarity - set all to 0.0.
    total = seal_controls['num_cells']
    grid = []

    # Create a grid of uniform cells with default x/y values = (0,0).
    for _ in range(total):
        grid.append(mod.Cell())

    return grid


def create_layout(seal_controls, grid):
    """Create grid (list of cells) with box coordinates on simple grid.

    Parameters
    ----------
    seal_controls = (dict) seal control parameters
    grid = (list of class) cells (empty)

    Returns
    -------
    grid = (list of class) updated cells - center coordinates

    Notes
    -----
    For simple box setup, coordinates assume tight rectangular arrangement
        and centers are at 1/2 of the relevant x/y dimension.
    """
    # Define variables for clarity.
    grid_y = seal_controls['grid_rows']
    grid_x = seal_controls['grid_cols']
    high = seal_controls['cell_height']
    wide = seal_controls['cell_width']

    # Create a grid of uniform-spaced cells with defined center offsets.
    side_shift = (wide * 0.5)
    up_shift = (high / 2.0)
    for i in range(grid_y):
        for j in range(grid_x):
            indx = j + (i * grid_x)
            y_center = (i * high) + up_shift
            x_center = (j * wide) + side_shift
            grid[indx].set_coord(x_center, y_center)

    return grid


def input_coordinates(seal_controls, grid, param_bounds):
    """Input or define x and y coordinate data for cells, as desired.

    Parameters
    ----------
    seal_controls = (dict) seal control parameters
    grid = (list of class) cells
    param_bounds = (dict) seal parameter bounds

    Returns
    -------
    grid = (list of class) cells
    """
    # Check if input is desired.
    if seal_controls['layout_approach']:
        grid = sup.define_coordinates(grid, param_bounds)
    else:
        # Assign general units if rectangular grid is assumed.
        if seal_controls['grid_approach']:
            grid = create_layout(seal_controls, grid)

    return grid


def convert_lognorm_terms(mean_val, std_dev):
    """Convert mean/std_dev into log-normal parameters for NumPy.

    Parameters
    ----------
    mean_val = (float) mean permeability in microdarcys
    std_dev = (float) std deviation of permeability

    Returns
    -------
    log_mu = (float) log-normal mean
    sigma = (float) log-normal variance

    Notes
    -----
    see https://en.wikipedia.org/wiki/Log-normal_distribution
    """
    # define common terms.
    variance = math.pow(std_dev, 2)
    mean_2 = math.pow(mean_val, 2)
    inside = 1.0 + (variance / mean_2)

    # Compute location - mean.
    term = math.sqrt(inside)
    log_mu = math.log(mean_val / term)

    # Compute sigma - variance.
    term = math.log(inside)
    sigma = math.sqrt(term)

    return log_mu, sigma


def preliminaries(alone, seal_controls):
    """Open control file and check data.

    Parameters
    ----------
    alone = (bool) status of code operation; True = if stand-alone
    seal_controls = (dict) seal control parameters

    Returns
    -------
    seal_controls = (dict) seal control parameters
    param_bounds = (dict) seal parameter bounds
    """
    # Check for any missing data.
    stats = check_missing_values(seal_controls)
    if stats < 0:
        sfile.opx_problem("Input Parameter(s) Outside of Range or Wrong!")

    # Identify input bounds and then check values.
    # --> param_bounds is a dict()
    param_bounds = define_input_limits()
    stats = check_input_parameters(seal_controls, param_bounds)
    stats += cross_check_input(seal_controls)
    stats += check_grid_area(seal_controls)
    if stats < 0:
        sfile.opx_problem("Input Parameter(s) Outside of Range or Wrong!")

    # Set output directory in config file.
    set_output_config(seal_controls)

    # Wipe output directory to remove old output files.
    if alone:
        sfile.echo_status(alone,
                          "DELETING OLD FILES & PERFORMING INPUT CHECK.")
    sfile.clean_output()

    # Check for inconsistent input with grid_approach.
    if seal_controls['area_approach']:
        seal_controls['grid_approach'] = False

    return seal_controls, param_bounds


def finish_definitions(seal_controls, param_bounds):
    """Finish input operations by defining values.

    Parameters
    ----------
    seal_controls = (dict) seal control parameters
    param_bounds = (dict) parameter bounds (dictionary)

    Returns
    -------
    grid = (list of class) cells
    num_cells = (int) = number of cells
    seal_controls = (dict) seal control parameters (revised)
    """
    # Create a grid of cells w default coordinates.
    grid = create_cells(seal_controls)

    # Define coordinates of grid from file if desired.
    grid = input_coordinates(seal_controls, grid, param_bounds)

    # Assign class variables to Cell Class from control file.
    mod.Cell.assign_controls(seal_controls)

    # Define number of lines & columns for reservoir input.
    seal_controls['r_lines'], seal_controls['r_columns'] = \
        sup.reservoir_check(seal_controls)

    # Define log-normal parameters to compute total permeabilities for grid.
    if not seal_controls['perm_input_approach']:
        seal_controls['perm_location'], seal_controls['perm_scale'] = \
            convert_lognorm_terms(seal_controls['perm_mean'],
                                  seal_controls['perm_std'])
    else:
        seal_controls['perm_location'] = 0.0
        seal_controls['perm_scale'] = 0.0

    # Define time-dependent permeability model and create list for results.
    mod.Cell.model = seal_controls['model']

    # Set control for reading permeabilities from file.
    #  - Note: In reading permeabilities from file, only one set of data is
    #          expected; therefore, the following control limits reading of
    #          input file only once.
    seal_controls['read_perm'] = True

    return grid, seal_controls


#
# -----------------------------------------------------------------------------
# - End of module
