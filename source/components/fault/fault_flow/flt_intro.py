#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions for checking and interpreting YAML file.

Author: Ernest N. Lindner
Date: 08/16/2022

Module Name
    flt_intro.py

Contents (12)
    esc(code)
    check_python(op_major, op_minor, version_name)
    define_model_limits(param_bounds)
    define_fault_core_limits(param_bounds)
    define_aperture_limits(param_bounds)
    define_let_limits(param_bounds)
    define_input_limits()
    check_blanks(fault_params)
    check_input_parameters(fault_params, param_bounds)
    cross_check_input(fault_params)
    flash(alone, op_major, op_minor, version_name)
    start(alone, yaml_file)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import sys                      # For program exit
import logging                  # For reporting errors
import time                     # For calculating runtime
from datetime import datetime   # For timing

import flt_decipher as deci     # For converting yaml data
import flt_file as fileop       # For file operations
import flt_message as mess      # For messages
import flt_warranty as war      # Print warranty

# Limit constants on flow
ERR_IN = "\n    > "             # Secondary Error messages indent
DUAL_CASE_DEPTH = 200.0         # m - Smallest near-surface distance.


def esc(code):
    """Define escape code.

    Parameters
    ----------
    code = (Int) print code to be precede by esc key.

    Returns
    -------
    code = escape code.
    """
    return f'\033[{code}m'


def check_python(op_major, op_minor, version_name):
    """Check that user is using proper Python version.

    Parameters
    ----------
    op_major = (int) Python version for code - major number (e.g., 3)
    op_minor = (int) Python version for code - minor number (e.g., 9)
    version_name = (str) Code version in string (e.g, "3.7")

    Returns
    -------
    N/A
    """
    # Get python version of system.
    logging.info('Checking Python Version')
    py_major, py_minor = sys.version_info[0:2]

    # Check if OK - they must use Version 3.9.
    # -> New versions are assumed backwards compatible
    if (py_major != op_major) or (py_minor < op_minor):
        msg = ('Python ' + version_name + ' is Required for this Code!'
               + ERR_IN
               + f'Python Version {py_major}.{py_minor} was Detected.')
        fileop.opx_problem(msg)

    # return None


def define_model_limits(param_bounds):
    """Define numeric limits/ranges of simulation run parameters.

    Parameters
    ----------
    pars_bounds = (dict) dictionary with min/max indices

    Returns
    -------
    param_bounds = (dict) updated dictionary
    """
    # model parameters.
    param_bounds['start_time'] = [0.0, 5000.0]
    param_bounds['end_time'] = [0.0, 10000.0]
    param_bounds['inject_end'] = [0.0, 10000.0]
    param_bounds['time_points'] = [0, 50]
    param_bounds['realizations'] = [0, 100000]

    return param_bounds


def define_fault_core_limits(param_bounds):
    """Define numeric limits/ranges of fault core.

    Parameters
    ----------
    pars_bounds = (dict) dictionary with min/max indices

    Returns
    -------
    param_bounds = (dict) updated dictionary
    """
    # FaultCore parameters:
    param_bounds['fault_probability'] = [0.0, 100.0]
    param_bounds['strike_mean'] = [0.0, 360.0]
    param_bounds['strike_sig'] = [0.0, 180.0]
    param_bounds['dip_mean'] = [10.0, 90.0]
    param_bounds['dip_std'] = [0.0, 90.0]
    param_bounds['length'] = [0.0, 1.0e+04]
    param_bounds['x_start'] = [-5.0e+07, 5.0e+07]
    param_bounds['y_start'] = [-5.0e+07, 5.0e+07]
    param_bounds['n_plates'] = [0, 100]
    param_bounds['sgr'] = [0, 100]
    param_bounds['state_correction'] = [0, 1.0]

    return param_bounds


def define_aperture_limits(param_bounds):
    """Define numeric limits/ranges of aperture parameters.

    Parameters
    ----------
    pars_bounds = (dict) dictionary with min/max indices

    Returns
    -------
    param_bounds = (dict) updated dictionary
    """
    # Aperture variability parameters.
    param_bounds['aperture_mean'] = [0.0, 101.0]
    param_bounds['aperture_std'] = [0.0, 20.0]
    param_bounds['aperture_min'] = [0.0, 1.0]
    param_bounds['aperture_max'] = [1.0E-03, 150.0]
    param_bounds['aperture_confines'] = [0.0, 50.0]

    return param_bounds


def define_let_limits(param_bounds):
    """Define numeric limits/ranges of LET model parameters.

    Parameters
    ----------
    pars_bounds = (dict) dictionary with min/max indices

    Returns
    -------
    param_bounds = (dict) updated dictionary
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
    """Define the numeric limits/ranges of fault input parameters.

    Parameters
    ----------
    N/A

    Returns
    -------
    param_bounds = (dict) dictionary with min/max indices
    """
    # Define dictionary of input bounds.
    param_bounds = {}

    # ModelParams:
    param_bounds = define_model_limits(param_bounds)

    # FaultCore:
    param_bounds = define_fault_core_limits(param_bounds)

    # Aperture variability parameters and limits:
    param_bounds = define_aperture_limits(param_bounds)

    # Field:
    param_bounds['aquifer_depth'] = [DUAL_CASE_DEPTH, 2.0e+04]
    param_bounds['aquifer_pressure'] = [1.0e+06, 6.0e+08]
    param_bounds['aquifer_temperature'] = [15.0, 180.0]
    param_bounds['inject_depth'] = [8.6e+02, 2.0e+04]   # see Model 1 break!
    param_bounds['field_pressure'] = [1.0e+05, 6.0e+07]
    param_bounds['final_pressure'] = [1.0e+05, 6.0e+07]
    param_bounds['inject_pressure'] = [7.0e+06, 6.0e+08]
    param_bounds['inject_temperature'] = [31.0, 180.0]
    param_bounds['inject_x'] = [-5.0e+07, 5.0e+07]
    param_bounds['inject_y'] = [-5.0e+07, 5.0e+07]

    # inject Conditions:
    param_bounds['salinity'] = [0.0, 8.0e+04]
    param_bounds['co2_density'] = [9.3E+01, 1.05e+03]
    param_bounds['co2_viscosity'] = [1.8E-05, 1.4e-04]
    param_bounds['brine_density'] = [8.8e+02, 1.08e+03]
    param_bounds['brine_viscosity'] = [1.5e-04, 1.6e-03]
    param_bounds['co2_solubility'] = [0.0, 2.0]

    # Aquifer Conditions:
    param_bounds['aqui_co2_density'] = [9.3E+01, 1.05e+03]
    param_bounds['aqui_co2_viscosity'] = [1.1E-05, 1.4e-04]
    param_bounds['aqui_brine_density'] = [8.8e+02, 1.08e+03]
    param_bounds['aqui_brine_viscosity'] = [1.5e-04, 1.6e-03]

    # Relative Permeability:
    param_bounds['resid_brine'] = [0.01, 0.35]
    param_bounds['resid_co2'] = [0.00, 0.35]
    param_bounds['perm_ratio'] = [0.00, 1.5]
    param_bounds['entry_pressure'] = [1.0E+02, 2.0E+06]
    param_bounds['zeta'] = [0.00, 5.0]

    # Stress:
    param_bounds['max_horizontal'] = [0.0, 5.0e+07]
    param_bounds['min_horizontal'] = [0.0, 5.0E+07]
    param_bounds['max_trend'] = [0.0, 180.0]

    # Define LET Parameters:
    param_bounds = define_let_limits(param_bounds)

    # Repository Input:
    param_bounds['reservoir_pressure'] = [1.0E+00, 1.0E+08]
    param_bounds['reservoir_saturation'] = [0.0, 1.0]
    param_bounds['coordinates'] = [-1.0E+05, 1.0E+05]

    return param_bounds


def check_blanks(fault_params):
    """Check whether if any input parameter is blank (= None).

    Parameters
    ----------
    fault_params = (dict) input parameters of fault component model

    Returns
    -------
    N/A

    Notes
    -----
    1. "logging" module still uses string format - added statements
         messages prior to logging.
    2. Error messages -> explicitly string outside of <logging>
        as <logging> still uses a string !
    """
    # Define error flag.
    err_flag = 0

    # Check all values for None.
    for key, val in fault_params.items():
        if val is None:
            msg = f'> Parameter <{key}> is Not Defined!'
            logging.error(msg)
            err_flag = -1

    # Report error if found. Missing data is Fatal!
    if err_flag == -1:
        fileop.opx_problem("Error: Input Parameter(s) are Missing!")

    # return None


def check_input_parameters(fault_params, param_bounds):
    """Check whether input parameters fall within specified limits.

    Parameters
    ----------
    fault_params = (dict) input parameters of fault component model
    param_bounds = (dict) bounds on fault parameters (w. 2 values)

    Returns
    -------
    stats = (int) error code, default = 0, if warning = -1

    Notes
    -----
    1. "logging" module still uses string format - add statements
        to string messages prior to logging.
    2. Some parameters are checked using a defined list.
    3. "title" is not checked.
    """
    # Define parameters for checking.
    stats = 0
    true_false_parameters = ['approach', 'switch', 'plot', 'time_input',
                             'use_si', 'vary_aperture', 'skip_output']
    truth_list = [True, False]
    relative_list = ["LET", "BC"]
    profile_parameters = ['profile']
    profile_list = [0, 1, 2]

    # Check if values for various conditions.
    for key, val in fault_params.items():
        # Check if relative permeability models are present.
        if 'relative' in key:
            if val not in relative_list:
                txt_not_model = (f'> Parameter <{key}> not Recognized as a '
                                 + 'Relative Permeability Model Type.')
                logging.error(txt_not_model)
                stats = -1

        # Check unit types.
        elif 'units' in key:
            if "MET" not in val.upper() and "ENG" not in val.upper():
                txt_units = '> Only Metric or English Units can be Specified!'
                logging.error(txt_units)
                stats = -1

        # Check profile types = 0,1,2.
        elif [ele for ele in profile_parameters if ele in key]:
            if val not in profile_list:
                txt_not_model = (f'Parameter <{key}> Value is not Correct '
                                 + f'and has a Value of <{val}>!')
                logging.error(txt_not_model)
                stats = -1

        # Check all True/False parameters for Yes/No.
        elif [ele for ele in true_false_parameters if ele in key]:
            if val not in truth_list:
                txt_wrong = (f'Parameter <{key}> is not Defined Correctly '
                             + f'and has a value of <{val}>!')
                logging.error(txt_wrong)
                stats = -1

        # Check other numeric values to be within bounds.
        elif key in param_bounds:
            if (val < param_bounds[key][0]) \
                    or (val > param_bounds[key][1]):
                txt_outside = (f'> Parameter <{key}> is Outside of the '
                               + f'Recommended Bounds with a Value of {val}!')
                logging.error(txt_outside)
                stats = -1

        else:
            # Do not check title.
            # Otherwise - everything else or missing in fault_controls.
            if'title' not in key:
                txt_not_fault = (f'> Parameter <{key}> Not Recognized as a '
                                 + 'Fault parameter.')
                logging.error(txt_not_fault)
                stats = -1

    return stats


def cross_check_input(fault_controls):
    """Check if input is numerically logical for variability.

    Parameters
    ----------
    fault_controls = (dict) input parameters of fault component

    Returns
    -------
    statx = (int) error code, default = 0, if warning = -1

    Notes
    -----
    "logging" module still uses string format -
        set string messages prior to calling <logging>.
    """
    # Set default return.
    statx = 0

    # Run check only if no file input for permeability.
    if not fault_controls['aperture_approach']:

        # Run checks on min and max values and mean.
        if fault_controls['aperture_min'] > fault_controls['aperture_max']:
            causitive = "aperture"
            txt_cross = f'> Minimum {causitive} Exceeds Maximum {causitive}!'
            logging.error(txt_cross)
            statx = -1

        data_mean = fault_controls['aperture_mean']
        if (data_mean > fault_controls['aperture_max'] or
                data_mean < fault_controls['aperture_min']):
            causitive = "aperture"
            txt_outside = f'> Mean {causitive} is Outside Min./Max. Bounds!'
            logging.error(txt_outside)
            statx = -1

    # Run check on near surface specifications.
    flow_approach = fault_controls['near_surface_approach']
    profile_type = fault_controls['profile_type']
    condition_1 = not flow_approach and profile_type > 0
    condition_2 = flow_approach and profile_type < 1

    if condition_1 or condition_2:
        msg = (f'> Profile = {profile_type} and Near-Surface Approach = '
               f'<{flow_approach}> are Inconsistent')
        logging.error(msg)
        statx = -1

    return statx


def flash(alone, op_major, op_minor, version_name):
    """Show warranty and first line.

    Parameters
    ----------
    alone = (bool) status of code operation; True = if stand-alone
    op_major = (int) Python version for code - major number (e.g., 3)
    op_minor = (int) Python version for code - minor number (e.g., 7)
    version_name = (str) Code version string (e.g, "1.4")

    Returns
    -------
    None
    """
    if alone:
        print('\n\n' + esc('31;1') + '  +*** ' + 'FAULT_FLO'
              + ' - Version: ' + version_name + " ***" + esc(0)
              + "                 Python "
              + str(op_major) + "." + str(op_minor))
        war.show_warranty()
        mess.echo_status(alone, "BEGIN")


def start(alone, yaml_file, op_major, op_minor, version_name):
    """Open the YAML control file, and establishes controls and checks.

    Parameters
    ----------
    alone = (bool) status of code operation; True = if stand-alone
    yaml_file = (str) name of input file
    op_major = (int) Python version for code - major number (e.g., 3)
    op_minor = (int) Python version for code - minor number (e.g., 7)
    version_name = (str) Code version in string (e.g, "1.4")

    Returns
    -------
    fault_controls = (dict) input parameters (dictionary)
    param_bounds = (dict) parameter bounds (dictionary)
    """
    # Set duration timer and set date.
    timer = time.monotonic()
    flag = datetime.now()

    # Show header statements on console - construct Python version.
    flash(alone, op_major, op_minor, version_name)

    # Check Python version.
    check_python(op_major, op_minor, version_name)

    # Initiate input control parameters from text file.
    fault_controls = {}
    fault_controls = deci.acquire_control_parameters(yaml_file, fault_controls)

    if alone:
        mess.echo_status(alone, "CLEAN")

    # Wipe output directory to remove old output files.
    fileop.clean_output()

    # Check for any missing data.
    check_blanks(fault_controls)

    # Identify input bounds, check values and report errors.
    param_bounds = define_input_limits()
    status = check_input_parameters(fault_controls, param_bounds)
    status2 = cross_check_input(fault_controls)

    # Consider all warnings fatal! - and exit program.
    if status < 0 or status2 < 0:
        fileop.opx_problem("")  # Do not re-write error message

    # Save time values in fault_controls for later use.
    fault_controls['code_start'] = timer
    fault_controls['record_date'] = flag
    fault_controls['seg_length'] = (fault_controls['length']
                                    / fault_controls['n_plates'])

    return fault_controls, param_bounds


#
# -----------------------------------------------------------------------------
# End of module
