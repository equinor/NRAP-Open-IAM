#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions for file operations including opening YAML file.

Author: Ernest N. Lindner
Date: 08/23/2022

Module Name
    flt_file

Contents (12)
    check_file_data(data_array, data_type, data_limits)
    opx_problem(msg, err)
    io_snag(err, subdirectory, file_selected)
    io_yml(err, subdirectory, file_selected)
    data_error(subdirectory, file_selected)
    acquire_yaml_file(yaml_file_selected)
    acquire_data_array(new_sub, file_selected, pluria)
    clean_output()
    get_path_name(subdirectory, file_selected, extension)
    find_local_path(file_selected)
    ---
    end_message(alone)
    terminate_code()

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import os                       # For paths
import sys                      # For program exit
import time                     # For pause
import logging                  # For reporting errors
from pathlib import Path        # To allow path definition - Windows/Unix
import yaml                     # For structured control file input
import numpy as np              # for array

import flt_config as scfg       # IO directory and file names

# Statement constants
COD_IN = "  --> "               # Message start - no new line
LED_IN = "\n  --> "             # Message start
RUN_IN = "\n  > "               # End message
ERR_IN = "\n    > "             # Secondary Error messages start
ECHO = False                    # Control to print for debugging
CLOSER = True                   # Control to add final line if used by IED.


def check_file_data(data_array, data_type, data_bounds):
    """Check data from a file if it is within assumed bounds.

    Parameters
    ----------
    data_array = (1D NumPy array) data to check
    data_type = (str) data description
    data_bounds = (1D array) minimum/maximum bounds

    Returns
    -------
    stats = error flag
          = 0 no error
          =-1 error in data
    -----------
    """
    # Define error flag.
    stats = 0

    # Check values to be inside min/max - 1D Data array.
    for indx in range(data_array.size):
        val = data_array[indx]
        if (val < data_bounds[0]) or (val > data_bounds[1]):
            msg = (data_type + ' from file is outside of defined '
                   f'bounds with value of {val}!')
            logging.error(msg)
            stats = (-1)

    return stats


def opx_problem(reason, err=''):
    """For operation errors, provide error message and stops.

    Parameters
    ----------
    reason = (str) function error string
    err = (str) OS error message

    Returns
    -------
    N/A
    """
    # Print error message of problem.
    print(flush=True)
    if reason != "":
        logging.error(reason)

    # Print details only if desired.
    if err != '' and ECHO:
        msg = 'OS Report: \n' + ERR_IN + f'{err}'
        logging.error(msg)

    # Stop code.
    terminate_code()

    # return None


def io_snag(err, subdirectory, file_selected):
    """For file IO errors, provide short io error message and stops.

    Parameters
    ----------
    reason (str) user-friendly error message
    err = (str) OS error message
    subdirectory = (str) directory where code was looking
    file_selected = (str) input file

    Returns
    -------
    N/A  (end program)
    """
    # Convert path to string and set default message.
    print(flush=True)
    msg = "OS Error - File Not Found/Accessed!"
    sub_directory = str(subdirectory)
    file_name = str(file_selected)

    # Show detailed message with location if possible.
    if sub_directory != "":
        msg += ERR_IN + "Looking for file: '" + file_name + "'"
        msg += (ERR_IN + "Looking in directory:" + ERR_IN + "    "
                + sub_directory)
        msg += ERR_IN + "Please correct issue and retry."
    logging.error(msg)

    # Print specific details, only if desired.
    if err != '' and ECHO:
        msg = 'OS Report: \n' + ERR_IN + f'{err}'
        logging.error(msg)

    # Stop code.
    terminate_code()

    # return None


def io_yml(err):
    """For file error in YAML file, provide short error message and stop.

    Parameters
    ----------
    err = (str) OS error message

    Returns
    -------
    N/A
    """
    # Show multi-line detailed message if position available.
    print(flush=True)
    if hasattr(err, 'problem_mark'):
        msg = "Problem While Parsing YAML Control File."
        mark = err.problem_mark
        msg += ERR_IN + "Problem location in file: "
        msg += f'line = {mark.line + 1}, column = {mark.column + 1}.'
        msg += (ERR_IN + "Cause: " + str(err.problem) + ' '
                + str(err.context))
        msg += ERR_IN + "Please correct control file and retry."
        logging.error(msg)

    # Otherwise, show general message.
    else:
        msg = "Problem While Parsing YAML Control File."
        msg += ERR_IN + "Undefined parse error in file."
        msg += ERR_IN + f'Error Report: {0}'
        msg += ERR_IN + "Please correct control file and retry."
        logging.error(msg)

    # Stop code.
    terminate_code()

    # return None


def data_error(reason, subdirectory, file_selected):
    """For missing or wrong data in file, provide an error message.

    Parameters
    ----------
    reason = (str) description of issue
    subdirectory = (str) directory where code was looking
    file_selected = (str) input file

    Returns
    -------
    N/A
    """
    # Show multi-line detailed message with location.
    print(flush=True)
    msg = reason
    msg += ERR_IN + "Looking in file: '" + file_selected + "'"
    msg += ERR_IN + "In directory: " + ERR_IN + subdirectory
    msg += ERR_IN + "Please correct issue and retry."
    logging.error(msg)

    # Stop code.
    terminate_code()

    # return None


def acquire_yaml_file(yaml_file_selected):
    """Open text file with YAML options and read data.

    Parameters
    ----------
    yaml_filename = (str) YAML file for fault operations

    Returns
    -------
    docs = YAML file
    """
    # Get full path for startup control file.
    subdirectory = ""
    subdirectory_path, destination = \
        get_path_name(subdirectory, yaml_file_selected, "yaml")

    # Check Open file.
    try:
        with open(destination, "r", encoding="utf8") as streamer:
            docs = yaml.safe_load(streamer)
    except OSError as err:
        io_snag(err, subdirectory_path, yaml_file_selected)
    except yaml.YAMLError as err:
        io_yml(err)
        docs = 0.0  # for inspection

    return docs


def acquire_data_array(new_sub, file_selected, pluria):
    """Read data from a *.csv file into a NumPy array.

    Parameters
    ----------
    new_sub = (str) subdirectory name where file resides
    file_select = (str) name of selected file
    pluria = (int) number of elements expected in array

    Returns
    -------
    out_array = (array) NumPy array (flat / 1D)

    Notes
    -----
    1. Data file assumes 2 header lines and is 2D!
    2. Extension assumed to be *.txt.
    """
    # Open the current directory and source file for input operations.
    sub_path, source_name = get_path_name(new_sub, file_selected,
                                          scfg.EXTENSION_TXT)
    try:
        # Read data as a NumPy array, excluding first two (2) header lines.
        out_array = np.genfromtxt(source_name, delimiter=",",
                                  autostrip=True, skip_header=2)
    except OSError as err:
        # File not found or accessible.
        io_snag(err, sub_path, file_selected)
        out_array = np.empty(0)  # for inspection
    except ValueError:
        # Value error.
        data_error("Value Error -> Check number of columns in file!",
                   sub_path, file_selected)
        out_array = np.empty(0)  # for inspection
    except EOFError:
        # End of File (EOF) interrupt or no data.
        data_error("EOF Error -> Check Data in File!",
                   sub_path, file_selected)
        out_array = np.empty(0)  # for inspection

    # Double check array if present.
    if out_array.size == 0:
        # No data.
        data_error("Data Error - Array size = 0! "
                   + "-> No data found in file.", sub_path, file_selected)

    # Double check if input array has the correct size.
    if pluria != out_array.size:
        # Too little data.
        if pluria > out_array.size:
            data_error("Data Error - "
                       + "Array size in file is too small!",
                       sub_path, file_selected)
        else:
            data_error("Data Error - "
                       + "Array size in file is too large!",
                       sub_path, file_selected)

    # If no error, flatten array from 2D to 1D - using "C" method.
    out_array = out_array.flatten('C')

    return out_array


def clean_output():
    """Delete output files in destination directory from prior run.

    Parameters
    ----------
    N/A

    Returns
    -------
    N/A

    Notes
    -----
    1. For Path details see:
    https://medium.com/@ageitgey/python-3-quick-tip-the-easy-way-to-deal
    -with-file-paths-on-windows-mac-and-linux-11a072b58d5f
    """
    # Construct full path from script location to directory file.
    directory_path = os.path.dirname(os.path.abspath(__file__))

    # Using Path, construct directory profile with correct "/" or "\".
    panther = Path(directory_path)
    subdirectory_path = panther / scfg.OUTPUT_DIR

    # Check directory - it may not exist.
    if not os.path.exists(subdirectory_path):
        # Create new directory, if none exists
        os.makedirs(subdirectory_path)
    else:
        # Create a list of all files in current output directory.
        old_file_list = os.listdir(subdirectory_path)

        # Delete all specific files in output directory. Note "." in list
        file_remove_list = [scfg.EXTENSION_CSV, scfg.EXTENSION_TXT,
                            scfg.EXTENSION_PNG]

        for file_name in old_file_list:
            # Create Path object for file.
            old_file_path = subdirectory_path / file_name

            # Select files with specific suffix.
            if old_file_path.suffix in file_remove_list:
                try:
                    old_file_path.unlink()
                except OSError as detail:
                    msg = 'File Error While Deleting File -> Is File Open?'
                    msg += COD_IN + f'FILE: {file_name}'
                    opx_problem(msg, str(detail))

    # return None


def get_path_name(subdirectory, file_selected, extension):
    """Provide file path.

    Parameters
    ----------
    subdirectory = (str) directory name where code is looking
    file_selected = (str) file name of object to be saved
    extension = (str) required file extension for file

    Returns
    -------
    subdirectory_path = (str) path to subdirectory (only)
    destination = (str) full path to file

    Notes
    -----
    Assumes destination directory is a subdirectory of executable.
    """
    # Define subdirectory path from script file, if defined.
    # directory_path = os.getcwd() # working directory.
    directory_path = os.path.dirname(os.path.abspath(__file__))

    # Define subdirectory path from script file, if it is to be defined.
    if subdirectory == "":
        subdirectory_path = Path(directory_path)
    else:
        pointer = Path(directory_path)
        subdirectory_path = pointer / subdirectory

    # Ensure that extension is in file name.
    if not file_selected.endswith(extension):
        file_selected = file_selected + extension

    # Define destination path with file name.
    destination = subdirectory_path / file_selected

    return subdirectory_path, str(destination)


def find_local_path(file_selected):
    """Provide file path for a file in source directory.

    Parameters
    ----------
    file_selected = (str) file name

    Returns
    -------
    destination = (str) full path to file
    """
    # Define directory path from script file, if defined.
    directory_path = os.path.dirname(os.path.abspath(__file__))
    panther = Path(directory_path)

    # Define destination path with file name.
    destination = panther / file_selected

    return destination


def end_message(alone):
    """Send a message to console on closing program.

    Parameters
    ----------
    alone = (bool) control parameter to show results to console

    Returns
    -------
    N/A
    """
    if alone:
        if CLOSER:
            print(LED_IN + "FAULT_FLO FINISHED.", end='', file=sys.stdout)
            print(" (Press <ENTER> to Close Program)", end='', file=sys.stdout)
            input()
        else:
            print(RUN_IN + "FAULT_FLO RUN COMPLETED.", file=sys.stdout,
                  flush=True)

    # return None


def terminate_code():
    """Send a message on termination of program due to error + stop.

    Parameters
    ----------
    N/A

    Returns
    -------
    N/A

    Notes
    -----
    1. Prints to screen at all times.
    """
    # Clear all warnings and pause briefly to allow logging to finish.
    logging.shutdown()
    sys.stdout.flush()
    time.sleep(0.1)

    if CLOSER:
        print(COD_IN + "TERMINATING PROGRAM DUE TO ERROR.", end='',
              file=sys.stdout)
        print(" (Press <ENTER> to Close Program)", end='', file=sys.stdout)
        input()
    else:
        print(LED_IN + "TERMINATING PROGRAM DUE TO ERROR.", file=sys.stdout)

    # Error exit.
    raise SystemExit(1)

    # return None


#
# -----------------------------------------------------------------------------
# End of module
