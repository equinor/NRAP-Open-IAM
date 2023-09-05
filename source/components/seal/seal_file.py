#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions for basic file operations for seal ROM.

Author: Ernest N. Lindner
Date: 08/19/2022

Module Name
    seal_file

Contents (16)
    check_python()
    check_file_floats(data_array, data_type, data_limits)
    check_file_values(data_array, data_type, data_limits)
    opx_problem(msg, err)
    io_snag(err)
    io_yml(err)
    issue_warning(reason)
    data_error(reason, subdirectory, file_selected)
    acquire_yaml_file(yaml_file_selected)
    acquire_data_array(new_sub, file_selected, extent)
        ----
    acquire_csv_list(destination, subdirectory_path, file_name)
    clean_output()
    get_path_name(subdirectory, file_selected, extension)
    find_local_path(file_selected)
    echo_status(alone, msg)
    terminate_code(kill, alive)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import os                           # For paths
import sys                          # For program exit
import time                         # For pause for warning to process
import logging                      # For reporting errors
import csv                          # For file input
from pathlib import Path            # To allow path definition - Windows/Unix
import yaml                         # For structured control file input
import numpy as np                  # for array

import seal_config as scfg           # IO directory and file names

# Print Constants
LED_IN = "\n  --> "                 # Message start
ERR_IN = "\n      -> "              # Secondary error messages start
LIX_IN = "        -> "              # Common error message
FIN_IN = "  >>> "                   # Last line start
SPX_IN = "      "                   # Spaced line start
ECHO = False                        # Control to print for debugging
END_TEXT = True                     # Add message start at last line


def check_python(code_major, code_minor):
    """Check that user is using correct version of Python.

    Parameters
    ----------
    code_major = (int) Python version for code - major number (e.g., 3)
    code_minor = (int) Python version for code - minor number (e.g., 7)

    Returns
    -------
    1. New versions are assumed backwards compatible with development version.
       py_minor >= code_minor => True
    """
    # Get python version of current system.
    logging.info('Checking Python Version')
    py_major, py_minor = sys.version_info[0:2]

    # Check if error - must use Python version consistent with code.
    if (py_major != code_major) or (py_minor < code_minor):
        msg = (f'Python {code_major}.{code_minor} is Required for this Code!'
               + '\n      -> '
               + f'Version Python {py_major}.{py_minor} was Detected')
        logging.error(msg)
        terminate_code(1, True)

    # return None


def check_file_floats(data_array, data_type, data_bounds):
    """Check float data from seal file if it is within assumed bounds.

    Parameters
    ----------
    data_array = (NumPy array) data to check
    data_type = (str) data description
    data_bounds = (dict) minimum/maximum bounds

    Returns
    -------
    statx = (int) error flag
        =0 no error
        =1 error in data
    """
    # Define error flag.
    statx = 0

    # Check values to be inside min/max.
    for key in data_array:
        if (key < data_bounds[0]) or (key > data_bounds[1]):
            msg = (data_type + ' From File is Outside of Defined ' +
                   f'Bounds with Value of {key}!')
            logging.error(msg)
            statx = -1

    return statx


def check_file_values(data_array, data_type, data_bounds):
    """Check integer data from seal file if it has defined values.

    Parameters
    ----------
    data_array = (array) data to check
    data_type = (str) data description
    data_bounds = (dict) minimum/maximum bounds

    Returns
    -------
    statx = (int) error flag
    =0 no error
    =1 error in data
    """
    # Define error flag.
    statx = 0

    # Check values to either equal to min or max.
    for key in data_array:
        if (int(key) != data_bounds[0]) and (int(key) != data_bounds[1]):
            msg = (data_type + ' From File is Incorrect ' +
                   f'with a Value of {key}!')
            logging.error(msg)
            statx = -1

    return statx


def opx_problem(reason, err=''):
    """For operation error, provides error message and stops.

    Parameters
    ----------
    reason = (str) function error description
    err = (str) OS error message

    Returns
    -------
    None
    """
    # Print error message.
    logging.error(reason)

    # Print details only if desired.
    if err != '' and ECHO:
        msg = 'OS Report: \n' + ERR_IN + f' {err}'
        logging.error(msg)

    terminate_code(1, True)

    # return None


def io_snag(err, subdirectory, file_selected):
    """For file IO errors, provide short error message and then stop.

    Parameters
    ----------
    err = (str) OS error message
    subdirectory = (str) directory where code was looking
    file_selected = (str) input file

    Returns
    -------
    None
    """
    # Convert directory to string.
    directory = str(subdirectory)

    # Show multi-line detailed message with location.
    msg = ERR_IN + "OS Error - File Not Found/Accessed!"
    msg += ERR_IN + "Looking for File: '" + file_selected + "'"
    msg += ERR_IN + "Looking in Directory: \n" + LIX_IN + directory
    msg += ERR_IN + "Please Correct Issue and Retry."

    # Send error message.
    logging.error(msg)

    # Print specific details only if desired.
    if err != '' and ECHO:
        msg = 'OS Report: \n' + ERR_IN + f' {err}'
        logging.error(msg)

    terminate_code(1, True)

    # return None


def io_yml(err):
    """For a file error in YAML file, provide error message and then stop.

    Parameters
    ----------
    err = (str) OS error message

    Returns
    -------
    None
    """
    # Show multi-line detailed message if position available.
    if hasattr(err, 'problem_mark'):
        msg = "Problem While Parsing YAML Control File."
        mark = err.problem_mark
        msg += (ERR_IN + "Problem Location in File: "
                f'line = {(mark.line+1)}, column = {(mark.column+1)}')
        msg += (ERR_IN + 'Cause: ' + str(err.problem) + ' '
                + str(err.context))
        msg += ERR_IN + "Please Correct Control File and Retry."
        logging.error(msg)

    # Otherwise, show general message.
    else:
        msg = "Problem While Parsing YAML Control File."
        msg += ERR_IN + "Undefined Parse Error in File."
        msg += ERR_IN + f'Error Report: {err}'
        msg += ERR_IN + "Please Correct Control File and Retry."
        logging.error(msg)

    # Exit.
    terminate_code(1, True)

    # return None


def issue_warning(reason):
    """Provide a warning for missing or wrong data.

    Parameters
    ----------
    reason = (str) description of issue
    subdirectory = (str) directory where code was looking
    file_selected = (str) input file

    Returns
    -------
    None

    Notes
    -----
    1. This function does not use standard logging warn method.
    2. This function does not stop operation.
    """
    # Show message formatted for printout.
    print(SPX_IN + reason, flush=True, file=sys.stderr)
    # return None


def data_error(reason, subdirectory, file_selected):
    """Provide an error message for missing or wrong data.

    Parameters
    ----------
    reason = (str) description of issue
    subdirectory = (str) directory where code was looking
    file_selected = (str) input file

    Returns
    -------
    None
    """
    # Show multi-line detailed message with location.
    msg = reason
    subdirectory = str(subdirectory)
    msg += ERR_IN + "Looking in File: '" + file_selected + "'"
    msg += ERR_IN + "In Directory: " + ERR_IN + subdirectory
    msg += ERR_IN + "Please Correct Issue and Retry."
    logging.error(msg)

    terminate_code(1, True)

    # return None


def acquire_yaml_file(yaml_file_selected):
    """Open text file with YAML options and read data.

    Parameters
    ----------
    yaml_filename = (str) YAML file for seal operations

    Returns
    -------
    docs = (stream) read in of YAML file
    """
    # Initialize tracker.
    docs = None

    # Get full path for yaml control file.
    subdirectory = ""
    subdirectory_path, destination = \
        get_path_name(subdirectory, yaml_file_selected, "yaml")

    # Check if file exists before download.
    if os.path.exists(destination):
        try:
            with open(destination, 'r', encoding="utf8") as y_stream:
                docs = yaml.safe_load(y_stream)
        except yaml.YAMLError as err:
            io_yml(err)
    else:
        # No directory found!
        io_snag(" ", subdirectory_path, yaml_file_selected)

    return docs


def acquire_data_array(new_sub, file_selected, extent):
    """Read data from a *.csv file into a NumPy array.

    Parameters
    ----------
    new_sub = (str) subdirectory name where file resides
    file_select = (str) name of selected file
    extent = (int) number of elements expected in array

    Returns
    -------
    out_array = (array) (flat / 1D) NumPy array from file

    Notes
    -----
    Reading file assumes data has 2 header lines and is in 2D format!
    out_array is a 1D (flattened) array for use in code.
    Extension assumed to be *.txt
    """
    # Provide default.
    out_array = None

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
    except ValueError:
        # Value error.
        data_error("Value Error -> Check Number of Columns in File!",
                   sub_path, file_selected)
    except EOFError:
        # End of File (EOF) interrupt or no data.
        data_error("EOF Error -> Check Data in File!",
                   sub_path, file_selected)

    # Double check array if present.
    if out_array.size == 0:
        # No data.
        data_error("Data Error - Array Size = 0! "
                   + "-> No Data Found in File.", sub_path, file_selected)

    # Double check if input array has the correct size.
    if extent != out_array.size:
        # Too little data.
        if extent > out_array.size:
            data_error("Data Error - "
                       + "Array Size in File is Too Small!",
                       sub_path, file_selected)
        else:
            data_error("Data Error - "
                       + "Array Size in File is Too Large!",
                       sub_path, file_selected)

    # If no error, flatten array from 2D to 1D - using "C" method.
    out_array = out_array.flatten('C')

    return out_array


def acquire_csv_list(destination, subdirectory_path, file_name):
    """Input data from CSV file as list.

    Parameters
    ----------
    file_name = (str) name of input file (w.o. destination directory)

    Returns
    -------
    user_list = (list) list of fracture data

    Notes
    -----
    1. Used only for fracture data input.
    """
    try:
        with open(destination, 'r', encoding="utf8") as csvfile:
            csv_reader = csv.reader(csvfile)

            # Skip 2 header lines - use first header for checking.
            header = next(csv_reader, None)
            next(csv_reader)

            # Do not read file without a header / empty file.
            if header is None:
                # No data.
                data_error("Data Error - Frac Array Header is 'None'! "
                           + "-> No data found in file.",
                           subdirectory_path, file_name)
            else:
                # Create list - iterate over each row after the header.
                user_list = []
                for row in csv_reader:
                    user_list.append(row)

    except OSError as err:
        io_snag(err, subdirectory_path, file_name)
    except ValueError:
        data_error("Value Error: Check Data in Frac Input File!",
                   subdirectory_path, file_name)
    except EOFError:
        # Data interrupt or no data.
        data_error("EOFError: -> Check Data in Frac Input File!",
                   subdirectory_path, file_name)

    # Double-check array if present.
    if not user_list:
        # No data.
        data_error("File Problem: No Data Found in File in <aquire_csv_list>!",
                   subdirectory_path, file_name)

    return user_list


def clean_output():
    """Delete files in destination/output directory.

    Parameters
    ----------
    N/A

    Returns
    -------
    None

    Notes
    -----
    1. For Path details see:
    https://medium.com/@ageitgey/python-3-quick-tip-the-easy-way-to-deal
    -with-file-paths-on-windows-mac-and-linux-11a072b58d5f
    """
    # Construct full path from script location to directory file.
    directory_path = os.path.dirname(os.path.abspath(__file__))

    # Using Path, construct directory path with correct "/" or "\".
    panther = Path(directory_path)
    subdirectory_path = panther / scfg.OUTPUT_DIRECTORY

    # Check subdirectory - it may not exist.
    if not os.path.exists(subdirectory_path):
        # Create new directory, if none exists.
        os.makedirs(subdirectory_path)
    else:
        # Create a list of all files in current output directory.
        old_file_list = os.listdir(subdirectory_path)

        # Define file types to be deleted in output dir. Note "." in list.
        file_remove_list = [".csv", ".txt", ".npy", ".png"]

        for file_name in old_file_list:
            # Create path object for each file.
            old_file_path = subdirectory_path / file_name

            # Select files with specific suffix and delete.
            if old_file_path.suffix in file_remove_list:
                try:
                    old_file_path.unlink()
                except OSError as detail:
                    msg = ('File Error While Deleting File' +
                           '- Is File Open?' + LED_IN + f'FILE: {file_name}')
                    opx_problem(msg, str(detail))

    # return None


def get_path_name(subdirectory, file_selected, extension):
    """Provide file path.

    Parameters
    ----------
    subdirectory = (str) directory name where code is looking
    file_selected = (str) file name of object to be saved
    extension = (str) file extension for file

    Returns
    -------
    subdirectory_path = (Path) path to subdirectory (only)
    destination = (str) full path to file

    Notes
    -----
    Assumes destination directory is a subdirectory of executable
    """
    # Directory is current location.
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
    """Provide file path for a file in the source directory.

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

    return str(destination)


def echo_status(alone, msg):
    """Echo status message to console, if in stand-alone mode.

    Parameters
    ----------
    alone = (bool) control parameter to show results
    msg = (str) text message

    Returns
    -------
    None
    """
    if alone:
        print(LED_IN + msg, flush=True)
        # logging.info(msg)

    # return None


def terminate_code(kill, alone):
    """Echo message to console on termination of program due to error.

    Parameters
    ----------
    kill = (int) flag to terminate code
         = 0 normal
         = 1 error
    alone = (bool) control parameter to show results to console.

    Returns
    -------
    None
    """
    # Clear all warnings first and pause briefly to allow logging to finish.
    logging.shutdown()
    sys.stdout.flush()
    time.sleep(0.1)

    # Provide final message - if code immediately ends, provide alternate end.
    if kill == 0:
        msg = LED_IN + "SEAL FLUX FINISHED."
    else:
        msg = LED_IN + "TERMINATING PROGRAM DUE TO ERROR."

    if END_TEXT and alone:
        # Add text and pause if needed.
        print(msg, end='', file=sys.stdout)
        print(" (Press <ENTER> to Continue)", end='', file=sys.stdout)
        input()
    else:
        # No pause at code end.
        print(msg, file=sys.stdout)

    # Error - use system to exit code if error.
    if kill == 1:
        raise SystemExit(1)

    # return None


#
# -----------------------------------------------------------------------------
# - End of module
