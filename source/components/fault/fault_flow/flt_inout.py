#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions to save output.

Author: Ernest N. Lindner
Date: 08/23/2022

Module Name
    flt_inout

Contents (10)
    save_csv_results(file_name, data_block, data_title)
    cache_step_results(sim, time_value, sim_title, mass_flow_co2
                       mass_flow_brine)
    cache_perm(fault_controls, sim_numbr, fault)
    save_fault_results(sim_numbr, entitle, fault)
    save_sim_list(file_name, sim_listing, title)
    cache_sim_results(sim_numbr, sim_listing, entitle)
    input_apertures(fault_controls, param_bounds)
    get_data_lines(destination, file_location, n_plates, skip_lines)
    input_reservoir_data(sequence, n_plates, param_bounds)
    list_data_head():

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
from itertools import repeat    # For line input
import csv                      # For saving a list of lists
import numpy as np              # For arrays & input

import flt_config as scfg       # Provides IO directory and file names
import flt_file as fileop       # For file operations

# Other constants
ECHO = False                    # Print header lines


def save_csv_results(file_name, data_block, run_title, data_title):
    """Save a data block to a *.csv file.

    Parameters
    ----------
    file_name = (str) Name of output file (string)
    data_block = (array) NumPy 2D data array
    data_title = (str) header for data file

    Returns
    -------
    N/A
    """
    # Construct full path for file.
    sub_path, destination = fileop.get_path_name(scfg.OUTPUT_DIR, file_name,
                                                 scfg.EXTENSION_CSV)

    # Write data - with two header lines.
    overall_title = " > Run:  " + run_title + "\n "
    overall_title += data_title
    try:
        np.savetxt(destination, data_block, fmt="%15.5e", delimiter=",",
                   comments=" ", header=overall_title)
    except OSError as err:
        fileop.io_snag(err, sub_path, file_name)

    # return None


def cache_step_results(simulation, time_value, sim_title, mass_flow_co2,
                       mass_flow_brine):
    """Control data storage of a single time-step - if debugging.

    Parameters
    ----------
    simulation = (int) realization number
    time_value = (float) time value of step (float)
    sim_title = (str) title of simulation
    mass_flow_co2 = (2D NumPy array) CO2 mass for realization
    mass_flow_brine = (2D NumPy array) brine mass for realization

    Returns
    -------
    N/A
    """
    # Construct file names.
    insert = str(simulation + 1) + "_at_Time_" + str(time_value)
    output1 = scfg.BRINE_INTRO + insert + scfg.EXTENSION_CSV
    output2 = scfg.CO2_INTRO + insert + scfg.EXTENSION_CSV

    # --> Store time-variant data for time X.
    title = " CO2 Flow Through Fault (tonnes)"
    save_csv_results(output2, mass_flow_co2, sim_title, title)

    title = " Brine Flow Through Fault (tonnes)"
    save_csv_results(output1, mass_flow_brine, sim_title, title)

    # return None


def cache_perm(fault_controls, sim_numbr, fault):
    """Store permeability values for a simulation - if debugging.

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    sim_numbr = (int) realization number
    fault = (class) group of fault segments

    Returns
    -------
    N/A
    """
    # set variables for clarity.
    entitle = fault_controls['title']
    segs = fault_controls['n_plates']

    # Construct unique file name for saving data.
    current_data = scfg.PERM_NAME + str(sim_numbr + 1)
    output_nomer = current_data + scfg.EXTENSION_CSV

    # Create perm array from results of each fault plate.
    perm_array = np.zeros(segs)
    for numbr, part in enumerate(fault):
        perm_array[numbr] = part.perm

    # Create header lines and save results to file - with average.
    perm_ave = np.average(perm_array)
    data_name = "Permeability Data (m^2) for Each Segment"
    data_name += f"   -  Average = {perm_ave:.2f}"
    save_csv_results(output_nomer, perm_array, entitle, data_name)

    # return None


def save_fault_results(sim_numbr, entitle, fault):
    """Save fault plate values for a simulation to a *.csv file.

    Parameters
    ----------
    sim_numbr = (int) realization number
    entitle = (str) title of run
    fault = (class) group of fault segments

    Returns
    -------
    N/A
    """
    # Construct full path for file and run title.
    file_name = scfg.PLATE_NAME + str(sim_numbr + 1) + scfg.EXTENSION_CSV
    overall_title = "  Run:  " + entitle

    # Construct full path to results directory.
    #  -- Check for extension "txt".
    sub_path, destination = fileop.get_path_name(scfg.OUTPUT_DIR, file_name,
                                                 scfg.EXTENSION_CSV)

    try:
        # Open file and write header [single element].
        with open(destination, "w", newline='', encoding="utf8") as csv_file:
            plate_writer = csv.writer(csv_file, delimiter=',')
            plate_writer.writerow([overall_title])

            # Write plate data.
            plate_writer.writerow(list_data_head())
            for numbr, part in enumerate(fault):
                plate_writer.writerow(part.list_data(numbr + 1))

    except OSError as err:
        fileop.io_snag(err, sub_path, file_name)

    # return None


def save_sim_list(file_name, sim_listing, entitle):
    """Save fault output from simulation in a *.csv file.

    Parameters
    ----------
    file_name = (str) Name of output file (string)
    sim_listing = (NumPy 2D array) data
    entitle = (str) header for data file

    Returns
    -------
    N/A
    """
    # Construct full path for file.
    sub_path, destination = fileop.get_path_name(scfg.OUTPUT_DIR, file_name,
                                                 scfg.EXTENSION_CSV)

    # Define header strings.
    overall_title = "  Run:  " + entitle
    field_names = [" Sim. No.", " Time (yrs)", " CO2 Flow (tonnes)",
                   " Brine Flow (tonnes)"]

    # Write data with 2 header lines.
    try:
        with open(destination, "w", newline='', encoding="utf8") as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            writer.writerow([overall_title])
            writer.writerow(field_names)
            for line in sim_listing:
                writer.writerow(line)
    except OSError as err:
        fileop.io_snag(err, sub_path, file_name)

    # return None


def cache_sim_results(sim_numbr, sim_listing, entitle, prevent):
    """Control the data storage of an entire simulation.

    Parameters
    ----------
    sim_numbr = (int) realization number
    sim_listing = (list of lists) results for a realization
    entitle = (str) title for run
    prevent = (bool) flag to not print

    Returns
    -------
    N/A
    """
    # Check if results should be printed.
    if not prevent:
        # Construct unique file name for saving a realization.
        current_sim = scfg.RESULTS_NAME + str(sim_numbr + 1)
        output_nomer = current_sim + scfg.EXTENSION_CSV

        # Save List of lists for time-history.
        save_sim_list(output_nomer, sim_listing, entitle)

    # return None


def input_apertures(fault_controls, param_bounds):
    """Acquire the aperture values from a file and checks them.

    Parameters
    ----------
    fault_controls = (dict) input parameters (dictionary)
    param_bounds = (dict) parameter bounds (dictionary)

    Returns
    -------
    fault_controls = (dict) input parameters (dictionary)
    """
    # Import data, if desired.
    if fault_controls['aperture_approach']:
        if ECHO:
            print("  --> Input from file: '" + scfg.APERTURE_NAME + "'")

        # Get data array from file.
        data_array = fileop.acquire_data_array(scfg.INPUT_DIR,
                                               scfg.APERTURE_NAME,
                                               fault_controls['n_plates'])

        # Check if numeric values are within limits.
        limiter = param_bounds['aperture_confines']
        stats = fileop.check_file_data(data_array, "Aperture value", limiter)

        # If error (statx < 0), halt program.
        err_flag = 0  # Added to eliminate inspection error on stats (why?).
        err_flag += stats
        if err_flag < 0:
            fileop.opx_problem("Error: Problem with Aperture Data " +
                               "Caused Program to Exit!")

        # Define parameters for next cycle as prior values.
        fault_controls['aperture_data'] = data_array

    return fault_controls


def get_data_lines(destination, file_location, n_plates, skip_lines):
    """Get a block of data from the repository file.

    Parameters
    ----------
    destination = (str) path to input file
    file_location = (list) subdirectory & file names
        [0] -> subdirectory
        [1] -> filename
    n_plate = (int) number of fault plates (columns)
    skip_lines = (int) present position count in data file

    Returns
    -------
    new_array = (NumPy array) data values
    """
    # Define list and adjust skip for single header line.
    data_list = []
    skip_lines += 1

    try:
        # Read data after skipping previous data.
        with open(destination, encoding="utf8") as line_input:
            # Skip to right position.
            for _ in repeat(None, skip_lines):
                line_input.readline()

            # Read one line/row of data.
            lined = line_input.readline()
            interim = lined.strip().split(',')
            divided = interim[0:n_plates]
            data_list.append(divided)

    except OSError as err:
        fileop.io_snag(err, file_location[0], file_location[1])

    # Convert to numpy array - possible error in reading.
    try:
        new_array = np.asarray(data_list, dtype=np.float64)
    except ValueError:
        # Value error found - exit.
        fileop.data_error("Value Error! -> Check Data Points in File!",
                          file_location[0], file_location[1])
        new_array = np.empty(0)  # for inspection

    # Ensure array is 1D - issue w. flatten.
    result = new_array.flatten('C')

    # Error check on shape >> dimension = n_plates.
    if result.size != n_plates:
        msg1 = "Error in Reservoir Data! "
        msg1 += f"-> Wrong Number of Data Points in {file_location[1]}."
        fileop.opx_problem(msg1)

    return result


def input_reservoir_data(sequence, n_plates, param_bounds):
    """Control the data input from reservoir files at a time point.

    Parameters
    ----------
    sequence = (int) current time point index (starts at 1)
    n_plate = (int) number of fault plates
    param_bounds = (dict) data bounds (dictionary)

    Returns
    -------
    pressure = (NumPy array) base CO2 pressure
    saturation = (NumPy array) base CO2 saturation

    Notes
    -----
    1. Files assumed to have extension *.txt,
        and are comma delimited, with NO header.
    2. Pressure and Saturation at start of time step are assumed to be
        constant for entire period - no averaging employed in this version.
    """
    # Setup data arrays.
    #    skip_lines -> position in file data.
    # saturation = np.zeros(n_plates)
    # pressure = np.zeros(n_plates)
    skip_lines = (sequence - 1)  # sequence starts at 1

    # Saturation - Open data array and read lines until start.
    if ECHO:
        print("  --> Input from File: '" + scfg.SATURATION_NAME + "'")

    sub_path, destination = \
        fileop.get_path_name(scfg.RESERVOIR_DIR, scfg.SATURATION_NAME,
                             scfg.EXTENSION_TXT)
    location = [sub_path, scfg.SATURATION_NAME]
    saturation = get_data_lines(destination, location, n_plates, skip_lines)

    # Reservoir - Open data array and read lines until start.
    if ECHO:
        print("  --> Input from File: '" + scfg.PRESSURE_NAME + "'")

    sub_path, destination = fileop.get_path_name(scfg.RESERVOIR_DIR,
                                                 scfg.PRESSURE_NAME,
                                                 scfg.EXTENSION_TXT)
    location = [sub_path, scfg.PRESSURE_NAME]
    pressure = get_data_lines(destination, location, n_plates, skip_lines)

    # Check that the 2 data arrays are within defined limits.
    err_flag = 0
    val_limits = param_bounds['reservoir_pressure']
    stats = fileop.check_file_data(pressure, "Reservoir Pressure Value",
                                   val_limits)
    err_flag += stats

    val_limits = param_bounds['reservoir_saturation']
    stats = fileop.check_file_data(saturation, "Reservoir Saturation Value",
                                   val_limits)
    err_flag += stats

    # If error, halt program.
    if err_flag < 0:
        fileop.opx_problem("Error: File value(s) Outside Defined Bounds " +
                           "Caused Program to Exit!")

    return pressure, saturation


def list_data_head():
    """Construct the header for instance properties of plate data.

    Parameters
    ----------
    N/A

    Returns
    -------
    txtr = (str) Header data string for output
    """
    # Construct header for position values.
    txtr = ('No.', 'Start X', 'Start Y', 'End X', 'End Y', 'Strike', 'Dip',
            'Aperture', 'Perm.', 'Threshold')

    return txtr


#
# -----------------------------------------------------------------------------
# End of module
