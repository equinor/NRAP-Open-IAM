#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions to get and define data from file input.

Author: Ernest N. Lindner
Date: 08/18/2022

Module Name
    seal_upload

Contents (15)
    random_seeder(switch=0)
    define_units_for_output(si_units)
    define_file_limits(param_bounds)
    define_influence_factor(seal_controls, grid, file_limits)
    define_coordinates(grid, param_bounds)
    define_thickness(seal_controls, grid, file_limits)
    define_active(seal_controls, grid, file_limits)
    define_area(seal_controls, grid, file_limits)
    define_depth(seal_controls, grid, file_limits)
    define_entry(seal_controls, grid, file_limits)
    ----
    input_various_data(input_controls, grid, param_bounds, rng)
    define_top_press_array(seal_controls, grid, param_bounds)
    get_data_lines(destination, file_location, lines, columns, skip_lines)
    reservoir_check(seal_controls)
    input_reservoir_data(sequence, seal_controls, param_bounds)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
from itertools import repeat        # For line input
import numpy as np                  # For arrays & input

import seal_config as scfg           # IO directory and file names
import seal_file as sfile           # For file operations
import seal_refresh as sref         # Functions for thickness and depth

# Other constants
ECHO = False                            # DEBUG printout - status of input.
STARTER = "  --> Input from file: "     # Line start


def random_seeder(switch=0):
    """Define random function seed.

    Parameters
    ----------
    switch = (int) flag to increase SEEDX

    Returns
    -------
    rng = (generator) random number generator
    """
    # Set seed for random numbers.
    if scfg.SEEDX is None:
        rng = np.random.default_rng()
    else:
        if switch > 0:
            scfg.SEEDX += 10
        randy = scfg.SEEDX
        rng = np.random.default_rng(randy)

    return rng


def define_units_for_output(si_units):
    """Define units for output - as strings.

    Parameters
    ----------
    si_units = (bool) control to convert English to SI units

    Returns
    -------
    uts = (list of str) parameter units for printing
    """
    if si_units:
        # Values are in SI; define units for printout.
        # set unit labels for printout in SI.
        uts = ['m', 'MPa', 'kg/m3', 'Pa-s', 'mol/kg', 'm2', 'metric tonne',
               'oC', 'mm']
    else:
        # values are in English/US units; define units for printout.
        uts = ['ft', 'psi', 'pcf', 'cP', 'mol/lb', 'ft2', 'short ton',
               'oF', 'ín']

    return uts


def define_file_limits(param_bounds):
    """Define numeric ranges of additional cell parameters.

    Parameters
    ----------
    param_bounds = (dict) min/max bounds for Yaml file

    Returns
    -------
    param_bounds2 = (dict) additional min/max limits for file input

    Notes
    -----
    Function modifies prior dictionary bounds for use with files.
    """
    # Define dictionary of input bounds.
    param_bounds2 = {}

    # Define area limits.
    param_bounds2['area'] = [param_bounds['area'][0],
                             param_bounds['area'][1]]

    # Define thickness limits.
    param_bounds2['thickness'] = [param_bounds['thickness_min'][0],
                                  param_bounds['thickness_max'][1]]

    # Define active cell limits.
    param_bounds2['status'] = [0, 1]

    # Define entry pressure limits.
    param_bounds2['entry'] = [param_bounds['entry_pressure'][0],
                              param_bounds['entry_pressure'][1]]

    # Define limits of depth to top of cell.
    param_bounds2['top_depth'] = [param_bounds['static_depth'][0],
                                  param_bounds['static_depth'][1]]

    return param_bounds2


def define_influence_factor(seal_controls, grid, file_limits):
    """Define influence model values for each cell.

    Parameters
    ----------
    seal_controls = (dict) seal parameters
    grid = (list of class) collection of cells
    file_limits = (dict) minimum / maximum bounds for input

    Returns
    -------
    statx = error flag
    """
    if seal_controls['initialize']:
        # Echo operation to user.
        if ECHO:
            print("   >> Input from file: '" + scfg.INFLUENCE_NAME + "'")

        # Get data array - ".txt" added by call.
        data_array = sfile.acquire_data_array(scfg.INPUT_DIRECTORY,
                                              scfg.INFLUENCE_NAME,
                                              len(grid))

        # Check if numeric values are within limits.
        val_limits = file_limits['influence']
        statx = sfile.check_file_floats(data_array, "Influence Factor value",
                                        val_limits)

        # Update influence factor for each cell.
        for numbr, cell in enumerate(grid):
            cell.influence = data_array[numbr]

    else:
        # Set default influence factor for each cell.
        statx = 0
        for cell in grid:
            cell.influence = seal_controls['influence']

    return statx


def define_coordinates(grid, param_bounds):
    """Input coordinate data from file.

    Parameters
    ----------
    grid = (list of class) cells
    param_bounds = (dict) seal parameter bounds

    Returns
    -------
    grid = (list of class) cells
    """
    # Write status of input process.
    if ECHO:
        print("   >> Input from Files: '" + scfg.X_COORD_NAME + "' and '"
              + scfg.Y_COORD_NAME + "' ")

    # Get data arrays - ".txt" extension is added in call.
    x_data_array = sfile.acquire_data_array(scfg.INPUT_DIRECTORY,
                                            scfg.X_COORD_NAME,
                                            len(grid))
    y_data_array = sfile.acquire_data_array(scfg.INPUT_DIRECTORY,
                                            scfg.Y_COORD_NAME,
                                            len(grid))

    # Check if numeric values are within typical limits.
    val_limits = param_bounds['coordinates']
    statx = sfile.check_file_floats(x_data_array, "X-Coordinate value",
                                    val_limits)
    statx += sfile.check_file_floats(y_data_array, "Y-Coordinate value",
                                     val_limits)
    if statx < 0:
        sfile.opx_problem("File Value(s) Outside Defined Bounds " +
                          "Caused Program to Exit!")

    # Update cell data for coordinate input.
    for numbr, cell in enumerate(grid):
        cell.set_coord(x_data_array[numbr], y_data_array[numbr])

    return grid


def define_thickness(seal_controls, grid, file_limits, rng):
    """Define thickness values for each cell.

    Parameters
    ----------
    seal_controls = (dict) seal parameters
    grid = (list of class) collection of cells
    file_limits = (dict) minimum / maximum bounds for input
    rng = (generator) random number generator

    Returns
    -------
    statx = error flag
    """
    if seal_controls['thickness_approach']:
        # echo operation to user.
        if ECHO:
            print(STARTER + scfg.THICKNESS_NAME + "'")

        # Define data array from file.
        new_thickness_array = sfile.acquire_data_array(scfg.INPUT_DIRECTORY,
                                                       scfg.THICKNESS_NAME,
                                                       len(grid))

        # Check if numeric values are within limits.
        val_limits = file_limits['thickness']
        statx = sfile.check_file_floats(new_thickness_array, "Thickness value",
                                        val_limits)
    else:
        # Generate thickness data array based on stochastic values (1D array).
        statx = 0
        new_thickness_array = sref.thickness_variability(seal_controls, rng)

    # Set cell thickness value of grid from 1D array.
    for numbr, cell in enumerate(grid):
        cell.thickness = new_thickness_array[numbr]

    return statx


def define_active(seal_controls, grid, file_limits):
    """Define active status values for each cell.

    Parameters
    ----------
    seal_controls = (dict) seal parameters
    grid = (list of class) collection of cells
    file_limits = (dict) minimum / maximum bounds for input

    Returns
    -------
    statx = error flag; 0 = no error
    """
    if seal_controls['active_approach']:
        # echo operation to user.
        if ECHO:
            print(STARTER + scfg.ACTIVE_NAME + "'")

        # Get data array from file.
        data_array = sfile.acquire_data_array(scfg.INPUT_DIRECTORY,
                                              scfg.ACTIVE_NAME,
                                              len(grid))

        # Check if numeric values are within limits.
        val_limits = file_limits['status']
        statx = sfile.check_file_values(data_array, "Cell Status value",
                                        val_limits)

        # Update cell values.
        for numbr, cell in enumerate(grid):
            cell.status = int(data_array[numbr])

    else:
        # Set cells active if not input from file.
        for cell in grid:
            cell.status = int(1)
        statx = 0

    return statx


def define_area(seal_controls, grid, file_limits):
    """Define area values for each cell.

    Parameters
    ----------
    seal_controls = (dict) seal parameters
    grid = (list of class) collection of cells
    file_limits = (dict) minimum / maximum bounds for input

    Returns
    -------
    statx = error flag; 0 = no error
    """
    if seal_controls['area_approach']:
        # echo operation to user.
        if ECHO:
            print(STARTER + scfg.AREA_NAME + "'")

        # Get data array from file.
        data_array = sfile.acquire_data_array(scfg.INPUT_DIRECTORY,
                                              scfg.AREA_NAME,
                                              len(grid))

        # Check if numeric values are within limits.
        val_limits = file_limits['area']
        statx = sfile.check_file_floats(data_array, "Area value",
                                        val_limits)

        # Update cell values.
        for numbr, cell in enumerate(grid):
            cell.area = data_array[numbr]

    else:
        # Set cells to default area parameter.
        statx = 0
        for cell in grid:
            cell.area = seal_controls['area']

    return statx


def define_depth(seal_controls, grid, file_limits):
    """Define depth values for each cell.

    Parameters
    ----------
    seal_controls = (dict) seal parameters
    grid = (list of class) collection of cells
    file_limits = (dict) minimum / maximum bounds for input

    Returns
    -------
    statx = (int) error flag; 0 = no error
    """
    if seal_controls['depth_approach']:
        # echo operation to user.
        if ECHO:
            print(STARTER + scfg.BASE_NAME + "'")

        # Get data array from file.
        data_array = sfile.acquire_data_array(scfg.INPUT_DIRECTORY,
                                              scfg.BASE_NAME,
                                              len(grid))

        # Check if numeric values are within limits.
        val_limits = file_limits['top_depth']
        statx = sfile.check_file_floats(data_array, "Depth value",
                                        val_limits)

        # Update cell values.
        for numbr, cell in enumerate(grid):
            cell.set_top_depth(data_array[numbr])

    else:
        # Set cells to default using base depth parameter for calculation.
        statx = 0
        for cell in grid:
            if cell.status > 0:
                cell.set_top_depth(seal_controls['ave_base_depth'])

    return statx


def define_entry(seal_controls, grid, file_limits):
    """Define threshold pressure values for each cell.

    Parameters
    ----------
    seal_controls = (dict) seal parameters
    grid = (list of class) collection of cells
    file_limits = (dict) minimum / maximum bounds for input

    Returns
    -------
    statx = (int) error flag; 0 = no error
    """
    if seal_controls['entry_approach']:
        # Get data array from file.
        if ECHO:
            print(STARTER + scfg.ENTRY_NAME + "'")
        data_array = sfile.acquire_data_array(scfg.INPUT_DIRECTORY,
                                              scfg.ENTRY_NAME,
                                              len(grid))

        # Check if numeric values are within limits.
        val_limits = file_limits['entry']
        statx = sfile.check_file_floats(data_array, "Entry pressure value",
                                        val_limits)

        # Update cell values.
        for numbr, cell in enumerate(grid):
            cell.entry = data_array[numbr]

    else:
        # Set cells to default entry pressure parameter.
        statx = 0
        for cell in grid:
            cell.entry = seal_controls['entry_pressure']

    return statx


def input_various_data(seal_controls, grid, param_bounds, rng):
    """Manage/check input of values for various data types.

    Parameters
    ----------
    seal_controls = (dict) seal controls
    grid = (list of class) collection of cells
    file_limits = (dict) minimum / maximum bounds for input
    rng = random number function

    Returns
    -------
    grid = (list) cells (class) - updated

    Notes
    -----
    input_controls: index for cell aspects, all are *.txt files
    """
    # Set error flag and define limits.
    err_flag = 0
    file_limits = define_file_limits(param_bounds)

    # 1. Initialize influence factor, as desired.
    result = define_influence_factor(seal_controls, grid, file_limits)
    err_flag += result

    # 2. Define thickness of grid.
    result = define_thickness(seal_controls, grid, file_limits, rng)
    err_flag += result

    # 3. Define active cells in grid; default is active.
    result = define_active(seal_controls, grid, file_limits)
    err_flag += result

    # 4. Define area of grid.
    result = define_area(seal_controls, grid, file_limits)
    err_flag += result

    # 5. Define depth to top of cell.
    result = define_depth(seal_controls, grid, file_limits)
    err_flag += result

    # 6. Define entry pressure of grid.
    result = define_entry(seal_controls, grid, file_limits)
    err_flag += result

    # Report error if found. Data outside bounds is deemed fatal!
    if err_flag < 0:
        sfile.opx_problem("File value(s) Outside Defined Bounds " +
                          "Caused Program to Exit!")

    return grid


def define_top_press_array(seal_controls, grid, param_bounds):
    """Define top pressure values.

    Parameters
    ----------
    seal_controls = (dict) seal parameters
    grid = (list of class) collection of cells
    param_bounds = (dict) minimum / maximum bounds for input

    Returns
    -------
    press_top = (1D array) Pressure at top of seal cells
    """
    # Get data depending on user choice.
    dimension = seal_controls['num_cells']
    press_top = np.zeros(dimension)

    if seal_controls['upper_approach']:
        # Input top pressures from file.
        if ECHO:
            print(STARTER + scfg.PRESS_NAME + "'")
        press_top = sfile.acquire_data_array(scfg.INPUT_DIRECTORY,
                                             scfg.PRESS_NAME,
                                             len(grid))

        # Check if numeric values are within typical pressure bounds.
        # - Use base pressure as representative.
        val_limits = param_bounds['ave_base_pressure']
        statx = sfile.check_file_floats(press_top, "Top Pressure Value",
                                        val_limits)

        if statx < 0:
            sfile.opx_problem("File value(s) Outside Defined Bounds " +
                              "Caused Program to Exit!")

    else:
        # Compute static pressure from reference depth and reference pressure.
        for numbr, cell in enumerate(grid):
            press_top[numbr] = \
                cell.compute_static_pressure(seal_controls['static_pressure'],
                                             seal_controls['static_depth'])

    return press_top


def get_data_lines(destination, file_location, lines, columns, skip_lines):
    """Read a block of data from reservoir file.

    Parameters
    ----------
    destination = (str) path to input file
    file_location = (list of str) subdirectory & file
          =0 -> subdirectory
          =1 -> filename
    lines = (int) lines of input
    columns = (int) columns of input
    skip_lines = (int) present position count

    Returns
    -------
    new_array = (array) pressure or saturation values in (1D) format
    """
    # Define list and add 1 to location for header line.
    data_list = []
    skip_lines += 1
    new_array = None

    try:
        # Read data after skipping previous data + format line.
        with open(destination, encoding="utf8") as line_input:
            # Skip to right position.
            for _ in repeat(None, skip_lines):
                line_input.readline()

            # Read data block as list.
            for _ in repeat(None, lines):
                lined = line_input.readline()
                interim = lined.strip().split(',')
                divided = interim[0:columns]
                data_list.append(divided)

    except OSError as err:
        sfile.io_snag(err, file_location[0], file_location[1])

    # Convert to array - possible error in reading.
    try:
        new_array = np.asarray(data_list, dtype=np.float64)
    except ValueError:
        # Value error found -> exit program.
        sfile.data_error("Value Error! -> Check Data Format "
                         + "in file!", file_location[0], file_location[1])

    # Error check on shape of 2D array.
    if new_array.shape != (lines, columns):
        msg1 = "Failure in Reservoir Data Input,"
        msg1 += "- Wrong Number of Lines/Columns!"
        sfile.opx_problem(msg1)

    # Flatten arrays for use in computations.
    new_array.shape = (lines*columns)

    return new_array


def reservoir_check(seal_controls):
    """Read the line/column format at start of the reservoir file.

    Parameters
    ----------
    seal_controls = (dict) seal parameters

    Returns
    -------
    lines = number of lines in reservoir data block
    columns = number of columns in reservoir data block

    Notes
    -----
    First entry in reservoir file must be line, column specification.

    """
    # Get filename and path.
    sub_path, destination = \
        sfile.get_path_name(scfg.RESERVOIR_DIR, scfg.PRESSURE_NAME,
                            scfg.EXTENSION_TXT)
    file_location = [sub_path, scfg.PRESSURE_NAME]

    try:
        # Read line/columns data from first line.
        with open(destination, encoding="utf8") as line_input:
            start = line_input.readline()
            item = start.strip().split(',')
            lines = int(item[0])
            columns = int(item[1])

    except OSError as err:
        sfile.io_snag(err, file_location[0], file_location[1])

    # Check that there is sufficient data for cells based on format.
    if (lines * columns) < seal_controls['num_cells']:
        sfile.opx_problem("Reservoir File has Insufficient Data " +
                          "for Specified Number of Cells!")

    return lines, columns


def input_reservoir_data(sequence, seal_controls, param_bounds):
    """Obtain data from reservoir files at a time point.

    Parameters
    ----------
    sequence = (float) current time point index (starts at 1)
    seal_controls = (dict) seal parameters
    param_bounds = (dict) data bounds

    Returns
    -------
    pressure = (array) base CO2 pressure - as NumPy 1D array
    saturation = (array) base CO2 saturation - as 1D NumPy array

    Notes
    -----
    1. Files assumed to have extension *.txt,
        and are comma delimited, with NO header.
    2. Pressure and Saturation at start of time step are assumed to be
        constant for entire period - no averaging employed in this version.
    """
    # Setup data arrays.
    lines = seal_controls['r_lines']
    columns = seal_controls['r_columns']
    skip_lines = lines * (sequence - 1)

    # Saturation - Open data array and read lines, out => 1D array.
    if ECHO:
        print(STARTER + scfg.SATURATION_NAME + "'")
    sub_path, destination = \
        sfile.get_path_name(scfg.RESERVOIR_DIR, scfg.SATURATION_NAME,
                            scfg.EXTENSION_TXT)
    location = [sub_path, scfg.SATURATION_NAME]
    saturation = get_data_lines(destination, location, lines, columns,
                                skip_lines)

    # Reservoir - Open data array and read lines, out => 1D array.
    if ECHO:
        print(STARTER + scfg.PRESSURE_NAME + "'")
    sub_path, destination = \
        sfile.get_path_name(scfg.RESERVOIR_DIR, scfg.PRESSURE_NAME,
                            scfg.EXTENSION_TXT)
    location = [sub_path, scfg.PRESSURE_NAME]
    pressure = get_data_lines(destination, location, lines, columns,
                              skip_lines)

    # Check that the data are within limits.
    err_flag = 0
    val_limits = param_bounds['reservoir_pressure']
    statx = sfile.check_file_floats(pressure, "Reservoir pressure value",
                                    val_limits)
    err_flag += statx

    val_limits = param_bounds['reservoir_saturation']
    statx = sfile.check_file_floats(saturation, "Reservoir saturation value",
                                    val_limits)
    err_flag += statx

    # If error, halt program.
    if err_flag < 0:
        sfile.opx_problem("File Value(s) Outside Defined Bounds " +
                          "Caused Program to Exit!")

    return pressure, saturation


#
# -----------------------------------------------------------------------------
# - End of module
