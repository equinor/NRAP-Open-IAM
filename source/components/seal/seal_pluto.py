#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains auxiliary functions to plot data.

Author: Ernest N. Lindner
Date: 08/18/2022

Module Name
    seal_pluto

Contents (17)
    line_space(start_value, end_value, n)
    set_exponent(z_data)
    test_contour(z_data)
    check_plot_data(z_data)
    construct_coordinate_data(y_vals, x_vals, grid)
    construct_permeability_data(y_vals, x_vals, grid)
    construct_straight_coords(grid, data_size)
    construct_irregular_permeability(grid)
    irregular_contour_data(grid, z_data)
    set_bar_arrays(y_vals, x_vals)
    ----
    round_exp(exp_value)
    setup_time_history_data(simulation, si_units)
    setup_co2_seepage_data(simulation):
    simulation_query(max_realizations, start_value=0, gapper=MIN_IN)
    create_history_query(query_stage, series)
    create_contour_query(trial)
    range_query(max_realizations)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import math                                     # For truncation
import copy                                     # For array adjustment
from scipy.interpolate import griddata          # For irregular points
import numpy as np                              # For array operations

import seal_config as scfg                       # IO directory names
import seal_file as sfile                       # For paths and error reports
import seal_units as sun                        # Unit conversion

# Constants for starting for statements
LED_IN = "      --> "                           # General message start
CMD_IN = "      >>> "                           # Command message start
ERR_IN = "      *** "                           # Error message start
QUE_IN = "      ?-> "                           # Question on console
QU2_IN = "\n      ?-> "                         # Repeat question
MIN_IN = "      ?-> "                           # Second question on console

# Controls and warnings
NEGATIVES = ('n', 'N', 'q', 'Q')                # Keys for negative answer
EXIT_COMMAND = ('x', 'X')                       # Keys that will exit
DOUBLES = ('nn', 'NN', 'qq', 'QQ', 'xx', 'XX', 'yy', 'YY')  # Common errors
POSITIVES = ('y', 'YY')
WAIT_FOR_IT = CMD_IN + "Please wait for plot ..."
AGAIN = "\n" + CMD_IN + "Please Try Again."
INPUT = ERR_IN + "You typed = {}"
BAR_SHIFT_X = 0.75          # Shift of 3D bars along X axis
BAR_SHIFT_Y = 0.50          # Shift of 3D bars along Y axis
NUM_POINTS_X = 100          # Number of x points for irregular plot
NUM_POINTS_Y = 100          # Number of y points for irregular plot


def line_space(start_value, end_value, num):
    """Compute a np array of equal spacing - a replacement for Numpy linspace.

    Parameters
    ----------
    start_value = (float) start value of array
    end_value = (float) end value of array
    num = (int) number of values

    Returns
    -------
    out_array = (numpy array) spaced array of values (1D)
    """
    # If num => then setup array; else default.
    if num >= 2:
        diff = (float(end_value) - start_value)/(num - 1)
        values = [diff * i + start_value for i in range(num)]
    else:
        values = end_value

    # Convert list to array.
    out_array = np.array(values)

    return out_array


def set_exponent(z_data):
    """Find exponent to fit data and then factor data array.

    Parameters
    ----------
    z_data = (array) data to plot - NumPy array

    Returns
    -------
    cite = (int) exponent for data
    show_units = (str) text for data magnitude
    zed_data = (array)
    """
    # Get exponent value for z data - Don't change source!
    z_max = np.max(z_data)
    cite = math.trunc(math.log10(abs(z_max)))
    if z_max < 1.0:
        cite -= 1
    zed_data = copy.deepcopy(z_data)
    zed_data /= pow(10.0, cite)
    show_units = r'$\times$10$^{%d}$' % cite

    return cite, show_units, zed_data


def test_contour(z_data):
    """Examine data to prevent a poor contour plot.

    Parameters
    ----------
    z_data = (array) data to plot - NumPy array

    Returns
    -------
    answer = (bool): True = data OK; False = error / do not plot
    """
    # Examine if spread in data - check min/max.
    answer = True
    z_min = np.min(z_data)
    z_max = np.max(z_data)

    if z_max == z_min:
        # No contour possible!
        print(CMD_IN + "Uniform Value Across Grid of " + str(z_max)
              + " Values!")
        print(CMD_IN + "Data is Non-Differentiated! " +
              "Contour Plot is Not Shown!")
        answer = False

    return answer


def check_plot_data(z_data):
    """Examine values to prevent an error in time series plot.

    Parameters
    ----------
    z_data = (array) data to plot - NumPy array

    Returns
    -------
    answer = (bool): True = data OK; False = error / do not plot
    """
    # Examine if spread in data - check min/max.
    answer = True
    z_min = np.min(z_data)
    z_max = np.max(z_data)

    if z_max <= 0.0 or z_min == z_max:
        # No plot possible!
        answer = False

    return answer


def construct_coordinate_data(rank, chain, grid):
    """Cast grid x/y coordinate data into NumPy coordinate arrays.

    Parameters
    ----------
    rank = (int) lines in plot
    chain = (int) columns in plot
    grid = (list) (class) a collection of cells (1D)

    Returns
    -------
    x_data, y_data = (array) coordinates for cell centers (2D)
    """
    # Establish x/y data for contour plots.
    x_data = np.zeros((rank, chain))
    y_data = np.zeros((rank, chain))

    for i in range(rank):
        for j in range(chain):
            indx = j + (i * chain)
            x_data[i, j] = grid[indx].x_center
            y_data[i, j] = grid[indx].y_center

    return x_data, y_data


def construct_permeability_data(rank, chain, grid):
    """Cast the permeability data into NumPy array.

    Parameters
    ----------
    rank = (int) lines in plot
    chain= (int) columns in plot
    grid = (list) (class) a collection of cells

    Returns
    -------
    z_data = (array) data to plot - NumPy array
    """
    # Create permeability data for plot.
    z_data = np.zeros((rank, chain))
    for i in range(rank):
        for j in range(chain):
            indx = j + (i * chain)
            # report only permeability for active cells.
            if grid[indx].status > 0:
                z_data[i, j] = grid[indx].permeability
            else:
                z_data[i, j] = 0.0

    return z_data


def construct_straight_coords(grid):
    """Cast irregular x/y coordinate data into NumPy coordinate arrays.

    Parameters
    ----------
    grid = (list) (class) a collection of cells (1D)

    Returns
    -------
    x_data, y_data = (array) coordinates for cell centers (1D)
    """
    # Establish x/y data for contour plots.
    data_size = len(grid)
    x_data = np.zeros(data_size)
    y_data = np.zeros(data_size)

    for indx in range(data_size):
        x_data[indx] = grid[indx].x_center
        y_data[indx] = grid[indx].y_center

    return x_data, y_data


def construct_irregular_permeability(grid):
    """Cast the permeability data into NumPy array for irregular contour.

    Parameters
    ----------
    grid = (list) (class) a collection of cells

    Returns
    -------
    z_data = (array) data to plot - NumPy array (1D)
    """
    # Create permeability data for irregular plot.
    data_size = len(grid)
    z_data = np.zeros(data_size)

    for indx in range(data_size):
        # report only permeability for active cells.
        if grid[indx].status > 0:
            z_data[indx] = grid[indx].permeability
        else:
            z_data[indx] = 0.0

    return z_data


def construct_irregular_data(grid, z_data):
    """Determine data arrays for contouring irregular spaced perm values.

    Parameters
    ----------
    grid = (list) (class) a collection of cells
    z_data - (array) data to be plotted (1D) array

    Returns
    -------
    x_new y_new = (array) regularly-spaced coordinates
    z_new = (array) regularly-spaced permeability values
    result = (bool) condition of array
              = True => array has size
              = False => array has no coordinates!

    Notes
    -----
    1. User may have forgotten to input coordinates needed here.
    """
    # Construct x/y 1D arrays from cells.
    x_data, y_data = construct_straight_coords(grid)

    # Find min/max of x/y coordinates to define plot extent.
    x_min = np.min(x_data)
    x_max = np.max(x_data)
    y_min = np.min(y_data)
    y_max = np.max(y_data)

    # Check if grid exists.
    if (x_min == x_max) or (y_min == y_max):
        # No coordinates for plot!!
        result = False

        # Create dummy arrays to prevent error on return.
        x_new = [0]
        y_new = [0]
        z_new = [0]

    else:
        result = True
        # Generate a regular x/y mesh to interpolate the data.
        x_line = np.linspace(x_min, x_max, num=NUM_POINTS_X, endpoint=True)
        y_line = np.linspace(y_min, y_max, num=NUM_POINTS_Y, endpoint=True)
        x_new, y_new = np.meshgrid(x_line, y_line)

        # Interpolate to find z values, using Delaunay triangularization.
        z_new = griddata((x_data, y_data), z_data, (x_new, y_new),
                         method='cubic')

    return x_new, y_new, z_new, result


def set_bar_arrays(y_vals, x_vals):
    """Define coordinates for 3D bar graph.

    Parameters
    ----------
    y_vals = (int) lines in plot
    x_vals = (int) columns in plot

    Returns
    -------
    x_loca = (array) flat array of x coordinates (1D)
    y_loca = (array) flat array of y coordinates (1D)
    """
    # create axes.
    x_indices = np.arange(x_vals)
    y_indices = np.arange(y_vals)

    # Create coordinate data for plot - shift for better image.
    x_coord, y_coord = np.meshgrid(BAR_SHIFT_X + x_indices,
                                   BAR_SHIFT_Y + y_indices)

    # Flatten data to 1D arrays.
    x_loca = x_coord.flatten()
    y_loca = y_coord.flatten()

    return x_loca, y_loca


def round_exp(exp_value):
    """Round-up exponential to desired degree.

    Parameters
    ----------
    exp_value = (float) general limit value

    Returns
    -------
    target = (float) rounded exponential
    """
    # Trap error if = 0.0.
    if exp_value != 0.0:
        # Make the number a decimal.
        core = abs(exp_value)
        level = math.trunc(math.log10(core))
        base = math.pow(10.0, level)
        root_exp = math.ceil(exp_value / base)

        # Round and then reconstitute the number.
        # -- adjust_number = round(root_exp, 1)
        target = root_exp * base  # round to zero decimal
    else:
        target = exp_value

    return target


def setup_time_history_data(simulation, use_si_units):
    """Get time history data from file.

    Parameters
    ----------
    simulation = (int) simulation number to plot
    use_si_units = (bool) use of metric units

    Returns
    -------
    leakage_data = (array) of co2 data - numpy
    time_data = (array) of time values - numpy
    """
    # Construct path name to file.
    sim_number = str(simulation)
    file_name = scfg.RESULTS_FILE + sim_number
    subdirectory_path, destination = sfile.get_path_name(scfg.OUTPUT_DIRECTORY,
                                                         file_name,
                                                         scfg.EXTENSION_CSV)
    # Get data list from destination.
    data_list = sfile.acquire_csv_list(destination, subdirectory_path,
                                       file_name)

    # Convert list to NumPy array.
    npa = np.asarray(data_list, dtype=np.float64)

    # Slice array for data.
    time_data = npa[:, 1]
    leakage_data = npa[:, 2]

    # Recast leakage units into English/US units, if desired.
    if not use_si_units:
        leakage_data = sun.tonne_to_ton() * leakage_data

    return leakage_data, time_data


def setup_co2_seepage_data(simulation):
    """Get the CO2 seepage data from NPY file.

    Parameters
    ----------
    simulation = (int) simulation number to plot

    Returns
    -------
    leakage_data = (array) CO2 data - NumPy
    """
    # Construct path name to file.
    sim_number = str(simulation)
    file_name = scfg.CO2_INTRO + sim_number
    subdirectory_path, destination = sfile.get_path_name(scfg.OUTPUT_DIRECTORY,
                                                         file_name,
                                                         scfg.EXTENSION_NPY)
    # To address Pycharm error report.
    leakage_data = 0.0

    try:
        # Load NumPy file and return.
        leakage_data = np.load(destination)

    except OSError as err_x:
        sfile.io_snag(err_x, subdirectory_path, file_name)

    return leakage_data


def simulation_query(max_realizations, start_value=0, gapper=MIN_IN):
    """Get simulation number to plot.

    Parameters
    ----------
    max_realizations = (int) maximum number of realizations
    start_value = (int) minimum/start number of realizations to be plotted
    gapper = (str) indent on question for Enter Realization

    Returns
    -------
    response = status for next plot:
        True = process plot
        False = exit
    realization = (int) simulation number to plot
    """
    # Default values.
    repeat_loop = True
    response = False
    realization = -1

    # Loop to get an appropriate response.
    while repeat_loop:
        try:
            # Get input with prompt.
            code = input(gapper + "Enter Realization Number "
                         + "to Plot ('Q'=quit): ")
            # Check for double letters.
            if code in DOUBLES:
                code = code[:1]
            # Check if user wishes to quit.
            if code in NEGATIVES or code in EXIT_COMMAND:
                print(QUE_IN + "Exiting Plot Option. \n")
                repeat_loop = False
                break
        except ValueError:
            # Parse fails.
            print(ERR_IN + "Invalid Input!  " + AGAIN)
            continue

        # Check if "entered" is a number.
        try:
            entered = int(code)
        except ValueError:
            # Not an integer.
            print(ERR_IN + "Input Is Not A Number!  " + AGAIN)
            continue

        # Check range of number
        if entered > (max_realizations + 1):
            print(INPUT.format(entered) + "\n"
                  + ERR_IN + "This Number Exceeds the Maximum "
                  + f' of {max_realizations}!', end='')
            print(AGAIN)
        elif entered <= 0:
            print(INPUT.format(entered) + "\n"
                  + ERR_IN + "The Input is Less than Zero!", end='')
            print(AGAIN)
        elif entered < start_value:
            print(INPUT.format(entered) + "\n"
                  + ERR_IN + "The Input is Less Than Starting Value"
                  + f' of {start_value}!', end='')
            print(AGAIN)
        else:
            # Data OK.
            realization = entered
            response = True
            # repeat_loop = False.
            break
    # end while

    return response, realization


def create_history_query(query_stage, series):
    """Check if user wants to plot a time series.

    Parameters
    ----------
    query_stage = (int) question queue
            0 = first time
            1 = plot again?
    Series = (bool) Type of plot
            True = series of time plots
            False = single time history

    Returns
    -------
    reply = (bool) answer plot question
    """
    # Change message depending on place in query queue.
    if series:
        # Time series plot.
        if query_stage == 0:
            code = input(QUE_IN
                         + "Create a Multiple Time-Series Plot? ('Q'=quit): ")
        else:
            code = input(QU2_IN
                         + "Create Another Multiple Series Plot? ('Q'=quit): ")
    else:
        # Single history plot.
        if query_stage == 0:
            code = input(QUE_IN
                         + "Create a Single Time-Series Plot? ('Q'=quit): ")
        else:
            code = input(QU2_IN
                         + "Create Another Single Time Plot? ('Q'=quit): ")

    # Get response.
    if code in POSITIVES:
        reply = True
    elif NEGATIVES or code in EXIT_COMMAND:
        reply = False
    else:
        reply = False   # Default is negative!

    return reply


def create_contour_query(trial):
    """Check if user wants to plot a time series.

    Parameters
    ----------
    trial = (int) question queue
            0 = first time
            1 = plot again (more plots?)

    Returns
    -------
    reply = (bool) answer plot question
    """
    # Ask question based on sequence.
    if trial == 0:
        code = input(QUE_IN + "Create a CO2 Contour Plot? ('Q'=quit): ")
    else:
        code = input(QU2_IN + "Create Another CO2 Contour Plot? ('Q'=quit): ")

    # Get response.
    if code in POSITIVES:
        reply = True
    elif NEGATIVES or code in EXIT_COMMAND:
        reply = False
    else:
        reply = False   # Default is negative!

    return reply


def range_query(max_realizations):
    """Get simulation numbers to plot.

    Parameters
    ----------
    max_realizations = (int) maximum number of realizations

    Returns
    -------
    response = (bool) code to plot or not
    start_run = (int) start realization to include
    end_run = (int) realization to include
    """
    # Default values.
    end_run = 0

    # Get start simulation for start of run with prompt.
    print(LED_IN + "START Realization:", end="")

    # Get start simulation number.
    starter = 0
    response, start_run = simulation_query(max_realizations, starter, "   ")

    # If OK, get end of range of simulations.
    if response:
        print(LED_IN + "END Realization:  ", end="")
        starter = start_run
        response, end_run = simulation_query(max_realizations, starter, "   ")

    return response, start_run, end_run


#
# -----------------------------------------------------------------------------
# - End of module
