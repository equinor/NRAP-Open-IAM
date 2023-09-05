#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions to compute seal thickness, layout, etc.

Author: Ernest N. Lindner
Date: 08/24/2022

Module Name
    seal_refresh

Contents (16)
    aggregate_results(co2_flow, brine_flow)
    thickness_variability(seal_controls, state)
    filter_thickness_data(big_array, min_val)
    create_time_values(seal_controls)
    establish_simulation_list(sim_step)
    update_simulation_list(sim_step, period, co2_results, brine_results,
                           sim_listing)
    refresh_data_arrays(seal_controls, sim_step, grid, press_top)
    create_sim_storage(seal_controls)
    create_reservoir_storage(seal_controls)
    clear_cycle_storage(co2_flow, brine_flow)
    ----
    clear_storage(co2_flow, brine_flow, simulation_list)
    print_data_arrays(grid_table, lines, columns, print_list, outlook)
    debug_check(checkup, grid, seal_controls)
    closeout(alone, start_code)
    echo_time_step(desired, alive, current)
    echo_sim_step(alive, sim_step)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import sys                          # For standard output
import time                         # For calculating time
from datetime import timedelta      # For calculating runtime
from scipy.stats import truncnorm   # For stochastic thickness generation
from scipy.ndimage import gaussian_filter    # For smoothing thickness data
import numpy as np                  # For arrays & input

import seal_config as scfg           # IO directory and filenames
import seal_file as sfile           # For file operations
import seal_units as sun            # For unit conversion

# Other operations constants
MUSE = 1.0                          # Level of Gaussian smoothing
OTH_IN = "   >> "                   # Info line start
LIN_IN = "  --> "                   # Line start

# Debug constants
# Note: For large grids, set DEBUG_TO_FILE to "True" as console cannot process
#       larger files and print to console shows only portion of arrays.
DEBUG_TO_FILE = False                    # Control to print Debug to file
OPTIONS = ["PERMEABILITY", "THICKNESS", "INFLUENCE", "AREA", "DEPTH",
           "STATUS", "X-COORDINATES", "Y-COORDINATES", "ENTRY PRESSURE",
           "HISTORY"]     # Debugging printout headers


def aggregate_results(co2_flow, brine_flow):
    """Compute cumulative flows for entire grid for a time step.

    Parameters
    ----------
    brine_flow_step = (array) brine flow during step
    co2_flow_step = (array) CO2 flow during step

    Returns
    -------
    sum_co2_flow = (float) total delta CO2 for grid
    sum_brine_flow = (float) total delta brine flow for grid
    """
    # Compute total CO2 flow and total brine flow.
    sum_co2_flow = co2_flow.sum()
    sum_brine_flow = brine_flow.sum()

    return sum_co2_flow, sum_brine_flow


def thickness_variability(seal_controls, state):
    """Compute an array of thickness using a truncated normal distribution.

    Parameters
    ----------
    seal_controls = (dict) seal parameters
    state = (int) random seed

    Returns
    -------
    out_array = (array) thickness (1D)
    """
    # Define variables and computation array.
    mean_val = seal_controls['thickness_ave']
    sigma = seal_controls['thickness_std']
    min_val = seal_controls['thickness_min']
    max_val = seal_controls['thickness_max']
    numbr = seal_controls['num_cells']
    out_array = np.zeros(numbr)

    # If variability is set small, do not vary thickness!
    if sigma > 0.01:
        # Define limits for truncnorm.
        lower_bound = (min_val - mean_val) / sigma
        upper_bound = (max_val - mean_val) / sigma

        # compute array from random distribution.
        big_array = truncnorm.rvs(lower_bound, upper_bound,
                                  loc=mean_val, scale=sigma,
                                  size=numbr, random_state=state)

        # Smooth data if in rectangular array.
        if seal_controls['grid_approach']:
            # Gaussian smoothing
            lines = seal_controls['grid_rows']
            columns = seal_controls['grid_cols']
            interim_array = big_array.reshape(lines, columns)

            out_array = filter_thickness_data(interim_array, min_val)
        else:
            # No smoothing of data.
            out_array = big_array
    else:
        # For low variability, set all to mean & no smoothing.
        out_array.fill(mean_val)

    # Return flat (1D) array.
    return out_array


def filter_thickness_data(rough_array, minimum):
    """Smooth the random independent values of thickness.

    Parameters
    ----------
    rough_array = (array of float) independent random thickness values
    minimum = (float) minimum thickness

    Returns
    -------
    smooth_array = (array of float) smoothed array 1D

    Notes
    -----
    No check on maximum as smoothing reduces peak values.
    Applicable to rectangular grid only.
    """
    # Define output array and standard deviation of current array.
    old_dev = rough_array.std()

    # Define gauss filter parameters - isotropic assumption.
    sigma = [MUSE, MUSE]

    # Apply filter – “reflect” values at boundary edge.
    filtered = gaussian_filter(rough_array, sigma,
                               mode='reflect')

    # Determine stochastic variables of new array.
    # Issue with array -> explicitly set to Numpy array.
    new_array = np.asarray(filtered)
    smooth_array = new_array.flatten()
    center = smooth_array.mean()
    new_dev = smooth_array.std()

    # Restore standard deviation of sample to create smooth array.
    #   -> Do not change for low deviation values.
    if new_dev > 0.001:
        npts = smooth_array.size
        ratio_dev = old_dev / new_dev

        for indx in range(npts):
            current = smooth_array[indx]
            smooth_array[indx] = ((current - center) * ratio_dev
                                  + center)

    # Ensure that all values are greater than the minimum.
    smooth_array[smooth_array < minimum] = minimum

    return smooth_array


def create_time_values(seal_controls):
    """Create a set of time steps for the analysis.

    Parameters
    ----------
    seal_controls = (dict) seal parameters

    Returns
    -------
    time_steps = (array of float) time (in yrs!)

    Notes
    -----
    Time steps include starting time of zero!
    """
    # Define major variables (in years)!
    numbr_steps = seal_controls['time_points']
    time_steps = np.zeros(seal_controls['time_points'])

    # Examine file input options.
    if seal_controls['time_input']:
        # Open the current directory and source file for input operations.
        sub_path, source_name = sfile.get_path_name(scfg.RESERVOIR_DIR,
                                                    scfg.TIME_STEP_NAME,
                                                    scfg.EXTENSION_TXT)
        try:
            # Read data as a NumPy array, excluding first header line.
            time_steps = np.genfromtxt(source_name, delimiter=",",
                                       autostrip=True, skip_header=1)
        except OSError as err:
            sfile.io_snag(err, sub_path, source_name)

        # Check for error in number of steps.
        if len(time_steps) != numbr_steps:
            msg1 = "Failure in Time Series - Wrong Number of Points!"
            sfile.opx_problem(msg1)

        # Correct data if in days.
        # time_steps *= sun.days_to_yrs()

    else:
        # Define variables.
        start = seal_controls['start_time']
        end = seal_controls['end_time']

        # Create a series of uniform time steps.
        interval = (end - start) / (numbr_steps - 1)
        for i in range(0, numbr_steps):
            time_steps[i] = (i * interval)

    return time_steps


def establish_simulation_list(sim_step):
    """Create a list for accumulating results of a simulation.

    Parameters
    ----------
    sim_step = (int) simulation/realization step number

    Returns
    -------
    sim_listing = (list) data list with header information
    """
    # Define header for list.
    sim_listing = []

    # Append first null step - start number at 1.
    new_list = [(sim_step + 1), "0.00", "0.00", "0.00"]
    sim_listing.append(new_list)

    # Return new list.
    return sim_listing


def update_simulation_list(sim_step, period, co2_results, brine_results,
                           sim_listing):
    """Append data to list of results from a simulation.

    Parameters
    ----------
    sim_step = (int) simulation/realization number
    period = (float) time value of step
    co2_results = (float) current CO2 cumulative flow
    brine_results = (float) current brine cumulative flow
    sim_listing = (list) accumulation list of results for simulation

    Returns
    -------
    sim_listing = (list) accumulation list of results for simulation
    """
    # Define data; steps start at 1 for output.
    track = str(sim_step + 1)
    clock = str(period)
    co2_value = str(co2_results)
    brine_value = str(brine_results)

    # Create new list and append.
    new_list = [track, clock, co2_value, brine_value]
    sim_listing.append(new_list)

    # Return new list.
    return sim_listing


def refresh_data_arrays(seal_controls, sim_step, grid, press_top, rng):
    """Refresh thickness and pressure arrays as needed.

    Parameters
    ----------
    seal_controls = (dict) seal parameters
    sim_step = (int) simulation/realization number
    grid = (list) (class) a collection of cells
    press_top = (float) pressure at top of cells
    rng = (generator) random number generator

    Returns
    -------
    press_top = (array) (original or updated)

    Notes
    -----
    1. Do not update on last simulation step - set Boolean.
    2. Do not adjust thickness, if it is from file.
    3. Pressure not changed, if thickness not changed.
    """
    # Set criterion for last step.
    step_allow = (sim_step < (seal_controls['realizations'] - 1))

    # Re-evaluate thickness if not from file and not on last step.
    if not seal_controls['thickness_approach'] and step_allow:

        # Get new stochastic thickness.
        new_thickness_array = thickness_variability(seal_controls, rng)

        # Re-evaluate pressure.
        if not seal_controls['upper_approach']:
            #  Update thickness and adjust pressure - if not from file.
            for numbr, cell in enumerate(grid):

                # Update cell thickness - save old value, set new.
                old_thickness = cell.thickness
                cell.thickness = new_thickness_array[numbr]

                # Compute changes relative to original cell data.
                depth_change = old_thickness - cell.thickness
                pressure_change = (cell.brineDensity * depth_change
                                   * sun.gravity())

                # Update both pressure and current depth of cell.
                press_top[numbr] += pressure_change
                cell.top_seal_depth += depth_change

        else:
            # Pressure from file - update thickness values only.
            for numbr, cell in enumerate(grid):
                cell.thickness = new_thickness_array[numbr]

    # else: No thickness or pressure change

    return press_top


def create_sim_storage(seal_controls):
    """Create NumPy arrays storage arrays for loops.

    Parameters
    ----------
    seal_controls = (dict) seal parameters

    Returns
    -------
    sim_co2_flow = (array) CO2 flow for simulation
    sim_brine_flow = (array) brine flow for simulation
    co2_flow_step = (array) CO2 flow for time step
    brine_flow_step = (array) brine flow for time step

    """
    # Define array structure.
    total = seal_controls['num_cells']

    # Create 1D storage arrays for a realization.
    sim_co2_flow = np.zeros(total)
    sim_brine_flow = np.zeros(total)

    # Create 1D storage arrays for each computation of delta time.
    co2_flow_step = np.zeros(total)
    brine_flow_step = np.zeros(total)

    return (sim_co2_flow, sim_brine_flow, co2_flow_step,
            brine_flow_step)


def create_reservoir_storage(seal_controls):
    """Create NumPy arrays for input from reservoir.

    Parameters
    ----------
    seal_controls = (dict) seal parameters

    Returns
    -------
    base_co2_saturation = (array) base CO2 saturation from input
    base_co2_pressure = (array) base CO2 pressure from input
    """
    # Define array structure.
    total = seal_controls['num_cells']

    # Create reservoir arrays.
    base_co2_saturation = np.zeros(total)
    base_co2_pressure = np.zeros(total)

    return base_co2_saturation, base_co2_pressure


def clear_cycle_storage(co2_flow, brine_flow):
    """Convert arrays to 1D and clear arrays for new loop.

    Parameters
    ----------
    co2_flow = (array) CO2 flow - NumPy array
    brine_flow = (array) brine flow - NumPy array
    total_cells = (float) total number of cells

    Returns
    -------
    co2_flow = (array) zeroed
    brine_flow = (array) zeroed
    """
    co2_flow.fill(0.0)
    brine_flow.fill(0.0)

    return co2_flow, brine_flow


def clear_storage(co2_flow, brine_flow, simulation_list):
    """Clear storage arrays for loop.

    Parameters
    ----------
    co2_flow = (array) CO2 flow - NumPy array
    brine_flow = (array) brine flow - NumPy array
    simulation_list = (list) sim data

    Returns
    -------
    co2_flow = (array) zeroed
    brine_flow = (array) zeroed
    simulation_list = (list) cleared
    """
    # Clear arrays for a time loop.
    co2_flow.fill(0.0)
    brine_flow.fill(0.0)
    simulation_list.clear()

    return co2_flow, brine_flow, simulation_list


def print_data_arrays(grid_table, lines, columns, print_list, outlook):
    """Print data values in controlled format / table.

    Parameters
    ----------
    grid_table = (list) (class) list of cells
    lines = (int) number of lines in return array
    columns = (int) number of columns in return array
    print_list  = (list on ints) list of selections to print:
        = 0: permeability
        = 1: thickness
        = 2: influence
        = 3: area
        = 4: depth
        = 5: status
        = 6: x-coordinate
        = 7: y-coordinate
        = 8: entry pressure
        = 9: history
    outlook = (str) file destination

    Returns
    -------
    None

    Notes
    -----
    Change "n_vals" for fewer columns in a line.
    Permeability set later than other input!
    """
    # If list exists, print each item in list.
    if print_list:
        for select in print_list:
            # Print header and setup for printout.
            print("\n  " + OPTIONS[select], file=outlook)

            # Define columns for printing.
            n_vals = columns - 1

            #  Print table.
            for i in range(lines):
                print("   ", end="", file=outlook)
                for j in range(columns):
                    stat_cell = j + (i * columns)

                    # Select print results as defined.
                    valu = grid_table[stat_cell].get_value(select)

                    # print results with correct format.
                    if select == 5:    # Status values are integers
                        print(f' {valu:2d}', end="", file=outlook)
                    else:
                        print(f' {valu:7.3e}', end="", file=outlook)

                    # check for new line.
                    if j % n_vals == 0 and j > 0:
                        print("", file=outlook)
                        print("", end="", file=outlook)

    # return None


def debug_check(checkup, grid, seal_controls):
    """Print arrays for debugging/checking.

    Parameters
    ----------
    checkup = (bool) flag to print list for debugging (set in seal_flux)
    grid = (list of class) = seal cells
    seal_controls = (dict) seal component values

    Returns
    -------
    None

    Notes
    -----
    1. Printout per code in print_list.
    2. Permeability is not input before computation loop start!
    3. See seal_model function: "get_value" for correlation eq. in code.
    4. Code: perm.=0; thickness=1; influence=2; area=3; depth=4;
         status=5; x-coordinates=6; y-coordinates=7;
         entry pressure=8; history=9.
    """
    # Print cell data, if desired.
    if checkup:
        # Define lists/parameters to be printed.
        print_list = [1, 2, 3, 4, 5, 6, 7, 8, 9]  # See "seal_model".

        # Define variables.
        lines = seal_controls['r_lines']
        columns = seal_controls['r_columns']

        # Print to file, if desired.
        if DEBUG_TO_FILE:
            # For file output, get path.
            sub_path, destination = sfile.get_path_name(scfg.OUTPUT_DIRECTORY,
                                                        scfg.DEBUG_FILE,
                                                        scfg.EXTENSION_TXT)

            # Write information to file, "outlook".
            try:
                with open(destination, 'w', encoding="utf8") as outlook:
                    # Print selected arrays.
                    print_data_arrays(grid, lines, columns, print_list,
                                      outlook)
            except OSError as err:
                sfile.io_snag(err, sub_path, scfg.DEBUG_FILE)

        else:
            # Print to console.
            outlook = sys.stdout
            print_data_arrays(grid, lines, columns, print_list, outlook)

        # Clear list for other uses.
        print_list.clear()

    # Print properties of a cell.
    # marker = 12
    # grid[marker].print_cell()  # ==> Commented out

    # return None


def closeout(alone, start_code):
    """Compute elapsed time and write to console.

    Parameters
    ----------
    alone = (bool) flag to indicate stand-alone operation
    start_code = (float) start time

    Returns
    -------
    elapse = (float) run time for seal code
    """
    # Compute run time and echo status to user.
    # Stop clock on runtime.
    end_code = time.monotonic()

    # Echo end, compute time analysis.
    sfile.echo_status(alone, "COMPLETED ANALYSIS COMPUTATIONS.")
    elapse = timedelta(seconds=(end_code - start_code))

    # Echo execution time.
    if alone:
        print(OTH_IN + f'Analysis Execution Time = {elapse}',
              flush=True)

    return elapse


def echo_time_step(desired, alive, current):
    """Echo time step progress to console.

    Parameters
    ----------
    desired = (bool) flag to indicate header is to be printed
    alive = (bool) stand-alone operation
    current = (int) current time step

    Returns
    -------
    None
    """
    # Echo progress line to user on console.
    if desired and alive:
        print(LIN_IN + f'Time Step = {current}')

    # return None


def echo_sim_step(alive, sim_step):
    """Echo simulation progress to console.

    Parameters
    ----------
    alive = (bool) stand-alone operation
    sim_step = (int) current simulation step

    Returns
    -------
    None
    """
    # Echo progress line to user on console.
    if alive:
        print(OTH_IN + f'Realization Loop = {(sim_step + 1)}')

    # return None


#
# -----------------------------------------------------------------------------
# - End of module
