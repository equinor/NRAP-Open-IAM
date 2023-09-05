#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions for defining time and storage allocations.

Author: Ernest N. Lindner
Date: 08/23/2022

Module Name
    flt_storage

Contents (7)
    create_sim_storage(n_segments)
    create_reservoir_storage(n_segments)
    create_simulation_list(sim_step)
    create_time_values(fault_controls)
    clear_cycle_storage(co2_flow, brine_flow, one_dimension)
    clear_storage(co2_flow, brine_flow, simulation_list)
    update_simulation_list(sim_step, period, co2_results, brine_results,
                           sim_listing)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import numpy as np              # For array of values

import flt_file as fileop       # For file operations
import flt_config as scfg       # IO directory and file names
# import flt_units as funit     # For unit conversion


def create_sim_storage(n_segments):
    """Create the NumPy arrays storage arrays for calculation loops.

    Parameters
    ----------
    n_segments = (int) number of plates in a fault

    Returns
    -------
    sim_co2_flow = (array) CO2 flow for entire simulation
    sim_brine_flow = (array) brine flow for entire simulation
    co2_flow_step = (array) CO2 flow for time step
    brine_flow_step = (array) brine flow for time step
    """
    # Create storage arrays for a realization.
    sim_co2_flow = np.zeros(n_segments)
    sim_brine_flow = np.zeros(n_segments)

    # Create 1D storage arrays for each computation of delta time.
    co2_flow_step = np.zeros(n_segments)
    brine_flow_step = np.zeros(n_segments)

    return (sim_co2_flow, sim_brine_flow, co2_flow_step,
            brine_flow_step)


def create_reservoir_storage(n_segments):
    """Create the NumPy arrays for input from reservoir.

    Parameters
    ----------
    n_segments = (int) number of plates in a fault

    Returns
    -------
    base_co2_saturation = (array) base CO2 saturation from input
    base_co2_pressure = (array) base CO2 pressure from input
    """
    # Create reservoir arrays.
    base_co2_saturation = np.zeros(n_segments)
    base_co2_pressure = np.zeros(n_segments)

    return base_co2_saturation, base_co2_pressure


def create_simulation_list(sim_step, sim_listing):
    """Update the list accumulating results of a simulation.

    Parameters
    ----------
    sim_step = (int) simulation/realization step number

    Returns
    -------
    sim_listing = (list) with header information

    Notes
    -----
    Start simulation numbers at 1 (=> sim_step + 1).
    """
    # Define header for list - clear old list.
    sim_listing.clear()

    # Append first null (0) step as header.
    new_list = [(sim_step + 1), "0.00", "0.00", "0.00"]
    sim_listing.append(new_list)

    # Return new list.
    return sim_listing


def create_time_values(fault_controls):
    """Create a set of time steps for the analysis.

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault control values

    Returns
    -------
    time_steps = (array) NumPy array of float values (in yrs!)

    Notes
    -----
    Time steps include starting time of zero!
    """
    # Define major variables (in years)!
    numbr_steps = fault_controls['time_points']
    time_steps = np.zeros(numbr_steps)

    # Examine file input options.
    if fault_controls['time_input']:
        # Open the current directory and source file for input operations.
        sub_path, source_name = fileop.get_path_name(scfg.RESERVOIR_DIR,
                                                     scfg.TIME_STEP_NAME,
                                                     scfg.EXTENSION_TXT)
        try:
            # Read data as a NumPy array, excluding first header line.
            time_steps = np.genfromtxt(source_name, delimiter=",",
                                       autostrip=True, skip_header=1)
        except OSError as err:
            fileop.io_snag(err, sub_path, source_name)

        # Check for error in number of steps.
        if len(time_steps) != numbr_steps:
            msg1 = "Error: Failure in Time Series - Wrong Number of Points!"
            fileop.opx_problem(msg1)

    else:
        # Uniform time steps - define variables.
        start = fault_controls['start_time']
        end = fault_controls['end_time']

        # Create a series of uniform time steps.
        interval = (end - start) / (numbr_steps - 1)
        for i in range(0, numbr_steps):
            time_steps[i] = (i * interval)

    return time_steps


def clear_cycle_storage(co2_flow, brine_flow, one_dimension):
    """Reset arrays to data arrays for new cycle.

    Parameters
    ----------
    co2_flow = (array) CO2 flow
    brine_flow = (array) brine flow
    one_dimension = (int) array dimension

    Returns
    -------
    co2_flow => (array) zeroed
    brine_flow => (array) zeroed
    """
    brine_flow.shape = one_dimension
    co2_flow.shape = one_dimension
    co2_flow.fill(0.0)
    brine_flow.fill(0.0)

    return co2_flow, brine_flow


def clear_storage(co2_flow, brine_flow, simulation_list):
    """Reset storage arrays for simulation calculation.

    Parameters
    ----------
    co2_flow = (array) CO2 flow - 1D NumPy array
    brine_flow = (array) brine flow - 1D NumPy array
    simulation_list = (list) data list

    Returns
    -------
    co2_flow => (array) zeroed
    brine_flow => (array) zeroed
    simulation_list => (list) clear
    """
    co2_flow.fill(0.0)
    brine_flow.fill(0.0)
    simulation_list.clear()

    return co2_flow, brine_flow, simulation_list


def update_simulation_list(sim_step, period, co2_results, brine_results,
                           sim_listing):
    """Append data to the list of results from a simulation.

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

    Notes
    -----
    Start simulation numbers at 1 (=> sim_step + 1)
    """
    # Define data.
    track = str(sim_step + 1)
    clock = str(period)
    co2_value = str(co2_results)
    brine_value = str(brine_results)

    # Create new list and append.
    new_list = [track, clock, co2_value, brine_value]
    sim_listing.append(new_list)

    # Return in new list.
    return sim_listing


#
# -----------------------------------------------------------------------------
# End of module
