#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions to compute second segment (tails) of profiles.

Author: Ernest N. Lindner
Date: 08/16/2022

Module Name
    flt_tail.py

Contents (6)
    trans_point()
    trans_point_temp()
    compute_poly_value(temp)
    compute_line_distance(fault_controls, start)
    compute_tail_length(fault_controls)
    compute_average_properties(fault_controls)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import math                     # For cosine

import flt_file as fileop       # For file operations
import flt_fluids as flu        # For T/P profile calculations
import flt_units as funit       # For unit conversion

# Constants for polynomial for saturation extension line.
POLY_A = 2.36E-05               # Saturation line extension = ax^3
POLY_B = 6.94E-04               # Saturation line extension = bx^2
POLY_C = 8.16E-02               # Saturation line extension = cx
POLY_D = 3.532                  # Saturation line extension = constant
BREAK_DEPTH = 854.7             # Depth for change in slope
BREAK_PRESSURE = 10.0E+06       # Pressure for change in slope


def trans_point():
    """Provide default transition point pressure and depth - Model #1.

    Parameters
    ----------
    N/A

    Returns
    -------
    pres = (float) Transition point pressure (Pa)

    Notes
    -----
    1. Assume two-line approximation profile - broken at press_break
    """
    press = BREAK_PRESSURE    # (Pa)
    depth = BREAK_DEPTH       # (m)
    return press, depth


def trans_point_temp():
    """Provide default transition point temperature - Model #1.

    Parameters
    ----------
    N/A

    Returns
    -------
    temp = (float) Transition point temperature (at 10 MPa)
    """
    temp = 42.25    # (oC)
    return temp


def compute_poly_value(temp):
    """Compute pressure along saturation line - complex model.

    Parameters
    ----------
    temp = (float) temperature of value (oC)

    Returns
    -------
    pressure = (float) current pressure (Pa)

    Notes
    -----
    1. Polynomial: P = ax^3 + bx^2 + cX + d, x=temperature

    """
    # Compute 3rd Order polynomial.
    term_1 = POLY_A * math.pow(temp, 3)
    term_2 = POLY_B * math.pow(temp, 2)
    term_3 = POLY_C * temp
    pressure = term_1 + term_2 + term_3 + POLY_D
    pressure *= funit.mpa_to_pa()

    return pressure


def compute_line_distance(fault_controls, start):
    """Compute distance of flow along dipping fault up to crux point.

    Parameters
    ----------
    fault_controls = (dict) with updated controls
    start = (float) depth of bottom of line (deeper) (i.e. start > x2)

    Returns
    -------
    length = (float) current length of flow profile tail (m)

    Notes
    -----
    1. Ensure that the result is positive.
    """
    # Define parameters; ensure positive value.
    term_1 = start  # deeper point
    term_2 = fault_controls['crux_depth']
    if term_2 > term_1:
        term_1 = term_2
        term_2 = start

    # Get fault inclination.
    dip = fault_controls.get('dip')
    dip = math.radians(dip)

    # Use trig. to compute distance on inclines fault.
    delta = term_1 - term_2
    length = math.sin(dip) * delta

    return length


def compute_tail_length(fault_controls):
    """Compute the travel profile length along fault for fault tail segment.

    Parameters
    ----------
    fault_controls = (dict) fault control values

    Returns
    -------
    extent = length (m)
    """
    # Establish default.
    base_depth = 0.0

    if fault_controls['profile_type'] == 1:
        # Model 1: Transition depth.
        base_depth = fault_controls['trans_point_depth']
    elif fault_controls['profile_type'] == 2:
        # Model 2: Supercritical depth is deepest point.
        base_depth = funit.supercrit_press_depth()
    else:
        msg = ("Code Error: " +
               "No Fault Model Identified in <compute_tail_length>")
        fileop.opx_problem(msg)

    # compute length from start depth to crux point.
    extent = compute_line_distance(fault_controls, base_depth)

    return extent


def compute_average_properties(fault_controls):
    """Compute average properties for tail at current time for models =1,2.

    Parameters
    ----------
    fault_controls = (dict) fault control values

    Returns
    -------
    co2_density = (float) density of CO2
    o2_viscosity = (float) viscosity of  CO2
    brine_density = (float) density of brine
    brine_viscosity = (float) viscosity of brine
    """
    # Default values.
    model_press = 0.0
    model_temp = 0.0

    # Establish start point for average properties, depending on model type.
    if fault_controls['profile_type'] == 1:
        model_press = fault_controls['trans_point_press']
        model_temp = trans_point_temp()
    elif fault_controls['profile_type'] == 2:
        model_temp = funit.supercrit_temperature()
        model_press = funit.supercrit_pressure()
    else:
        msg = ("Code Error: " +
               "No Fault Model Identified in <compute_tail_length>")
        fileop.opx_problem(msg)

    # Compute average conditions - pressure in MPa.
    ave_press = (model_press + fault_controls['crux_pressure']) / 2.0
    ave_press *= funit.pa_to_mpa()
    ave_temp = (model_temp + fault_controls['crux_temp']) / 2.0

    # Compute fluid properties for tail - temperature & pressure dependent.
    co2_density, co2_viscosity = \
        flu.co2_properties(ave_press, ave_temp)
    brine_density =  \
        flu.brine_density_property(ave_press, ave_temp)
    brine_viscosity = \
        flu.brine_viscosity_property(ave_press, ave_temp)

    return co2_density, co2_viscosity, brine_density, brine_viscosity

#
# -----------------------------------------------------------------------------
# End of module
