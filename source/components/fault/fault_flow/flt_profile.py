#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions for control of interpolation profiles - part #1.

Author: Ernest N. Lindner
Date: 08/17/2022

Module Name
    flt_profile

Contents (13)
    setup_transition_points(fault_controls)
    define_ave_interpolation_stages(alive, fault_controls, press_top,
    param_bounds)
    interpolate_initialize(fault_controls)
    compute_crux_point(timer)
    update_interpolate(timer, fault_controls)
    refresh_interpolate_controls(fault_controls)
    define_simple_ave(fault_controls)
    define_default_conditions(fault_controls)
    define_complex_ave(fault_controls)
    define_disjoint_ave(fault_controls)
    ---
    establish_crux_properties(pressure, temperature, fault_controls)
    ave_complex_properties(fault_controls)  NOT USED
    ave_disjoint_properties(fault_controls)  NOT USED

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import flt_file as fileop     # For file operations
import flt_fluids as flu      # For T/P profile calculations
import flt_message as mess    # For messages
import flt_tail as ftal       # Tail operations
import flt_units as funit     # For unit conversion

# Constants
SALINITY_DEFAULT = 0.0


def setup_transition_points(fault_controls):
    """Compute the profile type and conditions for near-surface transition.

    Parameters
    ----------
    fault_controls = (dict) fault control values

    Returns
    -------
    fault_controls = (dict) with updated control

    Notes
    -----
    1. Depths are positive values.
    2. Assume supercritical pressure is primary for finding transition.
    """
    # Condition if <aquifer> press./temp. are greater than supercritical.
    condition_below = (fault_controls['aquifer_pressure'] >=
                       funit.supercrit_pressure() and
                       fault_controls['aquifer_temperature']
                       >= funit.supercrit_temperature())

    # Check if near surface conditions are needed.
    if not fault_controls['near_surface_approach'] or condition_below:
        # Base Case Model #0.
        # --> "near-surface considerations" not selected,
        # --> OR profile is entirely supercritical.
        fault_controls['profile_type'] = 0  # default simple model
        fault_controls['simple_depth'] = fault_controls['aquifer_depth']
        fault_controls['simple_pressure'] = fault_controls['aquifer_pressure']
        fault_controls['simple_temperature'] = \
            fault_controls['aquifer_temperature']

    elif fault_controls['profile_type'] == 1:
        # Complex Model #1.
        # Use 10 MPa depth as break point.
        press, upper_depth = ftal.trans_point()

        # If fault starts at shallower depth, use injection point as start.
        if fault_controls['inject_depth'] >= upper_depth:
            # Use default values.
            fault_controls['trans_point_press'] = press
            fault_controls['trans_point_depth'] = upper_depth
            fault_controls['trans_point_temp'] = \
                funit.geothermal_temperature(upper_depth)
        else:
            # Use injection depth as reference as injection is too shallow.
            upper_depth = fault_controls["inject_depth"]
            fault_controls['trans_point_depth'] = upper_depth
            fault_controls['trans_point_press'] = \
                fault_controls['inject_pressure']
            fault_controls['trans_point_temp'] = \
                funit.geothermal_temperature(upper_depth)

        # Crux point starts at supercritical depth & function of time.
        fault_controls['crux_depth'] = funit.supercrit_press_depth()
        fault_controls['crux_pressure'] = funit.supercrit_pressure()
        fault_controls['crux_temp'] = funit.supercrit_temperature()

    elif fault_controls['profile_type'] == 2:
        # Disjoint Model #2.
        fault_controls['disjoint_depth'] = \
            funit.depth_for_pressure(funit.supercrit_pressure())
        fault_controls['disjoint_pressure'] = funit.supercrit_pressure()
        fault_controls['disjoint_temperature'] = \
            funit.geothermal_temperature(fault_controls['disjoint_depth'])

        fault_controls['end_depth'] = fault_controls['aquifer_depth']
        fault_controls['end_pressure'] = fault_controls['aquifer_pressure']
        fault_controls['end_temperature'] = \
            fault_controls['aquifer_temperature']

    else:
        msg = ("Code Error: " +
               "No Fault Model Defined in <setup_transition_points>")
        fileop.opx_problem(msg)

    return fault_controls


def define_ave_interpolation_stages(alive, fault_controls):
    """Develop average properties for brine & CO2 for segment #1.

    Parameters
    ----------
    alive = (bool) flag on stand-alone operations
    fault_controls = (dict) fault controls

    Returns
    -------
    fault_controls = (dict) fault controls (updated)
    """
    # Echo status to user.
    mess.echo_status(alive, "SURFACE")

    # Select conditions based on profile type.
    if fault_controls['profile_type'] == 0:
        # Simple model parameters. (profile=0).
        fault_controls = define_simple_ave(fault_controls)

    elif fault_controls['profile_type'] == 1:
        # Complex model parameters.(profile=1).
        fault_controls = define_complex_ave(fault_controls)

    elif fault_controls['profile_type'] == 2:
        # Disjoint model parameters. (profile=2).
        fault_controls = define_disjoint_ave(fault_controls)
    else:
        # Error.
        msg = 'Error: Model Number Wrong in >define_ave_interpolation_stages<'
        fileop.opx_problem(msg)

    return fault_controls


def interpolate_initialize(fault_controls):
    """Define fluid property interpolation for injection period.

    Parameters
    ----------
    fault_controls = (dict) fault controls

    Returns
    -------
    fault_controls = (updated) dictionary of fault controls
    """
    # Set default values for interpolation tracker and values during injection.
    if fault_controls['interpolate_approach']:
        fault_controls = \
            flu.reset_fluid_properties(fault_controls,
                                       fault_controls['injection_properties'])

        # Turn-on interpolate switch.
        fault_controls['update'] = True

    else:
        # Otherwise, turn-off interpolate switch.
        fault_controls['update'] = False

    return fault_controls


def compute_crux_point(timer, fault_controls):
    """Determine temperature and pressure at crux point.

    Parameters
    ----------
    timer = (float) time of simulation after injection start
    fault_controls = (dict) fault control values

    Returns
    -------
    temp = (float) temperature of crux point
    press = (float) pressure of crux point

    Notes
    -----
    1. Crux point is on the extended Saturation Line.
    """
    # Values per crux equation.
    if timer < 0.1:
        # If minor time, no change.
        press = 12.0
        temp = 49.0
    elif timer <= 49.0:
        # Curve fit.
        press = 8.0 * pow(timer, -0.175)
        temp = 24.0 * pow(timer, -0.350)
    else:
        # Values at larger times.
        press = 4.04
        temp = 6.10

    # Crux point pressure must be greater than or equal to aquifer.
    #  - early termination limit.
    if press < fault_controls['aquifer_pressure']:
        press = fault_controls['aquifer_pressure']
        temp = fault_controls['aquifer_temperature']

    return press, temp


def update_interpolate(timer, fault_controls):
    """Determine if fluid property update occurs after start loop.

    Parameters
    ----------
    timer = (float) current time in years
    fault_controls = (dict) fault controls

    Returns
    -------
    fault_controls = (dict) fault controls (updated)

    Notes
    -----
    1. Update fluid parameters after injection period.
    2. "update" prevents changes after the one at injection period end.
    """
    # Update crux values for complex profiles (type=1, =2).
    if fault_controls['profile_type'] > 0:
        crux_press, crux_temp = \
            compute_crux_point(timer, fault_controls)
        fault_controls = \
            establish_crux_properties(crux_press, crux_temp,
                                      fault_controls)

    # Interpolation changes after injection completed but only once.
    if fault_controls['update']:
        if timer > fault_controls['inject_end']:
            # Update fluid parameters.
            fault_controls = \
                flu.reset_fluid_properties(fault_controls,
                                           fault_controls[
                                               'long_term_properties'])
            # Turn-off interpolate control to stop further updates.
            fault_controls['update'] = False

    return fault_controls


def define_simple_ave(fault_controls):
    """Develop average properties for "Simple" profile segment.

    Parameters
    ----------
    fault_controls = (dict) fault controls

    Returns
    -------
    fault_controls = (dict) fault controls (updated)

    Notes
    -----
    1. Creates two lists for different time periods for single segment:
        (1) During injection, "injection_properties".
        (2) After injection, "long_term_properties".
    """
    # -------------------------------------------------------------------------
    # Simple model conditions  - segment #1 - Injection.
    if fault_controls['interpolate_approach']:

        # Compute properties - at Injection base for fluid flow.
        temp = fault_controls['inject_temperature']
        press = fault_controls['inject_pressure']
        saline = fault_controls['salinity']
        injection_list = flu.interpolate_fluid_properties(temp, press, saline)

        # Compute properties - At base of aquifer for fluid flow.
        temp = fault_controls['aquifer_temperature']
        press = fault_controls['aquifer_pressure']
        aquifer_list = flu.interpolate_fluid_properties(temp, press, saline)
        fault_controls['aquifer_properties'] = aquifer_list

        # Compute properties - At reservoir level at All other times.
        temp = fault_controls['inject_temperature']
        press = fault_controls['final_pressure']
        duration_list = flu.interpolate_fluid_properties(temp, press, saline)

        # Save/reset aquifer properties for use later.
        fault_controls = flu.reset_aqui_properties(fault_controls,
                                                   aquifer_list)

    else:
        # No interpolation - use default values from input.
        injection_list, aquifer_list, duration_list = \
            define_default_conditions(fault_controls)

    # -------------------------------------------------------------------------
    # --> Define average parameters during injection.
    fault_controls['injection_properties'] = \
        flu.ave_properties(injection_list, aquifer_list)

    # --> Define average parameters for long-term conditions.
    fault_controls['long_term_properties'] = \
        flu.ave_properties(duration_list, aquifer_list)

    return fault_controls


def define_default_conditions(fault_controls):
    """Develop default properties for profiles..

    Parameters
    ----------
    fault_controls = (dict) fault controls

    Returns
    -------
    injection_list = (list) injection properties
    aquifer_list = (list) aquifer properties
    duration_list = (list) long-term properties

    Notes
    -----
    1. Creates two lists for different time periods for single segment:
        (1) During injection, "injection_properties".
        (2) After injection, "long_term_properties".
    """
    # Define default values.
    injection_list = [fault_controls.get('co2_density'),
                      fault_controls.get('co2_viscosity'),
                      fault_controls.get('brine_density'),
                      fault_controls.get('brine_viscosity'),
                      fault_controls.get('co2_solubility')]

    aquifer_list = [fault_controls.get('aqui_co2_density'),
                    fault_controls.get('aqui_co2_viscosity'),
                    fault_controls.get('aqui_brine_density'),
                    fault_controls.get('aqui_brine_viscosity'),
                    fault_controls.get('co2_solubility')]

    duration_list = [fault_controls.get('co2_density'),
                     fault_controls.get('co2_viscosity'),
                     fault_controls.get('brine_density'),
                     fault_controls.get('brine_viscosity'),
                     fault_controls.get('co2_solubility')]

    return injection_list, aquifer_list, duration_list


def define_complex_ave(fault_controls):
    """Develop average properties for complex segment #1 (trunk).

    Parameters
    ----------
    fault_controls = (dict) fault controls

    Returns
    -------
    fault_controls = (dict) fault controls (updated)

    Notes
    -----
    1. Creates lists for segments:
        (1) During injection, "injection_properties".
        (2) After injection, "long_term_properties".
    """
    # -------------------------------------------------------------------------
    # Model 1 conditions  - segment #1 - Injection.
    # --> Compute properties - Injection base for fluid flow.
    if fault_controls['interpolate_approach']:

        # Compute properties - at Injection base for fluid flow.
        injection_list = flu.interpolate_fluid_properties(
            fault_controls['inject_temperature'],
            fault_controls['inject_pressure'],
            fault_controls['salinity'])

        divide_list = flu.interpolate_fluid_properties(
            fault_controls['trans_point_temp'],
            fault_controls['trans_point_press'],
            fault_controls['salinity'])

        duration_list = flu.interpolate_fluid_properties(
            fault_controls['inject_temperature'],
            fault_controls['final_pressure'],
            fault_controls['salinity'])

    else:
        # No interpolation - use default values from input.
        injection_list, aquifer_list, duration_list = \
            define_default_conditions(fault_controls)

        # Compute divide point list as average as compromise.
        divide_list = flu.ave_properties(aquifer_list, duration_list)

    # -------------------------------------------------------------------------
    # --> Define average parameters during injection for trunk.
    fault_controls['injection_properties'] = \
        flu.ave_properties(injection_list, divide_list)

    # --> Define average parameters for long-term conditions for trunk.
    fault_controls['long_term_properties'] = \
        flu.ave_properties(duration_list, divide_list)

    return fault_controls


def define_disjoint_ave(fault_controls):
    """Develop average properties for Disjoint Profile segment #1 (trunk).

    Parameters
    ----------
    fault_controls = (dict) fault controls

    Returns
    -------
    fault_controls = (dict) fault controls (updated)

    Notes
    -----
    1. Creates two lists for different time periods for single segment:
        (1) During injection, "injection_properties".
        (2) After injection, "long_term_properties".
    """
    # -------------------------------------------------------------------------
    if fault_controls['interpolate_approach']:
        # Disjoint model conditions  - segment #1 - Injection.
        # Compute properties - Injection base for fluid flow.
        injection_list = flu.interpolate_fluid_properties(
            fault_controls['inject_temperature'],
            fault_controls['inject_pressure'],
            fault_controls['salinity'])

        # Compute properties - At supercritical pressure boundary.
        depth = funit.supercrit_press_depth()
        temperature = funit.geothermal_temperature(depth)
        bound_list = flu.interpolate_fluid_properties(
            temperature,
            funit.supercrit_pressure(),
            fault_controls['salinity'])

        # Compute properties - At reservoir level at All other times.
        duration_list = flu.interpolate_fluid_properties(
            fault_controls['inject_temperature'],
            fault_controls['final_pressure'],
            fault_controls['salinity'])

    else:
        # No interpolation - use default values from input.
        injection_list, aquifer_list, duration_list = \
            define_default_conditions(fault_controls)

        # Compute divide point list as average as compromise.
        bound_list = flu.ave_properties(aquifer_list, duration_list)

    # -------------------------------------------------------------------------
    # Disjoint model conditions  - segment #1 - Long_term
    # --> Define average parameters during injection.
    fault_controls['injection_properties'] = \
        flu.ave_properties(injection_list, bound_list)

    # --> Define average parameters for long-term conditions.
    fault_controls['long_term_properties'] = \
        flu.ave_properties(duration_list, bound_list)

    return fault_controls


def establish_crux_properties(pressure, temperature, fault_controls):
    """Establish properties of the crux point.

    Parameters
    ----------
    pressure = (float) crux point pressure
    temperature = (float) crux point temperature
    fault_controls = (dict) fault control values

    Returns
    -------
    fault_controls = (dict) with updated controls

    Notes
    -----
    1. Aquifer depth is not deep enough for supercritical CO2 zone.
    2. crux_depth = depth where CO2 liquid phase ends.
    """
    # Define conditions at endpoint on Saturation Line.
    fault_controls['crux_pressure'] = pressure
    fault_controls['crux_temp'] = temperature
    fault_controls['crux_depth'] = funit.depth_for_pressure(pressure)

    return fault_controls


def refresh_interpolate_controls(fault_controls):
    """Refresh the interpolation parameters for a new sim.

    Parameters
    ----------
    fault_controls = (dict) fault controls

    Returns
    -------
    fault = (dict) fault controls (updated)
    """
    # ------------------------------------------
    # Reset default values for interpolation tracker and values (Pa).
    fault_controls['current_pressure'] = fault_controls['inject_pressure']

    return fault_controls


def ave_complex_properties(fault_controls):
    """Develop ave. properties for complex model, segment #2 - w. crux point.

    Parameters
    ----------
    fault_controls = (dict) fault controls

    Returns
    -------
    fault_controls = (dict) fault controls (updated)
    """
    # -------------------------------------------------------------------------
    # Complex model conditions  - segment #2.

    # Use interpolation values.
    if fault_controls['interpolate_approach']:
        divide_list = flu.interpolate_fluid_properties(
            fault_controls['trans_point_temp'],
            fault_controls['trans_point_press'],
            fault_controls['salinity'])
        crux_list = flu.interpolate_fluid_properties(
            fault_controls['crux_temp'],
            fault_controls['crux_pressure'],
            SALINITY_DEFAULT)     # Low salinity aquifer

    # Otherwise, use defined values as a compromise.
    else:
        crux_list = fault_controls['aquifer_properties']
        divide_list = flu.ave_properties(
            fault_controls['aquifer_properties'],
            fault_controls['long_term_properties'])

    # --> Define average parameters during injection.
    fault_controls['tail_properties'] = \
        flu.ave_properties(divide_list, crux_list)

    return fault_controls


def ave_disjoint_properties(fault_controls):
    """Develop average properties for disjoint segment #2 - w. crux point.

    Parameters
    ----------
    fault_controls = (dict) fault controls

    Returns
    -------
    fault_controls = (dict) fault controls (updated)
    """
    # -------------------------------------------------------------------------
    # Disjoint model conditions  - segment #2.

    # Use interpolation values.
    if fault_controls['interpolate_approach']:
        supercrit_list = \
            flu.interpolate_fluid_properties(
                                             funit.supercrit_temperature(),
                                             funit.supercrit_pressure(),
                                             fault_controls['salinity'])
        crux_list = \
            flu.interpolate_fluid_properties(
                                             fault_controls['crux_temp'],
                                             fault_controls['crux_pressure'],
                                             SALINITY_DEFAULT)

    # Otherwise use defined values as a compromise.
    else:
        crux_list = fault_controls['aquifer_properties']
        supercrit_list = flu.ave_properties(
            fault_controls['aquifer_properties'],
            fault_controls['long_term_properties'])

    # --> Define average parameters during injection.
    fault_controls['tail_properties'] = \
        flu.ave_properties(supercrit_list, crux_list)

    return fault_controls

#
# -----------------------------------------------------------------------------
# End of module
