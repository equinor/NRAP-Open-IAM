#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions for misc. fault setup and computations.

Author: Ernest N. Lindner
Date: 08/23/2022

Module Name
    flt_compute

Contents (19)
    random_seeder(switch=0)
    check_if_fault_exists(fault_controls, rng)
    compute_sgr_influence(fault_controls)
    setup_apertus(fault_controls)
    compute_new_gap(fault_controls, stress, gamma_ref, beta, vary)
    define_apertus(fault_controls)
    define_fracto(fault_controls)
    stress_factor(fault_controls, old_strike, new_strike)
    project_to_subsurface(fault_controls)
    create_features(fault_controls, fault, x_valu, y_valu)
    ---
    define_aspects(fault_controls, fault)
    define_amplius(fault_controls, fault, param_bounds, rng)
    compute_fault_length(fault_controls)
    fault_plane_eq(fault_controls)
    test_oblique(planar, focus, fault_controls)
    distance_vector(planar)
    is_almost_equal(x, y, epsilon=1.0e-08)
    check_inclination(fault_controls)
    fault_create(fault_controls, fault, param_bounds, rng)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import sys                      # For stdout
import math                     # For sin and other functions
import numpy as np              # For array of values

import flt_inout as fio         # File input/output operations
import flt_file as fileop       # For file operations
import flt_model as model       # For aperture variation
import flt_units as funit       # For unit conversion
import flt_basics as bcs        # For basic math operations
import flt_config as scfg       # For a fixed random seed

# SGR computation constants
SGR_SLOPE = 0.02                # Slope of linear pressure relation of SGR
SGR_PRESS = 0.4
SGR_MAXIMUM = 100.0             # Max. value of SGR for test
SGR_MINIMUM = 20.0              # Min. value of SGR for effect
SGR_CONSTANT = -4.0             # Constant for permeability relation

# Pressure-aperture (pa) constants
STIFFNESS = 1.8E+6              # (Pa) Fracture stiffness
LIMIT_STRESS = 20.0E+6          # (Pa) Maximum stress for relation / curve
MAX_APERTURE = 2                # (mm) Maximum aperture
RESIDUAL_APERTURE = 0.01        # (mm) Residual aperture

# Default constants
V_SMALL_ANGLE = 11.0            # Dip should be > 10.0
V_SMALL_DIFF = 1.0              # Very small difference in depth
V_SMALL_ZERO = 1.0E-08          # Very small number in division by zero
V_SMALL_ONE = 0.01
FAULT_LIMIT = 99.9              # Existence of fault probability limit
DIP_MIN = 10.0                  # minimum fault dip
DIP_MAX = 90.0                  # maximum dip

# Other constants
ECHO = False                    # Debug - Printout intermediate results
ECHO_2 = False                  # Debug - Printout plate results
EXTRA_D = 200                   # Extra distance on inclusion algorithm


def random_seeder(switch=0):
    """Define random function seed.

    Parameters
    ----------
    switch = (int) flag to increase SEEDX

    Returns
    -------
    rng = random number generator function
    """
    # Set seed for random numbers.
    if scfg.SEEDX is None:
        # Set undefined seed.
        rng = np.random.default_rng()
    else:
        # Set defined random seed.
        randy = scfg.SEEDX
        rng = np.random.default_rng(randy)

        # Vary seed with next simulation.
        if switch > 0:
            scfg.SEEDX += 100

    return rng


def check_if_fault_exists(fault_controls, rng):
    """Test if a fault exists for this simulation.

    Parameters
    ----------
    fault_controls = (dict) fault control values
    rng = (int) random function

    Returns
    -------
    status = (bool) fault status / existence of fault => True/False
    """
    # Check if the existence of fault is probabilistic.
    if fault_controls['fault_probability'] >= FAULT_LIMIT:
        status = True
    else:
        # Setup lists for determination.
        bool_list = ['True', 'False']
        present = fault_controls['fault_probability']
        not_present = 100.0 - present

        # Create list of >k=1< random Booleans.
        result_list = rng.random.choices(bool_list,
                                         weights=(present, not_present), k=1)
        # Choose element as result.
        status = (result_list[0] == "True")

    return status


def compute_sgr_influence(fault_controls):
    """Compute the effects of SGR on entry pressure & permeability.

    Parameters
    ----------
    fault_controls = (dict) fault control values

    Returns
    -------
    fault_controls = (dict) with updated control

    Notes
    -----
    1. Only for normal faulting.
    2. Simple models; assumes shale layers have significant clay content.
    3. Assume 20% SGR at start of blockage effects.
    4. SGR varies from 0.0 to 100.0 (percent).
    5. References:
        - Yielding, G.; Freeman, B.; Needham, D.T. (1997)
        - Manzocchi, T.; Walsh, J. J.; Nell, P.; Yielding, G. (1999)
        https://www.researchgate.net/publication/237231297
        _Fault_transmissibility_multipliers_for_flow_simulation_models
        - Benguigui, A.; Yin, Z.; MacBeth, C. (2013)
        https://doi.org/10.1088/1742-2132/11/2/025006
        - Bourg, I. C.; Ajo-Franklin, J. B. (2017)
               https://pubs.acs.org/doi/10.1021/acs.accounts.7b00261
    """
    # Define terms for clarity.
    sgr = fault_controls['sgr']

    # Set corrections if SGR is in range.
    if SGR_MINIMUM <= sgr <= SGR_MAXIMUM:
        # compute new entry pressure due to SGR (assigned to each plate).
        bot = math.pow(10, SGR_PRESS)
        top = math.pow(10, (SGR_SLOPE * sgr))
        press_factor = (top / bot)

        # Use log relationship to find correction factor on permeability.
        val = sgr / SGR_MINIMUM
        log_val = math.log10(val)
        log_val *= SGR_CONSTANT
        new_perm = math.pow(10, log_val)
    else:
        new_perm = 1.0
        press_factor = 1.0

    fault_controls['sgr_press_correction'] = press_factor
    fault_controls['sgr_perm_correction'] = new_perm

    return fault_controls


def setup_apertus(fault_controls):
    """Provide setup parameters for aperture computation.

    Parameters
    ----------
    fault_controls = (dict) input parameters

    Returns
    -------
    fault_controls = (dict) revised parameters
    """
    # Define model default/base aperture parameters.
    fault_controls['pa_limit_stress'] = LIMIT_STRESS
    fault_controls['pa_frac_stiffness'] = STIFFNESS
    fault_controls['pa_initial_max'] = MAX_APERTURE
    fault_controls['pa_initial_residual'] = RESIDUAL_APERTURE
    average_press = (fault_controls['field_pressure']
                     + fault_controls['aquifer_pressure']) / 2.0
    fault_controls['pa_initial_stress'] = average_press
    fault_controls['pa_initial_aperture'] = fault_controls['aperture_mean']

    return fault_controls


def compute_new_gap(fault_controls, stress, gamma_ref, beta, vary):
    """Compute new aperture using press-aperture relationship.

    Parameters
    ----------
    fault_controls = (dict) revised parameters
    stress = (float) current stress value for effect. normal stress / pressure
    gamma_ref = (float) reference gamma for specified stress/aperture
    beta = (float) beta value for current fracture
    vary = (float) variable aperture

    Returns
    -------
    new_aperture = (float) gamma function
    """
    stiffness = fault_controls['pa_frac_stiffness']
    gamma_current = 1.0 - beta * (stress / (stress + stiffness))
    gamma_delta = gamma_current - gamma_ref
    initial_aperture = fault_controls['pa_initial_aperture']
    new_aperture = gamma_delta * vary + initial_aperture

    return new_aperture


def define_apertus(fault_controls):
    """Provide setup parameters for pressure-aperture computation.

    Parameters
    ----------
    fault_controls = (dict) input parameters

    Returns
    -------
    fault_controls = (dict) revised parameters
    """
    # Define default parameters.
    fault_controls = setup_apertus(fault_controls)

    # Define shortcuts for use in equation.
    stress = fault_controls['pa_initial_stress']
    max_aperture = fault_controls['pa_initial_max']
    residual_aperture = fault_controls['pa_initial_residual']
    stiffness = fault_controls['pa_frac_stiffness']
    limit_stress = fault_controls['pa_limit_stress']

    # Define beta, extreme stress and aperture terms.
    beta = (limit_stress + stiffness) / limit_stress
    gamma_ref = (1.0 - beta * (stress / (stress + stiffness)))
    vary_aperture = max_aperture - residual_aperture

    # Compute new maximum and minimum apertures.
    # -- At max aperture, stress=0.0 and gamma_max = 1.0.
    new_max_gap = compute_new_gap(fault_controls, 0.0, gamma_ref, beta,
                                  vary_aperture)
    # -- At min aperture, stress=limit stress.
    new_residual_gap = compute_new_gap(fault_controls, limit_stress,
                                       gamma_ref, beta, vary_aperture)

    # Define fault_controls for future aperture calculations.
    fault_controls['pa_beta'] = beta
    fault_controls['pa_gamma_ref'] = gamma_ref
    fault_controls['pa_vary_aperture'] = vary_aperture
    fault_controls['pa_max_aperture'] = new_max_gap
    fault_controls['pa_residual_aperture'] = new_residual_gap

    return fault_controls


def define_fracto(fault_controls, param_bounds):
    """Provide the setup fault input values.

    Parameters
    ----------
    fault_controls = (dict) input parameters
    fault= (class) list of segments
    param_bounds = (dict) parameter bounds (dictionary)

    Returns
    -------
    fault_controls = (dict) revised
    """
    # Define log-normal parameters to compute apertures.
    if fault_controls['aperture_approach']:
        # No variability - data from file.
        fault_controls['aperture_mu'] = 0.0
        fault_controls['aperture_scale'] = 0.0
    else:
        # Setup parameters for variable aperture.
        cent, scale = \
            bcs.convert_lognorm_terms(fault_controls['aperture_mean'],
                                      fault_controls['aperture_std'])
        fault_controls['aperture_mu'] = cent
        fault_controls['aperture_scale'] = scale

    # Setup disperse for strike if needed.
    if fault_controls['strike_approach'] and \
            fault_controls['strike_sig'] > V_SMALL_ZERO:
        # Define variability.
        fault_controls['strike_disperse'] = \
            bcs.convert_kappa(fault_controls['strike_sig'])

    # Compute effects of SGR.
    fault_controls = compute_sgr_influence(fault_controls)

    # Input apertures from file - if desired.
    fault_controls = fio.input_apertures(fault_controls, param_bounds)

    # Define fault counter.
    fault_controls['simulations_with_fault'] = 0

    # Define aperture-pressure parameters, if they will be used.
    if fault_controls['pressure_approach']:
        fault_controls = define_apertus(fault_controls)

    return fault_controls


def project_to_subsurface(fault_controls):
    """Project fault start point on surface to subsurface.

    Parameters
    ----------
    fault_controls = (dict) fault control values

    Returns
    -------
    x_new = x coordinate of point in subsurface
    y_new = y coordinate of point in subsurface

    Notes
    -----
    "strike" may vary each simulation.
    """
    # Define variables and convert to radians.
    x_surface = fault_controls['x_start']
    y_surface = fault_controls['y_start']
    alpha = math.radians(fault_controls['dip'])
    gamma = math.radians(fault_controls['strike'])
    depth = fault_controls['inject_depth']

    # Incline length - use complement.
    beta = math.radians(90.0 - fault_controls['dip'])
    if beta < 90.0:
        # Compute inclination.
        incline = depth / math.cos(beta)
    else:
        # Error - division by zero!
        incline = 0.0
        fileop.opx_problem("Error: Dip Angle is Too Large!")

    # Along dip extent - local x coordinate.
    x_prime = incline * math.cos(alpha)

    # Rotate from local coordinates to N/E.
    delta_x = math.cos(gamma) * x_prime
    delta_y = -math.sin(gamma) * x_prime

    # Compute new coordinates.
    x_new = x_surface + delta_x
    y_new = y_surface + delta_y

    return x_new, y_new


def create_features(fault_controls, fault, x_start, y_start):
    """Create a fault (list of plates) with default values & coordinates.

    Parameters
    ----------
    fault_controls = (dict) fault control values
    fault = (class) list of plates
    x_start= (float) X coordinate of subsurface start point
    y_start= (float) Y coordinate of subsurface start point

    Returns
    -------
    fault = (list) fault segment/plates
    """
    # Define computation parameters.
    n_segs = fault_controls['n_plates']
    long = fault_controls['length'] / n_segs
    x_valu = x_start
    y_valu = y_start

    # Generate new points and new fault plates - save in "fault" list.
    for _ in range(n_segs):
        point_1 = model.Pointe(x_valu, y_valu)
        point_2 = point_1.arthro(fault_controls['strike'], long)

        new_plate = model.Plate(point_1, point_2,
                                fault_controls['strike'],
                                fault_controls['dip'])
        fault.append(new_plate)
        x_valu, y_valu = point_2.get_coords()

    # Define fault end point.
    fault_controls['fault_end'] = model.Pointe(x_valu, y_valu)

    # Debug.
    if ECHO_2:
        print("\n start sim")
        for part in fault:
            part.print_plate(sys.stdout)

    return fault


def stress_factor(fault_controls, old_strike, new_strike):
    """Compute the stress factor ratio for new strike based on old.

    Parameters
    ----------
    fault_controls = (dict) fault control values
    old_strike = (float) original strike of fault
    new_strike = (float) new strike of fault

    Returns
    -------
    stress ratio = (float) ratio of normal stresses to fault
    """
    # Rotate to get stress perpendicular to fault strike.
    # >> Note that is approximate as using horizontal components only!
    old_strike += 90.0
    new_strike += 90.0

    # Get ratio of new to old stress as correction factor.
    major = fault_controls['max_horizontal']
    minor = fault_controls['min_horizontal']
    stress_angle = fault_controls['max_trend']

    radius_old = bcs.radius_ellipse(major, minor, stress_angle, old_strike)
    radius_new = bcs.radius_ellipse(major, minor, stress_angle, new_strike)
    stress_ratio = radius_new / radius_old

    return stress_ratio


def define_aspects(fault_controls, fault, rng):
    """Define parameters for fault plates.

    Parameters
    ----------
    fault_controls = (dict) fault control values
    fault = (class) list of segments/plates
    file_gaps = (array) apertures from file (if desired)
    rng = (int) random function

    Returns
    -------
    fault_controls = (dict) fault control values revised
    """
    # Compute correction for stress field.
    normal_factor = stress_factor(fault_controls,
                                  fault_controls['strike_mean'],
                                  fault_controls['strike'])

    # Define terms for clarity.
    nsegs = fault_controls['n_plates']
    gaps = np.empty(nsegs)

    # Define apertures for plates - depending on option selected.
    if fault_controls['aperture_approach']:
        # Take apertures from file values - no correction for stress field.
        gaps = fault_controls['aperture_data']

    elif fault_controls['vary_aperture']:
        # Vary aperture and correct for stress field.
        for i in range(nsegs):
            gaps[i] = \
                bcs.evaluate_lognormal(fault_controls['aperture_mu'],
                                       fault_controls['aperture_scale'],
                                       fault_controls['aperture_min'],
                                       fault_controls['aperture_max'],
                                       rng)
            gaps[i] *= normal_factor
    else:
        # Use single value for all - mo correction.
        valu = fault_controls['aperture_mean']
        gaps.fill(valu)

    # -----------
    # Assign an aperture to each fault plate + remaining parameters.
    for index, part in enumerate(fault):

        # Set plate aperture in mm
        part.aperture = gaps[index]

        # Compute Permeability in meters.
        #   cubic law => k [permeability] - in m^2
        corrected = gaps[index] * funit.mm_to_m()
        squared = math.pow(corrected, 2)
        part.perm = squared / 12.0

        # Entry pressure.
        part.entry = fault_controls['entry_pressure']

        # Debug.
        if ECHO:
            # part.print_plate(sys.stdout)
            pass

    return fault


def define_amplius(fault_controls, fault, rng):
    """Complete fault input values and fault plates.

    Parameters
    ----------
    fault_controls = (dict) fault control values
    fault = (list) list of segments/plates
    rng = (int) random function

    Returns
    -------
    fault_controls = (revised)
    """
    # Define new strike, if desired.
    if fault_controls['strike_approach'] and \
            fault_controls['strike_sig'] > V_SMALL_ZERO:
        fault_controls['strike'] = \
            bcs.evaluate_orientation(fault_controls['strike_mean'],
                                     fault_controls['strike_disperse'],
                                     rng)
    else:
        fault_controls['strike'] = fault_controls['strike_mean']

    # Define new dip, if desired.
    if fault_controls['dip_approach'] and \
            fault_controls['dip_std'] > V_SMALL_ZERO:
        fault_controls['dip'] = \
            bcs.evaluate_norm(fault_controls['dip_mean'],
                              fault_controls['dip_std'],
                              DIP_MIN, DIP_MAX, rng)
    else:
        fault_controls['dip'] = fault_controls['dip_mean']

    # Compute subsurface coordinates.
    x_valu, y_valu = project_to_subsurface(fault_controls)
    fault_controls['fault_start'] = model.Pointe(x_valu, y_valu)

    # Create fault plates.
    fault = create_features(fault_controls, fault, x_valu, y_valu)

    # Assign aperture and other parameters to fault plates.
    fault = define_aspects(fault_controls, fault, rng)

    return fault_controls, fault


def compute_fault_length(fault_controls):
    """Compute the deeper travel profile length along fault for Darcy Flow.

    Parameters
    ----------
    fault_controls = (dict) fault control values

    Returns
    -------
    fault_controls = (dict) fault control values with updated controls

    Notes
    -----
    1. Assume no tortuosity / channelization on profile.
    2. Depths are positive values.
    3. top = supercritical transition depth for fault.
    """
    # Default lengths.
    start_depth = fault_controls['inject_depth']
    upper_depth = start_depth    # To eliminate error message

    # Transition is dependent on model type.
    if fault_controls['profile_type'] == 0:
        # Model 0: aquifer depth.
        upper_depth = fault_controls['aquifer_depth']

    elif fault_controls['profile_type'] == 1:
        # Model 1: Use 10 MPa depth as base point.
        upper_depth = fault_controls['trans_point_depth']
        if upper_depth > fault_controls['inject_depth']:
            max_depth = fault_controls['inject_depth']
            msg = 'Input Error: Model <1> Cutoff is Deeper than Injection '
            msg += f'Depth of {max_depth}'
            fileop.opx_problem(msg)

    elif fault_controls['profile_type'] == 2:
        # Model 2: Supercritical depth.
        upper_depth = funit.supercrit_press_depth()

    else:
        msg = ("Code Error: " +
               "No Fault Model Identified in <compute_fault_length>")
        fileop.opx_problem(msg)

    # Compute parameters for calculation.
    delta_depth = start_depth - upper_depth
    theta = fault_controls['dip']
    total_length = 0.0    # To eliminate error message

    # Find inclination and then length.
    if theta < V_SMALL_ANGLE:
        msg = 'Input Error: Fault is Near-Horizontal with a Dip of '
        msg += f'{fault_controls["dip"]}'
        fileop.opx_problem(msg)
    else:
        total_length = delta_depth / math.sin(math.radians(theta))

    # Set length as new fault_control.
    fault_controls['travel_distance'] = total_length

    return fault_controls


def fault_plane_eq(fault_controls):
    """Compute the vector equation for the normal to a fault.

    Parameters
    ----------
    fault_controls = (dict) fault control values

    Returns
    -------
    plane_new = (array) np_vector for plane with X-Y-Z components.

    Notes
    -----
    1. Direction cosines are converted from East-North-Down (NED)
        Use Right-Hand Rule (RHR)
        v_north = north component
        v_east = east component
        v_down = z (down) component
    2. XYZ is a RHR East-North-Down (ENU) system
    3. see:
        https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470099728.app3
    """
    # Define angles.
    phi = math.radians(fault_controls['strike'])
    theta = math.radians(fault_controls['dip'])

    # compute direction cosines (Ref #3).
    v_north = -math.sin(phi) * math.sin(theta)
    v_east = math.cos(phi) * math.sin(theta)
    v_down = -math.cos(theta)

    # Convert from NED to ENU (XYZ) RHS coordinate system (Ref. #1, #2).
    vector_x = v_east
    vector_y = v_north
    vector_z = -v_down

    # Compute "D" value using start point (z = 0.0).
    d_value = -(vector_x * fault_controls['x_start']
                + vector_y * fault_controls['y_start']
                + vector_z * 0.0)

    # Establish vector.
    plane_new = np.array([vector_x, vector_y, vector_z, d_value])

    return plane_new


def test_oblique(planar, focus, fault_controls):
    """Check if closest point to well is within fault segment limits.

    Parameters
    ----------
    planar = (array) vector values for plane
    focus = (array) well location coordinates (3D)
    fault_controls = (dict) fault control values

    Returns
    -------
    direct = (bool) flag if short point within fault extent.

    Notes
    -----
    1. Indicator that fault is facing/near well
        https://mathworld.wolfram.com/Point-PlaneDistance.html
    """
    # Define denominator of computation - value not squared!
    direct = False
    denominator = (planar[0] * planar[0]
                   + planar[1] * planar[1]
                   + planar[2] * planar[2])

    # Compute lambda (eta) factor and determine resulting point on plane.
    if denominator > 0.0:
        numerator = -(planar[0] * focus[0]
                      + planar[1] * focus[1]
                      + planar[2] * focus[2]
                      + planar[3])
        eta = numerator / denominator

        # Compute coordinates of near point on plane.
        normal_pt = np.zeros(3)
        for i in range(3):
            normal_pt[i] = focus[i] + eta * planar[i]

        # Define fault start and end points in subsurface.
        pt_start = fault_controls.get('fault_start')
        x_pt1, y_pt1 = pt_start.get_coords()
        pt_end = fault_controls.get('fault_end')
        x_pt2, y_pt2 = pt_end.get_coords()

        # Check if well location is within fault extent (in subsurface).
        # + EXTRA_D as a margin.
        if (min(x_pt1, x_pt2) - EXTRA_D) < normal_pt[0] \
                < (max(x_pt1, x_pt2) + EXTRA_D):
            if (min(y_pt1, y_pt2) - EXTRA_D) < normal_pt[1] \
                    < (max(y_pt1, y_pt2) + EXTRA_D):
                direct = True

        # Debug.
        if ECHO:
            print('  - Normal point = ', end='')
            print(normal_pt)

    return direct


def distance_vector(planar, well_pt):
    """Compute vector for shortest distance from injection well to fault.

    Parameters
    ----------
    planar = (array) vector values for plane (XYZ system)
    depth = (float) depth of injection point (XYZ System)

    Returns
    -------
    distance = (float) the shortest distance from well point to plane.
    new_vector = (array) vector from liner to plane (in XYZ)

    Notes
    -----
    1. Right Hand Rule, North-East-Down Convention
        https://mathworld.wolfram.com/Point-PlaneDistance.html
    """
    # Define denominator of distance computation without "d".
    denominator = math.sqrt(planar[0] * planar[0]
                            + planar[1] * planar[1]
                            + planar[2] * planar[2])

    # Check for zero divisor.
    if denominator > 0.0:
        numerator = (planar[0] * well_pt[0]
                     + planar[1] * well_pt[1]
                     + planar[2] * well_pt[2]
                     + planar[3])
        distance = abs(numerator / denominator)

        # Debug.
        if ECHO:
            val = numerator / denominator
            print(f'\n  - Distance = {val:.5e} m')

    else:
        distance = 0.00

    # ---------------------
    # Determine the vector from point to plane in XYZ space.
    vector_new = np.zeros(3)
    for i in range(3):
        vector_new[i] = -(planar[i] - well_pt[i])

    # Debug.
    if ECHO:
        print('  - Perpendicular Point Vector = ', end='')
        print(vector_new)

    return distance, vector_new


def is_almost_equal(x_val, y_val, epsilon=1.0e-08):
    """Return True if two float values are close within tolerance.

    Parameters
    ----------
    x_val = (float) one value for comparison (reference)
    y_val = (float) second value for comparison
    epsilon = (float) allowed difference

    Return
    ------
    result = (boolean) status if within tolerance

    Note
    ----
    1. Default tolerance is 1.0E-08.
    2. Values can be negative or positive.
    """
    result = math.isclose(x_val, y_val, abs_tol=epsilon)
    return result


def check_inclination(fault_controls):
    """Determine the fault inclination wrt to the well and set factor.

    Parameters
    ----------
    fault_controls = (dict) fault control values

    Returns
    -------
    fault_controls = (dict) fault control values with updated control
    """
    # Establish plane and well point in vector notation (NED).
    plane_vector = fault_plane_eq(fault_controls)
    well_pt = [fault_controls.get('inject_x'),
               fault_controls.get('inject_y'),
               -fault_controls.get('inject_depth')]

    # Determine distance and distance vector from well to plane.
    extent, well_vector = distance_vector(plane_vector, well_pt)
    fault_controls['shortest_distance'] = extent

    # Determine angle of plane normal wrt to well distance vector.
    factor, included_angle = bcs.angle_vector(plane_vector, well_vector)
    fault_controls['angle_between_vectors'] = included_angle

    # Determine if shortest distance point is on/near fault plate extent.
    nearby = test_oblique(plane_vector, well_pt, fault_controls)

    # Determine plane inclination condition and effect factor.
    #   Note for factor ~1.0, vector is perpendicular to plane
    if (is_almost_equal(factor, 1.0, V_SMALL_ONE)) and nearby:
        fault_controls['fault_inclined'] = \
            "Fault Normal Is Parallel to Projection Line."
        fault_controls['fault_orient_effect'] = 0.90
    elif factor > 0.0 and nearby:
        fault_controls['fault_inclined'] = \
            "Inclined Away From Injection Point."
        fault_controls['fault_orient_effect'] = 0.90
    elif factor < 0.0 and nearby:
        fault_controls['fault_inclined'] = \
            "Inclined Towards Injection Point."
        fault_controls['fault_orient_effect'] = 1.10
    else:
        fault_controls['fault_orient_effect'] = 1.00
        if not nearby:
            fault_controls['fault_inclined'] = \
                "Inclination is Oblique wrt Well."
        else:
            # Error.
            fault_controls['fault_inclined'] = \
                "Inclination Not Determined!"

    # Debug.
    if ECHO:
        print('  - Inclination Condition: '
              + fault_controls['fault_inclined'])

    return fault_controls


def fault_create(fault_controls, fault, rng):
    """Manage the fault setup/plate creation.

    Parameters
    ----------
    fault_controls = (dict) fault control values
    fault = (class) fault (with segments)
    rng = (int) random function

    Returns
    -------
    fault_controls = (dict) fault control values updated values
    fault = (class) fault plates
    """
    # Define faults.
    fault.clear()
    fault_controls, fault = \
        define_amplius(fault_controls, fault, rng)

    # Initialize: Compute travel length along fault - based on model type.
    fault_controls = compute_fault_length(fault_controls)

    # Determine inclination wrt to injection well.
    fault_controls = check_inclination(fault_controls)

    return fault_controls, fault


#
# -----------------------------------------------------------------------------
# End of module
