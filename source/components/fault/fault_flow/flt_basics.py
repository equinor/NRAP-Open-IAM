#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains basic math computations for fault analyses.

Author: Ernest N. Lindner
Date: 08/16/2022

Module Name
    flt_basics

Contents (9)
    radius_ellipse(major, minor, trend, theta):
    convert_lognorm_terms(mean_val, std_dev)
    evaluate_lognormal(ap_mu, ap_scale, minimal, maximal, rng)
    convert_kappa(fault_controls)
    evaluate_orientation(alpha, kappa, rng)
    find_line(point_x, point_y, zpoint_x, zpoint_y) -NOT USED
    find_cross(fault_controls) -NOT USED
    angle_vector(planar, liner)
    evaluate_norm(tn_mean, tn_sigma, tn_min, tn_max, rng)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import math                         # For sin and other functions
import numpy as np                  # For array of values
from scipy.stats import truncnorm   # For a truncated normal function

import flt_units as funit           # For unit conversion
import flt_file as fileop           # For file operations

# Other constants
DEG_TO_KAPPA = 0.00875          # Conversion factor for 2-sigma deg. to Kappa
V_SMALL_ZERO = 1.0E-08          # Very small number in division by zero
ECHO = False                    # Debugging


def radius_ellipse(major, minor, trend, theta):
    """Compute the radius for a stress ellipse at angle, in polar notation.

    Parameters
    ----------
    major = (float) value of major axis (Pa)
    minor = (float) value of minor axis (Pa)
    trend = (float) trend of major stress axis (degrees) clockwise from North
    theta = (float) strike of new radius (degrees), clockwise from North

    Returns
    -------
    radius = (float) length of radius

    Notes
    -----
    1) see polar form: https://en.wikipedia.org/wiki/Ellipse
    """
    # Get angle of strike relative to stress ellipse.
    diff_angle = theta - trend
    omega = math.radians(diff_angle)

    # polar equations of ellipse - use MPa for clarity.
    x_val = major * math.cos(omega) * funit.pa_to_mpa()
    y_val = minor * math.sin(omega) * funit.pa_to_mpa()
    radius = (major * minor) / math.hypot(x_val, y_val)

    return radius


def convert_lognorm_terms(mean_val, std_dev):
    """Convert mean/std_dev into log-normal parameters for NumPy.

    Parameters
    ----------
    mean_val = (float) mean permeability in microdarcys
    std_dev = (float) std deviation of permeability

    Returns
    -------
    log_mu = (float) log-normal mean
    sigma = (float) log-normal variance

    Notes
    -----
    See: https://en.wikipedia.org/wiki/Log-normal_distribution
    """
    # define common terms.
    variance = math.pow(std_dev, 2)
    mean_2 = math.pow(mean_val, 2)
    inside = 1.0 + (variance / mean_2)

    # Compute location - mean.
    term = math.sqrt(inside)
    log_mu = math.log(mean_val / term)

    # Compute sigma - variance.
    term = math.log(inside)
    sigma = math.sqrt(term)

    return log_mu, sigma


def evaluate_lognormal(ap_mu, ap_scale, minimal, maximal, rng):
    """Define stochastic aperture of a fault.

    Parameters
    ----------
    ap_mu = (float) location of log-distribution
    ap_scale = (float) scale of log-distribution
    minimal = (float) minimum aperture in mm
    maximal = (float) maximum aperture in mm
    rng = random number function

    Returns
    -------
    new_value = (float) aperture

    Notes
    -----
    Values are censored to be within the defined min. and max.
    """
    # Use log-normal distribution to get value within limits.
    while True:
        new_value = rng.lognormal(mean=ap_mu, sigma=ap_scale,
                                  size=None)
        if minimal <= new_value <= maximal:
            break

    return new_value


def convert_kappa(sigma):
    """Define kappa for a von Mises orientation distribution.

    Parameters
    ----------
    sigma = (float) the spread in strike - in degrees

    Returns
    -------
    disperse = (float) kappa value for strike

    Notes
    -----
    1. sigma is the spread in degrees - must be converted to kappa.
    2. Kappa:
        kappa is "concentration" -> larger k is a smaller spread.
        Rough empirical fit for kappa = (sigma*DEG_TO_KAPPA)^-2;
        sigma should be > 0.1 for numerical stability of kappa.
    """
    if sigma < 0.1:
        sigma = 0.1 * DEG_TO_KAPPA
    else:
        sigma *= DEG_TO_KAPPA

    disperse = math.pow(sigma, -2.0)

    return disperse


def evaluate_orientation(alpha, kappa, rng):
    """Define random orientation of a fault.

    Parameters
    ----------
    alpha = (float) average value (degrees)
    kappa = (float) term for dispersion of distribution of von Mises
    rng = (int) random function

    Returns
    -------
    results = (float) orientation (degrees) >> 0 to 360 from North

    Notes
    -----
    1. ""unction: random.vonmisis => on the interval [-pi, pi].
    2. Kappa:
        kappa is inverse of spread!; (Large kappa reduces spread of values);
        kappa = 0 is uniform across range (spread ~ across range);
        kappa ~ 1.0E+5 is spread of less than 1 deg => so no distribution.
        (kappa ~ 1.0E+4 for parallel lines (spread ~0) )
    """
    # Convert degrees to radians for math.
    trend = math.radians(alpha)
    limit = math.pi * 2.0

    # Use von Mises distribution to define orientation;
    #   --> Watch for limiting cases for low and high kappa values.
    if kappa <= 0.01:
        results = rng.uniform(0.0, limit)
    elif kappa > 1.0E+6:
        results = trend
    else:
        # Use Fisher-von Mises distribution.
        val_rad = rng.vonmises(trend, kappa)
        results = math.degrees(val_rad)

    return results


def find_line(point_x, point_y, zpoint_x, zpoint_y):
    """Find the slope and intercept "(m, b) of line.

    Parameters
    ----------
    point_x = (float) point-1 X Coordinate
    point_y = (float) point-1 Y Coordinate
    zpoint_x = (float) point-2 X Coordinate
    zpoint_y = (float) point-2 Y Coordinate

    Returns
    -------
    incline = (float) slope
    beta = (float) constant

    Notes
    -----
    1. Equation: <y = mx + b>; where m = slope and b = beta.
    """
    delta_x = zpoint_x - point_x
    delta_y = zpoint_y - point_y

    if abs(delta_x) < V_SMALL_ZERO:
        incline = math.inf
        beta = point_y
    else:
        incline = delta_y / delta_x
        if abs(incline) < V_SMALL_ZERO:
            incline = 0.0
        beta = point_y - incline * point_x

    return incline, beta


def find_cross(fault_controls):
    """Find the intersection of site gradient with Triple Line.

    Parameters
    ----------
    fault_controls = (dic) dictionary of control values

    Returns
    -------
    x_new = (float) intersection X-coordinate
    y_new = (float) intersection Y-coordinate

    Notes
    -----
    For theory, see URLs:
        1. https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
    """
    # Define variables for clarity from tuple.
    # --> Gradient line.
    pt_x1 = fault_controls['aquifer_temperature']
    pt_y1 = fault_controls['aquifer_pressure']
    pt_x2 = fault_controls['inject_temperature']
    pt_y2 = fault_controls['field_pressure']
    # --> Triple line.
    pt_x3 = funit.triple_temperature()
    pt_y3 = funit.triple_pressure()
    pt_x4 = funit.supercrit_temperature()
    pt_y4 = funit.supercrit_pressure()

    # Compute determinant.
    determinant = ((pt_x4 - pt_x3) * (pt_y1 - pt_y2)
                   - (pt_x1 - pt_x2) * (pt_y4 - pt_y3))

    # Compute extension numerators.
    extend_a = ((pt_y3 - pt_y4) * (pt_x1 - pt_x3) +
                (pt_x4 - pt_x3) * (pt_y1 - pt_y3))
    extend_b = ((pt_y1 - pt_y2) * (pt_x1 - pt_x3) +
                (pt_x2 - pt_x1) * (pt_y1 - pt_y3))

    # get extension but check for division by zero.
    if abs(determinant) < V_SMALL_ZERO:
        # Parallel - no intersection! - provide error code.
        msg1 = "Error - No Intersection Found!"
        fileop.opx_problem(msg1)
    else:
        # Get extensions for lines - divide by determinant.
        extend_a /= determinant
        extend_b /= determinant

    # Compute intersection point coordinates.
    x_new = pt_x1 + extend_a * (pt_x2 - pt_x1)
    y_new = pt_y1 + extend_a * (pt_y2 - pt_y1)

    return x_new, y_new


def angle_vector(planar, liner):
    """Compute the dot product of injection well to fault and fault normal.

    Parameters
    ----------
    planar = (list) list of vector values for plane - vector 1
    liner = (list) well vector (XYZ system) - vector 2

    Returns
    -------
    dot_product = (float) dot product of planar and liner
    angle_deg = (float) angle between vectors, in degrees

    Notes
    -----
    1. Takes a, b, c of plane as normal to plane.
    2. Right-Hand Rule, North-East-Down Convention.
    3. See URL:
        1. https://mathworld.wolfram.com/Point-PlaneDistance.html
    """
    # Create unit vectors.
    vector_1 = np.array([planar[0], planar[1], planar[2]])
    vector_2 = np.array([liner[0], liner[1], liner[2]])

    unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
    unit_vector_2 = vector_2 / np.linalg.norm(vector_2)

    # Find angle between two products in XYZ system.
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    angle_deg = math.degrees(np.arccos(dot_product))

    # Debug.
    if ECHO:
        print(f'  - Dot Product = {dot_product:.5e}')

    return dot_product, angle_deg


def evaluate_norm(tn_mean, tn_sigma, tn_min, tn_max, rng):
    """Get one random value from a truncated normal distribution.

    Parameters
    ----------
    tn_mean = (float) mean of distribution
    tn_sigma = (float) std deviation of distribution
    tn_min = (float) minimum of distribution
    tn_max =  (float) maximum of distribution
    rng = (int) random functions
    sizer = (int) number of values

    Returns
    -------
    value = (float, list) truncated value(s)
    """
    # Obtain value from truncated normal distribution.
    sizer = 1   # One value
    limit_1 = (tn_min - tn_mean) / tn_sigma
    limit_2 = (tn_max - tn_mean) / tn_sigma
    valu = truncnorm.rvs(a=limit_1, b=limit_2, loc=tn_mean, scale=tn_sigma,
                         size=sizer, random_state=rng)

    return valu


#
# -----------------------------------------------------------------------------
# End of module
