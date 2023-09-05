#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions for random stochastic computations.

Author: Ernest N. Lindner
Date: 08/18/2022

Module Name
    frac_random

Contents (10)
    convert_lognorm_terms(mean_val, std_dev)
    correlate_aperture(alpha, beta, frac_length, rng)
    evaluate_aperture(mu_d, ap_scale, minimal, maximal, rng)
    evaluate_density(ave_valu, min_valu, max_valu, rng)
    evaluate_power_length(eta, minimal, maximal, rng)
    evaluate_location(x_min, x_max, y_min, y_max, rng)
    evaluate_matrix_perm(mu_d, ap_scale, minimal, maximal, rng)
    evaluate_orientation(strike, kappa, rng)
    rand_power_law(eta, k_min, k_max, rng, nodes=1)
    evaluate_lognorm(mu_d, perm_scale, minimal, maximal, rng)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import math                         # For pi, sqrt, pow, radians
from scipy.stats import powerlaw    # For power law function, a>0
import numpy as np                  # For random functions


# Constants
LIMIT_UP = 1.0E-10                  # Just above zero for power law
LIMIT_DN = -1.0E-10                 # Just below zero for power law


def convert_lognorm_terms(mean_val, std_dev):
    """Convert mean/std_dev into Normal parameters for NumPy.

    Parameters
    ----------
    mean_val = (float) mean of data (log-normal)
    std_dev = (float) std deviation of data (log-normal)

    Returns
    -------
    mu_dist = (float) Normal mean
    sigma = (float) Normal standard deviation

    Notes
    -----
    1. variance = (standard deviation)^2
    2. mean of log-normal is <<not> typically the mode (peak value)
    3. For Equations:
        see: https://www.mathworks.com/help/stats/lognstat.html
        see: https://en.wikipedia.org/wiki/Log-normal_distribution
        (Estimation of parameters)
    """
    # define common terms.
    variance = math.pow(std_dev, 2)
    mean_2 = math.pow(mean_val, 2)
    inside = 1.0 + (variance / mean_2)

    # Compute location - mu - mean.
    mu_dist = math.log(mean_val / math.sqrt(inside))

    # Compute sigma - variance.
    sigma = math.sqrt(math.log(inside))

    return mu_dist, sigma


def correlate_aperture(alpha, beta, frac_length, rng):
    """Define aperture based on length correlation.

    Parameters
    ----------
    alpha = (float) exponent of correlation
    beta = (float) slope of correlation
    frac_length = (float) fracture length
    rng = (generator) random number generator

    Returns
    -------
    aperture_new = new aperture

    Notes
    -----
    1. Correlation based on simple power law on length (m)
        versus aperture (mm), b = beta * length^alpha.
    2. Due to high variability, values are considered
        to be within the defined minimum and maximum.
    """
    # Define limits for uncertainty on beta.
    minimum = 0.65 * beta
    maximum = 1.50 * beta

    # Compute new beta with uniform distribution of variability.
    beta_new = rng.uniform(minimum, maximum)

    # Define aperture using power law correlation.
    aperture_new = beta_new * math.pow(frac_length, alpha)

    return aperture_new


def evaluate_aperture(mu_d, ap_scale, minimal, maximal, rng):
    """Define aperture of a fracture.

    Parameters
    ----------
    mu_d = (float) location of log-distribution
    ap_scale = (float) scale of log-distribution
    minimal = (float) minimum aperture in mm
    maximal = (float) maximum aperture in mm
    rng = (generator) random number generator

    Returns
    -------
    new_value = (float) aperture

    Notes
    -----
    1. Values are censored to be within the defined min. and max.
    """
    # Use log-normal distribution to get values.
    while True:
        new_value = rng.lognormal(mean=mu_d, sigma=ap_scale, size=None)
        if minimal <= new_value <= maximal:
            break

    return new_value


def evaluate_density(ave_valu, min_valu, max_valu, rng):
    """Define general density of fractures over an area.

    Parameters
    ----------
    ave_valu = (float) average value
    min_valu = (float) minimum value
    max_valu = (float) maximum value
    rng = (generator) random number generator

    Returns
    -------
    density = density value
    """
    # Define density for analysis using triangular distribution.
    density = rng.triangular(min_valu, ave_valu, max_valu, size=1)

    return density


def evaluate_power_length(eta, minimal, maximal, rng):
    """Define line length using power law distribution.

    Parameters
    ----------
    eta = (float) power law exponent
    minimum = (float) minimum vertical permeability in μD
    maximum = (float) maximum vertical permeability in μD
    rng = (generator) random number generator

    Returns
    -------
    new_value = length value

    Notes
    -----
    1. Values are censored to be within the defined min. and max.
    """
    # Use power law distribution to get value.
    new_value = rand_power_law(eta, minimal, maximal, rng, nodes=1)

    return new_value


def evaluate_location(x_min, x_max, y_min, y_max, rng):
    """Define location coordinates within a rectangular region.

    Parameters
    ----------
    x_min = (float) minimum x / width bound of region
    x_max = (float) maximum x / width bound of region
    y_min = (float) minimum y / length bound of region
    y_max = (float) maximum y / length bound of region
    rng = (generator) random number generator

    Returns
    -------
    x_coord = (float) X coordinate of point in box
    y_coord = (float) Y coordinate of point in box
    """
    # Use random uniform distribution for defining coordinates.
    x_coord = rng.uniform(x_min, x_max)
    y_coord = rng.uniform(y_min, y_max)

    return x_coord, y_coord


def evaluate_matrix_perm(mu_d, ap_scale, minimal, maximal, rng):
    """Define matrix permeability of a grid box.

    Parameters
    ----------
    mu_d = (float) location of log-normal distribution
    ap_scale = (float) scale of log-normal distribution
    minimal = (float) minimum permeability in microdarcys
    maximal = (float) maximum permeability in microdarcys
    rng = (generator) random number generator

    Returns
    -------
    new_value = (float) permeability

    Notes
    -----
    1. Values are censored to be within the defined min. and max.
    """
    # Use log-normal distribution to get values + censor values.
    while True:
        new_value = rng.lognormal(mean=mu_d, sigma=ap_scale, size=None)
        if minimal <= new_value <= maximal:
            break

    return new_value


def evaluate_orientation(alpha, kappa, rng):
    """Define random orientation of a fracture.

    Parameters
    ----------
    alpha = (float) average value (degrees)
    kappa = (float) term for dispersion of distribution of von Mises.
    rng = (generator) random number generator

    Returns
    -------
    results = (float) orientation (radians) >> +/- Pi from North

    Notes
    -----
    1. kappa is inverse of spread; (Large kappa reduces spread of values);
        kappa = 0 is uniform across range (spread ~ across range);
        kappa ~ 1.0E+5 is spread of less than 1 deg => so no distribution.
        kappa ~ 1.0E+4 for parallel lines (spread ~0)
    """
    # Convert degrees to radians.
    trend = math.radians(alpha)
    limit = math.pi

    # Use von Mises distribution to define orientation;
    #   * Watch for limiting cases for low and high kappa values.
    if kappa <= 0.01:
        results = rng.uniform(-limit, limit)
    elif kappa > 1.0E+6:
        results = trend
    else:
        results = rng.vonmises(trend, kappa)

    return results


def rand_power_law(eta, k_min, k_max, rng, nodes=1):
    """Sample value from random power law function.

    Parameters
    ----------
    eta = (float) power law exponent
    k_min = (float) minimum value of variable
    k_max = (float) maximum value of variable
    rng = (generator) random number generator
    nodes = (float) number of values to be generated


    Returns
    -------
    valx = (array) NumPy array of random values

    Notes
    -----
    1. See https://stackoverflow.com/questions/31114330/
        python-generating-random-numbers-from-a-power-law-distribution
    2. Negative exponent: Use inverse method and uniform distribution.
        But check if conditions are met for approach.
    3. Check on distribution limits performed earlier.
    """
    # Set up NumPy array for results.
    valx = np.zeros(nodes, float)

    # Solution depends on exponent - use inverse method for negative.
    if eta <= 0.0:
        # Adjust power value for gamma = n+1 case.
        gamma = eta + 1.00

        if LIMIT_DN < gamma < LIMIT_UP:   # for gamma = zero
            # Use uniform distribution for exponent ~ 0.
            for idx in range(nodes):
                valx[idx] = np.random.uniform(k_min, k_max)
        else:
            # Use inverse method (see reference).
            # -- Define minimum and maximum values of equation.
            pl_min = pow(k_min, gamma)
            pl_max = pow(k_max, gamma)

            # Loop over all nodes.
            for idx in range(nodes):

                # Get value from uniform distribution.
                uniform_function = rng.uniform(0, 1)

                # Compute equation to get value.
                term1 = (pl_max - pl_min) * uniform_function + pl_min
                exponent = 1.0 / gamma
                valx[idx] = pow(term1, exponent)
    else:
        # For positive exponents <eta > 0>, use SciPy function.
        start = k_min
        mod = (k_max - k_min)  # scale
        for idx in range(nodes):
            valx[idx] = powerlaw.rvs(eta, loc=start, scale=mod, size=1,
                                     random_state=None)

    return valx


def evaluate_lognorm(mu_d, perm_scale, minimal, maximal, rng):
    """Define value from limited log normal distribution.

    Parameters
    ----------
    mu_d = (float) mean of log-normal distribution
    perm_scale = (float) scale of log-normal distribution
    minimum = (float) minimum vertical permeability in μD
    maximum = (float) maximum vertical permeability in μD
    rng = (generator) random number generator

    Returns
    -------
    new_value = (float) length value

    Notes
    -----
    1. Values are censored to be within the defined min. and max.
    """
    # Use log-normal distribution to get values + censor values.
    while True:
        new_value = rng.lognormal(mean=mu_d, sigma=perm_scale,
                                  size=None)
        if minimal <= new_value <= maximal:
            break

    return new_value


#
# -----------------------------------------------------------------------------
# - End of module
