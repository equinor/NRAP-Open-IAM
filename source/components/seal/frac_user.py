#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions to process user fracture data & compute matrix permeability.

Author: Ernest N. Lindner
Date: 07/27/2022

Module Name
    frac_user

Contents (6)
    select_fractures(frac_controls, line_list, grid, parameters)
    define_user_lines(user_data)
    process_user_lines(frac_controls, grid, user_data)
    compute_matrix_threshold(frac_controls, permeability)
    define_matrix_perm(frac_controls)
    sum_threshold_arrays(frac_controls, rand_array, user_array, matrix_array)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import math
import numpy as np                  # Array operations

import frac_model as fmod           # Class definitions
import frac_random as fran          # Random value generators
import frac_stochastic as fsch      # Basic functions for stochastic fractures

# Constants
P_MATIX = 0.25                      # Matrix permeability exponent
PARA = "    > "                     # Line start


def select_fractures(frac_controls, line_list, grid_box, parameters):
    """Create fractures within a specific grid box.

    Parameters
    ----------
    frac_controls = (dict) dictionary of fracture controls
    line_list = (list) list of lines
    grid_box = (rectangle class) grid box to consider
    parameters = (array) NumPy array of parameters for line being examined

    Returns
    -------
    frac_list = list of fractures
    """
    # Initialize list.
    frac_list = []

    # create 'n' lines.
    for indx, element in enumerate(line_list):

        # Clip Line.
        status, try_line = fmod.clipping_2(grid_box, element)

        # If clipping successful, add fracture to new list.
        if status:
            # Create a new fracture from clipped line.
            start, end = try_line.get_points()
            good_fracture = fmod.Fracture(start, end)

            # Compute strike from line data.
            good_fracture.strike = try_line.trend()

            # Define aperture from parameter input for line.
            aperture = parameters[indx, 4]

            # Compute permeability (in m^4) and transmissivity from data.
            good_fracture.trans = \
                fsch.compute_frac_transmissivity(frac_controls,
                                                 aperture,
                                                 try_line.length())

            # Compute threshold pressure and set aperture parameter.
            good_fracture.threshold = \
                fsch.compute_threshold(frac_controls['ref_user_aperture'],
                                       aperture)
            good_fracture.aperture = aperture

            # Store line in list.
            frac_list.append(good_fracture)

    # end loop

    return frac_list


def define_user_lines(user_data):
    """Construct data lines from user data.

    Parameters
    ----------
    user_data = (list fo line class) list of input of user fracture parameters

    Returns
    -------
    user_lines = (list) list of user lines
    user_array = (array) NumPy array of input data
    counter = (list) lines from file as read

    Notes
    -----
    First four columns of data have coordinates
    """
    # Convert list to NumPy array.
    user_array = np.asarray(user_data, dtype=np.float64)

    # Define user lines from input data array.
    counter = 0
    user_lines = []
    for indx in range(0, len(user_data)):

        # Define values for clarity; first four values are coordinates.
        xcoord_1 = user_array[indx, 0]
        ycoord_1 = user_array[indx, 1]
        xcoord_2 = user_array[indx, 2]
        ycoord_2 = user_array[indx, 3]

        # Construct points and line.
        point_1 = fmod.Pointe(xcoord_1, ycoord_1)
        point_2 = fmod.Pointe(xcoord_2, ycoord_2)
        new_line = fmod.Segment(point_1, point_2)

        user_lines.append(new_line)
        counter += 1

    return user_lines, user_array, counter


def process_user_lines(frac_controls, grid, user_data):
    """Determine user fracture permeability for each grid area.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters
    user_data = (list) list of user fracture parameters

    Returns
    -------
    bloc_perm = (array) permeability of each grid area
    bloc_threshold = (array) threshold of each grid area
    bloc_averages = (list) list of values of each grid area
    """
    # Define user lines for analysis.
    user_lines, user_parameters, counter = define_user_lines(user_data)

    # Set controls and output array.
    frac_controls['user_fractures'] = counter
    bloc_perm = np.zeros(frac_controls['num_cells'])
    bloc_thresh = np.zeros(frac_controls['num_cells'])
    bloc_averages = []

    # Iterate over each grid box and evaluate permeability - Reversed Order!.
    for indx, cell in enumerate(grid):
        # Define grid values.
        deltas = fsch.compute_grid_limits(cell.area, frac_controls)

        # Define grid box.
        grid_box = fsch.create_grid_box(deltas, cell.x_center, cell.y_center)

        # Create fractures from input.
        user_fractures = select_fractures(frac_controls, user_lines,
                                          grid_box, user_parameters)

        # Compute equivalent permeability / threshold for grid element.
        average_values = \
            fsch.compute_ave_terms(frac_controls, user_fractures,
                                   grid_box.area())
        bloc_perm[indx] = average_values[6]
        bloc_thresh[indx] = average_values[4]

        # Save average parameters for each block (for summary file).
        bloc_averages.append(average_values)

    return bloc_perm, bloc_thresh, bloc_averages


def compute_matrix_threshold(frac_controls, permeability):
    """Determine matrix threshold pressure of each grid area.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters
    permeability = (float) ave matrix permeability for area

    Returns
    -------
    threshold = (float) threshold pressure (Pa).

    Notes
    -----
    1. See Appendix E in user manual for logic.
    """
    # Define reference terms for clarity.
    ref_threshold = frac_controls['ref_matrix_threshold']
    ref_permeability = frac_controls['ref_matrix_perm']

    # compute threshold factor based on equation.
    term_1 = ref_permeability / permeability
    term_2 = math.pow((1.00 / term_1), P_MATIX)
    factor = math.sqrt(term_1*term_2)

    # Compute new threshold.
    threshold = factor * ref_threshold

    return threshold


def define_matrix_perm(frac_controls):
    """Determine matrix permeability and threshold of each grid area.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters

    Returns
    -------
    perm_array = (array) NumPy array of permeabilities for grid
    threshold_array = (array) NumPy array of threshold pressures for grid
    matrix_perm = (list) list of matrix values for grid.
    """
    # Define variables from frac_controls.
    numbr = frac_controls['num_cells']
    mu_rock = frac_controls['rock_perm_mu']
    scale_rock = frac_controls['rock_perm_scale']
    max_rock = frac_controls['rock_perm_max']
    min_rock = frac_controls['rock_perm_min']

    # Define lists for computations.
    matrix_threshold = []
    matrix_perm = []
    matrix_list = []

    # Define matrix permeability and threshold pressure for each grid area.
    for _ in range(0, numbr):
        new_permeability = fran.evaluate_matrix_perm(mu_rock, scale_rock,
                                                     min_rock, max_rock,
                                                     frac_controls['rng'])
        new_threshold = compute_matrix_threshold(frac_controls,
                                                 new_permeability)
        matrix_perm.append(new_permeability)
        matrix_threshold.append(new_threshold)
        matrix_list.append([new_permeability, new_threshold])

    # For file output, convert lists to NumPy array and resize.
    threshold_array = np.asarray(matrix_threshold)
    perm_array = np.asarray(matrix_perm)

    return perm_array, threshold_array, matrix_list


def sum_threshold_arrays(frac_controls, rand_array, user_array, matrix_array):
    """Sum threshold arrays of random, user and matrix cases.

    Parameters
    ----------
    frac_controls = (dict) dictionary of fracture controls
    rand_array = (array) random fracture threshold pressures for grid
    user_array = (array) user defined fracture threshold pressures for grid
    matrix_array = (array) matrix threshold pressures for grid

    Returns
    -------
    final_threshold = (array) resulting threshold pressures

    Notes
    -----
    1. Logic: Use of threshold pressures either from matrix <or> random
    fractures for final value:
        a) If fractures exist, use maximum value.
        b) If no fractures exist, use matrix values.
        c) Fracture approach = True, it does not mean that there are
          results for each grid element.
    2. Matrix value is considered default value.
    """
    # Define basic variables.
    numbr = frac_controls['num_cells']
    final_threshold = np.zeros(numbr)

    # Compute resulting pressure - based on logic of code.
    for indx in range(numbr):

        # If both fracture options were used, use random fractures.
        if (frac_controls['random_approach']
                and frac_controls['user_approach']):
            # Check if array values may be zero.
            if rand_array[indx] > 0.0:
                final_threshold[indx] = rand_array[indx]
            elif user_array[indx] > 0.0:
                final_threshold[indx] = user_array[indx]
            else:
                final_threshold[indx] = matrix_array[indx]

        # If only random fractures results exist, (i.e., no user) use.
        elif frac_controls['random_approach']:
            # If random value significant, use it.
            if rand_array[indx] > 0.0:
                final_threshold[indx] = rand_array[indx]
            else:
                final_threshold[indx] = matrix_array[indx]

        # If only user fractures exist, use these results.
        elif frac_controls['user_approach']:
            if user_array[indx] > 0.0:
                final_threshold[indx] = min(user_array[indx],
                                            matrix_array[indx])
            else:
                final_threshold[indx] = matrix_array[indx]

        else:
            # if either random or user fractures exist,
            #   * use matrix threshold pressures.
            final_threshold[indx] = matrix_array[indx]

    return final_threshold


#
# -----------------------------------------------------------------------------
# - End of module
# -------1---------2---------3---------4---------5---------6---------7---------8
