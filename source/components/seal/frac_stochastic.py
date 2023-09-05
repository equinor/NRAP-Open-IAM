#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions to create stochastic fracture characteristics and plots.

Author: Ernest N. Lindner
Date: 08/23/2022

Module Name
    frac_stochastic

Contents (20)
    compute_grid_limits(area, frac_controls))
    compute_margin(frac_controls, deltas)
    create_grid_box(deltas, x_coord, y_coord)
    create_scan_box(deltas, x_coord, y_coord, margin):
    compute_tally(areal_density, scan_area)
    determine_line_length(frac_controls, rng)
    create_line(frac_controls, box, trend)
    provide_line_list(frac_controls, num_lines, box)
    clip_lines(line_list, grid)
    identify_inclusive_lines(frac_controls, area, x_coord, y_coord, indx)
    ----
    debug_plot_lines(line_list, separate, indx)
    debug_plot_points(line_list, grid_box, indx)
    define_apertures(frac_controls, frac_lines)
    update_aperture(new_pressure, current_aperture, frac_controls)
    compute_frac_transmissivity(frac_controls, aperture, frac_length):
    compute_threshold(frac_controls, frac_perm)
    create_fractures(frac_controls, rand_apertures, good_lines)
    compute_ave_terms(frac_controls, rand_fractures, area)
    assign_zero_terms()
    compute_random_permeability(frac_controls, grid, alive)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import math                         # For sqrt, pow
import numpy as np                  # For arrays
import matplotlib.pyplot as plt     # For graphics
from matplotlib import collections as mcol   # For multiple-line plots

import seal_config as scfg           # IO directory and file names
import frac_model as fmod           # Class definitions
import frac_random as fran          # Random value generators
import seal_units as sun            # For unit conversion
import seal_file as sfile           # Error in file operations

# Fonts for plotting
LEGEND_FONT = 16                    # Font size for figure title (pt)
TYP_FONT = 12                       # Font size for axes numbers (pt)
REG_FONT = {'fontname': 'Arial',
            'color':    'red',
            'size':      TYP_FONT,
            }                       # X-axis, Y-axis text
TITLE_FONT = {'fontname': 'Arial',
              'color':    'red',
              'size':      LEGEND_FONT,
              'fontweight': 'bold'
              }                     # Figure title text

# Computation Constants
EXTRA = 0.10                        # Percentage of side used as margin

# String lead-in constants
INTRO = "  >>> "                    # Start - printing
TXTIN = "    * "                    # Start for informing on plot

# Debugging options
ECHO_1 = False                      # Provide fracture information
BOZO_1 = False                      # Show fracture lines for grid unit
BOZO_2 = False                      # Show fracture points for grid unit
BOZO_3 = False                      # Save figures
NOTA = False                        # Tell user if no random fractures in cell


def compute_grid_limits(area, frac_controls):
    """Compute grid dimensions for stochastic analysis.

    Parameters
    ----------
    area = (float) area of cell

    Returns
    -------
    side_list = (list)
        [0] = width of a section / grid box
        [1] = length of a section / grid box

    Notes
    -----
    1. Freeform area is arbitrary - so assume a square shape for analysis
    """
    # Use grid values if grid_approach is used.
    if frac_controls['grid_approach']:
        sect_length = frac_controls['cell_height']
        sect_width = frac_controls['cell_width']

    # Else define grid dimensions based on square assumption.
    else:
        sect_width = math.sqrt(area)
        sect_length = sect_width

    # Values for computations.
    side_list = (sect_width, sect_length)

    return side_list


def compute_margin(frac_controls, deltas):
    """Compute margin for stochastic fracture generation.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters
    deltas = (list)
        [0] = width of a grid element
        [1] = length of a grid element

    Returns
    -------
    margin = (float) margin around grid box to include in generation

    Notes
    -----
    1. Assumes a margin is 10% (EXTRA) of grid width/length <or>
       maximum fracture length, whichever is larger
    2. Margin is assumed to be equal on all side of grid box
    """
    # Define limits for selection.
    side = max(deltas[0], deltas[1]) * EXTRA
    extent = frac_controls['length_min']

    # Define margin around box.
    margin = max(side, extent)

    return margin


def create_grid_box(deltas, x_coord, y_coord):
    """Create a grid element from basic dimensions.

    Parameters
    ----------
    deltas => (list) basic x-dimension/ y-dimension  of box
    x_coord, y_coord = (float) coordinates of center of area

    Returns
    -------
    new_box = (rectangle class) rectangular box of one grid sector

    Notes
    -----
    1. No reference grid; use center coordinates for basis.
    """
    # Define left corner based on center coordinates.
    grid_x = x_coord - (deltas[0] / 2.0)
    grid_y = y_coord - (deltas[1] / 2.0)

    # Define point and new box.
    grid_point = fmod.Pointe(grid_x, grid_y)
    new_box = fmod.Rectangle(grid_point, deltas[0], deltas[1])

    return new_box


def create_scan_box(deltas, x_coord, y_coord, margin):
    """Create a region for generation of fractures.

    Parameters
    ----------
    deltas => (list of floats) basic x-dimension / y-dimension of box
    x_coord, y_coord = (float) coordinates of center of area
    margin = (float) extra space around box for line creation / scanning

    Returns
    -------
    new_box = (rectangle class) rectangular scan box for creating fractures
    scan_area = (float) scan area

    Notes
    -----
    1. Use Baecher approach to use a margin around grid box to
        mitigate edge effects.
    """
    # For scan box, define bounds with margin on both sides.
    x_width = deltas[0] + (2.0 * margin)
    y_length = deltas[1] + (2.0 * margin)

    # Compute new corner coordinates - offset with margin.
    new_x_coord = x_coord - (x_width / 2.0)
    new_y_coord = y_coord - (y_length / 2.0)

    # Create scan box - rectangle class.
    scan_point = fmod.Pointe(new_x_coord, new_y_coord)
    new_box = fmod.Rectangle(scan_point, x_width, y_length)

    # Compute area assuming a rectangular area.
    scan_area = x_width * y_length

    return new_box, scan_area


def compute_tally(areal_density, scan_area):
    """Compute the number of fracture centers to evaluate.

    Parameters
    ----------
    areal_density = (float) defined areal density of line centers
    scan area = (float) rectangle area

    Returns
    -------
    result = (int) number of lines to create
    """
    # Compute number of lines to be created.
    result = int(scan_area * areal_density)

    return result


def determine_line_length(frac_controls, rng):
    """Determine line length using log-normal or power law functions.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters
    rng = (generator) random number generator

    Returns
    -------
    line_length = (float) fracture length (m)
    """
    # Select method for generation - lognormal or power.
    if frac_controls['length_approach'] == "LOGNORM":
        line_length = fran.evaluate_lognorm(frac_controls['length_mu'],
                                            frac_controls['length_scale'],
                                            frac_controls['length_min'],
                                            frac_controls['length_max'],
                                            rng)
    else:
        line_length = fran.evaluate_power_length(frac_controls['length_eta'],
                                                 frac_controls['length_min'],
                                                 frac_controls['length_max'],
                                                 rng)
    return line_length


def create_line(frac_controls, box, trend):
    """Create a single line starting inside the scan box.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters
    box = scan box => (rectangle class) region for creation of fractures
    trend = (float) orientation of line

    Returns
    -------
    result_line = (line class) new line to evaluate
    """
    # Define random x,y coordinate of a point within scan box.
    center_x, center_y = fran.evaluate_location(box.pt_1.x_coord,
                                                box.pt_2.x_coord,
                                                box.pt_2.y_coord,
                                                box.pt_3.y_coord,
                                                frac_controls['rng'])

    # Define length of new line.
    line_length = determine_line_length(frac_controls, frac_controls['rng'])

    # Define line orientation, theta(in radians).
    theta = fran.evaluate_orientation(trend,
                                      frac_controls['orient_disperse'],
                                      frac_controls['rng'])

    # Create line from above values.
    result_line = fmod.line_from_center(line_length, theta,
                                        center_x, center_y)

    return result_line


def provide_line_list(frac_controls, num_lines, box):
    """Create a list of lines, all inside scan box.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters
    num_lines = (float) number of lines / corners of line
    box = (rectangle class) scan box for line creation

    Returns
    -------
    line_list = list of lines to evaluate

    Notes
    -----
    1. Assume two orthogonal sets of joints; maintain density as input.
    """
    # Initialize list.
    line_list = []
    line_list.clear()   # necessary???

    # create 'n' lines within box - assume 2 orthogonal joint sets.
    set_1 = int(num_lines / 2)
    for _ in range(set_1):
        new_line = create_line(frac_controls, box,
                               frac_controls['orient_mu'])
        line_list.append(new_line)

    set_2 = int(num_lines / 2)
    for _ in range(set_2):
        new_line = create_line(frac_controls, box,
                               frac_controls['orient_set_2'])
        line_list.append(new_line)

    return line_list


def clip_lines(line_list, grid_box):
    """Clip lines to be within grid box.

    Parameters
    ----------
    line_list = (list) list of lines
    grid_box = (int) grid box to consider

    Returns
    -------
    frac_list = (list) list of lines to convert to fractures
    """
    # Initialize list.
    frac_list = []
    frac_list.clear()   # necessary???

    # create 'n' lines.
    for element in line_list:

        # Clip Line.
        status, try_line = fmod.clipping_2(grid_box, element)

        # If clipping successful, add to new list.
        if status:
            frac_list.append(try_line)

    # end loop.

    return frac_list


def debug_plot_lines(line_list, separate, indx):
    """Plot the lines within a grid box number "indx".

    Parameters
    ----------
    line_list = (list) list of Lines (see class)
    separate = (rectangle class) box around grid being processed
    indx = (int) area/cell number

    Returns
    -------
    None (picture)

    Notes
    -----
    1. Method uses a collection for faster plot.
    """
    # Define window size for plot (inches).
    plt.figure(figsize=(10, 8), dpi=90, num='Seal_Flux line plot')

    # Construct a list of plot lines as required for matplotlib.
    frac_list = []
    for element in line_list:
        iotm = element.delist()
        frac_list.append(iotm)

    # Plot collection of lines.
    tints = np.array([(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)])
    line_collect = mcol.LineCollection(frac_list, colors=tints, linewidths=1)
    plt.gca().add_collection(line_collect)

    # Scale - necessary!
    # x_min, x_max, y_min, y_max = separate.limits().
    plt.xlim(separate.xlimits())
    plt.ylim(separate.ylimits())

    # Add axes and figure titles.
    plt.margins(0.1)
    head = ("\n Line Plot for Area No. " + str(indx + 1) + "\n")
    plt.title(head, fontdict=TITLE_FONT)
    plt.ylabel('Y Position (m)', fontdict=REG_FONT)
    plt.xlabel('X Position (m)', fontdict=REG_FONT)

    # Add line count at bottom left.
    head = ("\n Number of Lines = " + str(len(line_list)))
    plt.text((separate.pt_1.x_coord + 7.0), (separate.pt_1.y_coord + 7.0),
             head, ha='left', va='top', fontsize=10, color='magenta')

    # Save figure.
    if BOZO_3:
        # print cell number, save copy and define file name.
        print(TXTIN + f'Saving Plot #{(indx + 1)}')
        output_file = "fracture_figure-" + str(indx + 1) + scfg.EXTENSION_PNG
        _, destination = sfile.get_path_name(scfg.OUTPUT_DIRECTORY,
                                             output_file, scfg.EXTENSION_PNG)
        # Save figure to file.
        plt.savefig(destination, dpi=300, bbox_inches='tight')

    # Show figure on console.
    plt.show()

    # return None


def debug_plot_points(line_list, grid_box, indx):
    """Plot points of grid lines within a grid box.

    Parameters
    ----------
    line_list = (list) list of Lines (see class)
    grid_box = (rectangle class) box around grid being processed
    indx = (int) area/cell number

    Returns
    -------
    None (picture)
    """
    # Define window size for plot (inches).
    plt.plot(figsize=(10, 8), dpi=90, num='Seal_Flux point plot')

    # Construct a list of data points as required for matplotlib.
    data_x = []
    data_y = []
    for element in line_list:
        x_val, y_val = element.get_mid_point()
        data_x.append(x_val)
        data_y.append(y_val)

    # Plot points using scatter option.
    tons = np.arange(len(data_x))
    plt.scatter(data_x, data_y, c=tons, cmap='jet', s=10)

    # Add axes and figure titles.
    plt.margins(0.1)
    head = ("Point Plot for Block No. " + str(indx + 1) + "\n")
    plt.title(head, fontdict=TITLE_FONT)
    plt.ylabel('Y Position (m)', fontdict=REG_FONT)
    plt.xlabel('X Position (m)', fontdict=REG_FONT)

    # Add point count at bottom left.
    head = ("\n Number of Points = " + str(len(line_list)))
    loc_x = grid_box.pt_1.x_coord
    loc_y = grid_box.pt_1.y_coord
    plt.text(loc_x, loc_y, head, ha='left', va='top', fontsize=10,
             color='magenta')

    # Save figure, if desired.
    if BOZO_3:
        print(TXTIN + f'Saving Point Plot #{(indx + 1)}')
        output_file = "point_figure-" + str(indx + 1) + scfg.EXTENSION_PNG
        _, destination = sfile.get_path_name(scfg.OUTPUT_DIRECTORY,
                                             output_file, scfg.EXTENSION_PNG)
        plt.savefig(destination, dpi=300, bbox_inches='tight')

    # Show figure on console.
    plt.show()

    # return None


# noinspection PyTypeChecker
def identify_inclusive_lines(frac_controls, area, x_coord, y_coord, indx):
    """Determine fractures from lines within a single grid box.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters
    area = (float) area of grid_box under analyses
    x_coord = (float) x-coordinate of center
    y_coord = (float) x-coordinate of center
    indx = (int) block number

    Returns
    -------
    fracture_lines = list of fractures in box

    """
    # Define dimensions of a single area.
    deltas = compute_grid_limits(area, frac_controls)
    extra = compute_margin(frac_controls, deltas)

    # Get density of fractures in box.
    center_density = fran.evaluate_density(frac_controls['density_ave'],
                                           frac_controls['density_min'],
                                           frac_controls['density_max'],
                                           frac_controls['rng'])

    # Create rectangles for grid and scan zones.
    grid_box = create_grid_box(deltas, x_coord, y_coord)
    scan_box, scan_area = create_scan_box(deltas, x_coord, y_coord, extra)

    # Determine number of lines to generate.
    center_tally = compute_tally(center_density, scan_area)

    # Generate lines within scan box. -> list of lines.
    line_list = provide_line_list(frac_controls, center_tally, scan_box)

    # Clip line list to <grid> box to obtain lines for fractures.
    frac_lines = clip_lines(line_list, grid_box)

# -----------------------------------------------------------------------------
    # DEBUG checks.

    # Echo data to console.
    if ECHO_1:
        print(INTRO + f"Grid Number = {indx}")
        print(INTRO + f'The Number of Basic Lines in Block = {center_tally}')
        print(INTRO + f'Grid Block Area = {area}')
        print(INTRO + f'The Density = {center_density}')
        print(INTRO + f'Number of Fracture Lines = {len(frac_lines)}')

    # Show/save fracture line plot.
    if BOZO_1 or frac_controls['plot_fractures']:
        debug_plot_lines(frac_lines, scan_box, indx)

    # Plot center points.
    if BOZO_2:
        # noinspection PyTypeChecker
        debug_plot_points(frac_lines, grid_box, indx)

    # Print line if a plot was produced or text echo.
    if BOZO_1 or BOZO_2 or ECHO_1:
        print()
    # -------------------------------------------------------------------------

    return frac_lines


def define_apertures(frac_controls, frac_lines):
    """Determine aperture of each fracture line.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters
    frac_lines = (list) list of fractures

    Returns
    -------
    aperture_list = (list) list of apertures

    Notes
    -----
    1. Apertures are not correlated with line length.
    """
    # Get number of lines.
    aperture_list = []
    mu_ap = frac_controls['aperture_mu']
    scale = frac_controls['aperture_scale']
    min_ap = frac_controls['aperture_min']
    max_ap = frac_controls['aperture_max']
    alpha = frac_controls['aperture_alpha']
    beta = frac_controls['aperture_beta']

    # Create aperture list as user wants.
    for element in frac_lines:
        if frac_controls['correlate_approach']:
            extent = element.length()
            new_aperture = fran.correlate_aperture(alpha, beta, extent,
                                                   frac_controls['rng'])
        else:
            new_aperture = fran.evaluate_aperture(mu_ap, scale,
                                                  min_ap, max_ap,
                                                  frac_controls['rng'])
        # Add to list.
        aperture_list.append(new_aperture)

    return aperture_list


def update_aperture(new_pressure, current_aperture, frac_controls):
    """Re-evaluate fracture aperture based on expected pressure change.

    Parameters
    ----------
    new_pressure = (float) pressure at fracture (Pa)
    current_aperture = (float) current fracture aperture (mm)
    frac_controls = (dict) fracture control parameters

    Returns
    -------
    corrected = (float) adjusted aperture (mm)

    Notes
    -----
    1. Typical hydrostatic gradient = 10.45 KPa/m
    2. Typical lithostatic gradient = 24.0 KPa/m (defined in seal_units)
    3. Assume current aperture was identified at original static pressure.
    4. Equation w factor is to compute new aperture.
    5. Watch: - stresses in routine are in MPA!!!
    """
    # Correct aperture for pressure, if desired.
    if frac_controls['pressure_approach']:
        # Constants for calc
        brine_density = 1007.0  # kg/m^3
        new_pressure *= sun.pa_to_mpa()  # now in MPa!

        # Compute original hydrostatic/lithostatic pressure at seal base.
        pressure_change = (frac_controls['ave_base_depth']
                           - frac_controls['static_depth'])
        pressure_change *= brine_density * sun.gravity() * sun.pa_to_mpa()

        hydro_pressure = frac_controls['static_pressure'] * sun.pa_to_mpa()
        hydro_pressure += pressure_change

        litho_stress = sun.earth_pressure(frac_controls['ave_base_depth'])
        litho_stress *= sun.pa_to_mpa()
        effect_stress = litho_stress - hydro_pressure

        # Compute aperture at static pressure w. orig. conditions.
        limit_stress = frac_controls['stress_limit'] * sun.pa_to_mpa()
        beta = ((limit_stress + frac_controls['theta_aperture'])
                / limit_stress)
        stress_term = (effect_stress
                       / (effect_stress + frac_controls['theta_aperture']))
        ap_term = (frac_controls['wide_aperture']
                   - frac_controls['residual_aperture'])
        aperture_k = (1.0 - beta * stress_term) * ap_term
        aperture_k += frac_controls['residual_aperture']

        # Compute correction term.
        factor = current_aperture / aperture_k

        # Determine new aperture at new pressure.
        effect_stress = litho_stress - new_pressure
        if effect_stress <= 0.0:
            # If negative stress, use maximum aperture value.
            corrected = frac_controls['wide_aperture'] * factor
        elif effect_stress > frac_controls['stress_limit'] * sun.pa_to_mpa():
            # If larger stress, use residual aperture value.
            corrected = frac_controls['res_aperture'] * factor
        else:
            # Recompute aperture at current (new) stress.
            stress_term = (effect_stress
                           / (effect_stress + frac_controls['theta_aperture']))
            aperture_k = ((1.0 - beta * stress_term) * ap_term
                          + frac_controls['residual_aperture'])
            corrected = factor * aperture_k
    else:
        # No correction to aperture.
        corrected = current_aperture

    return corrected


def compute_frac_transmissivity(frac_controls, aperture, frac_length):
    """Compute the transmissivity of a fracture (in m^4).

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters
    aperture = (float) fracture aperture
    frac_length = (float) length of fracture

    Returns
    -------
    trans = (float) transmissivity (m^4)

    Notes
    -----
    1. transmissivity = <kf>*<Af>, where:
        Af = fracture area = length x aperture
        kf = (aperture^2) / 12
    2. transmissivity = length * (aperture^3)/12
    3. Final units need to be => m^4.
    """
    # Convert fracture aperture from millimeters to meters.
    aperture_meter = aperture * sun.mm_to_m()

    # Define connectivity, stress, roughness, tortuosity factor.
    geo_factor = frac_controls['connect_factor']

    # Define transmissivity (m^4).
    trans = frac_length * math.pow(aperture_meter, 3.0) / 12.0
    trans *= geo_factor

    return trans


def compute_threshold(threshold_factor, aperture):
    """Compute the threshold pressure of a fracture.

    Parameters
    ----------
    threshold_factor = (float) factor to convert aperture to threshold pressure
    aperture = (float) fracture aperture (mm)

    Returns
    -------
    result = (float) threshold pressure for fracture (Pa/mm)
    """
    result = threshold_factor / aperture

    return result


def create_fractures(frac_controls, rand_apertures, good_lines):
    """Create list of fractures for grid box for each good line.

    Parameters
    ----------
    frac_controls = (dict) fracture control parameters
    rand_apertures = (list) list of apertures
    good_lines = (list) list of lines within grid box

    Returns
    -------
    fracture_list = (list) list of fractures

    Notes
    -----
    1. Apertures not correlated with line length.
    """
    # Create a list of fractures.
    fracture_list = []

    # Define each fracture and store in new list.
    for indx, element in enumerate(good_lines):
        # Create a fracture.
        start, end = element.get_points()
        new_fracture = fmod.Fracture(start, end)

        # Define aperture for fracture - correct aperture for pressure.
        frac_width = update_aperture(frac_controls['ave_base_pressure'],
                                     rand_apertures[indx], frac_controls)
        new_fracture.aperture = frac_width

        # Define strike for fracture.
        new_fracture.strike = element.trend()

        # Define transmissivity for fracture (in terms of m^4).
        new_fracture.trans = \
            compute_frac_transmissivity(frac_controls,
                                        rand_apertures[indx],
                                        element.length())

        # Define Threshold pressure based on aperture.
        new_fracture.threshold = \
            compute_threshold(frac_controls['threshold_factor'],
                              rand_apertures[indx])

        # Store new fracture.
        fracture_list.append(new_fracture)

    return fracture_list


def compute_ave_terms(frac_controls, fracture_list, rect_area):
    """Compute average terms for a grid element.

    Parameters
    ----------
    frac_controls = (dict) fracture controls
    fracture_list = (list) list of fractures in grid element
    rect_area = (float) total area of area under consideration

    Returns
    -------
    parameter_list = average parameters for block:
        0 = average element.aperture (mm)
        1 = average element.extent (m)
        2 = average element.strike (degrees from N)
        3 = average element.trans (m^4)
        4 = average element.threshold (Pa)
        5 = Number of fractures within grid element
        6 = total permeability of grid element (microdarcys!)
    """
    # Zero values to store from list.
    aperture = 0.0
    length = 0.0
    trend = 0.0
    transmissivity = 0.0
    cap_pressure = 0.0
    vals = 0

    # enter loop only if fracture exists.
    if fracture_list:
        # Compute sums and averages of parameters.
        for element in fracture_list:
            aperture += element.aperture
            length += element.extent()
            trend += element.strike
            transmissivity += element.trans
            cap_pressure += element.threshold

        # Compute averages of selected terms.
        vals = len(fracture_list)
        aperture /= vals
        length /= vals
        trend /= vals
        transmissivity /= vals
        cap_pressure /= vals

        # ----------------------------------------
        # Compute equivalent element permeability.
        sum_perm = 0.0
        for element in fracture_list:
            frac_wide = element.aperture * sun.mm_to_m()
            perm = math.pow(frac_wide, 2.0) / 12.0
            frac_area = frac_wide * element.extent()
            sum_perm += perm * frac_area
        sum_perm /= rect_area

        # Convert permeability to mD, and add correction.
        sum_perm *= (sun.metersq_to_microd() *
                     frac_controls.get('connect_factor'))
    else:
        # No fracture data.
        sum_perm = 0.0

    # Create a list of average for block.
    parameter_list = [aperture, length, trend, transmissivity,
                      cap_pressure, vals, sum_perm]

    return parameter_list


def assign_zero_terms():
    """Assign zero terms for a grid element.

    Parameters
    ----------
    N/A

    Returns
    -------
    permeability = total fracture permeability for grid box
    ave_thresh = average threshold value for grid box
    parameter_list = average parameters for block:
        0 = average element.aperture (mm)
        1 = average element.extent (m)
        2 = average element.strike (degrees from N)
        3 = average element.trans (m^4)
        4 = average element.threshold (Pa)
        5 = Number of fractures within grid element
        6 = total permeability of grid element (microdarcys!)
    """
    # Zero values to store from list.
    aperture = 0.0
    length = 0.0
    trend = 0.0
    transmissivity = 0.0
    vals = 0.0
    total_permeability = 0.0
    cap_pressure = 0.0

    # Create a list of average for block.
    parameter_list = [aperture, length, trend, transmissivity,
                      cap_pressure, vals, total_permeability]

    return parameter_list


def compute_random_permeability(frac_controls, grid, alive):
    """Determine equivalent fracture permeability for random fractures.

    Parameters
    ----------
    frac_controls = (dict) fracture controls
    seal_controls = (dict) seal controls
    grid = (list) collection of seal cells
    alive = (bool) = stand-alone operation

    Returns
    -------
    equi_perm = (array) equivalent permeabilities (mD)
    equi_threshold = (array) equivalent threshold pressures
    block_params = (list) list of averages for each block
    """
    # list of parameters for grid blocks.
    numbr = frac_controls['num_cells']
    block_params = []
    equi_perm = np.zeros(numbr)
    equi_threshold = np.zeros(numbr)

    # Tell user about saving figures.
    if frac_controls['plot_fractures'] and alive:
        if BOZO_3:
            sfile.echo_status(alive, "SAVING FIGURES; PLEASE WAIT.")
        else:
            sfile.echo_status(alive, "PLOTTING FIGURES; PLEASE WAIT.")
        print("        " + "If Necessary - "
              + "Hit Exit at Top-Right of Plot Window to Proceed ...")

    # Loop over each grid block.
    for indx, cell in enumerate(grid):
        # Find lines that are within a block.
        area_bloc = cell.area
        good_lines = identify_inclusive_lines(frac_controls, area_bloc,
                                              cell.x_center, cell.y_center,
                                              indx)

        if len(good_lines) > 0:
            # Define apertures of lines in block.
            rand_apertures = define_apertures(frac_controls, good_lines)

            # Create fractures with orientation, aperture and permeability.
            rand_fractures = create_fractures(frac_controls, rand_apertures,
                                              good_lines)

            # Compute equivalent permeability / threshold for grid element.
            instance_values = \
                compute_ave_terms(frac_controls, rand_fractures, area_bloc)
            equi_perm[indx] = instance_values[6]
            equi_threshold[indx] = instance_values[4]

            # Store parameters for each block for summary file.
            block_params.append(instance_values)

        else:
            # No Fractures! (Density may be too low; warn user.)
            if alive and NOTA:
                reason = f'No Fractures Created in Cell #{indx}!'
                sfile.issue_warning(reason)

            # Assign zero to equivalent permeability / threshold.
            equi_perm[indx] = frac_controls['rock_perm_ave']
            equi_threshold[indx] = frac_controls['ref_matrix_perm']
            instance_values = assign_zero_terms()

            # Store parameters for each block for summary file.
            block_params.append(instance_values)

    # Stop plot of fractures a second time.
    frac_controls['plot_fractures'] = False

    return equi_perm, equi_threshold, block_params


# -----------------------------------------------------------------------------
# - End of module
# -------1---------2---------3---------4---------5---------6---------7---------8
