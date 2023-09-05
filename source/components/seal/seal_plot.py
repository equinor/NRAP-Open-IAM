#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions to manage and plot data for seal_flux.

Author: Ernest N. Lindner
Date: 08/18/2022

Module Name
    seal_plot

Contents (14)
    control_time_history(max_realizations, si_units, max_time)
    check_series(start_run, si_units)
    control_time_series(max_realizations, si_units, max_time)
    manage_bar_plot(seal_controls, z_data)
    manage_permeability_plot(x_data, y_data, z_data)
    control_permeability_plot(alive, seal_controls, grid)
    manage_co2_contour(seal_controls, grid)
    control_co2_contour(max_realizations, x_data, y_data)
    plot_contour(fig_legend, bar_title, x_data, y_data, z_data)
    plot_bar_chart(x_vals, y_vals, z_values, z_names)
    ----
    plot_time_history(realization, x_data, y_data, max_time)
    save_figure_to_file()
    plot_time_series(x_data, data_list, data_maximum, start_sim, si_units,
    max_time)
    plot_manager(alive, seal_controls, grid, si_units)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import math                         # For truncation
import copy                         # For array adjustment
import numpy as np                  # For array operations
import matplotlib.pyplot as plt     # For plots
from matplotlib import cm           # For color map

# For older matplotlib, import of Axes 3D is needed for projection=3d keyword.
# from mpl_toolkits.mplot3d import Axes3D
import seal_pluto as aux            # Support functions for plot
import seal_file as sfile           # For error report
import seal_config as scfg           # config variables for I/O
import seal_units as sun             # Unit conversion

# Start characters for statements
MAJ_IN = "\n  --> "                 # Page + general message start
CMD_IN = "      >>> "               # Command message start

# Statements for plotting and error messages
CAUTION = CMD_IN + "Hit Exit at Top-Right of Plot Window to Proceed ..."
WAIT_FOR_IT = CMD_IN + "Please wait for plot ..."

# Plot constants
N_DIVISIONS = 10                    # Subdivisions for contouring
PD_LINEWIDTH = 1                    # line width for contour
RESOLVE = 90                        # Figure display resolution (dpi)
COLOR_BAR = 'jet'                   # Color palette for fill in contouring
C_LINES = 'black'                   # Line color for contours
SHIFT = 20                          # Distance to shift y-label on color bar
BAR_SHIFT_X = 0.75                  # Shift of 3D bars along X axis
BAR_SHIFT_Y = 0.50                  # Shift of 3D bars along Y axis
PLOT_LINE_LIMIT = 25                # Limit for legend in figure
SAVE_FIGURE = True                  # Control to save figure to file

# Size definitions
HEIGHT = 8                          # Figure height (inches)
WIDTH = 12                          # Figure width (inches)
LEGEND_FONT = 14                    # Font size for figure title (pt)
TYP_FONT = 12                       # Font size for axes numbers (pt)

# Font definitions
REG_FONT = {'fontname': 'Arial',    # X-axis, Y-axis text
            'color':    'black',
            'size':      TYP_FONT,
            }
TITLE_FONT = {'fontname': 'Arial',
              'color':    'black',
              'size':      LEGEND_FONT,
              }
#               'fontweight': 'bold'
LEG_FONT = {'fontname': 'Arial',
            'color':    'black',
            'size':      LEGEND_FONT,
            }
BAR_FONT = {'fontname': 'Arial',
            'color':    'black',
            'size':      TYP_FONT,
            }


def control_time_history(max_realizations, si_units, max_time):
    """Provide overall control for plotting a simulation from output.

    Parameters
    ----------
    max_realizations = (int) maximum number of realizations in analysis

    Returns
    -------
    None
    """
    # Loop until finished plotting.
    while True:
        # Query for simulation number.
        response, realization = aux.simulation_query(max_realizations, 0)

        # Plot - if user still desires a plot.
        if response:
            # Get data & plot series. Recast units if necessary.
            leakage_data, time_data = \
                aux.setup_time_history_data(realization, si_units)

            # Check data for range.
            check_up = aux.check_plot_data(leakage_data)

            # Plot data if OK.
            if len(time_data) > 1 and check_up:
                plot_time_history(realization, time_data, leakage_data,
                                  si_units, max_time)
            elif len(time_data) <= 1:
                # Major error.
                reason = "No Time Data for Plot."
                sfile.opx_problem(reason, err='')
            else:
                # No separation.
                print(CMD_IN + "Plot Issue! Data Negative Or Series Has "
                      + "Maximum Equal To Minimum!")
                print(CMD_IN + "Plotting Is Discontinued!")

            # Ask about another plot; exit if not.
            plot_again = aux.create_history_query(1, False)
            if not plot_again:
                break
        else:
            break
    # end while

    # return None


def check_series(sim_valu, si_units):
    """Provide error check of series data.

    Parameters
    ----------
    specific_real= (int) one simulation from analysis
    si_units = (bool) use of metric units

    Returns
    -------
    match_err = (bool) error in time series
    data_err = (bool) error in data
    """
    # Error check time data (use one series to check)!
    # --> Check that the time string is greater than 1.
    leaker, time_match = \
        aux.setup_time_history_data(sim_valu, si_units)

    # --> Error check on length of time series
    if len(time_match) > 1:
        match_err = False
    else:
        match_err = True      # Error found

    # -- > Check data for range.
    err_check = aux.check_plot_data(leaker)
    if err_check:
        data_err = False
    else:
        data_err = True        # Error found

    return match_err, data_err


def control_time_series(max_realizations, si_units, max_time):
    """Provide the overall control to plot a series of simulations from output.

    Parameters
    ----------
    max_realizations = (int) upper limit of simulations from analysis
    si_units = (bool) use of metric units
    max_time = (float) maximum time limit on x-axis

    Returns
    -------
    None
    """
    # Loop until finished plotting.
    while True:
        # Query to get simulation numbers (start and end).
        response, start_run, end_run = aux.range_query(max_realizations)

        # Plot - if user desires a plot.
        if response:

            # Get data and store NumPy arrays in list.
            #   -- adjust counter to catch last run.
            data_list = []
            final_max = 0.0
            time_data = np.array(0.0)

            # Set up for series of plots.
            for realization in range(start_run, (end_run + 1)):
                leakage_data, time_data = \
                    aux.setup_time_history_data(realization, si_units)
                data_list.append(leakage_data)
                instance_max = leakage_data.max()
                if instance_max > final_max:
                    final_max = instance_max

            # Error check data (use one series to check)!
            match_err, data_err = check_series(start_run, si_units)

            # Plot data arrays on same plot diagram if data OK.
            if not match_err and not data_err:
                # Data is AOK.
                plot_time_series(time_data, data_list, final_max,
                                 start_run, si_units, max_time)
            elif match_err:
                # Time string error.
                sfile.opx_problem("No Time Data for Plot.", err='')
            else:
                # Value error.
                print(CMD_IN + "Plot Issue! Data Negative Or Series Has "
                      + "Maximum Equal To Minimum!")
                print(CMD_IN + "Plotting Is Discontinued!")

            # Ask about another plot; exit if not.
            plot_again = aux.create_history_query(1, True)
            if not plot_again:
                break
        else:
            break
        # end while

    # return None


def manage_bar_plot(seal_controls, z_data):
    """Manage setup for permeability bar chart plot.

    Parameters
    ----------
    seal_controls = (dict) seal controls
    z_data = (array of floats) z-axis data (2D)

    Returns
    -------
    None

    Notes
    -----
    1. 3D Bar chart requires grid approach!
    2. No variability check.
    """
    # Get exponent and factor z_data.
    cite, show_units, zed_data = aux.set_exponent(z_data)
    if cite == 0:
        z_title = "Permeability (mD)"
    else:
        z_title = "Permeability (" + show_units + " mD)"

    # Define configuration for plot, 1D array needed for plot.
    lines = seal_controls['grid_rows']
    columns = seal_controls['grid_cols']
    z_values = zed_data.flatten()

    # plot values.
    plot_bar_chart(columns, lines, z_values, z_title)

    # return None


def manage_permeability_plot(x_data, y_data, z_data):
    """Structure data and plot contour of permeability and bar chart.

    Parameters
    ----------
    x_data = (array of floats) x-axis data (2D)
    y_data = (array of floats) y-axis data (2D)
    z_data = (array of floats) z-axis data (2D)

    Returns
    -------
    None

    """
    # Test for variation in z values.
    variable = aux.test_contour(z_data)

    # If data varies, then draw contour plot.
    if variable:
        # Get exponent and factor z_data.
        cite, show_units, zed_data = aux.set_exponent(z_data)

        # Construct titles for plot; add returns to add space.
        fig_legend = "\n " + "Seal Permeability of Last Simulation " + "\n"

        # Add exponent value as needed to axis title.
        if cite == 0:
            bar_title = "Permeability (mD)"
        else:
            bar_title = "Permeability (" + show_units + " mD)"

        # Contour plot of permeabilities.
        plot_contour(fig_legend, bar_title, x_data, y_data, zed_data)

    # return None


def control_permeability_plot(seal_controls, grid):
    """Construct permeability data for plotting & send to plotting routines.

    Parameters
    ----------
    seal_controls = (dict) seal controls
    grid = (list) (class) a collection of cells

    Returns
    -------
    None
    """
    # Setup data for contour plot.
    if seal_controls['grid_approach']:
        # ----------------------------
        # Uniform Grid Approach:
        # Get 2D coordinate data and permeabilities.
        coord_exist = True
        lines = seal_controls['grid_rows']
        columns = seal_controls['grid_cols']
        x_data, y_data = aux.construct_coordinate_data(lines, columns, grid)
        z_data = aux.construct_permeability_data(lines, columns, grid)

    else:
        # ----------------------------
        # Irregular Data Approach:
        # Construct uniform arrays from irregular data;
        #   - data may not have coordinates! check => condition.
        z_line = aux.construct_irregular_permeability(grid)
        x_data, y_data, z_data, coord_exist = \
            aux.construct_irregular_data(grid, z_line)

    # Plotting
    # If coordinates exist, create contour plot of last permeability.
    print(MAJ_IN + "1. CONTOUR PLOT.")
    if coord_exist:
        manage_permeability_plot(x_data, y_data, z_data)
    else:
        # If no coordinates, no contour; print error message.
        print(CMD_IN + "Missing Coordinate Data; Contour Plot Not Possible!"
              + "\n" + CMD_IN + "No Contour Plot Constructed!",
              flush=True)

    # >> If grid exists, plot 3-D bar chart of last permeability.
    print(MAJ_IN + "2. 3D BAR CHART PLOT.")
    if seal_controls['grid_approach']:
        manage_bar_plot(seal_controls, z_data)
    else:
        # If no coordinates, restrict bar chart.
        print(CMD_IN + "Missing a Grid; 3D Chart Not Possible!"
              + "\n" + CMD_IN + "No Bar Chart Constructed!", flush=True)

    # return None


def manage_co2_contour(seal_controls, grid, leakage, simulation):
    """Structure data and plot contour of CO2 release.

    Parameters
    ----------
    seal_controls = (dict) seal controls
    grid = (list) (class) a collection of cells
    leakage = (array) CO2 flux through seal (1D)
    simulation = specific simulation

    Returns
    -------
    coord_exist = (bool) on existence of coordinate data
        True = data exist
        False = no coordinates
    """
    # Setup data for contour plot.
    if seal_controls['grid_approach']:
        # ----------------------------
        # Uniform Grid Approach.

        # Get 2D coordinate data and permeabilities.
        coord_exist = True

        z_data = leakage.reshape(seal_controls['grid_rows'],
                                 seal_controls['grid_cols'])
        x_data, y_data = \
            aux.construct_coordinate_data(seal_controls['grid_rows'],
                                          seal_controls['grid_cols'], grid)

    else:
        # ----------------------------
        # Irregular Data Approach.

        # Construct uniform arrays from irregular data;
        #   - data may not have coordinates! check for this condition.
        x_data, y_data, z_data, coord_exist = \
            aux.construct_irregular_data(grid, leakage)

    # Check if coordinates exist - if OK, then plot.
    if coord_exist:
        # Get exponent for bar and then factor z_data.
        cite, show_units, zed_data = aux.set_exponent(z_data)

        # Construct titles and plot contours; add space for larger margin.
        fig_legend = ("\n Seal Seepage for Realization #" + str(simulation)
                      + "\n")

        # Define units and magnitude for figure.
        if seal_controls['use_si_units']:
            cubit = ' tonnes'
        else:
            cubit = ' tons'

        if cite == 0:
            bar_title = "Seepage (" + cubit + ")"
        else:
            bar_title = "Seepage (" + show_units + cubit + ")"

        # Plot contour of CO2 flux.
        plot_contour(fig_legend, bar_title, x_data, y_data, zed_data)

    return coord_exist


def control_co2_contour(max_realizations, seal_controls, grid):
    """Provide overall control to plot a simulation from output.

    Parameters
    ----------
    max_realizations = (int) upper limit of simulations from analysis
    seal_controls = (dict) seal controls
    grid = (list) (class) a collection of cells

    Returns
    -------
    None
    """
    # Loop until finished plotting.
    while True:
        # Query for simulation number.
        response, simulation = aux.simulation_query(max_realizations, 0)

        # Plot - if user still desires a plot.
        if response:
            # Get z data.
            leakage_data = aux.setup_co2_seepage_data(simulation)

            # Recast leakage units into English/US units, if desired.
            if not seal_controls['use_si_units']:
                leakage_data = sun.tonne_to_ton() * leakage_data

            # Create contour plot - first check if data is variable!
            variable = aux.test_contour(leakage_data)
            if variable:
                # Plot Contour.
                coord_exist = manage_co2_contour(seal_controls, grid,
                                                 leakage_data, simulation)

                # Stop process if no coordinate data.
                if not coord_exist:
                    break

            # Ask about another plot; exit if not.
            plot_again = aux.create_contour_query(1)
            if not plot_again:
                break
        else:
            break
    # end while

    # return None


def plot_contour(fig_legend, bar_title, x_data, y_data, zed_data):
    """Plot contour data using matplotlib functions.

    Parameters
    ----------
    fig_legend = (str) legend of plot
    bar_title = (str) title of scale bar
    x_data, y_data = (arrays) coordinate values
    z_data = (array) - data to plot (Numpy)

    Returns
    -------
    None
    """
    # Notify user to wait.
    print(WAIT_FOR_IT)

    # Define window size for plot (inches).
    plt.figure(figsize=(WIDTH, HEIGHT), dpi=RESOLVE,
               num='Seal_Flux Contour Plot')

    # Draw legend and axis labels - X/Y = areal basis; add space at top.
    plt.title(fig_legend, fontdict=TITLE_FONT, y=1.02)
    plt.xlabel('X Axis (m)', fontdict=REG_FONT)
    plt.ylabel('Y Axis (m)', fontdict=REG_FONT)

    # Plot contours for both line values and fill controls.
    plt.contour(x_data, y_data, zed_data, N_DIVISIONS,
                colors=C_LINES, linewidths=PD_LINEWIDTH)
    contour_x = plt.contourf(x_data, y_data, zed_data, N_DIVISIONS,
                             cmap=COLOR_BAR)

    # plt.clabel(contour_x, inline=True, fontsize=8)  # plot line numbers.
    clb = plt.colorbar(contour_x, shrink=0.8, orientation='vertical')
    clb.ax.set_ylabel(bar_title, BAR_FONT, rotation=270,
                      labelpad=SHIFT)

    # Output figure to file or screen.
    if SAVE_FIGURE:
        save_figure_to_file()
    else:
        print(CAUTION)
        plt.show()

    # Clear plot to prevent over-writes.
    plt.clf()

    # return None


def plot_bar_chart(x_vals, y_vals, z_values, z_name):
    """Plot a 3D Bar chart of permeability.

    Parameters
    ----------
    y_vals = (int) lines in plot
    x_vals = (int) columns in plot
    z_values = (array) data to plot - Numpy
    z_name = (str) z-axis title

    Returns
    -------
    None
    """
    # Notify user to wait.
    print(WAIT_FOR_IT)

    # set up the figure and axes.
    fig = plt.figure(figsize=(WIDTH, HEIGHT), dpi=RESOLVE,
                     num='Seal_Flux Bar Chart')
    ax1 = fig.add_subplot(111, projection='3d')

    # Create coordinate data for plot - shift for better image.
    x_loca, y_loca = aux.set_bar_arrays(y_vals, x_vals)

    # Define a range of colors
    cmap = cm.get_cmap('Spectral')
    max_height = np.max(z_values)
    min_height = np.min(z_values)
    colors = [cmap((k-min_height)/max_height) for k in z_values]

    # Plot data.
    width = 1.0
    depth = 1.0
    base = np.zeros_like(z_values)
    ax1.bar3d(x_loca, y_loca, base, width, depth, z_values, shade=True,
              color=colors)

    # Draw legend and axis labels (add space for vertical title).
    ax1.set_title('Seal Permeability', loc='left', fontdict=TITLE_FONT)
    ax1.set_xlabel('Column Number', fontdict=REG_FONT)
    ax1.set_ylabel('Row Number', fontdict=REG_FONT)
    ax1.set_zlabel('\n' + z_name, fontdict=REG_FONT)

    # Rotate view.
    ax1.view_init(45, -70)

    if SAVE_FIGURE:
        save_figure_to_file()
    else:
        print(CAUTION)
        plt.show()

    # Clear plot to prevent over-write.
    plt.clf()

    # return None


def plot_time_history(realization, x_data, y_data, si_units, max_time):
    """Plot a leakage time history (results).

    Parameters
    ----------
    realization = (int) realization number selected for plot
    x_data = (array) x-axis data = time
    y_data = (array) y-axis data = leakage
    si_units = (bool) use of metric units
    max_time = (float) maximum for x-axis plot

    Returns
    -------
    None

    Notes
    -----
    1. Uses matplotlib functions to plot single realization.
    2. Uses LaTex to format text for superscript.
    """
    # Establish window size for plotting in inches.
    plt.figure(figsize=(WIDTH, HEIGHT), dpi=RESOLVE, num='Seal_Flux History')

    # Get exponent for y-axis and normalize data - keep source!
    y_min = 0.0
    y_max = aux.round_exp(np.max(y_data))
    cite = math.trunc(math.log10(abs(y_max)))
    yed_data = copy.deepcopy(y_data)
    yed_data /= pow(10.0, cite)
    y_max /= pow(10.0, cite)

    # Plot x/y data.
    plt.plot(x_data, yed_data, linestyle='solid', color='blue', linewidth=1.5,
             marker='o', markersize=4)

    # Set axes limits to limit white space along axes.
    # -->  old: x_max_plot = math.ceil(np.max(x_data))
    x_min_plot = math.floor(np.min(x_data))
    x_max_plot = max_time
    plt.xlim(x_min_plot, x_max_plot)
    plt.ylim(y_min, y_max)

    # Erase exponent value in plot.
    pivot = plt.gca()
    pivot.yaxis.offsetText.set_visible(False)

    # Construct plot title, axis labels and grid (add exponent to tile).
    fig_legend = ("Total CO2 Migration for Realization # " + str(realization)
                  + "\n")
    plt.title(fig_legend, fontdict=LEG_FONT)
    plt.xlabel("Time (years)", fontdict=REG_FONT)

    if si_units:
        new_label = 'tonnes)'
    else:
        new_label = 'tons)'
    y_label = r'Leakage ($\times$10$^{%d}$ ' % cite + new_label  # <== Tex
    plt.ylabel(y_label, fontdict=REG_FONT)
    plt.grid(which='major', linewidth='0.5')

    # Output figure to file or screen.
    if SAVE_FIGURE:
        save_figure_to_file()
    else:
        print(CAUTION)
        plt.show()

    # Clear plot to prevent over-write.
    plt.clf()

    # return None


def save_figure_to_file():
    """Save figure of results to output file.

    Parameters
    ----------
    N/A

    Returns
    -------
    N/A
    """
    # Create filename and then save plot to file.
    file_name = "seal_flux_plot_" + str(scfg.PLOT_NO)
    _, destination = sfile.get_path_name(scfg.OUTPUT_DIRECTORY,
                                         file_name,
                                         scfg.EXTENSION_PNG)
    plt.savefig(destination)

    # Increase figure number for next plot.
    scfg.PLOT_NO += 1

    # return None


def plot_time_series(x_data, data_list, data_maximum, start_sim, si_units,
                     max_time):
    """Plot a series of leakage time histories (results).

    Parameters
    ----------
    x_data = (array of floats) x-axis data
    data_list = (list) leakage data (NumPy)
    data_maximum = (float) maximum of <all> y data
    start_sim = (int) number of first simulation for plotting
    si_units = (bool) use of metric units
    max_time = (float) defined maximum of time for plot

    Returns
    -------
    None

    Notes
    -----
    1. Uses matplotlib functions to plot single realization
    2. Uses LaTex to format text
    3. Time array (x-axis) is assumed same for all data series
    """
    # Establish window size for plotting in inches.
    plt.figure(figsize=(WIDTH, HEIGHT), dpi=RESOLVE,
               num='Seal_Flux Simulations')

    # Get exponent for y-axis - to normalize data.
    y_max = aux.round_exp(data_maximum)
    cite = math.trunc(math.log10(abs(y_max)))
    y_max /= pow(10.0, cite)

    # Plot each "normalized" data set as line with label.
    for sim, results_array in enumerate(data_list):
        describe = "Simulation #" + str(sim + start_sim)
        y_data = results_array / pow(10.0, cite)
        plt.plot(x_data, y_data, linestyle='solid', linewidth=1.0,
                 label=describe)

    # Set axes limits to minimize white space along axes.
    #  -> use input value to control time max
    # --> old: plt.xlim(math.floor(np.min(x_data)), math.ceil(np.max(x_data)))
    plt.xlim(math.floor(np.min(x_data)), max_time)
    plt.ylim(0.0, y_max)

    # Hide axis title.
    pivot = plt.gca()
    pivot.yaxis.offsetText.set_visible(False)

    # Construct plot title and provide axis labels and grid.
    legend = "Total CO2 Migration for a Series of Realizations \n"
    plt.title(legend, fontdict=LEG_FONT)
    plt.xlabel('Time (years)', fontdict=REG_FONT)
    if si_units:
        new_label = 'tonnes)'
    else:
        new_label = 'tons)'
    plt.ylabel(r'Leakage ($\times$10$^{%d}$ ' % cite + new_label,
               fontdict=REG_FONT)

    plt.grid(which='major', linewidth='0.5')

    # Plot key at upper left.
    if len(data_list) <= PLOT_LINE_LIMIT:
        # Plot key at upper left - if enough space - say 25 lines is limit.
        plt.legend(loc=2, shadow=True, fancybox=True)
    else:
        # If large number of lines, plot number of lines in top left corner.
        new_label = f'Number of Simulations = {len(data_list)}'
        plt.text(150, 600, new_label, ha='left', va='center',
                 transform=None, fontdict=REG_FONT)

    # Save or show figure.
    #   Note: WINDOW BUG: interactive window and plot window
    #                     are not be active together.
    if SAVE_FIGURE:
        # Save to output file.
        save_figure_to_file()
    else:
        print(CAUTION)
        # plt.ion()       # Allow interactive mode for follow-up questions.
        plt.show()
        # plt.waitforbuttonpress(0)  # this will wait for indefinite time.

    # Clear plot to prevent over-write.
    plt.clf()

    # return None


def plot_manager(alive, seal_controls, grid, si_units):
    """Control all plots produced by seal_flux.

    Parameters
    ----------
    alive = (bool) status of run; alive = True if stand-alone
    seal_controls = (dict) seal controls
    grid = (list) (class) a collection of cells
    si_units = (bool) use of metric units

    Returns
    -------
    None
    """
    # Check first - if in stand-alone operation mode.
    if alive:
        # Print header, if plotting desired.
        check_plot = [seal_controls['plot_permeability'],
                      seal_controls['plot_time_history'],
                      seal_controls['plot_co2_contour']]
        if True in check_plot:
            print(MAJ_IN + "PLOT OPTIONS ENABLED.")

        # Limit selection of plot numbers in query.
        max_plots = (seal_controls['realizations'] - 1)

        # -----------------------------------------------------
        # A) Permeability initial plot option.
        if seal_controls['plot_permeability']:

            print(MAJ_IN + "LAST SIMULATION PERMEABILITY PLOTS CREATED.")
            control_permeability_plot(seal_controls, grid)

        # -----------------------------------------------------
        # B) Time history plot options.
        if seal_controls['plot_time_history']:

            # Single time history  *****
            # Print without return on print as "input" provides new line.
            print(MAJ_IN + "3. SEAL PLOT OF SINGLE TIME-HISTORY.")
            # Plot single time-history, if desired.
            initial_query = 0
            response = aux.create_history_query(initial_query, False)
            if response:
                control_time_history(max_plots, si_units,
                                     seal_controls['max_draw_time'])

            # Multiple time histories  *****
            # Print without return on print as "input" provides new line.
            print(MAJ_IN + "4. SEAL MULTIPLE TIME-SERIES PLOTS.")

            # Plot a series of time-histories, if desired.
            response = aux.create_history_query(initial_query, True)
            if response:
                control_time_series(max_plots, si_units,
                                    seal_controls['max_draw_time'])

        # -----------------------------------------------------
        # C) CO2 contour flux plot option.
        if seal_controls['plot_co2_contour']:

            # Print without return on print as "input" provides new line.
            print(MAJ_IN + "5. SEAL CO2 FLUX CONTOUR PLOT.", end="")
            initial_query = 0
            response = aux.create_contour_query(initial_query)
            if response:
                control_co2_contour(max_plots, seal_controls, grid)

    # return None


#
# -----------------------------------------------------------------------------
# - End of module
