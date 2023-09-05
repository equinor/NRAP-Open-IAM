#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions to plot fault history data.

Author: Ernest N. Lindner
Date: 08/23/2022

Module Name
    flt_plot

Contents (12)
    test_data(z_data)
    check_series(simulation, final_max)
    round_exp(exp_value)
    setup_time_history_data(simulation)
    simulation_query(max_realizations, start_value=0)
    create_history_query(query_stage, series)
    range_query(max_realizations)
    control_time_series(max_realizations, existence)
    run_plot (time_data, data_list, final_max, plotted, existing_run)
    save_figure_to_file()
    ----
    plot_time_series(x_data, data_list, data_maximum, plotted)
    plot_manager(alive, fault_controls, existence)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import math                     # For normal distribution
import csv                      # For saving a list of lists
import numpy as np              # For array operations
import matplotlib.pyplot as plt   # For plots

import flt_file as fileop       # For paths and error reports
import flt_config as scfg       # IO directory and file names


# Plot constants for figures
HEIGHT = 8                      # Figure height (inches)
WIDTH = 12                      # Figure width (inches)
LEGEND_FONT = 13                # Font size for figure title (pt)
TYP_FONT = 12                   # Font size for axes numbers (pt)
RESOLVE = 90                    # Figure display resolution (dpi)
VERY_SMALL = 1.0E-16            # limit on plot values

# Title Plot - Time history series
LEGEND_SERIES = "\n\n Fault-Flo Results for CO2 Migration through a Fault \n"

# Font definitions
REG_FONT = {'fontname': 'Arial',    # X-axis, Y-axis font
            'color':    'black',
            'size':      TYP_FONT}
LEG_FONT = {'fontname': 'Arial',    # Title font
            'color':    'black',
            'size':      LEGEND_FONT}

# General Statements
MAJ_IN = "\n  --> "                             # Page + general message start
LEAD_IN = "      --> "                          # General screen message start
CMD_IN = "      >>> "                           # Command screen message start
ERR_IN = "      *** "                           # Error message start
START_TXT = LEAD_IN + "START Realization >"     # Multi-history plot start
END_TXT = LEAD_IN + "END Realization >"         # Multi-history plot end
CAUTION = (CMD_IN + "If Paused, "
           + "Hit Exit at Top-Right of Plot Window to Proceed ...")
AGAIN = "Please Try Again."                     # Simple repeat request
AGAIN2 = "\n" + CMD_IN + AGAIN                  # Repeat request with new line
EMPTY = ("\n" + CMD_IN + "Y Value Essentially Zero For Selected Simulations!!"
         + "\n" + CMD_IN + "- No Plot; Try Again.")

# Controls and warnings
NEGATIVES = ('n', 'N', 'q', 'Q')                # Keys for negative answer
EXIT_COMMAND = ('x', 'X')                       # Keys that will exit
DOUBLES = ('nn', 'NN', 'qq', 'QQ', 'xx', 'XX', 'yy', 'YY')  # Common errors
POSITIVES = ('y', 'Y')
STAR_IN = "\n         ?-> "                     # Question start
CIRC_IN = "      --> "                          # Other question/statement
MORE_IN = "\n    ?-> "                          # Question start + newline
SERIES_START = MORE_IN + "Create a Time-Series Plot? ('Q'=quit): "
SERIES_END = MORE_IN + "Create Another Time-Series Plot? ('Q'=quit): "
SAVE_FIGURE = True                             # Save plots to file


def test_data(z_data, max_all_data):
    """Examine data to prevent an error in time series plot.

    Parameters
    ----------
    z_data = (Numpy array) data to plot
    max_all_data = maximum y value of all data

    Returns
    -------
    answer = (bool): True = data OK; False = error / do not plot
    """
    # Examine if spread in data - check min/max.
    answer = True
    z_min = np.min(z_data)
    z_max = np.max(z_data)

    if z_max <= VERY_SMALL or z_min == z_max:
        # No plot possible!
        answer = False
    elif max_all_data <= VERY_SMALL:
        # No good plot possible!
        answer = False

    return answer


def check_series(simulation, final_max):
    """Provide error check of series data.

    Parameters
    ----------
    simulation = (Numpy array) one simulation from analysis
    final_max = maximum value of all data

    Returns
    -------
    match_err = (bool) error in time series
    data_err = (bool) error in data
    """
    # Error check time data (use one series to check)!
    # --> Check that the time string is greater than 1.
    leaker, time_match = \
        setup_time_history_data(simulation)

    # --> Data range error - check on data errors.
    err_check = test_data(leaker, final_max)
    if err_check:
        data_err = False
    else:
        data_err = True        # Error found

    # --> Time series error - check on length of time series.
    if len(time_match) > 1:
        match_err = False
    else:
        match_err = True      # Error found

    return match_err, data_err


def round_exp(exp_value):
    """Round-up exponential to desired degree.

    Parameters
    ----------
    exp_value = (float) axis value

    Returns
    -------
    target = (float) rounded exponential
    """
    # Trap error if = 0.0.
    if exp_value != 0.0:
        # Make the number a decimal.
        core = abs(exp_value)
        level = math.trunc(math.log10(core))
        base = math.pow(10.0, level)
        root_exp = math.ceil(exp_value / base)

        # Round and then reconstitute the number.
        # adjust_number = round(root_exp, 1)
        target = root_exp * base  # round to zero decimal
    else:
        target = exp_value

    return target


def setup_time_history_data(simulation):
    """Get time history data from file.

    Parameters
    ----------
    simulation = (int) simulation number to plot

    Returns
    -------
    leakage_data = (NumPy array) CO2 data
    time_data = (NumPy array) time values
    """
    # Construct path name to file.
    sim_number = str(simulation)
    file_name = scfg.RESULTS_NAME + sim_number
    subdirectory_path, destination = fileop.get_path_name(scfg.OUTPUT_DIR,
                                                          file_name,
                                                          scfg.EXTENSION_CSV)
    try:
        with open(destination, 'r', encoding="utf8") as csvfile:
            csv_reader = csv.reader(csvfile)

            # Skip 2 header lines.
            header = next(csv_reader, None)
            next(csv_reader)

            # Do not read file without a header / empty file.
            if header is None:
                # No data.
                fileop.data_error("Data Error - Array Header is None! "
                                  + "-> No data found in file.",
                                  subdirectory_path, file_name)
            else:
                # Create data array - iterate over each row after the header.
                data_list = []
                for row in csv_reader:
                    data_list.append(row)

        # Convert list to NumPy array.
        npa = np.asarray(data_list, dtype=np.float64)

        # Slice array for data.
        time_data = npa[:, 1]
        leakage_data = npa[:, 2]

    except OSError as err:
        fileop.io_snag(err, subdirectory_path, file_name)
        time_data = np.array([])  # for inspection
        leakage_data = np.array([])   # for inspection

    return leakage_data, time_data


def simulation_query(max_realizations, start_value=0):
    """Get the simulation number to plot.

    Parameters
    ----------
    max_realizations = (int) maximum number of realizations
    start_value = (int) minimum number of realizations

    Returns
    -------
    simulation = (int) simulation number to plot - integer
    """
    # Default values.
    repeat_loop = True
    response = False
    realization = -1
    entered = ""

    # Loop to get an appropriate response.
    while repeat_loop:
        try:
            # Get input with prompt.
            code = input(STAR_IN + "Enter Realization Number "
                         + "to Plot ('Q'=quit): ")
            # Check if user wishes to quit (alpha character).
            if code in NEGATIVES or code in EXIT_COMMAND:
                print(CIRC_IN + "Exiting Plot Option. \n")
                repeat_loop = False
                break
        except ValueError:
            # Parse fails.
            print(ERR_IN + "Invalid Number.  " + AGAIN2)
            continue

        # Check if "entered" is a number and then within correct range.
        try:
            entered = int(code)
            if entered > max_realizations:
                print(ERR_IN + f'You typed = {entered}')
                print(ERR_IN + "This Number Exceeds the Maximum "
                      + f' of {max_realizations}! ', end='')
                print(AGAIN)
            elif entered <= 0:
                print(ERR_IN + f'You typed = {entered}')
                print(ERR_IN + "The Input is Equal to, or Less than Zero! ",
                      end='')
                print(AGAIN)
            elif entered <= start_value:
                print(ERR_IN + f'You typed = {entered}')
                print(ERR_IN + "The Input is Less Than Starting Value" +
                      f' of {start_value}!', end='')
                print(AGAIN)
            else:
                # Input OK!
                # print(CIRC_IN + f"Simulation Selected = {entered}")
                realization = entered
                response = True
                repeat_loop = False
                break
        except ValueError:
            print(ERR_IN + f'You typed = {entered}')
            print(ERR_IN + "This is Not a Number!! ", end='')
            print(AGAIN)
            continue

    # end while

    return response, realization


def create_history_query(query_stage):
    """Check if user wants to plot a time series.

    Parameters
    ----------
    query_stage = (int) question queue
        0 = first time
        1 = plot again?

    Returns
    -------
    reply = answer plot question
    """
    # Change message depending on place in query queue.
    if query_stage == 0:
        code = input(SERIES_START)
    else:
        code = input(SERIES_END)

    # correct typing error.
    if code in DOUBLES:
        code = code[0]

    # Get response.
    if code in POSITIVES:
        reply = True
    elif NEGATIVES or code in EXIT_COMMAND:
        reply = False
    else:
        reply = False   # Default is negative!

    return reply


def range_query(max_realizations):
    """Get simulation number to plot.

    Parameters
    ----------
    max_realizations = (int) maximum number of realizations

    Returns
    -------
    response = to plot or not - Boolean
    min_simulation = minimum simulation to plot - integer
    max_simulation = maximum simulation to plot - integer
    """
    # Default values.
    end_run = 0

    # Get start simulation for start of run with prompt.
    print(START_TXT, end='')

    # Get simulation number.
    starter = 0
    response, start_run = simulation_query(max_realizations, starter)

    # If OK, get end of range of simulations.
    if response:
        print(END_TXT, end="")
        starter = start_run
        response, end_run = simulation_query(max_realizations, starter)

    return response, start_run, end_run


def control_time_series(max_realizations, existence):
    """Provide overall control to plot a series of simulation from output.

    Parameters
    ----------
    max_realizations = (int) upper limit of simulations from analysis
    existence = (dict) existence of fault

    Returns
    -------
    N/A
    """
    # Loop until finished plotting.
    while True:
        # Query to get simulation numbers (start and end).
        response, start_run, end_run = range_query(max_realizations)

        # Plot - if user desires a plot.
        if response:
            # Get data and save NumPy arrays in list.
            #   -- adjust counter to catch last run.
            data_list = []
            time_data = []
            final_max = 0.0

            # Construct plot data of series.
            existing_run = start_run
            plotted = []

            for realization in range(start_run, (end_run + 1)):
                # Only include a plot if fault data exists.
                # --> Existence list starts at 0!
                if existence[realization - 1]:
                    leakage_data, time_data = \
                        setup_time_history_data(realization)
                    data_list.append(leakage_data)
                    instance_max = leakage_data.max()
                    existing_run = realization
                    if instance_max > final_max:
                        final_max = instance_max
                    plotted.append(realization)

            # Check and plot data.
            run_plot(time_data, data_list, final_max, plotted, existing_run)

            # Ask about another plot; exit if not.
            plot_again = create_history_query(1)
            if not plot_again:
                break
        else:
            break
    # end while loop

    # return None


def run_plot(time_data, data_list, final_max, plotted, existing_run):
    """Check data and plot a series of simulation from output.

    Parameters
    ----------
    time_data = (NumPy array) time values
     data_list = (list of arrays) a list of NumPy arrays of leakage
    final_max = (float) maximum of <all> y data
    plotted = (list) of runs
    existing_run = (NumPy array) of initial data array

    Returns
    -------
    N/A
    """
    # Check data!
    match_err, data_err = check_series(existing_run, final_max)

    if len(plotted) > 0:
        print(LEAD_IN + "Events Plotted with Faults = ", plotted)
    else:
        print(LEAD_IN + "No Simulations in Range with Faults!")

    # Plot data arrays on same plot diagram if data OK.
    if not match_err and not data_err:
        # Data is AOK - Plot!
        plot_time_series(time_data, data_list, final_max, plotted)

    elif match_err:
        # Time string error - terminate plotting.
        reason = "Error: No Time Data for Plot!"
        fileop.opx_problem(reason, err='')
    else:
        # Value error - continue logic.
        print(CMD_IN + "Plot Issue! Data Negative Or Series Has "
              + "Maximum Equal To Minimum!")
        print(CMD_IN + "Plot Attempt Is Discontinued!")

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
    file_name = "Fault_Flo-plot_" + str(scfg.PLOT_NO)
    _, destination = fileop.get_path_name(scfg.OUTPUT_DIR,
                                          file_name,
                                          scfg.EXTENSION_PNG)
    plt.savefig(destination)

    # Increase figure number for next plot.
    scfg.PLOT_NO += 1

    # return None


def plot_time_series(x_data, data_list, data_maximum, plotted):
    """Plot a series of leakage time histories (results).

    Parameters
    ----------
    x_data = (array) time data - assumed same for all
    data_list = (list of arrays) a list of NumPy arrays of leakage
    data_maximum = (float) maximum of <all> y data
    plotted = (list) of runs

    Returns
    -------
    N/A

    Notes
    -----
    1. Uses matplotlib functions to plot realizations.
    2. Uses LaTex to format text.
    """
    # Establish window size for plotting in inches.
    plt.figure(figsize=(WIDTH, HEIGHT), dpi=RESOLVE,
               num='Fault_Flux Simulations \n')

    # Get exponent for y-axis - to normalize data.
    y_max = round_exp(data_maximum)
    cite = math.trunc(math.log10(abs(y_max)))
    y_max /= math.pow(10.0, cite)

    # Plot each "normalized" data set as line with label.
    for sim, results_array in enumerate(data_list):
        run_number = plotted[sim]
        describe = "Simulation #" + str(run_number)
        y_data = results_array / math.pow(10.0, cite)
        plt.plot(x_data, y_data, linestyle='solid', linewidth=1.0,
                 label=describe)

    # Set axes limits to limit white space along axes.
    x_min_plot = math.floor(np.min(x_data))
    x_max_plot = math.ceil(np.max(x_data))
    plt.xlim(x_min_plot, x_max_plot)
    plt.ylim(0.0, y_max)

    # Hide axis title.
    pivot = plt.gca()
    pivot.yaxis.offsetText.set_visible(False)

    # Construct plot title w. Latex and provide axis labels and figure grid.
    new_label = r'Leakage ($\times$10$^{%d}$ tonnes)' % cite
    plt.title(LEGEND_SERIES, fontdict=LEG_FONT)
    plt.xlabel('Time (years)', fontdict=REG_FONT)
    plt.ylabel(new_label, fontdict=REG_FONT)
    plt.grid(which='major', linewidth='0.5')

    # Plot key of some sort.
    if len(data_list) <= 25:
        # Plot key at upper left - if enough space -  say 25 lines.
        plt.legend(loc=2, shadow=True, fancybox=True)
    else:
        # If large number of lines, plot number of lines in top left corner.
        new_label = f'Number of Simulations = {len(data_list)}'
        plt.text(150, 600, new_label, ha='left', va='center',
                 transform=None, fontdict=REG_FONT)

    # Option to save plot to file or show on console.
    if SAVE_FIGURE:
        # Save to output file.
        save_figure_to_file()
        plt.show()
    else:
        # Else show figure on console.
        #   Note: WINDOW BUG: interactive window and plot window
        #                     are not be active together.
        print(CAUTION)
        plt.show()
        # plt.waitforbuttonpress(0) # this will wait for indefinite time

    # Clear plot to prevent over-writes.
    plt.clf()

    # return None


def plot_manager(alive, fault_controls, existence):
    """Control plots produced by program.

    Parameters
    ----------
    alive = (bool) status of program; alive = True if stand-alone
    fault_controls = (dict) dictionary of fault controls
    existence = (dict) of fault existence

    Returns
    -------
    N/A
    """
    # Check if in stand-alone operation mode.
    if alive:
        # Print header, if plotting is desired.
        if fault_controls['plot_time_history']:
            print(MAJ_IN + "PLOT OPTIONS.", end='')

            # Inquire on plotting; limit selection of plot numbers in query.
            max_realizations = (fault_controls['realizations'])
            initial_query = 0
            response = create_history_query(initial_query)

            if response:
                control_time_series(max_realizations, existence)

    # return None


#
# -----------------------------------------------------------------------------
# End of module
