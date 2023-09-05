"""
Example illustrates the increase in processing time as a function of
additional generic reservoir components. The locations of the reservoirs,
as well as the distance from the reservoir where the pressure and
CO2 saturation are measured, are random.

The maximum number of reservoirs to simulate can be specified during setup
of the simulation as an integer value (n). If not provided, the default
value for the maximum number of simulated reservoirs is 10. Additionally,
the ceiling for this value has been set at 30 due to limited system
resources and extensive runtime, but can be raised in the program if
necessary.

Author: Paul Holcomb

Example of run:
$ python iam_sys_generic_reservoir_multilocs_sims.py (n)
"""

import sys
import os
import time

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import SystemModel, GenericReservoir
from matk import pyDOE

if __name__=='__main__':

    # For multiprocessing in Spyder
    __spec__ = None

    # For saving graphs
    to_save = True

    # If saving graphs, define output directory for plots
    if to_save:
        out_dir = os.path.join(os.getcwd(), 'generic_reservoir_plots\\')
        if not os.path.exists(out_dir):
            os.mkdir(os.path.dirname(out_dir))

    # Define keyword arguments of the system model
    num_years = 75   # in years
    time_array = 365.25*np.arange(0.0, num_years+1)  # convert to days
    sm_model_kwargs = {'time_array': time_array}   # time is given in days

    # Create a range of Generic Reservoirs (n)
    # If the user did not supply a number of reservoirs
    if len(sys.argv) == 1:
        # use default of 10
        n = 10

    else:  # otherwise
        # get the number of reservoirs specified
        n = sys.argv[1]
        # test to see if the user input is an integer
        try:
            n = int(n)
        # and if not, exit with error message
        except:
            sys.exit('The number of generic reservoirs must be specified as an integer.')

        # test if number of reservoirs exceeds the maximum allowed
        if n > 30:
            sys.exit('The number of generic reservoirs may not exceed 30 due to computational limits.')

    print("Number of Iterations: {}".format(n))

    # Set array to hold length of simulations
    times = []

    # Run multiple reservoir simulations, increasing the number of reservoirs by 1
    # for each simulation up to n
    for num_reservoirs in range(1, n+1):
        print("Iteration {}".format(num_reservoirs))

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)
        print("\tSystem model created.")

        # Create randomly located generic reservoir locations
        # within box defined by xmin, xmax, ymin, ymax
        xymins = np.array([100., 100.])
        xymaxs = np.array([15000., 15000.])
        reservoir_xys = xymins + pyDOE.lhs(2, samples=num_reservoirs)*(xymaxs-xymins)
        print("\tRandom reservoir locations created.")

        # Create arrays to hold reservoir component models and distances
        sres = []
        distances = []

        # For each reservoir
        for i, crds in enumerate(reservoir_xys):

            # Get random distance from reservoir to measure pressure and CO2 saturation
            curr_dist = np.random.randint(500, 2000)

            # Append random distance to list of distances
            distances.append(curr_dist)

            # Create reservoir component
            print("\t\tReservoir {}:".format(i+1))
            sres.append(sm.add_component_model_object(
                GenericReservoir(name='gres'+str(i), parent=sm,
                                 injX=crds[0], injY=crds[1],
                                 locX=crds[0]+curr_dist, locY=crds[1])))
            print("\t\t\tReservoir #{} created.".format(i+1))

            # Add parameters of reservoir component model
            sres[-1].add_par('reservoirDepth', value=2000.0, vary=False)
            sres[-1].add_par('reservoirThickness', value=50.0, vary=False)
            sres[-1].add_par('logResPerm', value=-14.0, vary=False)
            sres[-1].add_par('resTempGradient', value=30.0, vary=False)
            sres[-1].add_par('injRate', value=100, vary=False)
            sres[-1].add_par('initialSalinity', value=0.05, vary=False)
            sres[-1].add_par('wellRadius', value=0.05, vary=False)
            sres[-1].add_par('reservoirPorosity', value=0.1, vary=False)
            print("\t\t\tReservoir #{} parameters added.".format(i+1))

            # Add observations of reservoir component model to be used by the next component
            sres[-1].add_obs_to_be_linked('pressure')
            sres[-1].add_obs_to_be_linked('CO2saturation')
            print("\t\t\tReservoir #{} linked observations created.".format(i+1))

            # Add observations (output) from the reservoir model
            sres[-1].add_obs('pressure')
            sres[-1].add_obs('CO2saturation')
            print("\t\t\tReservoir #{} output observations created.".format(i+1))

        # Get start time for simulation
        t0 = time.time()

        # Run forward simulation
        sm.forward()

        # Get total time for simulation
        t_final = np.round(time.time() - t0, 3)

        # Append time for current simulation to list of simulation times
        times.append(t_final)
        print("\tGeneric Reservoir {} reservoir model running time: {} secs".format(
            num_reservoirs, t_final))

        # Collect results
        print("\tCollecting results...")

        # Create arrays for pressure and saturation for each reservoir
        pressure_at_reservoirs = []
        saturation_at_reservoirs = []

        # Collect pressure and CO2 saturation for each reservoir
        for wnum, p in enumerate(sres):
            print("\t\tReservoir #{} result collection...".format(wnum+1), end="")

            pressure_at_reservoirs.append(sm.collect_observations_as_time_series(p, 'pressure'))
            saturation_at_reservoirs.append(sm.collect_observations_as_time_series(p, 'CO2saturation'))

            print("complete.")

        # -------------------------------------------------------------------------
        # # Plot pressure and saturation at a random distance from each reservoir
        # -------------------------------------------------------------------------

        # Set the starting reservoir number
        reservoir_num = 1

        # Set graph font size
        font_size = 6
        label_size = 6

        # Set graph line widths
        grid_line_width = 0.1
        line_width = 0.5

        # Define time array for graphing
        time_array_yr = time_array / 365.25

        # Graph the pressure and saturation for each reservoir
        for pressure_at_reservoir, saturation_at_reservoir, distance in zip(
                pressure_at_reservoirs, saturation_at_reservoirs, distances):

            # convert pressure to MPa
            pressure_MPa = pressure_at_reservoir * 1e-6

            # create figure with two vertically stacked subplots
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3, 5), dpi=200)

            # PRESSURE PLOT
            # plot pressure vs time
            ax1.plot(time_array_yr, pressure_MPa, linewidth=line_width,
                     color='blue', label='at {:.1f} m'.format(distance))

            # create legend for graph in upper right
            ax1.legend(loc='upper right', fontsize=font_size)

            # set label for x-axis
            ax1.set_xlabel('Time,yr', fontsize=font_size)

            # set label for y-axis
            ax1.set_ylabel('Pressure, MPa', fontsize=font_size)

            # set min and max limits for y-axis based on pressure data
            ax1.set_ylim((np.max([np.min(pressure_MPa) - 3, 0]), np.max(pressure_MPa) + 3))

            # set min and max limits for x-axis based on time array
            ax1.set_xlim((np.min(time_array_yr), np.max(time_array_yr)))

            # format graph ticks
            ax1.tick_params(axis='both', which='both', length=0, labelsize=label_size)

            # format graph grid lines
            ax1.grid(color='grey', linestyle='-', linewidth=grid_line_width)

            # SATURATION PLOT
            # plot saturation vs time
            ax2.plot(time_array_yr, saturation_at_reservoir, linewidth=line_width, color='red')

            # set label for x-axis
            ax2.set_xlabel('Time,yr', fontsize=font_size)

            # set label for y-axis
            ax2.set_ylabel(r'CO$_2$ saturation', fontsize=font_size)

            # set min and max limits for y-axis
            ax2.set_ylim((-0.05, 1.05))

            # set min and max limits for x-axis based on time array
            ax2.set_xlim((np.min(time_array_yr), np.max(time_array_yr)))

            # format graph ticks
            ax2.tick_params(axis='both', which='both', length=0, labelsize=label_size)

            # format graph grid lines
            ax2.grid(color='grey', linestyle='-', linewidth=grid_line_width)

            # minimize white space in figure
            fig.tight_layout()

            # if saving figures
            if to_save:
                print("\tSaving plot for reservoir {}...".format(reservoir_num), end="")

                # save figure based on total number of reservoirs in current simulation (num_reservoirs)
                # and the label of the current reservoir (reservoir_num)
                fig.savefig(os.path.join(
                    out_dir, 'generic_reservoir_results_{}reservoirs_reservoir{}.png'.format(
                        num_reservoirs, reservoir_num)), dpi=500, bbox_inches='tight')

                print("complete.\n")

            # increment the current reservoir number
            reservoir_num += 1

    # GRAPH SIMULATION TIMES
    print("SIMULATION TIMES:")

    # For each simulation, print the number of reservoirs and time to completion
    for iter_num, sim_time in enumerate(times):
        print("Reservoirs: {}\tSimulation Time: {}".format(iter_num+1, sim_time))

    # Create figure for time vs number of reservoirs
    plt.figure()

    # Plot time vs number of reservoirs
    plt.plot(np.linspace(1, n, n), times, linewidth=line_width, color='blue')

    # Set x-axis label
    plt.xlabel('Number of Reservoirs in Simulation')

    # Set y-axis label
    plt.ylabel('Time (s)')

    # if saving figures
    if to_save:
        print("\tSaving plot for time series...", end="")

        # save figure to file in output directory
        plt.savefig(os.path.join(
            out_dir, 'generic_reservoir_multi-reservoir_results.png'),
            dpi=500, bbox_inches='tight')
        print("complete.")

    else:  # otherwise
        # show figure
        plt.show()
