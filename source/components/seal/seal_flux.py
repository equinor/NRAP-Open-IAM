#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Main module for Seal_Flux code and performs seal analysis for NRAP IAM.

Author: Ernest N. Lindner
Date: 08/25/2022

Module Name - Main Line
    seal_flux

Contents (1)
    main()

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import logging                      # For warnings

import seal_config as scfg          # config variables for I/O
import seal_decipher as dec         # Read yaml data
import seal_display as dis          # Save summary data and create plots
import seal_file as sfile           # Input/output related operations
import seal_fluids as fluid         # Functions for density & viscosity
import seal_intro as intro          # Gets YAML data and checks input
import seal_model as mod            # Functions for Class definitions
import seal_perm as perm            # Functions for permeability definition
import seal_plot as splot           # Plot results
import seal_refresh as sref         # Functions for thickness and depth
import seal_upload as sup           # Input file data
import frac_origin as fog           # Fracture setup & generation
import frac_view as view            # Save fracture summary data

# Define Details
__version__: str = '3.13.5'
__author__ = 'Ernest N. Lindner'
__contact__ = 'Ernest.Lindner@netl.doe.gov'
__copyright__ = 'Copyright (C) 2019-22 Ernest Lindner'
__date__ = '08/25/2022'
__deprecated__ = False
__license__ = 'Personal'
__maintainer__ = 'Ernest N. Lindner'
__status__ = 'Development'
"""
-------------------------------------------------------------------------------
Caution
-------
    The software is provided "As Is", without warranty of any kind, express
    or implied, including but not limited to the warranties of merchantability,
    fitness for a particular purpose and noninfringement. In no event shall
    the authors or copyright holders be liable for any claim, damages or other
    liability, whether in an action of contract, tort or otherwise, arising
    from, out of or in connection with the software or the use or other
    dealings in the software.

Warranty Disclaimer
-------------------
    This software development was funded by the United States Department
    of Energy, National Energy Technology Laboratory, in part, through a site
    support contract. Neither the United States Government nor any agency
    thereof, nor any of their employees, nor the support contractor, nor any of
    their employees, makes any warranty, express or implied, or assumes any
    legal liability or responsibility for the accuracy, completeness, or
    usefulness of any information, apparatus, product, or process disclosed,
    or represents that its use would not infringe privately owned rights.
    Reference herein to any specific commercial product, process, or service
    by trade name, trademark, manufacturer, or otherwise does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the
    United States Government or any agency thereof. The views and opinions of
    authors expressed herein do not necessarily state or reflect those of the
    United States Government or any agency thereof.
-------------------------------------------------------------------------------
"""
# Operation controls and debugging - also controlled by 'RUNNING'
DETAIL_PERM = False                 # Save permeability matrix to file
DETAIL_STEP = False                 # Save results of time step data to file
ECHO = False                        # Print progress banner of each time step
BOZO = False                        # DEBUG - Print data arrays and example
JEEZ = False                        # DEBUG - Print permeability

# Version constants
PYTHON_MAJOR = 3                    # Python major version for code
PYTHON_MINOR = 10                   # Python minor version for code
SEAL_VERSION = __version__          # Seal Flux version

RUNNING = (__name__ == "__main__")

logging.basicConfig(format='\n    > %(levelname)s: %(message)s',
                    level=logging.WARNING)


def main():
    """Manage stand-alone operation of the code.

    Parameters
    ----------
    N/A (input from yaml file)

    Returns
    -------
    None (file output)
    """
    # -------------------------------------------------------------------------
    # A. SETUP
    # -------------------------------------------------------------------------
    # Check version and read seal control YAML File.
    seal_controls = dec.boot_code(RUNNING, scfg.S_CONTROL, PYTHON_MAJOR,
                                  PYTHON_MINOR, SEAL_VERSION)

    # Establish bounds and check seal_controls data.
    seal_controls, param_bounds = intro.preliminaries(RUNNING, seal_controls)

    # Assign variables to class and complete other input tasks.
    grid, seal_controls = intro.finish_definitions(seal_controls, param_bounds)

    # Input control parameters from fracture YAML file and check.
    frac_controls = fog.frac_launch(RUNNING, scfg.F_CONTROL, seal_controls)

    # -------------------------------------------------------------------------
    # B. INITIALIZATION
    # -------------------------------------------------------------------------
    sfile.echo_status(RUNNING, "INITIALIZING PARAMETERS.")
    sim_flux = []

    # Set seed for random numbers.
    rng = sup.random_seeder(0)
    frac_controls['rng'] = rng

    # Input data from major files, as desired with error check.
    grid = sup.input_various_data(seal_controls, grid, param_bounds, rng)

    # Define pressure at top of cells in either of two ways.
    # * Static pressure requires depth to top of cell.
    press_top = sup.define_top_press_array(seal_controls, grid, param_bounds)

    # Interpolate fluid properties, if desired.
    seal_controls = fluid.manage_interpolation(RUNNING, seal_controls,
                                               press_top, param_bounds)

    # Detailed DEBUG here, if BOZO = True: print arrays.
    sref.debug_check(BOZO, grid, seal_controls)

    # -------------------------------------------------------------------------
    # C. COMPUTATION LOOPS
    # -------------------------------------------------------------------------
    sfile.echo_status(RUNNING, "STARTING ANALYSIS COMPUTATIONS.")

    # Establish time step array. In years!!
    time_period = sref.create_time_values(seal_controls)

    # Create storage arrays for time steps within a realization.
    sim_co2_flow, sim_brine_flow, rate_co2, rate_brine = \
        sref.create_sim_storage(seal_controls)

    # Create reservoir input arrays.
    #   base_co2_saturation, base_co2_pressure = \
    #       sref.create_reservoir_storage(seal_controls)

    # REALIZATION LOOP <<<<<<<<
    for sim_step in range(seal_controls['realizations']):

        # Write status to console.
        sref.echo_sim_step(RUNNING, sim_step)

        # Compute permeability and threshold values for every cell.
        if seal_controls['fracture_approach']:
            # Evaluate fractures to get permeability.
            grid = fog.evaluate_fracs(grid, frac_controls, sim_step, RUNNING)
        else:
            # Evaluate parameters to get permeability.
            grid = perm.obtain_permeability(seal_controls, grid, param_bounds,
                                            rng)
            grid = perm.correlate_entry_pressure(seal_controls, grid)

        if JEEZ:
            print("\n  >>> DEBUG: Permeabilities for Each Cell in Grid")
            for cell in grid:
                print("    " + str(cell.permeability))

        # Define permeability heterogeneities, if desired.
        grid = perm.evaluate_areal_heter(grid, seal_controls)

        # Establish realization list for output.
        total_co2 = 0.0
        simulation_list = sref.establish_simulation_list(sim_step)

        # START TIME STEP LOOP <<<<<<<<
        for chrono in range(1, seal_controls['time_points']):

            # Values computed at end of first step -> step = #1!
            current = time_period[chrono]

            # Provide header, if showing time steps.
            sref.echo_time_step(ECHO, RUNNING, current)

            # Input reservoir pressure and saturation arrays from files.
            base_co2_pressure, base_co2_saturation = \
                sup.input_reservoir_data(chrono, seal_controls, param_bounds)

            # Compute Darcy-flow for each cell over one increment.
            for cell_number in range(seal_controls['num_cells']):

                # Get reservoir data for this cell.
                co2_base_pressure = base_co2_pressure[cell_number]
                co2_base_saturation = base_co2_saturation[cell_number]

                # Compute flow rates for each cell - vol./sec.
                # NOTE: Top saturation is not used at this version;
                #   bottom saturations are assumed to be steady-state value.
                # NOTE: Dissolved CO2 is included in CO2 flow calculation as
                #       fluids despite assuming immiscible flow as theory.
                flow_rate_co2, flow_rate_brine = \
                    grid[cell_number].compute_flow(co2_base_pressure,
                                                   co2_base_saturation,
                                                   press_top[cell_number],
                                                   seal_controls)

                # Update history of cell.
                grid[cell_number].track_history(co2_base_pressure,
                                                flow_rate_co2, current)

                # Update influence factors for cell.
                grid[cell_number].compute_model_influence(flow_rate_co2)

                # Accumulate vol. flows of each cell for this time step.
                rate_brine[cell_number] = flow_rate_brine
                rate_co2[cell_number] = flow_rate_co2

            # Convert flows into mass flows - in tonnes for interval.
            rate_co2, rate_brine = \
                mod.convert_flows(rate_co2, rate_brine, current,
                                  time_period[chrono - 1])

            # Sum time step flows with prior flows for cumulative flows.
            sim_co2_flow += rate_co2
            sim_brine_flow += rate_brine

            # Write 2D data to a file for <current> time step, if debugging.
            if DETAIL_STEP:
                dis.cache_step_results(sim_step, current, seal_controls,
                                       sim_co2_flow, sim_brine_flow)

            # At last time step of sim - only save results for contour plot.
            if chrono == (seal_controls['time_points'] - 1):
                dis.store_grid_results(sim_step, sim_co2_flow, sim_brine_flow,
                                       RUNNING, seal_controls)

            # Sum flows over grid for single total flow for time step.
            total_co2, total_brine = sref.aggregate_results(sim_co2_flow,
                                                            sim_brine_flow)

            # Store current results in a list of total simulation results.
            simulation_list = \
                sref.update_simulation_list(sim_step, time_period[chrono],
                                            total_co2, total_brine,
                                            simulation_list)

            # Flatten and zero-out calculation arrays for next cycle.
            rate_co2, rate_brine = \
                sref.clear_cycle_storage(rate_co2, rate_brine)

            # END TIME STEP LOOP <<<<<<<<

        # Write the list of results for this simulation.
        dis.write_sim_results(sim_step, simulation_list,
                              seal_controls['title'])

        # Write permeability data for simulation, if debugging.
        if DETAIL_PERM:
            dis.cache_perm(sim_step, grid, seal_controls)

        # Save simulation CO2 flux results to a list.
        sim_flux.append(total_co2)

        # Zero-out computation arrays for next simulation loop.
        sim_co2_flow, sim_brine_flow, simulation_list = \
            sref.clear_storage(sim_co2_flow, sim_brine_flow, simulation_list)

        # Change random seed for next simulation.
        rng = sup.random_seeder(1)
        frac_controls['rng'] = rng

        # Compute new depth, thickness and pressure values for next loop.
        press_top = \
            sref.refresh_data_arrays(seal_controls, sim_step, grid, press_top,
                                     rng)
    # END OF REALIZATION LOOP <<<<<<<<

    # -------------------------------------------------------------------------
    # D. AT RUN END: WRITE STATUS & SUMMARY FILE.
    # -------------------------------------------------------------------------
    # Compute run time and write, if stand-alone.
    seal_controls['elapsed_time'] = sref.closeout(RUNNING,
                                                  seal_controls['code_start'])

    # Define units as strings.
    uts = sup.define_units_for_output(seal_controls['use_si_units'])

    # Write seal dictionary to seal summary file.
    dis.write_summary(seal_controls, sim_flux, SEAL_VERSION,
                      uts)
    sfile.echo_status(RUNNING, "CREATED SEAL SUMMARY FILE.")

    # Write fracture dictionary values to frac summary file.
    view.write_frac_summary(frac_controls, RUNNING,
                            seal_controls['fracture_approach'], uts,
                            seal_controls['use_si_units'])

    # Show plots, if desired.
    splot.plot_manager(RUNNING, seal_controls, grid,
                       seal_controls['use_si_units'])

    # Print Exit statement.
    sfile.terminate_code(0, RUNNING)

# -----------------------------------------------------------------------------
# Start MAIN Line of Code.
# -----------------------------------------------------------------------------


if __name__ == "__main__":
    main()


#
# -----------------------------------------------------------------------------
# - End of module
# -------1---------2---------3---------4---------5---------6---------7---------8
