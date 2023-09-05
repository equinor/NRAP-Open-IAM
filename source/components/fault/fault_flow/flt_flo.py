#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Main Module - Code performs a fault analysis for NRAP IAM.

Author: Ernest N. Lindner
Date: 08/25/2022

Module Name
    flt_flo

Contents (1)
    main()

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import logging  # For reporting errors

import flt_compute as fcp   # For basic calc.s and memory ops
import flt_file as fileop   # For file operations
import flt_inout as fio     # For reservoir input
import flt_intro as intro   # Gets YAML data and checks input
import flt_message as mess  # For messages
import flt_model as fmodl   # Fault model functions
import flt_plot as fplot    # For plotting data
import flt_profile as pro   # For T/P profile calculations
import flt_reveal as revl   # Provides summary file and plots
import flt_storage as stor  # For defining lists

# Details:
__version__ = "6.7.5"
__author__ = "Ernest N. Lindner"
__contact__ = "Ernest.Lindner@netl.doe.gov"
__date__ = "2022/08/25"
__deprecated__ = False
__license__ = "Personal"
__maintainer__ = "Ernest N. Lindner"
__status__ = "Development"

"""
-------------------------------------------------------------------------------
Note
----
    This software is provided for use as-is and no representations about the
    suitability or accuracy of this software is made. It is provided without
    warranty.

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
# Input control file name
YAML_FILE = "ControlFile_fault.yaml"  # Input file name in source directory

# Python version constants
PYTHON_MAJOR = 3   # Python major version of code, op_major
PYTHON_MINOR = 10  # Python minor version of code, op_minor

# Debugging controls - also controlled by "RUNNING"
DETAIL_PERM = False   # Save permeability to file
DETAIL_STEP = False   # Save results of time step data to file
DETAIL_PLATE = False  # Save fault plate data to file

RUNNING = (__name__ == "__main__")

logging.basicConfig(format='\n  --> %(levelname)s: %(message)s',
                    level=logging.WARNING)


def main():
    """Control entire operation of code.

    Parameters
    ----------
    See Control File

    Returns
    -------
    N/A
    """
    # -------------------------------------------------------------------------
    # A. SETUP
    # -------------------------------------------------------------------------
    # Read control file and check fault_controls.
    fault_controls, param_bounds = intro.start(RUNNING, YAML_FILE,
                                               PYTHON_MAJOR, PYTHON_MINOR,
                                               __version__)

    # Create additional input values.
    fault_controls = fcp.define_fracto(fault_controls, param_bounds)

    # Initialize: establish time step array - clock_values.
    clock_values = stor.create_time_values(fault_controls)  # In years!!

    # Create flux list, sum variable and fault class.
    all_sim_flux = []
    total_co2 = 0.0
    fault = []

    # Create storage arrays for time steps within a realization.
    sim_co2_flow, sim_brine_flow, co2_flow_step, brine_flow_step = \
        stor.create_sim_storage(fault_controls['n_plates'])

    # Create output list.
    simulation_list = []
    fault_life_list = []

    # Compute interpolated conditions points.
    fault_controls = pro.setup_transition_points(fault_controls)

    # Define fluid properties depending on case / profile type.
    fault_controls = pro.define_ave_interpolation_stages(RUNNING,
                                                         fault_controls)

    # -------------------------------------------------------------------------
    # B. COMPUTATION LOOPS
    # -------------------------------------------------------------------------
    mess.echo_status(RUNNING, "RUNNING")

    # >>>>>>> REALIZATION LOOP <<<<<<<<
    for simulation_numbr in range(fault_controls['realizations']):
        # ***** DEFINE VARIABLES FOR EACH SIMULATION *****
        # Set seed for random numbers on each simulation.
        rng = fcp.random_seeder(1)
        fault_controls['rng'] = rng   # Not needed

        # Define strike, aperture, length for fault plates.
        fault_controls, fault = fcp.fault_create(fault_controls, fault, rng)

        # Check if fault exists and add to counter.
        fault_status = fcp.check_if_fault_exists(fault_controls, rng)
        if fault_status:
            fault_controls['simulations_with_fault'] += 1

        # Write fault status to records file.
        fault_life_list.append(fault_status)

        # DEBUG: Save fault plate data to file, if desired.
        if DETAIL_PLATE:
            fio.save_fault_results(simulation_numbr, fault_controls['title'],
                                   fault)

        # Echo simulation number to console.
        mess.echo_simulation_step(RUNNING, simulation_numbr)

        # Create list for simulation output --> with chrono=0.
        simulation_list = stor.create_simulation_list(simulation_numbr,
                                                      simulation_list)

        # Reset Interpolate Values for start of a simulation.
        # --- Injection starts at time=0.
        fault_controls = pro.interpolate_initialize(fault_controls)

        # ---------------------------------------------------------------------
        #  ******* CONDUCT TIME STEP LOOP *******
        # - start at #1 --> as step #0 has zero accumulation.
        for chrono in range(1, fault_controls['time_points']):

            # Provide header, if showing time step data.
            mess.echo_time_step(DETAIL_STEP, RUNNING, chrono)

            # Input reservoir pressure and saturation arrays from files.
            # --  input for all parts of fault for this time step.
            base_co2_pressure, base_co2_saturation = \
                fio.input_reservoir_data(chrono, fault_controls['n_plates'],
                                         param_bounds)

            # Update interpolation values for time step.
            fault_controls = pro.update_interpolate(clock_values[chrono],
                                                    fault_controls)

            # ------ Loop Over Fault Plates ------
            # Compute Darcy flow for each plate over one increment.
            for seg_number in range(fault_controls['n_plates']):

                # Get reservoir data for plate.
                co2_base_pressure = base_co2_pressure[seg_number]
                co2_base_saturation = base_co2_saturation[seg_number]

                # Compute flow for each plate - vol./sec.- if fault exists.
                if fault_status:

                    # if fault_controls['profile_type'] == 0:
                    # NOTE: Dissolved CO2 is included in CO2 flow calc. as
                    #       fluid despite assuming immiscible flow in theory.
                    flow_rate_co2, flow_rate_brine = \
                        fault[seg_number].compute_flow(co2_base_pressure,
                                                       co2_base_saturation,
                                                       fault_controls)

                    # *********************************************************
                    # Compute additional part of profile for complex models.
                    if fault_controls['profile_type'] == 1:
                        # Compute complex profile.
                        flow_rate_co2, flow_rate_brine = \
                            fault[seg_number].complex_tail(co2_base_pressure,
                                                           co2_base_saturation,
                                                           fault_controls,
                                                           flow_rate_co2,
                                                           flow_rate_brine)
                    if fault_controls['profile_type'] == 2:
                        # Compute disjoint profile.
                        flow_rate_co2, flow_rate_brine = \
                            fault[seg_number].complex_tail(co2_base_pressure,
                                                           co2_base_saturation,
                                                           fault_controls,
                                                           flow_rate_co2,
                                                           flow_rate_brine)
                    # *********************************************************
                else:
                    # No fault in this simulation.
                    flow_rate_co2 = 0.0E-50
                    flow_rate_brine = 0.0E-50

                # Accumulate vol. flows of each for this time step.
                co2_flow_step[seg_number] = flow_rate_co2
                brine_flow_step[seg_number] = flow_rate_brine
            # ------- End Fault Plate Loop ------

            # Convert flows into mass flows in tonnes for time step.
            co2_flow_step, brine_flow_step = \
                fmodl.convert_flows(co2_flow_step, brine_flow_step,
                                    clock_values[chrono],
                                    clock_values[chrono - 1],
                                    fault_controls)

            # Sum time step flows for cumulative flows, and compute totals.
            sim_co2_flow += co2_flow_step
            sim_brine_flow += brine_flow_step
            total_co2 = sim_co2_flow.sum()
            total_brine = sim_brine_flow.sum()

            # DEBUG: Write data to a file for <current> time step.
            if DETAIL_STEP:
                fio.cache_step_results(simulation_numbr, clock_values[chrono],
                                       fault_controls["title"],
                                       sim_co2_flow, sim_brine_flow)

            # Save results in a list of total simulation results.
            simulation_list = \
                stor.update_simulation_list(simulation_numbr,
                                            clock_values[chrono],
                                            total_co2, total_brine,
                                            simulation_list)

            # Zero-out calculation arrays for next cycle.
            co2_flow_step, brine_flow_step = \
                stor.clear_cycle_storage(co2_flow_step, brine_flow_step,
                                         fault_controls['n_plates'])
        # ******* END TIME STEP LOOP *******
        # ---------------------------------------------------------------------

        # Write the list of results for this simulation.
        fio.cache_sim_results(simulation_numbr, simulation_list,
                              fault_controls['title'],
                              fault_controls['skip_output'])

        # Save simulation CO2 flux results to a list.
        all_sim_flux.append(total_co2)

        # DEBUG: Write permeability data for simulation, if debugging.
        if DETAIL_PERM:
            fio.cache_perm(fault_controls, simulation_numbr, fault)

        # Zero-out computation arrays for next simulation loop.
        sim_co2_flow, sim_brine_flow, simulation_list = \
            stor.clear_storage(sim_co2_flow, sim_brine_flow, simulation_list)

        # Re-set interpolation values
        fault_controls = pro.refresh_interpolate_controls(fault_controls)
    # >>>>>>> END OF REALIZATION LOOP  <<<<<<<

    # -------------------------------------------------------------------------
    # C. WRITE STATUS & SUMMARY FILE, SHOW PLOTS.
    # -------------------------------------------------------------------------
    mess.echo_status(RUNNING, "SUMMARY")

    # Write run time, if stand-alone.
    fault_controls = mess.closeout_time(RUNNING, fault_controls)

    # Write dictionary values to summary file.
    revl.write_summary(fault_controls, all_sim_flux, fault)

    # Write to file the fault status for each simulation.
    mess.write_fault_list(fault_controls, fault_life_list)

    # Create dictionary of real results and show plots, if desired.
    #   -> Convert "fault status list" into a dictionary for plots.
    # Must use dict() - this will not work with {} to create a dictionary!
    existence = dict(enumerate(fault_life_list))
    fplot.plot_manager(RUNNING, fault_controls, existence)

    # Print Exit statement.
    fileop.end_message(RUNNING)

# -----------------------------------------------------------------------------
# Start MAIN.
# -----------------------------------------------------------------------------


if __name__ == "__main__":
    main()


#
# -----------------------------------------------------------------------------
# End of module
# -------1---------2---------3---------4---------5---------6---------7---------8
