#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions to write to summary file and show code progress.

Author: Ernest N. Lindner
Date: 08/18/2022

Module Name
    seal_display

Contents (20)
    recreate_array_shape(seal_controls, source_array)
    save_csv_results(file_name, data_block, data_title)
    save_npy_results(file_name, data_block)
    store_grid_results(sim, mass_flow_co2, mass_flow_brine,
                       allow, seal_controls)
    cache_step_results(sim, time_value, seal_controls, mass_flow_co2,
                       mass_flow_brine)
    cache_perm(sim_numbr, grid, seal_controls)
    write_sim_results(sim_numbr, sim_listing, entitle)
    write_model_parameters(zum, seal_controls)
    write_description_parameters(zum, seal_controls, uts)
    write_fluid_parameters(zum, seal_controls, uts)
    ---
    write_depth_parameters(zum, seal_controls, uts)
    write_control_values(zum, seal_controls)
    write_other_control_values(zum, seal_controls)
    write_permeability(zum, seal_controls)
    write_relative_perm(zum, seal_controls, uts)
    write_reactivity(zum, seal_controls)
    write_thickness(zum, seal_controls, uts)
    write_co2_results(zum, sim_flux_list, uts, seal_controls, si_units)
    write_last(zum, seal_controls)
    write_summary(seal_controls, sim_flux_list, version_num, uts)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
from datetime import datetime       # For current time
import logging                      # For reporting errors
import csv                          # For saving a list of lists
import copy                         # For array adjustment
import numpy as np                  # For array operations

import seal_config as scfg           # IO directory and file names
import seal_file as sfile           # For paths and error reports
import seal_units as sun            # For unit conversion

# Lead-in characters for text - for uniformity
BLK_IN = "\n  > "                   # Block start
LIN_IN = "  --> "                   # Line start
OTH_IN = "   >> "                   # Indented
HED_IN = "  > "                     # Title
PERMIT = True                       # Control to print output to file
MAX_FLUX_LIMIT = 20.0               # Percentage limit of total for code


def recreate_array_shape(seal_controls, source_array):
    """Create a new array from source for output.

    Parameters
    ----------
    seal_controls = (dict) seal parameters
    source_array = (array) source array

    Returns
    -------
    results = (array) reshaped new array (1D or 2D)
    """
    # Create new array using np.copy(a).
    results = copy.deepcopy(source_array)

    # Define lines/columns for output, if grid assumed.
    if seal_controls['grid_approach']:
        n_cols = seal_controls['grid_cols']
        n_rows = seal_controls['grid_rows']
        results.shape = (n_rows, n_cols)

    else:
        # Change array only if 2D; otherwise can cause error.
        if len(results[0]) != 1:
            # Ensure 1D array for freeform data.
            results.flatten()

    return results


def save_csv_results(file_name, data_block, run_title, data_title):
    """Store seal output in a *.csv file.

    Parameters
    ----------
    file_name = (str) Name of output file (string)
    data_block = (2D array) data block
    run_title = (str) title of run
    data_title = (str) header for data file

    Returns
    -------
    None
    """
    # Construct full path for file.
    sub_path, destination = sfile.get_path_name(scfg.OUTPUT_DIRECTORY,
                                                file_name,
                                                scfg.EXTENSION_CSV)

    # Write data - with two header lines.
    overall_title = " > Run:  " + run_title + "\n "
    overall_title += data_title
    try:
        np.savetxt(destination, data_block, fmt="%15.5e", delimiter=",",
                   comments=" ", header=overall_title)
    except OSError as err:
        sfile.io_snag(err, sub_path, file_name)

    # return None


def save_npy_results(file_name, data_block):
    """Store seal output in a NumPy file.

    Parameters
    ----------
    file_name = (str) Name of output file
    data_block = (2D array) data array

    Returns
    -------
    None
    """
    # Construct full path for file.
    sub_path, destination = sfile.get_path_name(scfg.OUTPUT_DIRECTORY,
                                                file_name,
                                                scfg.EXTENSION_NPY)

    # Write data to path.
    try:
        np.save(str(destination), data_block)
    except OSError as err:
        sfile.io_snag(err, sub_path, file_name)

    # return None


def store_grid_results(simulation, mass_flow_co2, mass_flow_brine,
                       allow, seal_controls):
    """Store NumPy file of grid results from a simulation.

    Parameters
    ----------
    simulation = (int) realization number
    mass_flow_co2 = (2D array) CO2 mass for realization
    mass_flow_brine = (2D array) brine mass for realization
    allow = (bool) stand-alone operations of code
    seal_controls = (dict) seal parameters

    Returns
    -------
    None

    Notes
    -----
    Option to create brine flux file is disabled.
    """
    # Conduct operation only if stand-alone and contour is desired.
    if allow:
        if seal_controls['plot_co2_contour']:

            # Construct file names; eliminate "zero" simulation.
            insert = str(simulation + 1)
            output1 = scfg.CO2_INTRO + insert + scfg.EXTENSION_NPY
            output2 = scfg.BRINE_INTRO + insert + scfg.EXTENSION_NPY

            # Store time-variant data for each cell for time X.
            save_npy_results(output1, mass_flow_co2)
            save_npy_results(output2, mass_flow_brine)

    # return None


def cache_step_results(simulation, time_value, seal_controls, mass_flow_co2,
                       mass_flow_brine):
    """Develop file names and control data storage of flow for a time-step.

    Parameters
    ----------
    simulation = (int) realization number
    time_value = (float) time value of step
    seal_controls = (dict) seal parameters
    mass_flow_co2 = (array) (2D) CO2 mass for realization time step
    mass_flow_brine = (array) (2D) brine mass for realization time step

    Returns
    -------
    None
    """
    # Convert arrays to 2D (if necessary) as new objects for output.
    co2_step = recreate_array_shape(seal_controls, mass_flow_co2)
    brine_step = recreate_array_shape(seal_controls, mass_flow_brine,)

    # Construct file names.
    insert = str(simulation + 1) + "_at_time_" + str(time_value)
    output1 = scfg.BRINE_INTRO + insert + scfg.EXTENSION_CSV
    output2 = scfg.CO2_INTRO + insert + scfg.EXTENSION_CSV

    # Store data in csv files for both flows for time X.
    sim_title = seal_controls['title']
    second_title = " CO2 Flow Through Seal (tonnes)"
    save_csv_results(output2, co2_step, sim_title, second_title)

    second_title = " Brine Flow Through Seal (tonnes)"
    save_csv_results(output1, brine_step, sim_title, second_title)

    # return None


def cache_perm(sim_numbr, grid, seal_controls):
    """Develop file name and store permeability values of a simulation.

    Parameters
    ----------
    sim_numbr = (int) realization number
    grid = (list of class) cells
    seal_controls = (dict) seal parameters

    Returns
    -------
    None
    """
    # set variables for clarity.
    entitle = seal_controls['title']
    total = seal_controls['num_cells']

    # Construct unique file name for saving data.
    current_data = scfg.TOTO_NAME + str(sim_numbr)
    output_name = current_data + scfg.EXTENSION_CSV

    # Create perm array from grid results.
    perm_array = np.zeros(total)
    for numbr, cell in enumerate(grid):
        perm_array[numbr] = cell.permeability

    # Set array per grid approach or freeform.
    perm_array = recreate_array_shape(seal_controls, perm_array)

    # Create header lines and save 2D results to file - with average.
    perm_ave = np.average(perm_array)
    perm_txt = f'{perm_ave:.2f}'
    data_name = "Total Permeability Data (mD) for Each Cell"
    data_name += "   -  Average = " + perm_txt
    save_csv_results(output_name, perm_array, entitle, data_name)

# return None


def write_sim_results(sim_numbr, sim_listing, entitle):
    """Develop file names and control the data storage of a simulation.

    Parameters
    ----------
    sim_numbr = (int) realization number
    sim_listing = (list of str) list of results for a realization
    entitle = (str) title for run

    Returns
    -------
    None
    """
    # Printout control - for long runs, set this control to False (above).
    if PERMIT:
        # Construct unique file name for writing a realization.
        output_name = (scfg.RESULTS_FILE + str(sim_numbr + 1)
                       + scfg.EXTENSION_CSV)

        # Construct full path for new file.
        sub_path, destination = sfile.get_path_name(scfg.OUTPUT_DIRECTORY,
                                                    output_name,
                                                    scfg.EXTENSION_CSV)

        # Define header strings.
        overall_title = "  Run:  " + entitle
        fieldnames = [" Sim. No.", " Time (yrs)", " CO2 Flow (tonnes)",
                      " Brine Flow (tonnes)"]

        # Write data with 2 header lines.
        try:
            with open(destination, "w", newline='', encoding="utf8") \
                    as csv_file:
                writer = csv.writer(csv_file, delimiter=',')
                writer.writerow([overall_title])
                writer.writerow(fieldnames)
                for line in sim_listing:
                    writer.writerow(line)
        except OSError as err:
            sfile.io_snag(err, sub_path, output_name)

    # return None


def write_model_parameters(zum, seal_controls):
    """Print time model info to file named "zum".

    Parameters
    ----------
    zum = (str) file label
    seal_controls = (dict) seal parameters

    Returns
    -------
    None
    """
    # Print time stamp for printout.
    print("\n  SEAL FLUX PROGRAM RUN", file=zum)

    # Print run details.
    print(BLK_IN + "Model Parameters", file=zum)

    print(LIN_IN + f'Analysis Title = {seal_controls.get("title"):<}',
          file=zum)
    if seal_controls.get('use_si_units'):
        unit_type = "METRIC UNITS"
    else:
        unit_type = "ENGLISH UNITS"
    print(LIN_IN + 'Output Units               = '
          f'{unit_type:<}', file=zum)
    print(LIN_IN + 'Start Time of Analysis     = '
          f'{seal_controls.get("start_time"):<.2f} years', file=zum)
    print(LIN_IN + 'End Time of Analysis       = '
          f'{seal_controls.get("end_time"):<.2f} years', file=zum)
    print(LIN_IN + 'Time Steps for Run         = '
          f'{seal_controls.get("time_points"):<5}', file=zum)
    print(LIN_IN + 'Number of Realizations     = '
          f'{seal_controls.get("realizations"):<5}', file=zum)

    if seal_controls.get('time_input'):
        print(LIN_IN + "Time Steps Input From File.", file=zum)
    else:
        print(LIN_IN + "Uniform Time Steps Assumed.", file=zum)

    print(LIN_IN + 'Output Directory = '
          f'{seal_controls.get("output_directory"):<}', file=zum)

    # return None


def write_description_parameters(zum, seal_controls, uts):
    """Print description and grid data to file named "zum".

    Parameters
    ----------
    zum = (str) file label
    seal_controls = (dict) seal parameters
    uts = (list - str) = unit abbreviation

    Returns
    -------
    None
    """
    # Get numeric values for printout.
    temperature = seal_controls.get('temperature')
    height = seal_controls.get('cell_height')
    width = seal_controls.get('cell_width')
    area = seal_controls.get('area')
    inject_co2 = seal_controls.get('total_inject')

    # Recast units for English printout, if used.
    if not seal_controls['use_si_units']:
        # temperature is not a factor!
        temperature = sun.celsius_to_fahrenheit(temperature)
        height *= sun.meters_to_feet()
        width *= sun.meters_to_feet()
        area *= sun.meters2_to_ft2()
        inject_co2 *= sun.tonne_to_ton

    # ----------------------------------------
    # Continue from last group.
    print(LIN_IN + 'Total Injected CO2         = '
          f'{inject_co2:.1e} ' + uts[6], file=zum)

    # Start new group.
    print(BLK_IN + "Seal Description", file=zum)

    print(LIN_IN + 'Number of Cells in Grid    = '
          f'{seal_controls.get("num_cells"):<5}', file=zum)
    if seal_controls.get('area_approach'):
        print(LIN_IN + 'Area of a Cells are Taken from File')
    else:
        print(LIN_IN + 'Area of a Cells in a Grid  = '
              f'{area:.0f} ' + uts[5], file=zum)

    print(LIN_IN + 'Reference Temperature      = '
          f'{temperature:.0f} ' + uts[7], file=zum)
    print(LIN_IN + 'Reference Salinity         = '
          f'{seal_controls.get("salinity"):.0f} ppm', file=zum)

    # ----------------------------------------
    print(BLK_IN + "Grid Details", file=zum)

    if seal_controls['grid_approach']:
        print(LIN_IN + "Grid Approach: Cells Assumed in Box Layout.",
              file=zum)
        print(OTH_IN + 'Rows in Grid               = '
              f'{seal_controls.get("grid_rows"):<5}', file=zum)
        print(OTH_IN + 'Columns in Grid            = '
              f'{seal_controls.get("grid_cols"):<5}', file=zum)
        print(OTH_IN + 'Average Cell Height        = '
              f'{height:.0f} ' + uts[0], file=zum)
        print(OTH_IN + 'Average Cell Width         = '
              f'{width:.0f} ' + uts[0], file=zum)
    else:
        print(LIN_IN + "Grid Approach: Cells Assumed in Free Form Layout.",
              file=zum)

    # return None


def write_fluid_parameters(zum, seal_controls, uts):
    """Print fluid data to file named "zum".

    Parameters
    ----------
    zum = (str) file label
    seal_controls = (dict) seal parameters
    uts = (list - str) = unit abbreviation

    Returns
    -------
    None
    """
    # Get numeric values for print out.
    co2_density = seal_controls.get('co2_density')
    co2_viscosity = seal_controls.get('co2_viscosity')
    br_density = seal_controls.get('brine_density')
    br_viscosity = seal_controls.get('brine_viscosity')
    solubility = seal_controls.get('co2_solubility')

    # Recast units to English units, if desired.
    if not seal_controls['use_si_units']:
        co2_density *= sun.kilo_convert_density()
        co2_viscosity *= sun.pa_s_to_centipoise()
        br_density *= sun.kilo_convert_density()
        br_viscosity *= sun.pa_s_to_centipoise()

    # ----------------------------------------
    print(BLK_IN + "Fluid Parameters", file=zum)

    # Print interpolation control option.
    if seal_controls.get('interpolate_approach'):
        print(LIN_IN + "Fluid Properties Interpolated.", file=zum)
    else:
        print(LIN_IN + "Fluid Properties from Input Values.", file=zum)

    # Print calculated fluid parameters.
    print(LIN_IN + 'CO2 Density                = '
          f'{co2_density:.0f}' + uts[2], file=zum)
    print(LIN_IN + 'CO2 Viscosity              = '
          f'{co2_viscosity:.4e} ' + uts[3], file=zum)
    print(LIN_IN + 'Brine Density              = '
          f'{br_density:.0f} ' + uts[2], file=zum)
    print(LIN_IN + 'Brine Viscosity            = '
          f'{br_viscosity:.4e} ' + uts[3], file=zum)
    print(LIN_IN + 'CO2 Solubility             = '
          f'{solubility:.4e} ' + uts[4], file=zum)

    # return None


def write_depth_parameters(zum, seal_controls, uts):
    """Print fluid data to file named "zum".

    Parameters
    ----------
    zum = (str) file label
    seal_controls = (dict) seal parameters
    uts = (list - str) = unit abbreviation

    Returns
    -------
    None
    """
    # Get numeric values for print out.
    base_depth = seal_controls.get('ave_base_depth')
    base_pressure = seal_controls.get('ave_base_pressure')
    static_pressure = seal_controls.get('static_pressure')
    static_depth = seal_controls.get('static_depth')

    # Recast pressure units to MPa.
    if seal_controls['use_si_units']:
        static_pressure *= sun.pa_to_mpa()
        base_pressure *= sun.pa_to_mpa()
    else:
        # Recast all SI units to English units, if desired.
        static_pressure *= sun.pa_to_psi()
        static_depth *= sun.meters_to_feet()
        base_depth *= sun.meters_to_feet()
        base_pressure *= sun.mpa_to_psi()

    # ----------------------------------------
    print(BLK_IN + "At-Depth Parameters", file=zum)

    # Print base parameters.
    print(LIN_IN + 'Ave. Depth to Seal Base    = '
          f'{base_depth:.0f} ' + uts[0], file=zum)
    print(LIN_IN + 'Ave. Seal Base Pressure    = '
          f'{base_pressure:.4e} ' + uts[1], file=zum)

    # Print stress estimate parameters.
    print(LIN_IN + 'Reference Depth            = '
          f'{static_depth:.0f} ' + uts[0], file=zum)
    print(LIN_IN + 'Reference Pressure         = '
          f'{static_pressure:.4e} ' + uts[1], file=zum)

    # return None


def write_control_values(zum, seal_controls):
    """Print controls (boolean parameters) to file named "zum".

    Parameters
    ----------
    zum = (str) file label
    seal_controls = (dict) seal parameters
    uts = (str) unit abbreviation

    Returns
    -------
    None
    """
    # Controls
    print(BLK_IN + "Seal Approach Parameters", file=zum)

    # Depth Control.
    if seal_controls.get('depth_approach'):
        print(LIN_IN + "Depth Option: Depth to Repository from File.",
              file=zum)
    else:
        print(LIN_IN
              + "Depth Option: All Cells use Repository Depth as Basis.",
              file=zum)

    # Area Control.
    if seal_controls.get('area_approach'):
        print(LIN_IN + "Area Option: Areas Input from File.", file=zum)
    else:
        print(LIN_IN + "Area Option: All Cells have a Uniform Area.",
              file=zum)

    # Layout Control.
    if seal_controls.get('layout_approach'):
        print(LIN_IN
              + "Layout Option: Coordinates Input from File.", file=zum)
    else:
        print(LIN_IN + "Layout Option: Simple Grid Layout Assumed.", file=zum)

    # return None


def write_other_control_values(zum, seal_controls):
    """Print controls (boolean parameters) to file named "zum".

    Parameters
    ----------
    zum = (str) file label
    seal_controls = (dict) seal parameters
    uts = (str) unit abbreviation

    Returns
    -------
    None
    """
    # Print boundary control.
    if seal_controls.get('upper_approach'):
        print(LIN_IN
              + "Top Boundary Option: Input Boundary Conditions from File.",
              file=zum)
    else:
        print(LIN_IN
              + "Top Boundary Option: Used Static Pressure for Analysis.",
              file=zum)

    # Print active control.
    if seal_controls.get('active_approach'):
        print(LIN_IN + "Active Cell Option: Input from File.", file=zum)
    else:
        print(LIN_IN + "Active Cell Option: All Cells are Active.", file=zum)

    # Print influence control.
    if seal_controls.get('initialize'):
        print(LIN_IN + "Influence Array Initialized at Start From File.",
              file=zum)
    else:
        print(LIN_IN + "Influence Array Not Initialized.", file=zum)

    # return None


def write_permeability(zum, seal_controls):
    """Print permeability data to file named "zum".

    Parameters
    ----------
    zum = (str) file label
    seal_controls = (dict) seal parameters

    Returns
    -------
    None

    Notes
    -----
    1. All units are in microdarcys (mD) for both Metric and English/US output
    """
    # ----------------------------------------
    print(BLK_IN + "Equivalent Permeability", file=zum)
    print(LIN_IN + "Units are in microdarcys", file=zum)

    # print permeability variability parameters as prescribed.
    if seal_controls.get('perm_input_approach'):
        print(LIN_IN + "Permeability Option: Permeabilities Input from File.",
              file=zum)
    elif seal_controls.get('fracture_approach'):
        print(LIN_IN
              + "Fracture Option: Permeability Computed Using "
              + "Random Fractures.",
              file=zum)
    elif not seal_controls.get('vary_perm_choice'):
        print(LIN_IN + "Uniform Permeability - Set to Mean Value.", file=zum)
        print(LIN_IN + 'Mean Total Permeability    = '
              f'{seal_controls.get("perm_mean"):<.3e} mD', file=zum)
    else:
        print(LIN_IN + 'Computed with Censored Log-normal Distribution.',
              file=zum)
        print(LIN_IN + 'Mean Total Permeability    = '
              f'{seal_controls.get("perm_mean"):<.3e} mD', file=zum)
        print(LIN_IN + 'Std. Dev. Permeability     = '
              f'{seal_controls.get("perm_std"):<.4e} mD', file=zum)
        print(LIN_IN + 'Minimum Total Permeability = '
              f'{seal_controls.get("perm_min"):.3e} mD', file=zum)
        print(LIN_IN + 'Maximum Total Permeability = '
              f'{seal_controls.get("perm_max"):.3e} mD',  file=zum)

    # ----------------------------------------
    # Print heterogeneity response.
    if (seal_controls.get('perm_input_approach') and
            seal_controls.get('heterogeneity_approach')):
        print(LIN_IN + "Heterogeneity Option: Included.", file=zum)
        print(LIN_IN + 'Heterogeneity Factor       = '
              f'{seal_controls.get("perm_heter_factor"):.3f}', file=zum)
    else:
        print(LIN_IN + "Heterogeneity Option is NOT Included.", file=zum)

    # return None


def write_relative_perm(zum, seal_controls, uts):
    """Print relative permeability controls to a file named "zum".

    Parameters
    ----------
    zum = (str) file label
    seal_controls = (dict) seal parameters
    uts = (str) = unit abbreviation

    Returns
    -------
    None
    """
    # Define stress values in MPa.
    cap_pressure = seal_controls.get('entry_pressure')
    if seal_controls['use_si_units']:
        cap_pressure *= sun.pa_to_mpa()
    else:
        cap_pressure *= sun.psi_to_mpa()

    max_press = 0.0
    # Define pressure for LET model.
    if 'LET' in seal_controls.get('relative_model'):
        max_press = seal_controls.get('max_capillary')
        if seal_controls['use_si_units']:
            max_press *= sun.pa_to_mpa()
        else:
            max_press *= sun.pa_to_psi

    # ----------------------------------------
    # Print limits on relative permeability.
    print(BLK_IN + "Relative Permeability Flow Limits", file=zum)
    print(LIN_IN + 'Residual Brine Saturation  = '
                   f'{seal_controls.get("resid_brine"):.2f}', file=zum)
    print(LIN_IN + 'Residual CO2 Saturation    = '
                   f'{seal_controls.get("resid_co2"):.2f}', file=zum)
    print(LIN_IN + 'Nonwetting/Wetting Ratio   = '
                   f'{seal_controls.get("perm_ratio"):.2f}', file=zum)
    print(LIN_IN + 'Ave. Threshold Pressure    = '
                   f'{cap_pressure:.4e} ' + uts[1], file=zum)

    # ----------------------------------------
    print(BLK_IN + "Relative Permeability Model", file=zum)
    # Print relative permeability model type.
    m_type = seal_controls.get('relative_model')
    if "LET" in m_type:
        print(LIN_IN + "Relative Permeability = L-E-T Model",
              file=zum)
    else:
        print(LIN_IN + "Relative Permeability = " +
              "Brooks-Corey Model", file=zum)

    # ----------------------------------------
    # LET & BC relative permeability parameters.
    print(BLK_IN + "Relative Permeability Parameters", file=zum)

    if 'LET' in seal_controls.get('relative_model'):
        print(LIN_IN + '"L" Wetting Parameter      = '
              f'{seal_controls.get("l_wetting"):.2f}', file=zum)
        print(LIN_IN + '"E" Wetting Parameter      = '
              f'{seal_controls.get("e_wetting"):.2f}', file=zum)
        print(LIN_IN + '"T" Wetting Parameter      = '
              f'{seal_controls.get("t_wetting"):.2f}', file=zum)
        print(LIN_IN + '"L" Nonwetting Parameter   = '
              f'{seal_controls.get("l_nonwet"):.2f}', file=zum)
        print(LIN_IN + '"E" Nonwetting Parameter   = '
              f'{seal_controls.get("e_nonwet"):.2f}', file=zum)
        print(LIN_IN + '"T" Nonwetting Parameter   = '
              f'{seal_controls.get("t_nonwet"):.2f}', file=zum)
    else:
        print(LIN_IN + 'Lambda Parameter for BC    = '
              f'{seal_controls.get("zeta"):.2f}', file=zum)

    # ----------------------------------------
    # LET model - capillary parameters.
    if 'LET' in seal_controls.get('relative_model'):
        print("\n  > Capillary Pressure Parameters", file=zum)
        print(LIN_IN + '"L" Capillary Parameter    = '
              f'{seal_controls.get("l_capillary"):.2f}', file=zum)
        print(LIN_IN + '"E" Capillary Parameter    = '
              f'{seal_controls.get("e_capillary"):.2f}', file=zum)
        print(LIN_IN + '"T" Capillary Parameter    = '
              f'{seal_controls.get("t_capillary"):.2f}', file=zum)
        print(LIN_IN + 'Maximum Capillary Pressure = '
              f'{max_press:.4e} ' + uts[1], file=zum)

    # return None


def write_reactivity(zum, seal_controls):
    """Print reactivity parameters to file named "zum".

    Parameters
    ----------
    zum = (str) file label
    seal_controls = (dict) seal parameters
    uts = (str) unit abbreviation

    Returns
    -------
    None
    """
    # Print reactivity and thickness parameters.
    print(BLK_IN + "Reactivity & Alteration Parameters", file=zum)
    time_dependent_model = seal_controls.get('model')
    print(LIN_IN + 'Time-Dependent Model       = '
          f'{time_dependent_model:1d}', file=zum)

    # Material parameters - if model > 0.
    if time_dependent_model > 0:
        print(LIN_IN + 'Time Rate of Change        = '
              f'{seal_controls.get("rate_effect"):.4f}', file=zum)
        print(LIN_IN + 'Total Effect of Change     = '
              f'{seal_controls.get("total_effect"):.4f}', file=zum)

        # For model=2 only.
        if time_dependent_model > 1:
            print(LIN_IN + 'Seal Reactivity            = '
                  f'{seal_controls.get("reactivity"):.2f}', file=zum)
            print(LIN_IN + 'Total Clay Content         = '
                  f'{seal_controls.get("clay_content"):.2f}%', file=zum)
            print(LIN_IN + 'Total Carbonate            = '
                  f'{seal_controls.get("carbonate_content"):.2f}%', file=zum)
            print(LIN_IN + 'Clay Type                  = '
                  f'{seal_controls.get("clay_type"):<}', file=zum)

    # return None


def write_thickness(zum, seal_controls, uts):
    """Print thickness parameters to file named "zum".

    Parameters
    ----------
    zum = (str) file label
    seal_controls = (dict) seal parameters
    uts = (str) unit abbreviation

    Returns
    -------
    None
    """
    # Get numeric values for print out.
    ave_thick = seal_controls.get('thickness_ave')
    std_thick = seal_controls.get('thickness_std')
    min_thick = seal_controls.get('thickness_min')
    max_thick = seal_controls.get('thickness_max')

    # Recast units for English output.
    if not seal_controls['use_si_units']:
        ave_thick *= sun.meters_to_feet()
        std_thick *= sun.meters_to_feet()
        min_thick *= sun.meters_to_feet()
        max_thick *= sun.meters_to_feet()

    # ----------------------------------------
    # Print thickness variables - tailor response to selection.
    print(BLK_IN + "Seal Thickness", file=zum)
    if seal_controls.get('thickness_approach'):
        print(LIN_IN + "Thickness Option: Thickness Values from File",
              file=zum)
    else:
        print(LIN_IN + "Thickness Option: Compute Thickness", file=zum)
        print(LIN_IN + "Method: Truncated Normal Distribution", file=zum)
        print(LIN_IN + 'Mean Thickness             = '
              f'{ave_thick:.1f} ' + uts[0], file=zum)
        print(LIN_IN + 'Thickness Std. Dev.        = '
              f'{std_thick:.1f} ' + uts[0], file=zum)
        print(LIN_IN + 'Minimum Thickness          = '
              f'{min_thick:.1f} ' + uts[0], file=zum)
        print(LIN_IN + 'Maximum Thickness          = '
              f'{max_thick:.1f} ' + uts[0], file=zum)

    # return None


def write_co2_results(zum, sim_flux_list, uts, seal_controls, si_units):
    """Compute and print results of simulations to file "zum".

    Parameters
    ----------
    zum = (str) file label
    sim_flux_list = (list) CO2 results at the end of simulations
    uts = (array str) unit abbreviation
    seal_controls = (dict) seal parameters
    si_units = (bool) metric or English units

    Returns
    -------
    ave_results = (float) flux as percentage of total
    """
    # Convert list to array.
    simulation_number = len(sim_flux_list)
    flux_data = np.asarray(sim_flux_list)

    # compute statistics of the array.
    ave_flux_data = np.mean(flux_data)
    std_dev = np.std(flux_data)
    ave_percent = ave_flux_data / seal_controls.get('total_inject') * 100.00
    min_flux_data = np.min(flux_data)
    max_flux_data = np.max(flux_data)
    argument_min = np.argmin(flux_data) + 1
    argument_max = np.argmax(flux_data) + 1

    # Convert to English units as desired.
    if not si_units:
        ave_flux_data *= sun.tonne_to_ton()
        std_dev *= sun.tonne_to_ton()
        min_flux_data *= sun.tonne_to_ton()
        max_flux_data *= sun.tonne_to_ton()

    # ----------------------------------------
    # print statistics.
    print("\n  " + "*" * 60, file=zum)

    print(BLK_IN + 'Summary of Results', file=zum)
    print(LIN_IN + 'Number of Simulations               = '
          f'{simulation_number:,}', file=zum)
    print(LIN_IN + 'Average CO2 Flux of All Simulations = '
          f'{ave_flux_data:.3e} ' + uts[6], file=zum)
    print(LIN_IN + 'Percentage of Total Injected CO2    = '
                   f'{ave_percent:.3f} %', file=zum)

    # Error check - flux must be small!
    if ave_percent > MAX_FLUX_LIMIT:
        msg = 'ERROR - Flux Value Exceeds Code Assumptions!'
        print(BLK_IN + msg, flush=True, file=zum)

    else:
        # If no error, print remaining statistics.
        print(LIN_IN + 'Standard Deviation of CO2 flux      = '
              f'{std_dev:.4e} ' + uts[6], file=zum)

        print(LIN_IN + 'Minimum CO2 Flux of All Simulations = '
              f'{min_flux_data:<.3e} ' + uts[6], file=zum)
        print(OTH_IN + f'- at Index =      {argument_min:}', file=zum)

        print(LIN_IN + 'Maximum CO2 Flux of All Simulations = '
              f'{max_flux_data:.3e} ' + uts[6], file=zum)
        print(OTH_IN + f'- at Index =      {argument_max:}', file=zum)

    return ave_percent


def write_last(zum, seal_controls, seal_version):
    """Write version & current time value to summary file named "zum" and end.

    Parameters
    ----------
    zum = (str) file label
    seal_controls = (dict) dictionary of fault parameters
    seal_version = (str) seal_flux version

    Returns
    -------
    N/A
    """
    # Print version number.
    print('\n  ', end='', file=zum)
    print('*' * 60, file=zum)
    print(BLK_IN + "Seal Flux Version Number = " + seal_version, file=zum)

    # Print date.
    now = datetime.now()
    clump = HED_IN + "Simulation Date =  %a %d %b %Y, Time: %H:%M"
    stamp = now.strftime(clump)
    print(stamp, file=zum)

    # Print run time.
    print(HED_IN + 'Computation Time (h.mm.ss.xxxx) = '
          f'{seal_controls.get("elapsed_time")}', file=zum)

    # Print last line --> End.
    print('\n  End', file=zum)

    # return None


def write_summary(seal_controls, sim_flux_list, version_num, uts):
    """Control the printout of summary data to file.

    Parameters
    ----------
    seal_controls = (dict) seal parameters
    sim_flux_list = (list) CO2 flux at end of simulations
    version_num = (str) Seal_Flux version number
    uts[] =  (list - str) units for parameters

    Returns
    -------
    None

    Notes
    -----
    SEAL_SUM_FILE = Defined destination file name in output directory.
    """
    # Construct full path to results directory and summary file.
    #  -- Check for extension "txt".
    sub_path, destination = sfile.get_path_name(scfg.OUTPUT_DIRECTORY,
                                                scfg.SEAL_SUM_FILE,
                                                scfg.EXTENSION_TXT)

    # Write information to file on seal run.
    try:
        with open(destination, 'w', encoding="utf8") as zum:
            write_model_parameters(zum, seal_controls)
            write_description_parameters(zum, seal_controls, uts)
            write_fluid_parameters(zum, seal_controls, uts)
            write_depth_parameters(zum, seal_controls, uts)
            write_control_values(zum, seal_controls)
            write_other_control_values(zum, seal_controls)
            write_permeability(zum, seal_controls)
            write_relative_perm(zum, seal_controls, uts)
            write_reactivity(zum, seal_controls)
            write_thickness(zum, seal_controls, uts)
            ave_results = write_co2_results(zum, sim_flux_list, uts,
                                            seal_controls,
                                            seal_controls['use_si_units'])
            write_last(zum, seal_controls, version_num)

    except OSError as err:
        sfile.io_snag(err, sub_path, scfg.SEAL_SUM_FILE)

    # Error in assumptions - terminate run.
    if ave_results > MAX_FLUX_LIMIT:
        # Record error.
        msg = 'Resulting Flux Values Exceeds Code Assumptions!'
        logging.error(msg)

        # End program.
        sfile.terminate_code(1, True)

    # return None


#
# -----------------------------------------------------------------------------
# - End of module
