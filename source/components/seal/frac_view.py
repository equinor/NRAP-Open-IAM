#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions for writing fracture summary file and fracture debug files.

Author: Ernest N. Lindner
Date: 08/18/2022

Module Name
    frac_view

Contents (16)
    write_frac_run_parameters(cast, frac_controls)
    write_frac_controls(cast, frac_controls)
    write_stochastic_parameters(cast, frac_controls, uts, si_units)
    write_more_stochastic_parameters(cast, frac_controls, uts, si_units)
    write_pressure_parameters(cast, frac_controls, uts, si_units)
    write_stochastic_options(cast, frac_controls, uts, si_units)
    write_defined_fractures(cast, frac_controls)
    write_rock_data(cast, frac_controls, uts)
    write_translated_controls(cast, frac_controls)
    write_time(cast, elapsed_time)
    ----
    write_completed(cast, frac_controls)
    write_frac_summary(frac_controls, alive, x_method, uts, si_units)

    DEBUG:
    write_overall_results(cast, frac_controls, data_list)
    write_list(cast, header, data_list)
    write_sections(cast, frac_controls, data_list)
    write_data_results(frac_controls, data_list, sim_numbr)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
from datetime import datetime       # For time stamp
import numpy as np                  # For array operations

import seal_config as scfg           # IO directory and file names
import seal_file as sfile           # Error in file operations
import seal_units as sun            # Unit conversion

# Constants for file handling
OUTPUT_COLS = 8                     # Columns for output
CONFIG = 11.3                       # Output format

# Constants for writing
LINER_2 = ("  " + "*" * 60)         # Line across page
INTRO = "\n  >>> "                  # Line start
PARA = "\n  > "                     # Indent at paragraph start
IDNT = "  --> "                     # Indent at item start
SPEC = "      "                     # General indent

# Debug constants
ECHO = False                         # Print averages for each block


def write_frac_run_parameters(cast, frac_controls):
    """Print basic model info to file named "cast".

    Parameters
    ----------
    cast = (str) file label
    frac_controls = (dict) dictionary of fracture parameters

    Returns
    -------
    None
    """
    # Print header for printout.
    print("\n  FRACTURE GENERATION RUN", file=cast)

    # Print run details.
    print(PARA + "Model Parameters", file=cast)
    print(IDNT + f'Analysis Title = {frac_controls.get("title"):<}', file=cast)

    # return None


def write_frac_controls(cast, frac_controls):
    """Print region data and controls to file named "cast".

    Parameters
    ----------
    cast = (str) file label
    frac_controls = (dict) dictionary of fracture parameters

    Returns
    -------
    None
    """
    # Print header for input controls used.
    print(PARA + "Controls:", file=cast)

    # Print input random-fracture generation control.
    if frac_controls.get('random_approach'):
        print(IDNT + "Random Option: Random Fractures Generated", file=cast)
    else:
        print(IDNT + "Random Option: Random Fractures NOT Generated",
              file=cast)

    # Print user-fractures control.
    if frac_controls.get('user_approach'):
        print(IDNT + "Input Option: Fractures Read from File.",
              file=cast)
    else:
        print(IDNT + "Input Option: Fractures NOT Read from File.",
              file=cast)

    # Print aperture correlation control.
    if frac_controls.get('correlate_approach'):
        print(IDNT + "Aperture Option: Apertures from Length Correlation.",
              file=cast)
    else:
        print(IDNT + "Aperture Option: Apertures Generated Stochastically.",
              file=cast)

    # Print pressure-aperture control.
    if frac_controls.get('pressure_approach'):
        print(IDNT + "Pressure-Aperture Option: Correction Applied.",
              file=cast)
    else:
        print(IDNT + "Pressure-Aperture Option: Not Used.",
              file=cast)

    # Print header for connectivity factor.
    print(PARA + "General:", file=cast)
    # Print vertical connectivity parameter.
    print(IDNT + 'Fracture Connectivity Factor     = '
          f'{frac_controls.get("connect_factor"):.4} ', file=cast)

    # return None


def write_stochastic_parameters(cast, frac_controls, uts, si_units):
    """Print first part stochastic fracture data to file named "cast".

    Parameters
    ----------
    cast = (str) file label
    frac_controls = (dict) dictionary of fracture parameters
    uts[] = (list of str) units for parameters
    si_units = (bool) control to convert SI to English/US units

    Returns
    -------
    None
    """
    # Get numeric density & length values.
    ave_density = frac_controls.get('density_ave')
    min_density = frac_controls.get('density_min')
    max_density = frac_controls.get('density_max')
    ave_length = frac_controls.get('length_ave')
    std_length = frac_controls.get('length_dev')
    min_length = frac_controls.get('length_min')
    max_length = frac_controls.get('length_max')

    # Recast units for English/US output.
    if not si_units:
        ave_density *= sun.per_sqm_to_per_sft()
        min_density *= sun.per_sqm_to_per_sft()
        max_density *= sun.per_sqm_to_per_sft()
        ave_length *= sun.meters_to_feet()
        std_length *= sun.meters_to_feet()
        min_length *= sun.meters_to_feet()
        max_length *= sun.meters_to_feet()

    # --------------------------------------------------
    # Print density parameters.
    print(PARA + 'Density - Triangular Distribution', file=cast)
    print(IDNT + f'Average Fracture Density         = {ave_density:.5e} //'
          + uts[5], file=cast)
    print(IDNT + f'Minimum Fracture Density         = {min_density:.5e} //'
          + uts[5], file=cast)
    print(IDNT + f'Maximum Fracture Density         = {max_density:.5e} //'
          + uts[5], file=cast)

    # Print orientation parameters.
    print(PARA + 'Orientation - von Mises Distribution', file=cast)
    print(IDNT + 'Average Fracture Orientation     = '
          f'{frac_controls.get("orient_mu"):.2f} degrees', file=cast)
    print(IDNT + 'Orientation Spread - 2 Sigma     = '
          f'{frac_controls.get("orient_sigma"):.2f} degrees', file=cast)

    # Print length parameters - depending on function type.
    if frac_controls['length_approach'] == 'POWER':
        print(PARA + "Length - Power-Law Distribution", file=cast)
        print(IDNT + 'Power-Law Exponent               = '
              f'{frac_controls.get("length_eta"):.3f}', file=cast)
    else:
        print(PARA + 'Length - Censored Log-Normal Distribution', file=cast)
        print(IDNT + 'Average Fracture Length          = '
              f'{ave_length:.2f} ' + uts[0], file=cast)
        print(IDNT + 'Standard Deviation in Length     = '
              f'{std_length:.3f} ' + uts[0], file=cast)

    # Print remaining length parameters common to both functions.
    print(IDNT + 'Minimum Fracture Length          = '
          f'{min_length:.1f} ' + uts[0], file=cast)
    print(IDNT + 'Maximum Fracture Length          = '
          f'{max_length:.1f} ' + uts[0], file=cast)

    # return None


def write_more_stochastic_parameters(cast, frac_controls, uts, si_units):
    """Print part 2 stochastic fracture data to file named "cast".

    Parameters
    ----------
    cast = (str) file label
    frac_controls = (dict) dictionary of fracture parameters
    uts[] = (list of str) units for parameters
    si_units = (bool) control to convert SI to English/US units

    Returns
    -------
    None
    """
    # Get numeric values.
    ave_aperture = frac_controls.get('aperture_ave')
    dev_aperture = frac_controls.get('aperture_dev')
    min_aperture = frac_controls.get('aperture_min')
    max_aperture = frac_controls.get('aperture_max')
    entry_stress = frac_controls.get('entry_pressure')

    # Recast units for English/US printout.
    if si_units:
        entry_stress *= sun.pa_to_mpa()
    else:
        ave_aperture *= sun.mm_to_inch()
        dev_aperture *= sun.mm_to_inch()
        min_aperture *= sun.mm_to_inch()
        max_aperture *= sun.mm_to_inch()
        entry_stress *= sun.pa_to_psi()

    # ----------------------------------------
    # Print aperture parameters.
    if frac_controls['correlate_approach']:
        print(PARA + 'Aperture-Length Correlation Factors', file=cast)

        print(IDNT + 'Aperture Correlation - Alpha     = '
              f'{frac_controls.get("aperture_alpha"):.2f} ', file=cast)
        print(IDNT + 'Aperture Correlation - Beta      = '
              f'{frac_controls.get("aperture_beta"):.4f} '
              + uts[8], file=cast)

    else:
        print(PARA + 'Aperture - Censored Log-Normal Distribution',
              file=cast)
        print(IDNT + 'Average Fracture Aperture        = '
              f'{ave_aperture:.4e} ' + uts[8], file=cast)
        print(IDNT + 'Standard Deviation in Aperture   = '
              f'{dev_aperture:.4e} ' + uts[8], file=cast)
        print(IDNT + 'Minimum Fracture Aperture        = '
              f'{min_aperture:.4e} ' + uts[8], file=cast)
        print(IDNT + 'Maximum Fracture Aperture        = '
              f'{max_aperture:.4e} ' + uts[8], file=cast)

    # Threshold printout.
    print(PARA + 'Threshold References', file=cast)
    print(IDNT + 'Reference Threshold Pressure     = '
          f'{entry_stress:.3e} ' + uts[1], file=cast)

    # return None


def write_pressure_parameters(cast, frac_controls, uts, si_units):
    """Print pressure-approach data to file named "cast".

    Parameters
    ----------
    cast = (str) file label
    frac_controls = (dict) dictionary of fracture parameters
    uts[] = (list of str) units for parameters
    si_units = (bool) control to convert SI to English/US units

    Returns
    -------
    None
    """
    # Pressure-aperture correction parameters.
    if frac_controls['pressure_approach']:

        # Define variables for pressure-aperture.
        resid_aperture = frac_controls.get('residual_aperture')
        wide_aperture = frac_controls.get('wide_aperture')
        limit_pressure = frac_controls.get('stress_limit')
        theta_value = frac_controls.get('theta_aperture')

        # Define units for printout.
        if si_units:
            limit_pressure *= sun.pa_to_mpa()
        else:
            limit_pressure *= sun.pa_to_psi()
            resid_aperture *= sun.mm_to_inch()
            wide_aperture *= sun.mm_to_inch()

        # Print parameters.
        print(PARA + "Pressure Aperture Parameters",
              file=cast)
        print(IDNT + 'Residual Aperture                = '
              f'{resid_aperture:.3e} ' + uts[8], file=cast)
        print(IDNT + 'Maximum Aperture                 = '
              f'{wide_aperture:.3e} ' + uts[8], file=cast)
        print(IDNT + 'Stress Limit - Non-Linear Curve  = '
              f'{limit_pressure:.3e} ' + uts[1], file=cast)
        print(IDNT + 'Stiffness History Factor         = '
              f'{theta_value:.3f} ', file=cast)

    # return None


def write_stochastic_options(cast, frac_controls, uts, si_units):
    """Control the printout of summary fracture data to cast file.

    Parameters
    ----------
    cast = (str) file label
    frac_controls = (dict) dictionary of fracture parameters
    uts[] = (list - str) units for parameters
    si_units = (bool) control to convert English to SI units

    Returns
    -------
    None

    Notes
    -----
    1. FRAC_SUM_FILE = Destination file name
    """
    # Printer header for section.
    print("\n" + LINER_2, file=cast)
    print(PARA + "Stochastic Parameters", file=cast)

    if frac_controls.get("random_approach"):
        # Print note that stochastic data were generated.
        print(PARA + 'Stochastic Fractures Were Generated.', file=cast)

        # Print stochastic data.
        write_stochastic_parameters(cast, frac_controls, uts,
                                    si_units)
        write_more_stochastic_parameters(cast, frac_controls, uts,
                                         si_units)
        write_pressure_parameters(cast, frac_controls, uts,
                                  si_units)
    else:
        # No stochastic parameters.
        print(PARA + 'No Stochastic Fractures Were Generated.',
              file=cast)

    # return


def write_defined_fractures(cast, frac_controls):
    """Print user input fracture data to file named "cast".

    Parameters
    ----------
    cast = (str) file label
    frac_controls = (dict) dictionary of fracture parameters

    Returns
    -------
    None
    """
    # Print section header.
    print("\n" + LINER_2, file=cast)
    print(PARA + 'User Defined Fractures', file=cast)

    # Print User Data - if used.
    if frac_controls.get('user_approach'):

        if frac_controls.get('user_fractures') > 0:
            print(IDNT + "User Fractures Were Read from File.", file=cast)
            print(PARA + "User Defined Fractures:", file=cast)
            print(IDNT + 'Source File in Input = '
                  f'{frac_controls.get("input_source")}', file=cast)
            print(IDNT + 'Number of Fracture Lines Input     = '
                  f'{int(frac_controls.get("user_fractures")):6d}', file=cast)
        else:
            print(PARA + "No User Fractures Were Found in File.", file=cast)
    else:
        print(IDNT + "User Fractures Were Not Read from File.", file=cast)

    # return None


def write_rock_data(cast, frac_controls, uts, si_units):
    """Print matrix permeability values to file named "cast".

    Parameters
    ----------
    cast = (str) file label
    frac_controls = (dict) dictionary of fracture parameters
    uts[] = (list of str) units for parameters
    si_units = (bool) control to convert English to SI units

    Returns
    -------
    None
    """
    # Print matrix parameters in microdarcies.
    print("\n" + LINER_2, file=cast)
    print(PARA + "Matrix Permeability", file=cast)
    print(PARA + "Matrix Stochastic Parameters:", file=cast)
    print(IDNT + "Permeability Generated from Censored Normal Distribution.",
          file=cast)
    print(IDNT + 'Average Matrix Permeability        = '
          f'{frac_controls.get("rock_perm_ave"):.4e} mD', file=cast)
    print(IDNT + 'Standard Deviation in Permeability = '
          f'{frac_controls.get("rock_perm_dev"):.4e} mD', file=cast)
    print(IDNT + 'Minimum Matrix Permeability        = '
          f'{frac_controls.get("rock_perm_min"):.4e} mD', file=cast)
    print(IDNT + 'Maximum Matrix Permeability        = '
          f'{frac_controls.get("rock_perm_max"):.4e} mD', file=cast)

    # Print reference parameters.
    print(PARA + "Matrix Reference Parameters", file=cast)
    print(IDNT + 'Ref. Matrix Permeability           = '
          f'{frac_controls.get("ref_matrix_perm"):.4e} mD', file=cast)

    pressure = frac_controls.get('ref_matrix_threshold')
    if si_units:
        pressure *= sun.pa_to_mpa()
    else:
        pressure *= sun.pa_to_psi()
    print(IDNT + 'Ref. Matrix Threshold Pressure     = '
          f'{pressure:.1e} ' + uts[1], file=cast)

    # return None


def write_translated_controls(cast, frac_controls):
    """Print user translated lognormal data to file named "cast".

    Parameters
    ----------
    cast = (str) file label
    frac_controls = (dict) dictionary of fracture parameters

    Returns
    -------
    None
    """
    # Print stochastic parameters of Normal distributions.
    print("\n" + LINER_2, file=cast)
    print(PARA + "Translated Log-Normal Distribution Parameters", file=cast)
    print(SPEC + "- Values Computed Internally by Code.", file=cast)

    # Print length and aperture parameters - if generated.
    if frac_controls.get('random_approach'):

        # Print length parameters if log-normal distribution.
        if frac_controls['length_approach'] == 'LOGNORM':
            print(PARA + 'Length Log-Normal Distribution Parameters:',
                  file=cast)
            print(IDNT + 'Length - Mu Value                  = '
                  f'{frac_controls.get("length_mu"):.6f}', file=cast)
            print(IDNT + 'Length - Sigma Value               = '
                  f'{frac_controls.get("length_scale"):.6f}', file=cast)

        # Print aperture parameters - if used.
        if not frac_controls.get("correlate_approach"):
            print(PARA + 'Aperture Log-Normal Distribution Parameters:',
                  file=cast)
            print(IDNT + 'Aperture - Mu Value                = '
                  f'{frac_controls.get("aperture_mu"):.6f}', file=cast)
            print(IDNT + 'Aperture - Sigma Value             = '
                  f'{frac_controls.get("aperture_scale"):.6f}', file=cast)

    # Print matrix permeability parameters.
    print(PARA + 'Matrix Log-Normal Distribution Parameters:',
          file=cast)
    print(IDNT + 'Matrix Permeability - Mu Value     = '
          f'{frac_controls.get("rock_perm_mu"):.6f}', file=cast)
    print(IDNT + 'Matrix Permeability - Sigma Value  = '
          f'{frac_controls.get("rock_perm_scale"):.6f}', file=cast)

    # return None


def write_time(cast, elapsed_time):
    """Write run time in simple format.

    Parameters
    ----------
    cast = (str) file label for output - summary file
    elapsed_time = (float) delta time value for run

    Returns
    -------
    None
    """
    # compute time values for delta time in seconds.
    seconds = int(elapsed_time.total_seconds())
    hours, remainder = divmod(seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    milli_secs = int(elapsed_time.microseconds / 1000)

    # Write time in design format.
    print(INTRO + 'Fracture Execution Time = '
          f'{int(hours):02}:{int(minutes):02}:'
          f'{int(secs):02}.{int(milli_secs):03}  (hh:mm:ss:ms)',
          file=cast, flush=True)

    # return None


def write_completed(cast, frac_controls):
    """Print end parameters info to file named "cast".

    Parameters
    ----------
    cast = (str) file label
    frac_controls = (dict) dictionary of fracture parameters

    Returns
    -------
    None
    """
    # Print separator.
    print("\n" + LINER_2, file=cast)

    # Print version number.
    print(PARA + 'Seal Flux / Frac Version Number = '
          f'{frac_controls.get("version"):<}', file=cast)

    # Print date.
    now = datetime.now()
    clump = "  > " + "Run Date: %a %d %b %Y - Time: %H:%M"
    stamp = now.strftime(clump)
    print(stamp, file=cast)

    # Closing line.
    print("\n  End", file=cast)

    # return None


def write_frac_summary(frac_controls, alive, x_method, uts, si_units):
    """Control the printout of summary fracture data to cast file.

    Parameters
    ----------
    frac_controls = (dict) dictionary of fracture parameters
    alive = (bool) stand-alone run
    x_method = (bool) condition if fracture method was used
    uts[] = (list - str) units for parameters
    si_units = (bool) control to convert English to SI units

    Returns
    -------
    None

    Notes
    -----
    1. FRAC_SUM_FILE = Destination file name
    """
    if x_method:
        # Construct full path to results directory and summary file.
        sub_path, destination = sfile.get_path_name(scfg.OUTPUT_DIRECTORY,
                                                    scfg.FRAC_SUM_FILE,
                                                    scfg.EXTENSION_TXT)

        # Write information to file on fracture analysis.
        try:
            with open(destination, 'w', encoding="utf8") as cast:

                # Print run information.
                write_frac_run_parameters(cast, frac_controls)
                write_frac_controls(cast, frac_controls)

                # Write stochastic data.
                write_stochastic_options(cast, frac_controls, uts, si_units)

                # Write user and matrix parameters.
                write_defined_fractures(cast, frac_controls)
                write_rock_data(cast, frac_controls, uts, si_units)

                # Write derived parameters and time.
                write_translated_controls(cast, frac_controls)
                write_completed(cast, frac_controls)

        except OSError as err:
            sfile.io_snag(err, sub_path, scfg.FRAC_SUM_FILE)

        if alive:
            sfile.echo_status(alive, "CREATED FRACTURE SUMMARY FILE.")

    # return None


# **** DEBUGGING ROUTINES *******


def write_overall_results(cast, frac_controls, data_list):
    """For debugging, print averages for each block to the file "cast.

    Parameters
    ----------
    cast = (str) file label
    frac_controls = (dict) dictionary of fracture parameters
    data_list = (list) list results -> for each group:
        =[0] for random fractures
        =[1] for user fractures
        =[2] for matrix

    Returns
    -------
    None

    Notes
    -----
    1. See function "compute_ave_terms" for key to data_list
    """
    # Define loop control for each case (random and/or user fractures).
    ender = 2
    if data_list[0]:
        starter = 0
        if data_list[1]:
            ender = 2
        else:
            ender = 1
    else:
        starter = 1
        if not data_list[1]:
            starter = 2  # Will skip loop if start = end

    # Define header for each case.
    for numbr in range(starter, ender):
        # Control parameters for loop.
        if numbr == 0:
            header = PARA + "Average Parameters for Random Fracture per Grid"
            frac_info = np.asarray(data_list[0])
        else:
            print(file=cast)
            header = "n  > Average Parameters for User Fracture per Grid"
            frac_info = np.asarray(data_list[1])

        # Print block data for each case (reverse loop).
        print(header, file=cast)
        for indx in range(frac_controls['num_cells']):
            # Define block number.
            bloc_id = f'[ {indx} ]'
            bloc_id = PARA + "Average Data for Cell " + bloc_id
            print(bloc_id, file=cast)
            print(IDNT + 'Number of Fractures in Element    = '
                  f'{frac_info[indx, 5]:6d}', file=cast)
            print(IDNT + 'Average Aperture - All            = '
                  f'{frac_info[indx, 0]:.3g} mm', file=cast)
            print(IDNT + 'Average Length - All              = '
                  f'{frac_info[indx, 1]:.1f} m', file=cast)
            print(IDNT + 'Average Strike - All              = '
                  f'{frac_info[indx, 2]:.1f} deg', file=cast)
            print(IDNT + 'Average Transmissivity - All      = '
                  f'{frac_info[indx, 3]:.3g} m^4', file=cast)
            print(IDNT + 'Average Threshold - All           = '
                  f'{frac_info[indx, 4]:.3g} Pa', file=cast)
            print(IDNT + 'Permeability                      = '
                  f'{frac_info[indx, 6]:.3g} mD', file=cast)

    # return None


def write_list(cast, header, data_list):
    """Print list in arbitrary format to "cast" file.

    Parameters
    ----------
    cast = (str) file label for output
    header = (str) title of data
    data_list = (list) list for output to summary file

    Returns
    -------
    None

    Notes
    -----
    1. NumPy print functions will not print all elements or larger
        arrays vertically or laterally. Therefore, this compromise
        was developed for pretty printing.
    """
    # Compute the number of lines and the "remainder" columns on last line.
    columns = OUTPUT_COLS
    n_lines, remainder = divmod(len(data_list), columns)

    # Setup format and controls.
    pattern = f'  {{:{CONFIG}e}}'
    typ_line = '\n'.join(pattern * columns for _ in range(n_lines))
    last_line = pattern * remainder

    # Print data with formatting.
    print(IDNT + header, file=cast)
    print(typ_line.format(*data_list), file=cast)
    print(last_line.format(*data_list[n_lines*columns:]), file=cast)

    # return None


def write_sections(cast, frac_controls, data_list):
    """Print permeability and threshold of each grid to "cast" file.

    Parameters
    ----------
    cast = (str) file label for output - summary file
    frac_controls = (dict) dictionary of fracture parameters
    data_list = (list) list results -> for each group:
        =[0] for random fractures
        =[1] for user fractures
        =[2] for matrix

    Returns
    -------
    None
    """
    # Setup for lists.
    list_perm = []
    list_threshold = []
    total = frac_controls['num_cells']

    # Write random fracture data, if fracture data present.
    if data_list[0]:
        print(PARA + "Grid Parameters for Random Fractures:", file=cast)

        # Create sub-lists for printing threshold and permeability.
        for indx in range(0, total):
            list_threshold.append(data_list[0][indx][4])
            list_perm.append(data_list[0][indx][6])

        # Print sublists in row/column format.
        header = "Equivalent Permeability (microdarcys)"
        write_list(cast, header, list_perm)

        header = "Average Threshold Pressures (Pa)"
        write_list(cast, header, list_threshold)

    if data_list[1]:
        print(PARA + "Grid Parameters for User Fractures:", file=cast)

        # Create sub-lists.
        list_perm.clear()
        list_threshold.clear()
        for indx in range(0, total):
            list_threshold.append(data_list[1][indx][4])
            list_perm.append(data_list[1][indx][6])

        # Print sublists in row/column format.
        header = "Equivalent Permeability (microdarcys)"
        write_list(cast, header, list_perm)

        header = "Average Threshold Pressures (Pa)"
        write_list(cast, header, list_threshold)

    if data_list[2]:
        print(PARA + "Grid Parameters for Matrix:", file=cast)

        # Create sub-lists.
        list_perm.clear()
        list_threshold.clear()
        for indx in range(0, total):
            list_perm.append(data_list[2][indx][0])
            list_threshold.append(data_list[2][indx][1])

        # Print sublists in row/column format.
        header = "Equivalent Permeability (microdarcys)"
        write_list(cast, header, list_perm)

        header = "Average Threshold Pressures (Pa)"
        write_list(cast, header, list_threshold)

    # return None


def write_data_results(frac_controls, data_list, sim_numbr):
    """Write results to file.

    Parameters
    ----------
    frac_controls = (dict) dictionary of fracture parameters
    data_list = (list) stochastic parameters from each block for each case
    sim_numbr = (int) current simulation number

    Returns
    -------
    None
    """
    # Construct unique file name for saving data.
    output_file = scfg.FRAC_NAME + str(sim_numbr + 1) + scfg.EXTENSION_TXT

    # Construct full path to results directory and summary file.
    _, destination = sfile.get_path_name(scfg.OUTPUT_DIRECTORY, output_file,
                                         scfg.EXTENSION_TXT)

    # Write information to file on seal run.
    try:
        with open(destination, 'w', encoding="utf8") as cast:

            # Print current time stamp for printout.
            print("\n  Fracture Summary", file=cast)
            now = datetime.now()
            clump = IDNT + "Simulation Run Date: %a %d %b %Y - Time: %H:%M"
            stamp = now.strftime(clump)
            print(stamp, file=cast)

            # Print general information on run; add 1 to sim.
            print(IDNT + 'Analysis Title = '
                  f'{frac_controls.get("title"):<}', file=cast)
            print(IDNT + 'Simulation Number = '
                  f'{(sim_numbr + 1):4d}', file=cast)

            # Print emphasis for output.
            print("\n" + LINER_2, file=cast)
            print("  Results of Analysis  ---------------------", file=cast)
            print(LINER_2, file=cast)

            # Debug: Print fracture data for options selected, if desired.
            if ECHO:
                write_overall_results(cast, frac_controls, data_list)

            # Print results for each category.
            write_sections(cast, frac_controls, data_list)

            # Print current execution time at end using simple HH:MM:SS format.
            write_time(cast, frac_controls['elapsed'])

    except OSError as err:
        sfile.io_snag(err, scfg.OUTPUT_DIRECTORY, output_file)

    # return None


#
# -----------------------------------------------------------------------------
# - End of module
