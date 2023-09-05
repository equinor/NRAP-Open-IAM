#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions to model seal permeability aspects and history.

Author: Ernest N. Lindner
Date: 08/24/2022

Module Name
    seal_perm

Contents (16)
    evaluate_carbonate(react, carbonate)
    evaluate_clay(react, clay_type, clay_content)
    model_z_analysis(final_effect, delay_term, history)
    mineral_factor(reactivity, velocity, clay_content, clay_type, carbonate)
    compute_effective_saturation(co2_saturation, brine_residual, co2_residual)
    compute_capillary_pressure(zeta, beta, normal_saturation,
                               bubbling_pressure)
    compute_brine_pressure(base_co2_pressure, capillary_pressure)
    brine_relative_permeability(rel_model, effective_saturation, seal_controls)
    co2_relative_permeability(rel_model, effective_saturation, seal_controls)
    compute_current_permeability(rel_factor, total_permeability)
    ---
    define_permeability(grid, param_bounds)
    obtain_permeability(mid, scale, seal_controls, grid, param_bounds, rng)
    evaluate_areal_heter(grid_work, seal_controls)
    correlate_entry_pressure(seal_controls, grid)
    soluble_co2(seal_controls, brine_flow)
    print_tabular_data(title, data_list, columns, filename)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import random                       # For variability in values
import math                         # For Tanh + other functions

import seal_config as scfg          # IO directory and file names
import seal_file as sfile           # Input/output related operations
import seal_units as sun            # For unit conversion
import frac_random as fran          # Random value generators

# Seal model parameters
AXIS_SHIFT = 2.5                    # <S> - shift on x-axis for time model
HIGH_CLAY_SWELL = 1.0               # Model #2 - high swell degree
MID_CLAY_SWELL = 0.3                # Model #2 - med. swell degree
LOW_CLAY_SWELL = 0.0                # Model #2 - low. swell degree
V_CRITICAL = 0.5                    # Transition velocity (m/sec)
CARB_MINIMUM = 30.0                 # Minimum carb. content for reaction
CLAY_MINIMUM = 30.0                 # Minimum clay content for swell
CARB_CRITICAL = CARB_MINIMUM        # Critical carb. content for cementation
MAX_POSITIVE_CHANGE = 1.5           # Maximum + change due to carbonate
MAX_NEGATIVE_CHANGE = 0.9           # Maximum - change due to carbonate
MAX_NEGATIVE_SWELL = 0.99           # Maximum - change due to clay

# Permeability and heterogeneity controls
PERCENT = 8.0                       # Max. % cells for areal heterogeneity
MIN_SATURATION = 0.001              # Numerical minimum limit
MAX_SATURATION = 0.999              # Numerical maximum limit
V_SMALL_VARIANCE = 1.0e-05          # Small variability of log normal

# Other constants
BOZO = False                        # Control - print permeability ea loop
OUTPUT_COLS = 10                    # Columns for output in debugging
CONFIG = 11.3                       # Output format for debugging


def evaluate_carbonate(react, carbonate):
    """Compute carbonate factor for model #2.

    Parameters
    ----------
    react = (float) reactivity parameter (Range: 0 to 10)
    carbonate = (float) carbonate content as % (Range: 0.0 to 100.0)

    Returns
    -------
    change = (float) change factor on permeability for carbonate case
    """
    # Compute carbonate factor.
    change = (react / 10.0) * (carbonate / 100.0)

    return change


def evaluate_clay(react, clay_type, clay_content):
    """Compute clay factor for model #2.

    Parameters
    ----------
    react = (float) reactivity parameter (0 to 10)
    clay_type = (str) type of clay; must be:
        smectite = high
        illite = medium
        chlorite = low
    clay_content = (float) clay content as % (Range: 0.0 to 100.0)

    Returns
    -------
    change = (float) change factor on permeability for clay case
    """
    # Change based on clay type.
    if clay_type == "smectite":
        clay_factor = HIGH_CLAY_SWELL
    elif clay_type == "illite":
        clay_factor = MID_CLAY_SWELL
    else:
        clay_factor = LOW_CLAY_SWELL  # default = chlorite

    # Compute permeability change.
    change = clay_factor * (react / 10.0) * (clay_content / 100.0)
    return change


def model_z_analysis(final_effect, delay_term, history):
    """Construct logic to use time to compute influence factor.

    Parameters
    ----------
    final effect = (float) final ratio for change in permeability
    time_term = (float) time to see effect in time-line
    history = (float) cumulative time of alteration (yrs)

    Returns
    -------
    factor = (float) Influence factor (Range: 0 to 1.0)
    """
    # Compute equation terms.
    beta = math.tanh(AXIS_SHIFT)                # curve shift
    tau = (1.0 - final_effect) / (1.0 + beta)   # overall factor

    # Compute terms based on total flow time in cell.
    a_term = math.tanh(delay_term * history - AXIS_SHIFT)
    t_term = tau * (a_term + beta)

    # Compute influence as function of 1.0.
    factor = (1.0 - t_term)
    return factor


def mineral_factor(velocity, reactivity, clay_content, clay_type,
                   carbonate):
    """Construct the variables  to use Model #2 for influence factor.

    Parameters
    ----------
    velocity = (float) CO2 flow velocity (m/s)
    reactivity = (float) reactivity factor (1 to 10)
    clay_content = (float) clay content (%) (Range: 0.0 to 100.0)
    clay_type = (str) clay type name; must be:
        smectite = high
        illite = medium
        chlorite = low
    carbonate = (float) carbonate content (%) (Range: 0.0 to 100.0)

    Returns
    -------
    model_effect = (float) influence factor due to mineralogy
    """
    # Check if sufficient clay or carbonate content for reaction.
    if (carbonate >= CARB_MINIMUM) or (clay_content >= CLAY_MINIMUM):

        # - Check carbonate - if total carbonate content > 30%,
        #   ((carbonate model will dominate even if high clay)).
        if carbonate >= CLAY_MINIMUM:
            # Compute extent of change based on carbonate.
            change = evaluate_carbonate(reactivity, carbonate)

            # Determine if erosion or accumulation in this time step,
            # as a function of velocity.
            if velocity >= V_CRITICAL:
                model_effect = MAX_POSITIVE_CHANGE * change
            else:
                model_effect = 1.0 - (change * MAX_NEGATIVE_CHANGE)
        else:
            # Compute extent change based on clay model.
            change = evaluate_clay(reactivity, clay_type, clay_content)

            # Determine closure in this time step from swell;
            # - velocity has no effect.
            model_effect = 1.0 - (change * MAX_NEGATIVE_SWELL)
    else:
        # No change in influence.
        model_effect = 1.0

    return model_effect


def compute_effective_saturation(co2_saturation, brine_residual, co2_residual):
    """Compute the effective saturation value for two-phase flow.

    Parameters
    ----------
    co2_saturation = (float) current CO2 saturation (decimal)
    brine_residual = (float) brine residual saturation (decimal)
    co2_residual = (float) CO2 residual saturation (decimal)

    Returns
    -------
    effective_saturation => (float) normalized wet (brine) saturation
    """
    # Compute current brine saturation.
    wet_saturation = 1.0 - co2_saturation

    # Compute effective saturation within residual limits.
    if wet_saturation <= brine_residual:
        effective_saturation = 0.0
    elif wet_saturation > 1.0 - co2_residual:
        effective_saturation = 1.0
    else:
        effective_saturation = (wet_saturation - brine_residual) \
                                / (1.0 - brine_residual - co2_residual)
    return effective_saturation


def compute_brine_pressure(base_co2_pressure, capillary_pressure):
    """Compute the base brine pressures on seal horizon.

    Parameters
    ----------
    base_co2_pressure = (float) CO2 pressure at base from reservoir data
    capillary_pressure = (float) capillary pressure

    Returns
    -------
    base_pressure = (float) brine pressure at base of seal (MPa)
    """
    # Compute brine pressure from CO2 pressure.
    base_pressure = base_co2_pressure - capillary_pressure

    return base_pressure


def brine_relative_permeability(rel_model, effective_saturation,
                                seal_controls):
    """Compute the wetting relative permeabilities.

    Parameters
    ----------
    rel_model = (str) relative permeability model; either "LET" or "BC"
    effective_saturation = (float) effective saturation
    seal_controls = (dict) seal parameters

    Returns
    -------
    wet_perm = (float) wetting relative permeability (brine)

    Notes
    -----
    1. See Equation B-8 in manual for Brooks-Corey model.
    2. See Equation B-16 in manual for LET model.
    """
    # Define relative permeability based on model type.

    if 'LET' in rel_model:
        # Define wetting relative permeability using LET model.
        term_1 = pow(effective_saturation, seal_controls['l_wetting'])
        non_sat = 1.0 - effective_saturation
        term_2 = (seal_controls['e_wetting'] *
                  pow(non_sat, seal_controls['t_wetting']))
        wet_perm = term_1 / (term_1 + term_2)
    else:
        # Define nonwetting permeability using B_C, with zeta = lambda.
        term_1 = (2.0 + 3.0 * seal_controls['zeta']) / seal_controls['zeta']
        wet_perm = pow(effective_saturation, term_1)

    return wet_perm


def co2_relative_permeability(rel_model, effective_saturation, seal_controls):
    """Compute the nonwetting relative permeability.

    Parameters
    ----------
    rel_model = (str) relative permeability model; either "LET" or "BC"
    effective_saturation = (float) effective wetting saturation
    seal_controls = (dict) seal parameters

    Returns
    -------
    nonwet_perm = (float) nonwetting relative permeability (co2)

    Notes
    -----
    1. See Equation B-9 in manual for Brooks-Corey model.
    2. See Equation B-17 in manual for LET model.
    """
    # Define nonwetting saturation and ratio.
    non_saturation = 1.0 - effective_saturation
    p_ratio = seal_controls['perm_ratio']

    # Define relative permeability based on model type.

    if 'LET' in rel_model:
        # Define nonwetting permeability using LET model.
        term_1 = pow(non_saturation, seal_controls['l_nonwet'])
        term_2 = (seal_controls['e_nonwet'] *
                  pow(effective_saturation, seal_controls['t_nonwet']))
        fcap = term_1 / (term_1 + term_2)
        nonwet_perm = p_ratio * fcap
    else:
        # Define nonwetting permeability using B_C, with zeta = lambda.
        term_1 = pow(non_saturation, 2.0)
        term_2 = (2.0 + seal_controls['zeta']) / seal_controls['zeta']
        term_3 = pow(effective_saturation, term_2)
        nonwet_perm = p_ratio * term_1 * (1.0 - term_3)

    return nonwet_perm


def compute_capillary_pressure(rel_model, normal_saturation, seal_controls,
                               entry_pressure):
    """Compute the capillary pressure based on effective saturation.

    Parameters
    ----------
    rel_model = (str) relative permeability model = BC/LET
    normal_saturation = (float) normalized saturation between residual
           values = effective wetting saturation
    seal_controls = (dict) seal parameters
    entry_pressure = (float) limit pressure (Pa) for capillary pressure

    Returns
    -------
    capillary_pressure = (float) capillary pressure (Pa)

    Notes
    -----
    1. See Equation B-15a/b in manual for Brooks-Corey model.
    2. See Equations B-18 and B-19 in manual for LET model.
    """
    # Define nonwetting saturation.
    non_saturation = 1.0 - normal_saturation

    # Define capillary pressure based on model type.
    if 'LET' in rel_model:
        # Compute based on LET capillary terms.
        term_1 = pow(non_saturation, seal_controls['l_capillary'])
        term_2 = (seal_controls['e_capillary'] *
                  pow(normal_saturation, seal_controls['t_capillary']))
        fcap = term_1 / (term_1 + term_2)
        change = seal_controls['max_capillary'] - entry_pressure
        capillary_pressure = fcap * change + entry_pressure
    else:
        # Compute based on modified Brooks-Corey model.
        exponent = (1.0 / seal_controls['zeta'])
        if normal_saturation >= MAX_SATURATION:
            divisor = 1.0
        elif normal_saturation <= MIN_SATURATION:
            divisor = pow(MIN_SATURATION, exponent)
        else:
            divisor = pow(normal_saturation, exponent)

        capillary_pressure = (entry_pressure / divisor)

    return capillary_pressure  # in Pascals


def compute_current_permeability(rel_factor, affect, total_permeability):
    """Compute the current effective permeability of cell in m^2.

    Parameters
    ----------
    rel_factor = (float) relative permeability of cell (decimal)
    affect = (float) influence factor on permeability
    total permeability = (float) effective permeability of cell (mD)

    Returns
    -------
    effective_permeability = (float) permeability at base of seal (m2)

    Notes
    -----
    Input is originally in microdarcy but is converted to m2!
    """
    # Compute brine pressure from CO2 pressure; express in m^2.
    effective_permeability = (rel_factor * affect * total_permeability)
    effective_permeability *= sun.microd_to_metersq()

    return effective_permeability


def define_permeability(grid, param_bounds):
    """Define permeability values from file.

    Parameters
    ----------
    grid = (list) (class) a collection of cells
    param_bounds = (dict) input bounds

    Returns
    -------
    Data from file - (array) 1D NumPy data
    """
    # Get data array from file.
    data_array = sfile.acquire_data_array(scfg.INPUT_DIRECTORY, scfg.PERM_NAME,
                                          len(grid))

    # Check data to be within bounds.
    val_limits = [param_bounds['perm_min'][0],
                  param_bounds['perm_max'][1]]
    stage = sfile.check_file_floats(data_array, "Permeability input value",
                                    val_limits)
    if stage < 0:
        sfile.opx_problem("File Value(s) Outside Defined Bounds " +
                          "Caused Program to Exit!")

    return data_array


def obtain_permeability(seal_controls, grid, param_bounds, rng):
    """Define vertical permeability as variable or from file.

    Parameters
    ----------
    seal_controls = (dict) seal parameters
    grid = (list) (class) a collection of cells
    param_bounds = (dict) input bounds
    rng = (generator) random number generator

    Returns
    -------
    grid = grid = (list) (class) a collection of cells

    Notes
    -----
    1. Input values are in microdarcys.
    2. mid = location of log-distribution.
    3. scale = scale of log-distribution.
    """
    # Get data from file, if desired.
    if seal_controls['perm_input_approach']:
        # Read from file, but only once.
        if seal_controls['read_perm']:
            data_array = define_permeability(grid, param_bounds)
            for numbr, cell in enumerate(grid):
                cell.permeability = data_array[numbr]
            seal_controls['read_perm'] = False    # No additional reads!
    else:
        # Otherwise, set random values for use in loop.
        mid = seal_controls['perm_location']
        scale = seal_controls['perm_scale']

        # Define permeability value for every cell in grid.
        for cell in grid:
            if seal_controls['vary_perm_choice']:
                # For variable permeability, evaluate random value.
                perm_value = \
                    fran.evaluate_lognorm(mid, scale,
                                          seal_controls['perm_min'],
                                          seal_controls['perm_max'],
                                          rng)
            else:
                # For Uniform permeability, set to mean.
                perm_value = seal_controls['perm_mean']

            # Then set cell value.
            cell.permeability = perm_value

    # Debugging check - printout permeability for each simulation!
    if BOZO:
        # Create list of permeability data.
        print_list = []
        for cell in grid:
            print_list.append(cell.permeability)

        # Print data in tabular format.
        print_tabular_data("PERMEABILITY", print_list, OUTPUT_COLS,
                           "permeability_data")

    return grid


def evaluate_areal_heter(grid_work, seal_controls):
    """Generate random zones across seal grid of higher permeability.

    Parameters
    ----------
    grid_work = (list of class) a collection of cells
    seal_controls = (dict) control parameters

    Returns
    -------
    grid = (list) (class) a collection of cells
    """
    # Only implement if desired and if variability is used.
    if (seal_controls['heterogeneity_approach'] and
            seal_controls['vary_perm_choice']):

        # Define number of cells for inclusion.
        total_cells = seal_controls['num_cells']
        heter_limit = int(PERCENT * total_cells / 100.0)

        # Create list of cell index order and randomly shuffle.
        extent = len(grid_work)
        index_set = list(range(0, extent))
        random.shuffle(index_set)       # Eliminates duplicates!

        # Loop over grid until maximum number of random cells is adjusted.
        heter_count = 0
        for element in range(extent):
            heter_index = index_set[element]
            if grid_work[heter_index].status > 0:
                # Break loop if at limit.
                heter_count += 1
                if heter_count > heter_limit:
                    break

                # increase permeability of cell if still in loop.
                cell_perm = (grid_work[heter_index].permeability
                             * seal_controls['perm_heter_factor'])
                grid_work[heter_index].permeability = cell_perm

    return grid_work


def correlate_entry_pressure(seal_controls, grid):
    """Correlate the threshold value of cell with permeability.

    Parameters
    ----------
    seal_controls = (dict) control parameters
    grid = (list) (class) a collection of cells

    Returns
    -------
    grid = (list) (class) a collection of cells
    """
    # Change entry pressure only if desired and permeability is
    #    <not> from file.
    if seal_controls['correlate_entry_approach'] and \
            not seal_controls['entry_approach']:

        # Define parameters for clarity.
        entry_pressure = seal_controls['entry_pressure']
        reference_perm = seal_controls['perm_mean']

        # Define new threshold pressure for each cell.
        for cell in grid:
            if cell.permeability > 1.0e-20:
                ratio = math.sqrt(reference_perm / cell.permeability)
                cell.entry = entry_pressure * ratio

    return grid


def soluble_co2(seal_controls, brine_flow):
    """Compute amount of CO2 dissolved in brine (in m3/s).

    Parameters
    ----------
    seal_controls = (dict) dictionary of parameters
    brine_flow = (float) amount of brine (m^3/s)

    Returns
    -------
    result = soluble CO2 (m^3/s)

    Note
    ----
    1. Assumes brine is 100% saturated with CO2.
    2. Units defined in lookup table!
    """
    # Define variables for clarity.
    soluble = seal_controls['co2_solubility']
    brine_dense = seal_controls['brine_density']
    co2_dense = seal_controls['co2_density']

    # Define brine amount per time.
    # --> Define weight (kg) of brine per sec.
    brine_kg = brine_dense * brine_flow   # amount of brine flow - kg per sec.

    # Define solubility = find g-CO2 per g-brine,
    # soluble units: mol/kg
    # ---> soluble * sun.co2_molar_mass() / sun.kilogram_to_gram()
    co2_dissolved = soluble * sun.co2_molar_mass() / sun.kilogram_to_gram()

    # --> Compute volume of CO2 - in m^3/s. kg-CO2/kg-brine.mol/kg-brine.
    result = co2_dissolved * brine_kg / co2_dense

    return result


def print_tabular_data(title, data_list, columns, filename):
    """Print data values in arbitrary column/row to output directory.

    Parameters
    ----------
    title = (str) description of data
    data_list = (list) data values
    columns = (int) arbitrary number of columns for printing
    filename = output file path

    Returns
    -------
    None

    Notes
    -----
    1. Pattern => CONFIG = 11.3 (see constants).
    2. Columns is a printer-based choice and not the same as array dimension.
    5. See response at:
    https://stackoverflow.com/questions/9535954/printing-lists-as-tabular-data
    """
    # Get output file path.
    sub_path, destination = sfile.get_path_name(scfg.OUTPUT_DIRECTORY,
                                                filename,
                                                scfg.EXTENSION_TXT)

    # Compute the number of lines and the "remainder" columns on last line.
    n_lines, remainder = divmod(len(data_list), columns)

    # Setup format and controls.
    pattern = f'{{:{CONFIG}e}}'
    typ_line = '\n'.join(pattern * columns for _ in range(n_lines))
    last_line = pattern * remainder

    # Write information to output.
    try:
        with open(destination, 'w', encoding="utf8") as zum:
            # Print data with formatting.
            print("  " + title, file=zum)
            print(typ_line.format(*data_list), file=zum)
            print(last_line.format(*data_list[n_lines*columns:]), file=zum)
    except OSError as err:
        sfile.io_snag(err, sub_path, scfg.DEBUG_FILE)

    # return None


#
# -----------------------------------------------------------------------------
# - End of module
