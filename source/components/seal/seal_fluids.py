#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions for interpolation of density and viscosity.

Author: Ernest N. Lindner
Date: 08/18/2022

Module Name
    seal_fluids

Contents (13)
    conduct_2d_interp(vect_x, vect_y, interp_table, pressure, temperature)
    debug_shape(method, dotted, temperature, pressure, salinity)
    linear_interpolate(target, x1, x2, y1, y2)
    salinity_locate(target,vector)
    convert_salinity(salinity)
    define_lookup_array(option)
    co2_properties(pressure, temperature)
    brine_density_property(pressure, temperature, salinity=0.0)
    brine_viscosity_property(pressure, temperature, salinity=0.0)
    co2_solubility(pressure, temperature, salinity=0.0)
    ----
    reset_fluid_properties(seal_controls, property_list)
    interpolate_fluid_properties(seal_controls, pressure)
    manage_interpolation(alive, seal_controls, press_top, param_bounds)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import numpy as np
from scipy import interpolate

import data_arrays as sdata         # Data arrays indices
import seal_config as scfg           # IO directory and file names
import seal_file as sfile           # Error in file operations.
import seal_model as mod            # Functions for Class definitions
import seal_units as sun            # For unit conversion
import seal_upload as sup           # Input file data

# Other constants
ECHO = False                        # DEBUG Print variables


def conduct_2d_interp(vect_x, vect_y, interp_table, pressure, temperature):
    """Perform 2D interpolation on a density/viscosity table.

    Parameters
    ----------
    vectx = (list) vector of x values for table
    vecty = (list) vector of y values for table
    interp_table2d = (list) 2D table of values
    pressure = (float) pressure (MPa) at point
    temperature = (float) temperature (oC) at point

    Returns
    -------
    valu = (float) density or viscosity at given temperature and pressure
    """
    # Define vectors and array for interpolation.
    x_values = vect_x
    y_values = vect_y
    z_values = interp_table

    # Use SciPy function interp2d to obtain density value.
    condition = interpolate.interp2d(x_values, y_values, z_values,
                                     kind='cubic')
    valu = float(condition(pressure, temperature))

    # Function completed.
    return valu


def debug_shape(bmethod, dotted, temperature, pressure, salinity):
    """Write shape features for debugging.

    Parameters
    ----------
    bmethod = (int) interpolation method dimension
    dotted = (int) array position
    temperature = (array) temperature (oC)
    pressure = (array) reference pressure (Pa)
    salinity = (float) salinity (ppm)

    Returns
    -------
    None
    """
    # Print shapes without comma.
    print(f'   ** Interpolation Type: {bmethod}')
    print(f'   ** Salinity Index: {dotted}')

    temp_val = str(temperature.shape).replace(',', '')
    print(f'   ** Shape of Temperature Vector: {temp_val}')

    temp_val = str(pressure.shape).replace(',', '')
    print(f'   ** Shape of Pressure Vector {temp_val}')

    temp_val = str(salinity.shape).replace(',', '')
    print(f'   ** Shape of Salinity Array {temp_val}')

    # return None


def linear_interpolate(target, x_pt1, x_pt2, y_pt1, y_pt2):
    """Perform linear 1D interpolation of value from series.

    Parameters
    ----------
    target = (float) current salinity (x-value)
    x_pt1 = (float) x-value (salinity) at point 1
    x_pt2 = (float) x-value (salinity) at point 2
    y_pt1 = (float) function(density or viscosity) at point 1
    y_pt2 = (float) function(density or viscosity) at point 2

    Returns
    -------
    value = (float) interpolated value

    Notes
    -----
    See Wikipedia for equation.
    """
    # Define values for computation.
    del_x = x_pt2 - x_pt1
    del_y = y_pt2 - y_pt1

    # Return if denominator is OK.
    if del_x != 0.0:
        value = y_pt1 + ((del_y/del_x) * (target - x_pt1))
    else:
        value = 0.0
        msg1 = "Attempt to Interpolate Failed "
        msg1 += "during Input of Linear Interpolate."
        sfile.opx_problem(msg1)

    # return value.
    return float(value)


def salinity_locate(target, vector):
    """Search array to find appropriate location in array.

    Parameters
    ----------
    target = (float) specified salinity
    vector = (array) vector of salinity values controlling interpolation

    Returns
    -------
    answer = (int) index in vector where value is == or <= than target
    analysis = (int) type of analysis to be conducted (=2,3):
        2 = 2D cubic interpolation at constant salinity
        3 = 3D cubic/linear interpolations at different salinities
    """
    # Set search values.
    min_index = 0
    max_index = vector.size - 1
    min_of_array = vector[min_index]
    max_of_array = vector[max_index]
    analysis = 0
    answer = 0

    # Check if target is equal to any control points; use 2D analysis.
    if target in vector:
        answer = np.where(vector > target)[0][0] - 1  # zeros>tuple of indices
        analysis = 2
    # Check if target is less than interpolation range.
    elif target <= min_of_array:
        answer = min_index   # -> use lowest salinity array.
        analysis = 2
    # Check if target is greater than interpolation range.
    elif target >= max_of_array:
        answer = max_index   # -> use highest salinity array.
        analysis = 2
    # Otherwise, conduct basic search for values in-between controls points.
    else:
        for trak in range(0, max_index-1):
            if vector[trak] < target <= vector[trak+1]:
                answer = trak
                analysis = 3

    return answer, analysis


def convert_salinity(salinity):
    """Convert salinity in ppm to Molal for brine viscosity.

    Parameters
    ----------
    salinity = (float) NaCl in ppm

    Returns
    -------
    molal = (float) NaCl in Molal (mol/kg)

    Notes
    -----
    1. molar mass is in [g/kg]
    2. Check: 100000 ppm =. 0.1 weight% converts to 1.9 m in NACl
    """
    concentration = salinity * sun.ppm_convert()       # weight ratio
    solvent = (1 - concentration)                      # solvent alone
    molal = salinity / (sun.nacl_molar_mass() * sun.kilogram_to_gram())
    molal /= solvent

    return molal


def define_lookup_array(option):
    """Get full path to the lookup file and load data.

    Parameters
    ----------
    option = (str) lookup file name

    Returns
    -------
    NumPy array = (array) density or viscosity array
    """
    # Get file path.
    file_path = sfile.find_local_path(option)

    # Load NumPy file and return.
    lookup_array = np.load(str(file_path))

    return lookup_array


def co2_properties(pressure, temperature):
    """Evaluate density and viscosity parameters for CO2.

    Parameters
    ----------
    pressure = (float) pressure of CO2 (MPa)
    temperature = (float) temperature of CO2 (oC)

    Returns
    -------
    z_density = (float) CO2 density (kg/m3)
    viscosity = (float) CO2 viscosity (Pa-s)

    Notes
    -----
    1. CO2 Density from NIST REFPROP software, based on
        Span & Wagner (1996). Range 0.1 to 60 MPa and 1 to 180 oC.
    2. CO2 viscosity from NIST REFPROP software, based on
        Fenghour et al. (1998); Range 0.1 to 60 MPa and 1 to 180 oC.
    """
    # Load density and viscosity LUTs.
    co2_density_lut = define_lookup_array(scfg.CO2_DENSITY_FILE)
    co2_viscosity_lut = define_lookup_array(scfg.CO2_VISCOSITY_FILE)

    # Use 2D cubic interpolation to estimate values from table.
    z_density = conduct_2d_interp(sdata.CO2_DENSITY_PRESSURE,
                                  sdata.CO2_DENSITY_TEMP,
                                  co2_density_lut,
                                  pressure, temperature)

    z_viscosity = conduct_2d_interp(sdata.CO2_VISCOSITY_PRESSURE,
                                    sdata.CO2_VISCOSITY_TEMP,
                                    co2_viscosity_lut,
                                    pressure, temperature)
    return z_density, z_viscosity


def brine_density_property(pressure, temperature, salinity=0.0):
    """Estimate brine density at specific temperature, pressure, salinity.

    Parameters
    ----------
    pressure = (float) pressure of brine (MPa)
    temperature = (float) temperature of brine (oC)
    salinity = (float) NaCl in solution (ppm)

    Returns
    -------
    z_density = (float) brine density (kg/m3)

    Notes
    -----
    1. Temperature/pressure/salinity-dependent values based on
        Sun et al. (2008) + density of pure H2O provided by
        NIST (2010) (Tables from Sun; compiled/provided by J Morrell
    2. NIST data based on Wagner & Pruss [1995]).
    3. EOS valid to 100 MPa, 374 oC and salinity of
        ~80,000 mg/kg (ratio of mass of solute / mass of solution, or ppm).
    """
    # Identify interpolation method & index in 3D array:
    #   -> It is a function of (salinity, pressure, temperature)
    #   -> method = 2 for 2D interpolation - cubic
    #   -> method = 3 for 3D interpolation - cubic + linear (salinity)
    #   -> dot = position in array
    dot, method = salinity_locate(salinity, sdata.BRINE_DENSITY_SALINITY)

    # Catch error, if it occurs.
    if method not in (2, 3):
        sfile.opx_problem("Failure in Interpolation during Brine Density "
                          + "Calculation.")

    # Debug print.
    if ECHO:
        print("\n  For Brine Density")
        debug_shape(method, dot,
                    sdata.BRINE_DENSITY_TEMP,
                    sdata.BRINE_DENSITY_PRESSURE,
                    sdata.BRINE_DENSITY_SALINITY)
        print(f'   ** Salinity in ppm: = {salinity}')

    # STEP #1: 2D Interpolation  ----------------------------------------------

    # Establish variables and arrays for operation.
    # indx_1 = len(sdata.BRINE_DENSITY_TEMP)
    # indx_2 = len(sdata.BRINE_DENSITY_PRESSURE)

    # Load density LUT.
    brine_density_lut = define_lookup_array(scfg.BRINE_DENSITY_FILE)

    # Establish 2D subarray for density at low index salinity.
    # density_low = np.zeros((indx_1, indx_2))
    density_low = brine_density_lut[dot, :, :]

    # Conduct 2D interpolation with low array.
    z_density = conduct_2d_interp(sdata.BRINE_DENSITY_PRESSURE,
                                  sdata.BRINE_DENSITY_TEMP,
                                  density_low, pressure, temperature)

    # STEP #2: 3D Interpolation  ----------------------------------------------

    # Perform 3d interpolation across salinities, if required (method = 3).
    if method == 3:

        # Establish subarray for density at high index salinity.
        # density_high = np.zeros((indx_1, indx_2))
        density_high = brine_density_lut[(dot + 1), :, :]

        # Conduct 2D interpolation with low array.
        z_density_up = conduct_2d_interp(sdata.BRINE_DENSITY_PRESSURE,
                                         sdata.BRINE_DENSITY_TEMP,
                                         density_high, pressure, temperature)

        # Conduct linear interpolations between values.
        low_salinity = sdata.BRINE_DENSITY_SALINITY[dot]
        high_salinity = sdata.BRINE_DENSITY_SALINITY[dot + 1]
        z_density = linear_interpolate(salinity,
                                       low_salinity, high_salinity,
                                       z_density, z_density_up)

    return z_density


def brine_viscosity_property(pressure, temperature, salinity=0.0):
    """Establish viscosity at specific temperature, pressure, and salinity.

    Parameters
    ----------
    pressure = (float) pressure of brine (MPa)
    temperature = (float) temperature of brine (oC)
    salinity = (float) NaCl in ppm (to be converted to molal)

    Returns
    -------
    z_viscosity = (float) brine viscosity (Pa-s)

    Notes
    -----
    1. Temperature/pressure/salinity-dependent viscosity of a
        saline solution based on Mao and Duan (2009) in terms
        of molality + viscosity of pure H2O provided by
        NIST (2010).
    2. NIST data based on Huber et al. [2009]
        (Tables from Mao & Duan data provided by J Morrell).
    3. EOS valid to 100 MPa, 350 oC and salinity to 6 M.
    """
    # Identify interpolation method & index in 3D array:
    #   -> It is a function of (salinity, pressure, temperature)
    #   -> method = 2 for 2D interpolation - cubic
    #   -> method = 3 for 3D interpolation - cubic + linear (on saline)
    #   -> dot = position in array.

    # Note: lookup table is in Molal, not ppm!
    sal_molal = convert_salinity(salinity)
    dot, method = salinity_locate(sal_molal, sdata.BRINE_VISCOSITY_SALINITY)

    # Debug.
    if ECHO:
        print("\n  For Brine Viscosity")
        debug_shape(method, dot,
                    sdata.BRINE_VISCOSITY_TEMP,
                    sdata.BRINE_VISCOSITY_PRESSURE,
                    sdata.BRINE_VISCOSITY_SALINITY)
        print(f'   ** Salinity in Molal: = {sal_molal}')

    # Catch error, if it occurs.
    if method not in (2, 3):
        sfile.opx_problem("Failure in Interpolation during Brine Viscosity "
                          + "Calculation.")

    # STEP #1: 2D Interpolation  ----------------------------------------------

    # Establish variables for operation.
    # indx_1 = len(sdata.BRINE_VISCOSITY_TEMP)
    # indx_2 = len(sdata.BRINE_VISCOSITY_PRESSURE)

    # Load viscosity LUT.
    brine_viscosity_lut = define_lookup_array(scfg.BRINE_VISCOSITY_FILE)

    # Establish subarray for viscosity at low index salinity.
    # viscosity_low = np.zeros((indx_1, indx_2))
    viscosity_low = brine_viscosity_lut[dot, :, :]

    # Conduct 2D interpolation with low array.
    z_viscosity = conduct_2d_interp(sdata.BRINE_VISCOSITY_PRESSURE,
                                    sdata.BRINE_VISCOSITY_TEMP,
                                    viscosity_low, pressure, temperature)

    # STEP #2: 3D Interpolation  ----------------------------------------------

    # Perform 3d interpolation across salinities, if required (method = 3).
    if method == 3:

        # Establish subarray for viscosity at high index salinity.
        # viscosity_high = np.zeros((indx_1, indx_2))
        viscosity_high = brine_viscosity_lut[(dot + 1), :, :]

        # Conduct 2D interpolation with low array.
        z_viscosity_up = conduct_2d_interp(sdata.BRINE_VISCOSITY_PRESSURE,
                                           sdata.BRINE_VISCOSITY_TEMP,
                                           viscosity_high, pressure,
                                           temperature)

        # Conduct linear interpolations between values.
        low_salinity = sdata.BRINE_VISCOSITY_SALINITY[dot]
        high_salinity = sdata.BRINE_VISCOSITY_SALINITY[dot + 1]
        z_viscosity = linear_interpolate(sal_molal,
                                         low_salinity, high_salinity,
                                         z_viscosity, z_viscosity_up)

    return z_viscosity


def co2_solubility_property(pressure, temperature, salinity=0.0):
    """Estimate CO2 solubility in brine at temp., pressure and salinity.

    Parameters
    ----------
    pressure = (float) pressure of brine (MPa)
    temperature = (float) temperature of brine (oC)
    salinity = (float) brine (NaCl) in mol

    Returns
    -------
    z_solubility = (float) CO2 solubility (mol/kg)

    Notes
    -----
    Temperature/pressure/salinity-dependent values based on
    programming the equation of state correlations as reported
    by Duan et al. (2006) and computing required values by
    Jason Monnell, Research Assistant Professor, University of Pittsburgh.
    """
    # Convert ppm to mol for analysis.
    sal_molal = convert_salinity(salinity)

    # Identify interpolation method & index in 3D array:
    #   -> It is a function of (salinity, pressure, temperature)
    #   -> method = 2 for 2D interpolation - cubic
    #   -> method = 3 for 3D interpolation - cubic + linear (salinity)
    #   -> dot = position in array
    dot, method = salinity_locate(sal_molal, sdata.CO2_SOLUBILITY_SALINITY)

    # Catch error, if it occurs.
    if method not in (2, 3):
        sfile.opx_problem("Failure in Interpolation during CO2 Solubility "
                          + "Calculation.")

    # Debug print.
    if ECHO:
        print("\n  For CO2 solubility")
        debug_shape(method, dot,
                    sdata.CO2_SOLUBILITY_TEMP,
                    sdata.CO2_SOLUBILITY_PRESSURE,
                    sdata.CO2_SOLUBILITY_SALINITY)
        print(f'   ** Salinity in mol: = {sal_molal}')

    # STEP #1: 2D Interpolation  ----------------------------------------------

    # Establish variables and arrays for operation.
    # indx_1 = len(sdata.CO2_SOLUBILITY_TEMP)
    # indx_2 = len(sdata.CO2_SOLUBILITY_PRESSURE)

    # Read LUT and create NPY file:
    #  co2_solubility_lut = sol_data.CO2_SOLUBILITY_TABLE
    #  target = "data_co2_solubility"
    #  file_path = sfile.find_local_path(target)
    #  np.save(file_path, co2_solubility_lut).

    #  Load solubility LUT.
    co2_solubility_lut = define_lookup_array(scfg.CO2_SOLUBILITY_FILE)

    # Establish 2D subarray for solubility at low index salinity.
    # solubility_low = np.zeros((indx_1, indx_2))
    solubility_low = co2_solubility_lut[dot, :, :]

    # Conduct 2D interpolation with low array.
    z_solubility = conduct_2d_interp(sdata.CO2_SOLUBILITY_TEMP,
                                     sdata.CO2_SOLUBILITY_PRESSURE,
                                     solubility_low, temperature, pressure)

    # STEP #2: 3D Interpolation  ----------------------------------------------

    # Perform 3d interpolation across salinities, if required (method = 3).
    if method == 3:

        # Establish subarray for solubility at high index salinity.
        # solubility_high = np.zeros((indx_1, indx_2))
        solubility_high = co2_solubility_lut[(dot + 1), :, :]

        # Conduct 2D interpolation with low array.
        z_solubility_up = conduct_2d_interp(sdata.CO2_SOLUBILITY_TEMP,
                                            sdata.CO2_SOLUBILITY_PRESSURE,
                                            solubility_high, temperature,
                                            pressure)

        # Conduct linear interpolations between values.
        low_salinity = sdata.CO2_SOLUBILITY_SALINITY[dot]
        high_salinity = sdata.CO2_SOLUBILITY_SALINITY[dot + 1]
        z_solubility = linear_interpolate(sal_molal,
                                          low_salinity, high_salinity,
                                          z_solubility, z_solubility_up)

    return z_solubility


def reset_fluid_properties(seal_controls, property_list):
    """Update fluid parameters in Class and in seal_controls.

    Parameters
    ----------
    seal_controls = (dict) seal control parameters
    property_list (list):
        [0] = co2_density = new CO2 density
        [1] = co2_viscosity = new CO2 viscosity
        [2] = brine_density = new brine density
        [3] = brine_viscosity = new brine viscosity
        [4] = co2solubility = new CO2 solubility

    Returns
    -------
    seal_controls = (dict) seal controls (updated)
    """
    # Translate list.
    co2_density = property_list[0]
    co2_viscosity = property_list[1]
    brine_density = property_list[2]
    brine_viscosity = property_list[3]
    co2_solubility = property_list[4]

    # Reset Cell Class variables.
    mod.Cell.co2Density = co2_density
    mod.Cell.co2Viscosity = co2_viscosity
    mod.Cell.brineDensity = brine_density
    mod.Cell.brineViscosity = brine_viscosity
    mod.Cell.co2Solubility = co2_solubility

    # Reset seal controls for output.
    seal_controls['co2_density'] = co2_density
    seal_controls['co2_viscosity'] = co2_viscosity
    seal_controls['brine_density'] = brine_density
    seal_controls['brine_viscosity'] = brine_viscosity
    seal_controls['co2_solubility'] = co2_solubility

    return seal_controls


def interpolate_fluid_properties(seal_controls, pressure):
    """Interpolates to get fluid properties for brine & CO2.

    Parameters
    ----------
    seal_controls = (dict) seal controls (updated)
    pressure = (float) reference pressure (Pa)

    Returns
    -------
    seal_controls = (dict) seal controls (updated)
    """
    # Convert pressure to MPa and set other computation parameters.
    pressure *= sun.pa_to_mpa()
    temperature = seal_controls['temperature']
    salinity = seal_controls['salinity']

    # Lookup CO2 density, viscosity from lookup tables.
    co2_density, co2_viscosity = co2_properties(pressure, temperature)
    brine_density = brine_density_property(pressure, temperature, salinity)
    brine_viscosity = brine_viscosity_property(pressure, temperature, salinity)
    co2_solubility = co2_solubility_property(pressure, temperature, salinity)

    # Create list for transfer for reset.
    property_list = [co2_density, co2_viscosity, brine_density,
                     brine_viscosity, co2_solubility]

    # Reset controls within model and controls.
    reset_fluid_properties(seal_controls, property_list)

    # Debug print for check.
    if ECHO:
        # Debug - print variables for CO2 and brine.
        echo_txt = ('\n   ** CO2 density and viscosity = '
                    f'{co2_density:.3f} kg/m3 '
                    f'and {co2_viscosity:.5e} Pa-s')
        print(echo_txt)

        echo_txt = ('   ** Brine density and viscosity = '
                    f'{brine_density:.3f} kg/m3 '
                    f'and {brine_viscosity:.5e} Pa-s')
        print(echo_txt)

        print(f'   ** CO2 solubility = {co2_solubility:.3f} mol/kg')

    return seal_controls


def manage_interpolation(alone, seal_controls, press_top, param_bounds):
    """Manage interpolation of fluid properties for brine & CO2.

    Parameters
    ----------
    alone = (bool) flag on stand-alone operations
    seal_controls = (dict) seal controls
    press_top = (float) pressure at top of seal (NumPy array)
    param_bounds = (dict) parameter bounds

    Returns
    -------
    seal_controls = (dict) seal controls (updated)
    """
    # Generate new values as desired by user.
    if seal_controls['interpolate_approach']:

        # Echo status to user.
        sfile.echo_status(alone, "INTERPOLATING FLUID PROPERTIES.")

        # Compute average pressure top and bottom of seal.
        base_press_array, dummy = \
            sup.input_reservoir_data(1, seal_controls, param_bounds)
        base_pressure_ave = np.average(base_press_array)
        top_pressure_ave = np.average(press_top)

        # Compute average pressure at center of seal for interpolation.
        ave_pressure = (base_pressure_ave + top_pressure_ave) / 2.0

        # Interpolate data @ average pressure of seal center & reset values.
        seal_controls = \
            interpolate_fluid_properties(seal_controls, ave_pressure)

    return seal_controls


#
# -----------------------------------------------------------------------------
# - End of module
