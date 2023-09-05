#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions for interpolation of density and viscosity.

Author: Ernest N. Lindner
Date: 08/23/2022

Module Name
    flt_fluids

Contents (14)
    conduct_2d_interp(vect_x, vect_y, interp_table, pressure, temperature)
    linear_interpolate(target, x1, x2, y1, y2)
    salinity_locate(target,vector)
    convert_salinity(salinity)
    define_lookup_array(option)
    co2_properties(pressure, temperature)
    brine_density_property(pressure, temperature, salinity=0.0)
    brine_viscosity_property(pressure, temperature, salinity=0.0)
    co2_solubility(pressure, temperature, salinity=0.0)
    ave_properties(sub_list, aquifer_list)
    ---
    interpolate_fluid_properties(temperature, pressure, salinity)
    reset_aqui_properties(fault_controls, aqui_list)
    reset_fluid_properties(fault_controls, new_list)
    debug_shape(method, dotted, temperature, pressure, salinity)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import numpy as np
from scipy import interpolate

import data_arrays as fdata     # Data arrays indices
import flt_units as funit       # For unit conversion
import flt_file as fileop       # For file operations
import flt_config as scfg       # IO directory and file names

# Other constants
ECHO = False                    # DEBUG Print variables


def conduct_2d_interp(vect_x, vect_y, interp_table, pressure, temperature):
    """Perform 2D interpolation on a density/viscosity table.

    Parameters
    ----------
    vectx = (array) vector of x values for table
    vecty = (array) vector of y values for table
    interp_table2d = (array) 2D table of values
    pressure = (float) pressure (MPa) at point
    temperature = (float) temperature(oC) at point

    Returns
    -------
    value = (float) density or viscosity at given temperature and pressure
    """
    # Define vectors and array for interpolation.
    x_values = vect_x
    y_values = vect_y
    z_values = interp_table

    # Use SciPy function interp2d to obtain density value.
    condition = interpolate.interp2d(x_values, y_values, z_values,
                                     kind='cubic')
    value = float(condition(pressure, temperature))

    return value


def linear_interpolate(target, x_pt1, x_pt2, y_pt1, y_pt2):
    """Perform linear 1D interpolation of value from series.

    Parameters
    ----------
    target = (float) current salinity (x-value)
    x_pt1  = (float) x-value (salinity) at point 1
    x_pt2  = (float) x-value (salinity) at point 2
    y_pt1  = (float) function(density or viscosity) at point 1
    y_pt2  = (float) function(density or viscosity) at point 2

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
    value = 0.0

    # Return if denominator is OK.
    if del_x != 0.0:
        value = y_pt1 + ((del_y/del_x) * (target - x_pt1))
    else:
        msg1 = "Attempt to Linear Interpolate Failed "
        msg1 += "During Input of Linear Interpolate."
        fileop.opx_problem(msg1)

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
    analysis = (int) type of analysis to be conducted
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
        answer = np.where(vector > target)[0][0] - 1
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
    concentration = salinity * funit.ppm_convert()    # weight ratio
    solvent = (1 - concentration)                     # solvent alone
    molal = salinity / (funit.nacl_molar_mass() * funit.kilogram_to_gram())
    molal /= solvent

    return molal


def define_lookup_array(option):
    """Get the full path to the lookup file and load data.

    Parameters
    ----------
    option = (str) lookup file name

    Returns
    -------
    lookup_array = (array) density or viscosity array
    """
    # Get file path.
    file_path = fileop.find_local_path(option)

    # Load NumPy file and return.
    lookup_array = np.load(file_path)

    return lookup_array


def co2_properties(pressure, temperature):
    """Establish the density and viscosity parameters.

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
    1. References
        CO2 Density from NIST REFPROP software, based on
          Span & Wagner (1996). Range 0.1 to 60 MPa and 1 to 180 oC.
        CO2 viscosity from NIST REFPROP software, based on
          Fenghour et al (1998).  Range 0.1 to 60 MPa and 1 to 180 oC.
    """
    # Load density and viscosity LUTs.
    co2_density_lut = define_lookup_array(scfg.CO2_DENSITY_FILE)
    co2_viscosity_lut = define_lookup_array(scfg.CO2_VISCOSITY_FILE)

    # Use 2D cubic interpolation to estimate values from table.
    z_density = conduct_2d_interp(fdata.CO2_DENSITY_PRESSURE,
                                  fdata.CO2_DENSITY_TEMP,
                                  co2_density_lut,
                                  pressure, temperature)

    z_viscosity = conduct_2d_interp(fdata.CO2_VISCOSITY_PRESSURE,
                                    fdata.CO2_VISCOSITY_TEMP,
                                    co2_viscosity_lut,
                                    pressure, temperature)

    return z_density, z_viscosity


# noinspection DuplicatedCode,DuplicatedCode
def brine_density_property(pressure, temperature, salinity=0.0):
    """Estimate brine density using a 3D lookup table.

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
    Reference:
    1. Temperature/pressure/salinity-dependent values based on
        Sun et al. (2008) + density of pure H2O provided by
        NIST (2010) (Tables from Sun provided by J Morrell;
        NIST data based on Wagner & Pruss [1995]).
        EOS valid to 100 MPa, 374 oC and salinity of
        ~80,000 mg/kg (ratio of mass of solute / mass of solution,
        or ppm).
    """
    # Identify interpolation method & index in 3D array:
    #   -> It is a function of (salinity, pressure, temperature).
    #   -> method = 2 for 2D interpolation - cubic.
    #   -> method = 3 for 3D interpolation - cubic + linear (salinity).
    #   -> dot = position in array.
    dot, method = salinity_locate(salinity, fdata.BRINE_DENSITY_SALINITY)

    # Catch error, if it occurs.
    if method not in (2, 3):
        fileop.opx_problem("Failure in Interpolation During Brine Density "
                           + "Calculation.")

    # Debug print.
    if ECHO:
        print("\n  For Brine Density")
        debug_shape(method, dot,
                    fdata.BRINE_DENSITY_TEMP,
                    fdata.BRINE_DENSITY_PRESSURE,
                    fdata.BRINE_DENSITY_SALINITY)
        print(f'   ** Salinity in ppm: = {salinity}')

    # STEP #1: 2D Interpolation  ----------------------------------------------

    # Establish variables and arrays for operation.
    # i_rows = len(fdata.BRINE_DENSITY_PRESSURE).
    # i_cols = len(fdata.BRINE_DENSITY_TEMP).

    # Load density LUT.
    brine_density_lut = define_lookup_array(scfg.BRINE_DENSITY_FILE)

    # Establish 2D subarray for density at low index salinity.
    density_low = brine_density_lut[dot, :, :]

    # Conduct 2D interpolation with low array.
    z_density = conduct_2d_interp(fdata.BRINE_DENSITY_PRESSURE,
                                  fdata.BRINE_DENSITY_TEMP,
                                  density_low, pressure, temperature)

    # STEP #2: 3D Interpolation  ----------------------------------------------

    # Perform 3d interpolation across salinities, if required (method = 3).
    if method == 3:

        # Establish subarray for density at high index salinity.
        density_high = brine_density_lut[(dot + 1), :, :]

        # Conduct 2D interpolation with low array.
        z_density_up = conduct_2d_interp(fdata.BRINE_DENSITY_PRESSURE,
                                         fdata.BRINE_DENSITY_TEMP,
                                         density_high, pressure, temperature)

        # Conduct linear interpolations between values.
        low_salinity = fdata.BRINE_DENSITY_SALINITY[dot]
        high_salinity = fdata.BRINE_DENSITY_SALINITY[dot + 1]
        z_density = linear_interpolate(salinity,
                                       low_salinity, high_salinity,
                                       z_density, z_density_up)

    return z_density


# noinspection DuplicatedCode,DuplicatedCode
def brine_viscosity_property(pressure, temperature, salinity=0.0):
    """Establish viscosity using a 3D lookup table.

    Parameters
    ----------
    pressure = (float) pressure of brine (MPa)
    temperature = (float) temperature of brine (oC)
    salinity = (float) NaCl in ppm (to be converted to molal)

    Returns
    -------
    z_viscosity = brine viscosity (Pa-s)

    Notes
    -----
    1. Reference:
        Temperature/pressure/salinity-dependent viscosity of a
        saline solution based on Mao and Duan (2009) in terms
        of molality + viscosity of pure H2O provided by
        NIST (2010); NIST data based on Huber et al. [2009]
        (Tables from Mao & Duan data provided by J Morrell).
        EOS valid to 100 MPa, 350 oC and salinity to 6 M.
    """
    # Identify interpolation method & index in 3D array:
    #   -> It is a function of (salinity, pressure, temperature)
    #   -> method = 2 for 2D interpolation - cubic
    #   -> method = 3 for 3D interpolation - cubic + linear (on saline)
    #   -> dot = position in array

    #  Note: lookup table is in Molal, not ppm!
    sal_molal = convert_salinity(salinity)
    dot, method = salinity_locate(sal_molal, fdata.BRINE_VISCOSITY_SALINITY)

    # Debug.
    if ECHO:
        print("\n  For Brine Viscosity")
        debug_shape(method, dot,
                    fdata.BRINE_VISCOSITY_TEMP,
                    fdata.BRINE_VISCOSITY_PRESSURE,
                    fdata.BRINE_VISCOSITY_SALINITY)
        print(f"   ** Salinity in Molal: = {sal_molal}")

    # Catch error, if it occurs.
    if method not in (2, 3):
        fileop.opx_problem("Failure in Interpolation During Brine Viscosity "
                           + "Calculation.")

    # STEP #1: 2D Interpolation  ----------------------------------------------

    # Establish variables for operation.
    # i_rows = len(fdata.BRINE_VISCOSITY_PRESSURE).
    # i_cols = len(fdata.BRINE_VISCOSITY_TEMP).

    # Load viscosity LUT.
    brine_viscosity_lut = define_lookup_array(scfg.BRINE_VISCOSITY_FILE)

    # Establish subarray for viscosity at low index salinity.
    viscosity_low = brine_viscosity_lut[dot, :, :]

    # Conduct 2D interpolation with low array.
    z_viscosity = conduct_2d_interp(fdata.BRINE_VISCOSITY_PRESSURE,
                                    fdata.BRINE_VISCOSITY_TEMP,
                                    viscosity_low, pressure, temperature)

    # STEP #2: 3D Interpolation  ----------------------------------------------

    # Perform 3d interpolation across salinities, if required (method = 3).
    if method == 3:

        # Establish subarray for viscosity at high index salinity.
        # viscosity_high = np.zeros((i_cols, i_rows)).
        viscosity_high = brine_viscosity_lut[(dot + 1), :, :]

        # Conduct 2D interpolation with low array.
        z_viscosity_up = conduct_2d_interp(fdata.BRINE_VISCOSITY_PRESSURE,
                                           fdata.BRINE_VISCOSITY_TEMP,
                                           viscosity_high, pressure,
                                           temperature)

        # Conduct linear interpolations between values.
        low_salinity = fdata.BRINE_VISCOSITY_SALINITY[dot]
        high_salinity = fdata.BRINE_VISCOSITY_SALINITY[dot + 1]
        z_viscosity = linear_interpolate(sal_molal,
                                         low_salinity, high_salinity,
                                         z_viscosity, z_viscosity_up)

    return z_viscosity


# noinspection DuplicatedCode,DuplicatedCode
def co2_solubility_property(pressure, temperature, salinity=0.0):
    """Estimate CO2 solubility using a 3D lookup table.

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
    Reference:
        Temperature/pressure/salinity-dependent values based on
        programming the equation of state correlations as reported
        by Duan et al. (2006) and computing required values by
        Jason Monnell, Research Assistant Professor, University of Pittsburgh.
    """
    # Convert ppm to mol for analysis.
    sal_molal = convert_salinity(salinity)

    # Identify interpolation method & index in 3D array:
    #   -> It is a function of (salinity, pressure, temperature).
    #   -> method = 2 for 2D interpolation - cubic.
    #   -> method = 3 for 3D interpolation - cubic + linear (salinity).
    #   -> dot = position in array.
    dot, method = salinity_locate(sal_molal, fdata.CO2_SOLUBILITY_SALINITY)

    # Catch error, if it occurs.
    if method not in (2, 3):
        fileop.opx_problem("Failure in Interpolation During CO2 Solubility "
                           + "Calculation.")

    # Provide debug print.
    if ECHO:
        print("\n  For CO2 solubility")
        debug_shape(method, dot,
                    fdata.CO2_SOLUBILITY_TEMP,
                    fdata.CO2_SOLUBILITY_PRESSURE,
                    fdata.CO2_SOLUBILITY_SALINITY)
        print(f"   ** Salinity in mol: = {sal_molal}")

    # STEP #1: 2D Interpolation  ----------------------------------------------

    # Establish variables and arrays for operation.
    # i_rows = len(fdata.CO2_SOLUBILITY_PRESSURE).
    # i_cols = len(fdata.CO2_SOLUBILITY_TEMP).

    # Read LUT and create NPY file.
    #  co2_solubility_lut = sol_data.CO2_SOLUBILITY_TABLE.
    #  target = "data_co2_solubility".
    #  file_path = fileop.find_local_path(target).
    #  np.save(file_path, co2_solubility_lut).

    #  Load solubility LUT.
    co2_solubility_lut = define_lookup_array(scfg.CO2_SOLUBILITY_FILE)

    # Establish 2D subarray for solubility at low index salinity.
    solubility_low = co2_solubility_lut[dot, :, :]

    # Conduct 2D interpolation with low array.
    z_solubility = conduct_2d_interp(fdata.CO2_SOLUBILITY_TEMP,
                                     fdata.CO2_SOLUBILITY_PRESSURE,
                                     solubility_low, temperature, pressure)

    # STEP #2: 3D Interpolation  ----------------------------------------------

    # Perform 3d interpolation across salinities, if required (method = 3).
    if method == 3:

        # Establish subarray for solubility at high index salinity.
        solubility_high = co2_solubility_lut[(dot + 1), :, :]

        # Conduct 2D interpolation with low array.
        z_solubility_up = conduct_2d_interp(fdata.CO2_SOLUBILITY_TEMP,
                                            fdata.CO2_SOLUBILITY_PRESSURE,
                                            solubility_high, temperature,
                                            pressure)

        # Conduct linear interpolations between values.
        low_salinity = fdata.CO2_SOLUBILITY_SALINITY[dot]
        high_salinity = fdata.CO2_SOLUBILITY_SALINITY[dot + 1]
        z_solubility = linear_interpolate(sal_molal,
                                          low_salinity, high_salinity,
                                          z_solubility, z_solubility_up)

    return z_solubility


def ave_properties(reservoir_list, aquifer_list):
    """Compute average of two lists and update fluid parameters.

    Parameters
    ----------
    sub_list = (list) parameters at inject depth (see Notes)
    aquifer_list = (list) parameters at aquifer base depth

    Returns
    -------
    property_list = (list) averaged values

    Notes
    -----
    1. list properties (sub_ and aquifer_):
        list[0] = co2_density = new CO2 density
        list[1] = co2_viscosity = new CO2 viscosity
        list[2] = brine_density = new brine density
        list[3] = brine_viscosity = new brine viscosity
        list[4] = co2solubility = new CO2 solubility
    """
    # Average the elements of the two list values.
    zipped_lists = zip(reservoir_list, aquifer_list)
    property_list = [(x + y)/2.0 for (x, y) in zipped_lists]

    return property_list


def interpolate_fluid_properties(temperature, pressure, salinity):
    """Interpolate to get fluid properties for brine & CO2.

    Parameters
    ----------
    temperature = (float) local temperature (oC)
    pressure = (float) local pressure (Pa)
    salinity = (float) local salinity (ppm)

    Returns
    -------
    property_list = (list) interpolated properties

    Note
    ----
    1. Values in lookup tables are in MPa! need to convert pressure in Pa.
    """
    # Convert pressure to MPa and set other computation parameters.
    pressure = pressure * funit.pa_to_mpa()

    # Lookup CO2 density, viscosity from lookup tables.
    co2_density, co2_viscosity = co2_properties(pressure, temperature)
    brine_density = brine_density_property(pressure, temperature, salinity)
    brine_viscosity = brine_viscosity_property(pressure, temperature, salinity)
    co2_solubility = co2_solubility_property(pressure, temperature, salinity)

    # Create list for transfer for reset.
    property_list = [co2_density, co2_viscosity, brine_density,
                     brine_viscosity, co2_solubility]

    # Debug print for check.
    if ECHO:
        # Debug - print variables for CO2 and brine.
        echo_txt = "\n   ** CO2 density and viscosity = "
        echo_txt += f"{co2_density:.3f} kg/m3 "
        echo_txt += f"and {co2_viscosity:.5e} Pa-s"
        print(echo_txt)
        echo_txt = ("   ** brine density and viscosity = "
                    f"{brine_density:.3f} kg/m3 ")
        echo_txt += f"and {brine_viscosity:.5e} Pa-s"
        print(echo_txt)
        print(f"   ** CO2 solubility = {co2_solubility:.3f} kg/kg")

    return property_list


def reset_fluid_properties(fault_controls, new_list):
    """Update/reset fluid properties from new list.

    Parameters
    ----------
    fault_controls = (dict) fault controls
    new_list = (list) new fluid properties
    aqui_list = (list) new aquifer properties

    Returns
    -------
    fault_controls = (dict) dictionary of fault controls (updated)
    """
    # Reset fluid properties.
    fault_controls['co2_density'] = new_list[0]
    fault_controls['co2_viscosity'] = new_list[1]
    fault_controls['brine_density'] = new_list[2]
    fault_controls['brine_viscosity'] = new_list[3]
    fault_controls['co2_solubility'] = new_list[4]

    return fault_controls


def reset_aqui_properties(fault_controls, aqui_list):
    """Reset fluid properties for aquifer.

    Parameters
    ----------
    fault_controls = (dict) fault controls
    aqui_list = (list) new aquifer properties

    Returns
    -------
    fault_controls = (dict) dictionary of fault controls (updated)
    """
    fault_controls['aqui_co2_density'] = aqui_list[0]
    fault_controls['aqui_co2_viscosity'] = aqui_list[1]
    fault_controls['aqui_brine_density'] = aqui_list[2]
    fault_controls['aqui_brine_viscosity'] = aqui_list[3]

    return fault_controls


def debug_shape(bmethod, dotted, temperature, pressure, salinity):
    """Print shape features for debug.

    Parameters
    ----------
    bmethod = (int) interpolation method (=1,2,3)
    dotted = (int) array position
    temperature = (float) temperature (oC)
    pressure = (float) reference pressure (Pa)
    salinity = (float) salinity (ppm)

    Returns
    -------
    N/A
    """
    # Convert with units.
    print(f'   ** Interpolation Type: {bmethod}')
    print(f'   ** Salinity Index: {dotted}')

    temp_val = str(temperature.shape).replace(',', '')
    print(f'   ** Shape of Temperature Vector: {temp_val}')

    temp_val = str(pressure.shape).replace(',', '')
    print(f'   ** Shape of Pressure Vector: {temp_val}')

    temp_val = str(salinity.shape).replace(',', '')
    print(f'   ** Shape of Salinity Array: {temp_val}')

    # return None


#
# -----------------------------------------------------------------------------
# End of module
