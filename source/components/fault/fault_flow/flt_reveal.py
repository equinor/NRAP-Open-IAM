#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions to write data to summary file.

Author: Ernest N. Lindner
Date: 08/23/2022

Module Name
    flt_reveal

Contents (18)
    define_units_for_output(si_units)
    write_iam_parameters(fault_controls, zum)
    write_fault_controls(fault_controls, zum)
    write_description_parameters(fault_controls, uts, zum)
    write_field_parameters(fault_controls, uts, zum)
    write_aperture_parameters(fault_controls, uts, zum)
    write_aqui_fluid_values(fault_controls, uts, zum)
    write_deep_fluid_values(fault_controls, uts, zum)
    write_relative_perm(fault_controls, uts, zum)
    write_stress_parameters(fault_controls, uts, zum)
    ---
    write_profile_parameters(fault_controls, uts, zum)
    write_profile_simple(fault_controls, uts, zum)
    write_profile_complex(fault_controls, uts, zum)
    write_pressure_aperture_parameters(fault_controls, uts, zum)
    write_fault_parameters(fault, zum)
    write_co2_results(sim_flux_list, fault_controls, uts, zum)
    write_last(zum)
    write_summary(fault_controls, sim_flux_list, fault)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import numpy as np              # For array operations

import flt_file as fileop       # For paths and error reports
import flt_config as scfg       # IO directory and file names
import flt_units as funits      # For unit conversion

# Constants for file handling & writing to Summary File
ZIPX_IN = "  --> "              # Parameter start
SPEC_IN = "    >>> "            # Message start
BLK_IN = "\n  > "               # Block start
ECHO = False                    # Additional printout control
FAULT_LIMIT = 99.9              # Fault probability limit


def define_units_for_output(si_units):
    """Define units for output - as strings.

    Parameters
    ----------
    si_units = (bool) control to convert English to SI units

    Returns
    -------
    uts[] = (list of str) units for printing
    """
    if si_units:
        # Values are in SI; define units for printout.
        # set unit labels for printout in SI.
        uts = ['m', 'MPa', 'kg/m3', 'Pa-s', 'mol/kg', 'm2', 'tonne',
               'oC', 'mm', 'ppm']
    else:
        # values are in English/US units; define units for printout.
        uts = ['ft', 'psi', 'pcf', 'cP', 'mol/lb', 'ft2', 'short ton',
               'oF', 'ín', 'ppm']

    return uts


def write_iam_parameters(fault_controls, zum):
    """Write basic model info to file named "zum".

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    zum = (str) file label

    Returns
    -------
    N/A
    """
    # Write time stamp for start of the run.
    print('\n  Fault_Flo Summary', file=zum)
    now = fault_controls['record_date']
    clump = ZIPX_IN + 'Simulation Run Date: %a %d %b %Y - Time: %H:%M'
    stamp = now.strftime(clump)
    print(stamp, file=zum)

    # Write run details.
    print(BLK_IN + 'Model Parameters', file=zum)
    print(ZIPX_IN + f'Title = {fault_controls.get("title"):<} ', file=zum)

    print(ZIPX_IN + 'Start Time of Analysis     = '
          f'{fault_controls.get("start_time"):<.1f} years',
          file=zum)
    print(ZIPX_IN + 'End Time of Analysis       = '
          f'{fault_controls.get("end_time"):<.1f} years', file=zum)
    print(ZIPX_IN + 'End Time of Injection      = '
          f'{fault_controls.get("inject_end"):<.1f} years', file=zum)
    print(ZIPX_IN + 'Number of Realizations     = '
          f'{fault_controls.get("realizations"):<5}', file=zum)
    print(ZIPX_IN + 'Time Steps for Run         = '
          f'{fault_controls.get("time_points"):<5}', file=zum)

    if fault_controls.get('time_input'):
        print(SPEC_IN + 'Time Steps are Input From File', file=zum)
    else:
        print(SPEC_IN + 'Uniform Time Steps are Assumed', file=zum)

    input_form = (fault_controls.get('units')).upper()  # upper case
    print(ZIPX_IN + f'Units for Output = {input_form:<}', file=zum)

    # return None


def write_fault_controls(fault_controls, zum):
    """Write controls to file named "zum".

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    zum = (str) file label

    Returns
    -------
    N/A
    """
    # Write header for fault controls.
    print(BLK_IN + 'Fault Controls', file=zum)

    # profile type.
    input_form = int(fault_controls.get('profile_type'))
    if input_form == 0:
        print(SPEC_IN + 'Single Geothermal Profile is Used (profile Type = 0)',
              file=zum)
    else:
        print(SPEC_IN + 'Near Surface Profile Considered',
              file=zum)

    # File input statements.
    if fault_controls.get('aperture_approach'):
        print(SPEC_IN + 'Aperture Values are Input from File.', file=zum)
    else:
        print(SPEC_IN + 'Aperture Values are Computed.', file=zum)

    if fault_controls.get('strike_approach'):
        print(SPEC_IN + 'Fault Strike is Varied in Simulations', file=zum)
    else:
        print(SPEC_IN + 'Fault Strike is Constant', file=zum)

    if fault_controls.get('dip_approach'):
        print(SPEC_IN + 'Fault Dip is Varied in Simulations', file=zum)
    else:
        print(SPEC_IN + 'Fault Dip is Constant', file=zum)

    if fault_controls.get('pressure_approach'):
        print(SPEC_IN + 'Fault Aperture is Varied with Pressure', file=zum)
    else:
        print(SPEC_IN + 'Fault Aperture is Constant', file=zum)

    # return None


def write_description_parameters(fault_controls, uts, zum):
    """Write description and orientation data to file named "zum".

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    uts = (str, array) units for output
    zum = (str) file label

    Returns
    -------
    N/A
    """
    # Define parameters for either unit option.
    length = fault_controls.get('length')
    start_x = fault_controls.get('x_start')
    start_y = fault_controls.get('y_start')
    shortest = fault_controls.get('shortest_distance')

    if not fault_controls['use_si']:
        # US units.
        length *= funits.meters_to_feet()
        start_x *= funits.meters_to_feet()
        start_y *= funits.meters_to_feet()
        shortest *= funits.meters_to_feet()

    #  *******************************************************************
    # Write parameters on fault core to summary.
    print(BLK_IN + 'Fault Core Description', file=zum)

    # Write Yes or No responses for each approach in English.
    fault_exists = fault_controls.get('fault_probability')
    if fault_exists < FAULT_LIMIT:
        print(SPEC_IN + 'Fault Existence Varied in Simulations', file=zum)

        # Write fault probability.
        if fault_exists < FAULT_LIMIT:
            print(ZIPX_IN + 'Probability of a Fault     = '
                  + f'{fault_exists:.1f}%',
                  file=zum)
    else:
        print(SPEC_IN + 'A Fault Present in All Simulations', file=zum)

    # Write strike parameters.
    print(ZIPX_IN + 'Fault Strike - Mean        = '
          f'{fault_controls.get("strike_mean"):.1f} deg', file=zum)
    if fault_controls['strike_approach']:
        print(ZIPX_IN + 'Fault Strike - Spread      = '
              f'{fault_controls.get("strike_sig"):.1f} deg', file=zum)
        print(ZIPX_IN + 'Fault Strike - Kappa       = '
              f'{fault_controls.get("strike_disperse"):.1f}', file=zum)

    # Write dip parameters.
    print(ZIPX_IN + 'Fault Dip                  = '
          f'{fault_controls.get("dip"):.1f} deg', file=zum)
    if fault_controls['dip_approach']:
        print(ZIPX_IN + 'Fault Dip - Std. Deviation = '
              f'{fault_controls.get("dip_std"):.1f} deg', file=zum)
    print(ZIPX_IN + 'Fault Length               = '
          f'{length:.1f} {uts[0]}', file=zum)
    print(ZIPX_IN + 'At Surface - X Coordinate  = '
          f'{start_x:.1f} {uts[0]}', file=zum)
    print(ZIPX_IN + 'At Surface - Y Coordinate  = '
          f'{start_y:.1f} {uts[0]}', file=zum)
    print(ZIPX_IN + 'Number of Segments         = '
          f'{fault_controls.get("n_plates"):<5}', file=zum)

    # Write orientation parameters with header.
    print(BLK_IN + 'Orientation Factors - Computed', file=zum)
    print(ZIPX_IN + 'Fault Inclination          = ', end='', file=zum)
    print(fault_controls.get('fault_inclined'), file=zum)

    print(ZIPX_IN + 'Flow Correction            = '
          f'{fault_controls.get("fault_orient_effect"):.2f}', file=zum)
    print(ZIPX_IN + 'Shortest Distance to Well  = '
          f'{shortest:.5e} {uts[0]}', file=zum)
    print(ZIPX_IN + 'Angle Between Vectors      = '
          f'{fault_controls.get("angle_between_vectors"):.1f} deg.', file=zum)

    # return None


def write_field_parameters(fault_controls, uts, zum):
    """Write field data to file named "zum".

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    uts = (str, array) units for output
    zum = (str) file label

    Returns
    -------
    N/A
    """
    # Define parameters for units.
    aquifer_deep = fault_controls.get('aquifer_depth')
    aquifer_press = fault_controls.get('aquifer_pressure')
    aquifer_temp = fault_controls.get('aquifer_temperature')
    inject_deep = fault_controls.get('inject_depth')
    inject_temp = fault_controls.get('inject_temperature')
    field_press = fault_controls.get('field_pressure')
    inject_press = fault_controls.get('inject_pressure')
    final_press = fault_controls['final_pressure']
    x_inject = fault_controls.get('inject_x')
    y_inject = fault_controls.get('inject_y')
    if not fault_controls['use_si']:
        # US units.
        aquifer_deep *= funits.meters_to_feet()
        aquifer_press *= funits.pa_to_psi()
        aquifer_temp = funits.celsius_to_fahrenheit(aquifer_temp)
        inject_deep *= funits.meters_to_feet()
        inject_temp = funits.celsius_to_fahrenheit(inject_temp)
        field_press *= funits.pa_to_psi()
        inject_press *= funits.pa_to_psi()
        final_press *= funits.pa_to_psi()
        x_inject *= funits.meters_to_feet()
        y_inject *= funits.meters_to_feet()
    else:
        aquifer_press *= funits.pa_to_mpa()
        field_press *= funits.pa_to_mpa()
        inject_press *= funits.pa_to_mpa()
        final_press *= funits.pa_to_mpa()

    #  *********************************************************
    # Write header for field descriptors.
    print(BLK_IN + 'Field Conditions', file=zum)

    # Write numeric data for field conditions.
    print(ZIPX_IN + 'Aquifer Depth              = '
          f'{aquifer_deep:.0f} {uts[0]}', file=zum)
    print(ZIPX_IN + 'Aquifer Pressure           = '
          f'{aquifer_press:.4e} {uts[1]}', file=zum)
    print(ZIPX_IN + 'Aquifer Temperature        = '
          f'{aquifer_temp:.0f} {uts[7]}', file=zum)
    print(ZIPX_IN + 'Reservoir Depth            = '
          f'{inject_deep:.0f} {uts[0]}', file=zum)
    print(ZIPX_IN + 'Reservoir Temperature      = '
          f'{inject_temp:.0f} {uts[7]}', file=zum)
    print(ZIPX_IN + 'Initial Reservoir Pressure = '
          f'{field_press:.4e} {uts[1]}', file=zum)
    print(ZIPX_IN + 'Ave. Injection Pressure    = '
          f'{inject_press:.4e} {uts[1]}', file=zum)
    print(ZIPX_IN + 'Final Pressure             = '
          f'{final_press:.4e} {uts[1]}', file=zum)
    print(ZIPX_IN + 'X Coordinate of Injection  = '
          f'{x_inject:.1f} {uts[0]}', file=zum)
    print(ZIPX_IN + 'Y Coordinate of Injection  = '
          f'{y_inject:.1f} {uts[0]}', file=zum)

    # return None


def write_aperture_parameters(fault_controls, uts, zum):
    """Write aperture and correction data to file named "zum".

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    uts = (str, array) units for output
    zum = (str) file label

    Returns
    -------
    N/A
    """
    # Define parameters for units.
    mean = fault_controls.get('aperture_mean')
    standard = fault_controls.get('aperture_std')
    min_aperture = fault_controls.get('aperture_min')
    max_aperture = fault_controls.get('aperture_max')
    if not fault_controls['use_si']:
        # US units.
        mean *= funits.mm_to_inch()
        standard *= funits.mm_to_inch()
        min_aperture *= funits.mm_to_inch()
        max_aperture *= funits.mm_to_inch()

    #  ************************************************************************
    # print aperture variability and correction parameters.
    print(BLK_IN + 'Aperture & Correction Parameters', file=zum)

    if fault_controls.get('aperture_approach'):
        print(SPEC_IN + 'Aperture Option: Values Input from File.', file=zum)
    elif not fault_controls.get('vary_aperture'):
        print(SPEC_IN + 'Uniform Aperture - All Set to Mean Value.', file=zum)
        print(ZIPX_IN + 'Mean Total Aperture        = '
              f'{mean:<.3f} {uts[8]}', file=zum)
    else:
        print(SPEC_IN + 'Computed with Censored Log-normal Distribution.',
              file=zum)
        print(ZIPX_IN + 'Mean Aperture              = '
              f'{mean:.5f} {uts[8]}', file=zum)
        print(ZIPX_IN + 'Std. Dev. Aperture         = '
              f'{standard:.5f} {uts[8]}', file=zum)
        print(ZIPX_IN + 'Minimum Aperture           = '
              f'{min_aperture:.5f} {uts[8]}', file=zum)
        print(ZIPX_IN + 'Maximum Aperture           = '
              f'{max_aperture:.5f} {uts[8]}', file=zum)

        # Correction parameters.
        print(ZIPX_IN + 'Non-Isothermal Correction  = '
              f'{fault_controls.get("state_correction"):.2f}', file=zum)
        print(ZIPX_IN + 'SGR Value                  = '
              f'{fault_controls.get("sgr"):.0f}%', file=zum)
        print(ZIPX_IN + 'SGR Factor                 = '
              f'{fault_controls.get("sgr_perm_correction"):.2e}', file=zum)

    # return None


def write_aqui_fluid_values(fault_controls, uts, zum):
    """Write aquifer fluid parameters to file named "zum".

    Parameters
    ----------
    fault_controls = (dict) fault parameters
    uts = (list - str) = unit abbreviation
    zum = (str) file label

    Returns
    -------
    None
    """
    # Get numeric values for print out.
    aq_co2_density = fault_controls.get('aqui_co2_density')
    aq_co2_viscosity = fault_controls.get('aqui_co2_viscosity')
    aq_br_density = fault_controls.get('aqui_brine_density')
    aq_br_viscosity = fault_controls.get('aqui_brine_viscosity')

    if not fault_controls['use_si']:
        # US units.
        aq_co2_density *= funits.kilo_convert_density()
        aq_co2_viscosity *= funits.pa_s_to_centipoise()
        aq_br_density *= funits.kilo_convert_density()
        aq_br_viscosity *= funits.pa_s_to_centipoise()

    #  ***********************************************************
    # Write fluid parameters.

    # Print calculated fluid parameters at aquifer.
    print(BLK_IN + 'Aquifer Fluid Parameters', file=zum)
    print(ZIPX_IN + 'CO2 Density                = '
          f'{aq_co2_density:.0f} {uts[2]}', file=zum)
    print(ZIPX_IN + 'CO2 Viscosity              = '
          f'{aq_co2_viscosity:.4e} {uts[3]}', file=zum)
    print(ZIPX_IN + 'Brine Density              = '
          f'{aq_br_density:.0f} {uts[2]}', file=zum)
    print(ZIPX_IN + 'Brine Viscosity            = '
          f'{aq_br_viscosity:.4e} {uts[3]}', file=zum)

    # return None


def write_deep_fluid_values(fault_controls, uts, zum):
    """Write deep fluid parameters to file named "zum".

    Parameters
    ----------
    fault_controls = (dict) fault parameters
    uts = (list - str) = unit abbreviation
    zum = (str) file label

    Returns
    -------
    None
    """
    # Get numeric values for print out.
    co2_density = fault_controls.get('co2_density')
    co2_viscosity = fault_controls.get('co2_viscosity')
    br_density = fault_controls.get('brine_density')
    br_viscosity = fault_controls.get('brine_viscosity')
    solubility = fault_controls.get('co2_solubility')

    if not fault_controls['use_si']:
        # US units.
        co2_density *= funits.kilo_convert_density()
        co2_viscosity *= funits.pa_s_to_centipoise()
        br_density *= funits.kilo_convert_density()
        br_viscosity *= funits.pa_s_to_centipoise()
        solubility *= funits.mol_per_kg_to_mol_per_lb()

    #  ***********************************************************
    # Write fluid parameters.
    print(BLK_IN + 'Fluid Parameters', file=zum)
    if fault_controls.get('interpolate_approach'):
        print(SPEC_IN + 'Fluid Parameters are Interpolated.',
              file=zum)
    else:
        print(SPEC_IN + 'Default Fluid Parameters were Used.',
              file=zum)

    # Print calculated fluid parameters at depth.
    print(BLK_IN + 'Fluid Parameters at Injection Depth', file=zum)
    print(ZIPX_IN + 'Salinity                   = '
          f'{fault_controls.get("salinity"):.0f} {uts[9]}', file=zum)
    print(ZIPX_IN + 'CO2 Density                = '
          f'{co2_density:.0f} {uts[2]}', file=zum)
    print(ZIPX_IN + 'CO2 Viscosity              = '
          f'{co2_viscosity:.4e} {uts[3]}', file=zum)
    print(ZIPX_IN + 'Brine Density              = '
          f'{br_density:.0f} {uts[2]}', file=zum)
    print(ZIPX_IN + 'Brine Viscosity            = '
          f'{br_viscosity:.4e} {uts[3]}', file=zum)
    print(ZIPX_IN + 'CO2 Solubility             = '
          f'{solubility:.4e} {uts[4]}', file=zum)

    # return None


def write_relative_perm(fault_controls, uts, zum):
    """Write relative permeability controls to a file named "zum".

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    uts = (str, array) units for output
    zum = (str) file label

    Returns
    -------
    N/A

    Notes
    -----
    1. Models BC and LET.
    2. LET not implemented in this version.
    """
    # Define stress values in MPa.
    cap_pressure = fault_controls.get('entry_pressure')
    if fault_controls['use_si']:
        cap_pressure *= funits.psi_to_mpa()
    else:
        cap_pressure *= funits.pa_to_mpa()

    # Define pressure for LET model.
    max_press = 0.0
    if 'LET' in fault_controls.get('relative_model'):
        max_press = fault_controls.get('max_capillary')
        if fault_controls['use_si']:
            max_press *= funits.psi_to_mpa()
        else:
            max_press *= funits.pa_to_mpa()

    #  **********************************************************
    # print relative permeability parameters.
    print(BLK_IN + 'Relative Permeability Flow Limits', file=zum)

    # Write general limits on relative permeability model.
    print(ZIPX_IN + 'Residual Brine Saturation  = '
          f'{fault_controls.get("resid_brine"):.2f}', file=zum)
    print(ZIPX_IN + 'Residual CO2 Saturation    = '
          f'{fault_controls.get("resid_co2"):.2f}', file=zum)
    print(ZIPX_IN + 'Nonwetting/Wetting Ratio   = '
          f'{fault_controls.get("perm_ratio"):.2f}', file=zum)
    print(ZIPX_IN + 'Ave. Entry Pressure        = '
          f'{cap_pressure:.4e} {uts[1]}', file=zum)

    #  *********************************************************
    # Write relative permeability model and parameters.
    print(BLK_IN + 'Relative Permeability Model', file=zum)
    # Write model type.
    m_type = fault_controls.get('relative_model')
    if 'LET' in m_type:
        print(SPEC_IN + 'Relative Permeability = L-E-T Model',
              file=zum)
    else:
        print(SPEC_IN + 'Relative Permeability = ' +
              'Brooks-Corey Model', file=zum)

    # Write LET & BC relative permeability parameters.
    print(BLK_IN + 'Relative Permeability Parameters', file=zum)
    if "LET" in fault_controls.get('relative_model'):
        print(ZIPX_IN + '"L" Wetting Parameter      = '
              f'{fault_controls.get("l_wetting"):.2f}', file=zum)
        print(ZIPX_IN + '"E" Wetting Parameter      = '
              f'{fault_controls.get("e_wetting"):.2f}', file=zum)
        print(ZIPX_IN + '"T" Wetting Parameter      = '
              f'{fault_controls.get("t_wetting"):.2f}', file=zum)
        print(ZIPX_IN + '"L" Nonwetting Parameter   = '
              f'{fault_controls.get("l_nonwet"):.2f}', file=zum)
        print(ZIPX_IN + '"E" Nonwetting Parameter   = '
              f'{fault_controls.get("e_nonwet"):.2f}', file=zum)
        print(ZIPX_IN + '"T" Nonwetting Parameter   = '
              f'{fault_controls.get("t_nonwet"):.2f}', file=zum)
    else:
        print(ZIPX_IN + 'Lambda Parameter for BC    = '
              f'{fault_controls.get("zeta"):.2f}', file=zum)

    # Write LET model capillary parameters.
    if 'LET' in fault_controls.get('relative_model'):
        print(BLK_IN + 'Capillary Parameters', file=zum)
        print(ZIPX_IN + '"L" Capillary Parameter    = '
              f'{fault_controls.get("l_capillary"):.2f}', file=zum)
        print(ZIPX_IN + '"E" Capillary Parameter    = '
              f'{fault_controls.get("e_capillary"):.2f}', file=zum)
        print(ZIPX_IN + '"T" Capillary Parameter    = '
              f'{fault_controls.get("t_capillary"):.2f}', file=zum)
        print(ZIPX_IN + 'Maximum Capillary Pressure = '
              f'{max_press:.4e} {uts[1]}', file=zum)

    # return None


def write_stress_parameters(fault_controls, uts, zum):
    """Write in situ stress parameters to file named "zum".

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    uts = (list - str) = unit abbreviation
    zum = (str) file label

    Returns
    -------
    N/A
    """
    # Define default in situ stress values.
    maximal = fault_controls.get('max_horizontal')
    minimal = fault_controls.get('min_horizontal')
    if not fault_controls['use_si']:
        # US units.
        maximal *= funits.pa_to_psi()
        minimal *= funits.pa_to_psi()
    else:
        # Convert metric pressure values to MPa.
        maximal *= funits.pa_to_mpa()
        minimal *= funits.pa_to_mpa()

    #  ************************************************************************
    # Write values.
    print(BLK_IN + 'In-Situ, Secondary Principal Horizontal Stress Parameters',
          file=zum)
    print(ZIPX_IN + 'Maximum Horizontal Stress  = '
          f'{maximal:.2e} {uts[1]}', file=zum)
    print(ZIPX_IN + 'Minimum Horizontal Stress  = '
          f'{minimal:.2e} {uts[1]}', file=zum)
    print(ZIPX_IN + 'Strike of Maximum Stress   = '
          f'{fault_controls.get("max_trend"):.1f} deg', file=zum)

    # return None


def write_profile_parameters(fault_controls, uts, zum):
    """Write CO2 transition parameters to file named "zum".

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    uts = (list - str) = unit abbreviation
    zum = (str) file label

    Returns
    -------
    N/A
    """
    # print separator.
    print('\n  ', end='', file=zum)
    print('*' * 60, file=zum)

    # Write values for transition.
    print(BLK_IN + 'CO2 Transition Parameters', file=zum)

    if fault_controls['profile_type'] == 0:
        # For simple profile - model #0.
        print(SPEC_IN + 'Simple Supercritical Flow Model', file=zum)
        write_profile_simple(fault_controls, uts, zum)
    elif fault_controls['profile_type'] == 1:
        # For complex profile - model #1.
        print(SPEC_IN + 'Complex Profile Flow Model', file=zum)
        write_profile_complex(fault_controls, uts, zum)
    else:
        # For disjointed profile - model #2.
        print(SPEC_IN + 'Disjointed Profile Flow Model', file=zum)
        #  -- > Printout same as Model 1.
        write_profile_complex(fault_controls, uts, zum)

    # return None


def write_profile_simple(fault_controls, uts, zum):
    """Write simple CO2 transition parameters to file named "zum".

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    uts = (list - str) = unit abbreviation
    zum = (str) file label

    Returns
    -------
    N/A
    """
    # Define default values.
    base_press = fault_controls.get('simple_pressure')
    base_temp = fault_controls.get('simple_temperature')
    base_depth = fault_controls.get('simple_depth')

    if not fault_controls['use_si']:
        # US units.
        base_press *= funits.pa_to_psi()
        base_temp = funits.celsius_to_fahrenheit(base_temp)
        base_depth = funits.meters_to_feet()
    else:
        # In metric, convert pressure to MPa.
        base_press *= funits.pa_to_mpa()

    #  *******************************************************
    # Print simple model data.
    print(ZIPX_IN + 'Near-Surface Depth         = '
          f'{base_depth:.1f} {uts[0]}', file=zum)
    print(ZIPX_IN + 'Near-Surface Pressure      = '
          f'{base_press:.4e} {uts[1]}', file=zum)
    print(ZIPX_IN + 'Near-Surface Temperature   = '
          f'{base_temp:.0f} {uts[7]}', file=zum)

    # return None


def write_profile_complex(fault_controls, uts, zum):
    """Write complex CO2 transition parameters to file named "zum".

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    uts = (list - str) = unit abbreviation
    zum = (str) file label

    Returns
    -------
    N/A
    """
    # Define default values.
    crux_press = fault_controls.get('crux_pressure')
    crux_temp = fault_controls.get('crux_temp')
    crux_depth = fault_controls.get('crux_depth')

    if not fault_controls['use_si']:
        # US units.
        crux_press *= funits.pa_to_psi()
        crux_temp = funits.celsius_to_fahrenheit(crux_temp)
        crux_depth = funits.meters_to_feet()
    else:
        # In metric, convert pressure to MPa.
        crux_press *= funits.pa_to_mpa()

    #  *********************************************************
    # Print complex model data.
    print(ZIPX_IN + 'Final Crux Depth            = '
          f'{crux_depth:.1f} {uts[0]}', file=zum)
    print(ZIPX_IN + 'Final Crux Pressure         = '
          f'{crux_press:.4e} {uts[1]}', file=zum)
    print(ZIPX_IN + 'Final Crux Temperature      = '
          f'{crux_temp:.0f} {uts[7]}', file=zum)

    # return None


def write_pressure_aperture_parameters(fault_controls, uts, zum):
    """Write parameters for pressure-aperture correction to file named "zum".

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    uts = (list - str) = unit abbreviation
    zum = (str) file label

    Returns
    -------
    N/A
    """
    # Write header, convert terms for pressure, if used.
    if fault_controls.get('pressure_approach'):
        print(BLK_IN + 'Pressure-Aperture Parameters Used', file=zum)

        # Define terms.
        frac_stiffness = fault_controls.get('pa_frac_stiffness')
        limit = fault_controls.get('pa_limit_stress')
        initial = fault_controls.get('pa_initial_stress')
        vary_aperture = fault_controls.get('pa_vary_aperture')
        max_aperture = fault_controls.get('pa_max_aperture')
        residual_aperture = fault_controls.get('pa_residual_aperture')
        initial_residual = fault_controls.get('pa_initial_residual')
        initial_max = fault_controls.get('pa_initial_max')

        if not fault_controls['use_si']:
            # US units.
            frac_stiffness *= funits.pa_to_psi()
            limit *= funits.pa_to_psi()
            initial *= funits.pa_to_psi()
            vary_aperture *= funits.mm_to_inch()
            max_aperture *= funits.mm_to_inch()
            residual_aperture *= funits.mm_to_inch()
            initial_residual *= funits.mm_to_inch()
            initial_max *= funits.mm_to_inch()
        else:
            frac_stiffness *= funits.pa_to_mpa()
            limit *= funits.pa_to_mpa()
            initial *= funits.pa_to_mpa()

    #  *****************************************************************
        # Write values if approach was used.
        print(SPEC_IN + 'Pressure-Aperture Approach was Used', file=zum)

        print(SPEC_IN + 'Computed Terms for Approach:', file=zum)
        print(ZIPX_IN + 'Beta Factor                = '
              f'{fault_controls.get("pa_beta"):.5f}', file=zum)
        print(ZIPX_IN + 'Gamma Factor               = '
              f'{fault_controls.get("pa_gamma_ref"):.5f}', file=zum)
        print(ZIPX_IN + 'Stiffness Factor           = '
              f'{frac_stiffness:.3f} {uts[1]}', file=zum)
        print(ZIPX_IN + 'Limit Pressure             = '
              f'{limit:.1f} {uts[1]}', file=zum)
        print(ZIPX_IN + 'Initial Pressure           = '
              f'{initial:.1f} {uts[1]}', file=zum)
        print(ZIPX_IN + 'Variable Aperture          = '
              f'{vary_aperture:.5f} {uts[8]}', file=zum)
        print(ZIPX_IN + 'Updated Maximum Aperture   = '
              f'{max_aperture:.5f} {uts[8]}', file=zum)
        print(ZIPX_IN + 'Updated Residual Aperture  = '
              f'{residual_aperture:.5f} {uts[8]}', file=zum)
        print(ZIPX_IN + 'Initial Residual Aperture  = '
              f'{initial_residual:.5f} {uts[8]}', file=zum)
        print(ZIPX_IN + 'Initial Maximum Aperture   = '
              f'{initial_max:.5f} {uts[8]}', file=zum)

    else:
        # Write that approach was NOT used.
        print(SPEC_IN + 'A Pressure-Aperture Approach was Not Used', file=zum)

    # return None


def write_fault_parameters(fault, zum):
    """Write fault plate parameters to file named "zum".

    Parameters
    ----------
    fault = (class) collection of fault section parameters
    zum = (str) file label

    Returns
    -------
    N/A
    """
    # Write plate values for last simulation.
    print(BLK_IN + 'Details for Last Simulation', file=zum)
    for indx, part in enumerate(fault):
        # Write plate values.
        print(BLK_IN + f'Properties of Fault Section #{indx}', file=zum)
        part.print_plate(zum)

    # return None


def write_co2_results(sim_flux_list, fault_controls, uts, zum):
    """Compute and print results of simulations to file "zum".

    Parameters
    ----------
    sim_flux_list = (list) CO2 results at the end of simulations
    fault_controls = (dict) dictionary of fault parameters
    uts = (str, array) units for output
    zum = (str) file label

    Returns
    -------
    N/A
    """
    # Convert data to Numpy array.
    flux_data = np.asarray(sim_flux_list)
    simulation_number = len(flux_data)

    # Create positive simulation data and get length.
    positive_data = flux_data[np.where(flux_data > 0.0)]
    fault_number = len(positive_data)

    # compute statistics of the array.
    if len(positive_data) > 1:
        ave_flux_data = np.mean(positive_data)
        std_dev = np.std(positive_data)
        min_flux_data = np.min(positive_data)
        max_flux_data = np.max(positive_data)
        # argument_min = np.argmin(positive_data) + 1  # Watch zeros!
        low_index = sim_flux_list.index(min_flux_data) + 1
        argument_max = np.argmax(flux_data) + 1
    else:
        # ERROR!
        ave_flux_data = 0.0
        std_dev = -999.0
        min_flux_data = 00.0
        max_flux_data = 0.00
        low_index = 0
        argument_max = 0

    if not fault_controls['use_si']:
        # US units.
        ave_flux_data *= funits.tonne_to_ton()
        std_dev *= funits.tonne_to_ton()
        min_flux_data *= funits.tonne_to_ton()
        max_flux_data *= funits.tonne_to_ton()

    #  **************************************************************
    # print separator.
    print('\n  ', end='', file=zum)
    print('*' * 60, file=zum)

    # print statistics.
    print(BLK_IN + 'Summary of Results', file=zum)
    print(ZIPX_IN + 'Total Number of Simulations           = '
                    f'{simulation_number}', file=zum)

    # Print statistics when data exists.
    if len(positive_data) > 1:
        print(ZIPX_IN + 'Simulations with Faults               = '
                        f'{fault_number}', file=zum)
        print(ZIPX_IN + 'Average CO2 Flux of Fault Simulations = '
                        f'{ave_flux_data:.3E} {uts[6]}', file=zum)
        print(ZIPX_IN + 'Standard Deviation of CO2 Flux        = '
                        f'{std_dev:.3E} {uts[6]}', file=zum)
        # --> Minimum index may be incorrect if faults before index.
        print(ZIPX_IN + 'Minimum CO2 Flux of Fault Simulations = '
              f'{min_flux_data:<.3E} {uts[6]}', file=zum)
        print(f'       - at Index =      {low_index}', file=zum)
        print(ZIPX_IN + 'Maximum CO2 Flux of Fault Simulations = '
              f'{max_flux_data:.3E} {uts[6]}', file=zum)
        print(f'       - at Index =      {argument_max}', file=zum)
    else:
        # Error in data for statistics.
        print(ZIPX_IN + 'INPUT ERROR - ALL DATA iS NEAR_ZERO!', file=zum)

    # return None


def write_last(fault_controls, zum):
    """Write current time value to file named "zum" and end.

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    zum = (str) file label

    Returns
    -------
    N/A
    """
    # Print separator.
    print('\n  ', end='', file=zum)
    print('*' * 60, file=zum)

    # Print computation time.
    print(f'\n  Computation Time = {fault_controls.get("elapsed")} '
          '(h.mm.ss.xxxx)', file=zum)

    # Print last line.
    print('\n  End', file=zum)

    # return None


def write_summary(fault_controls, sim_flux_list, fault):
    """Control point for the printout of summary data to file.

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    sim_flux_list = (list) CO2 flux at end of simulations
    fault = (class) a collection of plates

    Returns
    -------
    N/A

    Notes
    -----
    SUMMARY_FILE = Defined destination file name.
    """
    # Construct full path to results directory and summary file.
    #  -- Check for extension "txt".
    sub_path, destination = fileop.get_path_name(scfg.OUTPUT_DIR,
                                                 scfg.SUMMARY_FILE,
                                                 scfg.EXTENSION_TXT)

    # Define US or metric unit array of names.
    uts = define_units_for_output(fault_controls['use_si'])

    # Write  information to file on Fault_Flo run.
    try:
        with open(destination, 'w', encoding="utf8") as zum:
            # Input parameters.
            write_iam_parameters(fault_controls, zum)
            write_fault_controls(fault_controls, zum)
            write_description_parameters(fault_controls, uts, zum)
            write_field_parameters(fault_controls, uts, zum)
            write_aperture_parameters(fault_controls, uts, zum)
            write_deep_fluid_values(fault_controls, uts, zum)
            write_aqui_fluid_values(fault_controls, uts, zum)
            write_relative_perm(fault_controls, uts, zum)
            write_stress_parameters(fault_controls, uts, zum)
            write_profile_parameters(fault_controls, uts, zum)
            write_pressure_aperture_parameters(fault_controls, uts, zum)

            # Write fault plate info.
            if ECHO:
                write_fault_parameters(fault, zum)

            # Write results.
            write_co2_results(sim_flux_list, fault_controls, uts, zum)
            write_last(fault_controls, zum)

    except OSError as err:
        fileop.io_snag(err, sub_path, scfg.SUMMARY_FILE)

    # return None


#
# -----------------------------------------------------------------------------
# End of module
