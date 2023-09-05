#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions for interpreting YAML file.

Author: Ernest N. Lindner
Date: 08/17/2022

Module Name
    flt_decipher

Contents (13)
    recast_unit_to_default(cubit, parameter_valu, term)
    translate_model_parameters(docs, fault_controls)
    translate_switches(docs, fault_controls)
    translate_fault_core(docs, fault_controls):
    translate_field_parameters(docs, fault_controls)
    translate_inject_parameters(docs, fault_controls)
    translate_aquifer_parameters(docs, fault_controls)
    translate_aperture_parameters(docs, fault_controls)
    translate_phase_parameters(docs, fault_controls)
    translate_stress_parameters(docs, fault_controls)
    ---
    translate_plots(docs, fault_controls)
    interpret_fault_parameters(docs, fault_controls)
    acquire_control_parameters(yaml_file_selected, fault_controls)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import flt_file as fileop  # For file operations
import flt_units as funits  # For units

# Misc constants
SMALL_DEV = 1.0e-6  # Minimum std deviation for aperture

# Unit dictionaries for conversion
LIST_OF_ALL = [
    ' ', 'm', 'ft', 'Pa', 'psi', 'm^2', 'ft^2', 'oF', 'oC', 'pcf', 'kg/m^3',
    'cP', 'Pa*s', 'yrs', 'ppm', 'mD', 'pm^2', 'mol/kg', 'mol/lbs',
    'mm', 'in', 'MPa', 'per_ftˆ2', 'per_mˆ2'
]

UNIT_DICT = {
    'ft': funits.feet_to_m(),
    'in': funits.inch_to_mm(),
    'ftˆ2': funits.ft2_to_m2(),
    'per_ft^2': funits.per_sqft_to_sqm(),
    'psi': funits.psi_to_pa(),
    'pcf': funits.pounds_convert_density(),
    'pm^2': funits.metersq_to_microd(),
    'mol/lbs': 1.00 / funits.pound_to_kilogram(),
    'cP': funits.centipoise_to_pa_s(),
    'MPa': funits.mpa_to_pa()
}


def recast_unit_to_default(cubit, parameter_valu, term):
    """Convert units to default units used in code.

    Parameters
    ----------
    cubit = (str) abbreviation for unit
    parameter_valu = (float) current parameter value
    term = (str) name of variable

    Returns
    -------
    new_valu = (float) converted parameter value in default units
    """
    # Ensure float; specify default for default values.
    try:
        real_valu = float(parameter_valu)
    except ValueError:
        reason = f"{parameter_valu} is not a valid number for input for {term}"
        real_valu = 0
        fileop.opx_problem(reason, err='')

    # Check for temperature first.
    if cubit == 'oF':
        real_valu = funits.fahrenheit_to_celsius(real_valu)

    else:
        # Check other values for conversion factor
        for key, val in UNIT_DICT.items():
            # Check each key for conversion match.
            if cubit == key:
                real_valu *= val
                break

    # Check if unit definition is a default.
    if cubit not in LIST_OF_ALL:
        txt_blank = f'Unit Abbreviation <{cubit}> used for '
        txt_blank += f'Parameter <{term}> in YAML File is Not Accepted!'
        fileop.opx_problem(txt_blank, err='')

    return real_valu


def translate_model_parameters(docs, fault_controls):
    """Convert model parameters into fault controls.

    Parameters
    ----------
    docs = (str) YAML array input
    fault_controls = (dict) fault control parameters

    Returns
    -------
    fault_controls = (dict) fault control parameters(modified)
    """
    # Set model parameters - run title.
    fault_controls['title'] = docs['ModelParams']['runTitle']
    if fault_controls['title'] == "":
        fault_controls['title'] = " Fault Flo Simulation Run"

    # Other ModelParams.
    fault_controls['start_time'] = docs['ModelParams']['startTime']
    fault_controls['end_time'] = docs['ModelParams']['endTime']
    fault_controls['inject_end'] = docs['ModelParams']['injectEndTime']
    fault_controls['time_points'] = int(docs['ModelParams']['timePoints'])
    fault_controls['realizations'] = int(docs['ModelParams']['realizations'])
    fault_controls['time_input'] = docs['ModelParams']['timeInput']

    # --> Set output units selection.
    choice_units = str(docs['ModelParams']['useSI'])
    select_units = choice_units.upper()
    fault_controls['units'] = select_units
    if "ENG" in select_units:
        fault_controls['use_si'] = False
    else:
        fault_controls['use_si'] = True

    return fault_controls


def translate_switches(docs, fault_controls):
    """Convert switches/control parameters into fault controls.

    Parameters
    ----------
    docs = (str) YAML array input
    fault_controls = (dict) fault control parameters

    Returns
    -------
    fault_controls = (dict) fault control parameters(modified)
    """
    # Define controls.
    fault_controls['aperture_approach'] = \
        bool(docs['Controls']['apertureApproach'])
    fault_controls['near_surface_approach'] = \
        bool(docs['Controls']['considerNearApproach'])
    fault_controls['profile_type'] = \
        int(docs['Controls']['profileType'])
    fault_controls['strike_approach'] = \
        bool(docs['Controls']['strikeApproach'])
    fault_controls['dip_approach'] = \
        bool(docs['Controls']['dipApproach'])
    fault_controls['pressure_approach'] = \
        bool(docs['Controls']['pressureApproach'])
    fault_controls['interpolate_approach'] = \
        docs['Controls']['interpolateApproach']

    return fault_controls


def translate_fault_core(docs, fault_controls):
    """Convert fault core parameters into values.

    Parameters
    ----------
    docs = (str) YAML array input
    fault_controls = (dict) fault control parameters

    Returns
    -------
    fault_controls = (dict) fault control parameters(modified)
    """
    # Define fault core variables.
    fault_controls['fault_probability'] = docs['FaultCore']['faultProbability']

    fault_controls['strike_mean'] = docs['FaultCore']['aveStrike']
    if fault_controls['strike_approach']:
        fault_controls['strike_sig'] = docs['FaultCore']['spreadStrike']

    fault_controls['dip_mean'] = docs['FaultCore']['aveDip']
    if fault_controls['dip_approach']:
        fault_controls['dip_std'] = docs['FaultCore']['stdDip']

    cubit = str(docs['FaultCore']['length']['unit'])
    worth = docs['FaultCore']['length']['valu']
    term = 'fault length'
    fault_controls['length'] = recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['FaultCore']['xStart']['unit'])
    worth = docs['FaultCore']['xStart']['valu']
    term = 'x start coordinate'
    fault_controls['x_start'] = recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['FaultCore']['yStart']['unit'])
    worth = docs['FaultCore']['yStart']['valu']
    term = 'y start coordinate'
    fault_controls['y_start'] = recast_unit_to_default(cubit, worth, term)

    fault_controls['n_plates'] = int(docs['FaultCore']['nSegments'])

    return fault_controls


def translate_field_parameters(docs, fault_controls):
    """Convert field parameters into fault controls.

    Parameters
    ----------
    docs = (str) YAML array input
    fault_controls = (dict) fault control parameters

    Returns
    -------
    fault_controls = (dict) fault control parameters(modified)
    """
    # Define aquifer properties.
    cubit = str(docs['Field']['aquiferDepth']['unit'])
    worth = docs['Field']['aquiferDepth']['valu']
    term = 'aquifer depth'
    fault_controls['aquifer_depth'] = recast_unit_to_default(cubit, worth,
                                                             term)

    cubit = str(docs['Field']['aquiferTemperature']['unit'])
    worth = docs['Field']['aquiferTemperature']['valu']
    term = 'aquifer temperature'
    fault_controls['aquifer_temperature'] = recast_unit_to_default(cubit,
                                                                   worth,
                                                                   term)

    cubit = str(docs['Field']['aquiferPressure']['unit'])
    worth = docs['Field']['aquiferPressure']['valu']
    term = 'aquifer pressure'
    fault_controls['aquifer_pressure'] = recast_unit_to_default(cubit, worth,
                                                                term)

    # Define injection properties.
    cubit = str(docs['Field']['injectDepth']['unit'])
    worth = docs['Field']['injectDepth']['valu']
    term = 'injection depth'
    fault_controls['inject_depth'] = recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Field']['injectTemperature']['unit'])
    worth = docs['Field']['injectTemperature']['valu']
    term = 'inject temperature'
    fault_controls['inject_temperature'] = recast_unit_to_default(cubit, worth,
                                                                  term)

    cubit = str(docs['Field']['fieldPressure']['unit'])
    worth = docs['Field']['fieldPressure']['valu']
    term = 'field pressure'
    fault_controls['field_pressure'] = recast_unit_to_default(cubit, worth,
                                                              term)

    cubit = str(docs['Field']['finalPressure']['unit'])
    worth = docs['Field']['finalPressure']['valu']
    term = 'final pressure'
    fault_controls['final_pressure'] = recast_unit_to_default(cubit, worth,
                                                              term)

    cubit = str(docs['Field']['injectPressure']['unit'])
    worth = docs['Field']['injectPressure']['valu']
    term = 'inject pressure'
    fault_controls['inject_pressure'] = recast_unit_to_default(cubit, worth,
                                                               term)

    cubit = str(docs['Field']['injectX']['unit'])
    worth = docs['Field']['injectX']['valu']
    term = 'inject x'
    fault_controls['inject_x'] = recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Field']['injectY']['unit'])
    worth = docs['Field']['injectY']['valu']
    term = 'inject x'
    fault_controls['inject_y'] = recast_unit_to_default(cubit, worth, term)

    return fault_controls


def translate_aperture_parameters(docs, fault_controls):
    """Convert aperture and correction parameters into fracture controls.

    Parameters
    ----------
    docs = (stream) YAML array input
    fault_controls = (dict) control parameter dictionary

    Returns
    -------
    frac_controls = (dict) (revised)
    """
    # Set VaryAperture parameters:
    cubit = str(docs['Aperture']['aveAperture']['unit'])
    worth = docs['Aperture']['aveAperture']['valu']
    term = 'average aperture'
    fault_controls['aperture_mean'] = \
        recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Aperture']['stdDevAperture']['unit'])
    worth = docs['Aperture']['stdDevAperture']['valu']
    term = 'aperture standard deviation'
    fault_controls['aperture_std'] = \
        recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Aperture']['minAperture']['unit'])
    worth = docs['Aperture']['minAperture']['valu']
    term = 'aperture minimum'
    fault_controls['aperture_min'] = \
        recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Aperture']['maxAperture']['unit'])
    worth = docs['Aperture']['maxAperture']['valu']
    term = 'aperture maximum'
    fault_controls['aperture_max'] = \
        recast_unit_to_default(cubit, worth, term)

    # Internal control on stochastic aperture - False = no variability.
    if ((fault_controls['aperture_std'] < SMALL_DEV) or
            (fault_controls['aperture_approach'])):
        fault_controls['vary_aperture'] = False
    else:
        fault_controls['vary_aperture'] = True

    # Corrections.
    fault_controls['sgr'] = docs['Aperture']['SGR']
    fault_controls['state_correction'] = docs['Aperture']['stateVariable']

    return fault_controls


def translate_aquifer_parameters(docs, fault_controls):
    """Convert aquifer fluid parameters into fault controls.

    Parameters
    ----------
    docs = (str) YAML array input
    fault_controls = (dict) fault control parameters

    Returns
    -------
    fault_controls = (dict) fault control parameters(modified)
    """
    # Define aquifer fluid properties:
    cubit = str(docs['AquiferConditions']['CO2Density']['unit'])
    worth = docs['AquiferConditions']['CO2Density']['valu']
    term = 'aquifer CO2 density'
    fault_controls['aqui_co2_density'] = \
        recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['AquiferConditions']['CO2Viscosity']['unit'])
    worth = docs['AquiferConditions']['CO2Viscosity']['valu']
    term = 'aquifer CO2 viscosity'
    fault_controls['aqui_co2_viscosity'] = \
        recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['AquiferConditions']['brineDensity']['unit'])
    worth = docs['AquiferConditions']['brineDensity']['valu']
    term = 'aquifer brine density'
    fault_controls['aqui_brine_density'] = \
        recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['AquiferConditions']['brineViscosity']['unit'])
    worth = docs['AquiferConditions']['brineViscosity']['valu']
    term = 'aquifer brine viscosity'
    fault_controls['aqui_brine_viscosity'] = \
        recast_unit_to_default(cubit, worth, term)

    return fault_controls


def translate_inject_parameters(docs, fault_controls):
    """Convert injection horizon fluid parameters into fault controls.

    Parameters
    ----------
    docs = (str) YAML array input
    fault_controls = (dict) fault control parameters

    Returns
    -------
    fault_controls = (dict) fault control parameters(modified)
    """
    # Set InjectConditions - fluid properties.
    fault_controls['salinity'] = docs['InjectConditions']['salinity']

    cubit = str(docs['InjectConditions']['CO2Density']['unit'])
    worth = docs['InjectConditions']['CO2Density']['valu']
    term = 'CO2 density'
    fault_controls['co2_density'] = \
        recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['InjectConditions']['CO2Viscosity']['unit'])
    worth = docs['InjectConditions']['CO2Viscosity']['valu']
    term = 'CO2 viscosity'
    fault_controls['co2_viscosity'] = \
        recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['InjectConditions']['brineDensity']['unit'])
    worth = docs['InjectConditions']['brineDensity']['valu']
    term = 'brine density'
    fault_controls['brine_density'] = \
        recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['InjectConditions']['brineViscosity']['unit'])
    worth = docs['InjectConditions']['brineViscosity']['valu']
    term = 'brine viscosity'
    fault_controls['brine_viscosity'] = \
        recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['InjectConditions']['CO2Solubility']['unit'])
    worth = docs['InjectConditions']['CO2Solubility']['valu']
    term = 'CO2 solubility'
    fault_controls['co2_solubility'] = \
        recast_unit_to_default(cubit, worth, term)

    return fault_controls


def translate_phase_parameters(docs, fault_controls):
    """Convert two_phase parameters into fault controls.

    Parameters
    ----------
    docs = (stream) YAML array input
    fault_controls = (dict) fault control parameters

    Returns
    -------
    fault_controls = (dict) fault control parameters (modified)
    """
    # Set relative flow limits.
    fault_controls['relative_model'] = \
        docs['RelativeFlowLimits']['relativeModel']
    fault_controls['resid_brine'] = \
        docs['RelativeFlowLimits']['brineResSaturation']
    fault_controls['resid_co2'] = \
        docs['RelativeFlowLimits']['CO2ResSaturation']
    fault_controls['perm_ratio'] = \
        docs['RelativeFlowLimits']['permRatio']

    # Set capillary pressure.
    cubit = str(docs['RelativeFlowLimits']['entryPressure']['unit'])
    worth = docs['RelativeFlowLimits']['entryPressure']['valu']
    term = 'capillary pressure'
    fault_controls['entry_pressure'] = \
        recast_unit_to_default(cubit, worth, term)

    # Set two_phase parameters depending on model.
    if 'BC' in fault_controls['relative_model']:
        # Brooks Corey model,
        #   -> Avoid using lambda in Python source -> use "zeta".
        fault_controls['zeta'] = docs['BrooksCoreyModel']['lambda']
    else:
        # Set LET model parameters.
        fault_controls['l_wetting'] = docs['LETModel']['wetting1']
        fault_controls['e_wetting'] = docs['LETModel']['wetting2']
        fault_controls['t_wetting'] = docs['LETModel']['wetting3']
        fault_controls['l_nonwet'] = docs['LETModel']['nonwet1']
        fault_controls['e_nonwet'] = docs['LETModel']['nonwet2']
        fault_controls['t_nonwet'] = docs['LETModel']['nonwet3']

        # Set LET Capillary pressure.
        fault_controls['l_capillary'] = docs['LETCapillaryModel'][
            'capillary1']
        fault_controls['e_capillary'] = docs['LETCapillaryModel'][
            'capillary2']
        fault_controls['t_capillary'] = docs['LETCapillaryModel'][
            'capillary3']

        cubit = str(docs['LETCapillaryModel']['maxCapillary']['unit'])
        worth = docs['LETCapillaryModel']['maxCapillary']['valu']
        term = 'maximum capillary pressure'
        fault_controls['max_capillary'] = \
            recast_unit_to_default(cubit, worth, term)

    return fault_controls


def translate_stress_parameters(docs, fault_controls):
    """Convert plot parameters into fault controls.

    Parameters
    ----------
    docs = (str) YAML array input
    fault_controls = (dict) fault control parameters

    Returns
    -------
    fault_controls = (dict) fault control parameters(modified)
    """
    # Define 3 stress components:
    cubit = str(docs['Stress']['maxHorizontal']['unit'])
    worth = docs['Stress']['maxHorizontal']['valu']
    term = 'maximum horizontal stress'
    fault_controls['max_horizontal'] = \
        recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Stress']['minHorizontal']['unit'])
    worth = docs['Stress']['minHorizontal']['valu']
    term = 'minimum horizontal stress'
    fault_controls['min_horizontal'] = \
        recast_unit_to_default(cubit, worth, term)

    fault_controls['max_trend'] = docs['Stress']['maxTrend']

    return fault_controls


def translate_plots(docs, fault_controls):
    """Convert YAML plot parameters into fault controls.

    Parameters
    ----------
    docs = (str) YAML array input
    fault_controls = (dict) fault control parameters

    Returns
    -------
    fault_controls = (dict) fault control parameters(modified)
    """
    # Set plot control parameters.
    fault_controls['plot_time_history'] = docs['FlowPlot']['timeSeriesPlot']
    fault_controls['skip_output'] = docs['FlowPlot']['skipOutput']

    return fault_controls


def interpret_fault_parameters(docs, fault_controls):
    """Translate YAML parameters into fault controls.

    Parameters
    ----------
    docs = (stream) YAML array input
    fault_controls = (dict) fault control parameters

    Returns
    -------
    fault_controls = (dict) fault control parameters - updated

    Notes
    -----
    1. Parameters that not needed are not set and therefore are not checked.
    2. Breakout in several routines is to avoid report of errors by Spyder.
    """
    # Set other header parameters.
    fault_controls = translate_model_parameters(docs, fault_controls)

    # Set control parameters
    fault_controls = translate_switches(docs, fault_controls)

    # Set fault core parameters.
    fault_controls = translate_fault_core(docs, fault_controls)

    # Set aperture parameters.
    fault_controls = translate_aperture_parameters(docs, fault_controls)

    # Set field parameters.
    fault_controls = translate_field_parameters(docs, fault_controls)

    # Set injection parameters.
    fault_controls = translate_inject_parameters(docs, fault_controls)

    # Set aquifer parameters.
    fault_controls = translate_aquifer_parameters(docs, fault_controls)

    # Set two-phase parameters.
    fault_controls = translate_phase_parameters(docs, fault_controls)

    # Set In Situ stress controls, if model > 0.
    fault_controls = translate_stress_parameters(docs, fault_controls)

    # Set fault plot controls.
    fault_controls = translate_plots(docs, fault_controls)

    return fault_controls


def acquire_control_parameters(yaml_file_selected, fault_controls):
    """Obtain the fault parameters from a YAML file.

    Parameters
    ----------
    yaml_file_selected = (str) file name of YAML file in same directory
    fault_controls = (dict) fault controls (empty)

    Returns
    -------
    fault_controls = (dict) fault parameters - with values
    """
    # Get YAML data from file.
    docs = fileop.acquire_yaml_file(yaml_file_selected)

    # Load parameters into fault_controls dictionary.
    try:
        interpret_fault_parameters(docs, fault_controls)
    except KeyError as errd:
        cause = errd.args[0]
        msg = f'Parameter <{cause}> is Missing from YAML file!'
        fileop.opx_problem(msg + "\n      " +
                           "- KeyError in Defining fault_controls Dictionary!")

    return fault_controls

#
# -----------------------------------------------------------------------------
# - End of module
