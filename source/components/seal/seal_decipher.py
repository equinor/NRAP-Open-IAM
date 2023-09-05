#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions for reading seal YAML data and units.

Author: Ernest N. Lindner
Date: 08/19/2022

Module Name
    seal_decipher

Contents (17)
    esc(code)
    set_output_config(seal_controls)
    translate_other(docs, seal_controls)
    translate_description(docs, seal_controls)
    translate_grid_parameters(docs, seal_controls)
    translate_fluid_parameters(docs, seal_controls)
    translate_switches(docs, seal_controls):
    translate_thickness(docs, seal_controls)
    translate_permeability(docs, seal_controls)
    translate_two_phase(docs, seal_controls)
    ----
    translate_time_params(docs, seal_controls)
    translate_plots(docs, seal_controls)
    acquire_seal_parameters(yaml_file_selected, seal_controls)
    interpret_yaml_parameters(docs, seal_controls)
    recast_unit_to_default(units, parameter_valu, term)
    preface(alone, python_major, python_minor, seal_version)
    boot_code(alone, yaml_file, python_major, python_minor, seal_version)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import logging                      # For reporting errors
import time                         # For calculating runtime

import seal_file as sfile           # For file operations
import seal_config as scfg           # IO directory and file names
import seal_units as sun           # Unit translations
import seal_warranty as war         # Print warranty

# Constants
MINIMUM_PERM = 1.0e-06              # Minimum std. dev. of permeability
MISSING = [None, "", "default"]     # Check for Yaml file name
ECHO = False                        # Control to print for debugging

# Unit dictionaries for conversion
LIST_OF_ALL = [
    ' ', 'm', 'ft', 'Pa', 'psi', 'm^2', 'ft^2', 'oF', 'oC', 'pcf', 'kg/m^3',
    'cP', 'Pa*s', 'yrs', 'ppm', 'mD', 'pm^2', 'mol/kg', 'mol/lbs',
    'mm', 'in', 'MPa', 'per_ft^2', 'per_m^2', 'ton', 'tonne'
    ]

UNIT_DICT = {
    'ft': sun.feet_to_m(),
    'in': sun.inch_to_mm(),
    'ftˆ2': sun.ft2_to_m2(),
    'per_ft^2': sun.per_sqft_to_sqm(),
    'psi': sun.psi_to_pa(),
    'pcf': sun.pounds_convert_density(),
    'pm^2': sun.metersq_to_microd(),
    'mol/lbs': 1.00 / sun.pound_to_kilogram(),
    'cP': sun.centipoise_to_pa_s(),
    'MPa': sun.mpa_to_pa(),
    'ton': sun.ton_to_tonne()
    }

logging.basicConfig(format='\n    > %(levelname)s: %(message)s',
                    level=logging.WARNING)


def esc(code):
    """Print escape code.

    Parameters
    ----------
    code = (Int) print code to be preceded by esc key.

    Returns
    -------
    code = escape code.

    Notes
    -----
    1. Used to change color of text on screen.
    """
    return f'\033[{code}m'


def set_output_config(seal_controls):
    """Reset output directory for code config file.

    Parameters
    ----------
    seal_controls = (dict) seal control parameters

    Returns
    -------
    None
    """
    test = seal_controls['output_directory']

    if test and test.strip():
        scfg.OUTPUT_DIRECTORY = seal_controls['output_directory']

    # return None


def translate_other(docs, seal_controls):
    """Convert YAML other header parameters into seal controls.

    Parameters
    ----------
    docs = (str) YAML array input
    seal_controls = (dict) seal control parameters

    Returns
    -------
    seal_controls = (dict) seal control parameters(modified)
    """
    # Get other ModelParams.
    cubit = str(docs['ModelParams']['startTime']['unit'])
    worth = docs['ModelParams']['startTime']['valu']
    term = 'start time'
    seal_controls['start_time'] = recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['ModelParams']['endTime']['unit'])
    worth = docs['ModelParams']['endTime']['valu']
    term = 'end time'
    seal_controls['end_time'] = recast_unit_to_default(cubit, worth, term)

    seal_controls['time_points'] = int(docs['ModelParams']['timePoints'])
    seal_controls['realizations'] = int(docs['ModelParams']['realizations'])
    seal_controls['time_input'] = docs['ModelParams']['timeInput']
    seal_controls['output_directory'] = docs['OutputDirectory']

    cubit = str(docs['ModelParams']['totalInject']['unit'])
    worth = docs['ModelParams']['totalInject']['valu']
    term = 'total injection weight'
    seal_controls['total_inject'] = \
        recast_unit_to_default(cubit, worth, term)

    if seal_controls['output_directory'] in MISSING:
        seal_controls['output_directory'] = scfg.OUTPUT_DIRECTORY

    return seal_controls


def translate_description(docs, seal_controls):
    """Convert YAML grid parameters into seal controls.

    Parameters
    ----------
    docs = (str) YAML array input
    seal_controls = (dict) seal control parameters

    Returns
    -------
    seal_controls = (dict) seal control parameters(modified)
    """
    # Get Description parameters.
    seal_controls['num_cells'] = int(docs['Description']['numCells'])

    cubit = str(docs['Description']['area']['unit'])
    worth = docs['Description']['area']['valu']
    term = 'area'
    seal_controls['area'] = recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Description']['salinity']['unit'])
    worth = docs['Description']['salinity']['valu']
    term = 'salinity'
    seal_controls['salinity'] = recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Description']['aveTemperature']['unit'])
    worth = docs['Description']['aveTemperature']['valu']
    term = 'temperature'
    seal_controls['temperature'] = recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Description']['aveBaseDepth']['unit'])
    worth = docs['Description']['aveBaseDepth']['valu']
    term = 'average base depth'
    seal_controls['ave_base_depth'] = \
        recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Description']['aveBasePressure']['unit'])
    worth = docs['Description']['aveBasePressure']['valu']
    term = 'average base pressure'
    seal_controls['ave_base_pressure'] = recast_unit_to_default(cubit, worth,
                                                                term)

    cubit = str(docs['Description']['staticDepth']['unit'])
    worth = docs['Description']['staticDepth']['valu']
    term = 'static depth'
    seal_controls['static_depth'] = recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Description']['staticPressure']['unit'])
    worth = docs['Description']['staticPressure']['valu']
    term = 'static pressure'
    seal_controls['static_pressure'] = recast_unit_to_default(cubit, worth,
                                                              term)

    return seal_controls


def translate_grid_parameters(docs, seal_controls):
    """Convert YAML grid parameters into seal controls.

    Parameters
    ----------
    docs = (str) YAML array input
    seal_controls = (dict) seal control parameters

    Returns
    -------
    seal_controls = (dict) seal control parameters(modified)
    """
    # Get grid parameters.
    seal_controls['grid_approach'] = docs['Grid']['gridApproach']

    if seal_controls['grid_approach']:
        seal_controls['grid_rows'] = int(docs['Grid']['gridRows'])
        seal_controls['grid_cols'] = int(docs['Grid']['gridColumns'])

        cubit = str(docs['Grid']['cellHeight']['unit'])
        worth = docs['Grid']['cellHeight']['valu']
        term = 'cell height'
        seal_controls['cell_height'] = recast_unit_to_default(cubit, worth,
                                                              term)

        cubit = str(docs['Grid']['cellWidth']['unit'])
        worth = docs['Grid']['cellWidth']['valu']
        term = 'cell width'
        seal_controls['cell_width'] = \
            recast_unit_to_default(cubit, worth, term)

    return seal_controls


def translate_fluid_parameters(docs, seal_controls):
    """Convert YAML fluid parameters into seal controls.

    Parameters
    ----------
    docs = (str) YAML array input
    seal_controls = (dict) seal control parameters

    Returns
    -------
    seal_controls = (dict) seal control parameters(modified)
    """
    # Get conditions - fluid properties.
    cubit = str(docs['Conditions']['CO2Density']['unit'])
    worth = docs['Conditions']['CO2Density']['valu']
    term = 'CO2 density'
    seal_controls['co2_density'] = recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Conditions']['CO2Viscosity']['unit'])
    worth = docs['Conditions']['CO2Viscosity']['valu']
    term = 'CO2 viscosity'
    seal_controls['co2_viscosity'] = recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Conditions']['brineDensity']['unit'])
    worth = docs['Conditions']['brineDensity']['valu']
    term = 'brine density'
    seal_controls['brine_density'] = recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['Conditions']['brineViscosity']['unit'])
    worth = docs['Conditions']['brineViscosity']['valu']
    term = 'brine viscosity'
    seal_controls['brine_viscosity'] = recast_unit_to_default(cubit, worth,
                                                              term)

    cubit = str(docs['Conditions']['CO2Solubility']['unit'])
    worth = docs['Conditions']['CO2Solubility']['valu']
    term = 'CO2 solubility'
    seal_controls['co2_solubility'] = \
        recast_unit_to_default(cubit, worth, term)

    return seal_controls


def translate_switches(docs, seal_controls):
    """Convert YAML switches/control parameters into seal controls.

    Parameters
    ----------
    docs = (str) YAML array input
    seal_controls = (dict) seal control parameters

    Returns
    -------
    seal_controls = (dict) seal control parameters(modified)
    """
    # Get file input flags:
    seal_controls['depth_approach'] = docs['FileInput']['depthApproach']
    seal_controls['layout_approach'] = docs['FileInput']['layoutApproach']
    seal_controls['area_approach'] = docs['FileInput']['areaApproach']
    seal_controls['upper_approach'] = docs['FileInput']['upperBoundApproach']
    seal_controls['active_approach'] = docs['FileInput']['activeCellsApproach']
    seal_controls['entry_approach'] = docs['FileInput']['entryApproach']
    seal_controls['perm_input_approach'] = docs['FileInput']['permApproach']
    seal_controls['thickness_approach'] = \
        docs['FileInput']['thicknessApproach']

    # Controls.
    seal_controls['fracture_approach'] = docs['Controls']['fractureApproach']
    seal_controls['correlate_entry_approach'] = \
        docs['Controls']['correlateApproach']
    seal_controls['interpolate_approach'] = \
        docs['Controls']['interpolateApproach']
    seal_controls['initialize'] = docs['Controls']['initializeApproach']

    return seal_controls


def translate_thickness(docs, seal_controls):
    """Convert YAML thickness parameters into seal controls.

    Parameters
    ----------
    docs = (str) YAML array input
    seal_controls = (dict) seal control parameters

    Returns
    -------
    seal_controls = (dict) seal control parameters(modified)
    """
    # Set thickness parameters if thickness <not> from file.
    if not seal_controls['thickness_approach']:

        # Mean thickness.
        cubit = str(docs['Thickness']['aveThickness']['unit'])
        worth = docs['Thickness']['aveThickness']['valu']
        term = 'average thickness'
        seal_controls['thickness_ave'] = recast_unit_to_default(cubit, worth,
                                                                term)
        # Standard deviation thickness.
        cubit = str(docs['Thickness']['stdDevThickness']['unit'])
        worth = docs['Thickness']['stdDevThickness']['valu']
        term = 'standard deviation of thickness'
        seal_controls['thickness_std'] = recast_unit_to_default(cubit, worth,
                                                                term)
        # Minimum thickness.
        cubit = str(docs['Thickness']['minThickness']['unit'])
        worth = docs['Thickness']['minThickness']['valu']
        term = 'minimum thickness'
        seal_controls['thickness_min'] = recast_unit_to_default(cubit, worth,
                                                                term)
        # Maximum thickness.
        cubit = str(docs['Thickness']['maxThickness']['unit'])
        worth = docs['Thickness']['maxThickness']['valu']
        term = 'maximum thickness'
        seal_controls['thickness_max'] = recast_unit_to_default(cubit, worth,
                                                                term)

    return seal_controls


def translate_permeability(docs, seal_controls):
    """Convert YAML permeability parameters into seal controls.

    Parameters
    ----------
    docs = (str) YAML array input
    seal_controls = (dict) seal control parameters

    Returns
    -------
    seal_controls = (dict) seal control parameters(modified)
    """
    # Set permeability parameters if permeability <not> from file.
    if not seal_controls['perm_input_approach']:

        # Mean permeability.
        cubit = str(docs['Permeability']['avePermeability']['unit'])
        worth = docs['Permeability']['avePermeability']['valu']
        term = 'average permeability'
        seal_controls['perm_mean'] = recast_unit_to_default(cubit, worth, term)

        # Standard Deviation permeability.
        cubit = str(docs['Permeability']['stdDevPermeability']['unit'])
        worth = docs['Permeability']['stdDevPermeability']['valu']
        term = 'standard deviation permeability'
        seal_controls['perm_std'] = recast_unit_to_default(cubit, worth, term)

        # Minimum Deviation permeability.
        cubit = str(docs['Permeability']['minPermeability']['unit'])
        worth = docs['Permeability']['minPermeability']['valu']
        term = 'minimum permeability'
        seal_controls['perm_min'] = recast_unit_to_default(cubit, worth, term)

        # Maximum Deviation permeability.
        cubit = str(docs['Permeability']['maxPermeability']['unit'])
        worth = docs['Permeability']['maxPermeability']['valu']
        term = 'maximum permeability'
        seal_controls['perm_max'] = recast_unit_to_default(cubit, worth, term)

        # Heterogeneity.
        seal_controls['heterogeneity_approach'] = \
            docs['Permeability']['heterogeneityApproach']
        seal_controls['perm_heter_factor'] = \
            docs['Permeability']['heterFactor']

        # Internal control on stochastic permeability - False = no variability.
        if ((seal_controls['perm_std'] < MINIMUM_PERM) or
                (seal_controls['perm_input_approach'])):
            seal_controls['vary_perm_choice'] = False
        else:
            seal_controls['vary_perm_choice'] = True
    else:
        # For file input, heterogeneity is not considered.
        seal_controls['heterogeneity_approach'] = False
        seal_controls['vary_perm_choice'] = False

    return seal_controls


def translate_two_phase(docs, seal_controls):
    """Convert YAML two_phase parameters into seal controls.

    Parameters
    ----------
    docs = (stream) YAML array input
    seal_controls = (dict) seal control parameters

    Returns
    -------
    seal_controls = (dict) seal control parameters (modified)
    """
    # Set relative flow limits.
    seal_controls['relative_model'] = \
        docs['RelativeFlowLimits']['relativeModel']
    seal_controls['resid_brine'] = \
        docs['RelativeFlowLimits']['brineResSaturation']
    seal_controls['resid_co2'] = \
        docs['RelativeFlowLimits']['CO2ResSaturation']
    seal_controls['perm_ratio'] = \
        docs['RelativeFlowLimits']['permRatio']

    # Set capillary pressure.
    cubit = str(docs['CapillaryPressure']['entryPressure']['unit'])
    worth = docs['CapillaryPressure']['entryPressure']['valu']
    term = 'capillary pressure'
    seal_controls['entry_pressure'] = \
        recast_unit_to_default(cubit, worth, term)

    # Set two_phase parameters depending on model.
    if 'BC' in seal_controls['relative_model']:
        # Brooks Corey model,
        #   -> Avoid using lambda in Python source -> use "zeta".
        seal_controls['zeta'] = docs['BrooksCoreyModel']['lambda']

    else:
        # Set LET model parameters.
        seal_controls['l_wetting'] = docs['LETModel']['wetting1']
        seal_controls['e_wetting'] = docs['LETModel']['wetting2']
        seal_controls['t_wetting'] = docs['LETModel']['wetting3']
        seal_controls['l_nonwet'] = docs['LETModel']['nonwet1']
        seal_controls['e_nonwet'] = docs['LETModel']['nonwet2']
        seal_controls['t_nonwet'] = docs['LETModel']['nonwet3']

        # Set LET Capillary pressure.
        seal_controls['l_capillary'] = docs['LETCapillaryModel']['capillary1']
        seal_controls['e_capillary'] = docs['LETCapillaryModel']['capillary2']
        seal_controls['t_capillary'] = docs['LETCapillaryModel']['capillary3']

        cubit = str(docs['LETCapillaryModel']['maxCapillary']['unit'])
        worth = docs['LETCapillaryModel']['maxCapillary']['valu']
        term = 'maximum capillary pressure'
        seal_controls['max_capillary'] = recast_unit_to_default(cubit, worth,
                                                                term)

    return seal_controls


def translate_time_params(docs, seal_controls):
    """Convert YAML Time Model parameters into seal controls.

    Parameters
    ----------
    docs = (stream) YAML array input
    seal_controls = (dict) seal control parameters

    Returns
    -------
    seal_controls = (dict) seal control parameters (modified)
    """
    # Determine model type & set initial regardless of model.
    seal_controls['model'] = docs['TimeModel']['influenceModel']
    seal_controls['influence'] = docs['TimeModel']['influence']

    # Set model parameters if model =1 or =2.
    if seal_controls['model'] > 0:
        seal_controls['rate_effect'] = docs['TimeModel']['rateEffect']
        seal_controls['total_effect'] = docs['TimeModel']['totalEffect']

        # For model=2 only.
        if seal_controls['model'] > 1:
            seal_controls['reactivity'] = docs['TimeModel']['reactivity']
            seal_controls['clay_content'] = docs['TimeModel']['clayContent']
            seal_controls['carbonate_content'] = \
                docs['TimeModel']['carbonateContent']
            seal_controls['clay_type'] = docs['TimeModel']['clayType']

    return seal_controls


def translate_plots(docs, seal_controls):
    """Convert YAML plot parameters into seal controls.

    Parameters
    ----------
    docs = (str) YAML array input
    seal_controls = (dict) seal control parameters

    Returns
    -------
    seal_controls = (dict) seal control parameters(modified)
    """
    # Set  plot control parameters.
    seal_controls['plot_permeability'] = docs['SealPlots']['permeabilityPlot']
    seal_controls['plot_time_history'] = docs['SealPlots']['timeSeriesPlot']
    seal_controls['plot_co2_contour'] = docs['SealPlots']['CO2ContourPlot']

    # Set max value on x-axis of plots.
    cubit = str(docs['SealPlots']['maxDrawTime']['unit'])
    worth = docs['SealPlots']['maxDrawTime']['valu']
    term = 'maximum plot time'
    seal_controls['max_draw_time'] = \
        recast_unit_to_default(cubit, worth, term)

    return seal_controls


def acquire_seal_parameters(yaml_file_selected, seal_controls):
    """Obtain seal parameters from a formatted YAML file.

    Parameters
    ----------
    yaml_file_selected = (str) file name of YAML file in same directory
    seal_controls = (dict) seal control parameters (empty)

    Returns
    -------
    seal_controls = (dict) seal control parameters - with values
    """
    # Get YAML data from file.
    docs = sfile.acquire_yaml_file(yaml_file_selected)

    # Translate parameters into seal_controls dictionary.
    try:
        interpret_yaml_parameters(docs, seal_controls)
    except KeyError as errd:
        cause = errd.args[0]
        msg = f'Parameter <{cause}> is Missing from YAML file!'
        sfile.opx_problem(msg + "\n      " +
                          "- KeyError in Defining Seal_Controls Dictionary!")

    return seal_controls


def interpret_yaml_parameters(docs, seal_controls):
    """Translate YAML parameters into seal controls.

    Parameters
    ----------
    docs = (stream) YAML array input
    seal_controls = (dict) seal control parameters

    Returns
    -------
    seal_controls = (dict) seal control parameters - updated

    Notes
    -----
    1. Parameters that not needed are not set and therefore are not checked.
    2. Breakout in several routines is to avoid report of errors by Spyder.
    """
    # Set model parameters - run title:
    seal_controls['title'] = docs['ModelParams']['runTitle']
    if seal_controls['title'] == "":
        seal_controls['title'] = " Seal Flux Code Run"

    # -- Set output units selection.
    choice_units = str(docs['ModelParams']['useSI'])
    select_units = choice_units.upper()
    if "ENG" in select_units:
        seal_controls['use_si_units'] = False
    else:
        seal_controls['use_si_units'] = True

    # Set other header parameters.
    seal_controls = translate_other(docs, seal_controls)

    # Set condition description parameters
    seal_controls = translate_description(docs, seal_controls)

    # Set grid parameters.
    seal_controls = translate_grid_parameters(docs, seal_controls)

    # Set fluid parameters.
    seal_controls = translate_fluid_parameters(docs, seal_controls)

    # Set controls and switches.
    seal_controls = translate_switches(docs, seal_controls)

    # Set thickness parameters if thickness <not> from file.
    seal_controls = translate_thickness(docs, seal_controls)

    # Set permeability variability parameters.
    seal_controls = translate_permeability(docs, seal_controls)

    # Set two-phase model controls, depending on model selection.
    seal_controls = translate_two_phase(docs, seal_controls)

    # Set Time Model controls, if model > 0.
    seal_controls = translate_time_params(docs, seal_controls)

    # Sey seal plot controls.
    seal_controls = translate_plots(docs, seal_controls)

    return seal_controls


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
        sfile.opx_problem(reason, err='')

    # Check for temperature first.
    if cubit == 'oF':
        real_valu = sun.fahrenheit_to_celsius(real_valu)
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
        sfile.opx_problem(txt_blank, err='')

    return real_valu


def preface(alone, python_major, python_minor, seal_version):
    """Write header and warranty.

    Parameters
    ----------
    alone = (bool) status of code operation; True = if stand-alone
    python_major = (int) major number of python
    python_minor = (int) minor number of python
    seal_version = (str) seal_flux version number

    Returns
    -------
    None
    """
    if alone:
        # Show header statements on console - name in red.
        print('\n  ' + esc('31;1') + 'SEAL_FLUX' + esc(0))
        py_version = str(python_major) + "." + str(python_minor)
        print("\n  *** SEAL FLUX MODEL - Version: " + seal_version +
              " ***" + "           Python " + py_version)

        # Echo warranty and code run.
        war.show_warranty()
        sfile.echo_status(alone, "READING SEAL CONTROL FILE.")

    # return


def boot_code(alone, yaml_file, python_major, python_minor, seal_version):
    """Open control file and check data.

    Parameters
    ----------
    alone = (bool) status of code operation; True = if stand-alone
    yaml_file = (str) name of input control file (yaml file)
    python_major = (int) major number of python
    python_minor = (int) minor number of python
    seal_version = (str) seal_flux version number

    Returns
    -------
    seal_controls = (dict) seal control parameters
    """
    # Start clock for current simulations.
    timer = time.monotonic()

    # Print warranty.
    preface(alone, python_major, python_minor, seal_version)

    # Check Python version being used.
    sfile.check_python(python_major, python_minor)

    # Save code start time in seal_controls for use later.
    seal_controls = {}
    seal_controls['code_start'] = timer
    seal_controls['version'] = seal_version

    # Get input control parameters from text file.
    seal_controls = acquire_seal_parameters(yaml_file, seal_controls)

    return seal_controls


#
# -----------------------------------------------------------------------------
# - End of module
