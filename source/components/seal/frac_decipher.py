#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions for getting frac YAML data.

Author: Ernest N. Lindner
Date: 08/19/2022

Module Name
    frac_decipher

Contents (4)
    acquire_frac_parameters(yaml_file_selected, frac_controls)
    translate_random_parameters(docs, frac_controls)
    translate_other_random_parameters(docs, frac_controls)
    interpret_frac_parameters(docs, frac_controls)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import seal_file as sfile           # For file operations
import seal_config as scfg          # IO directory and file names
import seal_decipher as yam         # For routine to check units

# Constants
MISSING = [None, "default"]         # Check for Yaml file


def acquire_frac_parameters(yaml_file_selected, frac_controls):
    """Obtain seal parameters from a formatted YAML file.

    Parameters
    ----------
    yaml_filename = (str) file name of YAML file in same directory
    frac_controls = (dict) dictionary of fracture controls (empty)

    Returns
    -------
    frac_controls = (dict) dictionary of fracture parameters - with values
    """
    # Get YAML data from file.
    docs = sfile.acquire_yaml_file(yaml_file_selected)

    # load parameters into frac_controls dictionary.
    frac_controls = interpret_frac_parameters(docs, frac_controls)

    return frac_controls


def translate_random_parameters(docs, frac_controls):
    """Convert random fracture YAML parameters into fracture controls.

    Parameters
    ----------
    docs = (stream) YAML array input
    frac_controls = (dict) control parameter dictionary

    Returns
    -------
    frac_controls = (dict) (revised)
    """
    # Density:
    cubit = str(docs['RandomFracs']['Density']['aveDensity']['unit'])
    worth = docs['RandomFracs']['Density']['aveDensity']['valu']
    term = 'average density'
    frac_controls['density_ave'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['RandomFracs']['Density']['minDensity']['unit'])
    worth = docs['RandomFracs']['Density']['minDensity']['valu']
    term = 'minimum density'
    frac_controls['density_min'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['RandomFracs']['Density']['maxDensity']['unit'])
    worth = docs['RandomFracs']['Density']['maxDensity']['valu']
    term = 'maximum density'
    frac_controls['density_max'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    # Orientation:
    frac_controls['orient_mu'] = \
        docs['RandomFracs']['Orientation']['aveTrend']
    frac_controls['orient_sigma'] = \
        docs['RandomFracs']['Orientation']['spreadTrend']

    # Setup for a function selection and input.
    # - Ensure uppercase text.
    frac_controls['length_approach'] = \
        docs['RandomFracs']['Length']['function']
    condition = frac_controls['length_approach'].upper()

    # Length:
    if "LOG" in condition:
        # - Lognormal Parameters.
        frac_controls['length_approach'] = 'LOGNORM'

        cubit = str(docs['RandomFracs']['Length']['aveLength']['unit'])
        worth = docs['RandomFracs']['Length']['aveLength']['valu']
        term = 'average fracture length'
        frac_controls['length_ave'] = \
            yam.recast_unit_to_default(cubit, worth, term)

        cubit = str(docs['RandomFracs']['Length']['stdDevLength']['unit'])
        worth = docs['RandomFracs']['Length']['stdDevLength']['valu']
        term = 'fracture length deviation'
        frac_controls['length_dev'] = \
            yam.recast_unit_to_default(cubit, worth, term)

    else:
        # - Power Parameters.
        frac_controls['length_approach'] = 'POWER'
        frac_controls['length_eta'] = docs['RandomFracs']['Length'][
            'expLength']

    # - For either lognormal or power Parameters.
    cubit = str(docs['RandomFracs']['Length']['minLength']['unit'])
    worth = docs['RandomFracs']['Length']['minLength']['valu']
    term = 'fracture minimum length'
    frac_controls['length_min'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['RandomFracs']['Length']['maxLength']['unit'])
    worth = docs['RandomFracs']['Length']['maxLength']['valu']
    term = 'fracture maximum length'
    frac_controls['length_max'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    return frac_controls


def translate_other_random_parameters(docs, frac_controls):
    """Convert more random fracture YAML parameters into fracture controls.

    Parameters
    ----------
    docs = (stream) YAML array input
    frac_controls = (dict) control parameter dictionary

    Returns
    -------
    frac_controls = (dict) (revised)
    """
    # VaryAperture:
    cubit = str(docs['RandomFracs']['VaryAperture']['aveAperture']['unit'])
    worth = docs['RandomFracs']['VaryAperture']['aveAperture']['valu']
    term = 'average aperture'
    frac_controls['aperture_ave'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['RandomFracs']['VaryAperture']['stdDevAperture']['unit'])
    worth = docs['RandomFracs']['VaryAperture']['stdDevAperture']['valu']
    term = 'aperture standard deviation'
    frac_controls['aperture_dev'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['RandomFracs']['VaryAperture']['minAperture']['unit'])
    worth = docs['RandomFracs']['VaryAperture']['minAperture']['valu']
    term = 'aperture minimum'
    frac_controls['aperture_min'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['RandomFracs']['VaryAperture']['maxAperture']['unit'])
    worth = docs['RandomFracs']['VaryAperture']['maxAperture']['valu']
    term = 'aperture maximum'
    frac_controls['aperture_max'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    # CorrelateAperture:
    frac_controls['aperture_alpha'] = \
        docs['RandomFracs']['CorrelateAperture']['alpha']

    cubit = str(docs['RandomFracs']['CorrelateAperture']['beta']['unit'])
    worth = docs['RandomFracs']['CorrelateAperture']['beta']['valu']
    term = 'beta'
    frac_controls['aperture_beta'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    # RandomThreshold (pressure):
    cubit = str(docs['RandomFracs']['Threshold']['refEntry']['unit'])
    worth = docs['RandomFracs']['Threshold']['refEntry']['valu']
    term = 'entry/threshold pressure'
    frac_controls['entry_pressure'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    # Pressure-Aperture correlation parameters.
    # -- Stress_limit is only defined under this option.
    if frac_controls['pressure_approach']:

        # Residual aperture.
        cubit = str(docs['RandomFracs']['PressureCorrection']['resAperture']
                    ['unit'])
        worth = (docs['RandomFracs']['PressureCorrection']['resAperture']
                 ['valu'])
        term = 'residual aperture'
        frac_controls['residual_aperture'] = \
            yam.recast_unit_to_default(cubit, worth, term)

        # maximum aperture.
        cubit = str(docs['RandomFracs']['PressureCorrection']['wideAperture']
                    ['unit'])
        worth = (docs['RandomFracs']['PressureCorrection']['wideAperture']
                 ['valu'])
        term = 'wide aperture'
        frac_controls['wide_aperture'] = \
            yam.recast_unit_to_default(cubit, worth, term)

        # Stress limit.
        cubit = str(docs['RandomFracs']['PressureCorrection']['stressLimit']
                    ['unit'])
        worth = (docs['RandomFracs']['PressureCorrection']['stressLimit']
                 ['valu'])
        term = 'stress limit'
        frac_controls['stress_limit'] = \
            yam.recast_unit_to_default(cubit, worth, term)

        # Theta (no unit).
        frac_controls['theta_aperture'] = \
            docs['RandomFracs']['PressureCorrection']['theta']

    return frac_controls


def interpret_frac_parameters(docs, frac_controls):
    """Control conversion of YAML parameters into fracture controls.

    Parameters
    ----------
    docs = (stream) YAML array input
    frac_controls = (dict) control parameter dictionary

    Returns
    -------
    frac_controls = (dict) (revised)
    """
    # Already defined: frac_controls[ = dict().

    # Define vertical connectivity.
    frac_controls['connect_factor'] = float(docs['Connectivity']['geometric'])

    # Define controls:
    frac_controls['random_approach'] = docs['Controls']['randomFracApproach']
    frac_controls['user_approach'] = docs['Controls']['userFracApproach']
    frac_controls['plot_fractures'] = docs['Controls']['plotFracApproach']
    frac_controls['correlate_approach'] = \
        docs['Controls']['apertureLengthApproach']
    frac_controls['pressure_approach'] = \
        docs['Controls']['pressureApproach']

    # Random fractures input:
    if frac_controls['random_approach']:
        frac_controls = translate_random_parameters(docs, frac_controls)
        frac_controls = translate_other_random_parameters(docs, frac_controls)

    # Input User file name - use default if missing.
    frac_controls['user_input_file'] = docs['InputFracs']['Filename']
    if frac_controls['user_input_file'] in MISSING:
        frac_controls['user_input_file'] = scfg.USER_NAME

    # Input User definitions.
    if frac_controls['user_approach']:
        cubit = str(docs['InputFracs']['Threshold']['refAperture']['unit'])
        worth = docs['InputFracs']['Threshold']['refAperture']['valu']
        term = 'user reference aperture'
        frac_controls['ref_user_aperture'] = \
            yam.recast_unit_to_default(cubit, worth, term)

        cubit = str(docs['InputFracs']['Threshold']['refPressure']['unit'])
        worth = docs['InputFracs']['Threshold']['refPressure']['valu']
        term = 'user reference pressure'
        frac_controls['ref_user_pressure'] = \
            yam.recast_unit_to_default(cubit, worth, term)

    # Input matrix-permeability:
    cubit = str(docs['RockMatrix']['Permeability']['aveMatrix']['unit'])
    worth = docs['RockMatrix']['Permeability']['aveMatrix']['valu']
    term = 'average matrix permeability'
    frac_controls['rock_perm_ave'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['RockMatrix']['Permeability']['stdDevMatrix']['unit'])
    worth = docs['RockMatrix']['Permeability']['stdDevMatrix']['valu']
    term = 'matrix standard deviation permeability'
    frac_controls['rock_perm_dev'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['RockMatrix']['Permeability']['minMatrix']['unit'])
    worth = docs['RockMatrix']['Permeability']['minMatrix']['valu']
    term = 'matrix minimum permeability'
    frac_controls['rock_perm_min'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['RockMatrix']['Permeability']['maxMatrix']['unit'])
    worth = docs['RockMatrix']['Permeability']['maxMatrix']['valu']
    term = 'matrix maximum permeability'
    frac_controls['rock_perm_max'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    # Matrix-Reference:
    cubit = str(docs['RockMatrix']['Threshold']['matrixPerm']['unit'])
    worth = docs['RockMatrix']['Threshold']['matrixPerm']['valu']
    term = 'reference threshold permeability'
    frac_controls['ref_matrix_perm'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    cubit = str(docs['RockMatrix']['Threshold']['matrixPress']['unit'])
    worth = docs['RockMatrix']['Threshold']['matrixPress']['valu']
    term = 'reference threshold pressure'
    frac_controls['ref_matrix_threshold'] = \
        yam.recast_unit_to_default(cubit, worth, term)

    return frac_controls


#
# -----------------------------------------------------------------------------
# - End of module
