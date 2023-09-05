#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions to compute fault permeability.

Author: Ernest N. Lindner
Date: 08/23/2022

Module Name
    flt_perm

Contents (11)
    compute_aperture(fault_controls, stress)
    aperture_correct(top_pressure, base_pressure, aperture)
    current_average_pressure(fault_controls, co2_base_pressure)
    saturation_temperature(injection_time)
    saturation_line_pressure(end_temperature)
    compute_effective_saturation(co2_saturation, brine_residual, co2_residual)
    compute_brine_pressure(base_co2_pressure, capillary_pressure)
    compute_capillary_pressure(zeta, beta, normal_saturation,
                               bubbling_pressure)
    brine_relative_perm(rel_model, effective_saturation, fault_controls)
    co2_relative_perm(rel_model, effective_saturation, fault_controls)
    ---
    soluble_co2(self, fault_controls, brine_flow)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import math                     # for power ops

import flt_units as funit         # For unit conversion
import flt_file as fileop       # For file operations

# Permeability and heterogeneity constants.
MIN_SATURATION = 0.001              # Numerical minimum limit
MAX_SATURATION = 0.999              # Numerical maximum limit

# Saturation Line constants
MIN_SATURATION_TEMP = 9.0           # oC - Minimum temperature along saturation
MAX_SATURATION_TEMP = 50.0          # oC - Maximum temperature along saturation
END_TIME = 50.0                     # years - end time of saturation curve
TEMPERATURE_END_POINT = 9.54        # oC - the lowest temperature on curve
EQ_POWER = -0.35                    # Time equation exponent
EQ_COEFF = 40.0                     # Time equation constant
COEFF_A = 2.3646200E-05             # 3rd order coefficient
COEFF_B = 6.9390337E-04             # 2nd order coefficient
COEFF_C = 8.1590611E-02             # 1st order coefficient
COEFF_D = 3.53230165                # constant coefficient


def ave_pressure(top_pressure, base_pressure):
    """Compute the average pressure.

    Parameters
    ----------
    top_pressure = (float) pressure at top of flow path
    base_pressure = (float) pressure at base of flow path

    Returns
    -------
    average_pressure_1 = average pressure
    """
    average_pressure = (top_pressure + base_pressure) / 2.0

    return average_pressure


def aperture_correct(fault_controls, current_pressure):
    """Compute the change in aperture due to pressure.

    Parameters
    ----------
    fault_controls = (dict) control parameter dictionary

    Returns
    -------
    new_permeability = (float) new permeability
    new_aperture = (float) changed aperture
    """
    # Compute new aperture.
    new_aperture = compute_aperture(fault_controls, current_pressure)

    # Compute Permeability in meters.
    #   cubic law => k [permeability] - in m^2
    interim = new_aperture * funit.mm_to_m()
    squared = math.pow(interim, 2)
    new_permeability = squared / 12.0

    return new_aperture, new_permeability


def compute_aperture(fault_controls, stress):
    """Compute aperture using offset model.

    Parameters
    ----------
    frac_controls = (dic) initial aperture values
    stress = (float) effective (external) normal stress on fracture

    Returns
    -------
    aperture = (float) aperture value for normal stress
    """
    if 0.0 <= stress <= fault_controls['pa_limit_stress']:
        # Compute  general case.
        stiffness = fault_controls['pa_frac_stiffness']
        stress_factor = stress / (stress + stiffness)
        gamma_current = 1.0 - fault_controls['pa_beta'] * stress_factor
        gamma_delta = gamma_current - fault_controls['pa_gamma_ref']
        result = (gamma_delta * fault_controls['pa_vary_aperture']
                  + fault_controls['pa_initial_aperture'])
    elif stress > fault_controls['pa_limit_stress']:
        # Minimum case.
        result = fault_controls['pa_residual_aperture']
    else:
        # Maximum case.
        result = fault_controls['pa_max_aperture']

    return result


def saturation_temperature(injection_time):
    """Compute the temperature along saturation curve based on injection time.

    Parameters
    ----------
    injection_time = (float) injection time (years)

    Returns
    -------
    end_temperature = (float) temperature value along saturation line

    Notes
    -----
    1. Assumes pressure < 12 MPa - i.e.along saturation line
    """
    # Check if timer has reached end.
    if injection_time >= END_TIME:
        end_temperature = TEMPERATURE_END_POINT   # = 9.54 oC
    elif injection_time < 0.0:
        end_temperature = EQ_COEFF
    else:
        end_temperature = EQ_COEFF * math.pow((injection_time + 1.0), EQ_POWER)

    return end_temperature


def saturation_line_pressure(end_temperature):
    """Compute the pressure along saturation value at a specific temperature.

    Parameters
    ----------
    end_temperature = (float) current saturation temperature

    Returns
    -------
    pressure = (float) pressure value for temperature point

    Notes
    -----
    1. Curve fit is for a temperature range of -20 to 50 oC.
    2. See code: <curve_fit_poly>
    """
    pressure = -999.9  # Error for debugging

    # Check if proper range for code - hydrate point is 9.9 oC.
    if (end_temperature < MIN_SATURATION_TEMP           # 9.0
            or end_temperature > MAX_SATURATION_TEMP):  # 50.0
        # Code Error! Value outside defined curve fit range.
        fileop.opx_problem("Temperature is Outside Saturation Line Range!")
    else:
        # Coefficients are constants for saturation line poly.
        # Obtain pressure value - use 3rd order polynomial on temperature.
        pressure = (COEFF_A * math.pow(end_temperature, 3)
                    + COEFF_B * math.pow(end_temperature, 2)
                    + COEFF_C * end_temperature
                    + COEFF_D)

    return pressure


def compute_effective_saturation(co2_saturation, brine_residual, co2_residual):
    """Compute the effective saturation value for two-phase flow.

    Parameters
    ----------
    co2_saturation = (float) current CO2 saturation (decimal)
    brine_residual = (float) brine residual saturation (decimal)
    co2_residual = (float) CO2 residual saturation (decimal)

    Returns
    -------
    effective_saturation => normalized wet (brine) saturation
    """
    # Compute current brine saturation.
    wet_saturation = 1.0 - co2_saturation

    #  Compute effective saturation within residual limits.
    if wet_saturation <= brine_residual:
        effective_saturation = 0.0
    elif wet_saturation > 1.0 - co2_residual:
        effective_saturation = 1.0
    else:
        effective_saturation = (wet_saturation - brine_residual) \
                                / (1.0 - brine_residual - co2_residual)
    return effective_saturation


def compute_brine_pressure(fluid_pressure, capillary_pressure):
    """Compute the base brine pressures on fault horizon.

    Parameters
    ----------
    fluid_pressure = (float) pressure from reservoir data
    capillary_pressure = (float) capillary pressure

    Returns
    -------
    base_pressure = brine pressure at base of fault (MPa)
    """
    # Compute brine pressure from co2 pressure.
    new_pressure = fluid_pressure - capillary_pressure

    return new_pressure


def compute_capillary_pressure(rel_model, normal_saturation, fault_controls,
                               entry_pressure):
    """Compute the capillary pressure based on effective saturation.

    Parameters
    ----------
    rel_model = (str) relative permeability model = BC/LET
    normal_saturation = (float) normalized saturation between residual
    values = effective wetting saturation
    fault_controls = (dict) dictionary of fault parameters
    entry_pressure = (float) limit pressure (Pa) for capillary pressure (Pa)

    Returns
    -------
    capillary_pressure = (float) capillary pressure (Pa)

    Notes
    -----
    1. See Equation B-15a/b in seal manual for Brooks-Corey model.
    2. See Equations B-18 and B-19 in seal manual for LET model.
    """
    # Define nonwetting saturation.
    non_saturation = 1.0 - normal_saturation

    # Define capillary pressure based on model type.

    if 'LET' in rel_model:
        # Compute based on LET capillary terms.
        term_1 = math.pow(non_saturation, fault_controls['l_capillary'])
        term_2 = (fault_controls['e_capillary'] *
                  math.pow(normal_saturation, fault_controls['t_capillary']))
        fcap = term_1 / (term_1 + term_2)
        change = fault_controls['max_capillary'] - entry_pressure

        capillary_pressure = fcap * change + entry_pressure
    else:
        # Compute based on modified Brooks-Corey model.
        exponent = (1.0 / fault_controls['zeta'])
        if normal_saturation >= MAX_SATURATION:
            divisor = 1.0
        elif normal_saturation <= MIN_SATURATION:
            divisor = math.pow(MIN_SATURATION, exponent)
        else:
            divisor = math.pow(normal_saturation, exponent)

        capillary_pressure = (entry_pressure / divisor)

    return capillary_pressure  # in Pascals


def brine_relative_perm(rel_model, effective_saturation, fault_controls):
    """Compute the wetting relative permeabilities.

    Parameters
    ----------
    rel_model = (str) relative permeability model - either"LET" or "BC"
    effective_saturation = (float) effective saturation
    fault_controls = (dict) dictionary of fault parameters

    Returns
    -------
    wet_perm = (float) wetting relative permeability (brine)

    Notes
    -----
    1. See Equation B-8 in seal manual for Brooks-Corey model.
    2. See Equation B-16 in seal manual for LET model.
    """
    # Define relative permeability based on model type.

    if 'LET' in rel_model:
        # Define wetting relative permeability using LET model.
        term_1 = math.pow(effective_saturation, fault_controls['l_wetting'])
        non_sat = 1.0 - effective_saturation
        term_2 = (fault_controls['e_wetting'] *
                  math.pow(non_sat, fault_controls['t_wetting']))
        wet_perm = term_1 / (term_1 + term_2)
    else:
        # Define nonwetting permeability using B-C, with zeta = lambda.
        term_1 = (2.0 + 3.0 * fault_controls['zeta']) / fault_controls['zeta']
        wet_perm = math.pow(effective_saturation, term_1)

    return wet_perm


def co2_relative_perm(rel_model, effective_saturation, fault_controls):
    """Compute the nonwetting relative permeability.

    Parameters
    ----------
    rel_model = (str) relative permeability model - either"LET" or "BC"
    effective_saturation = (float) effective wetting saturation
    fault_controls = (dict) dictionary of fault parameters

    Returns
    -------
    nonwet_perm = (float) nonwetting relative permeability (co2)

    Notes
    -----
    1. See Equation B-9 in seal manual for Brooks-Corey model.
    2. See Equation B-17 in seal manual for LET model.
    """
    # Define nonwetting saturation and ratio.
    non_saturation = 1.0 - effective_saturation
    p_ratio = fault_controls['perm_ratio']

    # Define relative permeability based on model type.

    if 'LET' in rel_model:
        # Define nonwetting permeability using LET model.
        term_1 = math.pow(non_saturation, fault_controls['l_nonwet'])
        term_2 = (fault_controls['e_nonwet'] *
                  math.pow(effective_saturation, fault_controls['t_nonwet']))
        fcap = term_1 / (term_1 + term_2)
        nonwet_perm = p_ratio * fcap
    else:
        # Define nonwetting permeability using B_C, with zeta = lambda.
        term_1 = math.pow(non_saturation, 2.0)
        term_2 = (2.0 + fault_controls['zeta']) / fault_controls['zeta']
        term_3 = math.pow(effective_saturation, term_2)
        nonwet_perm = p_ratio * term_1 * (1.0 - term_3)

    return nonwet_perm


def soluble_co2(fault_controls, brine_flow):
    """Compute amount of CO2 dissolved in brine (in m3/s).

    Parameters
    ----------
    fault_controls = (dict) dictionary of parameters
    brine_flow = (float) amount of brine (mˆ3/s)

    Returns
    -------
    result = soluble CO2 (mˆ3/s)

    Note
    ----
    1. Assumes brine is 100% saturated with CO2.
    """
    # Define variables for clarity.
    soluble = fault_controls['co2_solubility']
    brine_dense = fault_controls['brine_density']
    co2_dense = fault_controls['co2_density']

    # Define soluble amount per time.
    # --> Define weight (kg) )of brine per sec.
    brine_kg = brine_dense * brine_flow   # kg - amount of brine^per sec.

    # --> Define solubility = mol/kg-brine * g/mol / g/kg  = kg-CO2/kg-brine.
    co2_dissolved = soluble * funit.co2_molar_mass() / funit.kilogram_to_gram()

    # --> Compute volume of CO2 - in mˆ3/s.
    result = co2_dissolved * brine_kg / co2_dense

    return result

#
# -----------------------------------------------------------------------------
# End of module
