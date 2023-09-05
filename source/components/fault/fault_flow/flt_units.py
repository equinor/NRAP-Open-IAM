#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module provides conversion factors and standard values.

Author: Ernest N. Lindner
Date: 08/16/2022

Module Name
    flt_units

Contents (46)
    yrs_to_seconds()
    gravity()
    mpa_to_pa()
    nacl_molar_mass()
    co2_molar_mass():
    ppm_convert()
    kilogram_to_tonne()
    kilogram_to_gram()
    microd_to_metersq()
    metersq_to_microd()
    ----
    mm_to_m()
    m_to_mm()
    mm2_to_square_m()
    ft2_to_m2()
    fahrenheit_to_celsius(degrees_f)
    psi_to_pa()
    psi_to_mpa()
    feet_to_m()
    pounds_convert_density()
    centipoise_to_pa_s()
    ----
    per_sqft_to_sqm()
    inch_to_mm()
    meters2_to_ft2()
    celsius_to_fahrenheit(degrees_c)
    mpa_to_psi()
    pa_to_psi()
    pa_to_mpa()
    meters_to_feet()
    kilo_convert_density()
    pa_s_to_centipoise()
    ----
    per_sqm_to_per_sft()
    mm_to_inch()
    tonne_to_ton()
    pound_to_kilogram()
    hydrate_temperature()
    triple_pressure()
    triple_temperature()
    supercrit_pressure()
    supercrit_temperature()
    earth_pressure(depth)
    ----
    brine_pressure(depth)
    geothermal_temperature(depth)
    supercrit_press_depth()
    depth_for_pressure(press)
    depth_for_temp(temp)
    mol_per_kg_to_mol_per_lb()

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""


def yrs_to_seconds():
    """Convert years to seconds.

    Parameters
    ----------
    N/A (years)

    Returns
    -------
    value = (float) conversion factor in seconds (for Julian calendar)

    Notes
    -----
    1. Year => 365.25 days => 3.155 760 E+07 sec.
    2. Value of days from Wikipedia.
    3. Computed in Excel.
    """
    value = 3.155760E+07  # in seconds
    return value


def gravity():
    """Provide standard gravity number - in metric terms.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) standard acceleration (g) due to gravity (m/sec2)

    Notes
    -----
    Value from Wikipedia.
    """
    value = 9.80665  # m/s2
    return value


def mpa_to_pa():
    """Convert MPa to Pa.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in pascals (Pa)
    """
    value = 1.000000E+06  # Pa
    return value


def nacl_molar_mass():
    """Provide NaCl molar mass.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) molar mass NaCl (g/mol)

    Notes
    -----
    From: https://www.webqc.org/molecular-weight-of-NaCl.html
    """
    value = 58.4428  # g/mol
    return value


def co2_molar_mass():
    """Provide CO2 molar mass.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) molar mass CO2 (g/mol)

    Notes
    -----
    From: https://sciencetrends.com/molar-mass-of-co2-carbon-dioxide/
    """
    value = 44.00964  # g/mol
    return value


def ppm_convert():
    """Convert ppm to decimal number.

    Parameters
    ----------
    N/A (ppm)

    Returns
    -------
    value = (float)  conversion factor - decimal value
    """
    value = 1.000000E-06  # decimal
    return value


def kilogram_to_tonne():
    """Convert kg to metric tonnes.

    Parameters
    ----------
    N/A (Kg)

    Returns
    -------
    value = (float)  conversion factor in metric tonne
    """
    value = 1.000000E-03  # tonnes
    return value


def kilogram_to_gram():
    """Convert kg to grams.

    Parameters
    ----------
    N/A (Kg)

    Returns
    -------
    value = (float) conversion factor in grams

    """
    value = 1.000000E+03  # g
    return value


def microd_to_metersq():
    """Convert microdarcys to m2.

    Parameters
    ----------
    N/A (microdarcy)

    Returns
    -------
    value = (float) conversion factor in a m2

    Notes
    -----
    1. From Cardarelli, F. (1997)
      Scientific Unit Conversion: A Practical Guide to Metrication.
    """
    value = 9.86923266716E-19    # m2
    return value


def metersq_to_microd():
    """Convert m2 to microdarcys.

    Parameters
    ----------
    N/A (meters^2)

    Returns
    -------
    value = (float)  conversion factor in microdarcys

    Notes
    -----
    1. From https://en.wikipedia.org/wiki/Darcy_(unit)#Conversions
    """
    value = 1.01325E+18     # microdarcys
    return value


def mm_to_m():
    """Conversion factor for mm to meters.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor value in m
    """
    value = 1.000000E-03  # m
    return value


def m_to_mm():
    """Conversion factor for meters to mm.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor value in mm
    """
    value = 1.000000E+03  # mm
    return value


def mm2_to_square_m():
    """Conversion factor for mm^2 to m2.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor value in m2
    """
    value = mm_to_m() * mm_to_m()
    return value


def ft2_to_m2():
    """Convert area in ft2 to meters2.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor value in m2

    Note
    ----
    1. URL: https://www.unitconverters.net/area/square-feet-to-square-meter.htm
    """
    value = 0.09290304
    return value


def fahrenheit_to_celsius(degrees_f):
    """Convert temperature in fahrenheit to temperature in celsius.

    Parameters
    ----------
    degrees_f = (float) temperature in fahrenheit

    Returns
    -------
    value = (float) temperature in degrees celsius

    Note
    ----
    1. URL: https://en.wikipedia.org/wiki/Fahrenheit
    """
    value = (degrees_f - 32.0) * (5.0/9.0)
    return value


def psi_to_pa():
    """Convert pounds per square inch (psi) to pascals.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor value in pascals

    Note
    ----
    1. URL: https://www.unitconverters.net/pressure/psi-to-pascal.htm
    """
    value = 6894.7572931783
    return value


def psi_to_mpa():
    """Convert pounds per square inch (psi) to megapascals.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor value in megapascals

    Note
    ----
    1. URL: https://en.wikipedia.org/wiki/Pound_per_square_inch
    """
    value = 6.894757293168E-3
    return value


def feet_to_m():
    """Convert feet to meters.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in meters

    Note
    ----
    1. URL: https://en.wikipedia.org/wiki/Foot_(unit)
    """
    value = 0.3048
    return value


def pounds_convert_density():
    """Convert pound-mass per cubic foot to kilogram per cubic meter.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in kg/m3

    Note
    ----
    1. Value is approximate (different authors differ on last two digits).
    2. URL: https://www.unitconverters.net/density/pound-cubic-foot-to-kilogram
       -cubic-meter.htm
    """
    value = 16.018463374
    return value


def centipoise_to_pa_s():
    """Convert viscosity in centipoise to Pa-s.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in Pa-s

    Note
    ----
    1. URL: https://en.wikipedia.org/wiki/Poise_(unit)
    """
    value = 1.0E-3
    return value


def per_sqft_to_sqm():
    """Convert per square feet to per square meters.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in per meterˆ2

    Note
    ----
    1. URL: https://www.thecalculatorsite.com/conversions/area/square-meters
           -to-square-feet.php
    """
    value = 10.76391041671
    return value


def inch_to_mm():
    """Convert inches to millimeter.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in mm

    Note
    ----
    1. URL: https://en.wikipedia.org/wiki/Inch
    """
    value = 25.4
    return value


def meters2_to_ft2():
    """Convert area in metersˆ2 to ft2.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in ft2

    Note
    ----
    1. URL: https://www.convertunits.com/from/square+meters/to/square+feet
    """
    value = 10.76391041671
    return value


def celsius_to_fahrenheit(degrees_c):
    """Convert temperature in celsius to temperature in fahrenheit.

    Parameters
    ----------
    degrees_f = (float) temperature in celsius

    Returns
    -------
    value = (float) temperature in degrees fahrenheit

    Note
    ----
    1. URL: https://en.wikipedia.org/wiki/Conversion_of_units_of_temperature
    """
    value = (degrees_c * 9.0/5.0) + 32.0
    return value


def mpa_to_psi():
    """Convert MPa to pounds per square inch (psi).

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in psi

    Note
    ----
    1. URL: https://www.unitconverters.net/pressure/psi-to-pascal.htm
    """
    value = 145.03773773
    return value


def pa_to_psi():
    """Convert Pa to pounds per square inch (psi).

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in psi

    Note
    ----
    1. URL: https://www.unitconverters.net/pressure/psi-to-pascal.htm
    """
    value = 1.4503773773E-4
    return value


def pa_to_mpa():
    """Convert Pa to MPa.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in megapascals (MPa)
    """
    value = 1.000000E-06  # MPa
    return value


def meters_to_feet():
    """Convert meters to feet.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in feet

    Note
    ----
    1. URL: https://www.unitconverters.net/length/meters-to-feet.htm
    """
    value = 3.280839895
    return value


def kilo_convert_density():
    """Convert kilogram per cubic meter to pound-mass per cubic foot.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in pounds per cubic foot (pcf)

    Note
    ----
    1. Value is approximate (different authors differ on last few digits).
    2. URL: https://www.convertunits.com/from/kg/m3/to/lb/ft3
    """
    value = 0.062427960576145
    return value


def pa_s_to_centipoise():
    """Convert viscosity in Pa-s to centipoise.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in centipoise

    Note
    ----
    1. URL: https://www.convertunits.com/from/Pa-s/to/centipoise
    """
    value = 1.0E+3
    return value


def per_sqm_to_per_sft():
    """Convert per square meters to per square feet.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in per meterˆ2

    Note
    ----
    1. URL: https://www.thecalculatorsite.com/conversions/area/square-meters
           -to-square-feet.php
    """
    value = 0.09290304
    return value


def mm_to_inch():
    """Convert millimeters to inches.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in inches

    Note
    ----
    1. URL: https://www.convertunits.com/from/mm/to/inches
    """
    value = 0.03937007874015748
    return value


def tonne_to_ton():
    """Convert metric tonne to short ton (US).

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in US (short) tons

    Note
    ----
    1. URL: https://www.convertunits.com/from/metric+tonnes/to/tons
    2. One tonne is also equal to 1 megagram.
    """
    value = 1.1023113109244
    return value


def pound_to_kilogram():
    """Convert pound to kilogram.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in US (short) tons

    Note
    ----
    1. URL: https://www.convertworld.com/en/mass/pound-avoirdupois-us/
       lbs-to-kg.html
    """
    value = 0.45359237
    return value


def hydrate_temperature():
    """Provide temperature limit for P Quadruple Point - hydrate point.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) in degrees celsius

    Notes
    -----
    SLOAN, E. D.; KOH, C.; SUM, A. K., (2011)
        Natural gas hydrates in flow assurance.
        Gulf Professional Pub./Elsevier.
        (see Wikipedia, on carbon dioxide clathrate)
        = 283.0 K
    """
    value = 9.85  # oC
    return value


def triple_pressure():
    """Provide pressure limit for CO2 triple point.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) in Pa

    Notes
    -----
    From Span and Wagner (1994)
        www.nist.gov/system/files/documents/srd/jpcrd516.pdf
        = 0.51795 MPa
    """
    value = 517950.0  # Pa
    return value


def triple_temperature():
    """Provide temperature limit for CO2 triple point.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) in degrees celsius

    Notes
    -----
    From Span and Wagner (1994)
        www.nist.gov/system/files/documents/srd/jpcrd516.pdf
        = 216.592 K
    """
    value = -56.558  # oC
    return value


def supercrit_pressure():
    """Provide pressure limit for supercritcal CO2.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) in Pa

    Notes
    -----
    From Wikipedia:
        https://en.wikipedia.org/wiki/Supercritical_carbon_dioxide
    """
    value = 7377300.0  # Pa
    return value


def supercrit_temperature():
    """Provide temperature limit for supercritical CO2.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) in degrees celsius

    Notes
    -----
    NIST:
        https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=4
    """
    value = 31.03  # oC
    return value


def earth_pressure(depth):
    """Compute lithostatic pressure (Pa) at a depth.

    Parameters
    ----------
    depth = (float) depth in meters

    Returns
    -------
    value = (float) earth pressure value in Pa

    Note
    ----
    1. Gradient = Total subsurface mass gradient, kPa/m
    2. https://www.sciencedirect.com/topics/engineering/lithostatic-pressure
    """
    lithostatic_gradient = 24.0  # kPa/m
    value = depth * lithostatic_gradient * 1000.0  # Pa
    return value


def brine_pressure(depth):
    """Compute brine / hydrostatic pressure (Pa) at a depth.

    Parameters
    ----------
    depth = (float) depth in meters

    Returns
    -------
    value = (float) hydrostatic pressure value in Pa

    Note
    ----
    1. Gradient = Saline water gradient, ~10.53 kPa/m
    2. Salinity is about 80,000 ppm of dissolved solids at 25oC
    """
    hydrostatic_gradient = 10.53  # kPa/m
    value = depth * hydrostatic_gradient * 1000.0  # Pa
    return value


def geothermal_temperature(depth):
    """Compute geothermal (earth) temperature at a depth.

    Parameters
    ----------
    depth = (float) depth in meters

    Returns
    -------
    value = (float) subsurface temperature (oC)

    Note
    ----
    1. Geothermal gradient about 25.0 oC/km
    2. Linear equation includes surface temperature, assumed = 20.0 oC
    3. URL: https://www.sciencedirect.com/topics/earth-and-planetary-sciences/
            geothermal-gradient
    """
    surface_temperature = 20.0  # = 68 oF
    thermal_gradient = 25.0  # oC/km
    value = ((depth / 1000.0) * thermal_gradient) + surface_temperature
    return value


def supercrit_press_depth():
    """Compute depth of supercritical pressure based on gradient.

    Parameters
    ----------
    N/A

    Returns
    -------
    depth = (float) depth in meters.

    Note
    ----
    1. Gradient = Total subsurface mass gradient, kPa/m
    2. https://www.sciencedirect.com/topics/engineering/lithostatic-pressure
    """
    hydrostatic_gradient = 1.053E+04   # Pa/m

    # Compute value from gradient - supercritical pressure is in Pa
    depth = (supercrit_pressure() / hydrostatic_gradient)
    return depth


def depth_for_pressure(press):
    """Compute depth of a pressure value based on gradient.

    Parameters
    ----------
    press = (float) pressure value - Pa

    Returns
    -------
    depth = (float) depth in meters

    Note
    ----
    1. https://www.sciencedirect.com/topics/engineering/lithostatic-pressure
    """
    depth = press / 1.053E+04
    return depth


def depth_for_temp(temp):
    """Compute depth of a temperature (oC) based on gradient.

    Parameters
    ----------
    temp = (float) temperature in oC

    Returns
    -------
    depth = (float) depth in meters.

    Notes
    -----
    1. Geothermal gradient 25.0 oC/km (Pruess (2203) assumes 30 oC/km surface.
    2. Linear equation includes surface temperature, assumed = 15.0 oC
       Based on Pruess(2003).
    3. URL: https://www.sciencedirect.com/topics/earth-and-planetary-sciences/
       geothermal-gradient
    4. DiPietro, J. A. Chapter 20 - Keys to the Interpretation of Geological
       History. In: Landscape Evolution in the United States, An Introduction
       to the Geography, Geology, and Natural History. Part III, 1st ed.
       2013, 327-344.
    """
    surface_temp = 15.0  # oC (or 59 oF)
    thermal_gradient = 25.0     # oC/km

    # Compute value from gradient.
    depth = (temp - surface_temp) / thermal_gradient
    return depth


def mol_per_kg_to_mol_per_lb():
    """Convert solubility from metric to US units.

    Parameters
    ----------
     N/A

    Returns
    -------
    mol/lb = (float) solubility in US units.

    Note
    ----
    1. Inverse value at: https://www.convert-measurement-units.com/conversion
       - calculator.php?type=stoffmenge
    """
    value = 0.45359237
    return value


#
# -----------------------------------------------------------------------------
# - End of module
