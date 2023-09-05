#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module provides conversion factors and standard values.

Author: Ernest N. Lindner
Date: 08/11/2022

Module Name
    seal_units

Contents (43)
    yrs_to_seconds()
    yrs_to_days()
    seconds_to_yrs()
    days_to_yrs()
    gravity()
    mpa_to_pa()
    nacl_molar_mass()
    co2_molar_mass():
    ppm_convert()
    kilogram_to_tonne()
    ----
    kilogram_to_gram()
    microd_to_metersq()
    metersq_to_microd()
    mm_to_m()
    m_to_mm()
    mm2_to_square_m()
    earth_pressure(depth)
    brine_pressure(depth)
    geothermal_temperature(depth)
    acres_to_m2()
    ----
    ft2_to_m2()
    fahrenheit_to_celsius(degrees_f)
    psi_to_pa()
    psi_to_mpa()
    feet_to_m()
    pounds_convert_density()
    centipoise_to_pa_s()
    per_sqft_to_sqm()
    inch_to_mm()
    meters2_to_acres()
    ----
    meters2_to_ft2()
    celsius_to_fahrenheit(degrees_c)
    mpa_to_psi()
    pa_to_psi()
    pa_to_mpa()
    meters_to_feet()
    kilo_convert_density()
    pa_s_to_centipoise()
    per_sqm_to_per_sft()
    mm_to_inch()
    ----
    tonne_to_ton()
    ton_to_tonne()
    pound_to_kilogram()

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
    1. Year = 365.25 days => 3.155 760 E+07 sec.
    2. Day from Wikipedia.
    3. Computed in Excel.
    """
    value = 3.155760E+07  # in seconds
    return value


def yrs_to_days():
    """Convert years to days.

    Parameters
    ----------
    N/A (years)

    Returns
    -------
    value = (float) conversion factor in seconds (Julian calendar)

    Notes
    -----
    1. Year = 365.25 days.
    2. Day value is from Wikipedia.
    """
    value = 3.6525E+02  # in days
    return value


def seconds_to_yrs():
    """Convert seconds to years.

    Parameters
    ----------
    N/A (seconds)

    Returns
    -------
    value = (float) conversion factor in yrs (Julian calendar)

    Notes
    -----
    1. Year = 365.25 days => 3.155 760 E+07 sec.
    2. 1/year = 3.168 808 781 402 89 E-08 yr/sec.
    3. Day value from Wikipedia.
    4. Computed in Excel.
    """
    value = 3.16880878140289E-08  # in years
    return value


def days_to_yrs():
    """Convert days to years.

    Parameters
    ----------
    N/A (days)

    Returns
    -------
    value = (float) conversion factor in yrs (Julian calendar)

    Notes
    -----
    1. Year = 365.25 days.
    2. 1/year = 2.7378508.
    3. Day value from Wikipedia.
    4. Computed in Excel.
    """
    value = 2.7378507871321E-03  # in years
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
    value = (float) conversion factor - decimal value
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
    value = (float) conversion factor in metric tonne
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
    value = (float) conversion factor in microdarcys

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
    """Compute geothermal temperature at a depth.

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


def acres_to_m2():
    """Recast area in acres to meters^2.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor value in m2

    Note
    ----
    1. URL: https://www.thecalculatorsite.com/conversions/area/square-meters
       -to-acres.php
    """
    value = 4046.8564224
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
    value = (float) conversion factor in per meter^2

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


def meters2_to_acres():
    """Convert area in meters^2 to acres.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in m2

    Note
    ----
    1. URL: https://www.thecalculatorsite.com/conversions/area/square-meters
       -to-acres.php
    """
    value = 0.00024710538146717
    return value


def meters2_to_ft2():
    """Convert area in meters^2 to ft2.

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
    N/A (Pa)

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
    value = (float) conversion factor in per meter^2

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


def ton_to_tonne():
    """Convert short ton (US) to metric tonne.

    Parameters
    ----------
    N/A

    Returns
    -------
    value = (float) conversion factor in US (short) tons

    Note
    ----
    1. URL: https://www.checkyourmath.com/convert/weight_mass/
    short_ton_metric_ton.php
    2. One tonne is also equal to 1 megagram.
    """
    value = 0.90718474
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


#
# -----------------------------------------------------------------------------
# - End of module
