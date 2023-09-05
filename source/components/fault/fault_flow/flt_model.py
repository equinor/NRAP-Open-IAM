#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains definitions of Classes Pointe, Segment and Plate.

Author: Ernest N. Lindner
Date: 08/23/2022

Module Name
    flt_model

Total Routines (20 + 4 = 24)
Contents (6 + 6 + 8 = 20)
    ** Class Pointe **** (6)
        __init__(self, x_value, y_value) *
        __repr__(self)  - NOT USED
        replica(self)  - NOT USED
        distant(self, end_pt)  - NOT USED
        arthro(self, theta, radius) *
        get_coords(self) *
    ** Class Segment **** (6)
        __init__(self, start_pt, end_pt)
        replica(self)  - NOT USED
        length(self) *
        trend(self) - NOT USED
        slope_params(self) - NOT USED
        get_points(self)  - NOT USED
    ** Class Plate **** (8)
        __init__(self, point_1, point_2)
        extent(self) - NOT USED
        print_plate(self) *
        list_data(self, numbr)
        plate_array(self) - NOT USED
        get_value(self, select) - NOT USED
        compute_flow(self, base_co2_pressure, base_co2_saturation,
                     fault_controls)
        complex_tail(self, co2_base_pressure, co2_base_saturation,
                     fault_controls, trunk_co2_flow, trunk_brine_flow)

RELATED FUNCTIONS (4)
    round_half_up(number, decimals=0)
    convert_flows(co2_flow, brine_flow, current_time, past_time,
                  fault_controls)
    flow_terms(fault_controls, initial_permeability, relative_perm, viscosity)
    darcy(fault_controls, base_pressure, top_pressure, effective_terms,
          density, aperture, length)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import math                     # For trigonometry, infinite, hypotenuse
import copy                     # For <deepcopy>

import flt_perm as perm         # Permeability
import flt_tail as ftal         # Fault computations
import flt_units as funits      # For unit conversion

# Constants
V_SMALL_ZERO = 1.0E-20          # Very small number in division by zero
V_SMALL_ATAN = 1.0E-04          # Very small number in atan
SAT_LIMIT = 1.0e-01             # % - saturation transition for solubility
ECHO = False                    # Debugging printout of trunk and tail flows


#   ---------------------------------------------------------------------------
class Pointe:
    """A class of geometric points in 2D."""

    point_count = 0

    def __init__(self, x_value, y_value):
        """Create a point in 2D cartesian coordinates.

        Parameters
        ----------
        x_value = (float) coordinate along N-S axis (m)
        y_value = (float) coordinate along E-W axis (m)
        """
        self.x_coord = x_value
        self.y_coord = y_value
        Pointe.point_count += 1
        # end

    def __repr__(self):
        """Provide a command string for printing."""
        return f'Pointe(x={self.x_coord}, y={self.y_coord})'

    def replica(self):
        """Return a copy of a point.

        Parameters
        ----------
        N/A

        Returns
        -------
        new_point = (class) a new point
        """
        new_point = copy.copy(self)
        return new_point

    def distant(self, end_pt):
        """Find the distance from current point to another point.

        Parameters
        ----------
            end_pt = (class) another point

        Returns
        -------
        distance

        Notes
        -----
        Python: hypot = sqrt(x*x + y*y).
        """
        delta_x = end_pt.x_coord - self.x_coord
        delta_y = end_pt.y_coord - self.y_coord
        return math.hypot(delta_x, delta_y)

    def arthro(self, theta, radius):
        """Create a point from current point at a radial distance away.

        Parameters
        ----------
        theta = (float) angle from north (degrees) +/- 180
        radius = (float) straight distance from point to new

        Returns
        -------
        new_line = (class) line from point to radial distance.
        """
        alpha = math.radians(theta)
        x_add = radius * math.sin(alpha)
        y_add = radius * math.cos(alpha)

        new_point = Pointe((self.x_coord + x_add), (self.y_coord + y_add))
        return new_point

    def get_coords(self):
        """Get the X and Y coordinates of a point.

        Parameters
        ----------
        N/A

        Returns
        -------
        x_coord = (float) coordinate along N-S axis (m)
        y_coord = (float) coordinate along E-W axis (m)
        """
        return self.x_coord, self.y_coord


# -----------------------------------------------------------------------------
class Segment:
    """A class of line segments in 2D."""

    def __init__(self, start_pt, end_pt):
        """Create a line in 2D cartesian coordinates with two points.

        Parameters
        ----------
        start_pt = (class) start point (left)
        end_pt = (class) end point (right)

        Notes
        -----
        Convention is that left point is smaller value along X-axis.
        """
        self.pt_1 = start_pt
        self.pt_2 = end_pt
        # end

    def replica(self):
        """Provide a copy of a line.

        Parameters
        ----------
        N/A

        Returns
        -------
        new_point= (class) copy of line.
        """
        new_line = copy.deepcopy(self)
        return new_line

    def length(self):
        """Compute the length of line.

        Parameters
        ----------
        N/A

        Returns
        -------
        distance = (float) distance between two point

        Notes
        -----
        Python: hypot = sqrt(x*x + y*y).
        """
        delta_x = self.pt_2.x_coord - self.pt_1.x_coord
        delta_y = self.pt_2.y_coord - self.pt_1.y_coord
        distance = math.hypot(delta_x, delta_y)
        return distance

    def trend(self):
        """Compute strike of line from North.

        Parameters
        ----------
        N/A

        Returns
        -------
        angle = (float) angle from north (deg), -> must be within < +/-90 deg.
        """
        delta_x = self.pt_2.x_coord - self.pt_1.x_coord
        delta_y = self.pt_2.y_coord - self.pt_1.y_coord

        # keep angle positive, in range +/- 90 degrees.
        if delta_y < 0.0:
            delta_y = -delta_y
            delta_x = -delta_x

        # prevent crash if delta_y = 0.
        if delta_y < V_SMALL_ATAN:
            north_angle = 90.0
        else:
            north_angle = math.degrees(math.atan(delta_x / delta_y))
        return north_angle

    def slope_params(self):
        """Define slope and intercept (m, b) for 2D line.

        Parameters
        ----------
        N/A

        Returns
        -------
        m = (float) slope
        b = (float) line constant

        Notes
        -----
        1. Equation: < y = mx + b >; where m = slope and b = beta.
        """
        delta_x = self.pt_2.x_coord - self.pt_1.x_coord
        delta_y = self.pt_2.y_coord - self.pt_1.y_coord

        if abs(delta_x) < V_SMALL_ZERO:
            incline = math.inf
            beta = self.pt_1.y_coord
        else:
            incline = delta_y / delta_x
            if abs(incline) < V_SMALL_ZERO:
                incline = 0.0
            beta = self.pt_1.y_coord - incline * self.pt_1.x_coord
        return incline, beta

    def get_points(self):
        """Get the points of a 2D segment.

        Parameters
        ----------
        N/A

        Returns
        -------
        points = (class) Points defining line.
        """
        return self.pt_1, self.pt_2


# -----------------------------------------------------------------------------
class Plate:
    """A class of plate (parts/sections of a fault)."""

    def __init__(self, point_1, point_2, a_strike, a_dip):
        """Initialize fault plate.

        Parameters
        ----------
        profile => Segment
            pt_1 = Pointe - Left
            pt_2 = Pointe - Right
        aperture = fracture aperture (mm)
        strike = fracture strike (deg)
        dip = inclination from horizontal
        perm = fracture permeability
        entry = threshold pressure (Pa)
        """
        self.profile = Segment(point_1, point_2)
        self.strike = a_strike
        self.dip = a_dip
        self.aperture = 0.0
        self.perm = 0.0
        self.entry = 0.0
        # end

    def extent(self):
        """Get the length of a plate.

        Parameters
        ----------
        N/A

        Returns
        -------
        result = (float) length of line
        """
        result = self.profile.length()
        return result

    def print_plate(self, dest):
        """Print instance properties of a plate.

        Parameters
        ----------
        dest = destination file

        Returns
        -------
        N/A
        """
        # print(" Fracture Parameters for instance are:", file=dest)
        start = (
            f'\n   Start_point    = [{self.profile.pt_1.x_coord},'
            f' {self.profile.pt_1.y_coord}] m'
        )
        end = (
            f'\n   End_point       = [{self.profile.pt_2.x_coord},'
            f' {self.profile.pt_2.y_coord}] m'
        )

        print(start, file=dest)
        print(end, file=dest)
        print(f'   Strike         = {self.strike} deg', file=dest)
        print(f'   Dip            = {self.dip} deg', file=dest)
        print(f'   Aperture       = {self.aperture} mmm', file=dest)
        print(f'   Permeability   = {self.perm: 7.4e} m^2', file=dest)
        print(f'   Threshold      = {self.entry} Pa', file=dest)
        # return None

    def list_data(self, numbr):
        """Construct a list of plate data - tabular.

        Parameters
        ----------
        numbr = (int) line number

        Returns
        -------
        txtr = (list) parameter data
        """
        txtr = (numbr,
                self.profile.pt_1.x_coord,
                self.profile.pt_1.y_coord,
                self.profile.pt_2.x_coord,
                self.profile.pt_2.y_coord,
                self.strike,
                self.dip,
                self.aperture,
                self.perm,
                self.entry)
        return txtr

    def plate_array(self):
        """Provide plate parameters in a list.

        Parameters
        ----------
        N/A

        Returns
        -------
        result = (list) an instance of fracture/plate data
        """
        result = [self.profile.pt_1.x_coord, self.profile.pt_1.y_coord,
                  self.profile.pt_2.x_coord, self.profile.pt_2.y_coord,
                  self.aperture, self.strike, self.perm, self.entry]
        return result

    def get_value(self, select):
        """Return a parameter value of a plate depending on code number.

        Parameters
        ----------
        select = (int) parameter number code

        Returns
        -------
        result = a specific parameter value
        """
        if select == 0:
            result = self.profile.pt_1.x_coord
        elif select == 1:
            result = self.profile.pt_1.y_coord
        elif select == 2:
            result = self.strike
        elif select == 3:
            result = self.dip
        elif select == 4:
            result = self.aperture
        elif select == 5:
            result = self.perm
        elif select == 6:
            result = self.entry
        else:
            result = "ERROR!!"
        return result

    def compute_flow(self, base_co2_pressure, base_co2_saturation,
                     fault_controls):
        """Compute CO2 flow rate through each fault part / parallel plate.

        Parameters
        ----------
        base_co2_pressure = (float) CO2 pressure at base of fault (Pa)
        base_co2_saturation = (float) CO2 saturation at base of fault (--)
        top_brine_pressure = (float) brine pressure at top of fault (Pa)
        fault_controls = (dict) dictionary of parameters

        Return
        ------
        co2_flow = (float) CO2 flux (rate - mˆ3/sec)
        brine_flow = (float) brine flux (rate - m3/sec)
        """
        # ---------------------------------------------------------------------
        # CO2 computations:

        # Compute the effective wet saturation.
        effective_saturation = perm.compute_effective_saturation(
            base_co2_saturation, fault_controls['resid_brine'],
            fault_controls['resid_co2'])

        # Compute the capillary pressure (Pa).
        capillary_press = perm.compute_capillary_pressure(
            fault_controls['relative_model'], effective_saturation,
            fault_controls, self.entry)

        # ---------------------------------------------------------------------
        # CO2 computations:
        # -> Compute only if pressure is large enough.
        if base_co2_pressure > self.entry:

            # Compute the CO2 pressure at top (Pa).
            top_pressure = (capillary_press
                            + fault_controls['aquifer_pressure'])
            base_pressure = base_co2_pressure

            # Re-evaluate permeability for pressure_approach.
            if fault_controls['pressure_approach']:
                current_pressure = perm.ave_pressure(top_pressure,
                                                     base_co2_pressure)
                self.aperture, self.perm = \
                    perm.aperture_correct(fault_controls, current_pressure)

            # Get the CO2 relative permeability factor.
            relative_perm = \
                perm.co2_relative_perm(fault_controls['relative_model'],
                                       effective_saturation,
                                       fault_controls)

            # Compute effective terms - & convert to m2 for permeability.
            effective_terms = flow_terms(fault_controls, self.perm,
                                         relative_perm,
                                         fault_controls['co2_viscosity'])

            # Compute darcy velocity (v) and flow (Q) for CO2.
            co2_flow = darcy(fault_controls, base_pressure, top_pressure,
                             effective_terms, fault_controls['co2_density'],
                             self.aperture, fault_controls['travel_distance'])

            # Adjust flow based on fault orientation wrt injection well.
            co2_flow *= fault_controls['fault_orient_effect']

        else:
            # If too small pressure, no CO2/brine flow.
            co2_flow = 0.0

        # ------------------------------------------------------------------
        # Brine computations:
        #   -> Effective wet saturation - from above.
        #   -> Compute the capillary pressure - from above.

        # Compute effective brine pressures (Pa).
        base_pressure = perm.compute_brine_pressure(base_co2_pressure,
                                                    capillary_press)
        top_pressure =  \
            perm.compute_brine_pressure(fault_controls['aquifer_pressure'],
                                        capillary_press)

        if abs(base_pressure - top_pressure) > V_SMALL_ZERO:

            # Get brine relative permeability.
            relative_perm = \
                perm.brine_relative_perm(fault_controls['relative_model'],
                                         effective_saturation,
                                         fault_controls)

            # Compute effective terms - & convert to m2 for permeability.
            effective_terms = flow_terms(fault_controls, self.perm,
                                         relative_perm,
                                         fault_controls['brine_viscosity'])

            # Compute darcy velocity (v) and flow (Q) for brine -> m3/sec.
            brine_flow = darcy(fault_controls, base_pressure, top_pressure,
                               effective_terms,
                               fault_controls['brine_density'],
                               self.aperture,
                               fault_controls['travel_distance'])

            # Adjust flow based on fault inclination.
            brine_flow *= fault_controls['fault_orient_effect']

            # Add dissolved CO2 in brine to flow @ 100% brine saturated
            if base_co2_pressure > self.entry:
                if base_co2_saturation > SAT_LIMIT and brine_flow > 0.0:
                    co2_flow += perm.soluble_co2(fault_controls, brine_flow)

        else:
            # If too small pressure, no brine flow.
            brine_flow = 0.0

        return co2_flow, brine_flow

    def complex_tail(self, co2_base_pressure, co2_base_saturation,
                     fault_controls, trunk_co2_flow, trunk_brine_flow):
        """Compute CO2 flow rate of tail for each fault segment.

        Parameters
        ----------
        co2_base_pressure = (float) CO2 pressure at base of fault (Pa)
        co2_base_saturation = (float) CO2 saturation at base of fault (--)
        fault_controls = (dict) dictionary of parameters
        trunk_co2_flow = (float) CO2 flow rate for part 1 of profile
        trunk_brine_flow = (float) brine flow rate for part 1 of profile

        Return
        ------
        final_co2_flow = (float) final CO2 flux (rate - mˆ3/sec)
        final_brine_flow = (float) final brine flux (rate - m3/sec)
        """
        # ----------------------------------------------------------------
        # CO2 computations:

        # Compute average properties for tail.
        c_density, c_viscosity, b_density, b_viscosity = \
            ftal.compute_average_properties(fault_controls)

        # Compute flow profile length.
        complex_length = ftal.compute_tail_length(fault_controls)

        # Compute the effective wet saturation (assume same as base).
        effective_saturation = perm.compute_effective_saturation(
            co2_base_saturation, fault_controls['resid_brine'],
            fault_controls['resid_co2'])

        # Compute the capillary pressure (Pa).
        capillary_press = perm.compute_capillary_pressure(
            fault_controls['relative_model'], effective_saturation,
            fault_controls, self.entry)

        # Compute the CO2 pressure at top (Pa).
        top_pressure = (capillary_press
                        + fault_controls['crux_pressure'])
        base_pressure = 0.0
        if fault_controls['profile_type'] == 1:
            # Model 1: 12 MPa depth is deepest.
            base_pressure = fault_controls['trans_point_press']

        elif fault_controls['profile_type'] == 2:
            # Model 2: Supercritical depth is deepest.
            base_pressure = funits.supercrit_pressure()

        # Get the CO2 relative permeability factor.
        relative_perm = \
            perm.co2_relative_perm(fault_controls['relative_model'],
                                   effective_saturation,
                                   fault_controls)

        # Compute effective terms - & convert to m2 for permeability.
        effective_terms = flow_terms(fault_controls, self.perm,
                                     relative_perm,
                                     c_viscosity)

        # Compute darcy velocity (v) and flow (Q) for CO2.
        tail_co2_flow = darcy(fault_controls, base_pressure,
                              top_pressure,
                              effective_terms,
                              c_density,
                              self.aperture, complex_length)

        # --------------------------------------------------------------------
        # Brine computations:
        #   -> Effective wet saturation - from above.
        #   -> Compute the capillary pressure - from above

        # Compute effective brine pressures (Pa).
        base_pressure = perm.compute_brine_pressure(co2_base_pressure,
                                                    capillary_press)
        top_pressure = \
            perm.compute_brine_pressure(fault_controls['aquifer_pressure'],
                                        capillary_press)

        # Get brine relative permeability.
        relative_perm = \
            perm.brine_relative_perm(fault_controls['relative_model'],
                                     effective_saturation,
                                     fault_controls)

        # Compute effective terms - & convert to m2 for permeability.
        effective_terms = flow_terms(fault_controls, self.perm,
                                     relative_perm,
                                     b_viscosity)

        # Compute darcy velocity (v) and flow (Q) for brine -> m3/sec.
        tail_brine_flow = darcy(fault_controls, base_pressure,
                                top_pressure, effective_terms,
                                b_density,
                                self.aperture, complex_length)

        # Select the smaller values.
        if trunk_co2_flow <= 0.0:
            co2_flow = tail_co2_flow
        else:
            co2_flow = min(tail_co2_flow, trunk_co2_flow)

        if trunk_co2_flow <= 0.0:
            brine_flow = tail_brine_flow
        else:
            brine_flow = min(tail_brine_flow, trunk_brine_flow)

        if ECHO:
            print(f'Complex Tail Flow = {tail_co2_flow}')
            print(f'Complex Trunk Flow = {trunk_co2_flow} \n')

        return co2_flow, brine_flow

    # --------------------------------------------------------------------------
    #   RELATED FUNCTIONS


def convert_flows(co2_flow, brine_flow, current_time, past_time,
                  fault_controls):
    """Convert flows to mass flux.

    Parameters
    ----------
    co2_flow = (NumPy array) CO2 flow
    brine_flow = (NumPy array) brine flow
    current_time = (float) time for current step
    past_time = (float) time for previous step
    fault_controls = (dict) fault parameters (dictionary)

    Returns
    -------
    co2_flow = corrected CO2 values
    brine_flow = corrected brine values
    """
    # Convert time to years and volume rates into mass flows;
    # -- Intervals are in years, rate is in seconds.
    interval = current_time - past_time
    co2_flow *= interval * funits.yrs_to_seconds()  # kg per year
    co2_flow *= fault_controls['co2_density'] * funits.kilogram_to_tonne()

    brine_flow *= interval * funits.yrs_to_seconds()  # kg per year
    brine_flow *= fault_controls['brine_density'] * funits.kilogram_to_tonne()

    return co2_flow, brine_flow


def round_half_up(valu, decimals=0):
    """Round-up a floating point to a defined number of digits.

    Parameters
    ----------
    valu = (float) number to be rounded
    decimals = (int) number of digits per the decimal

    Returns
    -------
    results = rounded number
    """
    multiplier = math.pow(10.0, decimals)
    results = math.floor(valu * multiplier + 0.5) / multiplier

    return results


def flow_terms(fault_controls, initial_permeability, relative_perm, viscosity):
    """Compute relative permeability and effective terms.

    Parameters
    ----------
    fault_controls = (dict) fault parameters (dictionary)
    initial_permeability = (float) starting perm. before corrections (mD)
    relative_perm = (float) relative permeability of fluid (wetting or non).
    viscosity = (float) fluid viscosity

    Returns
    -------
    terms = terms for darcy equation (permeability in m^2)
    """
    # Get effective terms => k/u * corrections.
    sgr_effect = fault_controls['sgr_perm_correction']
    state_effect = fault_controls['state_correction']
    correct = sgr_effect * state_effect * relative_perm

    terms = correct * initial_permeability / viscosity
    terms *= funits.microd_to_metersq()

    return terms


def darcy(fault_controls, base_pressure, top_pressure, effective_terms,
          density, aperture, length):
    """Use darcy equation to get CO2 flow rate through each parallel plate.

    Parameters
    ----------
    fault_controls = (dict) fault parameters (dictionary)
    base_pressure = (float) base pressure (Pa)
    top_pressure = (float) top pressure (Pa)
    effective_terms = (float) Terms including viscosity
    density = (float) fluid density (kg/m^3)
    aperture = (float) aperture of fracture (mm)
    length = (float) travel length along fault

    Return
    ------
    darcy_flow = (float) flux (rate -> mˆ3/sec)
    """
    # Compute heads - length  must be > 0!
    if length > 0.0:
        pressure_head = (base_pressure - top_pressure) / length
        hydraulic_head = density * funits.gravity()
        potential = (pressure_head - hydraulic_head)

        # No flow downwards on fault.
        potential = max(potential, 0.0)

        # Compute Darcy velocity (q) and darcy flow (Q) -> m3/sec.
        velocity = effective_terms * potential

        area = aperture * funits.mm_to_m() * fault_controls['seg_length']
        darcy_flow = area * velocity

    else:
        darcy_flow = 0.0

    return darcy_flow


#
# -----------------------------------------------------------------------------
# End of module
# -------1---------2---------3---------4---------5---------6---------7---------8
