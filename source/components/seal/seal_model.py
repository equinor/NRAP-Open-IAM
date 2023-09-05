#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains Cell class for seal model computations.

Author: Ernest N. Lindner
Date: 05/27/2022

Module Name
    seal_model

Class Functions (8)
    ** Class Cell ****
    __init__(self, cell_area, cell_thickness, cell_permeability)
    set_top_depth(self, repository_depth)
    set_coord(self, x_value, y_value)
    get_value(self, select)
    compute_static_pressure(self)
    compute_model_influence(self, velocity)
    compute_flow(self, base_co2_pressure, base_co2_saturation,
                 top_brine_pressure, seal_controls)
    track_history(self, base_co2_pressure, co2_flow_rate, time_step)
Class Method (1):
    assign_controls(cls, seal_controls):
Other (1):
    convert_flows(co2_flow, brine_flow, current_time, past_time)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import seal_units as sun        # For unit conversion
import seal_perm as perm        # For permeability definitions/computations

# Other constants
LIMIT_SAT = 1.0e-01           # Small saturation - start solubility option

"""
-------------------------------------------------------------------------------
Caution
-------
    The software is provided "As Is", without warranty of any kind, express
    or implied, including but not limited to the warranties of merchantability,
    fitness for a particular purpose and noninfringement. In no event shall
    the authors or copyright holders be liable for any claim, damages or other
    liability, whether in an action of contract, tort or otherwise, arising
    from, out of or in connection with the software or the use or other
    dealings in the software.

Warranty Disclaimer
-------------------
    This software development was funded by the United States Department
    of Energy, National Energy Technology Laboratory, in part, through a site
    support contract. Neither the United States Government nor any agency
    thereof, nor any of their employees, nor the support contractor, nor any of
    their employees, makes any warranty, express or implied, or assumes any
    legal liability or responsibility for the accuracy, completeness, or
    usefulness of any information, apparatus, product, or process disclosed,
    or represents that its use would not infringe privately owned rights.
    Reference herein to any specific commercial product, process, or service
    by trade name, trademark, manufacturer, or otherwise does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the
    United States Government or any agency thereof. The views and opinions of
    authors expressed herein do not necessarily state or reflect those of the
    United States Government or any agency thereof.
-------------------------------------------------------------------------------
"""


class Cell:
    """A group of seal areal cells above the injection horizon."""

    # Class constants with default values.
    #  NOTE: Names are NOT capitalized as these will be added within IAM.

    # Fluid parameters (Could differ for each cell)
    brineDensity = 1007.0       # kg/m^3
    co2Density = 582.6          # kg/m^3
    brineViscosity = 5.596e-4   # Pa*s
    co2Viscosity = 4.387e-5     # Pa*s
    co2Solubility = 2.0E-03     # mol/kg

    # Two-phase model parameters for relative permeability
    relativeModel = 'BC'        # Relative permeability model
    brineResSaturation = 0.1    # (--)
    co2ResSaturation = 0.0      # (--)

    # Time-model parameters (Could differ for each cell)
    model = 0                   # Model type 0/1/2
    totalEffect = 0.5           # Final Perm = 50% of total
    rateEffect = 0.1            # Time function in years
    flowMinimum = 1.0E-09       # Minimum flow limit for aging
    reactivity = 8.0            # Value: 1 to 10
    carbonateContent = 5        # Shale
    clayContent = 60.0          # Clay - for a typical shale
    clayType = "smectite"       # Types: ["smectite", "illite", "chlorite"]

    def __init__(self, x_center=0.0, y_center=0.0):
        """Initialize attributes of a flow cell as part of areal grid.

        Parameters
        ----------
        x_center = x-coordinate of center (m)
        y_center = y-coordinate of center (m)

        Default Values
        --------------
        area = (float) horizontal area of cell (m^2)
        thickness = (float) vertical thickness of cell (m)
        top_seal_depth = (float) depth to top of cell (m)
        status = (int) status of cell - active=1; not used=0
        permeability = (float) total vertical permeability (microdarcys)
        history = (float) exposure time (yrs)
        influence = (float) permeability factor for exposure (1.0 = None)
        entry = (float) threshold pressure (Pa)
        """
        self.x_center = x_center        # cell center, x-coordinate
        self.y_center = y_center        # cell center, y-coordinate
        self.area = 10000.0             # area of cell (m2)
        self.thickness = 100.0          # thickness of cell (m)
        self.top_seal_depth = 1100      # depth to <top> of cell (m)
        self.status = int(1)            # created as active (=1)
        self.permeability = 1.0         # total perm. (microdarcys)
        self.history = 0.0              # exposure history (yrs)
        self.influence = 1.0            # factor on permeability
        self.entry = 5000.0             # entry pressure for cell (Pa)

    def set_top_depth(self, repository_depth):
        """Define depth to "top" of a cell.

        Parameters
        ----------
        repository_depth = (float) depth to top of repository (m)

        Returns
        -------
        top_seal_depth = (float) depth to top of cell

        Notes
        -----
        Define thickness of cell first!
        """
        self.top_seal_depth = repository_depth - self.thickness

    def set_coord(self, x_value, y_value):
        """Set new coordinates for class.

        Parameters
        ----------
        x_value = (float) x coordinate
        y_value = (float) y coordinate

        Returns
        -------
        None
        """
        self.x_center = x_value
        self.y_center = y_value

    def get_value(self, select):
        """Return a parameter value of a cell depending on code.

        Parameters
        ----------
        select = (int) number code for various parameters

        Returns
        -------
        result = (float) specific parameter value
        """
        if select == 0:
            result = self.permeability
        elif select == 1:
            result = self.thickness
        elif select == 2:
            result = self.influence
        elif select == 3:
            result = self.area
        elif select == 4:
            result = self.top_seal_depth
        elif select == 5:
            result = self.status
        elif select == 6:
            result = self.x_center
        elif select == 7:
            result = self.y_center
        elif select == 8:
            result = self.entry
        elif select == 9:
            result = self.history
        else:
            result = "ERROR!!"

        return result

    def compute_static_pressure(self, reference_pressure, reference_depth):
        """Compute static pressure at top of cell.

        Parameters
        ----------
        reference_pressure = (float) reference pressure (Pa)
        reference_depth = (float) reference depth (m)

        Returns
        -------
        pressure = (float) static pressure at top of cell (Pa)
        """
        depth_change = self.top_seal_depth - reference_depth
        pressure_change = self.brineDensity * depth_change * sun.gravity()
        pressure = reference_pressure + pressure_change

        return pressure

    def compute_model_influence(self, velocity):
        """Compute influence of alteration on permeability.

        Parameters
        ----------
        velocity = (float) CO2 flow velocity (m/s)

        Returns
        -------
        None
        """
        # Update only for active cells.
        if self.status > 0:
            if self.model == 1:          # Time model.
                # Change influence factor based on time alone.
                self.influence = perm.model_z_analysis(self.totalEffect,
                                                       self.rateEffect,
                                                       self.history)
            elif self.model == 2:        # Multi-variant model.
                prior = self.influence
                # Compute effects due to constitutive model on mineralogy.
                change_factor = perm.mineral_factor(velocity,
                                                    self.reactivity,
                                                    self.clayContent,
                                                    self.clayType,
                                                    self.carbonateContent)

                # compute z-factor from current time.
                new_influence = perm.model_z_analysis(change_factor,
                                                      self.rateEffect,
                                                      self.history)

                # Compute incremental change to influence.
                self.influence += new_influence - prior

            else:
                # No model used; no change to factor.
                pass
        else:
            # cell inactive; no change to influence.
            pass

        # return None

    def compute_flow(self, base_co2_pressure, base_co2_saturation,
                     top_brine_pressure, seal_controls):
        """Compute CO2 flow rate through a cell.

        Parameters
        ----------
        base_co2_pressure = (float) CO2 pressure at base of cell (Pa)
        base_co2_saturation = (float) CO2 saturation at base of cell (--)
        top_brine_pressure = (float) brine pressure at top of cell (Pa)
        seal_controls = (dict) seal control parameters

        Returns
        -------
        co2_flow = (float) CO2 flux (rate - kg/sec)
        brine_flow = (float) brine flux (rate - kg/sec)
        """
        # ---------------------------------------------------------------------
        # Work only on active cells.
        if self.status > 0:
            # CO2 computations:

            # -> Only pressure is large enough.
            if base_co2_pressure > self.entry:

                # Compute the effective wet saturation for cell.
                effective_saturation = perm.compute_effective_saturation(
                    base_co2_saturation, self.brineResSaturation,
                    self.co2ResSaturation)

                # Compute the capillary pressure at top of seal.
                capillary_press = perm.compute_capillary_pressure(
                    self.relativeModel, effective_saturation, seal_controls,
                    self.entry)

                # Get the CO2 relative permeability factor for the cell.
                relative_perm = \
                    perm.co2_relative_permeability(self.relativeModel,
                                                   effective_saturation,
                                                   seal_controls)

                # Compute the CO2 pressure at top of seal.
                top_pressure = capillary_press + top_brine_pressure

                # Get CO2 effective vertical permeability in (m2) for cell
                # -- (currently, self.permeability is in microdarcys!).
                effective_permeability = perm.compute_current_permeability(
                    relative_perm, self.influence, self.permeability)

                # Compute volume flux (m3/s) for cell.
                pressure_head = ((base_co2_pressure - top_pressure)
                                 / self.thickness)
                hydraulic_head = sun.gravity() * self.co2Density

                # Control reverse flows.
                delta = pressure_head - hydraulic_head
                if delta < 0:
                    print(f'Negative Head = {delta}')
                co2_flow = (self.area * effective_permeability
                            * delta / self.co2Viscosity)
            else:
                co2_flow = 0.0

            # -----------------------------------------------------------------
            # Brine computations:

            # Compute the effective wet saturation for cell.
            effective_saturation = perm.compute_effective_saturation(
                base_co2_saturation, self.brineResSaturation,
                self.co2ResSaturation)

            # Get base capillary pressure (Pa) for cell.
            capillary_press = perm.compute_capillary_pressure(
                self.relativeModel, effective_saturation, seal_controls,
                self.entry)

            # Get effective brine pressure (MPa) for cell.
            base_pressure = perm.compute_brine_pressure(base_co2_pressure,
                                                        capillary_press)
            top_pressure = perm.compute_brine_pressure(top_brine_pressure,
                                                       capillary_press)

            # Get brine relative permeability for cell.
            relative_perm = perm.brine_relative_permeability(
                self.relativeModel, effective_saturation, seal_controls)

            # Get effective permeability in m2.
            effective_permeability = \
                perm.compute_current_permeability(relative_perm,
                                                  self.influence,
                                                  self.permeability)

            # Compute brine volume flux (m3/s) for cell.
            pressure_head = ((base_pressure - top_pressure)
                             / self.thickness)
            hydraulic_head = sun.gravity() * self.brineDensity
            brine_flow = (self.area * effective_permeability
                          * (pressure_head - hydraulic_head)
                          / self.brineViscosity)

            # Add CO2 dissolved in Brine to flow - 100% brine saturated.
            if base_co2_saturation > LIMIT_SAT and brine_flow > 0.0:
                co2_flow += perm.soluble_co2(seal_controls, brine_flow)

        else:
            # If not active, no CO2/brine flow.
            co2_flow = 0.0
            brine_flow = 0.0

        return co2_flow, brine_flow

    def track_history(self, base_co2_pressure, co2_flow_rate, time_step):
        """Increase the time history of grid cells where flow occurs.

        Parameters
        ----------
        base_co2_pressure = (float) CO2 pressure at base of cell (Pa)
        co2_flow_rate = (float) CO2 flow through cell
        time_step = (float) current time step (yrs)

        Returns
        -------
        None
        """
        # Flow occurs only in active cells and when pressure is large enough.
        if self.status > 0 and base_co2_pressure > self.entry:
            # Alteration occurs with a minimum flow / flow_rate.
            if co2_flow_rate > self.flowMinimum:
                self.history += time_step

        # Return None

    @classmethod
    def assign_controls(cls, seal_controls):
        """Assign cell class parameters from dictionary.

        Parameters
        ----------
        seal_controls = (dict) seal control parameters

        Returns
        -------
        None
        """
        # Set fluid parameters.
        cls.brineDensity = seal_controls['brine_density']
        cls.brineViscosity = seal_controls['brine_viscosity']
        cls.co2Density = seal_controls['co2_density']
        cls.co2Viscosity = seal_controls['co2_viscosity']
        cls.co2Solubility = seal_controls['co2_solubility']

        # Set two-phase model parameters.
        cls.relativeModel = seal_controls['relative_model']
        cls.brineResSaturation = seal_controls['resid_brine']
        cls.co2ResSaturation = seal_controls['resid_co2']

        # Set time-model parameters.
        cls.model = seal_controls['model']
        if cls.model > 0:
            cls.totalEffect = seal_controls['total_effect']
            cls.rateEffect = seal_controls['rate_effect']
            cls.reactivity = seal_controls['reactivity']
            cls.clayContent = seal_controls['clay_content']
            cls.carbonateContent = seal_controls['carbonate_content']
            cls.clayType = seal_controls['clay_type']

        # return None


def convert_flows(co2_flow, brine_flow, current_time, past_time):
    """Convert flows from one step (1D array) in flux to 2D.

    Parameters
    ----------
    co2_flow = (array) CO2 flow (NumPy array) in m/sec
    brine_flow = (array) brine flow (Numpy Array)in m/sec
    current_time = (float) time for current step (yr)
    past_time = (float) time for previous step (yr)

    Returns
    -------
    co2_flow = (array) corrected CO2 values
    brine_flow = (array) corrected brine values
    """
    # For grid, convert seal vol. rates into mass flows for each cell.
    # -- Intervals are in years, rate is in seconds.

    interval = current_time - past_time
    co2_flow *= interval * sun.yrs_to_seconds() * Cell.co2Density
    brine_flow *= interval * sun.yrs_to_seconds() * Cell.brineDensity

    # Convert flows from kg to tonnes.
    co2_flow *= sun.kilogram_to_tonne()
    brine_flow *= sun.kilogram_to_tonne()

    return co2_flow, brine_flow


# -----------------------------------------------------------------------------
# - End of module
# -------1---------2---------3---------4---------5---------6---------7---------8
