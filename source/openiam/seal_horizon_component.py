# -*- coding: utf-8 -*-
# @author: Veronika Vasylkivska, Ernest Lindner
# September, 2022
import copy
import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt

source_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(source_folder)
sys.path.append(os.path.join(source_folder, 'components', 'seal'))

try:
    from openiam import SystemModel, ComponentModel, IAM_DIR, LookupTableReservoir
except ImportError as err:
    print('Unable to load IAM class module: '+ str(err))
import components.seal.seal_model as smod

from openiam.cfi.commons import process_dynamic_inputs

try:
    import components.seal.seal_config as sconf
    import components.seal.seal_decipher as sdeci
    import components.seal.seal_intro as sintro
    import components.seal.seal_fluids as sfluid
    import components.seal.seal_perm as perm
    import components.seal.seal_refresh as sref
    import components.seal.seal_units as sunit
    import components.seal.frac_random as fran
    from components.seal.seal_setup import SEAL_SETUP_DICT
except ImportError:
    print('\nERROR: Unable to load additional modules for Seal Horizon component\n')
    sys.exit()


SH_SCALAR_OBSERVATIONS = ['CO2_aquifer_total', 'brine_aquifer_total',
                          'mass_CO2_aquifer_total', 'mass_brine_aquifer_total']
SH_GRID_OBSERVATIONS = ['CO2_aquifer', 'brine_aquifer',
                        'mass_CO2_aquifer', 'mass_brine_aquifer']

SH_CFI_CONTROLS = {'correlate_entry_approach': ['correlateApproach', False],
                   'interpolate_approach': ['interpolateApproach', False],
                   'heterogeneity_approach': ['heterogeneityApproach', False],
                   # not available in standalone: its purpose is to provide
                   # option for variable pressure at the top of the seal
                   'vary_pressure_approach': ['varyPressureApproach', False]}


class SealHorizon(ComponentModel):
    """
    The Seal Horizon component model simulates the flow of |CO2| through a low
    permeability but fractured rock horizon (a "seal" formation) overlying the
    storage reservoir into which |CO2| is injected.

    The rock horizon is represented by a number of "cells" arranged
    (conceptually) in an arbitrary shape grid. A two-phase, relative
    permeability approach is used with Darcy’s law for one-dimensional (1D)
    flow computations of |CO2| through the horizon in the vertical
    direction. The code also allows the simulation of time-dependent processes
    that can influence such flow.

    The model is based on an earlier code, NSealR, created with GoldSim,
    and described in :cite:`Lindner2015`. A stand-alone version of this code in Python
    is also available on the NETL EDX system, described as Seal ROM.

    In the NRAP-Open-IAM control file, the type name for the component is SealHorizon.
    The following is a list of the component parameters, including the
    parameter names, units, accepted value range and the default value.

    Reference parameters for each cell:

    * **area** [|m^2|] (1 to 2.6e+5) - area of the cell (default: 10000.0)

    * **thickness** [|m|] (5 to 1000) - thickness of the cell (vertically)
      (default: 100)

    * **baseDepth** [|m|] (800 to 9500) - depth to the base of seal (default: 1100)

    * **permeability** [|m^2|] (1.0e-22 to 1.0e-15) - cell equivalent initial
      permeability (default: 1.0e-18)

    * **entryPressure** [|Pa|] (100 to 2.0e+6) - entry threshold pressure
      that controls flow into rock (default: 5000)

    Distribution parameters for thickness of the seal layer

    * **aveThickness** [|m|] (10 to 1000) - mean of the truncated normal
      distribution for thickness (default: 100)

    * **stdDevThickness** [|m|] (0 to 500) - standard deviation
      of the thickness distribution (default: 0)

    * **minThickness** [|m|] (5 to 1000) - minimum thickness; this value
      truncates the distribution and limits lower values (default: 75)

    * **maxThickness** [|m|] (10 to 1000) - maximum thickness; this value
      truncates the distribution and limits higher values (default: 125)

    Note: The setup of the four distribution parameters above is not yet implemented
    in the control file interface or GUI of NRAP-Open-IAM and available only in the script
    interface.

    Distribution parameters for permeability of the seal layer:

    * **avePermeability** [|m^2|] (1.0e-22 to 1.0e-16) - mean total vertical
      permeability of a lognormal distribution; equivalent value
      for fractured rock (default: 2.5e-16)

    * **stdDevPermeability** [|m^2|] (0 to 1.0e-17) - standard deviation
      of the total vertical permeability distribution (default: 0.0)

    * **minPermeability** [|m^2|] (1.0e-24 to 1.0e-17) - minimum total vertical
      permeability; this value truncates (censors) the vertical random
      distribution and limits lower values (default: 1.0e-18)

    * **maxPermeability** [|m^2|] (1.0e-21 to 1.0e-12) - maximum total vertical
      permeability; this value truncates (censors) the random distribution and
      limits higher values (default: 1.0e-15)

    Note: The setup of the four distribution parameters above is not yet implemented
    in the control file interface or GUI of NRAP-Open-IAM and available
    only in the script interface.

    * **heterFactor** [-] (1.0e-2 to 100) - increase factor of the permeability
      of cells selected for heterogeneity, if the heterogeneity approach is used
      (default: 0.5).

    Reference parameters for all cells:

    * **aveBaseDepth** [|m|] (800 to 9500) - average depth to base of
      cell/reservoir top; interpolation depth (default: 1100)

    * **aveBasePressure** [|Pa|] (1.0e+6 to 6.0e+7) - average pressure
      at seal base during injection (default: 3.3e+7)

    * **aveTemperature** [|C|] (31 to 180) - average temperature of seal
      (default: 50)

    * **salinity** [|ppm|] (0 to 80000) - average salinity of seal
      (default: 1.5e+4)

    * **staticDepth** [|m|] (800 to 9500) - reference depth for computing
      static pressure at top of seal (default: 1000)

    * **staticPressure** [|Pa|] (1.0e+6 to 6.0e+7) - pressure at static
      reference depth for computing pressure at the cell top (default: 1.0e+7).

    Fluid (conditions) parameters:

    * **brineDensity** [|kg/m^3|] (880 to 1080) - density of brine phase
      (default: 1004)

    * **CO2Density** [|kg/m^3|] (93 to 1050) - density of |CO2| phase
      (default:  597.8)

    * **brineViscosity** [|Pa*s|] (1.5e-4 to 1.6e-3) - viscosity of brine phase
      (default: 5.634e-4)

    * **CO2Viscosity** [|Pa*s|] (1.8e-5 to 1.4e-4) - viscosity of |CO2| phase
      (default: 4.452e-5)

    * **CO2Solubility** [|mol/kg|] (0 to 2) - solubility of |CO2| phase in brine
      (default: 0.035).

    Two-phase model parameters for LET model:

    * **wetting1** [-] (0.5 to 5) - wetting phase parameter |L| (default: 1)

    * **wetting2** [-] (0.1 to 30) - wetting phase parameter |E| (default: 10)

    * **wetting3** [-] (0 to 3) - wetting phase parameter |T| (default: 1.25)

    * **nonwet1** [-] (0.5 to 5) - nonwetting phase parameter |L| (default: 1.05)

    * **nonwet2** [-] (0.1 to 30) - nonwetting phase parameter |E| (default: 10)

    * **nonwet3** [-] (0 to 3) - nonwetting phase parameter |T| (default: 1.25)

    * **capillary1** [-] (0.01 to 5) - LET-model parameter |L| for capillary
      pressure (default: 0.2)

    * **capillary2** [-] (0.01 to 30) - LET-model parameter |E| for capillary
      pressure (default: 2.8)

    * **capillary3** [-] (0.01 to 3) - LET-model parameter |T| for capillary
      pressure (default: 0.43)

    * **maxCapillary** [|Pa|] (100 to 2.0e+8) - maximum capillary pressure
      for model (default: 1.0e+7)

    Note: Parameters **wetting1**, **wetting2**, **wetting3**, **nonwet1**,
    **nonwet2**, **nonwet3**, **capillary1**, **capillary2**, **capillary3**,
    and **maxCapillary** are used only if parameter **relativeModel** is set
    to *LET*.

    Parameters for BC model:

    * **lambda** [-] (0 to 5) - lambda parameter in Brooks-Corey model (default: 2.5)

    Note: Parameter **lambda** is used only if parameter/keyword argument
    **relativeModel** is set to *BC*.

    Additional parameters for two-phase flow:

    * **brineResSaturation** [-] (0.01 to 0.35) - residual brine saturation
      (default: 0.15)

    * **CO2ResSaturation** [-] (0 to 0.35) - residual |CO2| saturation
      (default: 0)

    * **relativeModel** [-] (LET or BC) - relative permeability model (default: LET)

    * **permRatio** [-] (0 to 1.5) - ratio of nonwetting to wetting permeability
      (default: 0.6).

    Time-model and rock type parameters:

    * **influenceModel** [-] (integer: 0, 1, 2) - time-dependent permeability model
      (default: 0); deterministic parameter, i.e. cannot be set to be random.
      Model type used to compute the influence factor of the fluid flow
      on permeability for time-dependent response:

      - 0: No influence factor used.

      - 1: Use a time-dependent model based on exposure time to |CO2|.
        Parameters **rateEffect** and **totalEffect** control the initial
        time delay and the maximum extent of effect.

      - 2: Use a multivariant model that considers **reactivity**, **clayType**,
        **clayContent** and **carbonateContent** values together with **rateEffect**
        and **totalEffect** parameters to establish the magnitude
        of the influence factor.

    * **influence** [-] (0 to 1) - initial permeability influence factor
      (default: 1)

    * **rateEffect** [-] (0.01 to 0.65) - time variance parameter; this
      parameter controls the initial time delay in the permeability effect
      of the model (default: 0.1)

    * **totalEffect** [-] (0.01 to 200) - time variance parameter; this
      parameter defines the total change in permeability of the model (as a factor)
      (default: 0.1)

    * **reactivity** [-] (0 to 10) - reactivity of time model; factor controls
      the magnitude of permeability change (default: 8)

    * **clayType** [-] (smectite, illite, or chlorite) - predominate
      clay mineral content in the seal horizon, defined as one of following categories:

      - smectite (high swelling material)

      - illite (moderate swelling material)

      - chlorite (low swelling material) (default: smectite)

    * **carbonateContent** [%] (0 to 100) - carbonate content in seal layer rock
      (default: 8)

    * **clayContent** [%] (0 to 100) - clay mineral content in seal layer rock
      (default: 60).

    Note: Parameters **rateEffect** and **totalEffect** are used only when parameter
    **influenceModel** is set to 1 or 2. These parameters control the initial time
    delay and the maximum extent of effect.

    Note: Parameters **reactivity**, **clayType**, **carbonateContent**, and
    **clayContent** are used only when parameter **influenceModel** is set to 2.

    The possible outputs from the Seal Horizon component are leakage rates
    of |CO2| and brine to aquifer through seal layer. The names
    of the observations are of the form:

    * **CO2_aquifer**, **brine_aquifer** [|kg/s|] - |CO2| and brine leakage rates to
      aquifer through seal layer (individual cells) into overlying aquifer

    * **mass_CO2_aquifer**, **mass_brine_aquifer** [|kg|] - mass of the |CO2|
      and brine leaked through seal layer (individual cells) into overlying aquifer

    * **CO2_aquifer_total**, **brine_aquifer_total** [|kg/s|] - cumulative
      (for all cells) |CO2| and brine leakage rates to aquifer through
      seal layer into overlying aquifer

    * **mass_CO2_aquifer_total**, **mass_brine_aquifer_total** [|kg|] - cumulative
      (for all cells) mass of the |CO2| and brine leaked through seal layer
      into overlying aquifer.
    """

    def __init__(self, name, parent, locX=None, locY=None, area=None, **kwargs):
        """
        Constructor method of SealHorizon class

        :param name: name of component model
        :type name: [str]

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel [object]

        :param locX: x-coordinates of each individual cell center; between
            -1.0e+5 and 1.0e+5
        :type locX: [float]

        :param locY: y-coordinates of each individual cell center; between
            -1.0e+5 and 1.0e+5
        :type locY: [float]

        :param area: area of each individual cell of interest;
            area can be a scalar if it is the same for all cells
        :type area: [float] or float

        :param kwargs: dictionary of additional arguments. Possible keys:
            'correlate_entry_approach', 'interpolate_approach',
            'heterogeneity_approach', 'vary_pressure_approach'
        :type kwargs: dict

        :returns: SealHorizon class object
        """
        model_kwargs = {'time_point': 365.25, 'time_step': 365.25}

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'SealHorizon'

        # Copy default setup
        self.seal_controls = {}
        self.seal_controls = sdeci.interpret_yaml_parameters(
            SEAL_SETUP_DICT, self.seal_controls)

        # Process keyword arguments related to controls: should be after reading
        # dictionary with seal setup
        self.process_additional_input(**kwargs)

        # Define cells grid
        self.grid = []

        # Define dictionary of input parameters and temporal input boundaries
        # Parameter bounds are to be defined first as the subsequent methods
        # call for the parameter bounds dictionary
        self.define_pars_bounds()

        # Check provided cell coordinates; if conditions are satisfied then save.
        if locX is not None and locY is not None:
            if len(locX) == len(locY):  # should have the same length
                # Create array with cell centers of shape (N,2)
                self.cell_xy_centers = np.array([locX, locY]).T
            else:
                msg = ''.join(['Length of array locX is not equal ',
                               'to the length of array locY.'])
                logging.error(msg)
                raise IndexError(msg)

            if area is None:
                msg = ''.join(
                    ['Argument area is not provided with cell centers locations.',
                     'Default value of 10000.0 will be used.'])
                logging.warning(msg)
                area = 10000.0

            # Create grid
            self.create_grid()

            # Check area
            self.check_area(area)

            # Update model_kwargs
            # Note: If area is provided to the constructor method then it will be
            # used to update cell parameters in the simulation_model method
            # If it's not provided then default value of 1.0e+4 is used.
            self.model_kwargs['area'] = self.area

            # Add accumulators
            self.add_mass_accumulators()

        else:
            # This is for the control file interface and GUI when we create
            # an instance object without knowing the locations of cells.
            # The locations will be assigned with the connect_with_system method
            self.cell_xy_centers = None
            self.area = None

        # Set default parameters
        self.setup_default_pars()

        # Setup potential names of individual cell parameters
        self.grid_pars_keys = ['area', 'baseDepth', 'influence', 'entryPressure',
                               'permeability', 'thickness']

        # Define gridded observations names
        self.grid_obs_keys = SH_GRID_OBSERVATIONS

    def process_additional_input(self, **kwargs):
        """
        Process additional input keyword arguments provided to the constructor method.

        Possible keys:
            'correlate_entry_approach' (default: False),
            'interpolate_approach' (default: False),
            'heterogeneity_approach' (default: False)

        Parameters
        ----------
        **kwargs : dict
            Dictionary containing additional input provided to the constructor.

        Returns
        -------
        None.
        """
        for key, val in SH_CFI_CONTROLS.items():
            self.seal_controls[key] = kwargs.get(key, val[1])

    def add_mass_accumulators(self):
        self.add_accumulator('mass_CO2_aquifer', sim=np.zeros(self.num_cells))
        self.add_accumulator('mass_brine_aquifer', sim=np.zeros(self.num_cells))

    def update_num_cells(self):
        self.num_cells = self.cell_xy_centers.shape[0]

    def create_grid(self):
        """
        If locations are provided, setup dependent attributes
        """
        # Determine and update number of cells
        self.update_num_cells()

        # Create grid of cells
        if not self.grid:
            for ind in range(self.num_cells):
                # Create a cell
                self.grid.append(
                    smod.Cell(x_center=self.cell_xy_centers[ind, 0],
                              y_center=self.cell_xy_centers[ind, 1]))

    def check_area(self, area):
        # Check whether the area argument is a scalar
        if isinstance(area, (int, float)):
            self.area = area*np.ones(self.num_cells)
        else:
            if len(area) == self.num_cells:
                self.area = np.array(area)
            else:
                msg = ''.join(['Length of area array is not equal ',
                               'to the number of cells.'])
                logging.error(msg)
                raise IndexError(msg)

        # Check whether area satisfies the model limits
        if np.any(self.area < self.pars_bounds['area'][0]) or (
                np.any(self.area > self.pars_bounds['area'][1])):
            msg = 'Cell areas do not satisfy model limits [{}, {}].'.format(
                self.pars_bounds['area'][0], self.pars_bounds['area'][1])
            logging.error(msg)
            raise ValueError(msg)

    def setup_default_pars(self):
        """
        Add parameters of the component with default values.
        """
        # HINT: Since parameters thickness and permeability might assume a particular
        # distribution defined by the standalone code there are no default
        # parameter setup for these parameters. We process them separately,
        # in a different way from the rest of parameters.
        # Parameter baseDepth is also not defined in the default pars
        # as it is assumed the default value of aveBaseDepth
        # Parameter entryPressure is also not defined in the default pars
        # as it is defined by user, or calculated based on correlation with permeability,
        # or assumes default value when not provided
        # These 3 parameters (thickness, permeability, and baseDepth) are not
        # present in the original setup of the standalone code. Here, in SH
        # component they are introduced for user convenience.
        # Description: all cells
        self.add_default_par(
            'area', value=SEAL_SETUP_DICT['Description']['area']['valu'])
        self.add_default_par(
            'salinity', value=SEAL_SETUP_DICT['Description']['salinity']['valu'])
        self.add_default_par(
            'aveTemperature',
            value=SEAL_SETUP_DICT['Description']['aveTemperature']['valu'])
        self.add_default_par(
            'aveBaseDepth',
            value=SEAL_SETUP_DICT['Description']['aveBaseDepth']['valu'])
        self.add_default_par(
            'aveBasePressure',
            value=SEAL_SETUP_DICT['Description']['aveBasePressure']['valu'])
        self.add_default_par(
            'staticDepth',
            value=SEAL_SETUP_DICT['Description']['staticDepth']['valu'])
        self.add_default_par(
            'staticPressure',
            value=SEAL_SETUP_DICT['Description']['staticPressure']['valu'])

        # Conditions parameters: all cells
        self.add_default_par(
            'CO2Density',
            value=SEAL_SETUP_DICT['Conditions']['CO2Density']['valu'])
        self.add_default_par(
            'CO2Viscosity',
            value=SEAL_SETUP_DICT['Conditions']['CO2Viscosity']['valu'])
        self.add_default_par(
            'brineDensity',
            value=SEAL_SETUP_DICT['Conditions']['brineDensity']['valu'])
        self.add_default_par(
            'brineViscosity',
            value=SEAL_SETUP_DICT['Conditions']['brineViscosity']['valu'])
        self.add_default_par(
            'CO2Solubility',
            value=SEAL_SETUP_DICT['Conditions']['CO2Solubility']['valu'])

        # Thickness: all cells
        self.add_default_par(
            'aveThickness',
            value=SEAL_SETUP_DICT['Thickness']['aveThickness']['valu'])
        self.add_default_par(
            'stdDevThickness',
            value=SEAL_SETUP_DICT['Thickness']['stdDevThickness']['valu'])
        self.add_default_par(
            'minThickness',
            value=SEAL_SETUP_DICT['Thickness']['minThickness']['valu'])
        self.add_default_par(
            'maxThickness',
            value=SEAL_SETUP_DICT['Thickness']['maxThickness']['valu'])

        # Permeability: all cells
        scalar = sunit.microd_to_metersq()
        self.add_default_par(
            'avePermeability',
            value=scalar*SEAL_SETUP_DICT['Permeability']['avePermeability']['valu'])
        self.add_default_par(
            'stdDevPermeability',
            value=scalar*SEAL_SETUP_DICT['Permeability']['stdDevPermeability']['valu'])
        self.add_default_par(
            'minPermeability',
            value=scalar*SEAL_SETUP_DICT['Permeability']['minPermeability']['valu'])
        self.add_default_par(
            'maxPermeability',
            value=scalar*SEAL_SETUP_DICT['Permeability']['maxPermeability']['valu'])
        self.add_default_par(
            'heterFactor',
            value=SEAL_SETUP_DICT['Permeability']['heterFactor'])

        # Relative flow: all cells
        self.add_default_par(
            'brineResSaturation',
            value=SEAL_SETUP_DICT['RelativeFlowLimits']['brineResSaturation'])
        self.add_default_par(
            'CO2ResSaturation',
            value=SEAL_SETUP_DICT['RelativeFlowLimits']['CO2ResSaturation'])
        self.add_default_par(
            'permRatio', value=SEAL_SETUP_DICT['RelativeFlowLimits']['permRatio'])

        # String variable cannot be parameter so the only way to pass
        # string type parameters is through keyword arguments to the model method
        # LET, default value
        self.model_kwargs['relativeModel'] = \
            SEAL_SETUP_DICT['RelativeFlowLimits']['relativeModel']

        # BC model: all cells
        self.add_default_par('lambda',
                             value=SEAL_SETUP_DICT['BrooksCoreyModel']['lambda'])

        # Two-phase model parameters for L-E-T model: all cells
        self.add_default_par('wetting1',
                             value=SEAL_SETUP_DICT['LETModel']['wetting1'])
        self.add_default_par('wetting2',
                             value=SEAL_SETUP_DICT['LETModel']['wetting2'])
        self.add_default_par('wetting3',
                             value=SEAL_SETUP_DICT['LETModel']['wetting3'])
        self.add_default_par('nonwet1',
                             value=SEAL_SETUP_DICT['LETModel']['nonwet1'])
        self.add_default_par('nonwet2',
                             value=SEAL_SETUP_DICT['LETModel']['nonwet2'])
        self.add_default_par('nonwet3',
                             value=SEAL_SETUP_DICT['LETModel']['nonwet3'])
        self.add_default_par(
            'capillary1', value=SEAL_SETUP_DICT['LETCapillaryModel']['capillary1'])
        self.add_default_par(
            'capillary2', value=SEAL_SETUP_DICT['LETCapillaryModel']['capillary2'])
        self.add_default_par(
            'capillary3', value=SEAL_SETUP_DICT['LETCapillaryModel']['capillary3'])
        self.add_default_par(
            'maxCapillary',
            value=SEAL_SETUP_DICT['LETCapillaryModel']['maxCapillary']['valu'])

        # Time-model parameters: all cells
        self.add_default_par('influenceModel',
                             value=SEAL_SETUP_DICT['TimeModel']['influenceModel'])
        self.add_default_par('influence',
                             value=SEAL_SETUP_DICT['TimeModel']['influence'])
        self.add_default_par('totalEffect',
                             value=SEAL_SETUP_DICT['TimeModel']['totalEffect'])
        self.add_default_par('rateEffect',
                             value=SEAL_SETUP_DICT['TimeModel']['rateEffect'])
        self.add_default_par('reactivity',
                             value=SEAL_SETUP_DICT['TimeModel']['reactivity'])
        self.add_default_par('carbonateContent',
                             value=SEAL_SETUP_DICT['TimeModel']['carbonateContent'])
        self.add_default_par('clayContent',
                             value=SEAL_SETUP_DICT['TimeModel']['clayContent'])
        # String variable cannot be parameter so the only way to pass
        # string type parameters is through keyword arguments to the model method
        # 'smectite' is a default value
        self.model_kwargs['clayType'] = SEAL_SETUP_DICT['TimeModel']['clayType']

    def define_pars_bounds(self):
        """ Define dictionaries for limits of mode parameters and temporal inputs."""
        param_bounds = sintro.define_input_limits()

        # Define dictionary of boundaries
        self.pars_bounds = dict()

        # Controls
        # We don't check (for now) the next 4 parameters as they relate to the model setup
        self.pars_bounds['startTime'] = param_bounds['start_time']
        self.pars_bounds['endTime'] = param_bounds['end_time']
        self.pars_bounds['timePoints'] = param_bounds['time_points']
        self.pars_bounds['realizations'] = param_bounds['realizations']

        # Grid
        self.pars_bounds['gridRows'] = param_bounds['grid_rows']
        self.pars_bounds['gridColumns'] = param_bounds['grid_cols']
        self.pars_bounds['cellHeight'] = param_bounds['cell_height']
        self.pars_bounds['cellWidth'] = param_bounds['cell_width']

        # Geometry
        self.pars_bounds['numCells'] = param_bounds['num_cells']
        self.pars_bounds['area'] = param_bounds['area']
        self.pars_bounds['baseDepth'] = param_bounds['ave_base_depth']
        self.pars_bounds['staticDepth'] = param_bounds['static_depth']
        self.pars_bounds['staticPressure'] = param_bounds['static_pressure']

        # Conditions: all cells
        self.pars_bounds['aveTemperature'] = param_bounds['temperature']
        self.pars_bounds['salinity'] = param_bounds['salinity']
        self.pars_bounds['aveBaseDepth'] = param_bounds['ave_base_depth']
        self.pars_bounds['aveBasePressure'] = param_bounds['ave_base_pressure']

        # Fluid parameters: all cells
        self.pars_bounds['brineDensity'] = param_bounds['brine_density']
        self.pars_bounds['CO2Density'] = param_bounds['co2_density']
        self.pars_bounds['brineViscosity'] = param_bounds['brine_viscosity']
        self.pars_bounds['CO2Viscosity'] = param_bounds['co2_viscosity']
        self.pars_bounds['CO2Solubility'] = param_bounds['co2_solubility']

        # Permeability: all cells
        scalar = sunit.microd_to_metersq()
        self.pars_bounds['permeability'] = [
            scalar*param_bounds['perm_mean'][0], scalar*param_bounds['perm_mean'][1]]
        self.pars_bounds['avePermeability'] = [
            scalar*param_bounds['perm_mean'][0], scalar*param_bounds['perm_mean'][1]]
        self.pars_bounds['stdDevPermeability'] = [
            scalar*param_bounds['perm_std'][0], scalar*param_bounds['perm_std'][1]]
        self.pars_bounds['minPermeability'] = [
            scalar*param_bounds['perm_min'][0], scalar*param_bounds['perm_min'][1]]
        self.pars_bounds['maxPermeability'] = [
            scalar*param_bounds['perm_max'][0], scalar*param_bounds['perm_max'][1]]
        self.pars_bounds['heterFactor'] = param_bounds['perm_heter_factor']

        # Thickness: all cells
        self.pars_bounds['thickness'] = param_bounds['thickness_ave']
        self.pars_bounds['aveThickness'] = param_bounds['thickness_ave']
        self.pars_bounds['stdDevThickness'] = param_bounds['thickness_std']
        self.pars_bounds['minThickness'] = param_bounds['thickness_min']
        self.pars_bounds['maxThickness'] = param_bounds['thickness_max']

        # Relative flow: all cells
        self.pars_bounds['brineResSaturation'] = param_bounds['resid_brine']
        self.pars_bounds['CO2ResSaturation'] = param_bounds['resid_co2']
        self.pars_bounds['permRatio'] = param_bounds['perm_ratio']
        self.pars_bounds['relativeModel'] = ['BC', 'LET']
        self.pars_bounds['entryPressure'] = param_bounds['entry_pressure']

        # Two-phase model parameters for L-E-T model: all cells
        self.pars_bounds['wetting1'] = param_bounds['l_wetting']
        self.pars_bounds['wetting2'] = param_bounds['e_wetting']
        self.pars_bounds['wetting3'] = param_bounds['t_wetting']
        self.pars_bounds['nonwet1'] = param_bounds['l_nonwet']
        self.pars_bounds['nonwet2'] = param_bounds['e_nonwet']
        self.pars_bounds['nonwet3'] = param_bounds['t_nonwet']
        self.pars_bounds['capillary1'] = param_bounds['l_capillary']
        self.pars_bounds['capillary2'] = param_bounds['e_capillary']
        self.pars_bounds['capillary3'] = param_bounds['t_capillary']
        self.pars_bounds['maxCapillary'] = param_bounds['max_capillary']

        # BC model: all cells
        self.pars_bounds['lambda'] = param_bounds['zeta']

        # Time-model parameters:all cells
        self.pars_bounds['influenceModel'] = [0, 1, 2]
        self.pars_bounds['influence'] = param_bounds['influence']
        self.pars_bounds['totalEffect'] = param_bounds['total_effect']
        self.pars_bounds['rateEffect'] = param_bounds['rate_effect']
        self.pars_bounds['reactivity'] = param_bounds['reactivity']
        self.pars_bounds['clayContent'] = param_bounds['clay_content']
        self.pars_bounds['carbonateContent'] = param_bounds['carbonate_content']
        self.pars_bounds['clayType'] = ["smectite", "illite", "chlorite"]

        # Define dictionary of temporal data limits
        self.temp_data_bounds = dict()
        self.temp_data_bounds['pressure'] = ['Pressure']+param_bounds['reservoir_pressure']
        self.temp_data_bounds['top_brine_pressure'] = ['Top brine pressure']+param_bounds['reservoir_pressure']
        self.temp_data_bounds['CO2saturation'] = (
            ['CO2 saturation']+param_bounds['reservoir_saturation'])

    def check_distribution_data(self, base_par_name, mean_val, min_val, max_val):
        """
        Provide additional check on the thickness and permeability
        distribution parameters.

        Code is partially based on parameters checks in the method
        seal_intro.cross_check_input of the standalone version of the code.
        """
        # Compare min and max of the distribution
        if min_val > max_val:
            err_msg = 'Minimum {par} exceeds maximum {par}.'.format(par=base_par_name)
            logging.error(err_msg)
            raise ValueError(err_msg)

        # Check average thickness against bounds.
        condition_a = mean_val < min_val
        condition_b = mean_val > max_val
        if condition_a or condition_b:
            err_msg = 'Mean {par} is outside min/max bounds.'.format(par=base_par_name)
            logging.error(err_msg)
            raise ValueError(err_msg)

    def update_seal_controls(self, p):
        """ Update seal_controls with values provided by user."""

        # Model parameters
        if self._parent.time_array is not None:
            self.seal_controls['start_time'] = self._parent.time_array[0] # first time point
            self.seal_controls['end_time'] = self._parent.time_array[-1]  # last time point
            self.seal_controls['time_points'] = len(self._parent.time_array)
        else:
            err_msg = 'Argument time_array has to be defined for system model.'
            logging.error(err_msg)

        # Description section parameters
        self.seal_controls['num_cells'] = self.num_cells
        self.seal_controls['area'] = p['area']
        self.seal_controls['salinity'] = p['salinity']
        self.seal_controls['temperature'] = p['aveTemperature']
        self.seal_controls['static_depth'] = p['staticDepth']
        self.seal_controls['static_pressure'] = p['staticPressure']
        self.seal_controls['ave_base_depth'] = p['aveBaseDepth']
        self.seal_controls['ave_base_pressure'] = p['aveBasePressure']
        p['depth_data'] = p.get(
            'baseDepth', self.seal_controls['ave_base_depth']*np.ones(self.num_cells))

        # Conditions section parameters
        self.seal_controls['co2_density'] = p['CO2Density']
        self.seal_controls['co2_viscosity'] = p['CO2Viscosity']
        self.seal_controls['brine_density'] = p['brineDensity']
        self.seal_controls['brine_viscosity'] = p['brineViscosity']
        self.seal_controls['co2_solubility'] = p['CO2Solubility']

        # Thickness parameter(s)
        if 'thickness' in p:   # user provided data
            p['thickness_data'] = p['thickness']  # array
            self.seal_controls['thickness_approach'] = True
        else: # random from distribution
            self.seal_controls['thickness_ave'] = p['aveThickness']
            self.seal_controls['thickness_std'] = p['stdDevThickness']
            self.seal_controls['thickness_min'] = p['minThickness']
            self.seal_controls['thickness_max'] = p['maxThickness']
            p['thickness_data'] = None
            self.seal_controls['thickness_approach'] = False
            self.check_distribution_data('thickness',
                                         self.seal_controls['thickness_ave'],
                                         self.seal_controls['thickness_min'],
                                         self.seal_controls['thickness_max'])

        # Permeability parameter(s)
        self.seal_controls['perm_heter_factor'] = p['heterFactor']
        scalar = sunit.metersq_to_microd()
        if 'permeability' in p:   # user provided data
            p['perm_data'] = scalar*p['permeability']
            self.seal_controls['perm_mean'] = np.mean(p['perm_data']) # might need for entry pressure
            self.seal_controls['heterogeneity_approach'] = False # not considered for user data
            self.seal_controls['vary_perm_choice'] = False       # variability is unknown
            self.seal_controls['perm_input_approach'] = True     # user provided
            # Define log-normal parameters
            self.seal_controls['perm_location'] = 0.0  # might not need these two
            self.seal_controls['perm_scale'] = 0.0 # keeping for consistency with standalone code
        else:  # random from distribution
            self.seal_controls['perm_mean'] = scalar*p['avePermeability']
            self.seal_controls['perm_std'] = scalar*p['stdDevPermeability']
            if self.seal_controls['perm_std'] > sdeci.MINIMUM_PERM:
                self.seal_controls['vary_perm_choice'] = True
            else:
                self.seal_controls['vary_perm_choice'] = False
            self.seal_controls['perm_min'] = scalar*p['minPermeability']
            self.seal_controls['perm_max'] = scalar*p['maxPermeability']
            p['perm_data'] = None
            self.seal_controls['perm_input_approach'] = False
            self.check_distribution_data('permeability',
                                         self.seal_controls['perm_mean'],
                                         self.seal_controls['perm_min'],
                                         self.seal_controls['perm_max'])
            # Define log-normal parameters
            self.seal_controls['perm_location'], self.seal_controls['perm_scale'] = \
                sintro.convert_lognorm_terms(self.seal_controls['perm_mean'],
                                            self.seal_controls['perm_std'])

        # Check type of model for relative permeability
        if p['relativeModel'] == 'BC':
            # BrooksCoreyModel:
            self.seal_controls['zeta'] = p['lambda']
        else:
            # LETModel parameters
            self.seal_controls['l_wetting'] = p['wetting1']
            self.seal_controls['e_wetting'] = p['wetting2']
            self.seal_controls['t_wetting'] = p['wetting3']
            self.seal_controls['l_nonwet'] = p['nonwet1']
            self.seal_controls['e_nonwet'] = p['nonwet2']
            self.seal_controls['t_nonwet'] = p['nonwet3']

            # LET CapillaryPressure:
            self.seal_controls['l_capillary'] = p['capillary1']
            self.seal_controls['e_capillary'] = p['capillary2']
            self.seal_controls['t_capillary'] = p['capillary3']
            self.seal_controls['max_capillary'] = p['maxCapillary']

        # Additional parameters for two-phase flow
        self.seal_controls['resid_brine'] = p['brineResSaturation']
        self.seal_controls['resid_co2'] = p['CO2ResSaturation']
        self.seal_controls['relative_model'] = p['relativeModel']
        self.seal_controls['perm_ratio'] = p['permRatio']
        if 'entryPressure' in p:  # user provided
            p['entry_pressure_data'] = p['entryPressure']
            self.seal_controls['entry_approach'] = True
        else:  # will be defined based on correlation or default value
            p['entry_pressure_data'] = self.seal_controls['entry_pressure']*np.ones(
                self.num_cells)
            self.seal_controls['entry_approach'] = False

        # Time parameters
        self.seal_controls['model'] = p['influenceModel']
        self.seal_controls['influence'] = p['influence']
        self.seal_controls['rate_effect'] = p['rateEffect']
        self.seal_controls['total_effect'] = p['totalEffect']
        self.seal_controls['reactivity'] = p['reactivity']
        self.seal_controls['clay_content'] = p['clayContent']
        self.seal_controls['carbonate_content'] = p['carbonateContent']
        self.seal_controls['clay_type'] = p['clayType']

        return p

    def update_all_parameters(self, scalar_p, **kwarg_p):
        """
        Transform all model parameters into array versions of themselves.
        """
        # Update major arrays
        combined_p = {}
        for key in self.grid_pars_keys:
            if key in kwarg_p:
                combined_p[key] = np.array(kwarg_p[key])
            elif key in scalar_p:
                combined_p[key] = scalar_p[key]*np.ones(self.num_cells)

        for key in ['clayType', 'relativeModel']:
            combined_p[key] = kwarg_p[key]

        for key in scalar_p:
            if key not in self.grid_pars_keys:
                combined_p[key] = scalar_p[key]

        return combined_p

    def update_cell_strata(self, p):
        """
        Update thickness and depth to the top of the cells.
        """
        # Loop over each cell.
        for ind in range(self.num_cells):
            # user provided or randomly generated
            self.grid[ind].thickness = p['thickness_data'][ind]
            # Depth to top is (depth to base - thickness)
            self.grid[ind].set_top_depth(p['depth_data'][ind])

    def update_cell_attributes(self, p):
        """
        Update parameters of all cells in the grid.
        """
        # Loop over each cell.
        for ind in range(self.num_cells):
            self.grid[ind].area = p['area'][ind]
            # user provided or randomly generated
            self.grid[ind].permeability = p['perm_data'][ind]
            self.grid[ind].influence = p['influence'][ind]
            # Default or provided by user
            self.grid[ind].entry = p['entry_pressure_data'][ind]

        # Update attributes common to all cells
        smod.Cell.assign_controls(self.seal_controls)

    def check_temporal_inputs(self, time, temp_inputs):
        """
        Check whether temporal data fall within specified boundaries.

        :param time: time point at which temporal input is compared with the
            minimum and maximum allowed values.
        :type time: float

        :param temp_inputs: temporal input data of component model
        :type temp_inputs: dict
        """
        debug_msg = 'Temporal inputs of component {} are {}'.format(
            self.name, temp_inputs)
        logging.debug(debug_msg)

        for key, val in temp_inputs.items():
            if val is not None:
                arr = np.array(val)
                if (np.any(arr < self.temp_data_bounds[key][1])) or (
                        np.any(arr > self.temp_data_bounds[key][2])):
                    err_msg = ''.join([
                        'Temporal input {} of SealHorizon component {} ',
                        'is outside the model range [{}, {}] at time t = {} days']).format(
                            self.temp_data_bounds[key][0].lower(), self.name,
                            self.temp_data_bounds[key][1],
                            self.temp_data_bounds[key][2], time)
                    logging.error(err_msg)
                    raise ValueError('Temporal inputs are outside the ROM limits.')

    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        if self.cell_xy_centers is None:
            msg = ''.join(['x- and y-coordinates of cell centers are not '
                           'provided.'])
            logging.error(msg)
            raise ValueError(msg)

        msg = ('Input parameters of {name} component are {p}'.
               format(name=self.name, p=p))
        logging.debug(msg)

        for key, val in p.items():
            if key == 'influenceModel':
                # self.pars_bounds[key] is list of possible values [0, 1, 2]
                if val not in self.pars_bounds['influenceModel']:
                    msg = 'Invalid value of influenceModel parameter: {}.'.format(val)
                    logging.error(msg)
                    raise ValueError(msg)
                continue

            if key not in self.pars_bounds:
                msg = ('Parameter {key} not recognized as a SealHorizon input ' +
                       'parameter.')
                msg = msg.format(key=key)
                logging.warning(msg)
            else:
                if ((val < self.pars_bounds[key][0]) or
                        (val > self.pars_bounds[key][1])):
                    msg = 'Parameter {} is out of bounds.'.format(key)
                    logging.warning(msg)

    def check_keyword_parameters(self, **kw_pars):
        """
        Check cell parameters provided as keyword arguments.

        :param kw_pars: input parameters of component model
        :type :
        """
        # Check the values of the gridded parameters for initial time point
        for key in self.grid_pars_keys:
            if key == 'area' and self.area is not None:
                continue
            else:
                if key in kw_pars:
                    if len(kw_pars[key]) == self.num_cells:

                        # Check that parameters satisfy boundaries
                        par_array = np.array(kw_pars[key])
                        if ((np.any(par_array < self.pars_bounds[key][0])) or
                                (np.any(par_array > self.pars_bounds[key][1]))):
                            msg = ('For some cells parameter ' +
                                   '{} is out of boundaries.'.format(key))
                            logging.warning(msg)
                    else:
                        msg = ('Length of keyword argument {} '.format(key) +
                               'is not equal to the number of cells.')
                        logging.error(msg)
                        raise IndexError(msg)

        # Check for (string) "clayType"
        if 'clayType' in kw_pars:
            if kw_pars['clayType'] not in ['smectite', 'illite', 'chlorite']:
                msg = "Argument clayType has an unrecognized value: {}.".format(
                    kw_pars['clayType'])
                logging.warning(msg)

                msg = ("Argument clayType will be set to 'smectite'.")
                logging.warning(msg)
                kw_pars['clayType'] = 'smectite'
        else:
            kw_pars['clayType'] = 'smectite'

        if 'relativeModel' in kw_pars:
            if kw_pars['relativeModel'] not in ['LET', 'BC']:
                msg = "Argument relativeModel has an unrecognized value: {}.".format(
                    kw_pars['relativeModel'])
                logging.warning(msg)

                msg = ("Argument relativeModel will be set to 'LET'.")
                logging.warning(msg)
                kw_pars['relativeModel'] = 'LET'
        else:
            kw_pars['relativeModel'] = 'LET'

        return kw_pars

    def manage_interpolation(self, top_pressure, base_pressure):
        """Manage interpolation of fluid properties for brine & CO2.

        Based on method manage_interpolation from seal_fluids module on the
        standalone code.
        """

        base_pressure_ave = np.average(base_pressure)
        top_pressure_ave = np.average(top_pressure)

        # Compute average pressure at center of seal for interpolation
        ave_pressure = (base_pressure_ave + top_pressure_ave)/2.0

        # Interpolate data at average pressure of seal center and reset values.
        self.seal_controls = \
            sfluid.interpolate_fluid_properties(self.seal_controls, ave_pressure)

    def obtain_permeability(self, rng):
        """
        Generate permeability using provided distribution parameters.
        """
        # Define permeability value for every cell in grid.
        if self.seal_controls['vary_perm_choice']:  # high enough value of std deviation parameter
            # Initialize permeability array
            perm_array = np.zeros(self.num_cells)

            # Get distribution parameters
            mid = self.seal_controls['perm_location']
            scale = self.seal_controls['perm_scale']

            for ind in range(self.num_cells):
                # For variable permeability, evaluate random value.
                perm_array[ind] = \
                    fran.evaluate_lognorm(mid, scale,
                                          self.seal_controls['perm_min'],
                                          self.seal_controls['perm_max'],
                                          rng)
        else:
            # For uniform permeability, set to mean
            perm_array = self.seal_controls['perm_mean']*np.ones(self.num_cells)

        return perm_array

    def calculate_static_pressure(self):
        """
        Calculate static pressure at top of each cell for cases when user does
        not provide this data.

        """
        # Initiate pressure array
        top_pressure = np.ones(self.num_cells)

        # Loop over all cells in the grid
        for ind, cell in enumerate(self.grid):
            top_pressure[ind] = \
                cell.compute_static_pressure(self.seal_controls['static_pressure'],
                                             self.seal_controls['static_depth'])
        return top_pressure

    def connect_with_system(self, component_data, name2obj_dict,
                            locations, system_adapters, **kwargs):
        """
        Code to add SealHorizon to system model for control file interface.

        :param component_data: Dictionary of component data including 'Connection'
        :type component_data: dict

        :returns: None
        """
        # Process controls
        if 'Controls' in component_data:
            controls_data = component_data['Controls']
            for key, val in SH_CFI_CONTROLS.items():
                controls_data[key] = controls_data.pop(val[0], val[1])
            self.process_additional_input(**controls_data)

        # Check cells setup data
        cell_data = component_data['Cells']

        # Get cell centers coordinates
        locX = locations[self.name]['coordx']
        locY = locations[self.name]['coordy']

        # Setup cell centers
        self.cell_xy_centers = np.array([locX, locY]).T

        # Create grid
        self.create_grid()

        # Read other parameters
        default_vals = {'area': 10000.0, 'thickness': 100.0, 'baseDepth': 1100.0,
                        'permeability': 2.5e-4, 'entryPressure': 5000.0,
                        'influence': 1.0}
        for key in default_vals:
            if key in cell_data:
                if isinstance(cell_data[key], str):
                    filename = os.path.join(IAM_DIR, cell_data[key])
                    # Using somewhat elaborate way to extract data from txt data file
                    # to make sure all sizes of data works, including size of 1
                    # Data file looks like a column, similar to input for
                    # locations in LUT
                    data = np.asarray([np.genfromtxt(filename)]).flatten()
                    default_vals[key] = data
                else:
                    default_vals[key] = cell_data[key]
                # Check what type of data is provided
                if isinstance(default_vals[key], (int, float)):
                    # Add as parameter if scalar is provided
                    self.add_par(key, value=default_vals[key], vary=False)
                else:
                    if len(default_vals[key]) == self.num_cells:
                        self.model_kwargs[key] = default_vals[key]
                    else:
                        msg = ''.join(['Length of {} array is not equal ',
                                       'to the number of cells.']).format(key)
                        logging.error(msg)
                        raise IndexError(msg)

        # Add accumulators
        self.add_mass_accumulators()

        # Add parameters
        if ('Parameters' in component_data) and (component_data['Parameters']):
            # Check for presence of keyword type parameters
            for kwarg_nm in ['relativeModel', 'clayType']:
                if kwarg_nm in component_data['Parameters']:
                    kw_data = component_data['Parameters'][kwarg_nm]
                    # Get the value and remove it from the dictionary
                    if not isinstance(kw_data, dict):
                        self.model_kwargs[kwarg_nm] = component_data[
                            'Parameters'].pop(kwarg_nm)
                    else:
                        self.model_kwargs[kwarg_nm] = component_data[
                            'Parameters'].pop(kwarg_nm)['value']

            for key in component_data['Parameters']:
                if not isinstance(component_data['Parameters'][key], dict):
                    component_data['Parameters'][key] = {
                        'value': component_data['Parameters'][key],
                        'vary': False}
                self.add_par(key, **component_data['Parameters'][key])

        # Check whether dynamic kwargs are provided
        process_dynamic_inputs(
            self, component_data, array_like=True, check_second_dim=True,
            size=len(locX), quantity_to_compare='Number of cells')

        if 'Connection' in component_data:
            # collectors is a dictionary containing inputs collected into arrays
            collectors = self._parent.collectors
            # Get references to the needed collectors
            cdict1 = collectors['pressure'][self.name]
            cdict2 = collectors['CO2saturation'][self.name]

            connections = np.array([component_data['Connection']]).flatten()
            for connect_nm in connections:
                try:
                    connection = name2obj_dict[connect_nm]
                except KeyError:
                    err_msg = ''.join([
                        'Component {} is not setup in the yaml ',
                        'dict data']).format(connect_nm)
                    logging.error(err_msg)
                    raise KeyError(err_msg)

                for obs_nm in ['pressure', 'CO2saturation']:
                    # Add to be linked observation of the reservoir component
                    connection.add_obs_to_be_linked(obs_nm)

                # Update list of links to the observations
                cdict1['data'].append(connection.linkobs['pressure'])
                cdict2['data'].append(connection.linkobs['CO2saturation'])

            # Add link to observations to the data of collectors
            self.add_kwarg_linked_to_collection('pressure', cdict1['data'])
            self.add_kwarg_linked_to_collection('CO2saturation', cdict2['data'])

        if 'Outputs' in component_data:
            comp_outputs = component_data['Outputs']
            new_comp_outputs = []
            for obs_nm in comp_outputs:
                if obs_nm in SH_SCALAR_OBSERVATIONS:
                    new_comp_outputs.append(obs_nm)
                    self.add_obs(obs_nm)
                elif obs_nm in SH_GRID_OBSERVATIONS:
                    # Get number of cells
                    num_locs = len(locX)
                    for ind in range(num_locs):
                        augm_obs_nm = obs_nm+'_cell_{}'.format(ind+1)
                        self.add_local_obs(augm_obs_nm, obs_nm, 'array', ind)
                        new_comp_outputs.append(augm_obs_nm)
                    # Add gridded observations
                    self.add_grid_obs(
                        obs_nm, constr_type='array', output_dir=kwargs['output_dir'])
                else:
                    warn_msg = ''.join([
                        '{} is not recognised as observation name ',
                        'of Seal Horizon component {}.']).format(obs_nm, self.name)
                    logging.warning(warn_msg)

            component_data['Outputs'] = new_comp_outputs

    def simulation_model(self, p, time_point=365.25, time_step=365.25,
                         pressure=None, CO2saturation=None,
                         top_brine_pressure=None, **kwargs):
        """
        :param p: input parameters of SealHorizon component model
        :type p: dict

        param time_point: time point in days at which the leakage rates are
            to be calculated; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param time_step: difference between the current and previous
            time points in days; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param pressure: pressure at the bottom of cells in the reservoir, in Pa;
            by default, its value is None
        :type pressure: [float]

        :param CO2saturation: saturation of |CO2| phase at the bottom
            of cells in the reservoir; by default, its value is None
        :type CO2saturation: [float]

        :param kwargs: additional keyword arguments of the method. Possible keys
            are 'depth', 'thickness', 'permeability', 'influence', 'entryPressure',
            'clayType', and 'relativeModel'. Values for keys 'depth', 'thickness',
            'permeability', 'influence', 'entryPressure' are array-like of size
            equal to the number of cells. Values for key 'clayType' are of
            type string with possible values 'smectite', 'illite', or 'chlorite'.
            Values for key 'relativeModel' are of type string with possible values
            'BC' or 'LET'.
        :type kwargs: dict

        :Returns: out - dictionary of observations of Seal Horizon model;
        outputs with keys ['CO2_aquifer', 'brine_aquifer',
        'mass_CO2_aquifer', 'mass_brine_aquifer'] are arrays;
        outputs with keys ['CO2_aquifer_total', 'brine_aquifer_total',
        'mass_CO2_aquifer_total', 'mass_brine_aquifer_total'] are scalars.
        """
        # Check whether seed parameter is given
        if 'seed' in p:
            # Set seed for random numbers on each simulation
            rng = np.random.default_rng(p['seed'])
        else:
            rng = np.random.default_rng(sconf.SEEDX)

        # Obtain the default values of the parameters from dictionary of
        # default parameters.
        actual_p = {k: v.value for k, v in self.default_pars.items()}

        # Update default values of parameters with the provided ones
        actual_p.update(p)

        self.check_temporal_inputs(time_point,
                                   {'pressure': pressure,
                                    'CO2saturation': CO2saturation,
                                    'top_brine_pressure': top_brine_pressure})

        # Check parameters provided by keyword arguments for initial time.
        msg = 'Provided keyword arguments: {}'.format(kwargs)
        kwarg_p = copy.deepcopy(kwargs)
        if time_point == 0:
            # Check whether top pressure is provided
            if top_brine_pressure is None:
                self.seal_controls['upper_approach'] = False
                self.seal_controls['vary_pressure_approach'] = False
            else:
                self.seal_controls['upper_approach'] = True
                # self.seal_controls['vary_pressure_approach'] can be True or False
                # If it's True new value will be provided for each time step (theoretically).
                # If it's False there will be a single value only for the time point 0.

            # Check parameters provided as arrays or strings
            kwarg_p = self.check_keyword_parameters(**kwargs)

            # Check whether any of the potentially array-type parameters
            # (e.g., thickness or permeability) are provided as scalars
            combined_p = self.update_all_parameters(actual_p, **kwarg_p)

            # Update seal_controls attribute with values provided by user
            # update_seal_controls also takes care of thickness, permeability
            # and baseDepth if they are provided by user
            # If thickness or permeability are not provided they will be generated
            # If baseDepth is not provided it assume default value equal
            # to default of aveBaseDepth
            combined_p = self.update_seal_controls(combined_p)

            # Check if thickness is to be generated
            if combined_p['thickness_data'] is None:
                # Generate thickness using distribution parameters
                combined_p['thickness_data'] = sref.thickness_variability(
                    self.seal_controls, rng)

            # Update thickness and depth to the top of cells
            self.update_cell_strata(combined_p)

            # Check if permeability is to be generated
            if combined_p['perm_data'] is None:
                # Generate permeability using distribution parameters
                combined_p['perm_data'] = self.obtain_permeability(rng)

            # Update the rest of cell attributes
            self.update_cell_attributes(combined_p)

            # If pressure at the top of seal is provided
            if self.seal_controls['upper_approach']:
                if len(top_brine_pressure) == self.num_cells:
                    self.press_top = top_brine_pressure
                else:
                    err_msg = ''.join([
                        'Length of argument top_brine_pressure does not equal ',
                        'to the number of cells defined for component {}.']).format(
                            self.name)
                    logging.error(err_msg)
                    raise IndexError(err_msg)

            else:
                # This relies on the depth to top of cell and brine density to be already defined
                self.press_top = self.calculate_static_pressure()

            # Adjust entry pressure if desired and assign to cells
            self.grid = perm.correlate_entry_pressure(self.seal_controls, self.grid)

            # Generate new fluid properties values based on interpolation as desired
            # Update cells and cell model class properties
            if self.seal_controls['interpolate_approach']:
                self.manage_interpolation(self.press_top, pressure)

            # Define permeability heterogeneities, if desired.
            self.grid = perm.evaluate_areal_heter(self.grid, self.seal_controls)


        # Initialize outputs.
        out = {}
        out['CO2_aquifer'] = np.zeros(self.num_cells)
        out['brine_aquifer'] = np.zeros(self.num_cells)
        out['mass_CO2_aquifer'] = np.zeros(self.num_cells)
        out['mass_brine_aquifer'] = np.zeros(self.num_cells)
        out['CO2_aquifer_total'] = 0.0
        out['brine_aquifer_total'] = 0.0
        out['mass_CO2_aquifer_total'] = 0.0
        out['mass_brine_aquifer_total'] = 0.0

        if time_point == 0.0:
            return out

        # Check if we don't have top pressure data
        if top_brine_pressure is None or not self.seal_controls['vary_pressure_approach']:
            top_brine_pressure = self.press_top  # can come from time point 0 if provided only at t=0

        # Go over all cells
        for ind in range(self.num_cells):

            # Compute flow rates for each cell - vol./sec.
            flow_rate_co2, flow_rate_brine = self.grid[ind].compute_flow(
                pressure[ind], CO2saturation[ind],
                top_brine_pressure[ind], self.seal_controls)

            if ind == 0:
                msg = ''.join([
                    'Volume flow rates:\n',
                    'CO2 rates: {}\n'.format(flow_rate_co2),
                    'Brine rates: {}\n'.format(flow_rate_brine)])
                logging.debug(msg)

            # Update history of cell
            self.grid[ind].track_history(pressure[ind], flow_rate_co2, time_point/365.25)

            # Update influence factors for cell
            self.grid[ind].compute_model_influence(flow_rate_co2)

            # Convert flows from m^3/sec to kg/sec
            out['CO2_aquifer'][ind] = flow_rate_co2 * self.grid[ind].co2Density
            out['brine_aquifer'][ind] = flow_rate_brine * self.grid[ind].brineDensity

        out['mass_CO2_aquifer'] = self.accumulators['mass_CO2_aquifer'].sim+\
            time_step*days_to_seconds()*out['CO2_aquifer']
        out['mass_brine_aquifer'] = self.accumulators['mass_brine_aquifer'].sim+\
            time_step*days_to_seconds()*out['brine_aquifer']

        self.accumulators['mass_CO2_aquifer'].sim = out['mass_CO2_aquifer']
        self.accumulators['mass_brine_aquifer'].sim = out['mass_brine_aquifer']

        # Calculate cumulative flow rate for all cells
        out['CO2_aquifer_total'] = np.sum(out['CO2_aquifer'])
        out['brine_aquifer_total'] = np.sum(out['brine_aquifer'])
        out['mass_CO2_aquifer_total'] = np.sum(out['mass_CO2_aquifer'])
        out['mass_brine_aquifer_total'] = np.sum(out['mass_brine_aquifer'])

        msg = ''.join(['Outputs at time{}:\n CO2 rate: {}\n Brine rate: {}\n',
                       'Total CO2 rate: {}\n Total brine rate: {}']).format(
                           time_point, out['CO2_aquifer'], out['brine_aquifer'],
                           out['CO2_aquifer_total'], out['brine_aquifer_total'])
        logging.debug(msg)

        return out


def days_to_seconds():
    """ Convert 1 day to seconds.
    """
    value = 86400.0  # in seconds
    return value


def test_scenario1():
    """ Run script example illustrating work of Seal Horizon component.

    Run script example illustrating work of Seal Horizon component with its input
    linked to the output of the Lookup Table Reservoir component.
    """
    # Setup location of lookup table data set
    file_directory = os.sep.join(['..', 'components', 'reservoir',
                                  'lookuptables', 'Kimb_54_sims'])

    if not os.path.exists(os.sep.join([file_directory, 'Reservoir_data_sim01.csv'])):
        msg = ''.join(['\nKimberlina data set can be downloaded ',
                       'from one of the following places:\n',
                       '1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam \n',
                       '2. https://gitlab.com/NRAP/Kimberlina_data  \n',
                       'Check this example description for more information.'])
        logging.error(msg)

    # Define keyword arguments of the system model.
    num_years = 5
    time_array = 365.25 * np.arange(0.0, num_years+1)
    seal_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model.
    sm = SystemModel(model_kwargs=seal_model_kwargs)

    # Add reservoir component.
    loc_X = [35315.8, 36036.4, 36757.1, 37477.7, 38198.4, 38919, 39639.7,
              40360.3, 41081, 41801.6, 42522.3, 43242.9, 43963.6, 44684.2, 45404.8]

    loc_Y = 15*[51110.2]

    ltres = sm.add_component_model_object(
        LookupTableReservoir(name='ltres', parent=sm, locX=loc_X, locY=loc_Y))

    ltres.build_and_link_interpolators(file_directory=file_directory,
                                       intr_family='reservoir',
                                       default_values={'salinity': 0.1,
                                                       'temperature': 50.0},
                                       recompute_triangulation=False,
                                       build_on_the_fly=False)

    ltres.add_par('index', value=5, vary=False)

    # Add observations of reservoir component model to be used as input for
    # SealHorizon component
    ltres.add_obs_to_be_linked('pressure', obs_type='grid',
                               constr_type='array')
    ltres.add_obs_to_be_linked('CO2saturation', obs_type='grid',
                               constr_type='array')
    # Add local observation: pressure at one of the cells
    ltres.add_local_obs('loc_pressure', grid_obs_name='pressure',
                        constr_type='array', loc_ind=3)

    cell_coord_x = loc_X
    cell_coord_y = loc_Y

    shc = sm.add_component_model_object(
        SealHorizon(name='shc', parent=sm,
                    locX=cell_coord_x, locY=cell_coord_y, area=6.3e+4))

    shc.add_par('thickness', value=100.0, vary=False)
    shc.add_par('permeability', value=2.5e-4*9.86923266716E-19, vary=False)
    shc.add_par('baseDepth', 1100.0, vary=False)

    # Add time varying keyword arguments of SealHorizon component
    shc.add_kwarg_linked_to_obs('pressure', ltres.linkobs['pressure'],
                                obs_type='grid', constr_type='array')
    shc.add_kwarg_linked_to_obs('CO2saturation',
                                ltres.linkobs['CO2saturation'],
                                obs_type='grid', constr_type='array')
    # Add local observation of the component
    shc.add_local_obs('loc_CO2_aquifer', grid_obs_name='CO2_aquifer',
                      constr_type='array', loc_ind=3)
    # Add scalar observations
    shc.add_obs('CO2_aquifer_total')
    shc.add_obs('brine_aquifer_total')

    print('Starting simulation...')
    sm.forward()
    print('Simulation is finished.')

    # Collect observations
    pressure = sm.collect_observations_as_time_series(ltres, 'loc_pressure')
    print('Pressure at location 3:', pressure, sep='\n')
    CO2_rate = sm.collect_observations_as_time_series(shc, 'loc_CO2_aquifer')
    print('Rate at location 3:', CO2_rate, sep='\n')
    total_CO2_rate = sm.collect_observations_as_time_series(
        shc, 'CO2_aquifer_total')
    print('Cumulative CO2 rate at all 15 locations:', total_CO2_rate, sep='\n')
    total_brine_rate = sm.collect_observations_as_time_series(
        shc, 'brine_aquifer_total')
    print('Cumulative brine rate at all 15 locations:', total_brine_rate, sep='\n')

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12, 8), sharex=True)
    ax[0, 0].plot(time_array/365.25, pressure/1.0e+6, label='Pressure')
    ax[0, 1].plot(time_array/365.25, CO2_rate, label='Local CO2 rate')
    ax[1, 0].plot(time_array/365.25, total_CO2_rate, label='Total CO2 rate')
    ax[1, 1].plot(time_array/365.25, total_brine_rate, label='Total brine rate')
    ax[1, 0].set_xlabel('Time, [years]')
    ax[1, 1].set_xlabel('Time, [years]')
    ax[0, 0].set_ylabel('Pressure, [MPa]')
    ax[0, 1].set_ylabel(r'CO$_2$ saturation, [-]')
    ax[1, 0].set_ylabel(r'CO$_2$ leakage rates, [kg/s]')
    ax[1, 1].set_ylabel('Brine leakage rates, [kg/s]')
    fig.subplots_adjust(left=0.1, right=0.98, top=0.95, bottom=0.1, wspace=0.2)

    # Output:
    # Pressure at location 3:
    # [26018000. 26827000. 27092000. 27226000. 27360000. 27494000.]
    # Rate at location 3:
    # [0. 0. 0. 0. 0. 0.]
    # Cumulative CO2 rate at all 15 locations:
    # [0. 0. 0. 0. 0. 0.]
    # Cumulative brine rate at all 15 locations:
    # [0.00000000e+00 4.98252694e-05 5.07288455e-05 5.12306784e-05
    #  5.17325113e-05 5.22343442e-05]


def test_scenario2():
    """
    Run script example illustrating work of Seal Horizon component.

    Simple example for one cell.
    """
    # Time constants.
    num_years = 1
    time_array = 365.25 * np.arange(0.0, num_years+1)
    seal_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model.
    sm = SystemModel(model_kwargs=seal_model_kwargs)
    shc = sm.add_component_model_object(
        SealHorizon(name='shc', parent=sm,
                    locX=[50.0], locY=[550.0], area=[100.0]))

    # Parameters to vary. Arrays size of keyword arguments coincide with the size
    # of arrays containing coordinates of the cell coordinates
    shc.model_kwargs['thickness'] = [92.89171013411794]
    shc.model_kwargs['permeability'] = [2.5e-2 * sunit.microd_to_metersq()]
    shc.model_kwargs['baseDepth'] = [1037.108289865882+92.89171013411794]

    shc.add_par('aveBaseDepth', value=1130.0, vary=False)
    shc.add_par('aveBasePressure', value=32000000.0, vary=False)
    shc.add_par('aveTemperature', value=50.0, vary=False)
    shc.add_par('salinity', value=15000.0, vary=False)
    shc.add_par('staticDepth', value=1000.0, vary=False)
    shc.add_par('staticPressure', value=10000000.0, vary=False)

    # Fluid parameters
    shc.add_par('brineDensity', value=1004.0, vary=False)
    shc.add_par('CO2Density', value=597.8, vary=False)
    shc.add_par('brineViscosity', value=5.634e-4, vary=False)
    shc.add_par('CO2Viscosity', value=4.52e-5, vary=False)
    shc.add_par('CO2Solubility', value=0.0, vary=False)

    # Two-phase model parameters for L-E-T model
    shc.add_par('wetting1', value=1.0, vary=False)
    shc.add_par('wetting2', value=10.0, vary=False)
    shc.add_par('wetting3', value=1.25, vary=False)
    shc.add_par('nonwet1', value=1.05, vary=False)
    shc.add_par('nonwet2', value=3.5, vary=False)
    shc.add_par('nonwet3', value=1.25, vary=False)
    shc.add_par('capillary1', value=0.2, vary=False)
    shc.add_par('capillary2', value=2.8, vary=False)
    shc.add_par('capillary3', value=0.43, vary=False)
    shc.add_par('maxCapillary', value=10000000.0, vary=False)

    # Additional parameters for two-phase
    shc.add_par('permRatio', value=0.6, vary=False)
    shc.add_par('entryPressure', value=5.0e+3, vary=False)
    shc.add_par('brineResSaturation', value=0.15)
    shc.add_par('CO2ResSaturation', value=0.0)

    # Time-model parameters
    # influenceModel cannot be random, so vary=False
    shc.add_par('influenceModel', value=0, vary=False)
    # The next parameters are only needed if influenceModel is 1 or 2
    # shc.add_par('totalEffect', value=0.5, vary=False)
    # shc.add_par('rateEffect', value=0.1, vary=False)
    # shc.add_par('reactivity', value=8.0, vary=False)
    # shc.add_par('carbonateContent', value=5.0, vary=False)
    # shc.add_par('clayContent', value=60.0, vary=False)
    # shc.model_kwargs['clayType'] = 'smectite'

    shc.add_dynamic_kwarg('CO2saturation', [[0.5], [0.5]])
    shc.add_dynamic_kwarg('pressure', [[3.20E+07], [3.20E+07]])

    shc.add_obs('CO2_aquifer_total')
    shc.add_obs('brine_aquifer_total')
    shc.add_obs('mass_CO2_aquifer_total')
    shc.add_obs('mass_brine_aquifer_total')

    sm.forward()

    print('Cumulative CO2 flow rate:',
          sm.collect_observations_as_time_series(shc, 'CO2_aquifer_total'))
    print('Cumulative brine flow rate:',
          sm.collect_observations_as_time_series(shc, 'brine_aquifer_total'))
    print('Cumulative CO2 mass:',
          sm.collect_observations_as_time_series(shc, 'mass_CO2_aquifer_total'))
    print('Cumulative brine mass:',
          sm.collect_observations_as_time_series(shc, 'mass_brine_aquifer_total'))

    # The results should be
    # Cumulative CO2 flow rate: [0.00000000e+00 1.25034391e-06]
    # Cumulative brine flow rate: [0.00000000e+00 7.25886042e-08]
    # Cumulative CO2 mass: [ 0.         39.45785296]
    # Cumulative brine mass: [0.         2.29072213]


if __name__ == "__main__":

    # Setup logging and constants.
    logging.basicConfig(level=logging.DEBUG)
    test_case = 1  # 2 is a maximum number of test cases

    if test_case == 1:
        test_scenario1()

    elif test_case == 2:
        test_scenario2()
