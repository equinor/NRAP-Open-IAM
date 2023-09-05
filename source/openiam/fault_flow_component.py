# -*- coding: utf-8 -*-
# @author: Veronika Vasylkivska, Ernest Lindner
# September, 2022
import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt


source_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(source_folder)
sys.path.append(os.path.join(source_folder, 'components', 'fault', 'fault_flow'))

try:
    from openiam import SystemModel, ComponentModel
    from openiam.cfi.commons import process_dynamic_inputs
except ImportError as err:
    print('Unable to load IAM class module: '+ str(err))

try:
    from components.fault.fault_flow.fault_setup import FAULT_SETUP_DICT
    import components.fault.fault_flow.flt_basics as bcs
    import components.fault.fault_flow.flt_compute as fcp
    import components.fault.fault_flow.flt_decipher as deci
    import components.fault.fault_flow.flt_intro as intro
    import components.fault.fault_flow.flt_profile as pro
    import components.fault.fault_flow.flt_units as fit
    import components.fault.fault_flow.flt_config as conf
except ImportError:
    print('\nERROR: Unable to load ROM for Fault Flow component\n')
    sys.exit()


FF_SCALAR_OBSERVATIONS = ['CO2_aquifer_total', 'brine_aquifer_total',
                          'mass_CO2_aquifer_total', 'mass_brine_aquifer_total']
FF_GRID_OBSERVATIONS = ['CO2_aquifer', 'brine_aquifer',
                        'mass_CO2_aquifer', 'mass_brine_aquifer']

FF_CFI_CONTROLS = {'profile_type': ['profileType', 0],
                   'near_surface_approach': ['considerNearApproach', False],
                   'pressure_approach': ['pressureApproach', False],
                   'interpolate_approach': ['interpolateApproach', False]}

class FaultFlow(ComponentModel):
    """
    The Fault Flow component model simulates the flow of carbon dioxide along
    a low permeability fault from an injection horizon (into which carbon dioxide
    is injected) up to a freshwater aquifer. The theoretical base is predicated
    on one-dimension (1D), steady-state, two-phase flow of |CO2| through
    a saturated discontinuity (parallel plates) under |CO2| supercritical
    conditions. The flow in the current implementation uses the near-surface
    |CO2| supercritical point to be the upper point of flow. The surrounding
    rock matrix is considered relatively impermeable.

    In the NRAP-Open-IAM control file, the type name for the Fault Flow component is
    ``FaultFlow``. The description of the component's parameters is provided below:

    Fault core setup parameters:

    * **strike** [|deg|] (0 to 360) - direction of fault: trend of fault strike
      taken clockwise from north (default: 30)

    * **dip** [|deg|] (10 to 90) - inclination of fault plane from strike, using
      right-hand rule from strike (default: 70)

    * **length** [|m|] (0 to 10,000) - length of fault trace at surface from start
      point (default: 100)

    * **xStart** [|m|] (-5.0e+07 to 5.0e+07) - x-coordinate of the fault start point
      taken as the left point on fault trace (default: 500)

    * **yStart** [|m|] (-5.0e+07 to 5.0e+07) - y-coordinate of the fault start point
      taken as the left point on fault trace (default: 500)

    * **nSegments** [-] (1 to 100) - number of separate fault divisions
      of the fault (default: 4)

    * **faultProbability** [%] (0 to 100) - probability of fault existence
      (default: 100)

    Fault aperture setup parameters:

    * **aperture** [|m|] (0 to 0.05) - effective aperture of fault
      (default: 2.5e-6)

    * **SGR** [-] (0 to 100) - shale gouge ratio for fault (default: 0)

    * **stateVariable** [-] (0 to 1) - correction factor for near-surface flow
      (default: 1)

    For variability of fault properties and orientation setup one can use the
    following eight distribution parameters.

    Note: The setup of the eight distribution parameters below is not yet implemented
    in the control file interface or GUI of NRAP-Open-IAM and available
    only in the script interface.

    Strike distribution parameters:

    * **aveStrike** [|deg|] (0 to 360) - average direction of fault: average trend
      of fault strike taken clockwise from north (default: 90); also default value for
      no variation in strike

    * **spreadStrike** [|deg|] (0 to 180) - spread in strike orientation
      (range of 2-sigma around average) (default: 0)

    Dip distribution parameters:

    * **aveDip** [|deg|] (10 to 90) - average inclination of fault plane from
      strike, using right-hand rule from strike (default: 90)

    * **stdDevDip** [|deg|] (0 to 90) - standard deviation of angle of dip (default: 0)

    Aperture distribution parameters:

    * **aveAperture** [|m|] (0 to 1.01e-1) - average effective aperture of fault
      (default: 1.0e-2)

    * **stdDevAperture** [|m|] (0 to 2.0e-2) - standard deviation of effective aperture
      (default: 0.0)

    * **minAperture** [|m|] (0 to 1.0e-3) - minimum aperture (default: 1.0e-7)

    * **maxAperture** [|m|] (0 to 5.0e-2) - maximum aperture (default: 2.0e-2)

    Field parameters:

    * **aquiferDepth** [|m|] (200 to 2,000) - depth to base of deepest aquifer
      along/above fault (default: 240)

    * **aquiferTemperature** [|C|] (15 to 180) - temperature of brine of deepest
      aquifer at base (default: 22)

    * **aquiferPressure** [|Pa|] (1.0e+6 to 6.0e+8) - pressure at base of aquifer
      (default: 1.42E+07)

    * **injectDepth** [|m|] (860 to 20,000) - reference depth positive below grade
      to top of injection horizon (default: 1880)

    * **injectTemperature** [|C|] (31 to 180) - average temperature of brine
      at injection depth in reservoir (default: 95)

    * **fieldPressure** [|Pa|] (1.0e+5 to 6.0e+7) - initial pressure at injection
      depth before injection starts (default: 1.9140e+07)

    * **injectPressure** [|Pa|] (7.0e+6 to 6.0e+8) - average pressure
      at base during injection period for interpolation of viscosity and density
      (default: 2.9290E+07)

    * **finalPressure** [|Pa|] (1.0e+5 to 6.0e+7) - final average pressure
      at injection depth for interpolation of viscosity and density
      (default: 1.9140e+07)

    * **injectX** [|m|] (-5.0e+07 to 5.0e+07) - x-coordinate of the location
      of injection well (default: 0)

    * **injectY** [|m|] (-5.0e+07 to 5.0e+07) - y-coordinate of the location
      of injection well (default: 0)

    * **injectEndTime** [years] (0 to 10000) - time when injection stops
        (default: 50)

    Reservoir conditions parameters:

    * **salinity** [|ppm|] (0 to 80000) - salinity of the brine (default: 0).
      The value is used to compute density and viscosity of the brine

    * **CO2Density** [|kg/m^3|] (93 to 1050) - average density of |CO2| phase
      for fault (default: 673.84). The value is used if interpolation
      is not conducted by code

    * **CO2Viscosity** [|Pa*s|] (1.8e-05 to 1.4e-04) - viscosity of |CO2| phase
      for fault (default: 5.5173e-05). The value is used if interpolation
      is not conducted by code

    * **brineDensity** [|kg/m^3|] (880 to 1080) - density of brine phase for fault
      (default: 974.895). The value is used if interpolation is not conducted by code

    * **brineViscosity** [|Pa*s|] (1.5e-04 to 1.6e-03) - viscosity of brine phase
      for fault (default: 3.0491e-04). The value is used if interpolation
      is not conducted by code

    * **CO2Solubility** [|mol/kg|] (0 to 2) - solubility of |CO2| phase in brine
      for fault (default: 0.035). The value is used if interpolation is not conducted
      by code

    Aquifer conditions parameters:

    * **aquiferCO2Density** [|kg/m^3|] (93 to 1050) - density of |CO2| phase
      in the aquifer (default: 886.44)

    * **aquiferCO2Viscosity** [|Pa*s|] (1.1e-05 to 1.4e-04) - viscosity of |CO2| phase
      in the aquifer (default: 8.8010e-05)

    * **aquiferBrineDensity** [|kg/m^3|] (880 to 1080) - density of brine phase
      in the aquifer (default: 1004.10)

    * **aquiferBrineViscosity** [|Pa*s|] (1.5e-04 to 1.6e-03) - viscosity of brine phase
      in the aquifer (default: 3.0221e-04)

    Relative flow parameters:

    * **brineResSaturation** [-] (0.01 to 0.35) - residual wetting brine saturation
      used in two-phase model (default: 0.15)

    * **CO2ResSaturation** [-] (0 to 0.35) - residual nonwetting |CO2| saturation
      used in two-phase model (default: 0.0)

    * **relativeModel** [-] (LET or BC) - relative permeability model (default: LET)

    * **permRatio** [-] (0 to 1.5) - ratio of maximum nonwetting permeability
      to the maximum wetting permeability (default: 0.6)

    * **entryPressure** [|Pa|] (100 to 2.0e+6) - entry/threshold/bubbling pressure
      that controls flow into rock (default: 5000)

    Two-phase model parameters for LET model:

    * **wetting1** [-] (0.5 to 5) - wetting phase parameter |L| (default: 1)

    * **wetting2** [-] (1 to 30) - wetting phase parameter |E| (default: 10)

    * **wetting3** [-] (0 to 3) - wetting phase parameter |T| (default: 1.25)

    * **nonwet1** [-] (0.5 to 5) - nonwetting phase parameter |L| (default: 1.05)

    * **nonwet2** [-] (1 to 30) - nonwetting phase parameter |E| (default: 10)

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
    and **maxCapillary** are used only if parameter **relativeModel**
    is set to *LET*.

    BC model parameters:

    * **lambda** [-] (0 to 5) - lambda term in Brooks-Corey model (default: 2.5)

    Note: Parameter **lambda** is used only if parameter **relativeModel**
    is set to *BC*.

    Stress parameters:

    * **maxHorizontal** [|Pa|] (0 to 5.0e+7) - secondary maximum horizontal
      principal stress at top of injection horizon (default: 3.0e+7)

    * **minHorizontal** [|Pa|] (0 to 5.0e+7) - secondary minimum horizontal
      principal stress at top of injection interval (default: 2.0e+7)

    * **maxTrend** [|deg|] (0 to 180) - strike of secondary maximum horizontal
      stress clockwise from north (default: 55)

    The possible outputs from the Fault Flow component are
    leakage rates of |CO2| and brine to aquifer through fault. The names
    of the observations are of the form:

    * **CO2_aquifer**, **brine_aquifer** [|kg/s|] - |CO2| and brine leakage rates to
      aquifer through fault (individual segments) into overlying aquifer

    * **mass_CO2_aquifer**, **mass_brine_aquifer** [|kg|] - mass of the |CO2|
      and brine through fault (individual segments) to overlying aquifer

    * **CO2_aquifer_total**, **brine_aquifer_total** [|kg/s|] - cumulative |CO2| and brine
      leakage rates to aquifer through fault into overlying aquifer

    * **mass_CO2_aquifer_total**, **mass_brine_aquifer_total** [|kg|] - cumulative mass
      of the |CO2| and brine through fault (individual cells) to overlying aquifer.

    Observations with names CO2_aquifer, brine_aquifer, mass_CO2_aquifer and
    mass_brine_aquifer are provided as arrays of values of length equal
    to the number of fault segments. To output observations corresponding to
    a particular fault segment (e.g., segment 1) one can add observations
    with names CO2_aquifer_segm# where # is an index of a segment of interest
    (e.g., CO2_aquifer_segm1) to the output of the Fault Flow component.
    """
    def __init__(self, name, parent, **kwargs):
        """
        Constructor method of FaultFlow class

        :param name: name of component model
        :type name: [str]

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel [object]

        :returns: FaultFlow class object
        """
        model_kwargs = {'time_point': 365.25, 'time_step': 365.25}

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'FaultFlow'

        # Copy default setup dictionary
        self.fault_controls = {}
        self.fault_controls = deci.interpret_fault_parameters(
            FAULT_SETUP_DICT, self.fault_controls)

        self.process_additional_input(**kwargs)

        # Define parameters and keyword arguments boundaries
        self.param_bounds = intro.define_input_limits()
        self.define_pars_bounds()

        # Set default parameters
        self.setup_default_pars()

        # Reserve space for fault variable
        self.fault = list()

        # Add accumulators for masses for the first segment
        self.add_accumulator('mass_CO2_aquifer_segm1', sim=0.0)
        self.add_accumulator('mass_brine_aquifer_segm1', sim=0.0)

        # Define gridded observations names
        self.grid_obs_keys = FF_GRID_OBSERVATIONS

    def process_additional_input(self, **kwargs):
        """
        Process additional input keyword arguments provided to the constructor method.

        Possible keys: 'profile_type' (default: 0), 'near_surface_approach'
        (default: False), 'pressure_approach' (default: False),
        'interpolate_approach' (default: False)

        Parameters
        ----------
        **kwargs : dict
            Dictionary containing additional input provided to the constructor.

        Returns
        -------
        None.
        """
        for key, val in FF_CFI_CONTROLS.items():
            self.fault_controls[key] = kwargs.get(key, val[1])

        if self.fault_controls['profile_type'] not in [0, 1, 2]:
            warn_msg = ''.join([
                'Input argument profile_type to component {} has invalid ',
                'value of {}. It will be replaced by default value of 0.']).format(
                    self.name, self.fault_controls['profile_type'])
            logging.warning(warn_msg)
            self.fault_controls['profile_type'] = 0

    def setup_default_pars(self):
        """
        Add parameters of the component with default values.
        """
        # HINT: Since parameters strike, dip and aperture might assume a particular
        # distribution defined by the standalone code there are no default
        # parameter setup for these parameters. We process them separately,
        # in a different way from the rest of parameters.

        # Fault core parameters
        self.add_default_par('faultProbability',
                             value=FAULT_SETUP_DICT['FaultCore']['faultProbability'])
        self.add_default_par(
            'aveStrike', value=FAULT_SETUP_DICT['FaultCore']['aveStrike'])
        self.add_default_par(
            'spreadStrike', value=FAULT_SETUP_DICT['FaultCore']['spreadStrike'])
        self.add_default_par(
            'aveDip', value=FAULT_SETUP_DICT['FaultCore']['aveDip'])
        # Parameter stdDevDip have different name in the FAULT_SETUP_DICT.
        # For consistency across parameter naming conventions we're using
        # stdDevDip instead of just stdDip
        self.add_default_par(
            'stdDevDip', value=FAULT_SETUP_DICT['FaultCore']['stdDip'])
        self.add_default_par(
            'length', value=FAULT_SETUP_DICT['FaultCore']['length']['valu'])
        self.add_default_par(
            'xStart', value=FAULT_SETUP_DICT['FaultCore']['xStart']['valu'])
        self.add_default_par(
            'yStart', value=FAULT_SETUP_DICT['FaultCore']['yStart']['valu'])
        self.add_default_par(
            'nSegments', value=FAULT_SETUP_DICT['FaultCore']['nSegments'])

        # Field parameters
        self.add_default_par(
            'aquiferDepth', value=FAULT_SETUP_DICT['Field']['aquiferDepth']['valu'])
        self.add_default_par(
            'aquiferTemperature',
            value=FAULT_SETUP_DICT['Field']['aquiferTemperature']['valu'])
        self.add_default_par(
            'aquiferPressure',
            value=FAULT_SETUP_DICT['Field']['aquiferPressure']['valu'])
        self.add_default_par(
            'injectDepth', value=FAULT_SETUP_DICT['Field']['injectDepth']['valu'])
        self.add_default_par(
            'injectTemperature',
            value=FAULT_SETUP_DICT['Field']['injectTemperature']['valu'])
        self.add_default_par(
            'fieldPressure',
            value=FAULT_SETUP_DICT['Field']['fieldPressure']['valu'])
        self.add_default_par(
            'injectPressure',
            value=FAULT_SETUP_DICT['Field']['injectPressure']['valu'])
        self.add_default_par(
            'finalPressure',
            value=FAULT_SETUP_DICT['Field']['finalPressure']['valu'])
        self.add_default_par(
            'injectX', value=FAULT_SETUP_DICT['Field']['injectX']['valu'])
        self.add_default_par(
            'injectY', value=FAULT_SETUP_DICT['Field']['injectY']['valu'])
        self.add_default_par(
            'injectEndTime', value=FAULT_SETUP_DICT['ModelParams']['injectEndTime'])

        # Aperture setup parameters
        self.add_default_par('SGR', value=FAULT_SETUP_DICT['Aperture']['SGR'])
        self.add_default_par(
            'stateVariable', value=FAULT_SETUP_DICT['Aperture']['stateVariable'])
        scalar = fit.mm_to_m()
        self.add_default_par(
            'aveAperture',
            value=scalar*FAULT_SETUP_DICT['Aperture']['aveAperture']['valu'])
        self.add_default_par(
            'stdDevAperture',
            value=scalar*FAULT_SETUP_DICT['Aperture']['stdDevAperture']['valu'])
        self.add_default_par(
            'minAperture',
            value=scalar*FAULT_SETUP_DICT['Aperture']['minAperture']['valu'])
        self.add_default_par(
            'maxAperture',
            value=scalar*FAULT_SETUP_DICT['Aperture']['maxAperture']['valu'])

        # Reservoir conditions parameters
        self.add_default_par(
            'salinity', value=FAULT_SETUP_DICT['InjectConditions']['salinity'])
        self.add_default_par(
            'CO2Density',
            value=FAULT_SETUP_DICT['InjectConditions']['CO2Density']['valu'])
        self.add_default_par(
            'CO2Viscosity',
            value=FAULT_SETUP_DICT['InjectConditions']['CO2Viscosity']['valu'])
        self.add_default_par(
            'brineDensity',
            value=FAULT_SETUP_DICT['InjectConditions']['brineDensity']['valu'])
        self.add_default_par(
            'brineViscosity',
            value=FAULT_SETUP_DICT['InjectConditions']['brineViscosity']['valu'])
        self.add_default_par(
            'CO2Solubility',
            value=FAULT_SETUP_DICT['InjectConditions']['CO2Solubility']['valu'])

        # Aquifer conditions parameters
        # Since there are similar parameters for reservoir conditions I've added
        # 'aquifer' at the beginning of each parameter to distinguish them from
        # reservoir related parameters
        self.add_default_par(
            'aquiferCO2Density',
            value=FAULT_SETUP_DICT['AquiferConditions']['CO2Density']['valu'])
        self.add_default_par(
            'aquiferCO2Viscosity',
            value=FAULT_SETUP_DICT['AquiferConditions']['CO2Viscosity']['valu'])
        self.add_default_par(
            'aquiferBrineDensity',
            value=FAULT_SETUP_DICT['AquiferConditions']['brineDensity']['valu'])
        self.add_default_par(
            'aquiferBrineViscosity',
            value=FAULT_SETUP_DICT['AquiferConditions']['brineViscosity']['valu'])

        # Relative flow parameters
        self.add_default_par(
            'brineResSaturation',
            value=FAULT_SETUP_DICT['RelativeFlowLimits']['brineResSaturation'])
        self.add_default_par(
            'CO2ResSaturation',
            value=FAULT_SETUP_DICT['RelativeFlowLimits']['CO2ResSaturation'])
        self.add_default_par(
            'permRatio',
            value=FAULT_SETUP_DICT['RelativeFlowLimits']['permRatio'])
        self.add_default_par(
            'entryPressure',
            value=FAULT_SETUP_DICT['RelativeFlowLimits']['entryPressure']['valu'])

        # String variable cannot be parameter so the only way to pass
        # string type parameters is through keyword arguments to the model method
        # BC, default value
        self.model_kwargs['relativeModel'] = \
            FAULT_SETUP_DICT['RelativeFlowLimits']['relativeModel']
        self.add_default_par(
            'lambda', value=FAULT_SETUP_DICT['BrooksCoreyModel']['lambda'])

        # Two-phase model parameters for L-E-T model
        self.add_default_par('wetting1',
                             value=FAULT_SETUP_DICT['LETModel']['wetting1'])
        self.add_default_par('wetting2',
                             value=FAULT_SETUP_DICT['LETModel']['wetting2'])
        self.add_default_par('wetting3',
                             value=FAULT_SETUP_DICT['LETModel']['wetting3'])
        self.add_default_par('nonwet1',
                             value=FAULT_SETUP_DICT['LETModel']['nonwet1'])
        self.add_default_par('nonwet2',
                             value=FAULT_SETUP_DICT['LETModel']['nonwet2'])
        self.add_default_par('nonwet3',
                             value=FAULT_SETUP_DICT['LETModel']['nonwet3'])
        self.add_default_par(
            'capillary1', value=FAULT_SETUP_DICT['LETCapillaryModel']['capillary1'])
        self.add_default_par(
            'capillary2', value=FAULT_SETUP_DICT['LETCapillaryModel']['capillary2'])
        self.add_default_par(
            'capillary3', value=FAULT_SETUP_DICT['LETCapillaryModel']['capillary3'])
        self.add_default_par(
            'maxCapillary',
            value=FAULT_SETUP_DICT['LETCapillaryModel']['maxCapillary']['valu'])

        # Stress parameters
        self.add_default_par(
            'maxHorizontal', value=FAULT_SETUP_DICT['Stress']['maxHorizontal']['valu'])
        self.add_default_par(
            'minHorizontal', value=FAULT_SETUP_DICT['Stress']['minHorizontal']['valu'])
        self.add_default_par(
            'maxTrend', value=FAULT_SETUP_DICT['Stress']['maxTrend'])

        # Seed parameter
        self.add_default_par('seed', value=conf.SEEDX)

    def define_pars_bounds(self):
        """ Define dictionaries for limits of parameters and temporal inputs."""
        # Define dictionary of boundaries
        self.pars_bounds = dict()

        # Controls
        # For now we don't check these limits: as they relate to the model setup
        self.pars_bounds['startTime'] = self.param_bounds['start_time']
        self.pars_bounds['endTime'] = self.param_bounds['end_time']
        self.pars_bounds['injectEndTime'] = self.param_bounds['inject_end']
        self.pars_bounds['timePoints'] = self.param_bounds['time_points']
        self.pars_bounds['realizations'] = self.param_bounds['realizations']

        # Fault core parameters boundaries
        self.pars_bounds['faultProbability'] = self.param_bounds['fault_probability']
        self.pars_bounds['aveStrike'] = self.param_bounds['strike_mean']
        self.pars_bounds['spreadStrike'] = self.param_bounds['strike_sig']
        self.pars_bounds['strike'] = self.param_bounds['strike_mean']
        self.pars_bounds['aveDip'] = self.param_bounds['dip_mean']
        self.pars_bounds['stdDevDip'] = self.param_bounds['dip_std']
        self.pars_bounds['dip'] = self.param_bounds['dip_mean']
        self.pars_bounds['length'] = self.param_bounds['length']
        self.pars_bounds['xStart'] = self.param_bounds['x_start']
        self.pars_bounds['yStart'] = self.param_bounds['y_start']
        self.pars_bounds['nSegments'] = self.param_bounds['n_plates']
        self.pars_bounds['SGR'] = self.param_bounds['sgr']
        self.pars_bounds['stateVariable'] = self.param_bounds['state_correction']

        # Aperture parameters
        scalar = fit.mm_to_m()
        self.pars_bounds['aveAperture'] = [
            scalar*val for val in self.param_bounds['aperture_mean']]
        self.pars_bounds['stdDevAperture'] = [
            scalar*val for val in self.param_bounds['aperture_std']]
        self.pars_bounds['minAperture'] = [
            scalar*val for val in self.param_bounds['aperture_min']]
        self.pars_bounds['maxAperture'] = [
            scalar*val for val in self.param_bounds['aperture_max']]
        # HINT: In the standalone code the aperture boundaries entered by user
        # are smaller than the ones generated from distribution. We might
        # consider changing that in the future and bringing them closer to the
        # latter.
        self.pars_bounds['aperture'] = [
            scalar*val for val in self.param_bounds['aperture_confines']]

        # Field parameters
        self.pars_bounds['aquiferDepth'] = self.param_bounds['aquifer_depth']
        self.pars_bounds['aquiferPressure'] = self.param_bounds['aquifer_pressure']
        self.pars_bounds['aquiferTemperature'] = self.param_bounds['aquifer_temperature']
        self.pars_bounds['injectDepth'] = self.param_bounds['inject_depth']
        self.pars_bounds['fieldPressure'] = self.param_bounds['field_pressure']
        self.pars_bounds['injectPressure'] = self.param_bounds['inject_pressure']
        self.pars_bounds['finalPressure'] = self.param_bounds['final_pressure']
        self.pars_bounds['injectTemperature'] = self.param_bounds['inject_temperature']
        self.pars_bounds['injectX'] = self.param_bounds['inject_x']
        self.pars_bounds['injectY'] = self.param_bounds['inject_y']

        # Reservoir conditions parameters
        self.pars_bounds['salinity'] = self.param_bounds['salinity']
        self.pars_bounds['CO2Density'] = self.param_bounds['co2_density']
        self.pars_bounds['CO2Viscosity'] = self.param_bounds['co2_viscosity']
        self.pars_bounds['brineDensity'] = self.param_bounds['brine_density']
        self.pars_bounds['brineViscosity'] = self.param_bounds['brine_viscosity']
        self.pars_bounds['CO2Solubility'] = self.param_bounds['co2_solubility']

        # Aquifer conditions parameters
        self.pars_bounds['aquiferCO2Density'] = self.param_bounds['aqui_co2_density']
        self.pars_bounds['aquiferCO2Viscosity'] = self.param_bounds['aqui_co2_viscosity']
        self.pars_bounds['aquiferBrineDensity'] = self.param_bounds['aqui_brine_density']
        self.pars_bounds['aquiferBrineViscosity'] = self.param_bounds['aqui_brine_viscosity']

        # Relative flow parameters
        self.pars_bounds['brineResSaturation'] = self.param_bounds['resid_brine']
        self.pars_bounds['CO2ResSaturation'] = self.param_bounds['resid_co2']
        self.pars_bounds['permRatio'] = self.param_bounds['perm_ratio']
        # HINT: Boundaries for parameter/keyword argument relativeModel is a list of
        # possible values, not [a lower boundary, an upper boundary] list
        # Thus, it needs different treatment when checking for valid values
        self.pars_bounds['relativeModel'] = ['BC', 'LET']
        self.pars_bounds['entryPressure'] = self.param_bounds['entry_pressure']
        self.pars_bounds['lambda'] = self.param_bounds['zeta']

        # Stress parameters
        self.pars_bounds['maxHorizontal'] = self.param_bounds['max_horizontal']
        self.pars_bounds['minHorizontal'] = self.param_bounds['min_horizontal']
        self.pars_bounds['maxTrend'] = self.param_bounds['max_trend']

        # Two-phase model parameters for L-E-T model
        self.pars_bounds['wetting1'] = self.param_bounds['l_wetting']
        self.pars_bounds['wetting2'] = self.param_bounds['e_wetting']
        self.pars_bounds['wetting3'] = self.param_bounds['t_wetting']
        self.pars_bounds['nonwet1'] = self.param_bounds['l_nonwet']
        self.pars_bounds['nonwet2'] = self.param_bounds['e_nonwet']
        self.pars_bounds['nonwet3'] = self.param_bounds['t_nonwet']
        self.pars_bounds['capillary1'] = self.param_bounds['l_capillary']
        self.pars_bounds['capillary2'] = self.param_bounds['e_capillary']
        self.pars_bounds['capillary3'] = self.param_bounds['t_capillary']
        self.pars_bounds['maxCapillary'] = self.param_bounds['max_capillary']

        self.pars_bounds['coordinates'] = self.param_bounds['coordinates']

        # Seed boundaries
        self.pars_bounds['seed'] = [1, np.inf]

        # Define dictionary of temporal data limits
        self.temp_data_bounds = dict()
        self.temp_data_bounds['pressure'] = ['Pressure']+self.param_bounds['reservoir_pressure']
        self.temp_data_bounds['CO2saturation'] = (
            ['CO2 saturation']+self.param_bounds['reservoir_saturation'])

    def update_fault_controls(self, p, **kwarg_p):
        """ Update fault_controls with values provided by user."""

        # Model parameters
        if self._parent.time_array is not None:
            self.fault_controls['start_time'] = self._parent.time_array[0] # first time point
            self.fault_controls['end_time'] = self._parent.time_array[-1]  # last time point
            self.fault_controls['time_points'] = len(self._parent.time_array)
        else:
            err_msg = 'Argument time_array has to be defined for system model.'
            logging.error(err_msg)

        # Injection period; tracked for interpolation purposes
        self.fault_controls['inject_end'] = p['injectEndTime']

        # Fault core parameters
        # Process strike data
        if 'strike' in p:   # constant value, no distribution
            self.fault_controls['strike_mean'] = p['strike']
            self.fault_controls['strike_sig'] = 0.0
            self.fault_controls['strike_approach'] = False
        else:               # variability in strike is assumed
            self.fault_controls['strike_mean'] = p['aveStrike']
            self.fault_controls['strike_sig'] = p['spreadStrike']
            self.fault_controls['strike_approach'] = True

        # Process dip data
        if 'dip' in p:   # constant value, no distribution
            self.fault_controls['dip_mean'] = p['dip']
            self.fault_controls['dip_std'] = 0.0
            self.fault_controls['dip_approach'] = False
        else:            # variability in dip is assumed
            self.fault_controls['dip_mean'] = p['aveDip']
            self.fault_controls['dip_std'] = p['stdDevDip']
            self.fault_controls['dip_approach'] = True

        self.fault_controls['length'] = p['length']
        self.fault_controls['x_start'] = p['xStart']
        self.fault_controls['y_start'] = p['yStart']
        self.fault_controls['n_plates'] = p['nSegments']
        self.fault_controls['sgr'] = p['SGR']
        self.fault_controls['state_variable'] = p['stateVariable']
        self.fault_controls['seg_length'] = (
            self.fault_controls['length']/self.fault_controls['n_plates'])

        # Field parameters
        self.fault_controls['aquifer_depth'] = p['aquiferDepth']
        self.fault_controls['inject_depth'] = p['injectDepth']
        self.fault_controls['aquifer_pressure'] = p['aquiferPressure']
        self.fault_controls['field_pressure'] = p['fieldPressure']
        self.fault_controls['final_pressure'] = p['finalPressure']
        self.fault_controls['inject_pressure'] = p['injectPressure']
        self.fault_controls['aquifer_temperature'] = p['aquiferTemperature']
        self.fault_controls['inject_temperature'] = p['injectTemperature']
        self.fault_controls['inject_x'] = p['injectX']
        self.fault_controls['inject_y'] = p['injectY']

        # Conditions parameters
        self.fault_controls['salinity'] = p['salinity']
        self.fault_controls['co2_density'] = p['CO2Density']
        self.fault_controls['co2_viscosity'] = p['CO2Viscosity']
        self.fault_controls['brine_density'] = p['brineDensity']
        self.fault_controls['brine_viscosity'] = p['brineViscosity']
        self.fault_controls['co2_solubility'] = p['CO2Solubility']

        # Check type of model for relative permeability
        if kwarg_p['relativeModel'] == 'BC':
            # BrooksCoreyModel:
            self.fault_controls['zeta'] = p['lambda']
        else:
            # LETModel parameters
            self.fault_controls['l_wetting'] = p['wetting1']
            self.fault_controls['e_wetting'] = p['wetting2']
            self.fault_controls['t_wetting'] = p['wetting3']
            self.fault_controls['l_nonwet'] = p['nonwet1']
            self.fault_controls['e_nonwet'] = p['nonwet2']
            self.fault_controls['t_nonwet'] = p['nonwet3']

            # LET CapillaryPressure:
            self.fault_controls['l_capillary'] = p['capillary1']
            self.fault_controls['e_capillary'] = p['capillary2']
            self.fault_controls['t_capillary'] = p['capillary3']
            self.fault_controls['max_capillary'] = p['maxCapillary']

        # Relative flow parameters
        self.fault_controls['resid_brine'] = p['brineResSaturation']
        self.fault_controls['resid_co2'] = p['CO2ResSaturation']
        self.fault_controls['relative_model'] = kwarg_p['relativeModel']
        self.fault_controls['perm_ratio'] = p['permRatio']
        self.fault_controls['entry_pressure'] = p['entryPressure']

        # Time parameters
        self.fault_controls['max_horizontal'] = p['maxHorizontal']
        self.fault_controls['min_horizontal'] = p['minHorizontal']
        self.fault_controls['max_trend'] = p['maxTrend']

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

        for key, val_array in temp_inputs.items():
            for val in val_array:
                if (val < self.temp_data_bounds[key][1]) or (
                        val > self.temp_data_bounds[key][2]):
                    err_msg = ''.join([
                        'Temporal input {} ({}) of FaultFlow component {} ',
                        'is outside the model range [{}, {}] at time t = {} days']).format(
                            self.temp_data_bounds[key][0], val, self.name,
                            self.temp_data_bounds[key][1],
                            self.temp_data_bounds[key][2], time)
                    logging.error(err_msg)
                    raise ValueError('Temporal inputs are outside the FaultFlow component limits.')

    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        msg = ('Input parameters of {name} component are {p}'.
               format(name=self.name, p=p))
        logging.debug(msg)

        for key, val in p.items():
            if key not in self.pars_bounds:
                msg = ('Parameter {key} not recognized as a FaultFlow input ' +
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
        Check fault parameters provided as keyword arguments.

        :param kw_pars: input parameters of component model
        :type :
        """
        if 'relativeModel' in kw_pars:
            if kw_pars['relativeModel'] not in self.pars_bounds['relativeModel']:
                msg = "Argument relativeModel has an unrecognized value: {}.".format(
                    kw_pars['relativeModel'])
                logging.warning(msg)

                msg = ("Argument relativeModel will be set to 'BC'.")
                logging.warning(msg)
                kw_pars['relativeModel'] = 'BC'
        else:
            kw_pars['relativeModel'] = 'BC'

        if 'aperture' in kw_pars:
            for ind, val in enumerate(kw_pars['aperture']):
                if ((val < self.pars_bounds['aperture'][0]) or
                        (val > self.pars_bounds['aperture'][1])):
                    msg = ''.join(['Parameter aperture is out of bounds ',
                                   'for fault segment {}.']).format(ind+1)
                    logging.warning(msg)

        return kw_pars

    def connect_with_system(self, component_data, name2obj_dict,
                            locations, system_adapters, **kwargs):
        """
        Code to add FaultFlow to system model for control file interface.

        :param component_data: Dictionary of component data including 'Connection'
        :type component_data: dict

        :returns: None
        """
        # Add number of segments
        self.add_par('nSegments', value=locations[self.name]['number'], vary=False)

        # Process controls
        if 'Controls' in component_data:
            controls_data = component_data['Controls']
            for key, val in FF_CFI_CONTROLS.items():
                controls_data[key] = controls_data.pop(val[0], val[1])
            self.process_additional_input(**controls_data)

        # Add parameters
        if ('Parameters' in component_data) and (component_data['Parameters']):
            # Check whether aperture is defined as list
            if 'aperture' in component_data['Parameters']:
                data = component_data['Parameters']['aperture']
                if isinstance(data, list): # if provided data is a list
                    self.model_kwargs['aperture'] = data
                else:
                    if not isinstance(data, dict):
                        component_data['Parameters']['aperture'] = {
                            'value': data,
                            'vary': False}
                    self.add_par('aperture',
                                 **component_data['Parameters']['aperture'])
                # Remove aperture data from component_data['Parameters'] dictionary
                # so it's not processed again in the next loop
                component_data['Parameters'].pop('aperture')

            for key in component_data['Parameters']:
                if not isinstance(component_data['Parameters'][key], dict):
                    component_data['Parameters'][key] = {
                        'value': component_data['Parameters'][key],
                        'vary': False}
                self.add_par(key, **component_data['Parameters'][key])

        # Check whether dynamic kwargs are provided
        process_dynamic_inputs(
            self, component_data, array_like=True, check_second_dim=True,
            size=locations[self.name]['number'],
            quantity_to_compare='Number of fault segments')

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
                    raise KeyError(err_msg) from None

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
                if obs_nm in FF_SCALAR_OBSERVATIONS:
                    new_comp_outputs.append(obs_nm)
                    self.add_obs(obs_nm)
                elif obs_nm in FF_GRID_OBSERVATIONS:
                    # Get number of segments
                    num_locs = locations[self.name]['number']
                    for ind in range(num_locs):
                        augm_obs_nm = obs_nm+'_segm_{}'.format(ind+1)
                        self.add_local_obs(augm_obs_nm, obs_nm, 'array', ind)
                        new_comp_outputs.append(augm_obs_nm)
                    # Add gridded observations
                    self.add_grid_obs(
                        obs_nm, constr_type='array', output_dir=kwargs['output_dir'])
                else:
                    warn_msg = ''.join([
                        '{} is not recognised as observation name ',
                        'of Fault Flow component {}.']).format(obs_nm, self.name)
                    logging.warning(warn_msg)

            component_data['Outputs'] = new_comp_outputs

    def define_fracto(self):
        """Provide the setup fault input values.

        This method is based on the method define_fracto from flt_compute module
        of the standalone seal flux rewritten without not applicable statements.
        """
        # Define log-normal parameters to compute apertures.
        if self.fault_controls['aperture_approach']:
            # No variability - data from file.
            self.fault_controls['aperture_mu'] = 0.0
            self.fault_controls['aperture_scale'] = 0.0
        else:
            # Setup parameters for variable aperture.
            cent, scale = \
                bcs.convert_lognorm_terms(self.fault_controls['aperture_mean'],
                                          self.fault_controls['aperture_std'])
            self.fault_controls['aperture_mu'] = cent
            self.fault_controls['aperture_scale'] = scale

        # Setup disperse for strike if needed.
        if self.fault_controls['strike_approach'] and \
                self.fault_controls['strike_sig'] > fcp.V_SMALL_ZERO:
            # Define variability.
            self.fault_controls['strike_disperse'] = \
                bcs.convert_kappa(self.fault_controls['strike_sig'])

        # Compute effects of SGR.
        self.fault_controls = fcp.compute_sgr_influence(self.fault_controls)

        # Define fault counter.
        self.fault_controls['simulations_with_fault'] = 0

        # Define aperture-pressure parameters, if they will be used.
        if self.fault_controls['pressure_approach']:
            self.fault_controls = fcp.define_apertus(self.fault_controls)

    def process_aperture_data(self, pars, **kwarg_p):

        # Aperture distribution parameters
        scalar = fit.m_to_mm()
        if 'aperture' in pars:   # constant value, no distribution
            self.fault_controls['aperture_mean'] = scalar*pars['aperture']
            self.fault_controls['aperture_std'] = 0.0
            aperture_array = scalar*pars['aperture']*np.ones(self.fault_controls['n_plates'])
            self.fault_controls['aperture_data'] = aperture_array
            self.fault_controls['aperture_approach'] = True
            self.fault_controls['vary_aperture'] = False  # defined by user so False
        elif 'aperture' in kwarg_p:  # array is given; individual value for each segment
            if len(kwarg_p['aperture']) != self.fault_controls['n_plates']:
                err_msg = ''.join([
                    'Size of provided aperture array for component {} is not ',
                    'equal to the number of fault segments']).format(self.name)
                logging.error(err_msg)
                raise ValueError(err_msg)

            self.fault_controls['aperture_data'] = scalar*np.array(kwarg_p['aperture'])
            self.fault_controls['aperture_approach'] = True
            self.fault_controls['vary_aperture'] = False  # defined by user so False

        else:  # variability in aperture is assumed
            self.fault_controls['aperture_mean'] = scalar*pars['aveAperture']
            self.fault_controls['aperture_std'] = scalar*pars['stdDevAperture']
            self.fault_controls['aperture_min'] = scalar*pars['minAperture']
            self.fault_controls['aperture_max'] = scalar*pars['maxAperture']
            # The following check comes from decipher.translate_aperture_parameters
            if self.fault_controls['aperture_std'] < deci.SMALL_DEV:
                self.fault_controls['vary_aperture'] = False
            else:
                self.fault_controls['vary_aperture'] = True
            self.fault_controls['aperture_approach'] = False

    def simulation_model(self, p, time_point=365.25, time_step=365.25,
                         pressure=None, CO2saturation=None, **kwargs):
        """
        :param p: input parameters of FaultFlow component model
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

        :Returns: out - dictionary of observations of FaultFlow model;
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
            rng = np.random.default_rng(conf.SEEDX)

        # Obtain the default values of the parameters from dictionary of
        # default parameters.
        actual_p = {k: v.value for k, v in self.default_pars.items()}

        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Check parameters provided by keyword arguments for initial time.
        msg = 'Keyword arguments provided to the fault component {}: {}'.format(
            self.name, kwargs)
        logging.debug(msg)

        if time_point == 0:
            # Update keyword parameters
            kwarg_p = self.check_keyword_parameters(**kwargs)
            # Update fault_controls attribute with values provided by user.
            # Due to the way component has to receive its parameters and nature
            # of user parameters processing in the seal flux standalone code
            # we have to update fault controls with user provided values
            self.update_fault_controls(actual_p, **kwarg_p)

            # Add extra accumulators after the number of fault segments became known
            num_extra_accumulators = self.fault_controls['n_plates']-1

            for i in range(num_extra_accumulators):
                self.add_accumulator(
                    'mass_CO2_aquifer_segm{}'.format(i+2), sim=0.0)
                self.add_accumulator(
                    'mass_brine_aquifer_segm{}'.format(i+2), sim=0.0)

            self.add_accumulator('mass_CO2_aquifer',
                                 sim=np.zeros(self.fault_controls['n_plates']))
            self.add_accumulator('mass_brine_aquifer',
                                 sim=np.zeros(self.fault_controls['n_plates']))

            # Get initial aperture data
            self.process_aperture_data(actual_p, **kwarg_p)

            # Define fault and update aperture data if needed
            self.define_fracto()

            # Compute interpolated conditions points
            self.fault_controls = pro.setup_transition_points(self.fault_controls)

            # Define fluid properties depending on case / profile type.
            self.fault_controls = pro.define_ave_interpolation_stages(
                False, self.fault_controls)

            # Define strike, aperture, length for fault plates
            self.fault_controls, self.fault = fcp.fault_create(
                self.fault_controls, self.fault, rng)

            self.fault_controls = pro.interpolate_initialize(self.fault_controls)

        # Initialize outputs.
        out = {}
        out['CO2_aquifer'] = np.zeros(self.fault_controls['n_plates'])
        out['brine_aquifer'] = np.zeros(self.fault_controls['n_plates'])
        out['mass_CO2_aquifer'] = np.zeros(self.fault_controls['n_plates'])
        out['mass_brine_aquifer'] = np.zeros(self.fault_controls['n_plates'])
        out['CO2_aquifer_total'] = 0.0
        out['brine_aquifer_total'] = 0.0
        out['mass_CO2_aquifer_total'] = 0.0
        out['mass_brine_aquifer_total'] = 0.0

        if time_point == 0.0:
            output_debug_message(out, time_point)
            return out

        # Get pressure/saturation data and check that it satisfies the predefined limits
        self.check_temporal_inputs(
            time_point, {'pressure': pressure, 'CO2saturation': CO2saturation})

        base_co2_pressure = pressure
        base_co2_saturation = CO2saturation

        co2_flow_step = np.zeros(self.fault_controls['n_plates'])
        brine_flow_step = np.zeros(self.fault_controls['n_plates'])

        self.fault_controls = pro.update_interpolate(
            time_point/365.25, self.fault_controls) # time is needed in years!!!

        for seg_number in range(self.fault_controls['n_plates']):

            # Get reservoir data for plate.
            co2_base_pressure = base_co2_pressure[seg_number]
            co2_base_saturation = base_co2_saturation[seg_number]

            # Compute flow rates for each plate - vol./sec.
            # Dissolved CO2 is included in CO2 flow calculation as
            # fluids despite assuming immiscible flow as theory.
            flow_rate_co2, flow_rate_brine = \
                self.fault[seg_number].compute_flow(
                    co2_base_pressure, co2_base_saturation, self.fault_controls)

            # Compute additional part of profile for complex models.
            if self.fault_controls['profile_type'] == 1:
                # Compute complex profile.
                flow_rate_co2, flow_rate_brine = \
                    self.fault[seg_number].complex_tail(
                        co2_base_pressure, co2_base_saturation,
                        self.fault_controls, flow_rate_co2, flow_rate_brine)

            if self.fault_controls['profile_type'] == 2:
                # Compute disjoint profile.
                flow_rate_co2, flow_rate_brine = \
                    self.fault[seg_number].complex_tail(
                        co2_base_pressure, co2_base_saturation,
                        self.fault_controls, flow_rate_co2, flow_rate_brine)

            co2_flow_step[seg_number] = flow_rate_co2        # m^3/sec
            brine_flow_step[seg_number] = flow_rate_brine    # m^3/sec

        # Convert flows vol/sec to kg/sec
        co2_flow_step = self.fault_controls['co2_density']*co2_flow_step
        brine_flow_step = self.fault_controls['brine_density']*brine_flow_step

        # Assign results to out dictionary
        for seg_number in range(self.fault_controls['n_plates']):
            out['CO2_aquifer'][seg_number] = co2_flow_step[seg_number]
            out['brine_aquifer'][seg_number] = brine_flow_step[seg_number]

            nm_end = '_segm{}'.format(seg_number+1)

            # Calculate individual masses for each fault segment
            out['mass_CO2_aquifer'][seg_number] = (
                self.accumulators['mass_CO2_aquifer'].sim[seg_number] +
                time_step*days_to_seconds()*out['CO2_aquifer'][seg_number])

            out['mass_brine_aquifer'][seg_number] = (
                self.accumulators['mass_brine_aquifer'].sim[seg_number] +
                time_step*days_to_seconds()*out['brine_aquifer'][seg_number])

            # Update accumulators
            self.accumulators['mass_CO2_aquifer'+nm_end].sim = (
                out['mass_CO2_aquifer'][seg_number])
            self.accumulators['mass_brine_aquifer'+nm_end].sim = (
                out['mass_brine_aquifer'][seg_number])

        # Calculate cumulative flow rate for all cells
        out['CO2_aquifer_total'] = np.sum(out['CO2_aquifer'])
        out['brine_aquifer_total'] = np.sum(out['brine_aquifer'])
        out['mass_CO2_aquifer_total'] = np.sum(out['mass_CO2_aquifer'])
        out['mass_brine_aquifer_total'] = np.sum(out['mass_brine_aquifer'])

        output_debug_message(out, time_point)

        return out


def output_debug_message(out, time_point):
    # Output useful debug message
    msg = ''.join([
        'Fault component outputs at time {}:\n CO2 rate: {}\n Brine rate: {}\n',
        'Total CO2 rate: {}\n Total brine rate: {}']).format(
            time_point, out['CO2_aquifer'], out['brine_aquifer'],
            out['CO2_aquifer_total'], out['brine_aquifer_total'])
    logging.debug(msg)


def days_to_seconds():
    """Convert days to seconds.

    Parameters
    ----------
    N/A (days)

    Returns
    -------
    value = (float) conversion factor in seconds (for Julian calendar)

    Notes
    -----
    1. 1 day = 24 hours * 60 minutes * 60 seconds =
    2. Value of days from Wikipedia.
    3. Computed in Excel.
    """
    value = 86400.0 # in seconds
    return value


def test_case1():
    """
    Runs the first example illustrating work of the Fault Flow component.

    Returns
    -------
    None.

    """
    # Define keyword arguments of the system model.
    # 46 time points
    time_array = 365.25 * np.array([
        0.00E+00, 8.49E-02, 1.64E-01, 2.49E-01, 3.31E-01, 4.16E-01,
        4.98E-01, 5.83E-01, 6.68E-01, 7.50E-01, 8.35E-01, 9.17E-01,
        1.00E+00, 2.00E+00, 4.00E+00, 5.00E+00, 6.00E+00, 7.00E+00,
        8.00E+00, 9.00E+00, 1.00E+01, 1.50E+01, 2.00E+01, 2.50E+01,
        3.00E+01, 3.50E+01, 4.00E+01, 4.50E+01, 5.00E+01, 5.50E+01,
        6.00E+01, 6.50E+01, 7.00E+01, 7.50E+01, 8.00E+01, 8.50E+01,
        9.00E+01, 1.10E+02, 1.30E+02, 1.40E+02, 1.50E+02, 1.60E+02,
        1.70E+02, 1.80E+02, 1.90E+02, 2.00E+02])
    sm_model_kwargs = {'time_array': time_array}     # time is given in days

    # Create system model.
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    pressure = np.array([
        3.035851266631769761e+07, 3.035851266631769761e+07, 3.035851266631769761e+07,
        3.035851266631769761e+07, 3.035851266631769761e+07, 3.035851266631769761e+07,
        3.035851266631769761e+07, 3.035851266631769761e+07, 3.035851266631769761e+07,
        3.035851266631769761e+07, 3.035851266631769761e+07, 3.035851266631769761e+07,
        3.035851266631769761e+07, 3.035858161389059946e+07, 3.036002951292150095e+07,
        3.036237373040009663e+07, 3.036664847991989926e+07, 3.037368113235569745e+07,
        3.038402326829069853e+07, 3.039843331102679670e+07, 3.041732494600139558e+07,
        3.059127967242810130e+07, 3.090643902815400064e+07, 3.134935823646359891e+07,
        3.188976931285379827e+07, 3.249181951941659674e+07, 3.312200033572259545e+07,
        3.375555958310069889e+07, 3.437615668677359819e+07, 3.497427688168109953e+07,
        3.554240488237709552e+07, 3.607599014905019850e+07, 3.657275741179469973e+07,
        3.703387877934990078e+07, 3.746025057016349584e+07, 3.785414805414119363e+07,
        3.821798439633449912e+07, 3.9500000e+07, 4.028282630954369903e+07,
        4.062384100510709733e+07, 4.091507555303669721e+07, 4.116363155334119499e+07,
        4.137640376331059635e+07, 4.155835640819369256e+07, 4.171431581809349358e+07,
        4.184800516194659472e+07]).reshape(46, 1)

    CO2saturation = np.array([
        2.000000000000000042e-03, 2.000000000000000042e-03, 2.000000000000000042e-03,
        2.000000000000000042e-03, 2.000000000000000042e-03, 2.000000000000000042e-03,
        2.000000000000000042e-03, 2.000000000000000042e-03, 2.000000000000000042e-03,
        2.000000000000000042e-03, 2.000000000000000042e-03, 2.000000000000000042e-03,
        2.000000000000000042e-03, 2.000000000000000042e-03, 1.999920000000000083e-03,
        1.999789999999999988e-03, 1.999550000000000112e-03, 1.999170000000000200e-03,
        1.998600000000000116e-03, 1.997810000000000037e-03, 1.996780000000000083e-03,
        1.987329999999999826e-03, 1.970539999999999792e-03, 1.947619999999999950e-03,
        1.920690000000000071e-03, 1.891960000000000039e-03, 1.863220000000000067e-03,
        1.835610000000000106e-03, 1.809749999999999944e-03, 1.785860000000000008e-03,
        1.764050000000000080e-03, 1.744310000000000001e-03, 1.726550000000000090e-03,
        1.710570000000000016e-03, 1.696199999999999921e-03, 1.683269999999999949e-03,
        1.671610000000000023e-03, 1.630E-03, 1.610080000000000061e-03,
        1.600620000000000081e-03, 1.592689999999999904e-03, 1.586020000000000051e-03,
        1.580390000000000093e-03, 1.575629999999999956e-03, 1.571579999999999939e-03,
        1.568149999999999926e-03]).reshape(46, 1)

    ffc = sm.add_component_model_object(FaultFlow(name='ffc', parent=sm))
    ffc.model_kwargs['relativeModel'] = 'BC'
    ffc.add_dynamic_kwarg('pressure', pressure)
    ffc.add_dynamic_kwarg('CO2saturation', CO2saturation)

    # Add gridded observations
    ffc.add_grid_obs('CO2_aquifer', constr_type='array',
                     output_dir=os.sep.join(['..', '..', 'output']))
    ffc.add_grid_obs('brine_aquifer', constr_type='array',
                     output_dir=os.sep.join(['..', '..', 'output']))
    # Add scalar observations
    ffc.add_obs('CO2_aquifer_total')
    ffc.add_obs('brine_aquifer_total')

    print('Starting simulation...')
    sm.forward()
    print('Simulation is finished.')

    plt.figure(1)

    # Collect scalar observations
    CO2_aquifer_total = sm.collect_observations_as_time_series(
        ffc, 'CO2_aquifer_total')
    print('CO2 aquifer total', CO2_aquifer_total, sep='\n')
    _, ax = plt.subplots()
    ax.plot(time_array/365.25, CO2_aquifer_total, '-g', label='Total CO2 rate')
    ax.legend()
    ax.set_xlabel('Time, [years]')
    ax.set_ylabel('Cumulative CO$_2$ Leakage\nRate to Aquifer, [kg/s]')
    ax.set_title('FaultFlow Component Test 1')

    plt.figure(2)

    brine_aquifer_total = sm.collect_observations_as_time_series(
        ffc, 'brine_aquifer_total')
    print('Brine aquifer total', brine_aquifer_total, sep='\n')
    _, ax = plt.subplots()
    ax.plot(time_array/365.25, brine_aquifer_total, '-b', label='Total brine rate')
    ax.legend()
    ax.set_xlabel('Time, [years]')
    ax.set_ylabel('Cumulative Brine Leakage\nRate to Aquifer, [kg/s]')
    ax.set_title('FaultFlow Component Test 1')


def test_case2():
    """
    Runs the first example illustrating work of the Fault Flow component.

    Returns
    -------
    None.

    """
    # Define keyword arguments of the system model.
    # 46 time points
    time_array = 365.25 * np.array([
        0.00E+00, 8.49E-02, 1.64E-01, 2.49E-01, 3.31E-01, 4.16E-01,
        4.98E-01, 5.83E-01, 6.68E-01, 7.50E-01, 8.35E-01, 9.17E-01,
        1.00E+00, 2.00E+00, 4.00E+00, 5.00E+00, 6.00E+00, 7.00E+00,
        8.00E+00, 9.00E+00, 1.00E+01, 1.50E+01, 2.00E+01, 2.50E+01,
        3.00E+01, 3.50E+01, 4.00E+01, 4.50E+01, 5.00E+01, 5.50E+01,
        6.00E+01, 6.50E+01, 7.00E+01, 7.50E+01, 8.00E+01, 8.50E+01,
        9.00E+01, 1.10E+02, 1.30E+02, 1.40E+02, 1.50E+02, 1.60E+02,
        1.70E+02, 1.80E+02, 1.90E+02, 2.00E+02])
    sm_model_kwargs = {'time_array': time_array}     # time is given in days

    # Create system model.
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    pressure = np.array([
        3.035851266631769761e+07, 3.035851266631769761e+07, 3.035851266631769761e+07,
        3.035851266631769761e+07, 3.035851266631769761e+07, 3.035851266631769761e+07,
        3.035851266631769761e+07, 3.035851266631769761e+07, 3.035851266631769761e+07,
        3.035851266631769761e+07, 3.035851266631769761e+07, 3.035851266631769761e+07,
        3.035851266631769761e+07, 3.035858161389059946e+07, 3.036002951292150095e+07,
        3.036237373040009663e+07, 3.036664847991989926e+07, 3.037368113235569745e+07,
        3.038402326829069853e+07, 3.039843331102679670e+07, 3.041732494600139558e+07,
        3.059127967242810130e+07, 3.090643902815400064e+07, 3.134935823646359891e+07,
        3.188976931285379827e+07, 3.249181951941659674e+07, 3.312200033572259545e+07,
        3.375555958310069889e+07, 3.437615668677359819e+07, 3.497427688168109953e+07,
        3.554240488237709552e+07, 3.607599014905019850e+07, 3.657275741179469973e+07,
        3.703387877934990078e+07, 3.746025057016349584e+07, 3.785414805414119363e+07,
        3.821798439633449912e+07, 3.9500000e+07, 4.028282630954369903e+07,
        4.062384100510709733e+07, 4.091507555303669721e+07, 4.116363155334119499e+07,
        4.137640376331059635e+07, 4.155835640819369256e+07, 4.171431581809349358e+07,
        4.184800516194659472e+07]).reshape(46, 1)

    CO2saturation = np.array([
        2.000000000000000042e-03, 2.000000000000000042e-03, 2.000000000000000042e-03,
        2.000000000000000042e-03, 2.000000000000000042e-03, 2.000000000000000042e-03,
        2.000000000000000042e-03, 2.000000000000000042e-03, 2.000000000000000042e-03,
        2.000000000000000042e-03, 2.000000000000000042e-03, 2.000000000000000042e-03,
        2.000000000000000042e-03, 2.000000000000000042e-03, 1.999920000000000083e-03,
        1.999789999999999988e-03, 1.999550000000000112e-03, 1.999170000000000200e-03,
        1.998600000000000116e-03, 1.997810000000000037e-03, 1.996780000000000083e-03,
        1.987329999999999826e-03, 1.970539999999999792e-03, 1.947619999999999950e-03,
        1.920690000000000071e-03, 1.891960000000000039e-03, 1.863220000000000067e-03,
        1.835610000000000106e-03, 1.809749999999999944e-03, 1.785860000000000008e-03,
        1.764050000000000080e-03, 1.744310000000000001e-03, 1.726550000000000090e-03,
        1.710570000000000016e-03, 1.696199999999999921e-03, 1.683269999999999949e-03,
        1.671610000000000023e-03, 1.630E-03, 1.610080000000000061e-03,
        1.600620000000000081e-03, 1.592689999999999904e-03, 1.586020000000000051e-03,
        1.580390000000000093e-03, 1.575629999999999956e-03, 1.571579999999999939e-03,
        1.568149999999999926e-03]).reshape(46, 1)

    ffc = sm.add_component_model_object(FaultFlow(name='ffc', parent=sm))
    ffc.add_dynamic_kwarg('pressure', pressure)
    ffc.add_dynamic_kwarg('CO2saturation', CO2saturation)

    # Add parameters of the FaultFlow component
    ffc.add_par('strike', value=95.4, vary=False)
    ffc.add_par('dip', value=90.0, vary=False)
    ffc.add_par('length', value=100.0, vary=False)
    ffc.add_par('xStart', value=388618.7, vary=False)
    ffc.add_par('yStart', value=3436341.1, vary=False)
    ffc.add_par('nSegments', value=1, vary=False)
    ffc.add_par('SGR', value=0.0, vary=False)
    ffc.add_par('stateVariable', value=1.0, vary=False)

    ffc.model_kwargs['aperture'] = [1.50E-05]

    ffc.add_par('aquiferDepth', value=500.0, vary=False)
    ffc.add_par('injectDepth', value=2884.31, vary=False)
    ffc.add_par('aquiferPressure', value=4.9e+06, vary=False)
    ffc.add_par('fieldPressure', value=2.85e+7, vary=False)
    ffc.add_par('injectPressure', value=3.0778e+7, vary=False)
    ffc.add_par('finalPressure', value=3.90e+7, vary=False)
    ffc.add_par('aquiferTemperature', value=30.0, vary=False)
    ffc.add_par('injectTemperature', value=98.9, vary=False)
    ffc.add_par('injectX', value=388505.9, vary=False)
    ffc.add_par('injectY', value=3434629.9, vary=False)

    ffc.add_par('salinity', value=50000.0, vary=False)
    ffc.add_par('CO2Density', value=430.0, vary=False)
    ffc.add_par('CO2Viscosity', value=3.72e-5, vary=False)
    ffc.add_par('brineDensity', value=988.00, vary=False)
    ffc.add_par('brineViscosity', value=4.36e-4, vary=False)
    ffc.add_par('CO2Solubility', value=2.0e-3, vary=False)

    ffc.add_par('brineResSaturation', value=0.15, vary=False)
    ffc.add_par('CO2ResSaturation', value=0.0, vary=False)
    ffc.model_kwargs['relativeModel'] = 'BC'
    ffc.add_par('permRatio', value=0.6, vary=False)
    ffc.add_par('entryPressure', value=1.0e+05, vary=False)
    ffc.add_par('lambda', value=2.5, vary=False)

    ffc.add_par('maxHorizontal', value=3.0e+07, vary=False)
    ffc.add_par('minHorizontal', value=2.0e+07, vary=False)
    ffc.add_par('maxTrend', value=55.0, vary=False)

    # Add gridded observations
    ffc.add_grid_obs('CO2_aquifer', constr_type='array',
                     output_dir=os.sep.join(['..', '..', 'output']))
    ffc.add_grid_obs('brine_aquifer', constr_type='array',
                     output_dir=os.sep.join(['..', '..', 'output']))

    # Add scalar observations
    ffc.add_obs('CO2_aquifer_total')
    ffc.add_obs('brine_aquifer_total')

    print('Starting simulation...')
    sm.forward()
    print('Simulation is finished.')

    plt.figure(3)

    # Collect scalar observations
    CO2_aquifer_total = sm.collect_observations_as_time_series(
        ffc, 'CO2_aquifer_total')
    print('CO2 aquifer total', CO2_aquifer_total, sep='\n')
    _, ax = plt.subplots()
    ax.plot(time_array/365.25, CO2_aquifer_total, '-g', label=r'Total CO$_2$ rate')
    ax.legend()
    ax.set_xlabel('Time, [years]')
    ax.set_ylabel('Cumulative CO$_2$ Leakage\nRate to Aquifer, [kg/s]')
    ax.set_title('FaultFlow Component Test 2')

    plt.figure(4)

    brine_aquifer_total = sm.collect_observations_as_time_series(
        ffc, 'brine_aquifer_total')
    print('Brine aquifer total', brine_aquifer_total, sep='\n')
    _, ax = plt.subplots()
    ax.plot(time_array/365.25, brine_aquifer_total, '-b', label='Total brine rate')
    ax.legend()
    ax.set_xlabel('Time, [years]')
    ax.set_ylabel('Cumulative Brine Leakage\nRate to Aquifer, [kg/s]')
    ax.set_title('FaultFlow Component Test 2')


if __name__ == "__main__":

    # Setup logging and constants.
    logging.basicConfig(level=logging.WARN)
    test_case = 1  # test case

    if test_case == 1:
        test_case1()
    elif test_case == 2:
        test_case2()
