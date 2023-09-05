# -*- coding: utf-8 -*-
import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import SystemModel, ComponentModel
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

from openiam.cfi.commons import process_parameters, process_dynamic_inputs

try:
    import components.wellbore.open as owmodel
    import components.wellbore.open.open_wellbore_ROM as owrom
except ImportError:
    print('\nERROR: Unable to load ROM for Open Wellbore component\n')
    sys.exit()

# Assumed gravitational acceleration for critical pressure calculation, m/(s^2)
GRAV_ACCEL = 9.8

# Assumed density of water for critical pressure calculation, kg/(m^3).
WATER_DENSITY = 1000

class OpenWellbore(ComponentModel):
    """
    The Open Wellbore model is a lookup table reduced order model based on the
    drift-flux approach, see :cite:`RN1899`. This model treats the leakage of |CO2|
    up an open wellbore or up an open (i.e., uncemented) casing/tubing. The lookup
    table is populated using T2Well/ECO2N Ver. 1.0 :cite:`RN1900`, which treats
    the non-isothermal flow of |CO2| and brine up an open wellbore, allows for the
    phase transition of |CO2| from supercritical to gaseous, with Joule-Thompson
    cooling, and considers exsolution of |CO2| from the brine phase.

    By default, when used within the control file interface the Open Wellbore
    is connected to the upper aquifer (e.g., aquifer 2 if there are 2 aquifers
    in the system). For user-defined scenarios the ``LeakTo`` keyword can be used
    to specify either the name of the aquifer (e.g., *aquifer1*) |CO2| leaks to
    or *atmosphere* for leakage to the atmosphere. The default value is
    *aquifer#* where # is an index of the uppermost aquifer.

    The Open Wellbore component can be used to calculate leakage rates following
    positive change in reservoir pressure or only from changes in reservoir
    pressure above a critical pressure. To use the latter approach, the
    argument ``crit_pressure_approach`` should be set to *True* for the setup
    of the component. Here is an example of this setup in a script application:

        OpenWellbore(name='ow', parent=sm, crit_pressure_approach=True)

    To set ``crit_pressure_approach`` to *True* in the control file interface,
    the ``Controls`` section in the .yaml entry for the Open Wellbore should be
    included with an additional entry ``critPressureApproach: True`` indented
    beneath ``Controls``. For an example setup, see control file example 31a.

    If ``crit_pressure_approach`` is set to *True*, the default approach is for
    critical pressure to be calculated as:

        Pcrit = (rho_w * g * d_aq) + (rho_br * g * (d_res - d_aq)),

    where rho_w and rho_br are the densities of water and brine (by default,
    1000 |kg/m^3|) defined by the brineDensity parameter, respectively,
    g is gravitational acceleration(9.8 |m/s^2|), d_aq is the depth to the
    bottom of the aquifer impacted by leakage (m) (defined by the wellTop
    parameter value; if wellTop is 0 m, then the atmosphere receives leakage),
    and d_res is the depth to the top of the reservoir (m). Higher brine densities
    generally produce lower leakage rates.

    Instead of calculating critical pressure in this manner, one can enforce a
    particular critical pressure (the **critPressure** parameter) by setting the
    argument ``enforce_crit_pressure`` to *True* for the setup of the component.
    Here is an example in a script setup:

        OpenWellbore(name='ow', parent=sm, crit_pressure_approach=True,
                     enforce_crit_pressure=True)

    To set ``enforce_crit_pressure`` to *True* in the control file interface,
    file interface, the ``Controls`` section in the .yaml entry for the Open
    Wellbore should be included with additional entry ``enforceCritPressure: True``
    indented beneath ``Controls``. If enforce_crit_pressure is not set to *True*,
    then the **critPressure** parameter will not be used.

    When a critical pressure is used, flow through an Open Wellbore can still
    occur at pressures beneath the critical pressure if a |CO2| plume is present
    in the reservoir at the base of the well (i.e., buoyancy effects from the
    |CO2|).

    Component model input definitions:

    * **logReservoirTransmissivity** [|log10| |m^3|] (-11.27 to -8.40) - reservoir
      transmissivity (default: -9.83)

    * **logAquiferTransmissivity** [|log10| |m^3|] (-11.27 to -8.40) - reservoir
      transmissivity (default: -9.83)

    * **brineSalinity** [-] (0 to 0.2) - brine salinity (mass fraction) (default: 0.1)

    * **brineDensity** [|kg/m^3|] (900 to 1200) - brine density (default: 1012)

    * **wellRadius** [|m|] (0.025 to 0.25) - radius of the wellbore (default: 0.05)

    * **wellTop** [|m|] (0 to 500) - depth of well top (default: 500); *linked to
      Stratigraphy*. Note that this parameter represents how far leakage can
      extend from the reservoir along the open wellbore. For example, the cement
      used to plug the well may be of poor quality or damaged, but the damage
      may not allow leakage to reach the surface at 0 m. When using control files,
      wellTop can be set to the bottom depth of an aquifer by entering
      'aquifer#Depth,' where # is the aquifer number. If the aquifer is too deep,
      however, the limits for this parameter would still be enforced. In control
      files, if the 'LeakTo' entry is provided as the name of the aquifer
      receiving leakage from the open wellbore (e.g., 'LeakTo: aquifer2')
      then the wellTop parameter will automatically be set to the bottom depth
      of the corresponding aquifer. If 'LeakTo' is set to 'atmosphere',
      wellTop will automatically be set to 0.

    * **critPressure** [|Pa|] (1.0e+5 to 9.0e+7) - pressure above which the model
      initiates leakage rates calculations. Default value of this parameter is
      not defined: either the user provides it through component setup or the
      value is calculated based on the value of **brineDensity** parameter.

    * **reservoirDepth** [|m|] (1000 to 4000) - depth of reservoir (well base)
      (default: 2000); *linked to Stratigraphy*. Note that if 'shale1Depth'
      is entered for this parameter in a control file, the parameter will be
      set to the bottom depth of shale 1 (which is also the top of the reservoir).

    The possible outputs from the Open Wellbore component are leakage rates
    of |CO2| and brine to aquifer and atmosphere. The names of the
    observations are of the form:

    * **CO2_aquifer** and **CO2_atm** [|kg/s|] - |CO2| leakage rates

    * **brine_aquifer** and **brine_atm** [|kg/s|] - brine leakage rates.

    """
    def __init__(self, name, parent, header_file_dir=None,
                 crit_pressure_approach=False, enforce_crit_pressure=False):
        """
        Constructor method of OpenWellbore class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param header_file_dir: path to the directory that contains
            the look up tables
        :type workdir: str

        :param crit_pressure_approach: flag indicating whether change in pressure
            exceeding critical pressure (calculated or provided by user) will be used
            to initiate the leakage versus non-zero pressure change. By default,
            non-zero pressure change is used.
        :type crit_pressure_approach: boolean

        :param enforce_crit_pressure: flag indicating whether user provided
        value of critical pressure will be used to model a leakage. Value of False
        means that the value of brineDensity will be used to estimate critical
        pressure value.
        :type enforce_crit_pressure: boolean

        :returns: OpenWellbore class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25}   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'OpenWellbore'

        # Define output dictionary labels
        self.output_labels = ['CO2_atm', 'CO2_aquifer', 'brine_atm', 'brine_aquifer']

        # Setup default observations of the component
        self.default_obs = {obs_nm: 0.0 for obs_nm in self.output_labels}

        # Set default parameters of the component model
        self.initPressure = 0.
        self.add_default_par('wellRadius', value=0.05)
        self.add_default_par('wellTop', value=500.0)
        self.add_default_par('reservoirDepth', value=2000.0)
        self.add_default_par('logReservoirTransmissivity', value=-9.83)
        self.add_default_par('logAquiferTransmissivity', value=-9.83)
        self.add_default_par('brineSalinity', value=0.1)
        self.add_default_par('brineDensity', value=1012.0)

        # Define dictionary of model parameters boundaries
        # Boundaries of some parameters are defined by the ROM
        self.pars_bounds = dict()
        self.pars_bounds['logReservoirTransmissivity'] = [-11.27, -8.4]
        self.pars_bounds['logAquiferTransmissivity'] = [-11.27, -8.4]
        self.pars_bounds['brineSalinity'] = [0, 0.2]
        self.pars_bounds['brineDensity'] = [900, 1200]
        self.pars_bounds['wellRadius'] = [0.025, 0.25]
        self.pars_bounds['wellTop'] = [0.0, 500.0]
        self.pars_bounds['critPressure'] = [1.0e+5, 9.0e+7]
        self.pars_bounds['reservoirDepth'] = [1000.0, 4000.0]

        # Instantiate solution
        if header_file_dir is None:
            header_file_dir = owmodel.__path__[0]
        self.sol = owrom.Solution(header_file_dir)

        # Save whether critical pressure is to be used
        self.use_crit_pressure = crit_pressure_approach

        # If use_crit_pressure is True, enforce_crit_pressure sets whether the
        # critical pressure is calculated automatically (False, default) or setup
        # by the critPressure parameter (True)
        self.enforce_crit_pressure = enforce_crit_pressure

        # Instantiate attribute keeping the index of aquifer layer to which
        # the leakage is estimated
        self.leak_layer = 0

        self.default_out = {}

        debug_msg = 'OpenWellbore component created with name {}'.format(name)
        logging.debug(debug_msg)

    def simulation_model(self, p, time_point=365.25,
                         pressure=0.0, CO2saturation=0.0):
        """
        Return CO2 and brine leakage rates corresponding to the provided input.

        :param p: dictionary of input parameters
        :type p: dict

        :param time_point: current time point in days (at which the output is
            to be provided); by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param pressure: pressure at the bottom of leaking well
            at the current time point, in Pa; by default, its value is 0.0
        :type pressure: float

        :param saturation: CO2 saturation at the bottom of leaking well;
            by default, its value is 0.0
        :type saturation: float

        :returns: dictionary of observations (leakage rates, etc.) of Open
            wellbore model; keys:
            ['CO2_rate','brine_rate']

        """
        # Return the initial state of the model if requested
        if time_point == 0.0:
            self.initPressure = pressure
            # For time point 0 all outputs of the model are zeros
            out = self.default_out.copy()
            out.update(dict(list(zip(self.output_labels, len(self.output_labels)*[0.0]))))
            return out

        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Set parameters
        wellTop = actual_p['wellTop']
        reservoirDepth = actual_p['reservoirDepth']
        logReservoirTransmissivity = actual_p['logReservoirTransmissivity']
        logAquiferTransmissivity = actual_p['logAquiferTransmissivity']
        brineSalinity = actual_p['brineSalinity']
        brineDensity = actual_p['brineDensity']

        # Calculate critical pressure delta (EPA method)
        if self.use_crit_pressure:
            if self.enforce_crit_pressure:
                try:
                    critP = actual_p['critPressure']
                except KeyError:
                    warn_msg = "".join([
                        "Parameter 'critPressure' was not setup for Open Wellbore ",
                        "component {} but arguments 'crit_pressure_approach' ",
                        "and 'enforce_crit_pressure' were set to True. ",
                        "critPressure will be calculated based on brineDensity ",
                        "parameter"]).format(self.name)
                    logging.warning(warn_msg)
                    critP = (wellTop * GRAV_ACCEL * WATER_DENSITY) + (
                        brineDensity * GRAV_ACCEL * (reservoirDepth - wellTop))
            else:
                critP = (wellTop * GRAV_ACCEL * WATER_DENSITY) + (
                    brineDensity * GRAV_ACCEL * (reservoirDepth - wellTop))

            # Calculate pressure change above critical
            deltaP = pressure - critP

            # If deltaP is negative - set to zero
            deltaP = max(deltaP, 0)

        else:  # otherwise calculate change from initial pressure
            deltaP = pressure-self.initPressure

        # Define array of input parameters
        inputArray = np.array([wellTop, reservoirDepth, logReservoirTransmissivity,
                               logAquiferTransmissivity, deltaP, CO2saturation,
                               brineSalinity])

        # Find solution
        self.sol.find(inputArray)

        # Rescale solution by constant dependent on the wellbore radius
        wellRadius = actual_p['wellRadius']
        originalRadius = 0.05
        solScalar = (wellRadius/originalRadius)**2
        CO2LeakageRates = solScalar*self.sol.CO2LeakageRates
        brineLeakageRates = solScalar*self.sol.brineLeakageRates

        # Return dictionary of leakage rates
        out = self.default_out.copy()
        out.update(dict(list(zip(
            self.output_labels, np.concatenate([CO2LeakageRates, brineLeakageRates])))))
        out['CO2_aquifer{}'.format(self.leak_layer)] = out['CO2_aquifer']
        out['brine_aquifer{}'.format(self.leak_layer)] = out['brine_aquifer']
        return out

    def connect_with_system(self, component_data, name2obj_dict, *args, **kwargs):
        """
        Code to add open wellbore to system model for control file interface.

        :param component_data: Dictionary of component data including 'Connection'
        :type component_data: dict

        :param name2obj_dict: Dictionary to map component names to objects
        :type name2obj: dict

        :returns: None
        """
        # Process parameters data if any
        process_parameters(self, component_data, name2obj_dict)

        # Process dynamic inputs if any
        process_dynamic_inputs(self, component_data)

        # Check whether critical pressure is to be used for leakage criteria
        if 'Controls' in component_data:
            self.use_crit_pressure = component_data['Controls'].get(
                'critPressureApproach', False)

            self.enforce_crit_pressure = component_data['Controls'].get(
                'enforceCritPressure', False)

        # Determine number of shale layers in the stratigraphy
        strata = name2obj_dict['strata']
        if 'numberOfShaleLayers' in strata.deterministic_pars:
            num_shale_layers = strata.deterministic_pars['numberOfShaleLayers'].value
        elif 'numberOfShaleLayers' in strata.default_pars:
            num_shale_layers = strata.default_pars['numberOfShaleLayers'].value
        else:
            num_shale_layers = 3

        # Determine to which aquifer the leakage is simulated
        if 'LeakTo' not in component_data:
            err_msg = ''.join(["Required argument 'LeakTo' is missing ",
                               "in the setup of Open Wellbore component. ",
                               "It should be setup to 'atmosphere' or ",
                               "to 'aquifer#' where # is an index of aquifer ",
                               "of interest."])
            logging.error(err_msg)
            raise KeyError(err_msg)
        else:
            leak_to = component_data['LeakTo'].lower()

        # Determine index of aquifer to which the leakage is simulated
        if 'aquifer' in leak_to:
            self.leak_layer = int(leak_to[7:])
        else:
            self.leak_layer = num_shale_layers - 1

        # Make model connections
        if 'Connection' in component_data:
            connection = None
            try:
                connection = name2obj_dict[component_data['Connection']]
            except KeyError:
                pass
            system_inputs = ['pressure', 'CO2saturation']
            for sinput in system_inputs:
                connection.add_obs_to_be_linked(sinput)
                self.add_kwarg_linked_to_obs(sinput, connection.linkobs[sinput])

        # Consider an option when reservoirDepth and wellTop are already added
        # as referring to spatially varying stratigraphy
        if 'reservoirDepth' not in self.pars and \
                'reservoirDepth' not in self.deterministic_pars:
            res_depth_expr = ' + '.join(['{}.shale{}Thickness'.format(
                strata.name, ind) for ind in range(1, num_shale_layers+1)])+' + '+\
                    ' + '.join(['{}.aquifer{}Thickness'.format(
                        strata.name, ind) for ind in range(1, num_shale_layers)])

            # Depth to the top of reservoir (usually)
            self.add_composite_par('reservoirDepth', res_depth_expr)

        if leak_to == 'atmosphere':
            self.add_par('wellTop', value=0.0, vary=False)
        else:
            if 'wellTop' not in self.pars and \
                    'wellTop' not in self.deterministic_pars:
                well_top_expr = ' + '.join([
                    '{nm}.shale{ind1}Thickness + {nm}.aquifer{ind2}Thickness'.format(
                        nm=strata.name, ind1=ind+1, ind2=ind) for ind in range(
                            self.leak_layer, num_shale_layers)])
                self.add_composite_par('wellTop', well_top_expr)

        for il in range(1, num_shale_layers):
            aq = 'aquifer{}'.format(il)
            self.default_out['CO2_' + aq] = 0.0
            self.default_out['brine_' + aq] = 0.0


if __name__ == "__main__":
    try:
        from openiam import SimpleReservoir
    except ImportError as err:
        print('Unable to load IAM class module: '+str(err))

    test_case = 1
    if test_case == 1:
        crit_pressure_approach = False
    else:
        crit_pressure_approach = True

    # Create system model
    num_years = 10.
    delta_t = 0.1
    time_array = 365.25*np.arange(0.0, num_years+delta_t, delta_t)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))

    # Add parameters of reservoir component model
    sres.add_par('numberOfShaleLayers', value=3, vary=False)
    sres.add_par('shale1Thickness', min=500.0, max=550., value=525.0)
    sres.add_par('shale2Thickness', min=510.0, max=550., value=525.0)
    # Shale 3 has a fixed thickness of 1.0 m
    sres.add_par('shale3Thickness', min=1.0, max=2.0, value=1.5)
    # Aquifer 1 (thief zone has a fixed thickness of 1.0)
    sres.add_par('aquifer1Thickness', value=1.0, vary=False)
    # Aquifer 2 (shallow aquifer) has a fixed thickness of 200
    sres.add_par('aquifer2Thickness', value=200.0, vary=False)
    # Reservoir has a fixed thickness of 51.2
    sres.add_par('reservoirThickness', value=51.2, vary=False)
    # Injection Rate
    sres.add_par('injRate', value=0.0025, vary=False)

    # Add observations of reservoir component model
    # When add_obs method is called, system model automatically
    # (if no other options of the method is used) adds a list of observations
    # with name sres.obsnm_0, sres.obsnm_1,.... The final index is determined by the
    # number of time points in system model time_array attribute.
    # For more information, see the docstring of add_obs of ComponentModel class.
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')
    sres.add_obs_to_be_linked('pressure')
    sres.add_obs_to_be_linked('CO2saturation')

    # Add open wellbore component
    ow = sm.add_component_model_object(
        OpenWellbore(name='ow', parent=sm,
                     crit_pressure_approach=crit_pressure_approach))

    # Add parameters of open wellbore component
    ow.add_par('logReservoirTransmissivity', min=-11.0, max=-9.0, value=-10.0)
    ow.add_par('logAquiferTransmissivity', min=-11.0, max=-9.0, value=-10.0)
    ow.add_par('brineSalinity', min=0.0, max=0.0, value=0.0)
    # The next parameter is not used if crit_pressure_approach is False (default)
    ow.add_par('brineDensity', value=1020, vary=False)

    # Add keyword arguments of the open wellbore component model
    ow.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
    ow.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

    # Add composite parameter of open wellbore component
    ow.add_composite_par('reservoirDepth',
                         expr='+'.join(['sres.shale1Thickness',
                                        'sres.shale2Thickness',
                                        'sres.shale3Thickness',
                                        'sres.aquifer1Thickness',
                                        'sres.aquifer2Thickness']))
    ow.add_composite_par(
        'wellTop', expr='sres.shale3Thickness + sres.aquifer2Thickness')

    # Add observations of the open wellbore component
    ow.add_obs('CO2_aquifer')
    ow.add_obs('CO2_atm')
    ow.add_obs('brine_aquifer')
    ow.add_obs('brine_atm')

    sm.forward()
    # Since the observations at the particular time points are different variables,
    # method collect_observations_as_time_series creates lists of
    # values of observations belonging to a given component (e.g. cw) and having the same
    # common name (e.g. 'CO2_aquifer1', etc) but differing in indices.
    # More details are given in the docstring and documentation to the method
    # collect_observations_as_time_series of SystemModel class.
    pressure = sm.collect_observations_as_time_series(sres, 'pressure')
    CO2saturation = sm.collect_observations_as_time_series(sres, 'CO2saturation')
    print('------------------------------------------------------------------')
    print('Pressure', pressure, sep='\n')
    print('------------------------------------------------------------------')
    print('CO2 saturation', CO2saturation, sep='\n')

    CO2leakrates_aq1 = sm.collect_observations_as_time_series(ow, 'CO2_aquifer')
    CO2leakrates_atm = sm.collect_observations_as_time_series(ow, 'CO2_atm')
    print('------------------------------------------------------------------')
    print('CO2 leakage rates to aquifer', CO2leakrates_aq1, sep='\n')
    print('------------------------------------------------------------------')
    print('CO2 leakage rates to atmosphere', CO2leakrates_atm, sep='\n')

    brine_leakrates_aq1 = sm.collect_observations_as_time_series(ow, 'brine_aquifer')
    brine_leakrates_atm = sm.collect_observations_as_time_series(ow, 'brine_atm')
    print('------------------------------------------------------------------')
    print('Brine leakage rates to aquifer', brine_leakrates_aq1, sep='\n')
    print('------------------------------------------------------------------')
    print('Brine leakage rates to atmosphere', brine_leakrates_atm, sep='\n')
    print('------------------------------------------------------------------')

    reservoirDepth = (sres.pars['shale1Thickness'].value \
                      + sres.pars['shale2Thickness'].value \
                      + sres.pars['shale3Thickness'].value \
                      + sres.deterministic_pars['aquifer1Thickness'].value \
                      + sres.deterministic_pars['aquifer2Thickness'].value)

    wellTop =  sres.pars['shale3Thickness'].value \
        + sres.deterministic_pars['aquifer2Thickness'].value
    brineDensity = ow.deterministic_pars['brineDensity'].value

    # Calculate the critical pressure and add to plot
    critical_pres = wellTop*9.8*1000 + brineDensity*9.8*(reservoirDepth-wellTop)

    print('------------------------------------------------------------------')
    print('Critical Pressure', critical_pres, sep='\n')

    plt.figure(1)
    plt.plot(sm.time_array/365.25, pressure, color='#000066', linewidth=1)
    plt.hlines(critical_pres, 0, sm.time_array[-1]/365.25, linestyle='--', color='k')
    plt.xlabel('Time, t (years)')
    plt.ylabel('Reservoir Pressure, P (Pa)')
    plt.title('Reservoir Pressure vs. time')
    plt.show()

    plt.figure(2)
    plt.plot(sm.time_array/365.25, CO2leakrates_aq1, color='#000066', linewidth=1)
    plt.xlabel('Time, t (years)')
    plt.ylabel('Leakage rates, q (kg/s)')
    plt.title(r'Leakage of CO$_2$: Shallow aquifer')
    plt.show()

    plt.figure(3)
    plt.plot(sm.time_array/365.25, brine_leakrates_aq1, color='#000066', linewidth=1)
    plt.xlabel('Time, t (years)')
    plt.ylabel('Leakage rates, q (kg/s)')
    plt.title('Leakage of brine: Shallow aquifer')
    plt.show()
