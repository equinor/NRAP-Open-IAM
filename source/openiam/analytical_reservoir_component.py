# -*- coding: utf-8 -*-
from math import sqrt
import sys
import os
import logging
import numpy as np
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import SystemModel, ComponentModel
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

try:
    import components.reservoir.analytical.analytical_reservoir_ROM as resrom
except ImportError:
    print('\nERROR: Unable to load ROM for Analytical Reservoir component\n')
    sys.exit()

class AnalyticalReservoir(ComponentModel):
    """
    The Analytical Reservoir component model is a semi-analytical model for the
    reservoir. It is focused on flow across relatively large distances
    and does not take into account discrete features of the flow paths such as
    fractures, cracks, etc. The model is based on work of Nordbotten et al.,
    :cite:`C2011`. Further reading can be found in :cite:`N2011a`, :cite:`BaekEtAl2021`.

    In the NRAP-Open-IAM control file, the type name for the analytical reservoir
    component is ``AnalyticalReservoir``. The description of the component's parameters
    is provided below:

    * **logResPerm** [|log10| |m^2|] (-15.3 to -12) - logarithm of reservoir
      permeability (default: -13.69897)

    * **reservoirPorosity** [-] (0.1 to 0.3) - porosity of reservoir (default: 0.15)

    * **reservoirRadius** [|m|] (500 to 100,000) - distance between injection well
      and outer reservoir boundary (default: 500)

    * **brineDensity** [|kg/m^3|] (965 to 1195) - density of brine phase
      (default: 1045)

    * **CO2Density** [|kg/m^3|] (450 to 976) - density of |CO2| phase
      (default: 479)

    * **brineViscosity** [|Pa*s|] (2.3e-4 to 15.9e-4) - viscosity of brine phase
      (default: 2.535e-4)

    * **CO2Viscosity** [|Pa*s|] (0.455e-6 to 1.043e-4)  - viscosity of |CO2| phase
      (default: 3.95e-5)

    * **brineResSaturation** [-] (0 to 0.25) - residual saturation of brine phase
      (default: 0)

    * **brineCompressibility** [|Pa^-1|] (3.63e-12 to 2.31e-11) - brine compressibility
      (default: 4.5e-12 = 3.1e-8 1/psi)

    * **injRate** [|m^3/s|] (0.0024 to 3.776) - |CO2| injection rate (default: 0.01)

    * **numberOfShaleLayers** [-] (3 to 30) - number of shale layers in the
      system (default: 3); *linked to Stratigraphy*. The shale units must be
      separated by an aquifer.

    * **shaleThickness** [|m|] (1 to 1600) - thickness of shale layers (default:
      250); *linked to Stratigraphy*. Thickness of shale layer 1, for example,
      can be defined by **shale1Thickness**; otherwise, shale layers for which
      the thickness is not defined will be assigned a default thickness.

    * **aquiferThickness** [|m|] (1 to 1600) - thickness of aquifers (default: 100);
      *linked to Stratigraphy*. Thickness of aquifer 1, for example, can be defined
      by **aquifer1Thickness**; otherwise, aquifers for which the thickness
      is not defined will be assigned a default thickness.

    * **reservoirThickness** [|m|] (15 to 500) - thickness of reservoir (default: 50);
      *linked to Stratigraphy*

    * **datumPressure** [|Pa|] (80,000 to 300,000) - pressure at the top of the
      system (default: 101,325); *linked to Stratigraphy*.

    Possible observations from the Analytical Reservoir component are:

    * **pressure** [|Pa|] - pressure at top of the reservoir at the user
      defined location(s)

    * **pressureAve** [|Pa|] - pressure vertically averaged at the user
      defined location(s)

    * **CO2saturation** [-] - |CO2| saturation vertically averaged at the user
      defined location(s)

    * **mass_CO2_reservoir** [|kg|] - (injected total) mass of the |CO2|
      in the reservoir.

    """
    def __init__(self, name, parent, injX=0., injY=0., locX=100., locY=0.):

        """
        Constructor method of AnalyticalReservoir class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param injX: x-coordinate of the injector location
        :type injX: float

        :param injY: y-coordinate of the injector location
        :type injY: float

        :param locX: x-coordinate of the location for pressure and CO2 saturation
            to be calculated
        :type locX: float

        :param locY: y-coordinate of the location for pressure and CO2 saturation
            to be calculated
        :type locY: float

        :returns: AnalyticalReservoir class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25, 'time_step': 365.25}

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'AnalyticalReservoir'

        # Set default parameters of the component model
        self.add_default_par('logResPerm', value=-13.69897)
        self.add_default_par('reservoirPorosity', value=0.15)
        self.add_default_par('brineDensity', value=1045.0)
        self.add_default_par('CO2Density', value=479.0)
        self.add_default_par('brineViscosity', value=2.535e-4)
        self.add_default_par('CO2Viscosity', value=3.95e-5)
        self.add_default_par('brineResSaturation', value=0.0)
        self.add_default_par('brineCompressibility', value=4.5e-12)
        self.add_default_par('reservoirRadius', value=500)
        self.add_default_par('reservoirThickness', value=50.0)
        self.add_default_par('numberOfShaleLayers', value=3)
        self.add_default_par('shaleThickness', value=250.0)
        self.add_default_par('aquiferThickness', value=100.0)
        self.add_default_par('datumPressure', value=101325.0)
        self.add_default_par('injRate', value=0.01)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['logResPerm'] = [-15.3, -12.0]
        self.pars_bounds['reservoirPorosity'] = [0.10, 0.30]
        self.pars_bounds['brineDensity'] = [965.0, 1195.0]
        self.pars_bounds['CO2Density'] = [450.0, 976.0]
        self.pars_bounds['brineViscosity'] = [2.3e-4, 15.9e-4]
        self.pars_bounds['CO2Viscosity'] = [0.455e-6, 1.043e-4]
        self.pars_bounds['brineResSaturation'] = [0.0, 0.25]
        self.pars_bounds['brineCompressibility'] = [3.63e-12, 2.31e-11]
        self.pars_bounds['reservoirRadius'] = [500, 100000]
        self.pars_bounds['reservoirThickness'] = [15, 500.0]
        self.pars_bounds['numberOfShaleLayers'] = [3, 30]
        self.pars_bounds['shaleThickness'] = [1.0, 3000.0]
        self.pars_bounds['aquiferThickness'] = [1.0, 1600.0]
        self.pars_bounds['datumPressure'] = [8.0e+4, 3.0e+5]
        self.pars_bounds['injRate'] = [0.0024, 3.776]

        # Specify default locations of injection and leaking wells
        self.injX = injX
        self.injY = injY
        self.locX = locX
        self.locY = locY

        # Setup default observations of the component
        self.default_obs = {'pressure': 101325.0,
                            'pressureAve': 101325.0,
                            'CO2saturation': 0.0,
                            'mass_CO2_reservoir': 0.0}

        # Add accumulator observation mass_CO2_reservoir used within the component
        self.add_accumulator('mass_CO2_reservoir', sim=0.0)
        self.add_accumulator('volume_CO2_reservoir', sim=0.0)
        # self.prevCO2Saturation = 0.0

        debug_msg = 'AnalyticalReservoir component created with name {}'.format(self.name)
        logging.debug(debug_msg)

    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        debug_msg = 'Input parameters of component {} are {}'.format(self.name, p)
        logging.debug(debug_msg)

        for key, val in p.items():
            warn_msg = ''.join([
                'Parameter {} of AnalyticalReservoir component {} ',
                'is out of boundaries.']).format(key, self.name)
            if key.startswith('shale') and key.endswith('Thickness'):
                if (val < self.pars_bounds['shaleThickness'][0]) or (
                        val > self.pars_bounds['shaleThickness'][1]):
                    logging.warning(warn_msg)
            elif key.startswith('aquifer') and key.endswith('Thickness'):
                if (val < self.pars_bounds['aquiferThickness'][0]) or (
                        val > self.pars_bounds['aquiferThickness'][1]):
                    logging.warning(warn_msg)
            elif key in self.pars_bounds:
                if (val < self.pars_bounds[key][0]) or (
                        val > self.pars_bounds[key][1]):
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of AnalyticalReservoir component {}.']).format(key, self.name)
                logging.warning(warn_msg)

    def simulation_model(self, p, time_point=365.25, time_step=365.25,
                         injX=None, injY=None, locX=None, locY=None):
        """
        Return pressure and CO2 saturation at the bottom of leaking well.
        Note that the names of parameters contained in the input dictionary
        coincide with those defined in this modules docstring.

        :param p: input parameters of analytical reservoir model
        :type p: dict

        :param time_point: time point (in days) for which the pressure and
            saturation are to be calculated; by default, its value is 365.25
            (1 year in days)
        :type time_point: float

        :param time_step: difference between the current and previous
            time points in days; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param mass_CO2_reservoir: mass of |CO2| accumulated in reservoir at
            time specified by time_point parameter; by default, its value is
            0.0; this parameter supplies value for an accumulator
            (needed for the next time point) if solution has to be calculated
            for several time steps.
        :type pressure: float

        :param injX: x-coordinate of the injector location
        :type injX: float

        :param injY: y-coordinate of the injector location
        :type injY: float

        :param locX: x-coordinate of the location for pressure and CO2 saturation
            to be calculated
        :type locX: float

        :param locY: y-coordinate of the location for pressure and CO2 saturation
            to be calculated
        :type locY: float

        :returns: out - dictionary of observations of analytical reservoir
            component model;
            keys: ['pressure','CO2saturation','mass_CO2_reservoir','pressureAve']

        """
        debug_msg = ''.join([
            'AnalyticalReservoir component {} model call input: ',
            'parameters {}, location xy ({}, {})']).format(
                self.name, p, self.locX, self.locY)
        logging.debug(debug_msg)

        # Set coordinates of injector and leak point/point of interest
        if injX is None:
            injX = self.injX
        if injY is None:
            injY = self.injY
        if locX is None:
            locX = self.locX
        if locY is None:
            locY = self.locY

        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}

        # Update default values of parameters with the provided ones
        actual_p.update(p)
        # Create instance of Parameters class
        inputParameters = resrom.Parameters()

        # Set up number of shale layers
        nSL = int(actual_p['numberOfShaleLayers'])
        inputParameters.numberOfShaleLayers = nSL

        # Initialize shale and aquifer thickness arrays
        inputParameters.shaleThickness = actual_p['shaleThickness']*np.ones(nSL)
        inputParameters.aquiferThickness = actual_p['aquiferThickness']*np.ones((nSL-1))

        # Set up shale and aquifer thickness
        for i in range(nSL):
            nm = 'shale{}Thickness'.format(i+1)
            if nm in p:
                inputParameters.shaleThickness[i] = p[nm]
        for i in range(nSL-1):
            nm = 'aquifer{}Thickness'.format(i+1)
            if nm in p:
                inputParameters.aquiferThickness[i] = p[nm]

        inputParameters.reservoirThickness = actual_p['reservoirThickness']
        if inputParameters.reservoirThickness < 0:
            warn_msg = 'Reservoir thickness {0} is less than expected lowest value of 0.'.format(
                inputParameters.reservoirThickness)
            logging.warning(warn_msg)

        # Set up reservoir permeability
        inputParameters.reservoirPermeability = 10**actual_p['logResPerm']

        # Set up reservoir porosity
        inputParameters.reservoirPorosity = actual_p['reservoirPorosity']

        # Set up land surface pressure
        inputParameters.datumPressure = actual_p['datumPressure']

        # Set up brine and CO2 density
        inputParameters.brineDensity = actual_p['brineDensity']
        inputParameters.CO2Density = actual_p['CO2Density']

        # Set up brine and CO2 viscosity
        inputParameters.brineViscosity = actual_p['brineViscosity']
        inputParameters.CO2Viscosity = actual_p['CO2Viscosity']

        # Set up residual saturation
        inputParameters.brineResidualSaturation = actual_p['brineResSaturation']

        # Rate at the injection well
        inputParameters.rate = actual_p['injRate']

        # Reservoir Outer Radius
        inputParameters.resOuterRadius = actual_p['reservoirRadius']

        # Set up distance between injection and leaking well
        inputParameters.distance = sqrt((locX-injX)**2+(locY-injY)**2)

        if inputParameters.distance > inputParameters.resOuterRadius*2:
            warn_msg = ''.join([
                'Distance {}m between injection and observation wells ',
                'is greater than the model reservoir size {}m.']).format(
                    inputParameters.distance, inputParameters.resOuterRadius*2)
            logging.warning(warn_msg)

        # Brine Compressibility
        inputParameters.brineCompressibility = actual_p['brineCompressibility']

        # Parameters from the class attributes
        inputParameters.timeStep = time_step       # in days
        inputParameters.timePoint = time_point     # in days
        inputParameters.prevCO2Volume = self.accumulators['volume_CO2_reservoir'].sim

        # Create solution object with defined input parameters
        sol = resrom.Solution(inputParameters)

        # Initiate dictionary of leakage rates
        out = dict()
        if time_point == 0.0:   #  return initial state of the reservoir component
            sol.setup_initial_conditions()
            # Save observations
            out['pressure'] = sol.initialTopPressure
            out['CO2saturation'] = sol.interface/sol.thicknessH
            out['pressureAve'] = (sol.initialTopPressure+sol.initialPressure)/2
        else:
            # Find solution corresponding to the inputParameters
            sol.find()
            # Save observations
            out['pressure'] = sol.pressureTop
            out['CO2saturation'] = sol.saturation
            out['pressureAve'] = sol.pressure

        out['mass_CO2_reservoir'] = sol.CO2Volume*inputParameters.CO2Density
        self.accumulators['mass_CO2_reservoir'].sim = sol.CO2Volume*inputParameters.CO2Density
        self.accumulators['volume_CO2_reservoir'].sim = sol.CO2Volume

        debug_msg = 'Analytical Reservoir component {} model call output: {}'.format(
            self.name, out)
        logging.debug(debug_msg)

        # Return dictionary of outputs
        return out

    # Attributes for system connections
    needs_locXY = True
    needs_injXY = True
    system_params = ['numberOfShaleLayers',
                     'shale1Thickness',
                     'shale2Thickness',
                     'shale3Thickness',
                     'aquifer1Thickness',
                     'aquifer2Thickness',
                     'reservoirThickness',
                     'datumPressure']

    def reset(self):
        pass

if __name__ == "__main__":
    __spec__ = None
    logging.basicConfig(level=logging.WARNING)
    # Define keyword arguments of the system model
    isUQAnalysisOn = 0 # 0: one forward run, 1: stochastic runs
    num_samples = 10   # for stochastic runs
    to_save = False    # figure save

    num_years = 40
    time_array = 365.25*np.arange(0, num_years+1)

    # delta_time = 5
    # time_array = np.arange(0,365.25*num_years,delta_time)

    sm_model_kwargs = {'time_array': time_array} # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    res = sm.add_component_model_object(AnalyticalReservoir(
        name='res', parent=sm, injX=100., injY=0., locX=0., locY=0.))

    # Add parameters
    res.add_par('reservoirThickness', value=30.0, vary=False)
    res.add_par('shaleThickness', value=970.0, vary=False)
    res.add_par('aquiferThickness', value=30.0, vary=False)

    # Add parameters of reservoir component model
    # Ebigbo 2007 setting (validated)
    res.add_par('logResPerm', value=-13.69897, vary=False)
    res.add_par('reservoirPorosity', value=0.15, vary=False)
    res.add_par('datumPressure', value=101325.0, vary=False)
    res.add_par('brineDensity', value=1045.0, vary=False)
    res.add_par('CO2Density', value=479.0, vary=False)
    res.add_par('brineViscosity', value=2.535e-4, vary=False)
    res.add_par('CO2Viscosity', value=3.95e-5, vary=False)
    res.add_par('brineResSaturation', value=0.0, vary=False)
    res.add_par('injRate', value=0.01, vary=False)
    res.add_par('reservoirRadius', value=500, vary=False)
    res.add_par('brineCompressibility', value=1.e-11, vary=False)

    # UQ Analysis
    if isUQAnalysisOn:
        res.add_par('injRate', min=8.87/479.0*0.2, max=8.87/479.0*5, value=8.87/479.0)

    # Add observations of reservoir component model
    # When add_obs method is called, system model automatically
    # (if no other options of the method is used) adds a list of observations
    # with names res.obsnm_0, res.obsnm_1,.... The final index is determined by the
    # number of time points in system model time_array attribute.
    # For more information, see the docstring of add_obs of ComponentModel class.
    res.add_obs('pressure')
    res.add_obs('pressureAve')
    res.add_obs('CO2saturation')
    res.add_obs('mass_CO2_reservoir')

    if isUQAnalysisOn != 1:

        # Run system model using current values of its parameters
        sm.forward()

        print('------------------------------------------------------------------')
        print('                  Forward method illustration ')
        print('------------------------------------------------------------------')

        # Plot pressure and saturation
        # Since the observations at the particular time points are different variables,
        # method collect_observations_as_time_series creates lists of
        # values of observations belonging to a given component (e.g. cw) and having the same
        # common name (e.g. 'pressure', 'CO2saturation', etc) but differing in indices.
        # More details are given in the docstring and documentation to the method
        # collect_observations_as_time_series of SystemModel class.

        pressureX = sm.collect_observations_as_time_series(res, 'pressure')
        CO2Sat = sm.collect_observations_as_time_series(res, 'CO2saturation')
        pressureAve = sm.collect_observations_as_time_series(res, 'pressureAve')
        CO2mass = sm.collect_observations_as_time_series(res, 'mass_CO2_reservoir')

        # Plot results
        # Setup plot parameters
        line_width = 1
        xy_label_size = 14
        title_size = 18
        xlims = [0, num_years]
        fig = plt.figure(figsize=(14, 6))

        # pressure
        ax = fig.add_subplot(121)
        plt.plot(time_array/365.25, pressureX/1.0e+6,
                 color="maroon", linewidth=line_width, label="Ptop")
        plt.plot(time_array/365.25, pressureAve/1.0e+6,
                 color="blue", linewidth=line_width, label="Pave")
        plt.xlabel('Time, t (years)', fontsize=xy_label_size)
        plt.ylabel('Pressure, P (MPa)', fontsize=xy_label_size)
        plt.title('Pressure at leaking well', fontsize=title_size)
        plt.tight_layout()
        plt.tick_params(labelsize=12)
        plt.xlim(xlims)
        plt.legend()
        # plt.ylim([20, 100])
        ax.get_yaxis().set_label_coords(-0.081, 0.5)

        # saturation
        ax = fig.add_subplot(122)
        plt.plot(time_array/365.25, CO2Sat,
                 color="maroon", linewidth=line_width)
        plt.xlabel('Time, t (years)', fontsize=xy_label_size)
        plt.ylabel('CO2 Saturation', fontsize=xy_label_size)
        plt.title('CO2 Saturation at leaking well', fontsize=title_size)
        plt.tight_layout()
        plt.tick_params(labelsize=12)
        plt.xlim(xlims)
        plt.ylim([0, 1.2])
        ax.get_yaxis().set_label_coords(-0.08, 0.5)

        print('Forward simulation is done. Check the figure. ')

        if to_save:
            plt.savefig('ReservoirComponentCombinedPlot1.png', dpi=300)
            print('"ReservoirComponentCombinedPlot1.png" successfully saved.')

    # UQ Analysis
    if isUQAnalysisOn:
        print('------------------------------------------------------------------')
        print('                          UQ illustration ')
        print('------------------------------------------------------------------')

        import random

        # Draw Latin hypercube samples of parameter values
        seed = random.randint(500, 1100)
        s = sm.lhs(siz=num_samples, seed=seed) # create sample set

        # Run model using values in samples for parameter values
        ncpus = 5
        results = s.run(cpus=ncpus, verbose=False)

        # Extract results from stochastic simulations
        pressure = np.ones((num_samples, len(time_array)))
        CO2_saturation = np.ones((num_samples, len(time_array)))
        mass_CO2_reservoir = np.ones((num_samples, len(time_array)))

        for ind in range(len(time_array)):
            pressure[:, ind] = s.recarray['res.pressure_{}'.format(ind)]
            CO2_saturation[:, ind] = s.recarray['res.CO2saturation_{}'.format(ind)]
            mass_CO2_reservoir[:, ind] = s.recarray['res.mass_CO2_reservoir_{}'.format(ind)]

        # Plot results
        # Setup plot parameters
        line_width = 1
        xy_label_size = 14
        title_size = 18
        xlims = [0, num_years]
        fig = plt.figure(figsize=(18, 4))
        ax = fig.add_subplot(131)
        for j in range(num_samples):
            plt.plot(time_array/365.25, pressure[j]/1.0e+6,
                     color="maroon", linewidth=line_width)
        plt.xlabel('Time, t (years)', fontsize=xy_label_size)
        plt.ylabel('Pressure, P (MPa)', fontsize=xy_label_size)
        plt.title('Pressure: leaking well', fontsize=title_size)
        plt.tight_layout()
        plt.tick_params(labelsize=12)
        plt.xlim(xlims)
        plt.ylim([20, 50])
        ax.get_yaxis().set_label_coords(-0.12, 0.5)

        ax = fig.add_subplot(132)
        for j in range(num_samples):
            plt.plot(time_array/365.25, CO2_saturation[j],
                     color="green", linewidth=line_width)
        plt.xlabel('Time, t (years)', fontsize=xy_label_size)
        plt.ylabel('Saturation, S (-)', fontsize=xy_label_size)
        plt.title(r'CO$_2$ saturation: leaking well', fontsize=title_size)
        plt.tight_layout()
        plt.tick_params(labelsize=12)
        plt.xlim(xlims)
        plt.ylim([0.0, 1.1])
        ax.get_yaxis().set_label_coords(-0.12, 0.5)

        ax = fig.add_subplot(133)
        for j in range(num_samples):
            plt.plot(time_array/365.25, 479.0*mass_CO2_reservoir[j],
                     color="darkblue", linewidth=line_width)
        plt.xlabel('Time, t (years)', fontsize=xy_label_size)
        plt.ylabel(r'Mass of CO$_2$, m (kg)', fontsize=xy_label_size)
        plt.title(r'Amount of CO$_2$ in reservoir', fontsize=title_size)
        plt.tight_layout()
        plt.tick_params(labelsize=12)
        plt.xlim(xlims)
        plt.ylim([0.0, 6.0e+09])
        ax.get_yaxis().set_label_coords(-0.12, 0.5)

        print('UQ run is done. Check the figure. ')
        if to_save:
            plt.savefig('ReservoirComponentCombinedPlot2.png', dpi=300)
            print('"ReservoirComponentCombinedPlot2.png" successfully saved.')
