# -*- coding: utf-8 -*-
import logging
import sys
import os
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import SystemModel, ComponentModel
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

try:
    import components.wellbore.multisegmented.multisegmented_wellbore_ROM as mswrom
except ImportError:
    print('\nERROR: Unable to load ROM for Multisegmented Wellbore component\n')
    sys.exit()

class MultisegmentedWellbore(ComponentModel):
    """
    The Multisegmented Wellbore component estimates the leakage rates of brine and
    |CO2| along wells in the presence of overlying aquifers or thief zones.
    The model is based on work of Nordbotten et al.,
    :cite:`N2009`. Further reading can be found in :cite:`N2011a`.

    The model is focused on flow across relatively large distances and does not take
    into account discrete features of the flow paths such as fractures, cracks,
    etc. It assumes that leakage is occurring in the annulus between the
    outside of the casing and borehole. This area is assigned an “effective”
    permeability of the flow path. The permeability is applied over a length
    along the well that corresponds to the thickness of a shale formation.
    Each well is characterized by an effective permeability assigned to each
    segment of the well that crosses an individual formation. For example, if
    a well crosses N permeable formations, then it is characterized by
    N different permeability values. The model utilizes the one-dimensional
    multiphase version of Darcy’s law to represent flow along a leaky well.

    In the NRAP-Open-IAM control file, the type name for the Multisegmented
    Wellbore component is ``MultisegmentedWellbore``. The description
    of the component's parameters are provided below. Names of the component
    parameters coincide with those used by ``model`` method of the
    ``MultisegmentedWellbore`` class.

    * **logWellPerm** [|log10| |m^2|] (-101 to -9) - logarithm of well permeability
      along shale layer (default: -13). Logarithm of well permeability along shale 3,
      for example, can be defined by **logWell3Perm**. Permeability of the well
      along the shale layers not defined by user will be assigned a default value.

    * **logAquPerm** [|log10| |m^2|] (-17 to -9) - logarithm of aquifer
      permeability (default: -12). Logarithm of aquifer 1 permeability, for example,
      can be defined by **logAqu1Perm**. Aquifer permeability not defined by user will
      be assigned a default value.

    * **brineDensity** [|kg/m^3|] (900 to 1500) - density of brine phase
      (default: 1000)

    * **CO2Density** [|kg/m^3|] (100 to 1000) - density of |CO2| phase
      (default: 479)

    * **brineViscosity** [|Pa*s|] (1.0e-4 to 5.0e-3) - viscosity of brine phase
      (default: 2.535e-3)

    * **CO2Viscosity** [|Pa*s|] (1.0e-6 to 1.0e-4)  - viscosity of |CO2| phase
      (default: 3.95e-5)

    * **aquBrineResSaturation** [-] (0 to 0.99) - residual saturation of brine phase
      in each aquifer (default: 0.0). For example, the residual brine saturation
      of aquifer2 can be defined by **aqu2BrineResSaturation**; otherwise, aquifer
      layers for which the residual brine saturation is not defined will be
      assigned a default value.

    * **compressibility** [|Pa^-1|] (1.0e-13 to 1.0e-8) - compressibility of brine
      and |CO2| phases (assumed to be the same for both phases) (default: 5.1e-11)

    * **wellRadius** [|m|] (0.01 to 0.5) - radius of leaking well (default: 0.05)

    * **numberOfShaleLayers** [-] (3 to 30) - number of shale layers in the
      system (default: 3); *linked to Stratigraphy*. The shale units must be
      separated by an aquifer.

    * **shaleThickness** [|m|] (1 to 3000) - thickness of shale layers (default:
      250); *linked to Stratigraphy*. Thickness of shale layer 1, for example,
      can be defined by **shale1Thickness**; otherwise, shale layers for which
      the thickness is not defined will be assigned a default thickness.

    * **aquiferThickness** [|m|] (1 to 1600) - thickness of aquifers (default: 100);
      *linked to Stratigraphy*. Thickness of aquifer 1, for example, can be defined
      by **aquifer1Thickness**; otherwise, aquifers for which the thickness
      is not defined will be assigned a default thickness.

    * **reservoirThickness** [|m|] (1 to 1600) - thickness of reservoir (default: 30);
      *linked to Stratigraphy*.

    * **datumPressure** [|Pa|] (80,000 to 300,000) - pressure at the top of the
      system (default: 101,325); *linked to Stratigraphy*

    The possible outputs from the Multisegmented Wellbore component are
    leakage rates of |CO2| and brine to each of the aquifers in the system and
    atmosphere. The names of the observations are of the form:

    * **CO2_aquifer1**, **CO2_aquifer2**,..., **CO2_atm** [|kg/s|] -
      |CO2| leakage rates

    * **brine_aquifer1**, **brine_aquifer2**,..., **brine_atm** [|kg/s|] -
      brine leakage rates

    * **mass_CO2_aquifer1**, **mass_CO2_aquifer2**,..., **mass_CO2_aquiferN** [|kg|]
      - mass of the |CO2| leaked into the aquifer.

    """
    def __init__(self, name, parent):
        """
        Constructor method of MultisegmentedWellbore class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :returns: MultisegmentedWellbore class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25, 'time_step': 365.25}   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'MultisegmentedWellbore'

        # Set default parameters of the component model
        self.add_default_par('numberOfShaleLayers', value=3)
        self.add_default_par('shaleThickness', value=250.0)
        self.add_default_par('aquiferThickness', value=100.0)
        self.add_default_par('reservoirThickness', value=30.0)
        self.add_default_par('logWellPerm', value=-13.0)
        self.add_default_par('logAquPerm', value=-12.0)
        self.add_default_par('datumPressure', value=101325.0)
        self.add_default_par('brineDensity', value=1000.0)
        self.add_default_par('CO2Density', value=479.0)
        self.add_default_par('brineViscosity', value=2.535e-3)
        self.add_default_par('CO2Viscosity', value=3.95e-5)
        self.add_default_par('aquBrineResSaturation', value=0.0)
        self.add_default_par('compressibility', value=5.1e-11)
        self.add_default_par('wellRadius', value=0.05)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['numberOfShaleLayers'] = [3, 30]
        self.pars_bounds['shaleThickness'] = [1.0, 3000.0]
        self.pars_bounds['aquiferThickness'] = [1.0, 1600.0]
        self.pars_bounds['reservoirThickness'] = [1.0, 1600.0]
        self.pars_bounds['logWellPerm'] = [-101.0, -9.0]
        self.pars_bounds['logAquPerm'] = [-17.0, -9.0]
        self.pars_bounds['datumPressure'] = [8.0e+4, 3.0e+5]
        self.pars_bounds['brineDensity'] = [900.0, 1500.0]
        self.pars_bounds['CO2Density'] = [100.0, 1000.0]
        self.pars_bounds['brineViscosity'] = [1.0e-4, 5.0e-3]
        self.pars_bounds['CO2Viscosity'] = [1.0e-6, 1.0e-4]
        self.pars_bounds['aquBrineResSaturation'] = [0.0, 0.99]
        self.pars_bounds['compressibility'] = [1.0e-13, 1.0e-8]
        self.pars_bounds['wellRadius'] = [0.01, 0.5]

        # By default, the smallest number of aquifers the system can have is 2,
        # so we add two accumulators by default. Extra will be added as needed,
        # once the system knows how many there are
        for i in range(2):
            self.add_accumulator('mass_CO2_aquifer'+str(i+1), sim=0.0)
            self.add_accumulator('volume_CO2_aquifer'+str(i+1), sim=0.0)
            self.add_accumulator('CO2_saturation_aquifer'+str(i+1), sim=0.0)
        self.num_accumulators = 2

        # Define output dictionary labels
        self.output_labels = ['CO2_aquifer1', 'CO2_aquifer2', 'CO2_atm',
                              'brine_aquifer1', 'brine_aquifer2', 'brine_atm',
                              'mass_CO2_aquifer1', 'mass_CO2_aquifer2']

        # Setup default observations of the component
        self.default_obs = {obs_nm: 0.0 for obs_nm in self.output_labels}

        debug_msg = 'MultisegmentedWellbore component created with name {}'.format(name)
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
                'Parameter {} of MultisegmentedWellbore component {} ',
                'is out of boundaries.']).format(key, self.name)
            if key.startswith('shale') and key.endswith('Thickness'):
                if (val < self.pars_bounds['shaleThickness'][0]) or (
                        val > self.pars_bounds['shaleThickness'][1]):
                    logging.warning(warn_msg)
                continue
            if key.startswith('aquifer') and key.endswith('Thickness'):
                if (val < self.pars_bounds['aquiferThickness'][0]) or (
                        val > self.pars_bounds['aquiferThickness'][1]):
                    logging.warning(warn_msg)
                continue
            if key.startswith('logWell') and key.endswith('Perm'):
                if (val < self.pars_bounds['logWellPerm'][0]) or (
                        val > self.pars_bounds['logWellPerm'][1]):
                    logging.warning(warn_msg)
                continue
            if key.startswith('logAqu') and key.endswith('Perm'):
                if ((val < self.pars_bounds['logAquPerm'][0]) or
                        (val > self.pars_bounds['logAquPerm'][1])):
                    logging.warning(warn_msg)
                continue
            if key.startswith('aqu') and key.endswith('BrineResSaturation'):
                if ((val < self.pars_bounds['aquBrineResSaturation'][0]) or
                        (val > self.pars_bounds['aquBrineResSaturation'][1])):
                    logging.warning(warn_msg)
                continue
            if key in self.pars_bounds:
                if ((val < self.pars_bounds[key][0]) or
                        (val > self.pars_bounds[key][1])):
                    logging.warning(warn_msg)

    def simulation_model(self, p, time_point=365.25, time_step=365.25,
                         pressure=0.0, CO2saturation=0.0):
        """
        Return |CO2| and brine leakage rates corresponding to the provided input.
        Note that the names of parameters contained in the input dictionary
        coincide with those defined in this module's docstring.

        :param p: input parameters of multisegmented wellbore model
        :type p: dict

        :param pressure: pressure at the bottom of leaking well, in Pa;
            by default, its value is 0.0
        :type pressure: float

        :param CO2saturation: saturation of |CO2| phase at the bottom
            of leaking well; by default, its value is 0.0
        :type CO2saturation: float

        :param time_point: time point in days at which the leakage rates are
            to be calculated; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param time_step: difference between the current and previous
            time points in days; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :returns: out - dictionary of observations of multisegmented wellbore
            model; keys:
            ['CO2_aquifer1','CO2_aquifer2',...,'CO2_atm',
            'brine_aquifer1','brine_aquifer2',...,'brine_atm',
            'mass_CO2_aquifer1','mass_CO2_aquifer2',...]
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)
        inputParameters = mswrom.Parameters()

        nSL = int(actual_p['numberOfShaleLayers'])
        inputParameters.numberOfShaleLayers = nSL

        # Add extra accumulators after the number of aquifers became known
        num_extra_accumulators = nSL-1 - self.num_accumulators
        if num_extra_accumulators != 0: # if component has no enough accumulators
            for i in range(num_extra_accumulators):
                self.add_accumulator('mass_CO2_aquifer{}'.format(i+3), sim=0.0)
                self.add_accumulator('volume_CO2_aquifer{}'.format(i+3), sim=0.0)
                self.add_accumulator('CO2_saturation_aquifer{}'.format(i+3), sim=0.0)

            self.num_accumulators = nSL - 1

        # Check whether the initial state of the system is requested
        # Then if yes, return the state without proceeding further
        if time_point == 0.0:
            # Create dictionary of leakage rates
            out = dict()
            for i in range(nSL-1):
                out['CO2_aquifer{}'.format(i+1)] = 0.0
                out['brine_aquifer{}'.format(i+1)] = 0.0
                out['mass_CO2_aquifer{}'.format(i+1)] = 0.0
                self.accumulators['mass_CO2_aquifer{}'.format(i+1)].sim = 0.0
                self.accumulators['volume_CO2_aquifer{}'.format(i+1)].sim = 0.0
                self.accumulators['CO2_saturation_aquifer{}'.format(i+1)].sim = 0.0
            out['CO2_atm'] = 0.0
            out['brine_atm'] = 0.0
            return out

        inputParameters.shaleThickness = actual_p['shaleThickness']*np.ones(nSL)
        inputParameters.shalePermeability = 10**actual_p['logWellPerm']*np.ones(nSL)
        inputParameters.aquiferThickness = actual_p['aquiferThickness']*np.ones((nSL-1))
        inputParameters.aquiferPermeability = 10**actual_p['logAquPerm']*np.ones((nSL-1))
        if 'brineResSaturation' in actual_p:
            inputParameters.aquBrineResSaturation = \
                actual_p['brineResSaturation']*np.ones((nSL))
        else:
            inputParameters.aquBrineResSaturation = \
                actual_p['aquBrineResSaturation']*np.ones((nSL))

        # Set up shale, aquifer and reservoir thickness parameters
        for i in range(nSL):
            nm = 'shale{}Thickness'.format(i+1)
            if nm in p:
                inputParameters.shaleThickness[i] = p[nm]
        for i in range(nSL-1):
            nm = 'aquifer{}Thickness'.format(i+1)
            if nm in p:
                inputParameters.aquiferThickness[i] = p[nm]
        inputParameters.reservoirThickness = actual_p['reservoirThickness']

        # Set up permeability of the well segment along shale
        for i in range(nSL):
            nm = 'logWell{}Perm'.format(i+1)
            if nm in p:
                inputParameters.shalePermeability[i] = 10**p[nm]

        # Set up permeability of the wellsegment along aquifer
        for i in range(nSL-1):
            nm = 'logAqu{}Perm'.format(i+1)
            if nm in p:
                inputParameters.aquiferPermeability[i] = 10**p[nm]

        # Set up residual saturation of aquifers
        for i in range(nSL):
            nm = 'aqu{}BrineResSaturation'.format(i) # aqu0BrineResSaturation is reservoir
            if nm in p:
                inputParameters.aquBrineResSaturation[i] = p[nm]

        # Set up land surface pressure
        inputParameters.datumPressure = actual_p['datumPressure']

        # Set up brine and CO2 density
        inputParameters.brineDensity = actual_p['brineDensity']
        inputParameters.CO2Density = actual_p['CO2Density']

        # Set up brine and CO2 viscosity
        inputParameters.brineViscosity = actual_p['brineViscosity']
        inputParameters.CO2Viscosity = actual_p['CO2Viscosity']

        # Set up residual saturation and compressibility
        inputParameters.compressibility = actual_p['compressibility']

        # Set up well radius parameter
        inputParameters.wellRadius = actual_p['wellRadius']
        inputParameters.flowArea = np.pi*inputParameters.wellRadius**2

        # Get parameters from keyword arguments of the 'model' method
        inputParameters.timeStep = time_step         # in days
        inputParameters.timePoint = time_point       # in days
        inputParameters.pressure = pressure

        inputParameters.prevCO2Volume = np.zeros(nSL-1)
        inputParameters.CO2saturation = np.zeros(nSL)
        inputParameters.CO2saturation[0] = CO2saturation

        for i in range(nSL-1):
            inputParameters.prevCO2Volume[i] = (
                self.accumulators['volume_CO2_aquifer{}'.format(i+1)].sim)
            inputParameters.CO2saturation[i+1] = (
                self.accumulators['CO2_saturation_aquifer{}'.format(i+1)].sim)

        # Create solution object with defined input parameters
        sol = mswrom.Solution(inputParameters)

        # Find solution corresponding to the inputParameters
        sol.find()

        # Create dictionary of leakage rates
        out = dict()
        for i in range(nSL-1):
            out['CO2_aquifer{}'.format(i+1)] = sol.CO2LeakageRates[i]
            out['brine_aquifer{}'.format(i+1)] = sol.brineLeakageRates[i]
            # Convert mass in cubic meters to kg
            out['mass_CO2_aquifer{}'.format(i+1)] = \
                sol.CO2Volume[i]*inputParameters.CO2Density
            # Keep mass in cubic meters in accumulators
            self.accumulators['mass_CO2_aquifer{}'.format(i+1)].sim = \
                sol.CO2Volume[i]*inputParameters.CO2Density
            self.accumulators['volume_CO2_aquifer{}'.format(i+1)].sim = sol.CO2Volume[i]
            self.accumulators['CO2_saturation_aquifer{}'.format(i+1)].sim = \
                sol.CO2SaturationAq[i]
        out['CO2_atm'] = sol.CO2LeakageRates[inputParameters.numberOfShaleLayers-1]
        out['brine_atm'] = sol.brineLeakageRates[inputParameters.numberOfShaleLayers-1]

        # Return dictionary of outputs
        return out

    def reset(self):
        return

    # Attributes for system connections
    system_inputs = ['pressure',
                     'CO2saturation']
    system_params = ['numberOfShaleLayers',
                     'shale1Thickness',
                     'shale2Thickness',
                     'shale3Thickness',
                     'aquifer1Thickness',
                     'aquifer2Thickness',
                     'reservoirThickness',
                     'datumPressure']

def read_data(filename):
    """
    Routine used for reading data files.

    Read data from file and create numpy.array if file exists.
    """
    # Check whether the file with given name exists
    if os.path.isfile(filename):
        data = np.genfromtxt(filename)
        return data

    return None


if __name__ == "__main__":
    try:
        from openiam import AnalyticalReservoir
    except ImportError as err:
        print('Unable to load IAM class module: {}'.format(err))
    import matplotlib.pyplot as plt
    __spec__ = None

    logging.basicConfig(level=logging.WARNING)
    # Define keyword arguments of the system model
    isUQAnalysisOn = 0 # 0: one forward run, 1: stochastic runs
    isPlottingOn = 1
    to_save_png = False

    num_years = 40
    time_array = 365.25*np.arange(0, num_years+1)

    # delta_time = 5
    # time_array = np.arange(0,365.25*num_years,delta_time)

    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    res = sm.add_component_model_object(AnalyticalReservoir(
        name='res', parent=sm, injX=100., injY=0., locX=0., locY=0.))

    # Add parameters of reservoir component model

    res.add_par('injRate', min=1e-5, max=1e+2, value=0.0185, vary=False)
    res.add_par('reservoirRadius', min=500.0, max=1000000, value=500.0, vary=False)
    res.add_par('reservoirThickness', min=10.0, max=150., value=50.0, vary=False)
    res.add_par('logResPerm', min=-16.0, max=-9.0, value=-13.69897, vary=False)
    res.add_par('reservoirPorosity', min=0.01, max=0.5, value=0.15, vary=False)

    res.add_par('shaleThickness', value=970.0, vary=False)
    res.add_par('aquiferThickness', value=30.0, vary=False)

    res.add_par('brineDensity', min=900.0, max=1300., value=1045.0, vary=False)
    res.add_par('CO2Density', min=200.0, max=900., value=479.0, vary=False)
    res.add_par('brineViscosity', min=1.0e-5, max=1.0e-3, value=2.535e-4, vary=False)
    res.add_par('CO2Viscosity', min=1.0e-5, max=1.0e-4, value=3.95e-5, vary=False)
    res.add_par('brineResSaturation', min=0.0, max=0.99, value=0.0, vary=False)
    res.add_par('brineCompressibility', min=1e-9, max=1e-13, value=1.e-11, vary=False)

    res.add_par('datumPressure', value=101325.0, vary=False)

    # Add observations of reservoir component model
    res.add_obs_to_be_linked('pressure')
    res.add_obs_to_be_linked('CO2saturation')

    res.add_obs('pressure')
    res.add_obs('CO2saturation')
    res.add_obs('mass_CO2_reservoir')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms',
                                                              parent=sm))
    ms.add_par('wellRadius', min=0.01, max=0.2, value=0.15)

    ms.add_par('numberOfShaleLayers', value=5, vary=False)
    ms.add_par('shale1Thickness', min=1.0, max=200, value=100.0)
    ms.add_par('shale2Thickness', min=1.0, max=200, value=100.0)
    ms.add_par('shale3Thickness', min=1.0, max=200, value=100.0)
    ms.add_par('shale4Thickness', min=1.0, max=200, value=100.0)
    ms.add_par('shale5Thickness', min=1.0, max=2500, value=2450.0)

    ms.add_par('logWell1Perm', min=-16, max=-9, value=-12)
    ms.add_par('logWell2Perm', min=-16, max=-9, value=-12)
    ms.add_par('logWell3Perm', min=-16, max=-9, value=-12)
    ms.add_par('logWell4Perm', min=-16, max=-9, value=-12)
    ms.add_par('logWell5Perm', min=-101, max=-9, value=-100)

    ms.add_par('aquifer1Thickness', min=1.0, max=200., value=30.0)
    ms.add_par('aquifer2Thickness', min=1.0, max=200., value=30.0)
    ms.add_par('aquifer3Thickness', min=1.0, max=200., value=30.0)
    ms.add_par('aquifer4Thickness', min=1.0, max=200., value=30.0)

    ms.add_par('logAqu1Perm', min=-16, max=-9, value=-12)
    ms.add_par('logAqu2Perm', min=-16, max=-9, value=-12)
    ms.add_par('logAqu3Perm', min=-16, max=-9, value=-12)
    ms.add_par('logAqu4Perm', min=-16, max=-9, value=-12)

    ms.add_par('aqu1BrineResSaturation', min=0.0, max=0.99, value=0.18)
    ms.add_par('aqu2BrineResSaturation', min=0.0, max=0.99, value=0.4)
    ms.add_par('aqu3BrineResSaturation', min=0.0, max=0.99, value=0.5)
    ms.add_par('aqu4BrineResSaturation', min=0.0, max=0.99, value=0.0)

    # Add linked parameters: common to both components
    ms.add_par_linked_to_par('reservoirThickness',
                             res.default_pars['reservoirThickness'])
    ms.add_par_linked_to_par('aqu0BrineResSaturation',
                             res.default_pars['brineResSaturation']) # reservoir
    ms.add_par_linked_to_par('datumPressure',
                             res.default_pars['datumPressure'])

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', res.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', res.linkobs['CO2saturation'])

    # Add observations of multisegmented wellbore component model
    # When add_obs method is called, system model automatically
    # (if no other options of the method is used) adds a list of observations
    # with names ms.obsnm_0, ms.obsnm_1,.... The final index is determined by the
    # number of time points in system model time_array attribute.
    # For more information, see the docstring of add_obs of ComponentModel class.
    ms.add_obs('CO2_aquifer1')
    ms.add_obs('CO2_aquifer2')
    ms.add_obs('CO2_aquifer3')
    ms.add_obs('CO2_aquifer4')
    ms.add_obs('CO2_atm')

    ms.add_obs('mass_CO2_aquifer1')
    ms.add_obs('mass_CO2_aquifer2')
    ms.add_obs('mass_CO2_aquifer3')
    ms.add_obs('mass_CO2_aquifer4')

    ms.add_obs('brine_aquifer1')
    ms.add_obs('brine_aquifer2')
    ms.add_obs('brine_aquifer3')
    ms.add_obs('brine_aquifer4')

    if isUQAnalysisOn == 0:
        # Run system model using current values of its parameters
        sm.forward(save_output=True)  # system model is run deterministically

        print('------------------------------------------------------------------')
        print('                  Forward method illustration ')
        print('------------------------------------------------------------------')
        # Since the observations at the particular time points are different variables,
        # method collect_observations_as_time_series creates lists of
        # values of observations belonging to a given component (e.g. cw) and having the same
        # common name (e.g. 'CO2_aquifer1', etc) but differing in indices.
        # More details are given in the docstring and documentation to the method
        # collect_observations_as_time_series of SystemModel class.

        pressure = sm.collect_observations_as_time_series(res, 'pressure')
        CO2saturation = sm.collect_observations_as_time_series(res, 'CO2saturation')
        mass_CO2_reservoir = sm.collect_observations_as_time_series(res, 'mass_CO2_reservoir')
        CO2_aquifer1 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer1')
        CO2_aquifer2 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer2')
        CO2_aquifer3 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer3')
        CO2_aquifer4 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer4')
        CO2_atm = sm.collect_observations_as_time_series(ms, 'CO2_atm')

        mass_CO2_aquifer1 = sm.collect_observations_as_time_series(ms, 'mass_CO2_aquifer1')
        mass_CO2_aquifer2 = sm.collect_observations_as_time_series(ms, 'mass_CO2_aquifer2')
        mass_CO2_aquifer3 = sm.collect_observations_as_time_series(ms, 'mass_CO2_aquifer3')
        mass_CO2_aquifer4 = sm.collect_observations_as_time_series(ms, 'mass_CO2_aquifer4')

        brine_aquifer1 = sm.collect_observations_as_time_series(ms, 'brine_aquifer1')
        brine_aquifer2 = sm.collect_observations_as_time_series(ms, 'brine_aquifer2')
        brine_aquifer3 = sm.collect_observations_as_time_series(ms, 'brine_aquifer3')
        brine_aquifer4 = sm.collect_observations_as_time_series(ms, 'brine_aquifer4')

        num_samples = 1
        print('Forward run is done. ')

        if isPlottingOn:
            # Plot results
            label_size = 13
            font_size = 16
            ticks_size = 12
            line_width = 1
            fig = plt.figure(figsize=(13, 12))

            ax = fig.add_subplot(321)
            plt.plot(time_array/365.25, CO2_aquifer1,
                     color='steelblue', linewidth=line_width)
            plt.xlabel('Time, years', fontsize=label_size)
            plt.ylabel('Leakage rates, kg/s', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 1', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, 50])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(322)
            plt.plot(time_array/365.25, CO2_aquifer2,
                     color='steelblue', linewidth=line_width)
            plt.xlabel('Time, years', fontsize=label_size)
            plt.ylabel('Leakage rates, kg/s', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 2', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, 50])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(323)
            plt.plot(time_array/365.25, brine_aquifer1,
                     color='rosybrown', linewidth=line_width)
            plt.xlabel('Time, years', fontsize=label_size)
            plt.ylabel('Leakage rates, kg/s', fontsize=label_size)
            plt.title(r'Leakage of brine: aquifer 1', fontsize=font_size)
            plt.tight_layout()
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, 50])
            ax.get_yaxis().set_label_coords(-0.14, 0.5)

            ax = fig.add_subplot(324)
            plt.plot(time_array/365.25, brine_aquifer2,
                     color='rosybrown', linewidth=line_width)
            plt.xlabel('Time, years', fontsize=label_size)
            plt.ylabel('Leakage rates, kg/s', fontsize=label_size)
            plt.title(r'Leakage of brine: aquifer 2', fontsize=font_size)
            plt.tight_layout()
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, 50])
            plt.ylim([1.0e-16, 1.0e-6])
            ax.get_yaxis().set_label_coords(-0.14, 0.5)

            fig.subplots_adjust(
                left=0.1, top=0.95, right=0.95, bottom=0.05, wspace=0.22)

            # to_save_png = False
            if to_save_png:
                plt.savefig('MSWObservationsCombinedPlotVsTime.png', dpi=300)
                print('"MSWObservationsCombinedPlotVsTime.png" successfully saved.')

    else:

        print('------------------------------------------------------------------')
        print('                          UQ illustration ')
        print('------------------------------------------------------------------')

        import random
        num_samples = 50
        ncpus = 1
        # Draw Latin hypercube samples of parameter values
        seed = random.randint(500, 1100)
    #    s = sm.lhs(siz=num_samples,seed=345)
        s = sm.lhs(siz=num_samples, seed=seed)

        # Run model using values in samples for parameter values
        results = s.run(cpus=ncpus, verbose=False)

        # Extract results from stochastic simulations
        outputs = s.collect_observations_as_time_series()
        CO2_aquifer1 = outputs['ms.CO2_aquifer1']
        CO2_aquifer2 = outputs['ms.CO2_aquifer2']
        CO2_aquifer3 = outputs['ms.CO2_aquifer3']
        CO2_aquifer4 = outputs['ms.CO2_aquifer4']
        brine_aquifer1 = outputs['ms.brine_aquifer1']
        brine_aquifer2 = outputs['ms.brine_aquifer2']
        brine_aquifer3 = outputs['ms.brine_aquifer3']
        brine_aquifer4 = outputs['ms.brine_aquifer4']
        mass_CO2_aquifer1 = outputs['ms.mass_CO2_aquifer1']
        mass_CO2_aquifer2 = outputs['ms.mass_CO2_aquifer2']
        mass_CO2_aquifer3 = outputs['ms.mass_CO2_aquifer3']
        mass_CO2_aquifer4 = outputs['ms.mass_CO2_aquifer4']

        mass_CO2_reservoir = outputs['res.mass_CO2_reservoir']

        print('UQ run is done. ')

        if isPlottingOn:
            # Plot result 1 : CO2 mass rate, kg/s
            label_size = 8
            font_size = 10
            ticks_size = 6
            line_width = 1
            fig = plt.figure(1, figsize=(13, 3.5))

            ax = fig.add_subplot(141)
            for j in range(num_samples):
                plt.plot(time_array/365.25, CO2_aquifer1[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 1', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(142)
            for j in range(num_samples):
                plt.plot(time_array/365.25, CO2_aquifer2[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 2', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(143)
            for j in range(num_samples):
                plt.plot(time_array/365.25, CO2_aquifer3[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 3', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(144)
            for j in range(num_samples):
                plt.plot(time_array/365.25, CO2_aquifer4[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 4', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            # to_save_png = False
            if to_save_png:
                plt.savefig('MSWObservationsCombinedPlot1VsTime.png', dpi=300)
                print('"MSWObservationsCombinedPlot1VsTime.png" successfully saved.')

            # Plot result 2 : CO2 mass, kg
            label_size = 8
            font_size = 10
            ticks_size = 6
            line_width = 1
            fig = plt.figure(2, figsize=(13, 3.5))

            # CO2 mass
            ax = fig.add_subplot(141)
            for j in range(num_samples):
                plt.plot(time_array/365.25, mass_CO2_aquifer1[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leaked mass, kg', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 1', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            ax.yaxis.get_offset_text().set_fontsize(ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(142)
            for j in range(num_samples):
                plt.plot(time_array/365.25, mass_CO2_aquifer2[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leaked mass, kg', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 2', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            ax.yaxis.get_offset_text().set_fontsize(ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(143)
            for j in range(num_samples):
                plt.plot(time_array/365.25, mass_CO2_aquifer3[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leaked mass, kg', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 3', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            ax.yaxis.get_offset_text().set_fontsize(ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(144)
            for j in range(num_samples):
                plt.plot(time_array/365.25, mass_CO2_aquifer4[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leaked mass, kg', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 4', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            ax.yaxis.get_offset_text().set_fontsize(ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            # to_save_png = False
            if to_save_png:
                plt.savefig('MSWObservationsCombinedPlot2VsTime.png', dpi=300)
                print('"MSWObservationsCombinedPlot2VsTime.png" successfully saved.')

            # Plot result 3 : brine mass rate, kg/s
            label_size = 8
            font_size = 10
            ticks_size = 6
            line_width = 1
            fig = plt.figure(3, figsize=(13, 3.5))

            ax = fig.add_subplot(141)
            for j in range(num_samples):
                plt.plot(time_array/365.25, brine_aquifer1[j],
                         color='rosybrown', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of brine: aquifer 1', fontsize=font_size)
            plt.tight_layout()
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.14, 0.5)

            ax = fig.add_subplot(142)
            for j in range(num_samples):
                plt.plot(time_array/365.25, brine_aquifer2[j],
                         color='rosybrown', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of brine: aquifer 2', fontsize=font_size)
            plt.tight_layout()
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            # plt.ylim([1.0e-16, 1.0e-6])
            ax.get_yaxis().set_label_coords(-0.14, 0.5)

            ax = fig.add_subplot(143)
            for j in range(num_samples):
                plt.plot(time_array/365.25, brine_aquifer3[j],
                         color='rosybrown', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of brine: aquifer 3', fontsize=font_size)
            plt.tight_layout()
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            # plt.ylim([1.0e-16, 1.0e-6])
            ax.get_yaxis().set_label_coords(-0.14, 0.5)

            ax = fig.add_subplot(144)
            for j in range(num_samples):
                plt.plot(time_array/365.25, brine_aquifer4[j],
                         color='rosybrown', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of brine: aquifer 4', fontsize=font_size)
            plt.tight_layout()
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            # plt.ylim([1.0e-16, 1.0e-6])
            ax.get_yaxis().set_label_coords(-0.14, 0.5)

            # to_save_png = False
            if to_save_png:
                plt.savefig('MSWObservationsCombinedPlot3VsTime.png', dpi=300)
                print('"MSWObservationsCombinedPlot3VsTime.png" successfully saved.')
