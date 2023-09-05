# -*- coding: utf-8 -*-
# Created on Aug 11, 2022
# @author: Seunghwan Baek
# seunghwan.baek@pnnl.gov

from math import sqrt
import sys
import os
import logging
import warnings
from pickle import load
import numpy as np
import matplotlib.pyplot as plt
import joblib
from scipy import interpolate
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# For now ignoring warning related to differences between versions of skikit
# used when saving and loading the ML-based models
warnings.filterwarnings("ignore", category=UserWarning)

try:
    from openiam import IAM_DIR, SystemModel, ComponentModel
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

try:
    import components.reservoir.generic.generic_reservoir_ROM as resrom
except ImportError:
    print('\nERROR: Unable to load ROM for Generic Reservoir component\n')
    sys.exit()


class GenericReservoir(ComponentModel):
    """
    The Generic Reservoir component model is a reduced-order-model to predict
    pressure and |CO2| saturation of the top part of the reservoir during |CO2| injection
    (25 years) and post-injection period (50 years).

    The model is a machine learning regression model fitted to the results of
    STOMP-CO2E multiphase flow transport simulations using Random Forest and
    scikit-learn library. Total ~2,614,400 data from ~4,000 numerical simulations
    were used to develop the model. The model predicts pressure and saturation
    only for the top part of the storage reservoir. Homogeneous reservoir model
    with radius 150 km was used. Input parameters were sampled using
    Latin Hypercube Sampling across wide ranges.

    In the NRAP-Open-IAM control file, the type name for the generic reservoir
    component is ``GenericReservoir``. The description of the component's parameters
    is provided below:

    * **reservoirDepth** [|m|] (1000 to 3500) - depth to the base of reservoir
      (default: 2000); *linked to Stratigraphy*

    * **logResPerm** [|log10| |m^2|] (-15 to -12) - logarithm of reservoir
      permeability (default: -14.0)

    * **reservoirThickness** [|m|] (15 to 200) - reservoir thickness
      (default: 50); *linked to Stratigraphy*

    * **resTempGradient** [|C/km|] (18 to 32) - reservoir temperature gradient
      (default: 30)

    * **injRate** [|kg/s|] (29 to 179) - |CO2| injection rate (default: 100)

    * **initialSalinity** [-] (0.001 to 0.05) - reservoir initial salinity
      (default: 0.05)

    * **reservoirPorosity** [-] (0.08 to 0.40) - reservoir porosity (default: 0.1).

    Possible observations from the Generic Reservoir component are:

    * **pressure** [|Pa|] - pressure at top of the reservoir at the user
      defined location(s)

    * **CO2saturation** [-] - |CO2| saturation at top of the reservoir at the
      user defined location(s).

    """
    def __init__(self, name, parent, injX=0., injY=0., locX=1000., locY=0.):

        """
        Constructor method of GenericReservoir class

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

        :returns: GenericReservoir class object
        """
        # Check that component model data is present
        if not self.model_data_check():
            pressure_folder = os.sep.join([IAM_DIR, 'source', 'components',
                                           'reservoir', 'generic', 'pressure'])
            saturation_folder = os.sep.join([IAM_DIR, 'source', 'components',
                                             'reservoir', 'generic', 'saturation'])
            error_msg = ''.join(['Generic Reservoir component {} cannot be created ',
                                 'as required model files are missing ',
                                 'in the folders\n {}\n and/or\n {}.']).format(
                                     name, pressure_folder, saturation_folder)
            # TODO Add instructions on where to download required model files once
            # we put them on EDX
            logging.error(error_msg)
            raise FileNotFoundError(error_msg)

        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25, 'time_step': 365.25}

        super().__init__(
            name, parent, model=self.simulation_model, model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'GenericReservoir'

        # Check if another generic reservoir component is already part of the system model
        orig_res_cmpnt = self._parent.ml_models_components.get(
            'GenericReservoir', None)

        sm_time_array = self._parent.time_array
        if sm_time_array is not None:
            if sm_time_array[-1]/365.25 > 75.0:
                warn_msg = ''.join([
                    'The last simulation time point ({:.1f} yr) of a system ',
                    'model goes beyond the GenericReservoir model capability ',
                    '(max. 75 yr)']).format(sm_time_array[-1]/365.25)
                logging.warning(warn_msg)

        # Set default parameters of the component model
        self.add_default_par('reservoirDepth', value=2000.0)
        self.add_default_par('reservoirThickness', value=50.0)
        self.add_default_par('logResPerm', value=-14.0)
        self.add_default_par('resTempGradient', value=30.0)
        self.add_default_par('injRate', value=100)
        self.add_default_par('initialSalinity', value=0.05)
        self.add_default_par('reservoirPorosity', value=0.1)

        self.add_default_par('surfaceTemperature', value=15.0)
        self.add_default_par('datumPressure', value=101325.0)
        self.add_default_par('resPressGradient', value=9.785602)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['reservoirDepth'] = [1000.0, 3500.0]
        self.pars_bounds['reservoirThickness'] = [15.0, 200.0]
        self.pars_bounds['logResPerm'] = [-15.0, -12.0]
        self.pars_bounds['resTempGradient'] = [18.0, 32.0]
        self.pars_bounds['injRate'] = [29.0, 179.0]
        self.pars_bounds['initialSalinity'] = [1e-3, 5e-2]
        self.pars_bounds['wellRadius'] = [0.017, 0.073]
        self.pars_bounds['reservoirPorosity'] = [0.08, 0.40]

        # Specify default locations of injection and leaking wells
        self.assign_coordinates(injX, injY, locX, locY)

        # Initiate solution object
        self.sol = resrom.Solution()

        componentPath0 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        componentPath = os.path.join(componentPath0, 'components', 'reservoir')

        TargetVars = ['pressure', 'saturation']

        if orig_res_cmpnt is None:
            # Register the component as the one using ml models
            # The component will be registered only if no other generic reservoir
            # components were added to the same system model
            self._parent.register_ml_model_component(self, 'GenericReservoir')

            # Follow the originally developed design: loading models
            for target_var in TargetVars:

                ScalerXfilename = os.path.join(componentPath, 'generic', target_var,
                                               'scaler_features_{}.pkl'.format(target_var))
                ScalerYfilename = os.path.join(componentPath, 'generic', target_var,
                                               'scaler_targets_{}.pkl'.format(target_var))
                Modelfilename = os.path.join(componentPath, 'generic', target_var,
                                             '{}RegModel1.joblib'.format(target_var))

                if target_var == 'pressure':
                    with open(ScalerXfilename, "rb") as f:
                        self.x_scaler_p = load(f)

                    with open(ScalerYfilename, "rb") as f:
                        self.y_scaler_p = load(f)

                    with open(Modelfilename, "rb") as f:
                        self.model_p = joblib.load(f)
                else:
                    with open(ScalerXfilename, "rb") as f:
                        self.x_scaler_s = load(f)

                    with open(ScalerYfilename, "rb") as f:
                        self.y_scaler_s = load(f)

                    with open(Modelfilename, "rb") as f:
                        self.model_s = joblib.load(f)

        else: # If another generic reservoir component was already added
            # save references to the models loaded in that component
            self.x_scaler_p = orig_res_cmpnt.x_scaler_p
            self.y_scaler_p = orig_res_cmpnt.y_scaler_p
            self.model_p = orig_res_cmpnt.model_p
            self.x_scaler_s = orig_res_cmpnt.x_scaler_s
            self.y_scaler_s = orig_res_cmpnt.y_scaler_s
            self.model_s = orig_res_cmpnt.model_s

        debug_msg = 'GenericReservoir component created with name {}'.format(self.name)
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
            if key in self.pars_bounds:
                if (val < self.pars_bounds[key][0]) or (val > self.pars_bounds[key][1]):
                    warn_msg = ''.join([
                        'Parameter {} (value {}) of GenericReservoir component {} ',
                        'is out of boundaries.']).format(key, val, self.name)
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of GenericReservoir component {}.']).format(key, self.name)
                logging.warning(warn_msg)

        if self.distance >= 150000:
            warn_msg = 'The distance between the injection well and leaky well \
            is greater than the boundary size of GenericReservoir component'
            logging.warning(warn_msg)

    def assign_coordinates(self, injX=None, injY=None, locX=None, locY=None):
        """
        Save locations of injector and leaking well and calculate distance
        between them.

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
        """
        if injX is not None:
            self.injX = injX
        if injY is not None:
            self.injY = injY
        if locX is not None:
            self.locX = locX
        if locY is not None:
            self.locY = locY

        self.distance = sqrt((self.locX-self.injX)**2+(self.locY-self.injY)**2)

    def simulation_model(self, p, time_point=365.25, time_step=365.25,
                         injX=None, injY=None, locX=None, locY=None):
        """
        Return pressure and CO2 saturation at the bottom of leaking well.
        Note that the names of parameters contained in the input dictionary
        coincide with those defined in this modules docstring.

        :param p: input parameters of generic reservoir model
        :type p: dict

        :returns: out - dictionary of observations of generic reservoir
            component model;
            keys: ['pressure', 'CO2saturation']

        """
        debug_msg = ''.join([
            'GenericReservoir component {} model call input: ',
            'parameters {}, location xy ({}, {})']).format(
                self.name, p, self.locX, self.locY)
        logging.debug(debug_msg)

        # Assign default values
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        actual_p.update(p)

        # Default value not allowed to be changed by user
        # well radius w/ porosity 1.0
        actual_p['wellRadius'] = 0.05

        # This statement should take care of the setup of locations through
        # control file interface, updates of distance between locations and injection
        # wells, as well as cases when locations for observations are
        # generated randomly for each realization
        if time_point == 0.0:
            self.assign_coordinates(injX, injY, locX, locY)

        inputArray = np.array([
            actual_p['reservoirDepth'], actual_p['reservoirThickness'],
            actual_p['logResPerm'], actual_p['resTempGradient'],
            actual_p['injRate'], actual_p['initialSalinity'],
            actual_p['wellRadius'], actual_p['reservoirPorosity'],
            actual_p['surfaceTemperature'], actual_p['datumPressure'],
            actual_p['resPressGradient'], self.distance])

        out = {}

        if time_point == 0.0:
            # Solve once
            self.sol.find(inputArray, self.x_scaler_p, self.y_scaler_p, self.model_p,
                          self.x_scaler_s, self.y_scaler_s, self.model_s)

            # Create interpolation functions
            self.func_p = interpolate.interp1d(
                list(range(76)), self.sol.Outputs['pressure'], kind='cubic')

            self.func_s = interpolate.interp1d(
                list(range(76)), self.sol.Outputs['saturation'], kind='cubic')

            out['pressure'] = self.func_p(time_point/365.25)[0]

            out['CO2saturation'] = 0.0

        else:
            if time_point/365.25 > 75.0:
                warn_msg = ''.join([
                    'The last simulation time point ({:.1f} yr) goes beyond ',
                    'the model capability (max. 75 yr)']).format(time_point/365.25)
                logging.warning(warn_msg)

            # Apply the interpolation function
            out['pressure'] = self.func_p(time_point/365.25)[0]

            out['CO2saturation'] = np.clip(self.func_s(time_point/365.25), 0.0, 1.0)[0]

        return out

    # Attributes for system connections
    needs_locXY = True
    needs_injXY = True
    system_params = ['reservoirThickness', 'reservoirDepth']

    @staticmethod
    def model_data_check():
        """
        Check whether model data files for Generic Reservoir component was downloaded
        and placed to the right place.
        """
        model_data_dir = os.sep.join([
            IAM_DIR, 'source', 'components', 'reservoir', 'generic'])
        outputs = ['pressure', 'saturation']

        model_files = {'pressure': ['PressureRegModel1.joblib',
                                    'scaler_features_Pressure.pkl',
                                    'scaler_targets_Pressure.pkl'],
                       'saturation': ['SaturationRegModel1.joblib',
                                      'scaler_features_Saturation.pkl',
                                      'scaler_targets_Saturation.pkl']}

        check_flag = 0   # no issues: data is present

        for output_nm in outputs:
            for file in model_files[output_nm]:
                file_path = os.sep.join([model_data_dir, output_nm, file])
                if not os.path.isfile(file_path):
                    check_flag = check_flag + 1

        return (check_flag == 0)


if __name__ == "__main__":
    __spec__ = None
    logging.basicConfig(level=logging.WARNING)
    # Define keyword arguments of the system model
    to_save = False    # figure save

    distance = 2000    # in meters

    num_years = 75
    time_array = 365.25*np.arange(0, num_years+1, 5)

    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    res = sm.add_component_model_object(GenericReservoir(
        name='res', parent=sm, injX=0., injY=0., locX=distance, locY=0.))

    # Add parameters of reservoir component model
    res.add_par('reservoirDepth', value=2000.0, vary=False)
    res.add_par('reservoirThickness', value=50.0, vary=False)
    res.add_par('logResPerm', value=-14.0, vary=False)
    res.add_par('resTempGradient', value=30.0, vary=False)
    res.add_par('injRate', value=100, vary=False)
    res.add_par('initialSalinity', value=0.05, vary=False)
    res.add_par('wellRadius', value=0.05, vary=False)
    res.add_par('reservoirPorosity', value=0.1, vary=False)

    # Add observations (output) from the reservoir model
    res.add_obs('pressure')
    res.add_obs('CO2saturation')

    # Run system model using current values of its parameters
    import time
    t0 = time.time()

    sm.forward()

    print("   model running time %.3f secs"%(time.time() - t0))

    # Collect results
    # pressure (defined time) at defined location
    pressure_at_well = sm.collect_observations_as_time_series(res, 'pressure')

    # saturation (defined time) at defined location
    saturation_at_well = sm.collect_observations_as_time_series(res, 'CO2saturation')

    # # Expected
    # pressure_at_well = [
    #     1.91833e+07, 4.92151e+07, 5.48526e+07, 5.77966e+07, 5.96888e+07,
    #     6.10506e+07, 4.28154e+07, 3.70276e+07, 3.37106e+07, 3.16047e+07,
    #     3.00346e+07, 2.88335e+07, 2.78929e+07, 2.71535e+07, 2.653e+07,
    #     2.60014e+07]

     # saturation_at_well = [
     #    0.00000e+00, 8.20977e-01, 9.80491e-01, 9.92484e-01, 9.94269e-01,
     #    9.95659e-01, 9.97513e-01, 9.97742e-01, 9.97767e-01, 9.97757e-01,
     #    9.97733e-01, 9.97702e-01, 9.97664e-01, 9.97621e-01, 9.97576e-01,
     #    9.97530e-01]

    # -------------------------------------------------------------------------
    # # Plot pressure and saturation at a leaky well
    # -------------------------------------------------------------------------
    font_size = 6
    grid_line_width = 0.1
    line_width = 0.5
    label_size = 6
    Time_array_yr = time_array/365.25
    Pressure_MPa = pressure_at_well*1e-6
    fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(3, 5), dpi=200)
    ax1.plot(Time_array_yr, Pressure_MPa, linewidth=line_width,
             color='blue', label='at {:.1f} m'.format(distance))
    ax1.legend(loc='upper right', fontsize=font_size)
    ax1.set_xlabel('Time,yr', fontsize=font_size)
    ax1.set_ylabel('Pressure, MPa', fontsize=font_size)
    ax1.set_ylim((np.max([np.min(Pressure_MPa)-3, 0]), np.max(Pressure_MPa)+3))
    ax1.set_xlim((np.min(Time_array_yr), np.max(Time_array_yr)))
    ax1.tick_params(axis='both', which='both', length=0, labelsize=label_size)
    ax1.grid(color='grey', linestyle='-', linewidth=grid_line_width)

    ax2.plot(Time_array_yr, saturation_at_well, linewidth=line_width, color='red')
    ax2.set_xlabel('Time,yr', fontsize=font_size)
    ax2.set_ylabel(r'CO$_2$ saturation', fontsize=font_size)
    ax2.set_ylim((-0.05, 1.05))
    ax2.set_xlim((np.min(Time_array_yr), np.max(Time_array_yr)))
    ax2.tick_params(axis='both', which='both', length=0, labelsize=label_size)
    ax2.grid(color='grey', linestyle='-', linewidth=grid_line_width)
    fig.tight_layout()

    if to_save:
        fig.savefig(os.path.join(os.getcwd(), 'generic_reservoir_results.png'),
                    dpi=500, bbox_inches='tight')
    plt.show()
