# -*- coding: utf-8 -*-
import os
import sys
import logging
# import tensorflow as tf
# import joblib
import numpy as np
import pandas as pd

# Common command to create a new component based on ComponentModel class
# so that it will be using its inherent basic methods and attributes
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Since this component is based on ComponentModel we need to import it
try:
    from openiam import SystemModel, ComponentModel
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

from openiam.cfi.commons import process_parameters

# try:
#     import components.fault.fault_leakage.fault_leakage_model as flmod
# except ImportError:
#     print('\nERROR: Unable to load models for Fault Leakage component\n')
#     sys.exit()


class FaultLeakage(ComponentModel):
    """
    The Fault Leakage Model component uses deep neural networks to estimate
    the flow of brine carbon dioxide along a fault. The dynamics of multiphase flow
    in the storage reservoir and shallow aquifer into which the fault leaks are
    both taken into consideration in the estimation of leakage rates. The |CO2| is
    treated as supercritical throughout the leakage process, including within the
    shallow aquifer. No phase change occurs. The fault is modeled as a continuous,
    homogeneous, isotropic porous medium with Darcy flow. This component assumes
    the following:

    * The storage reservoir is 4000 m long (distance away from the fault),
      2400 m wide, and 200 m thick.

    * The shallow aquifer has similar dimensions to the storage reservoir, but
      is located on the opposite side of the fault.

    * The fault has a thickness of 3 m and a width of 2400 m, contacting the entire
      width of both the storage reservoir and the shallow aquifer.

    * The fault dip angle varies.

    * There is a caprock on both sides of the fault with a thickness of 100 m
      perpendicular to the fault. The caprock has a permeability of 1.0e-18
      (impermeable) and a porosity of 0.1.

    * The storage reservoir has no-flow boundaries.

    * The shallow aquifer has a constant pressure boundary at its side boundary.
      All other boundaries are no-flow.

    * The domain is 1750 m in depth from the top of the shallow aquifer
      to the bottom of the storage reservoir.

    * The pressure and temperature at the top of the shallow aquifer are
      8.15 MPa and 50 |C|, respectively.

    * The |CO2| injection temperature is 32 |C|.

    * The thermal conductivity of rock is 3 W/(m*K).

    * The specific heat capacity of rock is 920 J/(kg*K).

    * Well is located about 178-190 m above the bottom of the storage reservoir.

    The description of the component's parameters is provided below:

    * **damage_zone_perm** [|log10| |m^2|] (-15 to -12) - the permeability
      of the fault, including both the fault core and damage zone (default: -13.5)

    * **damage_zone_por** [-] (0.001 to 0.1) - the porosity of the fault,
      including both the fault core and damage zone (default: 0.01)

    * **shallow_aquifer_perm** [|log10| |m^2|] (-14 to -12) - the permeability
      of the shallow aquifer (default: -13.0)

    * **deep_aquifer_perm** [|log10| |m^2|] (-14 to -12) - the permeability
      of the storage aquifer (e.g., reservoir) (default: -13.0)

    * **shallow_aquifer_por** [-] (0.05 to 0.5) - the porosity of the shallow
      aquifer (default: 0.25)

    * **deep_aquifer_por** [-] (0.05 to 0.35) - the porosity of the deep aquifer
      (default: 0.2)

    * **well_index** [-] (integer: 0, 1, 2) - a proxy for the horizontal distance
      of the well from the fault; value of 0 means that the well is about 200 m from
      the fault, value of 1 means that the well is about 400 m from the fault;
      value of 2 means that the well is about 600 m from the fault (default: 0)

    * **well_rate** [|kg/s|] (0.5 to 25) - injection rate of the well
      for the aquifer (default: 15.8)

    * **dip_angle** [|deg|] (integer: 40, 60, 80, 100, 120, 140) - dip angle
      of the fault, measured from the horizontal plane (default: 60)

    * **injection_time** [|years|] (10 to 50) - duration of injection (default: 30)

    * **geothermal_gradient** [|C/km|] (8 to 44) - the geothermal gradient
      in the formation (default: 30)

    The possible outputs from the Fault Leakage component are the leakage rates
    and cumulative leakage amounts of brine and |CO2| to the shallow aquifer
    through the fault. The names of the observations are:

    * **brine_aquifer**, **CO2_aquifer** [|kg/s|] - brine and |CO2|
      leakage rates to the shallow aquifer, respectively.

    * **mass_brine_aquifer**, **mass_CO2_aquifer** [|kg|] - cumulative
      brine and |CO2| leakage into the shallow aquifer, respectively.

    Notes:

    * Due to the use of the trapezoidal rule in integrating instantaneous
      leakage rates, the cumulative leakage values might not conserve mass
      for |CO2| depending on the time step size used.

    * Leakage estimates are available only up to 100 years.

    * After extensive sensitivity analysis the key parameters for the model
      have been found to be: damage zone permeability, deep aquifer permeability,
      injection rate of the well, injection duration, and dip angle. The time
      at which the leakage rate or amount is to be estimated (i.e., simulation
      time) is also a key parameter. Uncertainties in these variables are,
      therefore, most important.

    """
    def __init__(self, name, parent):
        """
        Constructor method of FaultLeakageModel class

        :param name: name of component model
        :type name: [str]

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel [object]

        :returns: FaultLeakageModel class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        # Note that time_point is the exact time currently, and time_step is the current time step.
        # Units of time are in days here.
        model_kwargs = {'time_point': 365.25, 'time_step': 365.25}

        super().__init__(
            name, parent, model=self.simulation_model, model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'FaultLeakage'

        # Set default parameters of the component model
        self.add_default_par('damage_zone_perm', value=-13.5)
        self.add_default_par('damage_zone_por', value=0.01)
        self.add_default_par('shallow_aquifer_perm', value=-13.0)
        self.add_default_par('deep_aquifer_perm', value=-13.0)
        self.add_default_par('shallow_aquifer_por', value=0.25)
        self.add_default_par('deep_aquifer_por', value=0.2)
        self.add_default_par('well_index', value=0)
        self.add_default_par('well_rate', value=15.8)
        self.add_default_par('dip_angle', value=60)
        self.add_default_par('injection_time', value=30.0)
        self.add_default_par('geothermal_gradient', value=30.0)

        # Define dictionary of parameters boundaries
        # The following lines define the dictionary attribute pars_bounds containing
        # the upper and lower boundaries for each parameter of the model.
        # This attribute is used in the method check_input_parameters of the
        # ComponentModel class. The method check_input_parameters is called
        # before the start of each simulation to check whether the provided
        # parameters satisfy the defined boundaries
        # Each entry of pars_bounds dictionary is a list of lower and upper
        # boundary, respectively, for the correspoding parameter
        self.pars_bounds = dict()
        self.pars_bounds['damage_zone_perm'] = [-15.0, -12.0]     # log_10 m^2
        self.pars_bounds['damage_zone_por'] = [0.001, 0.1]
        self.pars_bounds['shallow_aquifer_perm'] = [-14.0, -12.0] # log_10 m^2
        self.pars_bounds['deep_aquifer_perm'] = [-14.0, 12.0]     # log_10 m^2
        self.pars_bounds['shallow_aquifer_por'] = [0.05, 0.5]
        self.pars_bounds['deep_aquifer_por'] = [0.05, 0.35]
        self.pars_bounds['well_index'] = [0, 1, 2]
        self.pars_bounds['well_rate'] = [0.5, 25.0]               # kg/s, full aquifer. not half
        self.pars_bounds['dip_angle'] = np.linspace(40, 140, 6, dtype=int)
        self.pars_bounds['injection_time'] = [10.0, 50.0]         # yr
        self.pars_bounds['geothermal_gradient'] = [8.0, 44.0]     # degrees C/km

        # Define dictionary of temporal data limits
        self.temp_data_bounds = dict()
        self.temp_data_bounds['time'] = ['Time', 0.0, 100.0*365.25] # in days

        # Define accumulators and their initial values.
        # The accumulators are used to track cumulative brine and CO2 leakage
        # and the previous instantaneous brine and CO2 leakage rates (for
        # use in calculating the cumulative leakage values using the trapezoidal rule).
        self.add_accumulator('cumulative_brine_leakage', sim=0.0)
        self.add_accumulator('cumulative_CO2_leakage', sim=0.0)
        self.add_accumulator('previous_brine_leakage_rate', sim=0.0)
        self.add_accumulator('previous_CO2_leakage_rate', sim=0.0)

    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        debug_msg = 'Input parameters of {name} component are {p}.'.format(
            name=self.name, p=p)
        logging.debug(debug_msg)

        for key, val in p.items():
            if key == 'well_index':
                if int(p['well_index']) not in self.pars_bounds['well_index']:
                    msg = ''.join(['Invalid value {} of well_index parameter: ',
                                   'not in the list [0, 1, 2].']).format(val)
                    logging.error(msg)
                    raise ValueError(msg)
                continue
            if key == 'dip_angle':
                if int(p['dip_angle']) not in self.pars_bounds['dip_angle']:
                    msg = ''.join(['Invalid value {} of dip_angle parameter: ',
                                   'not in the list ',
                                   '[40, 60, 80, 100, 120, 140].']).format(val)
                    logging.error(msg)
                    raise ValueError(msg)
                continue
            if key in self.pars_bounds:
                if ((val < self.pars_bounds[key][0])or(val > self.pars_bounds[key][1])):
                    warn_msg = 'Parameter {key} is out of boundaries.'.format(key=key)
                    logging.warning(warn_msg)
            else:
                warn_msg = ''.join([
                    'Parameter {key} is not recognized as component ',
                    '{name} input parameter.']).format(key=key, name=self.name)
                logging.warning(warn_msg)

    def connect_with_system(self, component_data, name2obj_dict, *args, **kwargs):
        """
        Code to add fault leakage component to system model for control file interface.

        :param component_data: Dictionary of component data
        :type component_data: dict

        :param name2obj_dict: Dictionary to map component names to objects
        :type name2obj: dict

        :returns: None
        """
        # Process parameters data if any
        process_parameters(self, component_data, name2obj_dict)

        if 'Outputs' in component_data:
            comp_outputs = component_data['Outputs']
            for obs_nm in comp_outputs:
                self.add_obs(obs_nm)

    def simulation_model(self, p, time_point=365.25, time_step=365.25):
        """
        :param p: input parameters of FaultLeakage model
        :type p: dict

        param time_point: time point in days at which the leakage rates are
            to be calculated; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param time_step: difference between the current and previous
            time points in days; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        The possible outputs from the Fault Leakage component are leakage rates
        and cumulative masses of |CO2| and brine leaked. The names of the
        observations are of the form:

        * **CO2_aquifer**, **brine_aquifer** [|kg/s|] - instantaneous |CO2| and
          brine leakage rates

        * **mass_CO2_aquifer**, **mass_brine_aquifer** [|kg|] - cumulative mass
          of |CO2| and brine.
        """
        if time_point/365.25 > 100.0:
            err_msg = 'FaultLeakage Component model is valid only up to 100 years.'
            raise ValueError(err_msg)
        # Obtain the default values of the parameters from dictionary
        # of default parameters
        # Note: the default parameters are out of bounds for the parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # This component assumes that observations are equal to zero at the
        # initial time point.
        tol = 1e-8
        if time_point < tol:
            # Define initial values of the model observations
            # These are the initial conditions (no leakage)
            out = {}
            out['brine_aquifer'] = 0.0       # metric tons/yr, instantaneous brine leakage rate
            out['mass_brine_aquifer'] = 0.0  # metric tons, cumulative brine leakage
            out['CO2_aquifer'] = 0.0         # metric tons/yr, instantaneous CO2 leakage rate
            out['mass_CO2_aquifer'] = 0.0    # metric tons, cumulative CO2 leakage
            # Exit method
            return out

        # Check whether the temporal inputs satisfy the model requirements
        assumptions_satisfied = 0
        if (self.temp_data_bounds['time'][1] <= time_point <= self.temp_data_bounds['time'][2]):
            assumptions_satisfied = 1

        if assumptions_satisfied:
            # Calculate output of the component using parameters and
            # temporal keyword arguments.
            # This is where the deep-neural network ROM for fault leakage
            # is called by this Python script.
            X = np.expand_dims(np.array(
                [actual_p['damage_zone_perm'], actual_p['damage_zone_por'],
                 actual_p['shallow_aquifer_perm'], actual_p['deep_aquifer_perm'],
                 actual_p['shallow_aquifer_por'], actual_p['deep_aquifer_por'],
                 # The well rate taken is assumed to be for the FULL aquifer
                 # and is multiplied by 0.5 to account for the fact
                 # that the ROM requires HALF of this injection rate.
                 actual_p['well_index'], actual_p['well_rate']/2.0,
                 actual_p['dip_angle'], actual_p['injection_time'],
                 # time needs to have units of years before inputting it to the ROM
                 actual_p['geothermal_gradient'], time_point/365.25]), axis=0)
            # This line is for suppressing the feature naming warnings with the scalers.
            X = pd.DataFrame(X, columns=['damage_zone_perm', 'damage_zone_por',
                                         'shallow_aquifer_perm', 'deep_aquifer_perm',
                                         'shallow_aquifer_por', 'deep_aquifer_por',
                                         'well_index', 'well_rate', 'dip_angle',
                                         'injection_time', 'geothermal_gradient',
                                         'time'])
            # Run the model.
            import components.fault.fault_leakage.fault_leakage_model as flmod
            X_scaled_brine = flmod.fl_scaler('brine').transform(X)
            X_scaled_CO2 = flmod.fl_scaler('CO2').transform(X)
            # t/yr, multiplied by 2.0 to account for the FULL aquifer
            # (ROM is for a HALF aquifer with symmetry)
            brine_aquifer = 2.0*flmod.fl_model('brine').predict(X_scaled_brine)[0][0]
            # t/yr, multiplied by 2.0 to account for the FULL aquifer
            # (ROM is for a HALF aquifer with symmetry)
            CO2_aquifer = 2.0*flmod.fl_model('CO2').predict(X_scaled_CO2)[0][0]
            # Assign values to the component accumulators.
            # The trapezoidal rule is used to integrate the brine and CO2
            # cumulative leakage. This can lead to cumulative leakage values
            # that do not follow the conservation of mass (e.g., cumulative CO2
            # leakage mass exceeding injected CO2 mass) due to approximation error.
            prev_brine_aquifer = self.accumulators['previous_brine_leakage_rate'].sim
            prev_CO2_aquifer   = self.accumulators['previous_CO2_leakage_rate'].sim
            mass_brine_aquifer = self.accumulators['cumulative_brine_leakage'].sim \
                + 0.5*(prev_brine_aquifer+brine_aquifer)*time_step/365.25
            mass_CO2_aquifer   = self.accumulators['cumulative_CO2_leakage'].sim \
                + 0.5*(prev_CO2_aquifer+CO2_aquifer)*time_step/365.25

            self.accumulators['cumulative_brine_leakage'].sim = mass_brine_aquifer
            self.accumulators['cumulative_CO2_leakage'].sim = mass_CO2_aquifer
            self.accumulators['previous_brine_leakage_rate'].sim = brine_aquifer
            self.accumulators['previous_CO2_leakage_rate'].sim = CO2_aquifer

            # Assign model observations
            out = {}
            # Outputs need conversion from tonnes/year to kg/s
            ty_to_kgs = 1.0e+3/(365.25*24*3600)
            tonne_to_kg = 1.0e+3
            # kg/s, instantaneous brine leakage rate
            out['brine_aquifer'] = ty_to_kgs * brine_aquifer
            # kgs, cumulative brine leakage
            out['mass_brine_aquifer'] = tonne_to_kg * mass_brine_aquifer
            # kg/s, instantaneous CO2 leakage rate
            out['CO2_aquifer'] = ty_to_kgs * CO2_aquifer
            # kg, cumulative CO2 leakage
            out['mass_CO2_aquifer'] = tonne_to_kg * mass_CO2_aquifer

            return out

# Add test information here.
if __name__ == "__main__":
    # Define keyword arguments of the system model.

    times = np.linspace(0.0, 30.0, num=20)*365.25 # time in days
    sm_model_kwargs = {'time_array': times} # time must be given in days

    # Create the system model.
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # well_index is the well location index (not to be confused with Peaceman's well index)
    params = {'damage_zone_perm': -13.0, 'damage_zone_por': 0.05,
              'shallow_aquifer_perm': -13.0, 'deep_aquifer_perm': -12.0,
              'shallow_aquifer_por': 0.30, 'deep_aquifer_por': 0.20,
              'well_index': 0.0, 'well_rate': 20.0, 'dip_angle': 40.0,
              'injection_time': 20.0, 'geothermal_gradient': 35.0}

    ffc = sm.add_component_model_object(FaultLeakage(name='ffc', parent=sm))

    # Add parameters.
    for key in params:
        ffc.add_par(key, value=params[key], vary=False)

    # Add scalar observations
    ffc.add_obs('brine_aquifer')
    ffc.add_obs('mass_brine_aquifer')
    ffc.add_obs('CO2_aquifer')
    ffc.add_obs('mass_CO2_aquifer')

    print('Starting simulation...')
    sm.forward()
    print('Simulation is finished...')

    # Collect observations and convert from kg/sec to tonne/year
    kg_to_tonne = 1.0e-3
    kgs_to_ty = kg_to_tonne*365.25*24*3600

    brine_aquifer_obs = kgs_to_ty * sm.collect_observations_as_time_series(
        ffc, 'brine_aquifer')
    mass_brine_aquifer_obs = kg_to_tonne * sm.collect_observations_as_time_series(
        ffc, 'mass_brine_aquifer')
    CO2_aquifer_obs = kgs_to_ty * sm.collect_observations_as_time_series(
        ffc, 'CO2_aquifer')
    mass_CO2_aquifer_obs = kg_to_tonne * sm.collect_observations_as_time_series(
        ffc, 'mass_CO2_aquifer')

    # Print the results.
    print('Brine leakage rate (t/yr)', brine_aquifer_obs, sep='\n')
    print('Cumulative brine leakage (t)', mass_brine_aquifer_obs, sep='\n')
    print('CO2 leakage rate (t/yr)', CO2_aquifer_obs, sep='\n')
    print('Cumulative CO2 leakage (t)', mass_CO2_aquifer_obs, sep='\n')

    # Perform an accuracy test with a known case.
    true_brine_aquifer =  np.array([
        0., 175979.9, 145872.62, 129209.84, 125567.67,
        127366.84, 134765.28, 140630.7, 146373.23, 149061.16,
        150104.33, 148167.31, 146717.6, 124379.67, 84806.8,
        68673.25, 58678.89, 48936.89, 40774.6, 33600.992])   # t/yr

    true_CO2_aquifer = np.array([
        0., 239869.08, 445098.28, 467009.75, 476176.1,
        478896.06, 479526.47, 482197.6, 489237., 504251.44,
        529424.25, 557741.75, 576606.5, 414577.22, 214732.17,
        110940.52, 69257.9, 51961.707, 46253.23, 42923.3])   # t/yr

    true_mass_brine_aquifer = np.array([
        0., 138931.50493421, 393025.60855263, 610195.97861842,
        811336.12253289, 1011021.26644737, 1217967.68092105, 1435385.57565789,
        1661967.63157895, 1895205.29605263, 2131388.58552632, 2366866.18421053,
        2599670.05756579, 2813694.20230263, 2978841.41447368, 3100009.87253289,
        3200551.03618421, 3285510.86348684, 3356335.72574013, 3415053.29975329])    # t/yr

    true_mass_CO2_aquifer = np.array([
        0., 189370.32483553, 730134.04194079, 1450219.30509868,
        2194839.73273026, 2948844.04194079, 3705493.38404605, 4464749.22286184,
        5231671.29523026, 6016004.27220395, 6832064.02549342, 7690352.97286184,
        8585891.06496711, 9368404.55180921, 9865227.74259869, 10122337.75904606,
        10264599.67105263, 10360299.36266448, 10437837.47121711, 10508239.99588816]) # t/yr

    # Set tolerance to have a stringent tolerance while accounting
    # for rounding in the true values.
    tol = 1e-4
    delta_brine_aquifer = brine_aquifer_obs-true_brine_aquifer
    delta_CO2_aquifer = CO2_aquifer_obs-true_CO2_aquifer
    delta_brine_mass = mass_brine_aquifer_obs-true_mass_brine_aquifer
    delta_CO2_mass = mass_CO2_aquifer_obs-true_mass_CO2_aquifer
    assert np.max(abs(delta_brine_aquifer[1:])/true_brine_aquifer[1:]) < tol,\
        'Brine leakage rates failed the test for the fault leakage component.'
    assert np.max(abs(delta_CO2_aquifer[1:])/true_CO2_aquifer[1:]) < tol, \
        'CO2 leakage rates failed the test for the fault leakage component.'
    assert np.max(abs(delta_brine_mass[1:])/true_mass_brine_aquifer[1:]) < tol, \
        'Total brine leakage failed the test for the fault leakage component.'
    assert np.max(abs(delta_CO2_mass[1:])/true_mass_CO2_aquifer[1:]) < tol, \
        'Total CO2 leakage failed the test for the fault leakage component.'

    # Test the zeros separately
    assert brine_aquifer_obs[0] < tol, \
        'Brine leakage rates failed the test for the fault leakage component.'
    assert CO2_aquifer_obs[0] < tol, \
        'CO2 leakage rates failed the test for the fault leakage component.'
    assert mass_brine_aquifer_obs[0] < tol, \
        'Total brine leakage failed the test for the fault leakage component.'
    assert mass_CO2_aquifer_obs[0] < tol, \
        'Total CO2 leakage failed the test for the fault leakage component.'
