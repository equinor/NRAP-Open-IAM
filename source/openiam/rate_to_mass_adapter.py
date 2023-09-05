"""
This module contains RateToMassAdapter class.

Created in September 2017
Last modified in November 2019

@author: Veronika Vasylkivska, Seth King

"""
import sys
import os
import logging
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from openiam import SystemModel, ComponentModel
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))


class RateToMassAdapter(ComponentModel):
    """ NRAP-Open-IAM RateToMassAdapter component class. """
    def __init__(self, name, parent, calc_method=2):
        """
        Constructor method of RateToMassAdapter class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param calc_method: method used to calculate the accumulated mass.
            Can be 1 or 2. 1 corresponds to using only the incoming rates for the
            calculation of accumulated masses. 2 which is a default value corresponds to
            using both the incoming and outcoming rates.
        :type calc_method: int

        :returns: RateToMassAdapter class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_step': 365.25}   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'RateToMassAdapter'

        self.calc_method = calc_method

    def simulation_model(self, p, time_step=365.25, **leakage_rates):
        """
        Return CO2 and/or brine accumulated mass depending on provided
        parameters.

        :param p: dictionary of input parameters
        :type p: dict

        :param time_step: difference between current and previous time points
            in days; by default, its value is 365.25 (1 year in days)
        :type time_step: float

        :returns: dictionary of observations (accumulated masses of CO2 and brine)
            where the keys will be each leakage_rates key name prepended by "mass" and underscore
        """
        # Check whether the accumulators were already added to the component (self)
        if len(self.accumulators) == 0:
            # For every parameter in dictionary leakage_rates add a corresponding accumulator
            for leak_rates_nm in leakage_rates:
                self.add_accumulator('mass_{}'.format(leak_rates_nm), sim=0.0)

        scalar = 24*3600  # number of seconds in day
        out = dict()
        for leak_rates_nm, leak_rates_vals in leakage_rates.items():
            # Setup the observation name
            obs_nm = 'mass_{}'.format(leak_rates_nm)
            # Leakage rates are in kg/s, time_step is in days, scalar is present
            # to convert to the correct units
            if self.calc_method == 1:
                # Only entering leakage rates are considered
                self.accumulators[obs_nm].sim = max(
                    0.0,
                    self.accumulators[obs_nm].sim + leak_rates_vals[0]*time_step*scalar)
            elif self.calc_method == 2:
                # Setup the new value in the accumulator considering the difference
                # in the entering and leaving leakage rates and time step
                # during which the leakage occurs.
                self.accumulators[obs_nm].sim = max(
                    0.0, (self.accumulators[obs_nm].sim +
                          (leak_rates_vals[0]-leak_rates_vals[1])*time_step*scalar))
            out[obs_nm] = self.accumulators[obs_nm].sim

        # Return dictionary of accumulated masses
        return out

    def connect_with_system(self, component_data, name2obj_dict, *args, **kwargs):
        """
        Code to add RateToMassAdapters to system model for control file interface.

        :param component_data: Dictionary of component data including 'Connection'
        :type component_data: dict

        :param name2obj_dict: Dictionary to map component names to objects
        :type name2obj_dict: dict

        :returns: None
        """
        strata = name2obj_dict['strata']
        if 'numberOfShaleLayers' in strata.deterministic_pars:
            num_shale_layers = strata.deterministic_pars['numberOfShaleLayers'].value
        elif 'numberOfShaleLayers' in strata.default_pars:
            num_shale_layers = strata.default_pars['numberOfShaleLayers'].value
        else:
            num_shale_layers = 3

        well = name2obj_dict[component_data['Connection']]
        well.add_obs_to_be_linked('CO2_aquifer1')
        well.add_obs_to_be_linked('brine_aquifer1')
        self.add_obs_to_be_linked('mass_CO2_aquifer1')
        self.add_obs_to_be_linked('mass_brine_aquifer1')

        for il in range(2, num_shale_layers):
            well.add_obs_to_be_linked('CO2_aquifer{il}'.format(il=il))
            well.add_obs_to_be_linked('brine_aquifer{il}'.format(il=il))
            self.add_kwarg_linked_to_collection(
                'CO2_aquifer{ilm1}'.format(ilm1=il-1),
                [well.linkobs['CO2_aquifer{ilm1}'.format(ilm1=il-1)],
                 well.linkobs['CO2_aquifer{il}'.format(il=il)]])
            self.add_kwarg_linked_to_collection(
                'brine_aquifer{ilm1}'.format(ilm1=il-1),
                [well.linkobs['brine_aquifer{ilm1}'.format(ilm1=il-1)],
                 well.linkobs['brine_aquifer{il}'.format(il=il)]])
            self.add_obs_to_be_linked('mass_CO2_aquifer{il}'.format(il=il))
            self.add_obs_to_be_linked('mass_brine_aquifer{il}'.format(il=il))

        well.add_obs_to_be_linked('CO2_atm')
        well.add_obs_to_be_linked('brine_atm')
        il = num_shale_layers - 1
        self.add_kwarg_linked_to_collection(
            'CO2_aquifer{0}'.format(il),
            [well.linkobs['CO2_aquifer{0}'.format(il)], well.linkobs['CO2_atm']])
        self.add_kwarg_linked_to_collection(
            'brine_aquifer{0}'.format(il),
            [well.linkobs['brine_aquifer{0}'.format(il)], well.linkobs['brine_atm']])


if __name__ == "__main__":
    # Define logging level
    logging.basicConfig(level=logging.WARNING)

    try:
        from openiam import SimpleReservoir, MultisegmentedWellbore
    except ImportError as err:
        print('Unable to load IAM class module: '+str(err))

    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))

    # Add parameters of reservoir component model
    sres.add_par('numberOfShaleLayers', value=3, vary=False)
    sres.add_par('shale1Thickness', min=30.0, max=150., value=45.0)

    # Add observations of reservoir component model
    sres.add_obs_to_be_linked('pressure')
    sres.add_obs_to_be_linked('CO2saturation')
    # When add_obs method is called, system model automatically
    # (if no other options of the method is used) adds a list of observations
    # with name sres.obsnm_0, sres.obsnm_1,.... The final index is determined by the
    # number of time points in system model time_array attribute.
    # For more information, see the docstring of add_obs of ComponentModel class.
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')
    sres.add_obs('mass_CO2_reservoir')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))
    ms.add_par('wellRadius', min=0.01, max=0.02, value=0.015)

    # Add linked parameters: common to both components
    ms.add_par_linked_to_par('numberOfShaleLayers',
                             sres.deterministic_pars['numberOfShaleLayers'])
    ms.add_par_linked_to_par('shale1Thickness', sres.pars['shale1Thickness'])
    ms.add_par_linked_to_par('shale2Thickness',
                             sres.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('shale3Thickness',
                             sres.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('aquifer1Thickness',
                             sres.default_pars['aquiferThickness'])
    ms.add_par_linked_to_par('aquifer2Thickness',
                             sres.default_pars['aquiferThickness'])
    ms.add_par_linked_to_par('reservoirThickness',
                             sres.default_pars['reservoirThickness'])
    ms.add_par_linked_to_par('datumPressure',
                             sres.default_pars['datumPressure'])

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

    # Add observations of multisegmented wellbore component model
    ms.add_obs_to_be_linked('CO2_aquifer1')
    ms.add_obs_to_be_linked('CO2_aquifer2')
    ms.add_obs_to_be_linked('CO2_atm')
    ms.add_obs_to_be_linked('brine_aquifer1')
    ms.add_obs_to_be_linked('brine_aquifer2')
    ms.add_obs_to_be_linked('brine_atm')
    ms.add_obs('CO2_aquifer1')
    ms.add_obs('CO2_aquifer2')
    ms.add_obs('mass_CO2_aquifer1')
    ms.add_obs('mass_CO2_aquifer2')

    adapt = sm.add_component_model_object(RateToMassAdapter(name='adapt', parent=sm))
    adapt.add_kwarg_linked_to_collection(
        'CO2_aquifer1', [ms.linkobs['CO2_aquifer1'], ms.linkobs['CO2_aquifer2']])
    adapt.add_kwarg_linked_to_collection(
        'CO2_aquifer2', [ms.linkobs['CO2_aquifer2'], ms.linkobs['CO2_atm']])
    adapt.add_kwarg_linked_to_collection(
        'brine_aquifer1', [ms.linkobs['brine_aquifer1'], ms.linkobs['brine_aquifer2']])
    adapt.add_kwarg_linked_to_collection(
        'brine_aquifer2', [ms.linkobs['brine_aquifer2'], ms.linkobs['brine_atm']])
    adapt.add_obs('mass_CO2_aquifer1')
    adapt.add_obs('mass_CO2_aquifer2')
    adapt.add_obs('mass_brine_aquifer1')
    adapt.add_obs('mass_brine_aquifer2')

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically
    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

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

    CO2leakrates_aq1 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer1')
    CO2leakrates_aq2 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer2')
    print('------------------------------------------------------------------')
    print('CO2_aquifer1', CO2leakrates_aq1, sep='\n')
    print('------------------------------------------------------------------')
    print('CO2_aquifer2', CO2leakrates_aq2, sep='\n')

    CO2mass_aq1 = sm.collect_observations_as_time_series(ms, 'mass_CO2_aquifer1')
    CO2mass_aq2 = sm.collect_observations_as_time_series(ms, 'mass_CO2_aquifer2')
    adapt_CO2mass_aq1 = sm.collect_observations_as_time_series(adapt, 'mass_CO2_aquifer1')
    adapt_CO2mass_aq2 = sm.collect_observations_as_time_series(adapt, 'mass_CO2_aquifer2')
    print('------------------------------------------------------------------')
    print('Compare accumulated masses returned from wellbore component and adapter')
    print('------------------------------------------------------------------')
    print('mass_CO2_aquifer1 returned by wellbore component', CO2mass_aq1, sep='\n')
    print('mass_CO2_aquifer1 returned by adapter component', adapt_CO2mass_aq1, sep='\n')
    print('------------------------------------------------------------------')
    print('mass_CO2_aquifer2 returned by wellbore component', CO2mass_aq2, sep='\n')
    print('mass_CO2_aquifer2 returned by adapter component', adapt_CO2mass_aq2, sep='\n')
