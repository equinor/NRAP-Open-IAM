# -*- coding: utf-8 -*-
import logging
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import SystemModel, ComponentModel
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

from openiam.cfi.commons import process_parameters

class GeneralizedFlowRate(ComponentModel):
    """
    The Generalized Flow Rate component model is a model representing wide range of
    carbon dioxide (|CO2|) and brine leakage flow rates and created based
    on the results of multiple wellbore simulations. The generalized models
    facilitate the implementation of flow rates in an uncertainty quantification
    (UQ) framework since the relevant leakage rate and time parameters
    can be generated randomly. The basic shape for these models were constructed
    from the results of numerical wellbore simulations based on pressure and
    saturation profiles derived from the Kimberlina reservoir model
    :cite:`Wainwright2012` coupled with wellbore permeability to yield |CO2| and
    complimentary brine leakage functions. More details covering derivation and
    application of the model can be found in :cite:`MansoorEtAl2014`.

    In the IAM control file, the type name for the Generalized Flow Rate component is
    ``GeneralizedFlowRate``. The description of the possible component's parameters are
    provided below.

    * **numberOfShaleLayers** [-] (3 to 30) - number of shale layers in the
      system (default: 3); *linked to Stratigraphy*. The shale units must be
      separated by an aquifer.

    * **logPeakCO2Rate** [|log10| |kg/s|] (-inf to 5) - logarithm of the largest |CO2| flow rate
      (default: -5)

    * **timePeakCO2Rate** [years] (0 to 1000) - time to reach the largest |CO2| flow rate
      from initial time (default: 5)

    * **durationPeakCO2Rate** [years] (0 to 1000) - length of time period during which |CO2|
      flow rate was the largest (default: 10)

    * **durationPeakZeroCO2Rate** [years] (0 to 1000) - length of time period during which |CO2|
      flow rate decreased from the largest rate to zero (default: 100)

    * **logInitBrineRate** [|log10| |kg/s|] (-inf to 5) - logarithm of the initial brine flow
      rate (default: -10)

    * **logFinalBrineRate** [|log10| |kg/s|] (-inf to 5) - logarithm of the final brine
      flow rate (default: -11.5). Ratio of initial brine rate over final brine rate
      is recommended to be between 0.2 and 0.3

    * **durationInitBrineRate** [years] (0 to 1000) - length of initial brine flow rate
      time period (default: 2)

    * **durationInitFinalBrineRate** [years] (0 to 1000) - length of time period during
      which brine flow rate decreased from initial to final rate (default: 10)

    * **mitigationTime** [years] (0 to inf) - time at which the leakage
      was remediated (default: 10000)

    The possible outputs from the Generalized Flow Rate component are leakage rates
    of |CO2| and brine to the aquifer specified by user. The names of the
    observations are of the form:

    * **CO2_aquifer#** [|kg/s|] - |CO2| leakage rates where # is an aquifer index

    * **brine_aquifer#** [|kg/s|] - brine leakage rates

    * **mass_CO2_aquifer#** [|kg|] - mass of |CO2| leaked into the specified aquifer.

    """
    def __init__(self, name, parent, leak_to='aquifer1'):
        """
        Constructor method of GeneralizedFlowRate class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param leak_to: name of the aquifer leakage to which is to be
            simulated with the component. By default, the leakage is modeled
            to be to aquifer 1. Format of the aquifer name is 'aquifer#' where
            # is an index of aquifer, 1 is an index of the deepest aquifer.
        :type leak_to: str

        :returns: GeneralizedFlowRate class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25}   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'GeneralizedFlowRate'

        # Set default parameters of the component model
        self.add_default_par('numberOfShaleLayers', value=3)
        self.add_default_par('logPeakCO2Rate', value=-5)
        self.add_default_par('timePeakCO2Rate', value=5.0) # in years
        self.add_default_par('durationPeakCO2Rate', value=10.0)
        self.add_default_par('durationPeakZeroCO2Rate', value=100.0)
        self.add_default_par('logInitBrineRate', value=-10.0)
        self.add_default_par('logFinalBrineRate', value=-11.5)
        self.add_default_par('durationInitBrineRate', value=2.0)
        self.add_default_par('durationInitFinalBrineRate', value=10.0)
        self.add_default_par('mitigationTime', value=10000.0)

        # Define dictionary of boundaries
        # Time and duration boundaries are in days
        self.pars_bounds = dict()
        self.pars_bounds['numberOfShaleLayers'] = [3, 30]
        self.pars_bounds['logPeakCO2Rate'] = [-np.inf, 5.0]
        self.pars_bounds['timePeakCO2Rate'] = [0.0, 1000.0]    # 1000 years
        self.pars_bounds['durationPeakCO2Rate'] = [0.0, 1000.0]
        self.pars_bounds['durationPeakZeroCO2Rate'] = [0.0, 1000.0]
        self.pars_bounds['logInitBrineRate'] = [-np.inf, 5.0]
        self.pars_bounds['logFinalBrineRate'] = [-np.inf, 5.0]
        self.pars_bounds['durationInitBrineRate'] = [0.0, 1000.0]
        self.pars_bounds['durationInitFinalBrineRate'] = [0.0, 1000.0]
        self.pars_bounds['mitigationTime'] = [0.0, np.inf]

        self.sec_per_year = 365.25*24*60*60
        self.sec_per_min = 24*60*60.0

        self.leak_to = leak_to
        self.leak_layer = int(leak_to[7:])

        # Add accumulators that would keep the mass of CO2 and brine leaked
        self.add_accumulator('mass_CO2_{}'.format(self.leak_to), sim=0.0)
        self.add_accumulator('mass_brine_{}'.format(self.leak_to), sim=0.0)

        debug_msg = 'GeneralizedFlowRate component created with name {}'.format(name)
        logging.debug(debug_msg)

    def flow_rate(self, t, max_rate, t1, dt2, dt3, lmbda, mitigation_time):
        """
        Calculate flow rate of CO2 and brine based on the provided parameters.
        """
        t2 = t1 + dt2
        t3 = t1 + dt2 + dt3

        if t1 == 0:
            a0 = 0.
        else:
            a0 = max_rate/t1

        a1 = max_rate*(1.0-lmbda)/(t2-t3)
        a2 = max_rate*(lmbda*t2-t3)/(t2-t3)

        if t >= mitigation_time:
            rate = 0.0
        else:
            if t < t1:
                rate = a0 * t
            elif t1 <= t < t2:
                rate = max_rate
            elif t2 <= t < t3:
                rate = a1 * t + a2
            elif t >= t3:
                rate = lmbda * max_rate
            rate = max(rate, 0.0)

        return rate

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
                if (val < self.pars_bounds[key][0]) or (
                        val > self.pars_bounds[key][1]):
                    error_msg = ''.join([
                        'Parameter {} of GeneralizedFlowRate component {} ',
                        'is out of boundaries.']).format(key, self.name)
                    logging.error(error_msg)
                    raise ValueError(''.join([
                        "Entered parameter values do not satisfy assumptions",
                        " of the GeneralizedFlowRate component. Parameter {}",
                        " with value {} is out of boundaries. "]).format(key, val))
            else:
                warn_msg = ''.join([
                    'Parameter {} is not recognized as an input parameter ',
                    'of GeneralizedFlowRate component {}.']).format(key, self.name)
                logging.warning(warn_msg)

    def simulation_model(self, p, time_point=365.25, time_step=365.25):
        """
        Return |CO2| and brine leakage rates calculated as generalized rates.
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        out = {}
        for ind in range(1, actual_p['numberOfShaleLayers']):
            out['CO2_aquifer{}'.format(ind)] = 0.0
            out['brine_aquifer{}'.format(ind)] = 0.0
            out['mass_CO2_aquifer{}'.format(ind)] = 0.0
            out['mass_brine_aquifer{}'.format(ind)] = 0.0
        out['CO2_atm'] = 0.0
        out['brine_atm'] = 0.0

        # Calculate CO2 leakage rates
        out['CO2_{}'.format(self.leak_to)] = self.flow_rate(
            time_point/365.25,
            10**actual_p['logPeakCO2Rate']/self.sec_per_year,
            actual_p['timePeakCO2Rate'], actual_p['durationPeakCO2Rate'],
            actual_p['durationPeakZeroCO2Rate'], 0.0,
            actual_p['mitigationTime'])*self.sec_per_year
        out['brine_{}'.format(self.leak_to)] = self.flow_rate(
            time_point/365.25,
            10**actual_p['logInitBrineRate']/self.sec_per_year, 0.0,
            actual_p['durationInitBrineRate'],
            actual_p['durationInitFinalBrineRate'],
            10**(actual_p['logFinalBrineRate']-actual_p['logInitBrineRate']),
            actual_p['mitigationTime'])*self.sec_per_year

        out['mass_CO2_{}'.format(self.leak_to)] = self.accumulators[
            'mass_CO2_{}'.format(self.leak_to)].sim
        out['mass_brine_{}'.format(self.leak_to)] = self.accumulators[
            'mass_brine_{}'.format(self.leak_to)].sim

        self.accumulators['mass_CO2_{}'.format(self.leak_to)].sim = (
            self.accumulators['mass_CO2_{}'.format(self.leak_to)].sim
            + time_step*self.sec_per_min*out['CO2_{}'.format(self.leak_to)])
        self.accumulators['mass_brine_{}'.format(self.leak_to)].sim = (
            self.accumulators['mass_brine_{}'.format(self.leak_to)].sim
            + time_step*self.sec_per_min*out['brine_{}'.format(self.leak_to)])

        return out

    def connect_with_system(self, component_data, name2obj_dict,
                            *args, **kwargs):
        """
        Code to add the component to system model for control file interface.

        :param component_data: Dictionary of component data including 'Connection'
        :type component_data: dict

        :param name2obj_dict: Dictionary to map component names to objects
        :type name2obj: dict

        :returns: None
        """
        # Flag whether numberOfShaleLayers is provided as parameter
        par_present = False
        if ('Parameters' in component_data) and (component_data['Parameters']):
            if 'numberOfShaleLayers' in component_data['Parameters']:
                par_present = True

        # Process parameters data if any
        process_parameters(self, component_data, name2obj_dict)

        if not par_present:
            if 'numberOfShaleLayers' in component_data:
                self.add_par('numberOfShaleLayers',
                             value=component_data['numberOfShaleLayers'],
                             vary=False)
            else:
                strata = name2obj_dict['strata']
                if 'numberOfShaleLayers' in strata.pars:
                    self.add_par_linked_to_par(
                        'numberOfShaleLayers',
                        strata.pars['numberOfShaleLayers'])

                elif 'numberOfShaleLayers' in strata.deterministic_pars:
                    self.add_par_linked_to_par(
                        'numberOfShaleLayers',
                        strata.deterministic_pars['numberOfShaleLayers'])

                elif 'numberOfShaleLayers' in strata.default_pars:
                    self.add_par_linked_to_par(
                        'numberOfShaleLayers',
                        strata.default_pars['numberOfShaleLayers'])

                else:
                    self.add_par('numberOfShaleLayers', value=3, vary=False)

        if 'LeakTo' in component_data:
            leak_to = component_data['LeakTo'].lower()
        else:
            leak_to = 'aquifer1'
        self.leak_to = leak_to

        if 'aquifer' in leak_to:
            self.leak_layer = int(leak_to[7:])

        self.add_accumulator('mass_CO2_{}'.format(self.leak_to), sim=0.0)
        self.add_accumulator('mass_brine_{}'.format(self.leak_to), sim=0.0)

    def reset(self):
        pass

    system_params = ['numberOfShaleLayers']


if __name__ == "__main__":
    __spec__ = None

    example = 1
    if example == 1:
        logging.basicConfig(level=logging.DEBUG)
        # Define keyword arguments of the system model
        num_years = 50
        time_array = 365.25*np.arange(0.0, num_years+1)
        sm_model_kwargs = {'time_array': time_array} # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        gfr = sm.add_component_model_object(GeneralizedFlowRate(
            name='gfr', parent=sm))

        # Add parameters of generalized flow rate component model
        gfr.add_par('logPeakCO2Rate', value=-4.0, vary=False)
        gfr.add_par('timePeakCO2Rate', value=10.0, vary=False)
        gfr.add_par('durationPeakCO2Rate', value=15.0, vary=False)
        gfr.add_par('durationPeakZeroCO2Rate', value=100.0, vary=False)
        gfr.add_par('logInitBrineRate', value=-10.0, vary=False)
        gfr.add_par('logFinalBrineRate', value=-11.5, vary=False)
        gfr.add_par('durationInitBrineRate', value=2.0, vary=False)
        gfr.add_par('durationInitFinalBrineRate', value=10.0, vary=False)
        gfr.add_par('mitigationTime', value=45.0, vary=False)

        # Add observations of generalized flow rate component model
        gfr.add_obs('CO2_aquifer1')
        gfr.add_obs('brine_aquifer1')
        gfr.add_obs('mass_CO2_aquifer1')
        gfr.add_obs('mass_brine_aquifer1')

        sm.forward()  # system model is run deterministically

        print('------------------------------------------------------------------')
        print('                  Forward method illustration ')
        print('------------------------------------------------------------------')

        CO2_aquifer = sm.collect_observations_as_time_series(gfr, 'CO2_aquifer1')
        brine_aquifer = sm.collect_observations_as_time_series(gfr, 'brine_aquifer1')

        print('CO2_aquifer1', CO2_aquifer, sep='\n')
        print('------------------------------------------------------------------')
        print('brine_aquifer1', brine_aquifer, sep='\n')
        print('------------------------------------------------------------------')
        print('mass_CO2_aquifer1',
              sm.collect_observations_as_time_series(gfr, 'mass_CO2_aquifer1'),
              sep='\n')
        print('------------------------------------------------------------------')
        print('mass_brine_aquifer1',
              sm.collect_observations_as_time_series(gfr, 'mass_brine_aquifer1'),
              sep='\n')

        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(121)
        plt.plot(time_array/365.25, CO2_aquifer, color='red', linewidth=1)
        plt.xlabel('Time, t [years]', fontsize=14)
        plt.ylabel(r'Flow rate, q [kg/s]', fontsize=14)
        plt.title(r'Leakage of CO$_2$ to aquifer 1', fontsize=18)
        plt.tick_params(labelsize=12)
        plt.xlim([0, 50])
    #    ax.get_yaxis().set_label_coords(-0.13, 0.5)

        ax = fig.add_subplot(122)
        plt.plot(time_array/365.25, brine_aquifer, color='blue', linewidth=1)
        plt.xlabel('Time, t [years]', fontsize=14)
        plt.ylabel('Flow rate, q [kg/s]', fontsize=14)
        plt.title(r'Leakage of brine to aquifer 1', fontsize=18)
        plt.tick_params(labelsize=12)
        plt.xlim([0, 50])
    #    ax.get_yaxis().set_label_coords(-0.13, 0.5)
