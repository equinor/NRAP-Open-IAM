# -*- coding: utf-8 -*-
"""
Methods developed for the implementation of different setups
of the risk configurer component.

Last modified: June, 2023

@author: Veronika Vasylkivska, Greg Lackey
"""
import sys
import os

import numpy as np
import pandas as pd
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# Lists of probability options
DETECTION_PROBABILITY_OPTIONS = ['always', 'never', 'likely', 'possibly']

LEAKAGE_PROBABILITY_OPTIONS = ['always', 'never', 'likely', 'possibly',
                               'score_based', 'simple_score_based']


def create_details(num_leak_types):
    """
    Create configuration of risk details associated with each well

    """
    details = {'risk_type': 'TBD',
               'location': None,
               # Code for the conditions met; 1st entry is for leak type and exist path,
               # all other for other leakage types
               'num_leak_types': num_leak_types,
               'cond_indicator': np.zeros(num_leak_types),
               'run_freq_init': 0,              # 0 - well is not leaking, 2 - well is leaking
               # Leakage risk details
               'leak_prob': None,
               'leak_start_time': None,
               # Blowout condition details
               # pressure value threshold at which blowout can occur
               'blowout_pres': None,
               # probability that blowout will occur if the condition is met
               'blowout_prob': None,
               # time at which blowout condition is met for a first time
               'blowout_cond_time': None,
               # Corrosion condition details
               # CO2 saturation value threshold at which corrosion starts
               'corrosion_sat': None,
               'corrosion_prob': None,
               # time point (in days) at which CO2 reached the well
               'corrosion_cond_time': None,
               # Mechanical stress conditions
               'stress_pres': None,
               'stress_prob': None,
               'stress_cond_time': None,
               # Leakage detection details
               # probability that leakage will be detected
               'leak_detection_prob': None,
               # time of first leakage detected
               'leak_detected_time': None,
               # Remediation details
               # probability of the well being remediated after the leakage was detected
               'remediation_prob': None,
               # time at which well was remediated after leakage was detected
               'remediation_time': None,
               # Well status
               'well_status':None,
               # Leakage risk scores
               'deviation_score': None,
               'scvf_gmdeviation_score': None,
               'suspension_status_score': None,
               'cement_top_score': None,
               'casing_failure_score': None,
               'age_score': None,
               'depth_score': None,
               'type_score': None,
               'leak_score': None,
               'normalized_leak_score': None,
               'logWellPerm': None
              }

    return details


def reset_details(cmpnt):
    """
    Reset details of the wellbore component to the initial values.

    :param cmpnt: system model component for which details are to be reset
    :type cmpnt: ComponentModel class object
    """
    # Set all times relevant to non-leaking wells to None
    times = ['leak_start_time',
             'leak_detected_time',
             'blowout_cond_time',
             'corrosion_cond_time',
             'stress_cond_time']

    for t in times:
        cmpnt.details[t] = None

    cmpnt.details['risk_type'] = 'fixed'
    cmpnt.details['cond_indicator'] = np.zeros(cmpnt.details['num_leak_types'])

    # The wellore component under control will not be run
    cmpnt.run_frequency = 0


def update_risk_type(cmpnt, risk_type, leak_start_time):
    """
    Update relevant details of the wellbore component when it starts to leak.

    :param cmpnt: instance of component class for which details are to be updated
    :type cmpnt: ComponentModel class object
    """
    cmpnt.details['risk_type'] = risk_type
    cmpnt.details['leak_start_time'] = leak_start_time
    cmpnt.details['remediation_time'] = None
    cmpnt.run_frequency = 2


def always_detected(*args, **kwargs):
    """
    Return probability of detection of 1.
    """
    return 1.0


def never_detected(*args, **kwargs):
    """
    Return probability of detection of 0.
    """
    return 0.0


def possibly_detected(cmpnt, time, *args, **kwargs):
    """
    Return probability of detection which value depends on additional parameters
    and is based on predetermined formula.
    """
    # Get value from the dictionary kwargs; default 0.1
    scalar = kwargs.get('scalar', 0.1)

    if cmpnt.details['risk_type'] == 'blowout':
        return 1.0

    if cmpnt.details['leak_detected_time'] is not None:
        return 1.0

    if cmpnt.details['leak_start_time'] is not None:
        delta_t = time - cmpnt.details['leak_start_time']
    else:
        delta_t = 0.0

    return 1. - np.exp(-delta_t/365.25*scalar)


def likely_detected(*args, **kwargs):
    # Get value from the dictionary kwargs; default 0.75
    p = kwargs.get('value', 0.75)

    return p


def detection_prob(option, **keyword_args):
    """
    Method returns function evaluating probability of leakage detection

    The general signature of the returned function is f(*args, **kwargs).
    It is assumed that cmpnt and time are passed to this function as arguments
    and all other parameters should be sent as a list of arguments or
    # keyword arguments, i.e., the returned function can also have
    # a signature f(cmpnt, time, *args, **kwargs).
    Configurer calling the returned function will assume this structure.
    """

    prob_function = {'always': always_detected,
                     'never': never_detected,
                     'possibly': possibly_detected,
                     'likely': likely_detected}

    # Check which of the options were chosen
    if option not in DETECTION_PROBABILITY_OPTIONS:
        raise ValueError('Unknown option was passed to the detection_prob method.')

    return prob_function[option]


def remediation_prob(cmpnt, time, **keyword_args):
    """
    Evaluate probability of well remediation.
    """
    prob = 0.99
    return prob


def always_leaking(*args, **kwargs):
    return 1.0


def never_leaking(*args, **kwargs):
    return 0.0


def likely_leaking(*args, **kwargs):
    # Get value from the dictionary kwargs; default 0.75
    p = kwargs.get('value', 0.75)
    return p


def possibly_leaking(cmpnt, time, *args, **kwargs):
    """
    Return probability of leakage calculated based on predetermined formula:
        prob_of_leakage = 1 - exp(-dt/365.25*scalar)

        dt is difference between the current time and time at which the leak is
        started.
        scalar is a scaling parameter controling the behaviour of the probability
        function.
    """
    # Get value from the dictionary kwargs; default 0.1
    scalar = kwargs.get('scalar', 0.1)

    if cmpnt.details['leak_start_time'] is not None:
        delta_t = time - cmpnt.details['leak_start_time']
    else:
        delta_t = 0.0

    return 1. - np.exp(-delta_t/365.25*scalar)


def simple_score_based(*args, **kwargs):
    """
    Return probability of leakage calculated based on leakage score and initial
    probability using linear interpolation.
    """
    base_p = kwargs.get('value', 0.75)

    score = kwargs.get('score', 0.5)

    return np.interp(score, [0.2, 0.5, 1], [base_p/2, base_p, base_p*2])


def score_based(cmpnt, time, pressure, *args, **kwargs):
    """
    Return probability of leakage calculated based on leakage score and some
    additional parameters using one of the available functions.
    """
    # Check whether the component's parent defines the time array
    if cmpnt._parent.time_array is not None:
        time_array = cmpnt._parent.time_array

    base_prob = kwargs.get('base_prob', 0.5)

    adj_factor = kwargs.get('adj_factor', 2.0)

    tot_prob = cmpnt.details['normalized_leak_score']*(base_prob*adj_factor)

    if tot_prob >= 1:
        tot_prob = 0.99

    num_time_points = len(time_array)
    if num_time_points > 1:
        fun_type = kwargs.get('fun_type', 'exponential')

        if fun_type == 'point':
            prob = 1. - np.exp(np.log(1-base_prob)/num_time_points)

        elif fun_type == 'constant':
            prob = tot_prob/num_time_points

        elif fun_type == 'linear':
            intercept = kwargs.get('intercept', 0)
            rate = 2*(tot_prob-intercept*time_array[-1])/time_array[-1]**2
            prob = rate*time + intercept

        elif fun_type == 'exponential':
            rate = kwargs.get('rate', 1.0e-6)
            # Calculate the constant appropriate for tot_prob
            const = tot_prob/(time_array[-1] + 1/rate*(np.exp(-rate*time_array[-1])-1))
            prob = const*(1. - np.exp(-rate*time))
    else:
        prob = tot_prob

    return prob


def leak_prob(option, **keyword_args):
    """
    Method returns function evaluating probability of leakage
    after the right condition is met.

    The general signature of the returned function is f(*args, **kwargs).
    It is assumed that cmpnt, time, pressure, saturation will be passed
    to this function as arguments and all other parameters will be sent
    as a list of arguments and/or keyword arguments, i.e.
    the returned function can also have a signature
    f(cmpnt, time, pressure, saturation, *args, **kwargs).
    Configurer calling the returned function will assume this structure.
    """
    prob_function = {'always': always_leaking,
                     'never': never_leaking,
                     'possibly': possibly_leaking,
                     'likely': likely_leaking,
                     'score_based': score_based,
                     'simple_score_based': simple_score_based}

    # Check which of the options were chosen
    if option not in LEAKAGE_PROBABILITY_OPTIONS:
        raise ValueError('Unknown option was passed to the leak_prob method.')

    return prob_function[option]


def leakage_risk_configuration(num_wells, num_leak_types, leak_prob, well_data=None):
    """
    Return dictionary of well properties describing risk of leakage.

    Method illustrates one of the many ways to assign well properties.
    """
    details = []

    decision_var = np.random.rand(num_wells)

    # Setup an example of hypothetical well leakage ditribution
    for ind in range(num_wells):
        new_entry = create_details(num_leak_types)

        # The following setup can depend on the properties of a particular well
        # Some of the wells might have their own values of pressure/saturation thresholds
        # For other wells the thresholds are defined by configurer model parameter
        if decision_var[ind] < 0.25:
            new_entry['blowout_pres'] = 39.5e+6
            new_entry['corrosion_sat'] = 0.2
            new_entry['stress_pres'] = 37.0e+6

        # The probability functions might also be defined differently for each well
        new_entry['leak_prob'] = leak_prob()
        new_entry['blowout_prob'] = leak_prob('blowout')
        new_entry['corrosion_prob'] = leak_prob('corrosion')
        new_entry['stress_prob'] = leak_prob('stress')
        new_entry['leak_detection_prob'] = detection_prob('always')
        new_entry['remediation_prob'] = remediation_prob(1)

        if isinstance(well_data, pd.DataFrame):
            new_entry['well_status'] = well_data['lateral_status'][ind]
            new_entry['deviation_score'] = well_data['deviation_score'][ind]
            new_entry['scvf_gmdeviation_score'] = well_data['scvf_gm_score'][ind]
            new_entry['suspension_status_score'] = well_data['suspension_status_score'][ind]
            new_entry['cement_top_score'] = well_data['cement_top_score'][ind]
            new_entry['casing_failure_score'] = well_data['casing_failure_score'][ind]
            new_entry['age_score'] = well_data['age_score'][ind]
            new_entry['depth_score'] = well_data['depth_score'][ind]
            new_entry['type_score'] = well_data['type_score'][ind]
            new_entry['leak_score'] = well_data['leak_score'][ind]
            new_entry['normalized_leak_score'] = well_data['normalized_leak_score'][ind]
            new_entry['logWellPerm'] = -17.0

        details.append(new_entry)

    return details


if __name__ == "__main__":
    # Check detection probabilities
    prob1 = detection_prob(option='always')
    print(prob1())

    prob2 = detection_prob(option='never')
    print(prob2())

    prob3 = detection_prob(option='likely')
    print(prob3(value=0.6))
