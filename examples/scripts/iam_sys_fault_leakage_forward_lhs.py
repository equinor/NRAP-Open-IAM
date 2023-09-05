# -*- coding: utf-8 -*-
"""
Example illustrates setup of hydrocarbon leakage component.
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.sep.join(['..', '..', 'source']))

try:
    from openiam import SystemModel, FaultLeakage
except ImportError as err:
    print('Unable to load NRAP-Open-IAM modules: {}'.format(err))


if __name__ == "__main__":
    test_case = 2

    # Define keyword arguments of the system model.
    time_array = np.linspace(0.0, 30.0, num=20)*365.25 # time in days
    sm_model_kwargs = {'time_array': time_array} # time must be given in days

    # Create the system model.
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # well_index is the well location index (not to be confused with Peaceman's well index)
    params = {'damage_zone_perm': [-13.0],
              'damage_zone_por': [0.05],
              'shallow_aquifer_perm': [-13.0],
              'deep_aquifer_perm': [-12.0],
              'shallow_aquifer_por': [0.30],
              'deep_aquifer_por': [0.20],
              'well_index': [0.0],
              'well_rate': [20.0, 7.0, 25.0],
              'dip_angle': [40.0],
              'injection_time': [20.0],
              'geothermal_gradient': [35.0]}

    flc = sm.add_component_model_object(FaultLeakage(name='flc', parent=sm))

    # Add parameters.
    for key in params:
        if key == 'well_rate':
            flc.add_par(key, value=params[key][0],
                        min=params[key][1], max=params[key][2])
        else:
            flc.add_par(key, value=params[key][0], vary=False)

    # Add observations of fault leakage component
    flc_observations = ['CO2_aquifer', 'mass_CO2_aquifer',
                        'brine_aquifer', 'mass_brine_aquifer']
    for obs_nm in flc_observations:
        flc.add_obs(obs_nm)

    if test_case == 1:
        print('--------------------------------------------------------------')
        print('Started forward simulation...')
        print('--------------------------------------------------------------')
        # Run forward simulation
        sm.forward()

        # Collect ouputs into a dictionary
        outputs = {obs_nm: sm.collect_observations_as_time_series(flc, obs_nm)
                   for obs_nm in flc_observations}

        # Setup plot parameters
        label_size = 13
        font_size = 16
        ticks_size = 12
        line_width = 3

        # Plot results
        # First group of plots
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        for ind, obs_nm in enumerate(['CO2_aquifer', 'brine_aquifer']):
            ax[ind].plot(time_array/365.25, outputs[obs_nm], '-',
                         color='steelblue', linewidth=line_width)
            ax[ind].set_xlabel('Time, years', fontsize=label_size)
        ax[0].set_ylabel(r'CO$_2$ leakage rates into aquifer, [kg/s]',
                         fontsize=label_size)
        ax[1].set_ylabel(r'Brine leakage rates into aquifer, [kg/s]',
                         fontsize=label_size)
        fig.suptitle('Results of simulation', fontsize=font_size)

        # Second group of plots
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        for ind, obs_nm in enumerate(['mass_CO2_aquifer', 'mass_brine_aquifer']):
            ax[ind].plot(time_array/365.25, outputs[obs_nm]/1.0e+9, '-',
                         color='steelblue', linewidth=line_width)
            ax[ind].set_xlabel('Time, years', fontsize=label_size)
        ax[0].set_ylabel(r'Mass of CO$_2$ in aquifer, [Mt]', fontsize=label_size)
        ax[1].set_ylabel(r'Mass of brine in aquifer, [Mt]', fontsize=label_size)
        fig.suptitle('Results of simulation', fontsize=font_size)

    else:
        num_samples = 25
        ncpus = 5
        s = sm.lhs(siz=num_samples, seed=111)

        # Run model using values in samples for parameter values
        print('--------------------------------------------------------------')
        print('Started LHS based simulations...')
        print('--------------------------------------------------------------')
        s.run(cpus=ncpus, verbose=False)

        # Extract results from stochastic simulations
        out = s.collect_observations_as_time_series()

        # Setup plot parameters
        label_size = 13
        font_size = 16
        ticks_size = 12
        line_width = 3

        # Plot results
        # First group of plots
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        for ind, obs_nm in enumerate(['CO2_aquifer', 'brine_aquifer']):
            for sind in range(num_samples):
                ax[ind].plot(time_array/365.25, out['flc.'+obs_nm][sind], '-',
                             color='steelblue', linewidth=line_width)
            ax[ind].set_xlabel('Time, years', fontsize=label_size)
        ax[0].set_ylabel(r'CO$_2$ leakage rates into aquifer, [t/s]',
                         fontsize=label_size)
        ax[1].set_ylabel(r'Brine leakage rates into aquifer, [t/s]',
                         fontsize=label_size)
        fig.suptitle('Results of simulation', fontsize=font_size)

        # Second group of plots
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        for ind, obs_nm in enumerate(['mass_CO2_aquifer', 'mass_brine_aquifer']):
            for sind in range(num_samples):
                ax[ind].plot(time_array/365.25, out['flc.'+obs_nm][sind]/1.0e+9, '-',
                         color='steelblue', linewidth=line_width)
            ax[ind].set_xlabel('Time, years', fontsize=label_size)
        ax[0].set_ylabel(r'Mass of CO$_2$ in aquifer, [Mt]', fontsize=label_size)
        ax[1].set_ylabel(r'Mass of brine in aquifer, [Mt]', fontsize=label_size)
        fig.suptitle('Results of simulation', fontsize=font_size)
