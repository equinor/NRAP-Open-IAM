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
    from openiam import SystemModel, HydrocarbonLeakage
except ImportError as err:
    print('Unable to load NRAP-Open-IAM modules: {}'.format(err))


if __name__ == "__main__":
    test_case = 2
    # Define keyword arguments of the system model
    num_years = 410
    t0 = 0.0  # initial time point
    time_array = 365.25*np.arange(t0, t0+num_years+1, 10)

    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add hydrocarbon leakage component
    hcl_comp = sm.add_component_model_object(HydrocarbonLeakage(
        name='hcl', parent=sm))

    # Add parameters of hydrocarbon component model
    hcl_comp.add_par('reservoirDepth', value=1046.88, vary=False)
    hcl_comp.add_par('NTG', value=0.444, vary=False)
    hcl_comp.add_par('logResPerm', value=-13.2, min=-13.9, max=-13.1)
    hcl_comp.add_par('reservoirPressureMult', value=1.13, vary=False)
    hcl_comp.add_par('logWellPerm', value=-12.44, vary=False)
    hcl_comp.add_par('avgWaterSaturation', value=0.719, vary=False)
    hcl_comp.add_par('FCO2', value=0.647, vary=False)
    hcl_comp.add_par('FC1', value=0.062, vary=False)
    hcl_comp.add_par('FC4', value=0.082, vary=False)
    hcl_comp.add_par('FC7Plus', value=0.209, vary=False)

    # Add observations (output) of the hydrocarbon leakage component
    hcl_observations = ['mass_CO2_aquifer', 'mass_CO2_gas_aquifer',
                        'mass_methane_oil_aquifer', 'mass_methane_gas_aquifer',
                        'mass_oil_aquifer', 'mass_gas_aquifer']
    for obs_nm in hcl_observations:
        hcl_comp.add_obs(obs_nm)

    if test_case == 1:
        print('--------------------------------------------------------------')
        print('Started forward simulation...')
        print('--------------------------------------------------------------')
        # Run forward simulation
        sm.forward()

        # Collect ouputs into a dictionary
        outputs = {obs_nm: sm.collect_observations_as_time_series(hcl_comp, obs_nm)
                   for obs_nm in hcl_observations}

        # Setup plot parameters
        label_size = 13
        font_size = 16
        ticks_size = 12
        line_width = 3

        # Plot results
        # First group of plots
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        for ind, obs_nm in enumerate(['mass_CO2_aquifer', 'mass_CO2_gas_aquifer']):
            ax[ind].plot(time_array/365.25, outputs[obs_nm]/1.0e+6, '-',
                         color='steelblue', linewidth=line_width)
            ax[ind].set_xlabel('Time, years', fontsize=label_size)
        ax[0].set_ylabel(r'Mass of liquid CO$_2$ in aquifer, [Kt]', fontsize=label_size)
        ax[1].set_ylabel(r'Mass of gas CO$_2$ in aquifer, [Kt]', fontsize=label_size)
        fig.suptitle('Results of simulation', fontsize=font_size)

        # Second group of plots
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        for ind, obs_nm in enumerate(['mass_methane_oil_aquifer', 'mass_methane_gas_aquifer']):
            ax[ind].plot(time_array/365.25, outputs[obs_nm]/1.0e+6, '-',
                         color='steelblue', linewidth=line_width)
            ax[ind].set_xlabel('Time, years', fontsize=label_size)
        ax[0].set_ylabel('Mass of methane in oil phase in aquifer, [Kt]', fontsize=label_size)
        ax[1].set_ylabel('Mass of methane in aquifer, [Kt]', fontsize=label_size)
        fig.suptitle('Results of simulation', fontsize=font_size)

        # Third group of plots
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        for ind, obs_nm in enumerate(['mass_oil_aquifer', 'mass_gas_aquifer']):
            ax[ind].plot(time_array/365.25, outputs[obs_nm]/1.0e+6, '-',
                         color='steelblue', linewidth=line_width)
            ax[ind].set_xlabel('Time, years', fontsize=label_size)
        ax[0].set_ylabel('Mass of oil leaked in aquifer, [Kt]', fontsize=label_size)
        ax[1].set_ylabel('Mass of gas leaked in aquifer, [Kt]', fontsize=label_size)
        fig.suptitle('Results of simulation', fontsize=font_size)

    else:
        num_samples = 25
        ncpus = 4
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
        for ind, obs_nm in enumerate(['mass_CO2_aquifer', 'mass_CO2_gas_aquifer']):
            for sind in range(num_samples):
                ax[ind].plot(time_array/365.25, out['hcl.'+obs_nm][sind]/1.0e+6, '-',
                             color='steelblue', linewidth=line_width)
            ax[ind].set_xlabel('Time, years', fontsize=label_size)
        ax[0].set_ylabel(r'Mass of liquid CO$_2$ in aquifer, [Kt]', fontsize=label_size)
        ax[1].set_ylabel(r'Mass of gas CO$_2$ in aquifer, [Kt]', fontsize=label_size)
        fig.suptitle('Results of simulation', fontsize=font_size)

        # Second group of plots
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        for ind, obs_nm in enumerate(['mass_methane_oil_aquifer', 'mass_methane_gas_aquifer']):
            for sind in range(num_samples):
                ax[ind].plot(time_array/365.25, out['hcl.'+obs_nm][sind]/1.0e+6, '-',
                         color='steelblue', linewidth=line_width)
            ax[ind].set_xlabel('Time, years', fontsize=label_size)
        ax[0].set_ylabel('Mass of methane in oil phase in aquifer, [Kt]', fontsize=label_size)
        ax[1].set_ylabel('Mass of methane in aquifer, [Kt]', fontsize=label_size)
        fig.suptitle('Results of simulation', fontsize=font_size)

        # Third group of plots
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        for ind, obs_nm in enumerate(['mass_oil_aquifer', 'mass_gas_aquifer']):
            for sind in range(num_samples):
                ax[ind].plot(time_array/365.25, out['hcl.'+obs_nm][sind]/1.0e+6, '-',
                         color='steelblue', linewidth=line_width)
            ax[ind].set_xlabel('Time, years', fontsize=label_size)
        ax[0].set_ylabel('Mass of oil leaked in aquifer, [Kt]', fontsize=label_size)
        ax[1].set_ylabel('Mass of gas leaked in aquifer, [Kt]', fontsize=label_size)
        fig.suptitle('Results of simulation', fontsize=font_size)
