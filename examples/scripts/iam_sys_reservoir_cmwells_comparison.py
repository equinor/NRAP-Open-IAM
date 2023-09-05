'''
This example illustrate comparison of results obtained from linking reservoir
component to two versions of the cemented wellbore component. Due to the narrower
ranges of input parameters of CementedWellbore component class, the parameters are
setup to stay within those ranges for both versions of the wellbore components.

Examples of run:
$ python iam_sys_reservoir_cmwells_comparison.py
'''

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, SimpleReservoir, CementedWellbore, CementedWellboreWR


if __name__=='__main__':
    # Create system model
    num_years = 50.
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}   # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))

    # Add parameters of reservoir component model
    sres.add_par('numberOfShaleLayers', value=3, vary=False)
    sres.add_par('shale1Thickness', min=500.0, max=550., value=535.0)
    sres.add_par('shale2Thickness', min=450.0, max=500., value=475.0)
    # Shale 3 has a fixed thickness of 11.2 m
    sres.add_par('shale3Thickness', value=11.2, vary=False)
    # Aquifer 1 (thief zone has a fixed thickness of 22.4)
    sres.add_par('aquifer1Thickness', value=22.4, vary=False)
    # Aquifer 2 (shallow aquifer) has a fixed thickness of 19.2
    sres.add_par('aquifer2Thickness', value=19.2, vary=False)
    # Reservoir has a fixed thickness of 51.2
    sres.add_par('reservoirThickness', value=51.2, vary=False)

    # Add observations of reservoir component model
    for key in ['pressure', 'CO2saturation']:
        sres.add_obs(key)
        sres.add_obs_to_be_linked(key)

    # Add cemented wellbore components
    cw = {1: sm.add_component_model_object(CementedWellbore(name='cw1',
                                                            parent=sm)),
          2: sm.add_component_model_object(CementedWellboreWR(name='cw2',
                                                              parent=sm))}

    # Add similar parameters of cemented wellbore components
    for ind in [1, 2]:
        cw[ind].add_par('logWellPerm', min=-13.5, max=-11.0, value=-12.5)
        cw[ind].add_par('logThiefPerm', min=-13.5, max=-12.1, value=-13.0)

    # Add parameters specific to CementedWellboreWR
    cw[2].add_par_linked_to_par('thiefZoneThickness',
                                sres.deterministic_pars['aquifer2Thickness'])
    cw[2].add_par_linked_to_par('aquiferThickness',
                                sres.deterministic_pars['aquifer1Thickness'])
    cw[2].add_par_linked_to_par('reservoirThickness',
                                sres.deterministic_pars['reservoirThickness'])

    # Add keyword arguments of the cemented wellbore components
    for ind in [1, 2]:
        for key in ['pressure', 'CO2saturation']:
            cw[ind].add_kwarg_linked_to_obs(key, sres.linkobs[key])

    # Add composite parameters of cemented wellbore component
        cw[ind].add_composite_par(
            'wellDepth', expr=' + '.join([
                'sres.shale1Thickness', 'sres.shale2Thickness',
                'sres.shale3Thickness', 'sres.aquifer1Thickness',
                'sres.aquifer2Thickness']))
        cw[ind].add_composite_par('depthRatio',
            expr=''.join([
                '(sres.shale2Thickness + sres.shale3Thickness',
                '+ sres.aquifer2Thickness + sres.aquifer1Thickness/2)',
                '/cw{}.wellDepth']).format(ind))
        cw[ind].add_composite_par('initPressure',
            expr=''.join([
                'sres.datumPressure + cw{ind}.wellDepth*cw{ind}.g',
                '*sres.brineDensity']).format(ind=ind))

        # Add observations of the cemented wellbore components
        for key in ['CO2_aquifer1', 'CO2_aquifer2', 'CO2_atm',
                    'brine_aquifer1', 'brine_aquifer2']:
            cw[ind].add_obs(key)

    # Run forward simulation
    sm.forward()

    # Collect observations
    res_obs = {}
    for key in ['pressure', 'CO2saturation']:
        res_obs[key] = sm.collect_observations_as_time_series(sres, key)

    CO2_leakrates_aq1 = {}
    CO2_leakrates_aq2 = {}
    brine_leakrates_aq1 = {}
    brine_leakrates_aq2 = {}
    for ind in [1, 2]:
        CO2_leakrates_aq1[ind] = sm.collect_observations_as_time_series(
            cw[ind], 'CO2_aquifer1')
        CO2_leakrates_aq2[ind] = sm.collect_observations_as_time_series(
            cw[ind], 'CO2_aquifer2')

        brine_leakrates_aq1[ind] = sm.collect_observations_as_time_series(
            cw[ind], 'brine_aquifer1')
        brine_leakrates_aq2[ind] = sm.collect_observations_as_time_series(
            cw[ind], 'brine_aquifer2')

    colors = ['b', 'g']
    # Plot CO2 and brine leakage rates along the wellbore
    plt.figure(1)
    for ind in [1, 2]:
        plt.semilogy(sm.time_array/365.25, CO2_leakrates_aq1[ind],
                     color=colors[ind-1], linewidth=2, label="cw {}".format(ind))
    plt.legend()
    plt.xlabel('Time, t [years]')
    plt.ylabel('Leakage rates, q [kg/s]')
    plt.title(r'Leakage of CO$_2$ to thief zone')

    plt.figure(2)
    for ind in [1, 2]:
        plt.semilogy(sm.time_array/365.25, CO2_leakrates_aq2[ind],
                     color=colors[ind-1], linewidth=2, label="cw {}".format(ind))
    plt.legend()
    plt.xlabel('Time, t [years]')
    plt.ylabel('Leakage rates, q [kg/s]')
    plt.title(r'Leakage of CO$_2$ to aquifer')

    # Plot  brine leakage rates along the wellbore
    plt.figure(3)
    for ind in [1, 2]:
        plt.semilogy(sm.time_array/365.25, brine_leakrates_aq1[ind],
                     color=colors[ind-1], linewidth=2, label="cw {}".format(ind))
    plt.legend()
    plt.xlabel('Time, t [years]')
    plt.ylabel('Leakage rates, q [kg/s]')
    plt.title(r'Leakage of brine to thief zone')

    plt.figure(4)
    for ind in [1, 2]:
        plt.semilogy(sm.time_array/365.25, brine_leakrates_aq2[ind],
                     color=colors[ind-1], linewidth=2, label="cw {}".format(ind))
    plt.legend()
    plt.xlabel('Time, t [years]')
    plt.ylabel('Leakage rates, q [kg/s]')
    plt.title(r'Leakage of brine to aquifer')

    # Plot pressure
    plt.figure(5)
    plt.plot(sm.time_array/365.25, res_obs['pressure']/1.0e+6,
             color='k', linewidth=2)
    plt.xlabel('Time, t [years]')
    plt.ylabel('Pressure, P [MPa]')

    # Plot saturation
    plt.figure(6)
    plt.plot(sm.time_array/365.25, res_obs['CO2saturation'],
             color='k', linewidth=2)
    plt.xlabel('Time, t [years]')
    plt.ylabel(r'CO$_2$ saturation, S [-]')

    # Plot comparison chart
    # CO2 leakage rates
    plt.figure(7)
    plt.loglog(CO2_leakrates_aq1[1], CO2_leakrates_aq1[2], 'ok', linewidth=2)
    plt.xlabel('Leakage rates (cw1), q [kg/s]')
    plt.ylabel('Leakage rates (cw2), q [kg/s]')
    plt.title(r'Leakage of CO$_2$ to thief zone')

    plt.figure(8)
    plt.loglog(CO2_leakrates_aq2[1], CO2_leakrates_aq2[2], 'ok', linewidth=2)
    plt.xlabel('Leakage rates (cw1), q [kg/s]')
    plt.ylabel('Leakage rates (cw2), q [kg/s]')
    plt.title(r'Leakage of CO$_2$ to aquifer')

    # Brine leakage rates
    plt.figure(9)
    plt.loglog(brine_leakrates_aq1[1], brine_leakrates_aq1[2], 'ok', linewidth=2)
    plt.xlabel('Leakage rates (cw1), q [kg/s]')
    plt.ylabel('Leakage rates (cw2), q [kg/s]')
    plt.title('Leakage of brine to thief zone')

    plt.figure(10)
    plt.loglog(brine_leakrates_aq2[1], brine_leakrates_aq2[2], 'ok', linewidth=2)
    plt.xlabel('Leakage rates (cw1), q [kg/s]')
    plt.ylabel('Leakage rates (cw2), q [kg/s]')
    plt.title('Leakage of brine to aquifer')
