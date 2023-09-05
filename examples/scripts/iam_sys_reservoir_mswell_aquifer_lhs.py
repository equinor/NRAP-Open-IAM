'''
This example couples the simple reservoir, multisegmented wellbore and
carbonate aquifer models. The saturation/pressure output produced by simple
reservoir model is used to drive leakage from a single multisegmented wellbore
model, which is passed to the input of an adapter that provides well
coordinates, |CO2| and brine leakage rates and cumulative mass fluxes to the
carbonate aquifer model. Plots of relevant observations are created.

Example of run:
$ python iam_sys_reservoir_mswell_aquifer_lhs.py
'''

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, SimpleReservoir, MultisegmentedWellbore,
                     CarbonateAquifer, RateToMassAdapter)



if __name__ == "__main__":
    # For multiprocessing in Spyder
    __spec__ = None
    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}   # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))

    # Add parameters of reservoir component model
    sres.add_par('numberOfShaleLayers', value=3, vary=False)
    sres.add_par('shale1Thickness', min=30.0, max=150., value=45.0)
    sres.add_par('aquifer1Thickness', min=100.0, max=500.0, value=300.0)
    sres.add_par('injRate', min=0.1, max=10.0, value=1.0)

    # Add observations of reservoir component model
    sres.add_obs_to_be_linked('pressure')
    sres.add_obs_to_be_linked('CO2saturation')
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')
    sres.add_obs('mass_CO2_reservoir')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))
    ms.add_par('wellRadius', min=0.01, max=0.02, value=0.015)
    ms.add_par('logWellPerm', min=-14.0, max=-12.0, value=-13.0)

    # Add linked parameters: common to reservoir and wellbore components
    ms.add_par_linked_to_par('numberOfShaleLayers',
                             sres.deterministic_pars['numberOfShaleLayers'])
    ms.add_par_linked_to_par('shale1Thickness', sres.pars['shale1Thickness'])
    ms.add_par_linked_to_par('shale2Thickness',
                             sres.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('shale3Thickness',
                             sres.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('aquifer1Thickness',
                             sres.pars['aquifer1Thickness'])
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
    ms.add_obs_to_be_linked('brine_aquifer1')
    ms.add_obs_to_be_linked('brine_aquifer2')
    ms.add_obs_to_be_linked('mass_CO2_aquifer1')
    ms.add_obs_to_be_linked('mass_CO2_aquifer2')
    ms.add_obs_to_be_linked('brine_atm')
    ms.add_obs_to_be_linked('CO2_atm')
    ms.add_obs('brine_aquifer1')
    ms.add_obs('CO2_aquifer1')

    # Add adapter that transforms leakage rates to accumulated mass
    adapt = sm.add_component_model_object(RateToMassAdapter(name='adapt', parent=sm))
    adapt.add_kwarg_linked_to_collection('CO2_aquifer1',
        [ms.linkobs['CO2_aquifer1'], ms.linkobs['CO2_aquifer2']])
    adapt.add_kwarg_linked_to_collection('CO2_aquifer2',
        [ms.linkobs['CO2_aquifer2'], ms.linkobs['CO2_atm']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer1',
        [ms.linkobs['brine_aquifer1'], ms.linkobs['brine_aquifer2']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer2',
        [ms.linkobs['brine_aquifer2'], ms.linkobs['brine_atm']])
    adapt.add_obs_to_be_linked('mass_CO2_aquifer1')
    adapt.add_obs_to_be_linked('mass_CO2_aquifer2')
    adapt.add_obs_to_be_linked('mass_brine_aquifer1')
    adapt.add_obs_to_be_linked('mass_brine_aquifer2')
    adapt.add_obs('mass_CO2_aquifer1')
    adapt.add_obs('mass_brine_aquifer1')

    # Add carbonate aquifer model object and define parameters
    ca = sm.add_component_model_object(CarbonateAquifer(name='ca', parent=sm))
    ca.add_par('perm_var', min=0.017, max=1.89, value=0.9535)
    ca.add_par('corr_len', min=1.0, max=3.95, value=2.475)
    ca.add_par('aniso', min=1.1, max=49.1, value=25.1)
    ca.add_par('mean_perm', min=-13.8, max=-10.3, value=-12.05)
    ca.add_par_linked_to_par('aqu_thick', sres.pars['aquifer1Thickness'])
    ca.add_par('hyd_grad', min=2.88e-4, max=1.89e-2, value=9.59e-3)
    ca.add_par('calcite_ssa', min=0, max=1.0e-2, value=5.5e-03)
    ca.add_par('organic_carbon', min=0, max=1.0e-2, value=5.5e-03)
    ca.add_par('benzene_kd', min=1.49, max=1.73, value=1.61)
    ca.add_par('benzene_decay', min=0.15, max=2.84, value=1.5)
    ca.add_par('nap_kd', min=2.78, max=3.18, value=2.98)
    ca.add_par('nap_decay', min=-0.85, max=2.04, value=0.595)
    ca.add_par('phenol_kd', min=1.21, max=1.48, value=1.35)
    ca.add_par('phenol_decay', min=-1.22, max=2.06, value=0.42)
    ca.add_par('cl', min=0.1, max=6.025, value=0.776)
    ca.model_kwargs['x'] = [0.]
    ca.model_kwargs['y'] = [0.]

    CO2_rate_obs_list = []
    brine_rate_obs_list = []
    CO2_mass_obs_list = []
    brine_mass_obs_list = []
    CO2_rate_obs_list.append(ms.linkobs['CO2_aquifer1'])
    brine_rate_obs_list.append(ms.linkobs['brine_aquifer1'])
    CO2_mass_obs_list.append(adapt.linkobs['mass_CO2_aquifer1'])
    brine_mass_obs_list.append(adapt.linkobs['mass_brine_aquifer1'])

    # Add aquifer component's keyword argument co2_rate linked to the collection created above
    ca.add_kwarg_linked_to_collection('co2_rate', CO2_rate_obs_list)
    ca.add_kwarg_linked_to_collection('brine_rate', brine_rate_obs_list)
    ca.add_kwarg_linked_to_collection('co2_mass', CO2_mass_obs_list)
    ca.add_kwarg_linked_to_collection('brine_mass', brine_mass_obs_list)

    # Add observations (output) from the carbonate aquifer model
    ca.add_obs('pH_volume')
    ca.add_obs('TDS_volume')

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically
    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    print('Pressure',
          sm.collect_observations_as_time_series(sres, 'pressure'), sep='\n')
    print('------------------------------------------------------------------')
    print('CO2saturation',
          sm.collect_observations_as_time_series(sres, 'CO2saturation'), sep='\n')
    print('------------------------------------------------------------------')
    print('CO2_aquifer1',
          sm.collect_observations_as_time_series(ms, 'CO2_aquifer1'), sep='\n')
    print('------------------------------------------------------------------')
    print('brine_aquifer1',
          sm.collect_observations_as_time_series(ms, 'brine_aquifer1'), sep='\n')
    print('------------------------------------------------------------------')
    print('mass_CO2_aquifer1',
          sm.collect_observations_as_time_series(adapt, 'mass_CO2_aquifer1'), sep='\n')
    print('------------------------------------------------------------------')
    print('mass_brine_aquifer1',
          sm.collect_observations_as_time_series(adapt, 'mass_brine_aquifer1'), sep='\n')
    print('------------------------------------------------------------------')
    print('pH',
          sm.collect_observations_as_time_series(ca, 'pH_volume'), sep='\n')
    print('------------------------------------------------------------------')
    print('TDS',
          sm.collect_observations_as_time_series(ca, 'TDS_volume'), sep='\n')

    print('------------------------------------------------------------------')
    print('                          UQ illustration ')
    print('------------------------------------------------------------------')

    import random
    num_samples = 25
    ncpus = 1
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=random.randint(500, 1100))

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)

    # sm.obs_base_names:
    # Observations: 'sres.pressure', 'sres.CO2saturation', 'sres.mass_CO2_reservoir',
    # 'ms.brine_aquifer1', 'ms.CO2_aquifer1', 'adapt.mass_CO2_aquifer1',
    # 'adapt.mass_brine_aquifer1', 'ca.pH', 'ca.TDS'

    # Plot results
    num_obs = len(time_array)
    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(221)
    for j in range(num_samples):
        plt.semilogy(time_array[0:]/365.25, results[j, (5*num_obs):(6*num_obs)],
                     color='#1B4F72', linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel(r'Mass of CO$_2$ leaked (kg)', fontsize=14)
    plt.title(r'Mass of CO$_2$ leaked into aquifer 1', fontsize=18)
    plt.tick_params(labelsize=12)

    ax = fig.add_subplot(223)
    for j in range(num_samples):
        plt.semilogy(time_array[0:]/365.25, results[j, (6*num_obs):(7*num_obs)],
                     color='#000066', linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel(r'Mass of brine leaked (kg)', fontsize=14)
    plt.title(r'Mass of brine leaked into aquifer 1', fontsize=18)
    plt.tick_params(labelsize=12)

    ax = fig.add_subplot(222)
    for j in range(num_samples):
        plt.semilogy(time_array[0:]/365.25, results[j, (7*num_obs):(8*num_obs)],
                     color='#663300', linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel('pH Volume (m$^3$)', fontsize=14)
    plt.title(r'Volume of aquifer 1 below pH threshold', fontsize=18)
    plt.tight_layout()
    plt.tick_params(labelsize=12)

    ax = fig.add_subplot(224)
    for j in range(num_samples):
        plt.plot(time_array[0:]/365.25, results[j, (8*num_obs):(9*num_obs)],
                 color='#CC9900', linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel('TDS Volume (m$^3$)', fontsize=14)
    plt.title(r'Volume of aquifer 1 above TDS threshold', fontsize=18)
    plt.tight_layout()
    plt.tick_params(labelsize=12)

    to_save = False
    if to_save:
        plt.savefig('aquifer_example01_plot.png', dpi=300)
