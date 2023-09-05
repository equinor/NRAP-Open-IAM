'''
This example couples the simple reservoir, multisegmented wellbore and
alluvium aquifer models. The saturation/pressure output produced by simple
reservoir model is used to drive leakage from a single multisegmented wellbore
model, which is passed to the input of an adapter that provides well
coordinates, |CO2| and brine leakage rates and cumulative mass fluxes to the
alluvium aquifer model. Plots of relevant observations are created.

Example of run:
$ python iam_sys_reservoir_mswell_aquifer_lhs.py
'''

import sys,os
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
import numpy as np
from openiam import SystemModel, SimpleReservoir
from openiam import MultisegmentedWellbore, AlluviumAquifer, RateToMassAdapter
import matplotlib.pyplot as plt


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
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm,
                                                         locX=400, locY=400))

    # Add parameters of reservoir component model
    sres.add_par('numberOfShaleLayers', value=3, vary=False)
    sres.add_par('shale1Thickness', value=35.0, vary=False)
    sres.add_par('injRate', value=8.5, min=8.5, max=9)

    # Add observations of reservoir component model
    sres.add_obs_to_be_linked('pressure')
    sres.add_obs_to_be_linked('CO2saturation')
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')
    sres.add_obs('mass_CO2_reservoir')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))
    ms.add_par('wellRadius', value=0.02, vary=False)
    ms.add_par('logWellPerm', value=-12.0, vary=False)

    # Add linked parameters: common to reservoir and wellbore components
    ms.add_par_linked_to_par('numberOfShaleLayers',
                             sres.deterministic_pars['numberOfShaleLayers'])
    ms.add_par_linked_to_par('shale1Thickness',
                             sres.deterministic_pars['shale1Thickness'])
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

    # Create alluvium aquifer component
    aa = sm.add_component_model_object(AlluviumAquifer(name='aa', parent=sm))
    # Add parameters of the alluvium aquifer component
    aa.add_par('sandFraction', value=0.5, min=0.5, max=0.55)
    aa.add_par('correlationLengthX', value=360.050, min=360, max=370)
    aa.add_par('correlationLengthZ', value=18.000, vary=False)
    aa.add_par('permeabilityClay', value=-17.000, vary=False)
    aa.add_par('NaMolality', value=0.100, vary=False)
    aa.add_par('PbMolality', value=-6.000, vary=False)
    aa.add_par('benzeneMolality', value=0.200, vary=False)
    aa.add_par('tMitigation', value=87.914, vary=False)
    aa.add_par('CEC', value=32.073, vary=False)
    aa.add_par('Asbrine', value=-5.397, vary=False)
    aa.add_par('Babrine', value=-3.397, vary=False)
    aa.add_par('Cdbrine', value=-8.574, vary=False)
    aa.add_par('Pbbrine', value=-7.719, vary=False)
    aa.add_par('Benzene_brine', value=-8.610, vary=False)
    aa.add_par('Benzene_kd', value=-3.571, vary=False)
    aa.add_par('Benzene_decay', value=-2.732, vary=False)
    aa.add_par('PAH_brine', value=-7.118, vary=False)
    aa.add_par('PAH_kd', value=-0.985, vary=False)
    aa.add_par('PAH_decay', value=-3.371, vary=False)
    aa.add_par('phenol_brine', value=-6.666, vary=False)
    aa.add_par('phenol_kd', value=-1.342, vary=False)
    aa.add_par('phenol_decay', value=-3.546, vary=False)
    aa.add_par('porositySand', value=0.468, vary=False)
    aa.add_par('densitySand', value=2165.953, vary=False)
    aa.add_par('VG_mSand', value=0.627, vary=False)
    aa.add_par('VG_alphaSand', value=-4.341, vary=False)
    aa.add_par('permeabilitySand', value=-12.430, vary=False)
    aa.add_par('Clbrine', value=-0.339, vary=False)
    aa.add_par('calcitevol', value=0.165, vary=False)
    aa.add_par('V_goethite', value=0.004, vary=False)
    aa.add_par('V_illite', value=0.006, vary=False)
    aa.add_par('V_kaolinite', value=0.004, vary=False)
    aa.add_par('V_smectite', value=0.010, vary=False)

    aa.add_kwarg_linked_to_obs('co2_rate', ms.linkobs['CO2_aquifer1'])
    aa.add_kwarg_linked_to_obs('co2_mass', adapt.linkobs['mass_CO2_aquifer1'])
    aa.add_kwarg_linked_to_obs('brine_rate', ms.linkobs['brine_aquifer1'])
    aa.add_kwarg_linked_to_obs('brine_mass', adapt.linkobs['mass_brine_aquifer1'])

    # Add observations (output) from the alluvium aquifer model
    aa.add_obs('pH_volume')
    aa.add_obs('TDS_volume')

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically
    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    print('pH:',
          sm.collect_observations_as_time_series(aa, 'pH_volume'), sep='\n')
    print('------------------------------------------------------------------')
    print('TDS:',
          sm.collect_observations_as_time_series(aa, 'TDS_volume'), sep='\n')

    print('------------------------------------------------------------------')
    print('                          UQ illustration ')
    print('------------------------------------------------------------------')

    import random
    num_samples = 50
    ncpus = 1
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples,seed=random.randint(500, 1100))

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)

    # sm.obs_base_names:
    # Observations: 'sres.pressure', 'sres.CO2saturation', 'sres.mass_CO2_reservoir',
    # 'ms.brine_aquifer1', 'ms.CO2_aquifer1', 'adapt.mass_CO2_aquifer1',
    # 'adapt.mass_brine_aquifer1', 'ca.pH_volume', 'ca.TDS_volume'

    # Plot results
    num_obs = len(time_array)
    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(221)
    for j in range(num_samples):
        plt.semilogy(time_array[0:]/365.25, results[j, (5*num_obs):(6*num_obs)],
                     color='#1B4F72', linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel(r'Mass of CO$_2$ leaked (kg)', fontsize=14)
    plt.title(r'Mass of CO$_2$ leaked into aquifer 1', fontsize=14)
    plt.tick_params(labelsize=12)

    ax = fig.add_subplot(223)
    for j in range(num_samples):
        plt.semilogy(time_array[0:]/365.25, results[j, (6*num_obs):(7*num_obs)],
                     color='#000066', linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel(r'Mass of brine leaked (kg)', fontsize=14)
    plt.title(r'Mass of brine leaked into aquifer 1', fontsize=14)
    plt.tick_params(labelsize=12)

    ax = fig.add_subplot(222)
    for j in range(num_samples):
        plt.semilogy(time_array[0:]/365.25, results[j, (7*num_obs):(8*num_obs)],
                     color='#663300', linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel('pH Volume (m$^3$)', fontsize=14)
    plt.title(r'Volume of aquifer 1 below pH threshold', fontsize=14)
    plt.tight_layout()
    plt.tick_params(labelsize=12)

    ax = fig.add_subplot(224)
    for j in range(num_samples):
        plt.plot(time_array[0:]/365.25, results[j, (8*num_obs):(9*num_obs)],
                 color='#CC9900', linewidth=1)
    plt.xlabel('Time, t (years)', fontsize=14)
    plt.ylabel('TDS Volume (m$^3$)', fontsize=14)
    plt.title(r'Volume of aquifer 1 above TDS threshold', fontsize=14)
    plt.tight_layout()
    plt.tick_params(labelsize=12)

    to_save = False
    if to_save:
        plt.savefig('aluv_aquifer_example01_plot.png', dpi=300)
