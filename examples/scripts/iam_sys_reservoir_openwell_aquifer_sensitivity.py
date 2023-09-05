'''
This example couples the simple reservoir, open wellbore and
carbonate aquifer models. The saturation/pressure output produced by simple
reservoir model is used to drive leakage from a single open wellbore
model, which is passed to the input of an adapter that provides well
coordinates, |CO2| and brine leakage rates and cumulative mass fluxes to the
carbonate aquifer model. The sensitivity of outputs to the model parameters
is investigated and relevant plots are created.

Usage examples:
$ python iam_sys_reservoir_openwell_aquifer_sensitivity.py
'''
import sys
import os
import random
import datetime

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, SimpleReservoir, OpenWellbore,
                     CarbonateAquifer, RateToMassAdapter)
from openiam.visualize import (correlations_at_time, time_series_sensitivities,
                               multi_sensitivities_barplot,
                               simple_sensitivities_barplot)


if __name__ == "__main__":
    # For multiprocessing in Spyder
    __spec__ = None
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
    sres.add_par('shale1Thickness', min=900.0, max=1100., value=1000.0)
    sres.add_par('shale2Thickness', min=900.0, max=1100., value=1000.0)
    # Shale 3 has a fixed thickness of 11.2 m
    sres.add_par('shale3Thickness', value=11.2, vary=False)
    # Aquifer 1 (thief zone has a fixed thickness of 22.4)
    sres.add_par('aquifer1Thickness', value=22.4, vary=False)
    # Aquifer 2 (shallow aquifer) has a fixed thickness of 19.2
    sres.add_par('aquifer2Thickness', value=400, vary=False)
    # Reservoir has a fixed thickness of 51.2
    sres.add_par('reservoirThickness', value=51.2, vary=False)

    # Add observations of reservoir component model
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')
    sres.add_obs_to_be_linked('pressure')
    sres.add_obs_to_be_linked('CO2saturation')

    # Add open wellbore component
    ow = sm.add_component_model_object(OpenWellbore(name='ow', parent=sm))

    # Add parameters of open wellbore component
    ow.add_par('wellRadius', min=0.01, max=0.02, value=0.015)
    ow.add_par('logReservoirTransmissivity', min=-11.0, max=-9.0, value=-10.0)
    ow.add_par('logAquiferTransmissivity', min=-11.0, max=-9.0, value=-10.0)
    ow.add_par('brineSalinity', value=0.1, vary=False)

    # Add keyword arguments of open wellbore component
    ow.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
    ow.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

    # Add composite parameter of open wellbore component
    ow.add_composite_par('reservoirDepth',
                         expr='+'.join(['sres.shale1Thickness',
                                        'sres.shale2Thickness',
                                        'sres.shale3Thickness',
                                        'sres.aquifer1Thickness',
                                        'sres.aquifer2Thickness']))
    ow.add_composite_par(
        'wellTop', expr='sres.shale3Thickness + sres.aquifer2Thickness')

    # Add observations of open wellbore component
    ow.add_obs_to_be_linked('CO2_aquifer')
    ow.add_obs_to_be_linked('brine_aquifer')
    ow.add_obs_to_be_linked('brine_atm')
    ow.add_obs_to_be_linked('CO2_atm')
    ow.add_obs('CO2_aquifer')
    ow.add_obs('brine_aquifer')
    ow.add_obs('CO2_atm') # zero since well top is in aquifer
    ow.add_obs('brine_atm') # zero since well top is in aquifer

    # Add adapter that transforms leakage rates to accumulated mass
    adapt = sm.add_component_model_object(
        RateToMassAdapter(name='adapt', parent=sm))
    adapt.add_kwarg_linked_to_collection(
        'CO2_aquifer', [ow.linkobs['CO2_aquifer'], ow.linkobs['CO2_atm']])
    adapt.add_kwarg_linked_to_collection(
        'brine_aquifer', [ow.linkobs['brine_aquifer'], ow.linkobs['brine_atm']])
    adapt.add_obs_to_be_linked('mass_CO2_aquifer')
    adapt.add_obs_to_be_linked('mass_brine_aquifer')
    adapt.add_obs('mass_CO2_aquifer')
    adapt.add_obs('mass_brine_aquifer')

    # Add carbonate aquifer model object and define parameters
    ca = sm.add_component_model_object(CarbonateAquifer(name='ca', parent=sm))
    ca.add_par('perm_var', min=0.017, max=1.89, value=0.9535)
    ca.add_par('corr_len', min=1.0, max=3.95, value=2.475)
    ca.add_par('aniso', min=1.1, max=49.1, value=25.1)
    ca.add_par('mean_perm', min=-13.8, max=-10.3, value=-12.05)
    ca.add_par('aqu_thick', min=100., max=500., value=300.)
    ca.add_par('hyd_grad', min=2.88e-4, max=1.89e-2, value=9.59e-3)
    ca.add_par('calcite_ssa', min=0, max=1.e-2, value=5.5e-03)
    ca.add_par('organic_carbon', min=0, max=1.e-2, value=5.5e-03)
    ca.add_par('benzene_kd', min=1.49, max=1.73, value=1.61)
    ca.add_par('benzene_decay', min=0.15, max=2.84, value=1.5)
    ca.add_par('nap_kd', min=2.78, max=3.18, value=2.98)
    ca.add_par('nap_decay', min=-0.85, max=2.04, value=0.595)
    ca.add_par('phenol_kd', min=1.21, max=1.48, value=1.35)
    ca.add_par('phenol_decay', min=-1.22, max=2.06, value=0.42)
    ca.add_par('cl', min=0.1, max=6.025, value=0.776)
    ca.model_kwargs['x'] = [100.]
    ca.model_kwargs['y'] = [100.]

    CO2_rate_obs_list = []
    brine_rate_obs_list = []
    CO2_mass_obs_list = []
    brine_mass_obs_list = []
    CO2_rate_obs_list.append(ow.linkobs['CO2_aquifer'])
    brine_rate_obs_list.append(ow.linkobs['brine_aquifer'])
    CO2_mass_obs_list.append(adapt.linkobs['mass_CO2_aquifer'])
    brine_mass_obs_list.append(adapt.linkobs['mass_brine_aquifer'])

    # Add aquifer component's keyword argument co2_rate linked to the collection created above
    ca.add_kwarg_linked_to_collection('co2_rate', CO2_rate_obs_list)
    ca.add_kwarg_linked_to_collection('brine_rate', brine_rate_obs_list)
    ca.add_kwarg_linked_to_collection('co2_mass', CO2_mass_obs_list)
    ca.add_kwarg_linked_to_collection('brine_mass', brine_mass_obs_list)

    # Add observations (output) from the carbonate aquifer model
    ca.add_obs('pH_volume')
    ca.add_obs('TDS_volume')

    # Get sampling for sensitivity analysis
    saltelli_sample = sm.saltelli(50, calc_second_order=True)

    # Run model on sampling
    saltelli_sample.run(cpus=10, verbose=False)

    # Run sensitivity analysis on mass of CO2 leaked into aquifer
    mass_CO2_sensitivity = saltelli_sample.sobol(
        obsname='adapt.mass_CO2_aquifer_50', calc_second_order=True)

    output_directory = os.sep.join(
        ['..', '..', 'Output', 'Sensitivity_{date_time_stamp}'.format(
            date_time_stamp=datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))])

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    plt.figure(figsize=(15, 10))
    ax = plt.subplot(1, 2, 1)
    x_vals = list(range(len(sm.parnames)))
    plt.bar(x_vals, mass_CO2_sensitivity['S1'])
    ax.set_xticks(x_vals)
    ax.set_xticklabels(sm.parnames, rotation=90)
    plt.subplots_adjust(bottom=0.3)
    plt.title('First Order Sensitivity\nfor Mass of $CO_2$ Leaked into the Aquifer')
    plt.ylabel('First Order Sensitivity\n(Calculated using Sobol Method)')

    ax = plt.subplot(1, 2, 2)
    x_vals = list(range(len(sm.parnames)))
    plt.bar(x_vals, mass_CO2_sensitivity['ST'])
    ax.set_xticks(x_vals)
    ax.set_xticklabels(sm.parnames, rotation=90)
    plt.subplots_adjust(bottom=0.3, wspace=0.3)
    plt.title('Total Sensitivity\nfor Mass of $CO_2$ Leaked into the Aquifer')
    plt.ylabel('Total Sensitivity\n(Calculated using Sobol Method)')

    to_save = True
    if to_save:
        plt.savefig(os.sep.join([output_directory, 'sobol_co2_leak_plot.png']),
                    dpi=300)

#    # Calculate the Pearson correlation coefficient for mass of CO2 leaked into the aquifer
#    correlation_coeff, p_tail = calc_pearson_correlation_coeff(s, sm, 'adapt.mass_CO2_aquifer_50')
#    plt.figure(figsize=(10, 10))
#    ax = plt.subplot()
#    x_values = range(len(sm.parnames))
#    plt.bar(x_values, correlation_coeff)
#    ax.set_xticks(x_values)
#    ax.set_xticklabels(sm.parnames, rotation=90)
#    plt.subplots_adjust(bottom=0.3)
#    plt.title('Pearson Correlation Coefficient\nfor Mass of $CO_2$ Leaked into the Aquifer')
#    plt.ylabel('Pearson Correlation Coefficient')

#    # Now look at the aquifer plume size where pH has been reduced
#    ind_list = range(len(time_array))
#    ph_values = np.array([s.recarray['ca.pH_volume_'+str(indd)] for indd in ind_list])
#    ph = sm.collect_observations_as_time_series(ca, 'pH_volume')
#    plt.figure()
#    plt.plot(time_array/365.35, ph_values, 'b')
#    plt.xlabel('Time $(y)$')
#    plt.ylabel('Size of pH plume $(m^3)$')

    # Add capture point so sensitivity is calculated at time other than ending
    capture_point = 24   # in years
    cp_index = np.argmin(np.abs(time_array - capture_point*365.25))
    sensitivity_obs = 'ca.pH_volume'
    sen_obs_name = sensitivity_obs + '_{0}'.format(cp_index)
    # Run sensitivity analysis on aquifer pH plume size at capture point
    ph_sensitivity = saltelli_sample.sobol(obsname=sen_obs_name, calc_second_order=True)

    plt.figure(figsize=(15, 10))
    ax = plt.subplot(1, 2, 1)
    x_vals = list(range(len(sm.parnames)))
    plt.bar(x_vals, ph_sensitivity['S1'])
    ax.set_xticks(x_vals)
    ax.set_xticklabels(sm.parnames, rotation=90)
    plt.subplots_adjust(bottom=0.3)
    plt.title('First Order Sensitivity\nfor pH Plume in Aquifer')
    plt.ylabel('First Order Sensitivity\n(Calculated using Sobol Method)')

    ax = plt.subplot(1, 2, 2)
    x_vals = list(range(len(sm.parnames)))
    plt.bar(x_vals, ph_sensitivity['ST'])
    ax.set_xticks(x_vals)
    ax.set_xticklabels(sm.parnames, rotation=90)
    plt.subplots_adjust(bottom=0.3, wspace=0.3)
    plt.title('Total Sensitivity\nfor pH Plume in Aquifer')
    plt.ylabel('Total Sensitivity\n(Calculated using Sobol Method)')

    if to_save:
        plt.savefig(os.sep.join([output_directory, 'sobol_plot.png']), dpi=300)

    num_samples = 100
    ncpus = 10
    # Draw Latin hypercube samples of parameter values
    lhs_sample = sm.lhs(siz=num_samples, seed=random.randint(500, 1100))

    lhs_sample.run(cpus=ncpus, verbose=False)

    # Run RBD_fast sensitivity analysis on mass of CO2 leaked into aquifer
    mass_CO2_sensitivity = lhs_sample.rbd_fast(obsname='adapt.mass_CO2_aquifer_50')
    simple_sensitivities_barplot(
        mass_CO2_sensitivity, sm,
        title='RBD-Fast Sensitivity\nfor Mass of $CO_2$ Leaked into the Aquifer',
        ylabel='First Order Sensitivity\n(Calculated using RBD-Fast Method)',
        savefig=os.sep.join([output_directory, 'RBD-fast_co2_leak_plot.png']))

    # Add capture point so sensitivity is calculated at time other than ending
    capture_point = 24   # in years
    cp_index = np.argmin(np.abs(time_array - capture_point*365.25))
    sensitivity_obs = 'ca.pH_volume'
    sen_obs_name = sensitivity_obs + '_{0}'.format(cp_index)

    # Run sensitivity analysis on aquifer pH plume size at capture point
    ph_sensitivity = lhs_sample.rbd_fast(obsname=sen_obs_name)

    simple_sensitivities_barplot(
        ph_sensitivity, sm,
        title='First Order Sensitivity\nfor pH Plume in Aquifer',
        ylabel='First Order Sensitivity\n(Calculated using RBD-Fast Method)',
        savefig=os.sep.join([output_directory, 'RBD-fast_ph_Plume_plot.png']))

    time_series_sensitivities(
        'adapt.mass_CO2_aquifer', sm, lhs_sample, time_array,
        savefig=os.sep.join([output_directory, 'RBD-fast_MassCO2_time_series.png']))

    time_series_sensitivities(
        'ca.pH_volume', sm, lhs_sample, time_array,
        savefig=os.sep.join([output_directory,
                             'RBD-fast_pH_plume_vol_time_series.png']),
        capture_point=15)

    multi_sensitivities_barplot(
        ['sres.pressure_50', 'ow.CO2_aquifer_50', 'ow.brine_aquifer_50',
         'adapt.mass_CO2_aquifer_50', 'adapt.mass_brine_aquifer_50',
         'ca.pH_volume_50', 'ca.TDS_volume_50'], sm, lhs_sample,
        savefig=os.sep.join([output_directory, 'multi-bar-sensitivities.png']))

    corr_coeffs = correlations_at_time(
        lhs_sample, capture_point=50, excludes=['ow.CO2_atm', 'ow.brine_atm'],
        plot=True, figsize=(15, 15), printout=False, xrotation=90,
        title='Pearson Correlation Coefficients at 50 years',
        savefig=os.sep.join([output_directory, 'Corr_coeff_at_time_50.png']))
