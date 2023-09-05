'''
This example demonstrates the use of the time_series_plot() function from the
script time_series.py. A variety of figures are created with time_series_plot(),
with this variety intended to demonstrate different aspects of the function. For
example, some graphs have only one subplot while others have multiple subplots
(with the same or different types of metrics). Additionally, these examples
show how to assign a subplot title corresponding with a specific output name
(e.g., 'ca4.TDS_volume': 'Aquifer 4' within the subplot variable when plotting
 TDS_volume).

Last modified: February 9th, 2023

Author: Nate Mitchell
LRST (Battelle|Leidos) supporting NETL
Nathaniel.Mitchell@netl.doe.gov

Contributor: Veronika Vasylkivska
Veronika.Vasylkivska@netl.doe.gov
LRST (Battelle|Leidos) supporting NETL

Example of run:
$ python iam_sys_reservoir_mswell_4aquifers_timeseries.py
'''

import sys
import os
import datetime
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import MultisegmentedWellbore, RateToMassAdapter, SystemModel, \
    SimpleReservoir, CarbonateAquifer, Stratigraphy
import openiam.visualize as vis

if __name__ == "__main__":

    # SET UP DIRECTORY
    output_location = os.path.join(os.getcwd(), '..', '..', 'Output')

    if os.path.isdir(output_location) is False:
        os.mkdir(output_location)

    # The output folder name uses the current date and time
    example_output_name = 'timeseries_plot_example_' + str(datetime.date.today())
    example_output_location = os.path.join(output_location, example_output_name)

    if os.path.isdir(example_output_location) is False:
        os.mkdir(example_output_location)

    # SET UP MODEL
    # Number of years assessed in each model run
    num_years = 100
    # time is given in days
    time_array = 365.25 * np.arange(0.0, num_years + 1)
    sm_model_kwargs = {'time_array': time_array}

    # Create system model
    sm = SystemModel(model_kwargs = sm_model_kwargs)

    # Add stratigraphy component
    strata = sm.add_component_model_object(Stratigraphy(name='strata', parent=sm))

    strata.add_par('numberOfShaleLayers', value=5, vary=False)
    strata.add_par('reservoirThickness', value=100, vary=False)
    strata.add_par('shale1Thickness', value=250, vary=False)
    strata.add_par('aquifer1Thickness', value=150, vary=False)
    strata.add_par('shale2Thickness', value=400, vary=False)
    strata.add_par('aquifer2Thickness', value=200, vary=False)
    strata.add_par('shale3Thickness', value=150, vary=False)
    strata.add_par('aquifer3Thickness', value=200, vary=False)
    strata.add_par('shale4Thickness', value=200, vary=False)
    strata.add_par('aquifer4Thickness', value=100, vary=False)
    strata.add_par('shale5Thickness', value=100, vary=False)
    strata.add_par('datumPressure', value=101325, vary=False)

    # Add reservoir component
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm,
                                                         locX=100, locY=100))

    sres.add_par('injRate', value=0.25, vary=False)
    sres.add_par('logResPerm', value=-11.5, vary=False)
    sres.add_par_linked_to_par('numberOfShaleLayers',
                               strata.deterministic_pars['numberOfShaleLayers'])
    sres.add_par_linked_to_par('shale1Thickness',
                               strata.deterministic_pars['shale1Thickness'])
    sres.add_par_linked_to_par('shale2Thickness',
                               strata.deterministic_pars['shale2Thickness'])
    sres.add_par_linked_to_par('shale3Thickness',
                               strata.deterministic_pars['shale3Thickness'])
    sres.add_par_linked_to_par('shale4Thickness',
                               strata.deterministic_pars['shale4Thickness'])
    sres.add_par_linked_to_par('shale5Thickness',
                               strata.deterministic_pars['shale5Thickness'])
    sres.add_par_linked_to_par('aquifer1Thickness',
                               strata.deterministic_pars['aquifer1Thickness'])
    sres.add_par_linked_to_par('aquifer2Thickness',
                               strata.deterministic_pars['aquifer2Thickness'])
    sres.add_par_linked_to_par('aquifer3Thickness',
                               strata.deterministic_pars['aquifer3Thickness'])
    sres.add_par_linked_to_par('aquifer4Thickness',
                               strata.deterministic_pars['aquifer4Thickness'])
    sres.add_par_linked_to_par('reservoirThickness',
                               strata.deterministic_pars['reservoirThickness'])
    sres.add_par_linked_to_par('datumPressure',
                               strata.default_pars['datumPressure'])

    # Add observations of reservoir component model
    sres.add_obs_to_be_linked('pressure')
    sres.add_obs_to_be_linked('CO2saturation')
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')
    sres.add_obs('mass_CO2_reservoir')

    # Add multisegmented wellbore component.
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))

    ms.add_par('logWellPerm', value=-12, vary=False)
    ms.add_par('logAquPerm', value=-10, vary=False)
    ms.add_par('brineDensity', value=1100, vary=False)
    ms.add_par('wellRadius', value=0.1, vary=False)

    # Add linked parameters: common to reservoir and wellbore components
    ms.add_par_linked_to_par('numberOfShaleLayers',
                             strata.deterministic_pars['numberOfShaleLayers'])
    ms.add_par_linked_to_par('shale1Thickness',
                             strata.deterministic_pars['shale1Thickness'])
    ms.add_par_linked_to_par('shale2Thickness',
                             strata.deterministic_pars['shale2Thickness'])
    ms.add_par_linked_to_par('shale3Thickness',
                             strata.deterministic_pars['shale3Thickness'])
    ms.add_par_linked_to_par('shale4Thickness',
                             strata.deterministic_pars['shale4Thickness'])
    ms.add_par_linked_to_par('shale5Thickness',
                             strata.deterministic_pars['shale5Thickness'])
    ms.add_par_linked_to_par('aquifer1Thickness',
                             strata.deterministic_pars['aquifer1Thickness'])
    ms.add_par_linked_to_par('aquifer2Thickness',
                             strata.deterministic_pars['aquifer2Thickness'])
    ms.add_par_linked_to_par('aquifer3Thickness',
                             strata.deterministic_pars['aquifer3Thickness'])
    ms.add_par_linked_to_par('aquifer4Thickness',
                             strata.deterministic_pars['aquifer4Thickness'])
    ms.add_par_linked_to_par('reservoirThickness',
                             strata.deterministic_pars['reservoirThickness'])
    ms.add_par_linked_to_par('datumPressure',
                             strata.default_pars['datumPressure'])

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

    # Add observations of multisegmented wellbore component model
    ms.add_obs_to_be_linked('CO2_aquifer1')
    ms.add_obs_to_be_linked('CO2_aquifer2')
    ms.add_obs_to_be_linked('CO2_aquifer3')
    ms.add_obs_to_be_linked('CO2_aquifer4')
    ms.add_obs_to_be_linked('brine_aquifer1')
    ms.add_obs_to_be_linked('brine_aquifer2')
    ms.add_obs_to_be_linked('brine_aquifer3')
    ms.add_obs_to_be_linked('brine_aquifer4')
    ms.add_obs_to_be_linked('mass_CO2_aquifer1')
    ms.add_obs_to_be_linked('mass_CO2_aquifer2')
    ms.add_obs_to_be_linked('mass_CO2_aquifer3')
    ms.add_obs_to_be_linked('mass_CO2_aquifer4')
    ms.add_obs_to_be_linked('brine_atm')
    ms.add_obs_to_be_linked('CO2_atm')
    ms.add_obs('CO2_aquifer1')
    ms.add_obs('CO2_aquifer2')
    ms.add_obs('CO2_aquifer3')
    ms.add_obs('CO2_aquifer4')
    ms.add_obs('CO2_atm')
    ms.add_obs('brine_aquifer1')
    ms.add_obs('brine_aquifer2')
    ms.add_obs('brine_aquifer3')
    ms.add_obs('brine_aquifer4')
    ms.add_obs('brine_atm')

    # Add adapter that transforms leakage rates to accumulated mass
    adapt = sm.add_component_model_object(RateToMassAdapter(name='adapt',
                                                            parent=sm))
    adapt.add_kwarg_linked_to_collection('CO2_aquifer1',
                                         [ms.linkobs['CO2_aquifer1'],
                                          ms.linkobs['CO2_aquifer2']])
    adapt.add_kwarg_linked_to_collection('CO2_aquifer2',
                                         [ms.linkobs['CO2_aquifer2'],
                                          ms.linkobs['CO2_aquifer3']])
    adapt.add_kwarg_linked_to_collection('CO2_aquifer3',
                                         [ms.linkobs['CO2_aquifer3'],
                                          ms.linkobs['CO2_aquifer4']])
    adapt.add_kwarg_linked_to_collection('CO2_aquifer4',
                                         [ms.linkobs['CO2_aquifer4'],
                                          ms.linkobs['CO2_atm']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer1',
                                         [ms.linkobs['brine_aquifer1'],
                                          ms.linkobs['brine_aquifer2']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer2',
                                         [ms.linkobs['brine_aquifer2'],
                                          ms.linkobs['brine_aquifer3']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer3',
                                         [ms.linkobs['brine_aquifer3'],
                                          ms.linkobs['brine_aquifer4']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer4',
                                         [ms.linkobs['brine_aquifer4'],
                                          ms.linkobs['brine_atm']])
    adapt.add_obs_to_be_linked('mass_CO2_aquifer1')
    adapt.add_obs_to_be_linked('mass_CO2_aquifer2')
    adapt.add_obs_to_be_linked('mass_CO2_aquifer3')
    adapt.add_obs_to_be_linked('mass_CO2_aquifer4')
    adapt.add_obs_to_be_linked('mass_brine_aquifer1')
    adapt.add_obs_to_be_linked('mass_brine_aquifer2')
    adapt.add_obs_to_be_linked('mass_brine_aquifer3')
    adapt.add_obs_to_be_linked('mass_brine_aquifer4')
    adapt.add_obs('mass_CO2_aquifer1')
    adapt.add_obs('mass_CO2_aquifer2')
    adapt.add_obs('mass_CO2_aquifer3')
    adapt.add_obs('mass_CO2_aquifer4')
    adapt.add_obs('mass_brine_aquifer1')
    adapt.add_obs('mass_brine_aquifer2')
    adapt.add_obs('mass_brine_aquifer3')
    adapt.add_obs('mass_brine_aquifer4')

    # Add 4 carbonate aquifer model objects. They use the default parameters
    cas = []
    for ind in range(4):
        cas.append(sm.add_component_model_object(
            CarbonateAquifer(name='ca' + str(ind + 1), parent=sm)))

        cas[-1].add_par_linked_to_par(
            'aqu_thick', strata.deterministic_pars['aquifer{}Thickness'.format(ind + 1)])
        cas[-1].model_kwargs['x'] = [0.]
        cas[-1].model_kwargs['y'] = [0.]

        CO2_rate_obs_list = []
        brine_rate_obs_list = []
        CO2_mass_obs_list = []
        brine_mass_obs_list = []
        CO2_rate_obs_list.append(ms.linkobs['CO2_aquifer' + str(ind + 1)])
        brine_rate_obs_list.append(ms.linkobs['brine_aquifer' + str(ind + 1)])
        CO2_mass_obs_list.append(adapt.linkobs['mass_CO2_aquifer' + str(ind + 1)])
        brine_mass_obs_list.append(adapt.linkobs['mass_brine_aquifer' + str(ind + 1)])

        # Add aquifer component's keyword argument co2_rate linked to the collection created above
        cas[-1].add_kwarg_linked_to_collection('co2_rate', CO2_rate_obs_list)
        cas[-1].add_kwarg_linked_to_collection('brine_rate', brine_rate_obs_list)
        cas[-1].add_kwarg_linked_to_collection('co2_mass', CO2_mass_obs_list)
        cas[-1].add_kwarg_linked_to_collection('brine_mass', brine_mass_obs_list)

        # Add observations (output) from the carbonate aquifer model
        cas[-1].add_obs('pH_volume')
        cas[-1].add_obs('TDS_volume')

    # RUN SYSTEM MODEL
    sm.forward()  # system model is run deterministically

    # FIGURES
    os.chdir(example_output_location)
    fig_num = 0

    plot_data = {'TimeSeries': None,
                 'UseMarkers': False,       # default is False
                 'UseLines': True,          # default is True
                 'VaryLineStyles': False,   # default is False
                 'FigureDPI': 100,
                 'subplot': {'use': False}}

    # FIGURE 1: PRESSURE
    fig_num += 1
    output_names = ['pressure']
    output_list = {sres: 'pressure'}

    plot_data['TimeSeries'] = output_names
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis='forward', subplot=None, plot_type='real',
                         savefig='Pressure')

    # FIGURE 2: CO2 SATURATION AND CO2 MASS IN RESERVOIR
    fig_num += 1
    output_names = ['CO2saturation', 'mass_CO2_reservoir']
    output_list = {sres: ['CO2saturation', 'mass_CO2_reservoir']}
    subplot = {'use': True, 'ncols': 1}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis='forward', subplot=subplot,
                         plot_type='real', figsize=(6, 6),
                         savefig='CO2_Reservoir')

    # FIGURE 3: CO2 LEAKAGE RATES
    fig_num += 1
    output_names = ['CO2_aquifer1', 'CO2_aquifer2', 'CO2_aquifer3', 'CO2_aquifer4']
    output_list = {ms: ['CO2_aquifer1', 'CO2_aquifer2', 'CO2_aquifer3', 'CO2_aquifer4']}
    subplot = {'use': True, 'ncols': 2}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis='forward', subplot=subplot, plot_type='real',
                         savefig='CO2_Leakage_Rates',
                         title='CO$_2$ Leakage to Aquifers 1-4')

    # FIGURE 4: BRINE LEAKAGE RATES
    fig_num += 1
    output_names = ['brine_aquifer1', 'brine_aquifer2',
                    'brine_aquifer3', 'brine_aquifer4']
    output_list = {ms: ['brine_aquifer1', 'brine_aquifer2',
                        'brine_aquifer3', 'brine_aquifer4']}
    subplot = {'use': False}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis='forward', subplot=subplot, plot_type='real',
                         savefig='Brine_Leakage_Rates')

    # FIGURE 5: CO2 AND BRINE LEAKAGE TO ATMOSPHERE
    fig_num += 1
    output_names = ['CO2_atm', 'brine_atm']
    output_list = {ms: ['CO2_atm', 'brine_atm']}
    subplot = {'use': True, 'ncols': 1}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis='forward', subplot=subplot,
                         plot_type='real', figsize=(6, 6),
                         savefig='Leakage_to_Atmosphere',
                         title='Leakage to Atmosphere')

    # FIGURE 6: CO2 LEAKAGE RATES AND CO2 MASSES
    fig_num += 1
    output_names = ['CO2_aquifer1', 'CO2_aquifer2',
                    'mass_CO2_aquifer1', 'mass_CO2_aquifer2']
    output_list = {ms: ['CO2_aquifer1', 'CO2_aquifer2'],
                   adapt: {'mass_CO2_aquifer1', 'mass_CO2_aquifer2'}}
    subplot = {'use': True, 'ncols': 2}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis='forward', subplot=subplot, plot_type='real',
                         savefig='CO2_Leakage_and_Mass_Aqs1and2')

    # FIGURE 7: CO2 LEAKAGE RATES AND CO2 MASSES
    fig_num += 1
    output_names = ['CO2_aquifer3', 'CO2_aquifer4',
                    'mass_CO2_aquifer3', 'mass_CO2_aquifer4']
    output_list = {ms: ['CO2_aquifer3', 'CO2_aquifer4'],
                   adapt: {'mass_CO2_aquifer3', 'mass_CO2_aquifer4'}}
    subplot = {'use': True, 'ncols': 2}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis='forward', subplot=subplot, plot_type='real',
                         savefig='CO2_Leakage_and_Mass_Aqs3and4')

    # FIGURE 8: TDS PLUME VOLUMES
    fig_num += 1
    output_names = ['TDS_volume']
    output_list = {cas[0]: ['TDS_volume'], cas[1]: {'TDS_volume'},
                   cas[2]: ['TDS_volume'], cas[3]: {'TDS_volume'}}
    subplot = {'use': True, 'ncols': 2,
               'ca1.TDS_volume': 'Aquifer 1', 'ca2.TDS_volume': 'Aquifer 2',
               'ca3.TDS_volume': 'Aquifer 3', 'ca4.TDS_volume': 'Aquifer 4'}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis='forward', subplot=subplot, plot_type='real',
                         savefig='TDS_Plume_Volumes')
