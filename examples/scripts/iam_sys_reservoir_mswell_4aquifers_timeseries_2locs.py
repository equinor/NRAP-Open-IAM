'''
This example demonstrates the use of the time_series_plot() function from the
script time_series.py. A variety of figures are created with time_series_plot(),
with this variety intended to demonstrate different aspects of the function. For
example, some graphs have only one subplot while others have multiple subplots
(with the same or different types of metrics). Additionally, these examples
show how to assign a subplot title corresponding with a specific output name
(e.g., 'ca4_001.TDS_volume': 'Aquifer 4, Location 2' within the subplot variable
 when plotting TDS_volume).

This example is similar to iam_sys_reservoir_mswell_4aquifers_timeseries.py, but
it utilizes multiple locations. The time_series_plot() function can handle
multiple locations by location details in subplot titles and legends. Note that
the component names must include the location designation '_###', where ### is
the location number (e.g., 'sres_000' for the reservoir at location 1 and
'ca4_001' for Aquifer 4 at location 2). The function uses this approach for
location designations because it is used by the control file interface.

Created: 8-10-2022
Last modified: February 9th, 2023

Author: Nate Mitchell
LRST (Battelle|Leidos) supporting NETL
Nathaniel.Mitchell@netl.doe.gov

Contributor: Veronika Vasylkivska
Veronika.Vasylkivska@netl.doe.gov
LRST (Battelle|Leidos) supporting NETL

Example of run:
$ python iam_sys_reservoir_mswell_4aquifers_timeseries_2locs.py
'''

import sys
import os
import datetime
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (MultisegmentedWellbore, RateToMassAdapter, SystemModel,
                     SimpleReservoir, CarbonateAquifer, Stratigraphy)
import openiam.visualize as vis

if __name__ == "__main__":

    # SET UP DIRECTORY
    output_location = os.path.join(os.getcwd(), '..', '..', 'Output')

    if os.path.isdir(output_location) is False:
        os.mkdir(output_location)

    # The output folder name uses the current date and time
    example_output_name = 'timeseries_plot_2locs_example_' + str(datetime.date.today())
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
    sm = SystemModel(model_kwargs=sm_model_kwargs)

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

    x_vals = [100, 500]
    y_vals = [100, 500]

    sress = []
    mss = []
    adapts = []
    cas = []

    for loc, (x_val, y_val) in enumerate(zip(x_vals, y_vals)):
        # Add simple reservoir object
        sress.append(sm.add_component_model_object(
            SimpleReservoir(name='sres_{}'.format(loc+1), parent=sm,
                            locX=x_val, locY=y_val)))

        sress[-1].add_par('injRate', value=0.25, vary=False)
        sress[-1].add_par('logResPerm', value=-12, vary=False)

        sress[-1].add_par_linked_to_par(
            'numberOfShaleLayers', strata.deterministic_pars['numberOfShaleLayers'])
        sress[-1].add_par_linked_to_par(
            'shale1Thickness', strata.deterministic_pars['shale1Thickness'])
        sress[-1].add_par_linked_to_par(
            'shale2Thickness', strata.deterministic_pars['shale2Thickness'])
        sress[-1].add_par_linked_to_par(
            'shale3Thickness', strata.deterministic_pars['shale3Thickness'])
        sress[-1].add_par_linked_to_par(
            'shale4Thickness', strata.deterministic_pars['shale4Thickness'])
        sress[-1].add_par_linked_to_par(
            'shale5Thickness', strata.deterministic_pars['shale5Thickness'])
        sress[-1].add_par_linked_to_par(
            'aquifer1Thickness', strata.deterministic_pars['aquifer1Thickness'])
        sress[-1].add_par_linked_to_par(
            'aquifer2Thickness', strata.deterministic_pars['aquifer2Thickness'])
        sress[-1].add_par_linked_to_par(
            'aquifer3Thickness', strata.deterministic_pars['aquifer3Thickness'])
        sress[-1].add_par_linked_to_par(
            'aquifer4Thickness', strata.deterministic_pars['aquifer4Thickness'])
        sress[-1].add_par_linked_to_par(
            'reservoirThickness', strata.deterministic_pars['reservoirThickness'])
        sress[-1].add_par_linked_to_par(
            'datumPressure', strata.default_pars['datumPressure'])

        # Add observations of reservoir component model
        sress[-1].add_obs_to_be_linked('pressure')
        sress[-1].add_obs_to_be_linked('CO2saturation')
        sress[-1].add_obs('pressure')
        sress[-1].add_obs('CO2saturation')
        sress[-1].add_obs('mass_CO2_reservoir')

        # Add multisegmented wellbore component.
        mss.append(sm.add_component_model_object(
            MultisegmentedWellbore(name='ms_{}'.format(loc+1), parent=sm)))

        mss[-1].add_par('logWellPerm', value=-12, vary=False)
        mss[-1].add_par('logAquPerm', value=-10, vary=False)
        mss[-1].add_par('brineDensity', value=1100, vary=False)
        mss[-1].add_par('wellRadius', value=0.1, vary=False)

        # Add linked parameters: common to reservoir and wellbore components
        mss[-1].add_par_linked_to_par(
            'numberOfShaleLayers', strata.deterministic_pars['numberOfShaleLayers'])
        mss[-1].add_par_linked_to_par(
            'shale1Thickness', strata.deterministic_pars['shale1Thickness'])
        mss[-1].add_par_linked_to_par(
            'shale2Thickness', strata.deterministic_pars['shale2Thickness'])
        mss[-1].add_par_linked_to_par(
            'shale3Thickness', strata.deterministic_pars['shale3Thickness'])
        mss[-1].add_par_linked_to_par(
            'shale4Thickness', strata.deterministic_pars['shale4Thickness'])
        mss[-1].add_par_linked_to_par(
            'shale5Thickness', strata.deterministic_pars['shale5Thickness'])
        mss[-1].add_par_linked_to_par(
            'aquifer1Thickness', strata.deterministic_pars['aquifer1Thickness'])
        mss[-1].add_par_linked_to_par(
            'aquifer2Thickness', strata.deterministic_pars['aquifer2Thickness'])
        mss[-1].add_par_linked_to_par(
            'aquifer3Thickness', strata.deterministic_pars['aquifer3Thickness'])
        mss[-1].add_par_linked_to_par(
            'aquifer4Thickness', strata.deterministic_pars['aquifer4Thickness'])
        mss[-1].add_par_linked_to_par(
            'reservoirThickness', strata.deterministic_pars['reservoirThickness'])
        mss[-1].add_par_linked_to_par(
            'datumPressure', strata.default_pars['datumPressure'])

        # Add keyword arguments linked to the output provided by reservoir model
        mss[-1].add_kwarg_linked_to_obs('pressure', sress[-1].linkobs['pressure'])
        mss[-1].add_kwarg_linked_to_obs('CO2saturation', sress[-1].linkobs['CO2saturation'])

        # Add observations of multisegmented wellbore component model
        mss[-1].add_obs_to_be_linked('CO2_aquifer1')
        mss[-1].add_obs_to_be_linked('CO2_aquifer2')
        mss[-1].add_obs_to_be_linked('CO2_aquifer3')
        mss[-1].add_obs_to_be_linked('CO2_aquifer4')
        mss[-1].add_obs_to_be_linked('brine_aquifer1')
        mss[-1].add_obs_to_be_linked('brine_aquifer2')
        mss[-1].add_obs_to_be_linked('brine_aquifer3')
        mss[-1].add_obs_to_be_linked('brine_aquifer4')
        mss[-1].add_obs_to_be_linked('mass_CO2_aquifer1')
        mss[-1].add_obs_to_be_linked('mass_CO2_aquifer2')
        mss[-1].add_obs_to_be_linked('mass_CO2_aquifer3')
        mss[-1].add_obs_to_be_linked('mass_CO2_aquifer4')
        mss[-1].add_obs_to_be_linked('brine_atm')
        mss[-1].add_obs_to_be_linked('CO2_atm')
        mss[-1].add_obs('CO2_aquifer1')
        mss[-1].add_obs('CO2_aquifer2')
        mss[-1].add_obs('CO2_aquifer3')
        mss[-1].add_obs('CO2_aquifer4')
        mss[-1].add_obs('CO2_atm')
        mss[-1].add_obs('brine_aquifer1')
        mss[-1].add_obs('brine_aquifer2')
        mss[-1].add_obs('brine_aquifer3')
        mss[-1].add_obs('brine_aquifer4')
        mss[-1].add_obs('brine_atm')

        # Add adapter that transforms leakage rates to accumulated mass
        adapts.append(sm.add_component_model_object(
            RateToMassAdapter(name='adapt_{}'.format(loc+1), parent=sm)))
        adapts[-1].add_kwarg_linked_to_collection('CO2_aquifer1',
                                             [mss[-1].linkobs['CO2_aquifer1'],
                                              mss[-1].linkobs['CO2_aquifer2']])
        adapts[-1].add_kwarg_linked_to_collection('CO2_aquifer2',
                                             [mss[-1].linkobs['CO2_aquifer2'],
                                              mss[-1].linkobs['CO2_aquifer3']])
        adapts[-1].add_kwarg_linked_to_collection('CO2_aquifer3',
                                             [mss[-1].linkobs['CO2_aquifer3'],
                                              mss[-1].linkobs['CO2_aquifer4']])
        adapts[-1].add_kwarg_linked_to_collection('CO2_aquifer4',
                                             [mss[-1].linkobs['CO2_aquifer4'],
                                              mss[-1].linkobs['CO2_atm']])
        adapts[-1].add_kwarg_linked_to_collection('brine_aquifer1',
                                             [mss[-1].linkobs['brine_aquifer1'],
                                              mss[-1].linkobs['brine_aquifer2']])
        adapts[-1].add_kwarg_linked_to_collection('brine_aquifer2',
                                             [mss[-1].linkobs['brine_aquifer2'],
                                              mss[-1].linkobs['brine_aquifer3']])
        adapts[-1].add_kwarg_linked_to_collection('brine_aquifer3',
                                             [mss[-1].linkobs['brine_aquifer3'],
                                              mss[-1].linkobs['brine_aquifer4']])
        adapts[-1].add_kwarg_linked_to_collection('brine_aquifer4',
                                             [mss[-1].linkobs['brine_aquifer4'],
                                              mss[-1].linkobs['brine_atm']])
        adapts[-1].add_obs_to_be_linked('mass_CO2_aquifer1')
        adapts[-1].add_obs_to_be_linked('mass_CO2_aquifer2')
        adapts[-1].add_obs_to_be_linked('mass_CO2_aquifer3')
        adapts[-1].add_obs_to_be_linked('mass_CO2_aquifer4')
        adapts[-1].add_obs_to_be_linked('mass_brine_aquifer1')
        adapts[-1].add_obs_to_be_linked('mass_brine_aquifer2')
        adapts[-1].add_obs_to_be_linked('mass_brine_aquifer3')
        adapts[-1].add_obs_to_be_linked('mass_brine_aquifer4')
        adapts[-1].add_obs('mass_CO2_aquifer1')
        adapts[-1].add_obs('mass_CO2_aquifer2')
        adapts[-1].add_obs('mass_CO2_aquifer3')
        adapts[-1].add_obs('mass_CO2_aquifer4')
        adapts[-1].add_obs('mass_brine_aquifer1')
        adapts[-1].add_obs('mass_brine_aquifer2')
        adapts[-1].add_obs('mass_brine_aquifer3')
        adapts[-1].add_obs('mass_brine_aquifer4')

        # Add the carbonate aquifer model objects. They use the default parameters
        for ind in range(4):
            cas.append(sm.add_component_model_object(
                CarbonateAquifer(name='ca{}_{}'.format(ind + 1, loc+1), parent=sm)))

            cas[-1].add_par_linked_to_par(
                'aqu_thick', strata.deterministic_pars['aquifer{}Thickness'.format(ind + 1)])

            cas[-1].model_kwargs['x'] = [0.]
            cas[-1].model_kwargs['y'] = [0.]

            CO2_rate_obs_list = []
            brine_rate_obs_list = []
            CO2_mass_obs_list = []
            brine_mass_obs_list = []
            CO2_rate_obs_list.append(mss[-1].linkobs['CO2_aquifer' + str(ind + 1)])
            brine_rate_obs_list.append(mss[-1].linkobs['brine_aquifer' + str(ind + 1)])
            CO2_mass_obs_list.append(adapts[-1].linkobs['mass_CO2_aquifer' + str(ind + 1)])
            brine_mass_obs_list.append(adapts[-1].linkobs['mass_brine_aquifer' + str(ind + 1)])

            # Add aquifer component's keyword argument co2_rate linked
            # to the collection created above
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
    analysis_type = 'forward'
    plot_type = 'real'
    plot_data = {'TimeSeries': None,
                 'UseMarkers': False,       # default is False
                 'UseLines': True,          # default is True
                 'VaryLineStyles': False,   # default is False
                 'FigureDPI': 100,
                 'subplot': {'use': False}}

    # FIGURE 1: PRESSURE
    fig_num += 1
    output_names = ['pressure']
    output_list = {sress[0]: 'pressure', sress[1]: 'pressure'}

    plot_data['TimeSeries'] = output_names
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=None, plot_type=plot_type,
                         savefig='Pressure')

    # FIGURE 2: CO2 SATURATION AND CO2 MASS IN RESERVOIR
    fig_num += 1
    output_names = ['CO2saturation', 'mass_CO2_reservoir']
    output_list = {sress[0]: ['CO2saturation', 'mass_CO2_reservoir'],
                    sress[1]: ['CO2saturation', 'mass_CO2_reservoir']}
    subplot = {'use': True, 'ncols': 2}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=subplot,
                         plot_type=plot_type, savefig='CO2_Reservoir')

    # FIGURE 3: CO2 LEAKAGE RATES, AQUIFERS 1-2
    fig_num += 1
    output_names = ['CO2_aquifer1', 'CO2_aquifer2']
    output_list = {mss[0]: ['CO2_aquifer1', 'CO2_aquifer2'],
                    mss[1]: ['CO2_aquifer1', 'CO2_aquifer2']}
    subplot = {'use': True, 'ncols': 2}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=subplot,
                         plot_type=plot_type,
                         savefig='CO2_Leakage_Rates_Aqs1and2',
                         title='CO$_2$ Leakage to Aquifers 1-2')

    # FIGURE 4: CO2 LEAKAGE RATES, AQUIFERS 3-4
    fig_num += 1
    output_names = ['CO2_aquifer3', 'CO2_aquifer4']
    output_list = {mss[0]: ['CO2_aquifer3', 'CO2_aquifer4'],
                    mss[1]: ['CO2_aquifer3', 'CO2_aquifer4']}
    subplot = {'use': True, 'ncols': 2}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=subplot,
                         plot_type=plot_type,
                         savefig='CO2_Leakage_Rates_Aqs3and4',
                         title='CO$_2$ Leakage to Aquifers 3-4')

    # FIGURE 5: BRINE LEAKAGE RATES, LOCATION 1
    fig_num += 1
    output_names = ['brine_aquifer1', 'brine_aquifer2',
                    'brine_aquifer3', 'brine_aquifer4']
    output_list = {mss[0]: ['brine_aquifer1', 'brine_aquifer2',
                            'brine_aquifer3', 'brine_aquifer4']}
    subplot = {'use': False}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=subplot,
                         plot_type=plot_type, title='Location 1',
                         savefig='Brine_Leakage_Rates_Loc1')

    # FIGURE 6: BRINE LEAKAGE RATES, LOCATION 2
    fig_num += 1
    output_names = ['brine_aquifer1', 'brine_aquifer2',
                    'brine_aquifer3', 'brine_aquifer4']
    output_list = {mss[1]: ['brine_aquifer1', 'brine_aquifer2',
                            'brine_aquifer3', 'brine_aquifer4']}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=subplot,
                         plot_type=plot_type,
                         savefig='Brine_Leakage_Rates_Loc2', title='Location 2')

    # FIGURE 7: CO2 AND BRINE LEAKAGE TO ATMOSPHERE
    fig_num += 1
    output_names = ['CO2_atm']
    output_list = {mss[0]: 'CO2_atm', mss[1]: 'CO2_atm'}
    subplot = {'use': True, 'ncols': 1}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=subplot,
                         plot_type=plot_type, figsize=(6, 6),
                         savefig='CO2_Leakage_to_Atmosphere')

    # FIGURE 8: BRINE LEAKAGE TO ATMOSPHERE
    fig_num += 1
    output_names = ['brine_atm']
    output_list = {mss[0]: 'brine_atm', mss[1]: 'brine_atm'}
    subplot = {'use': True, 'ncols': 1}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=subplot,
                         plot_type=plot_type, figsize=(6, 6),
                         savefig='Brine_Leakage_to_Atmosphere')

    # FIGURE 9: CO2 MASSES, LOCATION 1
    fig_num += 1
    output_names = ['mass_CO2_aquifer1', 'mass_CO2_aquifer2',
                    'mass_CO2_aquifer3', 'mass_CO2_aquifer4']
    output_list = {adapts[0]: ['mass_CO2_aquifer1', 'mass_CO2_aquifer2',
                               'mass_CO2_aquifer3', 'mass_CO2_aquifer4']}
    subplot = {'use': False}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=subplot,
                         plot_type=plot_type,
                         savefig='CO2_Masses_Loc1', title='Location 1')


    # FIGURE 10: CO2 MASSES, LOCATION 2
    fig_num += 1
    output_names = ['mass_CO2_aquifer1', 'mass_CO2_aquifer2',
                    'mass_CO2_aquifer3', 'mass_CO2_aquifer4']
    output_list = {adapts[1]: ['mass_CO2_aquifer1', 'mass_CO2_aquifer2',
                               'mass_CO2_aquifer3', 'mass_CO2_aquifer4']}
    subplot = {'use': False}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=subplot,
                         plot_type=plot_type, title='Location 2',
                         savefig='CO2_Masses_Loc2')

    # FIGURE 11: BRINE MASSES, LOCATION 1
    fig_num += 1
    output_names = ['mass_brine_aquifer1', 'mass_brine_aquifer2',
                    'mass_brine_aquifer3', 'mass_brine_aquifer4']
    output_list = {adapts[0]: ['mass_brine_aquifer1', 'mass_brine_aquifer2',
                               'mass_brine_aquifer3', 'mass_brine_aquifer4']}
    subplot = {'use': False}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=subplot,
                         plot_type=plot_type,
                         savefig='Brine_Masses_Loc1', title='Location 1')

    # FIGURE 12: BRINE MASSES, LOCATION 2
    fig_num += 1
    output_names = ['mass_brine_aquifer1', 'mass_brine_aquifer2',
                    'mass_brine_aquifer3', 'mass_brine_aquifer4']
    output_list = {adapts[1]: ['mass_brine_aquifer1', 'mass_brine_aquifer2',
                               'mass_brine_aquifer3', 'mass_brine_aquifer4']}
    subplot = {'use': False}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=subplot,
                         plot_type=plot_type, title='Location 2',
                         savefig='Brine_Masses_Loc2')

    # FIGURE 11: TDS Plume, Location 1
    fig_num += 1
    output_names = ['TDS_volume']
    output_list = {cas[0]: 'TDS_volume', cas[1]: 'TDS_volume',
                    cas[2]: 'TDS_volume', cas[3]: 'TDS_volume'}

    subplot = {'use': True, 'ncols': 2,
               'ca1_000.TDS_volume': 'Aquifer 1, Location 1',
               'ca2_000.TDS_volume': 'Aquifer 2, Location 1',
               'ca3_000.TDS_volume': 'Aquifer 3, Location 1',
               'ca4_000.TDS_volume': 'Aquifer 4, Location 1'}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=subplot,
                         plot_type=plot_type,
                         savefig='TDS_Volumes_Location1')

    # FIGURE 11: TDS Plume, Location 2
    fig_num += 1
    output_names = ['TDS_volume']
    output_list = {cas[4]: 'TDS_volume', cas[5]: 'TDS_volume',
                    cas[6]: 'TDS_volume', cas[7]: 'TDS_volume'}

    subplot = {'use': True, 'ncols': 2,
               'ca1_001.TDS_volume': 'Aquifer 1, Location 2',
               'ca2_001.TDS_volume': 'Aquifer 2, Location 2',
               'ca3_001.TDS_volume': 'Aquifer 3, Location 2',
               'ca4_001.TDS_volume': 'Aquifer 4, Location 2'}

    plot_data['TimeSeries'] = output_names
    plot_data['subplot'] = subplot
    vis.time_series_plot(output_names, sm, None, plot_data, output_list,
                         name='Figure {}'.format(fig_num),
                         analysis=analysis_type, subplot=subplot,
                         plot_type=plot_type,
                         savefig='TDS_Volumes_Location2')
