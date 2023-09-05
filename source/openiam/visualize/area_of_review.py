# -*- coding: utf-8 -*-
"""
Code to create a map-view figure showing an area of review (AoR) based on
metrics like pH_volume and TDS_volume. This approach is based on the work of
Bacon et al. (2020), "Probabilistic risk-based Area of Review (AoR) determination
for a deep-saline carbon storage site."

Examples illustrating applications or setup of area_of_review_plot method:
    ControlFile_ex31a.yaml
    ControlFile_ex31b.yaml
    ControlFile_ex31c.yaml
    ControlFile_ex31d.yaml
    ControlFile_ex32a.yaml
    ControlFile_ex32b.yaml
    ControlFile_ex32c.yaml

Created: August 25th, 2022
Last Modified: June, 2023

@author: Nate Mitchell (Nathaniel.Mitchell@NETL.DOE.GOV)
@contributor: Veronika Vasylkivska (Veronika.Vasylkivska@NETL.DOE.GOV)
LRST (Battelle|Leidos) supporting NETL
"""
import csv
import os
import sys
import logging

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import cm
from matplotlib.lines import Line2D
from matk.sampleset import percentile

SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(SOURCE_DIR)

try:
    import openiam as iam
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

import openiam.cfi.commons as iamcommons
import openiam.cfi.strata as iam_strata

AOR_RESERVOIR_COMPONENTS = ['LookupTableReservoir', 'SimpleReservoir',
                        'AnalyticalReservoir', 'GenericReservoir']

BACKGROUND_COLOR = [0.67, 0.67, 0.67]

COLORMAP_OPTIONS = {'pressure': 'viridis',
                    'CO2saturation': 'plasma',
                    'pH_volume': 'YlGnBu',
                    'TDS_volume': 'YlOrRd',
                    'Dissolved_CO2_volume': 'YlGnBu',
                    'Dissolved_salt_volume': 'YlOrRd'}

DEFAULT_COLORMAP = 'cividis'

RC_FONT = {'family': 'Arial', 'weight': 'normal', 'size': None}

CBAR_LABELS = {'pressure': 'Pressure (MPa)',
               'CO2saturation': 'CO$_2$ Saturation [-]',
               'pH_volume': 'pH Plume Volume (m$^3$)',
               'TDS_volume': 'TDS Plume Volume (m$^3$)',
               'Dissolved_CO2_volume': 'CO$_2$ Plume Volume (m$^3$)',
               'Dissolved_salt_volume': 'Salt Plume Volume (m$^3$)'}

TITLE_OPTIONS_LHS = {
    'pressure': ''.join(['Maximum Reservoir Pressure{} Over Time ',
                         'Across {} LHS\nSimulations{}{}(Gray: 0 MPa){}']),
    'CO2saturation': ''.join(['Maximum Reservoir CO$_2$ Saturation{} Over Time ',
                              'Across {} LHS\nSimulations{}{}(Gray: 0){}']),
    'pH_volume': ''.join(['Maximum Aquifer {} pH Plume Volume Over Time ',
                          'Across {} LHS\nSimulations{}{}(Gray: 0 m$^3$){}']),
    'TDS_volume': ''.join(['Maximum Aquifer {} TDS Plume Volume Over Time ',
                           'Across {} LHS\nSimulations{}{}(Gray: 0 m$^3$){}']),
    'Dissolved_CO2_volume': ''.join(['Maximum Aquifer {} CO$_2$ Plume ',
                                     'Volume Over Time Across {} LHS\n',
                                     'Simulations{}{}(Gray: 0 m$^3$){}']),
    'Dissolved_salt_volume': ''.join(['Maximum Aquifer {} Salt Plume Volume ',
                                      'Over Time Across {} LHS\nSimulations',
                                      '{}{}(Gray: 0 m$^3$){}'])}

TITLE_OPTIONS_FORWARD = {
    'pressure': ''.join(['Maximum Reservoir Pressure{} Over Time{}\n',
                         '{}(Gray: 0 MPa){}']),
    'CO2saturation': ''.join(['Maximum Reservoir CO$_2$ Saturation{} Over Time{}\n',
                              '{}(Gray: 0){}']),
    'pH_volume': ''.join(['Maximum Aquifer {} pH Plume Volume Over Time{}\n',
                          '{}(Gray: 0 m$^3$){}']),
    'TDS_volume': ''.join(['Maximum Aquifer {} TDS Plume Volume Over Time{}\n',
                           '{}(Gray: 0 m$^3$){}']),
    'Dissolved_CO2_volume': ''.join(['Maximum Aquifer {} CO$_2$ Plume Volume Over ',
                                     'Time{}\n{}(Gray: 0 m$^3$){}']),
    'Dissolved_salt_volume': ''.join(['Maximum Aquifer {} Salt Plume Volume Over ',
                                      'Time{}\n{}(Gray: 0 m$^3$){}'])}

TITLE_OPTIONS_LHS_T_INDEX = {
    'pressure': ''.join(['Maximum Reservoir Pressure{} at t = {} years ',
                         'Across {} LHS\nSimulations{}{}(Gray: 0 MPa){}']),
    'CO2saturation': ''.join(['Maximum Reservoir CO$_2$ Saturation{} at t = {} years ',
                              'Across {} LHS\nSimulations{}{}(Gray: 0){}']),
    'pH_volume': ''.join(['Maximum Aquifer {} pH Plume Volume at t = {} years ',
                          'Across {} LHS\nSimulations{}{}(Gray: 0 m$^3$){}']),
    'TDS_volume': ''.join(['Maximum Aquifer {} TDS Plume Volume at t = {} years ',
                           'Across {} LHS\nSimulations{}{}(Gray: 0 m$^3$){}']),
    'Dissolved_CO2_volume': ''.join(['Maximum Aquifer {} CO$_2$ Plume Volume at ',
                                     't = {} years Across {} LHS\n',
                                     'Simulations{}{}(Gray: 0 m$^3$){}']),
    'Dissolved_salt_volume': ''.join(['Maximum Aquifer {} Salt Plume Volume at ',
                                      't = {} years Across {} LHS\n',
                                      'Simulations{}{}(Gray: 0 m$^3$){}'])}

TITLE_OPTIONS_FORWARD_T_INDEX = {
    'pressure': ''.join(['Maximum Reservoir Pressure{} at ',
                         't = {} years{}\n{}(Gray: 0 MPa){}']),
    'CO2saturation': ''.join(['Maximum Reservoir CO$_2$ Saturation{} at ',
                              't = {} years{}\n{}(Gray: 0){}']),
    'pH_volume': ''.join(['Maximum Aquifer {} pH Plume Volume at t = {} years{}\n',
                          '{}(Gray: 0 m$^3$){}']),
    'TDS_volume': ''.join(['Maximum Aquifer {} TDS Plume Volume at t = {} years{}\n',
                           '{}(Gray: 0 m$^3$){}']),
    'Dissolved_CO2_volume': ''.join(['Maximum Aquifer {} CO$_2$ Plume Volume ',
                                     'at t = {} years{}\n{}(Gray: 0 m$^3$){}']),
    'Dissolved_salt_volume': ''.join(['Maximum Aquifer {} Salt Plume Volume ',
                                      'at t = {} years{}\n{}(Gray: 0 m$^3$){}'])}

TITLE_RANGE = {
    'pressure': 'Range: {} MPa to {} MPa ',
    'CO2saturation': 'Range: {} to {} ',
    'pH_volume': 'Range: {} m$^3$ to {} m$^3$ ',
    'TDS_volume': 'Range: {} m$^3$ to {} m$^3$ ',
    'Dissolved_CO2_volume': 'Range: {} m$^3$ to {} m$^3$ ',
    'Dissolved_salt_volume': 'Range: {} m$^3$ to {} m$^3$ ',}

TITLE_SINGLE_VALUE = {
    'pressure': 'Value: {} MPa ',
    'CO2saturation': 'Value: {} ',
    'pH_volume': 'Value: {} m$^3$ ',
    'TDS_volume': 'Value: {} m$^3$ ',
    'Dissolved_CO2_volume': 'Value: {} m$^3$ ',
    'Dissolved_salt_volume': 'Value: {} m$^3$ ',}

CSV_FILE_COLUMNS = {'pressure': 'Max pressure [MPa]',
                    'CO2saturation': 'Max CO2 saturation [-]',
                    'pH_volume': 'Max pH plume volume [m^3]',
                    'TDS_volume': 'Max TDS plume volume [m^3]',
                    'Dissolved_CO2_volume': 'Max CO2 plume volume [m^3]',
                    'Dissolved_salt_volume': 'Max salt plume volume [m^3]'}

CSV_FILE_NAME_TAGS = {'pressure': 'Pressure', 'CO2saturation': 'CO2_Saturation',
                      'pH_volume': 'pH_Volume', 'TDS_volume': 'TDS_Volume',
                      'Dissolved_CO2_volume': 'Dissolved_CO2_Volume',
                      'Dissolved_salt_volume': 'Dissolved_Salt_Volume'}

MAX_VAL_ADJUST = 1.1

# Assumed gravitational acceleration for critical pressure calculation, m/(s^2)
GRAV_ACCEL = 9.8

# Assumed density of water for critical pressure calculation, kg/(m^3).
WATER_DENSITY = 1000

def area_of_review_plot(yaml_data, model_data, output_names, sm, s,
                        output_list, locations, name='AoR_Figure1', analysis='forward',
                        savefig=None, title=None, figsize=(10, 8), figdpi=100,
                        fontname='Arial', gen_font_size=12, axis_label_font_size=14,
                        title_font_size=14, colormap=None, bold_labels=True,
                        save_results=False, enforce_levels=True):
    """
    Makes a map-view figure showing the maximum values of a given metric (e.g.,
    pressure, CO2saturation, pH_volume, or TDS_volume) at each location used for
    for an OpenWellbore and aquifer component (e.g., FutureGen2Aquifer). Some of
    these plots (e.g., pH_volume and TDS_volume) are meant to help define an Area
    of Review (AoR).

    :param yaml_data: Dictionary of input values from the entire .yaml file
    :type yaml_data: dict

    :param model_data: Input from the 'ModelParams' section of the .yaml file
    :type model_data: dict

    :param output_names: List of observation names to match with component models and plot.
    :type output_names: list

    :param sm: OpenIAM System model for which plots are created
    :type sm: openiam.SystemModel object

    :param s: SampleSet with results of analysis performed on the system model.
        None for forward analysis.
    :type s: matk.SampleSet object

    :param output_list: Dictionary mapping component models to observations
    :type output_list: dict

    :param locations: dictionary of locations assigned to each applicable component
    :type locations: dict

    :param name: Figure Name to be used/created.
    :type name: str

    :param analysis: Type of OpenIAM system analysis performed ('forward',
        'lhs', or 'parstudy')
    :type analysis: str

    :param savefig: Filename to save resulting figure to. No file saved if None.
    :type savefig: str

    :param title: Optional Title for figure
    :type title: str

    :param figsize: width and height of the figure (width, height), in inches.
        Default value is (10, 8).
    :type figsize: tuple

    :param figdpi: dpi (dots-per-inch) for the figure
    :type figdpi: float or int

    :param fontname: name of the font type to be used
    :type fontname: str

    :param gen_font_size: fontsize for tick labels, etc.
    :type gen_font_size: float or int

    :param axis_label_font_size: fontsize for x and y axis labels
    :type axis_label_font_size: float or int

    :param title_font_size: fontsize for the title
    :type title_font_size: float or int

    :param colormap: string designation for a particular colormap (e.g., 'viridis')
    :type bold_labels: str

    :param bold_labels: option to use bold x and y labels and bold titles. Set to
        True for bold labels, False for normal labels.
    :type bold_labels: bool

    :param save_results: option to save the maximum values at each x and y value
        as a .csv file
    :type save_results: bool

    :param enforce_levels: option to enforce minimum and maximum levels for colorbars
        on plots
    :type enforce_levels: bool

    :return: None
    """
    # When a hypothetical wellbore is placed on the injection site itself, the
    # pressure results can include nan or Inf values. The code is set up to
    # ignore those results, but they can still cause numpy to print warning
    # statements. This statement supresses those warning statements.
    np.seterr(invalid='ignore')

    time = sm.time_array / 365.25

    yaml_input = get_AoR_yaml_input(yaml_data, name)

    if yaml_input['SaveCSVFiles'] is not None:
        save_results = yaml_input['SaveCSVFiles']

    if yaml_input['dpi_input'] is not None:
        figdpi = yaml_input['dpi_input']

    critPressureInput = yaml_input['CriticalPressureMPa']

    # Get the stratigraphy information from the .yaml file
    strata_var_info = iam_strata.get_strata_var_info_from_yaml(yaml_data)
    var_type = strata_var_info['var_type']

    # This option specifies whether to evaluate the max. values over all times
    # (False) or evaluate the max. values for specific times (True).
    time_option = False
    time_index_list = None
    if yaml_input['TimeList'] is not None:
        time_option = True
        time_list = yaml_input['TimeList']

        if time_list == 'All':
            time_index_list = range(len(time))
        else:
            time_index_list = get_t_indices(time_list, time)

    components = list(sm.component_models.values())
    # If injection sites need to be plotted, get the injection sites
    if yaml_input['plot_injection_sites']:
        InjectionCoordx = []
        InjectionCoordy = []

        for comp in components:
            if comp.class_type in AOR_RESERVOIR_COMPONENTS:
                if comp.class_type != 'LookupTableReservoir':
                    InjectionCoordx.append(comp.injX)
                    InjectionCoordy.append(comp.injY)

        InjectionCoordx = np.unique(InjectionCoordx).tolist()
        InjectionCoordy = np.unique(InjectionCoordy).tolist()

    else:
        InjectionCoordx = None
        InjectionCoordy = None

    aq_number = None

    # Get the OpenWellbore data required for process_wellbore_locations()
    for comp_model in model_data['Components']:
        comp_data = yaml_data[comp_model]

        ow_cmpnt_found = False
        if 'Type' in comp_data:
            if comp_data['Type'] == 'OpenWellbore':
                comp_model_ow = comp_model
                comp_data_ow = comp_data

                if 'LeakTo' in comp_data:
                    leakTo = comp_data['LeakTo']

                    if leakTo[0:7] == 'aquifer':
                        aq_number = int(leakTo[7:None])

                ow_cmpnt_found = True
                grid_option = 'grid' in comp_data_ow['Locations']

                break

    if not ow_cmpnt_found:
        err_msg = "".join(["'Type: OpenWellbore' was not found in the components ",
                           "specified in the .yaml file. This component type is ",
                           "required for the creation of an AoR plot."])
        logging.error(err_msg)
        raise KeyError(err_msg)

    # Get locations associated with given open wellbore component
    locations_ow = locations[comp_model_ow]

    x_loc = np.array(locations_ow['coordx'])
    y_loc = np.array(locations_ow['coordy'])

    # This variable is used to check if the current figure uses pressure and if
    # a critical pressure was given
    pressure_critP_check = False

    if 'pressure' in output_names and critPressureInput is not None:
        pressure_critP_check = True

        critP_warning_msg = "".join([
            'The calculated critical pressure for an OpenWellbore component ',
            'was found to be 0 Pa. This result can occur if the simulation uses ',
            'the Latin Hypercube Sampling (lhs) or Parameter Study (parstudy) analysis ',
            'types and the wellTop and/or reservoirDepth parameters are composite ',
            'parameters. These parameters will be composite parameters if they are ',
            'set only through the .yaml entry LeakTo: aquifer#, where # is a unit number. ',
            'To avoid this problem, explicitly set the wellTop and reservoirDepth ',
            'parameters. You can use string inputs, like wellTop: aquifer2Depth and ',
            'reservoirDepth: shale1Depth. The critical pressure option will not be ',
            'used for the AoR plot.'])

        # This is used to check if the error was already logged. That way, it
        # can be logged only once.
        critP_error_print_check = False

    elif not 'pressure' in output_names and critPressureInput is not None:
        # If the metric is not pressure but the .yaml plot entry included a
        # critical pressure entry, get rid of that entry.
        critPressureInput = None

    # Make the figures
    if not time_option:
        # One figure with max values from all model times
        results, critPressure = get_AoR_results(
            x_loc, output_names, sm, s, output_list, yaml_data, analysis=analysis,
            time_option=time_option, critPressureInput=critPressureInput,
            var_type=var_type)

        if pressure_critP_check:
            if np.max(critPressure) == 0 and critPressureInput == 'Calculated':
                if not critP_error_print_check:
                    logging.warning(critP_warning_msg)
                    critP_error_print_check = True

                critPressureInput = None

            elif critPressure is None and critPressureInput == 'Calculated':
                critPressureInput = None

            elif critPressureInput == 'Calculated':
                # Convert to MPa
                critPressureInput = critPressure / 1.0e+6

        plot_AoR_results(aq_number, x_loc, y_loc, results,
                         yaml_data, model_data, output_names, sm, name=name,
                         analysis=analysis, savefig=savefig, title=title,
                         figsize=figsize, figdpi=figdpi,
                         gen_font_size=gen_font_size,
                         axis_label_font_size=axis_label_font_size,
                         title_font_size=title_font_size, colormap=colormap,
                         bold_labels=bold_labels, save_results=save_results,
                         InjectionCoordx=InjectionCoordx,
                         InjectionCoordy=InjectionCoordy,
                         grid_option=grid_option, critPressureInput=critPressureInput,
                         var_type=var_type)
    else:
        # Get the min and max values over time, so the colorbar can use the
        # same limits in each figure.
        min_value = None
        max_value = None

        if enforce_levels:
            min_value = 9.9e+99
            max_value = -9.9e+99

            for time_index in time_index_list:
                results, critPressure = get_AoR_results(
                    x_loc, output_names, sm, s, output_list, yaml_data,
                    analysis=analysis, time_option=time_option,
                    time_index=time_index, var_type=var_type)

                if len(results[results > 0].tolist()) > 0:
                    if np.min(results[results > 0]) < min_value:
                        min_value = np.min(results[results > 0])

                if len(results[results > 0].tolist()) > 0:
                    if np.max(results[results > 0]) > max_value:
                        results_temp = np.max(results[results > 0])
                        max_value = np.min(results_temp[results_temp < np.Inf])

            if min_value == 9.9e+99:
                min_value = 0

            if max_value == -9.9e+99:
                max_value = 1

        # Make a figure for each time
        for time_index in time_index_list:
            results, critPressure = get_AoR_results(
                x_loc, output_names, sm, s, output_list, yaml_data, analysis=analysis,
                time_option=time_option, time_index=time_index,
                critPressureInput=critPressureInput, var_type=var_type)

            if pressure_critP_check:
                if critPressureInput == 'Calculated':
                    if np.max(critPressure) == 0:
                        if not critP_error_print_check:
                            logging.warning(critP_warning_msg)
                            critP_error_print_check = True

                        critPressureInput = None

                    elif critPressure is None:
                        critPressureInput = None

                    else:
                        # Convert to MPa
                        critPressureInput = critPressure / 1.0e+6

            plot_AoR_results(aq_number, x_loc, y_loc, results, yaml_data,
                             model_data, output_names, sm, name=name,
                             analysis=analysis, savefig=savefig, title=title,
                             figsize=figsize, figdpi=figdpi,
                             gen_font_size=gen_font_size,
                             axis_label_font_size=axis_label_font_size,
                             title_font_size=title_font_size, colormap=colormap,
                             bold_labels=bold_labels, save_results=save_results,
                             time_option=time_option, time_index=time_index,
                             InjectionCoordx=InjectionCoordx, InjectionCoordy=InjectionCoordy,
                             grid_option=grid_option, enforce_levels=enforce_levels,
                             min_value=min_value, max_value=max_value,
                             critPressureInput=critPressureInput, var_type=var_type)


def get_AoR_results(x_loc, output_names, sm, s, output_list, yaml_data,
                    analysis='forward',time_option=False, time_index=None,
                    critPressureInput=None, var_type='noVariation'):
    """
    Evaluates and returns the maximum values of a metric for all locations.
    These maximum values are then used in a plot that is meant to inform the
    boundaries of an area of review (AoR). If time_option is set to True, the
    results will only be evaluated at the time_index provided. Otherwise, the
    results returned are the maximum values across all times.
    """
    time = sm.time_array / 365.25

    # This is used to store the maximum value of a metric at each location
    results = np.zeros((len(x_loc), 1))

    if var_type == 'noVariation':
        critPressure = None
    elif var_type in ['strikeAndDip', 'LookupTable']:
        critPressure = np.zeros((len(x_loc), 1))

    for output_nm in output_names:
        # output_list is all the components with augmented names
        for output_component in list(output_list.keys()):

            # It it's an OpenWellbore and critPressureInput is 'Calculated',
            # get the critical pressure
            if isinstance(output_component, iam.OpenWellbore) and \
                    critPressureInput == 'Calculated':
                if var_type == 'noVariation' and critPressure is None:
                    # If using uniform stratigraphy, only do this once
                    critPressureVal = get_crit_pressure(
                        output_component, sm=sm, yaml_data=yaml_data)

                    critPressure = critPressureVal

                elif var_type in ['strikeAndDip', 'LookupTable']:
                    critPressureVal = get_crit_pressure(
                        output_component, sm=sm, yaml_data=yaml_data)

                    # Find the location index
                    loc_ref = int(output_component.name.split('_')[-1])

                    critPressure[loc_ref] = critPressureVal

            if output_nm in output_list[output_component]:
                full_obs_nm = '.'.join([output_component.name, output_nm])

                aq_comp_check = False
                res_comp_check = False

                if isinstance(output_component, (iam.FutureGen2Aquifer,
                                                 iam.FutureGen2AZMI,
                                                 iam.GenericAquifer,
                                                 iam.DeepAlluviumAquifer)):
                    aq_comp_check = True

                if isinstance(output_component, (iam.SimpleReservoir,
                                                 iam.AnalyticalReservoir,
                                                 iam.LookupTableReservoir)):
                    res_comp_check = True

                # Get maximum values at each location
                if analysis == 'forward':
                    values = sm.collect_observations_as_time_series(
                        output_component, output_nm)

                    if aq_comp_check or res_comp_check:
                        # Find the location index
                        loc_ref = int(output_component.name.split('_')[-1])

                        if not time_option:
                            if output_nm == 'pressure':
                                # Get the maximum pressure in MPa.
                                results[loc_ref] = max(values) / 1e+6
                            else:
                                results[loc_ref] = max(values)
                        else:
                            if output_nm == 'pressure':
                                # Get the maximum pressure in MPa
                                results[loc_ref] = values[time_index] / 1e+6
                            else:
                                results[loc_ref] = values[time_index]

                elif analysis in ['lhs', 'parstudy']:
                    ind_list = list(range(len(time)))

                    obs_names = [full_obs_nm + '_{0}'.format(indd)
                                 for indd in ind_list]
                    obs_percentiles = percentile(s.recarray[obs_names],
                                                 [0, 25, 50, 75, 100])

                    obs_t0 = [full_obs_nm + '_0']
                    obs_percentiles_t0 = percentile(s.recarray[obs_t0],
                                                    [0, 25, 50, 75, 100])

                    if aq_comp_check or res_comp_check:
                        # Find the location index
                        loc_ref = int(output_component.name.split('_')[-1])

                        if not time_option:
                            if output_nm == 'pressure':
                                # Get the maximum pressure in MPa
                                results[loc_ref] = max(obs_percentiles[4, :]) / 1e+6
                            else:
                                results[loc_ref] = max(obs_percentiles[4, :])
                        else:
                            if output_nm == 'pressure':
                                # Get the maximum pressure in MPa
                                results[loc_ref] = obs_percentiles[4, time_index] / 1e+6
                            else:
                                results[loc_ref] = obs_percentiles[4, time_index]

    return results, critPressure


def plot_AoR_results(aq_number, x_loc, y_loc, results, yaml_data, model_data,
                     output_names, sm, name='AoR_Figure1',
                     analysis='forward', savefig=None, title=None, figsize=(10, 8),
                     figdpi=100, gen_font_size=12,
                     axis_label_font_size=14, title_font_size=14, colormap=None,
                     bold_labels=True, save_results=False, time_option=False,
                     time_index=None, InjectionCoordx=None, InjectionCoordy=None,
                     grid_option=True, enforce_levels=False, min_value=None,
                     max_value=None, critPressureInput=None, var_type='noVariation'):
    """
    Plots maximum results across all x and y values (x_loc and y_loc) for either
    all times (time_option is False) or a specific time (time option is True).
    """
    time = sm.time_array / 365.25

    if bold_labels:
        selected_labelfontweight = 'bold'
    else:
        selected_labelfontweight = 'normal'

    yaml_input = get_AoR_yaml_input(yaml_data, name)

    output_nm = output_names[0]
    # If no colormap is specified, select a colormap. This selection
    # is done so that the default option of colormap = None produces
    # distinct colormaps for different metrics (e.g., pH_volume vs. TDS_volume).
    if colormap is None:
        colormap = COLORMAP_OPTIONS.get(output_nm, DEFAULT_COLORMAP)

    # Figures
    font = RC_FONT
    font['size'] = gen_font_size
    plt.rc('font', **font)

    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    # AoR Figure
    fig = plt.figure(figsize=figsize, dpi=figdpi)

    ax = plt.gca()
    ax.set_facecolor(BACKGROUND_COLOR)
    plt.plot(x_loc / 1000, y_loc / 1000, linestyle='none', marker='o',
             color='k', markeredgewidth=1, markersize=5, markerfacecolor='none',
             label='Hypothetical Open Wellbore for Area of Review')

    # This is the number of columns in the legend
    ncol_number = 1

    cmap = plt.cm.get_cmap(colormap)

    # This is used to specify if the critical pressure was exceeded, if the metric
    # is pressure and a critical pressure was given in the .yaml plot entry.
    title_pressure = ''

    # This specifies whether the critical pressure is included in the domain.
    # If it is included in the domain, a red line is shown for it.
    Pcrit_Included = False

    if np.max(results[:, 0]) != 0:
        # I use np.max(results[:, 0]) * MAX_VAL_ADJUST in levels because having
        # the maximum level equal to the maximum value can cause the area with
        # the highest values to be left out (i.e., there would be an uncolored
        # area there).
        if np.Inf not in results[:, 0]:
            min_val = np.min(results[:, 0][results[:, 0] > 0])
            a, b = '{:.2e}'.format(min_val).split('e')
            b = int(b)
            min_val_str = r'${}\times10^{{{}}}$'.format(a, b)

            # Having > 0 here prevents the inclusion of nan values
            max_val = np.max(results[:, 0][results[:, 0] > 0])
            a, b = '{:.2e}'.format(max_val).split('e')
            b = int(b)
            max_val_str = r'${}\times10^{{{}}}$'.format(a, b)

            if min_val == max_val:
                range_str = TITLE_SINGLE_VALUE.get(
                    output_nm, 'Value: {} ')
                range_str = range_str.format(min_val_str)

            else:
                range_str = TITLE_RANGE.get(
                    output_nm, 'Range: {} to {} ')
                range_str = range_str.format(min_val_str, max_val_str)

            if enforce_levels and time_option:
                min_val = min_value
                max_val = max_value

            interval = ((max_val * MAX_VAL_ADJUST) - min_val) / 100
            levels = np.arange(min_val, (max_val * MAX_VAL_ADJUST) + interval,
                               interval)

        elif np.Inf in results[:, 0]:
            # If an OpenWellbore is placed on the injection site, you can get
            # an infinite pressure value. I include "< np.Inf" to avoid that case.
            warning_msg = ''.join([
                'The results used for the AoR plot included an infinite value. ',
                'Infinite pressures can occur if an OpenWellbore is placed ',
                'at the injection location itself. Infinite values will be ',
                'excluded from the AoR plot.'])
            logging.warning(warning_msg)

            min_val = np.min(results[:, 0][results[:, 0] > 0])
            a, b = '{:.2e}'.format(min_val).split('e')
            b = int(b)
            min_val_str = r'${}\times10^{{{}}}$'.format(a, b)

            max_val = np.max(results[:, 0][results[:, 0] < np.Inf])
            a, b = '{:.2e}'.format(max_val).split('e')
            b = int(b)
            max_val_str = r'${}\times10^{{{}}}$'.format(a, b)

            if min_val == max_val:
                range_str = TITLE_SINGLE_VALUE.get(
                    output_nm, 'Value: {} ')
                range_str = range_str.format(min_val_str)

            else:
                range_str = TITLE_RANGE.get(
                    output_nm, 'Range: {} to {} ')
                range_str = range_str.format(min_val_str, max_val_str)

            if enforce_levels and time_option:
                min_val = min_value
                max_val = max_value

            interval = ((max_val * MAX_VAL_ADJUST) - min_val) / 100
            levels = np.arange(min_val, (max_val * MAX_VAL_ADJUST) + interval,
                               interval)

        try:
            # Gets rid of any nan values (from a wellbore being on the injection site)
            x_loc_temp = x_loc[results[:, 0] > 0]
            y_loc_temp = y_loc[results[:, 0] > 0]
            results_temp = results[results[:, 0] > 0]

            # Gets rid of any Inf values (from a wellbore being on the injection site)
            x_loc_temp = x_loc_temp[results_temp[:, 0] < np.Inf]
            y_loc_temp = y_loc_temp[results_temp[:, 0] < np.Inf]
            results_temp = results_temp[results_temp[:, 0] < np.Inf]

            plt.tricontourf(
                x_loc_temp / 1000.0, y_loc_temp / 1000.0,
                results_temp[:, 0], levels, cmap=colormap,
                locator=ticker.MaxNLocator(nbins=100, prune='lower'))

            cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                                values=levels, format=ticker.FuncFormatter(fmt))

            cbar_ticks = cbar.ax.get_yticks()
            cbar_ticks = cbar_ticks[cbar_ticks > 0].tolist()
            cbar.ax.set_yticks(cbar_ticks)
            cbar.ax.set_ylim([np.min(levels), np.max(levels)])

        except:
            cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                                values=levels, format=ticker.FuncFormatter(fmt))

            cbar_ticks = cbar.ax.get_yticks()
            cbar_ticks = cbar_ticks[cbar_ticks > 0].tolist()
            cbar.ax.set_yticks(cbar_ticks)
            cbar.ax.set_ylim([np.min(levels), np.max(levels)])

        # Plot colors for individual points so there is less ambiguity
        lgnd_check = False
        for loc_ref in range(len(results[:, 0])):
            if results[loc_ref, 0] > 0 and results[loc_ref, 0] < np.Inf:
                if not lgnd_check:
                    rgba = cmap(0)
                    # This does not have any color, it's just used for the legend
                    plt.plot(x_loc[loc_ref] / 1000, y_loc[loc_ref] / 1000,
                             marker='o', markerfacecolor=rgba[0:3],
                             markeredgecolor='k', markeredgewidth=1.5,
                             markersize=12, linestyle='none',
                             label='Wellbore with Nonzero Result')
                    lgnd_check = True
                    ncol_number += 1

                # Get the color by the value's proximity to the upper and lower
                # limits used for levels (which set the colorbar limits above)
                rgba = cmap((results[loc_ref, 0] - np.min(levels))
                             / (np.max(levels) - np.min(levels)))
                plt.plot(x_loc[loc_ref] / 1000, y_loc[loc_ref] / 1000,
                         marker='o', markerfacecolor=rgba[0:3],
                         markeredgecolor='k', markeredgewidth=1.5,
                         markersize=12, linestyle='none')

        if output_nm == 'pressure' and critPressureInput is not None:
            if var_type == 'noVariation':
                pressure_levels = np.array([critPressureInput])

                a, b = '{:.2e}'.format(critPressureInput).split('e')
                b = int(b)
                critPressure_str = r'${}\times10^{{{}}}$'.format(a, b)

                if np.max(results_temp[:, 0]) < critPressureInput:
                    title_pressure = ',\nNever Exceeded the P$_{crit}$ of ' + \
                        '{} MPa'.format(critPressure_str)

                elif np.min(results_temp[:, 0]) > critPressureInput:
                    title_pressure = ',\nAll Pressures Exceeded the P$_{crit}$ of ' + \
                        '{} MPa'.format(critPressure_str)

                elif np.min(results_temp[:, 0]) < critPressureInput and \
                        np.max(results_temp[:, 0]) >= critPressureInput:
                    title_pressure = ',\nCertain Pressures Exceeded the P$_{crit}$ of ' + \
                        '{} MPa'.format(critPressure_str)
                    Pcrit_Included = True

                    # The handle and label for this are created manually below
                    plt.tricontour(
                        x_loc_temp / 1000.0, y_loc_temp / 1000.0,
                        results_temp[:, 0], pressure_levels, colors = 'r')
                    ncol_number += 1

            elif var_type in ['strikeAndDip', 'LookupTable']:
                # Gets rid of any nan values (from a wellbore being on the injection site)
                x_loc_temp = x_loc[results[:, 0] > 0]
                y_loc_temp = y_loc[results[:, 0] > 0]
                critPressure_temp = critPressureInput[results[:, 0] > 0]
                results_temp = results[results[:, 0] > 0]

                # Gets rid of any Inf values (from a wellbore being on the injection site)
                x_loc_temp = x_loc_temp[results_temp[:, 0] < np.Inf]
                y_loc_temp = y_loc_temp[results_temp[:, 0] < np.Inf]
                critPressure_temp = critPressure_temp[results_temp[:, 0] > 0]
                results_temp = results_temp[results_temp[:, 0] < np.Inf]

                critPressureDiff = results_temp - critPressure_temp

                # Set the levels to zero, which will show where the pressure
                # equals the local critical pressure
                critPressureDiff_levels = np.array([0])

                if np.max(critPressureDiff[:, 0]) < 0:
                    title_pressure = ',\nNever Exceeded the Local P$_{crit}$ Values'

                elif np.min(critPressureDiff[:, 0]) > 0:
                    title_pressure = ',\nAll Pressures Exceeded the Local P$_{crit}$ Values'

                elif np.min(critPressureDiff[:, 0]) < 0 and \
                        np.max(critPressureDiff[:, 0]) >= 0:
                    title_pressure = ',\nCertain Pressures Exceeded the Local P$_{crit}$ Values'

                    Pcrit_Included = True

                    # The handle and label for this are created manually below
                    plt.tricontour(
                        x_loc_temp / 1000.0, y_loc_temp / 1000.0,
                        critPressureDiff[:, 0], critPressureDiff_levels, colors = 'r')
                    ncol_number += 1

    else:
        if enforce_levels and time_option:
            min_val = min_value
            max_val = max_value

            interval = ((max_val * MAX_VAL_ADJUST) - min_val) / 100
            levels = np.arange(min_val, (max_val * MAX_VAL_ADJUST) + interval,
                               interval)

            cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                                values=levels, format=ticker.FuncFormatter(fmt))

            cbar_ticks = cbar.ax.get_yticks()
            cbar_ticks = cbar_ticks[cbar_ticks > 0].tolist()
            cbar.ax.set_yticks(cbar_ticks)
            cbar.ax.set_ylim([np.min(levels), np.max(levels)])

            range_str = ''

        else:
            # If there are no results and no enforced colorbar, make the
            # colorbar go from 0 to 1.
            levels = np.arange(0, 1.01, 0.01)
            cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                                values=levels, format=ticker.FuncFormatter(fmt))

            cbar_ticks = cbar.ax.get_yticks()
            cbar_ticks = cbar_ticks[cbar_ticks > 0].tolist()
            cbar.ax.set_yticks(cbar_ticks)
            cbar.ax.set_ylim([np.min(levels), np.max(levels)])

            range_str = ''

    # Default cbar label is output name
    cbar_label = CBAR_LABELS.get(output_nm, output_nm)
    cbar.set_label(
        cbar_label, rotation=90, fontsize=axis_label_font_size,
        fontweight=selected_labelfontweight)

    if yaml_input['plot_injection_sites']:
        if len(InjectionCoordx) == 0 and yaml_input['InjectionCoordx'] is not None:
            # The lists will be empty for LookupTableReservoirs, use the .yaml input
            InjectionCoordx = yaml_input['InjectionCoordx']
            InjectionCoordy = yaml_input['InjectionCoordy']

        if isinstance(InjectionCoordx, list) and len(InjectionCoordx) > 0:
            for injRef, (xcoord_val, ycoord_val) in enumerate(
                    zip(InjectionCoordx, InjectionCoordy)):
                if injRef == 0:
                    plt.plot(xcoord_val / 1000, ycoord_val / 1000,
                             marker='s', color='k', linestyle='none',
                             markeredgewidth=2, markersize=6,
                             markerfacecolor='none', zorder=1e6,
                             label='Injection Site')
                    ncol_number += 1
                else:
                    plt.plot(xcoord_val / 1000, ycoord_val / 1000,
                             marker='s', color='k', linestyle='none',
                             markeredgewidth=2, markersize=6,
                             markerfacecolor='none', zorder=1e6)
        else:
            plt.plot(InjectionCoordx / 1000, InjectionCoordy / 1000,
                     marker='s', color='k', linestyle = 'none',
                     markeredgewidth=2, markersize=6,
                     markerfacecolor='none', zorder=1e6,
                     label='Injection Site')
            ncol_number += 1

    x_loc_copy = x_loc.tolist()
    y_loc_copy = y_loc.tolist()
    if yaml_input['plot_injection_sites']:
        x_loc_copy = x_loc.tolist()
        y_loc_copy = y_loc.tolist()
        if isinstance(InjectionCoordx, list) and len(InjectionCoordx) > 0:
            for (xcoord_val, ycoord_val) in zip(InjectionCoordx, InjectionCoordy):
                x_loc_copy.append(xcoord_val)
                y_loc_copy.append(ycoord_val)
        else:
            x_loc_copy.append(InjectionCoordx)
            y_loc_copy.append(InjectionCoordy)

    # These contain the OpenWellbore locations and the injection locations,
    # if the injection locations are being plotted.
    x_loc_copy = np.array(x_loc_copy)
    y_loc_copy = np.array(y_loc_copy)

    if grid_option:
        # These are used for the buffer room b/c they will not include the
        # injection location(s), only the grid-based OpenWellbore locations.
        x_vals_temp = np.unique(x_loc)
        y_vals_temp = np.unique(y_loc)
        plt.xlim((np.min(x_loc_copy) - ((x_vals_temp[1]
                                         - x_vals_temp[0]))) / 1000.0,
                  (np.max(x_loc_copy) + ((x_vals_temp[1]
                                          - x_vals_temp[0]))) / 1000.0)
        plt.ylim((np.min(y_loc_copy) - ((y_vals_temp[1]
                                         - y_vals_temp[0]))) / 1000.0,
                  (np.max(y_loc_copy) + ((y_vals_temp[1]
                                          - y_vals_temp[0]))) / 1000.0)
    else:
        xlim_adjust_val = (np.max(x_loc_copy) - np.min(x_loc_copy)) / 20
        ylim_adjust_val = (np.max(y_loc_copy) - np.min(y_loc_copy)) / 20
        plt.xlim((np.min(x_loc_copy) - xlim_adjust_val) / 1000.0,
                 (np.max(x_loc_copy) + xlim_adjust_val) / 1000.0)
        plt.ylim((np.min(y_loc_copy) - ylim_adjust_val) / 1000.0,
                 (np.max(y_loc_copy) + ylim_adjust_val) / 1000.0)

    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])

    ax.set_aspect('equal')
    ax.set_xlabel('Easting (km)', fontsize=axis_label_font_size,
                  fontweight=selected_labelfontweight)
    ax.set_ylabel('Northing (km)', fontsize=axis_label_font_size,
                  fontweight=selected_labelfontweight)

    # Create legend
    handle_list = []
    label_list = []
    handles, labels = ax.get_legend_handles_labels()

    for handle, label in zip(handles, labels):
        if label not in label_list:
            handle_list.append(handle)
            label_list.append(label)

    # If the metric is pressure and a critical pressure was given, include it in the legend
    if Pcrit_Included:
        if var_type == 'noVariation':
            critPressureLabel = 'P$_{crit}$'
        elif var_type in ['strikeAndDip', 'LookupTable']:
            critPressureLabel = 'P > Local P$_{crit}$'
        legend_element_critPressure = Line2D([0], [0], color='r',
                                             lw=2, label=critPressureLabel)

        handle_list.append(legend_element_critPressure)
        label_list.append(critPressureLabel)

    if ncol_number <= 2:
        bbox_val = (0.5, -0.1)
    elif ncol_number == 3:
        bbox_val = (0.55, -0.1)
    elif ncol_number >= 4:
        bbox_val = (0.6, -0.1)

    ax.legend(handle_list, label_list, fancybox=False, fontsize=gen_font_size - 2,
              ncol=ncol_number, edgecolor=[0, 0, 0], loc='upper center',
              bbox_to_anchor=bbox_val, framealpha=0.67)

    # aq_number is used in the figure title
    if output_nm in ['pressure', 'CO2saturation']:
        aq_number = ''

    # Get title corresponding to the output name
    if analysis == 'forward':
        if not time_option:
            # Need a space afterwards for the case where fig_title is set to
            # '{} at Each Point{}{}(Gray: Zero)'
            if range_str != '':
                comma_str = ', '
            else:
                comma_str = ' '

            fig_title = TITLE_OPTIONS_FORWARD.get(
                output_nm, '{} at Each Point{}{}(Gray: Zero){}'.format(
                    output_nm, comma_str, range_str, title_pressure))

            # No space afterwards b/c the comma_str goes right before a '\n'
            if range_str != '':
                comma_str = ','
            else:
                comma_str = ''

            if output_nm in TITLE_OPTIONS_FORWARD:
                fig_title = fig_title.format(aq_number, comma_str, range_str,
                                             title_pressure)
        else:
            # No space afterwards b/c the comma_str goes right before a '\n'
            if range_str != '':
                comma_str = ','
            else:
                comma_str = ''

            fig_title = TITLE_OPTIONS_FORWARD_T_INDEX.get(
                output_nm, '{} at Each Point at Time\nt = {} years{}{}(Gray: Zero){}'.format(
                    output_nm, time[time_index], comma_str, range_str, title_pressure))

            if output_nm in TITLE_OPTIONS_FORWARD_T_INDEX:
                fig_title = fig_title.format(aq_number, time[time_index],
                                             comma_str, range_str, title_pressure)

    elif analysis in ['lhs', 'parstudy']:
        realization_number = yaml_data['ModelParams']['Analysis']['siz']

        if range_str != '':
            comma_str = ', '
        else:
            comma_str = ' '

        if not time_option:
            fig_title = TITLE_OPTIONS_LHS.get(
                output_nm,
                ''.join(['Maximum {} at Each Point\nAcross {} ',
                         'LHS Simulations{}{}(Gray: Zero){}']).format(
                    output_nm, realization_number, comma_str, range_str, title_pressure))

            if output_nm in TITLE_OPTIONS_LHS:
                fig_title = fig_title.format(aq_number, realization_number,
                                             comma_str, range_str, title_pressure)
        else:
            fig_title = TITLE_OPTIONS_LHS_T_INDEX.get(
                output_nm,
                ''.join(['Maximum {} at Each Point at Time t = {} years\n',
                         'Across {} LHS Simulations{}{}(Gray: Zero){}']).format(
                    output_nm, time[time_index], realization_number,
                    comma_str, range_str, title_pressure))

            if output_nm in TITLE_OPTIONS_LHS_T_INDEX:
                fig_title = fig_title.format(aq_number, time[time_index],
                                             realization_number,
                                             comma_str, range_str, title_pressure)

    # Set figure title
    plt.title(fig_title, fontsize=title_font_size,
              fontweight=selected_labelfontweight)

    plt.subplots_adjust(left=0.001, bottom=0.15, right=0.901,
                        top=0.875, wspace=0.1, hspace=0.1)

    if title:
        plt.suptitle(title, fontweight=selected_labelfontweight,
                     fontsize=title_font_size)

    if savefig:
        output_dir = model_data['OutputDirectory']

        if save_results:
            # Set up data for .csv file
            if output_nm == 'pressure' and not critPressureInput is None:
                results_formatted = np.empty(((len(x_loc) + 1), 5), dtype=list)
            else:
                results_formatted = np.empty(((len(x_loc) + 1), 3), dtype=list)

            results_formatted[0, 0] = 'x (km)'
            results_formatted[0, 1] = 'y (km)'

            for row_ref, (x_val, y_val) in enumerate(zip(x_loc, y_loc)):
                results_formatted[row_ref + 1, 0] = str(x_val / 1000)
                results_formatted[row_ref + 1, 1] = str(y_val / 1000)

            # Get the column name depending on the output
            # If name does not match known observations use it as a column label
            results_formatted[0, 2] = CSV_FILE_COLUMNS.get(output_nm, output_nm)

            results_formatted[1:None, 2] = results[:, 0]

            if output_nm == 'pressure' and critPressureInput is not None:
                results_formatted[0, 3] = 'Critical Pressure [MPa]'

                if var_type == 'noVariation':
                    results_formatted[1:None, 3] = np.ones(len(results)) * critPressureInput

                elif var_type in ['strikeAndDip', 'LookupTable']:
                    results_formatted[1:None, 3] = critPressureInput[:, 0]

                results_formatted[0, 4] = 'Critical Pressure Exceeded'
                critPressureExceeded = np.zeros(len(results))

                for locRef in range(len(results)):
                    if var_type == 'noVariation':
                        if results[locRef, 0] >= critPressureInput:
                            critPressureExceeded[locRef] = 1

                    elif var_type in ['stirkeAndDip', 'LookupTable']:
                        if results[locRef, 0] >= critPressureInput[locRef, 0]:
                            critPressureExceeded[locRef] = 1

                results_formatted[1:None, 4] = critPressureExceeded

            if not os.path.exists(os.path.join(output_dir, 'csv_files')):
                os.mkdir(os.path.join(output_dir, 'csv_files'))

            file_name_addition = CSV_FILE_NAME_TAGS.get(output_nm, output_nm)

            if not time_option:
                filename = os.path.join(output_dir, 'csv_files',
                                        'AoR_{}.csv'.format(file_name_addition))
            else:
                filename = os.path.join(output_dir, 'csv_files',
                                        'AoR_{}_tIndex_{:.0f}.csv'.format(
                                            file_name_addition, time_index))

            # Save the ouput for the simulation to a .csv file
            with open(filename, 'w', newline='') as f:
                writer = csv.writer(f)
                for row_ref in range(len(x_loc) + 1):
                    writer.writerow(results_formatted[row_ref, :])

            f.close()

        if '.' in name:
            name_main = name[0:name.index('.')]
            name_main += '_tIndex_{:.0f}'.format(time_index)
            name_extension = name[name.index('.'):]
        else:
            name_main = name
            name_extension = '.png'

        if time_option:
            name_main += '_tIndex_{:.0f}'.format(time_index)

        plt.savefig(os.sep.join([output_dir,
                                 name_main + name_extension]), dpi=figdpi)
        plt.close()
    else:
        plt.show()


def not_boolean_debug_message_AoR(input_name, name, default_value):
    """
    Returns string delivering debug message regarding a variable not being
    of boolean type and setting it to the default value (True or False).

    input_name: string
    default_value: True or False
    """
    msg = ''.join(['Input provided for {} within the AoR plot {} was not of ',
                   'boolean type. Using the default value of {}.']).format(
                       input_name, name, default_value)
    return msg


def get_AoR_yaml_input(yaml_data, name, workflow_figure=False):
    """
    Function that reads the AoR plot input provided in the .yaml file
    and returns a dictionary containing the input. Note that InjectionCoordx
    and InjectionCoordy are only required if a LookupTableReservoir is being
    used. Other reservoir components will have locX and locY values that can
    be used to plot injection locations, but LookupTableReservoir components do
    not have locX and locY values.

    :param yaml_data: Dictionary of input values from the entire .yaml file
    :type yaml_data: dict

    :param name: Name of the AoR plot provided in the Plots section of the
        .yaml file
    :type name: str

    :returns: yaml_input
    """

    yaml_input_keys = [
        'dpi_input', 'plot_injection_sites', 'InjectionCoordx',
        'InjectionCoordy', 'SaveCSVFiles', 'TimeList', 'CriticalPressureMPa']

    InjectionCoord_debug_msg = ''.join([
        'InjectionCoord{} was provided for the AoR plot ', name,
        ', but not InjectionCoord{}. Check your input. Injection ',
        'sites will not be displayed.'])

    # Initialize output
    yaml_input = {key: None for key in yaml_input_keys}

    # Get shortcut to data to be analyzed
    if not workflow_figure:
        AoR_plot_data = yaml_data['Plots'][name]
    else:
        AoR_plot_data = yaml_data['Workflow']['Options']

    if AoR_plot_data is not None:
        if 'FigureDPI' in AoR_plot_data:
            yaml_input['dpi_input'] = AoR_plot_data['FigureDPI']

        if 'CriticalPressureMPa' in AoR_plot_data:
            if AoR_plot_data['CriticalPressureMPa'] == 'Calculated':
                yaml_input['CriticalPressureMPa'] = AoR_plot_data['CriticalPressureMPa']
            else:
                try:
                    yaml_input['CriticalPressureMPa'] = float(AoR_plot_data['CriticalPressureMPa'])
                except:
                    debug_msg = ''.join([
                        'The input provided for CriticalPressureMPa under the AoR plot ',
                        name, ' was not given as "Calculated" and it could not ',
                        'be turned into a float numeric value. Check your input. ',
                        'The CriticalPressureMPa entry will not be used.'])
                    logging.debug(debug_msg)
                    yaml_input['CriticalPressureMPa'] = None

        if 'TimeList' in AoR_plot_data:
            yaml_input['TimeList'] = AoR_plot_data['TimeList']

        if 'PlotInjectionSites' in AoR_plot_data:
            yaml_input['plot_injection_sites'] = AoR_plot_data[
                'PlotInjectionSites']
            if not isinstance(yaml_input['plot_injection_sites'], bool):
                debug_msg = not_boolean_debug_message_AoR(
                    'PlotInjectionSites', name, False)
                logging.debug(debug_msg)
                yaml_input['plot_injection_sites'] = False

        if 'InjectionCoordx' in AoR_plot_data:
            try:
                yaml_input['InjectionCoordx'] = float(AoR_plot_data['InjectionCoordx'])
            except TypeError:
                yaml_input['InjectionCoordx'] = list(AoR_plot_data['InjectionCoordx'])

        if 'InjectionCoordy' in AoR_plot_data:
            try:
                yaml_input['InjectionCoordy'] = float(AoR_plot_data['InjectionCoordy'])
            except TypeError:
                yaml_input['InjectionCoordy'] = list(AoR_plot_data['InjectionCoordy'])

            if yaml_input['InjectionCoordx'] is None:
                debug_msg = InjectionCoord_debug_msg.format('y', 'x')
                logging.debug(debug_msg)
                yaml_input['plot_injection_sites'] = False

            elif isinstance(yaml_input['InjectionCoordx'], list) and \
                    isinstance(yaml_input['InjectionCoordy'], list):
                if len(yaml_input['InjectionCoordx']) != len(
                        yaml_input['InjectionCoordy']):
                    debug_msg = ''.join([
                        'The InjectionCoordy provided for the AoR plot ', name,
                        ' was not of the same length as the InjectionCoordx ',
                        'provided. Check your input. Injection sites will not ',
                        'be displayed.'])
                    logging.debug(debug_msg)
                    yaml_input['plot_injection_sites'] = False

        if yaml_input['InjectionCoordx'] is not None and yaml_input[
                'InjectionCoordy'] is None:
            debug_msg = InjectionCoord_debug_msg.format('x', 'y')
            logging.debug(debug_msg)
            yaml_input['plot_injection_sites'] = False

        if 'SaveCSVFiles' in AoR_plot_data:
            yaml_input['SaveCSVFiles'] = AoR_plot_data['SaveCSVFiles']

            if not isinstance(yaml_input['SaveCSVFiles'], bool):
                debug_msg = not_boolean_debug_message_AoR(
                    'SaveCSVFiles', name, True)
                logging.debug(debug_msg)
                yaml_input['SaveCSVFiles'] = None

    return yaml_input


def get_t_indices(time_list, time_array):
    """
    Returns the time index corresponding to the time_array value closest to
    the times in time_list. Both time_list and time_array have units of years.
    """
    time_index_list = []

    corr_t_index = np.arange(0, len(time_array))

    if not isinstance(time_list, list):
        time_list = [time_list]

    for time in time_list:
        abs_diff = np.zeros(len(time_array))

        for t_ref, t_val in enumerate(time_array):
            abs_diff[t_ref] = np.abs(time - (t_val))

        closest_t_index = corr_t_index[abs_diff == np.min(abs_diff)]

        # If more than one time_array value had the same distance to the time
        # in question (e.g., one value is 0.5 yrs before and another is 0.5 yrs
        # after), then just pick one.
        if len(closest_t_index) > 1:
            closest_t_index = closest_t_index[-1]

        if isinstance(closest_t_index, list):
            closest_t_index = closest_t_index[0]

        time_index_list.append(int(closest_t_index))

    return time_index_list


def get_crit_pressure(output_component, sm=None, yaml_data=None):
    """
    This function calculates the critical pressure for an OpenWellbore component.
    """
    wellTop = iamcommons.get_parameter_val(output_component, 'wellTop',
                                           sm=sm, yaml_data=yaml_data)

    reservoirDepth = iamcommons.get_parameter_val(output_component, 'reservoirDepth',
                                                  sm=sm, yaml_data=yaml_data)

    brineDensity = iamcommons.get_parameter_val(output_component, 'brineDensity',
                                                sm=sm, yaml_data=yaml_data)

    critPressureVal = (wellTop * GRAV_ACCEL * WATER_DENSITY) + (
        brineDensity * GRAV_ACCEL * (reservoirDepth - wellTop))

    return critPressureVal


def plot_aor_workflow_results(yaml_data, sm, All_x_points_km, All_y_points_km,
                              AoR_point_included, time_index=None, analysis='forward',
                              pressure_included=True, figsize=(10, 8),
                              bold_labels=True, figdpi=100, gen_font_size=12,
                              axis_label_font_size=14, title_font_size=14):
    """
    Makes a plot for the AoR Workflow. The AoR shown in this plot reflects
    multiple metrics (e.g., pressure, CO2saturation, pH_Volume, and TDS_volume).
    """
    output_dir = yaml_data['ModelParams']['OutputDirectory']

    if bold_labels:
        selected_labelfontweight = 'bold'
    else:
        selected_labelfontweight = 'normal'

    if 'AoRFigureName' in yaml_data['Workflow']['Options']:
        name_main = yaml_data['Workflow']['Options']['AoRFigureName']
    else:
        name_main = 'AoR_Workflow_Plot'

    if '.' in name_main:
        name_extension = name_main[name_main.index('.'):]
        name_main = name_main[0:name_main.index('.')]
    else:
        name_extension = '.png'

    if time_index is None:
        fig_title = 'Area of Review Based on {}{}'
        name_addition = ''
    else:
        time = sm.time_array / 365.25
        fig_title = 'Area of Review Based on {}{},' + ' t = {} years'.format(
            time[time_index])
        name_addition = '_tIndex_{}'.format(time_index)

    analysis = yaml_data['ModelParams']['Analysis']

    if yaml_data['s'] is not None: # lhs or parstudy was run
        # Find number of run realizations
        realizations = yaml_data['s'].samples.values.shape[0]
    else: # for forward simulations
        realizations = 1

    if 'type' in analysis:
        analysis = analysis['type']

    if pressure_included:
        pressure_inclusion_str = 'All Metrics'
    else:
        pressure_inclusion_str = 'All Metrics but Pressure'

    if analysis in ['lhs', 'parstudy']:
        fig_title = fig_title.format(
            pressure_inclusion_str, '\nAcross {} Realizations'.format(realizations))
    else:
        fig_title = fig_title.format(pressure_inclusion_str, '')

    yaml_input = get_AoR_yaml_input(yaml_data, name_main, workflow_figure=True)

    components = list(sm.component_models.values())
    # If injection sites need to be plotted, get the injection sites
    if yaml_input['plot_injection_sites']:
        InjectionCoordx = []
        InjectionCoordy = []

        for comp in components:
            if comp.class_type in AOR_RESERVOIR_COMPONENTS:
                if comp.class_type != 'LookupTableReservoir':
                    InjectionCoordx.append(comp.injX)
                    InjectionCoordy.append(comp.injY)

        InjectionCoordx = np.unique(InjectionCoordx).tolist()
        InjectionCoordy = np.unique(InjectionCoordy).tolist()

    else:
        InjectionCoordx = None
        InjectionCoordy = None

    if not yaml_input['dpi_input'] is None:
        figdpi = yaml_input['dpi_input']

    # Figures
    font = RC_FONT
    font['size'] = gen_font_size
    plt.rc('font', **font)

    ncol_number = 0

    # AoR Figure
    _ = plt.figure(figsize=figsize, dpi=figdpi)

    ax = plt.gca()
    ax.set_facecolor(BACKGROUND_COLOR)

    plt.plot(All_x_points_km, All_y_points_km, linestyle='none',
             marker='o', color='k', markeredgewidth=1, markersize=5,
             markerfacecolor='none', label='All Points Considered', zorder=1)
    ncol_number += 1

    plt.plot(All_x_points_km[AoR_point_included == 1],
             All_y_points_km[AoR_point_included == 1],
             marker='o', markerfacecolor='red', markeredgecolor='k',
             markeredgewidth=1.5, markersize=12, linestyle='none',
             label='Points included in the Area of Review', zorder=2)
    ncol_number += 1

    if yaml_input['plot_injection_sites']:
        if len(InjectionCoordx) == 0 and yaml_input['InjectionCoordx'] is not None:
            # The lists will be empty for LookupTableReservoirs, use the .yaml input
            InjectionCoordx = yaml_input['InjectionCoordx']
            InjectionCoordy = yaml_input['InjectionCoordy']

        if isinstance(InjectionCoordx, list) and len(InjectionCoordx) > 0:
            for injRef, (xcoord_val, ycoord_val) in enumerate(
                    zip(InjectionCoordx, InjectionCoordy)):
                if injRef == 0:
                    plt.plot(xcoord_val / 1000, ycoord_val / 1000,
                             marker='s', color='k', linestyle='none',
                             markeredgewidth=2, markersize=6,
                             markerfacecolor='none', zorder=10,
                             label='Injection Site')
                    ncol_number += 1
                else:
                    plt.plot(xcoord_val / 1000, ycoord_val / 1000,
                             marker='s', color='k', linestyle='none',
                             markeredgewidth=2, markersize=6,
                             markerfacecolor='none', zorder=10)
        else:
            plt.plot(InjectionCoordx / 1000, InjectionCoordy / 1000,
                     marker='s', color='k', linestyle = 'none',
                     markeredgewidth=2, markersize=6,
                     markerfacecolor='none', zorder=10,
                     label='Injection Site')
            ncol_number += 1

    grid_option = False
    if 'OpenWellbore1' in yaml_data:
        if 'Locations' in yaml_data['OpenWellbore1']:
            if 'grid' in yaml_data['OpenWellbore1']['Locations']:
                grid_option = True

    if grid_option:
        x_vals_temp = np.unique(All_x_points_km)
        y_vals_temp = np.unique(All_y_points_km)

        plt.xlim(np.min(All_x_points_km) - ((x_vals_temp[1] - x_vals_temp[0])),
                 np.max(All_x_points_km) + ((x_vals_temp[1] - x_vals_temp[0])))

        plt.ylim(np.min(All_y_points_km) - ((y_vals_temp[1] - y_vals_temp[0])),
                 np.max(All_y_points_km) + ((y_vals_temp[1] - y_vals_temp[0])))
    else:
        xlim_adjust_val = (np.max(All_x_points_km) - np.min(All_x_points_km)) / 20
        ylim_adjust_val = (np.max(All_y_points_km) - np.min(All_y_points_km)) / 20

        plt.xlim((np.min(All_x_points_km) - xlim_adjust_val),
                 (np.max(All_x_points_km) + xlim_adjust_val))
        plt.ylim((np.min(All_y_points_km) - ylim_adjust_val),
                 (np.max(All_y_points_km) + ylim_adjust_val))

    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])

    ax.set_aspect('equal')
    ax.set_xlabel('Easting (km)', fontsize=axis_label_font_size,
                  fontweight=selected_labelfontweight)
    ax.set_ylabel('Northing (km)', fontsize=axis_label_font_size,
                  fontweight=selected_labelfontweight)

    ax.legend(fancybox=False, fontsize=gen_font_size - 2, ncol=ncol_number,
              edgecolor=[0, 0, 0], loc='upper center', bbox_to_anchor=(0.5, -0.1),
              framealpha=0.67)

    # Set figure title
    plt.title(fig_title, fontsize=title_font_size,
              fontweight=selected_labelfontweight)

    plt.subplots_adjust(left=0.051, bottom=0.15, right=0.951,
                        top=0.875, wspace=0.1, hspace=0.1)

    plt.savefig(os.sep.join([output_dir, name_main + name_addition + name_extension]),
                dpi=figdpi)
    plt.close()
