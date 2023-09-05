"""
Code to create figures illustrating the gridded observation produced by components
such as Seal Horizon and Fault Flow.

Examples illustrating setup of GriddedMetric plot:
    ControlFile_ex18.yaml
    ControlFile_ex19.yaml
    ControlFile_ex23.yaml

@author: Nate Mitchell (Nathaniel.Mitchell@NETL.DOE.GOV)
"""

import os
import csv
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import ticker

from .label_setup import Y_LABEL_DICT
from openiam.visualize import time_series
import openiam.cfi.commons as iamcommons


RC_FONT = {'family': 'Arial', 'weight': 'normal', 'size': None}

# Applicable components
GRIDDED_METRIC_COMPONENTS = ['SealHorizon', 'FaultFlow']

# Components where the gridded output represents cells that each have a specified area
AREA_COMPONENTS = ['SealHorizon']

# Components where the gridded output represents fault segments that each have
# a starting position, strike, dip, and length
FAULT_COMPONENTS = ['FaultFlow']

RESERVOIR_COMPONENTS = ['LookupTableReservoir',
                        'SimpleReservoir',
                        'AnalyticalReservoir']

BACKGROUND_COLOR = [0.67, 0.67, 0.67]

# The formats for each metric type is currently the same, but using this dictionary
# allows the formats for other metric types to be added more easily in future updates.
FILENAME_OPTIONS = {
    'CO2_aquifer': '{comp_name}_{metric_name}_sim_{sim_index}_time_{t_index}.npz',
    'brine_aquifer': '{comp_name}_{metric_name}_sim_{sim_index}_time_{t_index}.npz',
    'mass_CO2_aquifer': '{comp_name}_{metric_name}_sim_{sim_index}_time_{t_index}.npz',
    'mass_brine_aquifer': '{comp_name}_{metric_name}_sim_{sim_index}_time_{t_index}.npz',
    }

DEFAULT_FILENAME = '{comp_name}_{metric_name}_sim_{sim_index}_time_{t_index}.npz'

COLORMAP_OPTIONS = {'CO2_aquifer': 'plasma',
                    'mass_CO2_aquifer': 'plasma',
                    'brine_aquifer': 'viridis',
                    'mass_brine_aquifer': 'viridis'}

DEFAULT_COLORMAP = 'cividis'

# The units need a space before them for formatting purposes
UNIT_STR_OPTIONS = {
    'CO2_aquifer': ' kg/s',
    'mass_CO2_aquifer': ' kg',
    'brine_aquifer': ' kg/s',
    'mass_brine_aquifer': ' kg'
    }

MASS_METRICS = ['mass_CO2_aquifer', 'mass_brine_aquifer']
RATE_METRICS = ['CO2_aquifer', 'brine_aquifer']

CSV_FILE_COLUMNS = {
    'CO2_aquifer': 'CO2 leakage rate to aquifer (kg/s)',
    'mass_CO2_aquifer': 'Mass of CO2 leaked to aquifer (kg)',
    'brine_aquifer': 'Brine leakage rate to aquifer (kg/s)',
    'mass_brine_aquifer': 'Mass of brine leaked to aquifer (kg)'}

DEFAULT_MIN = 9.0e99
DEFAULT_MAX = -9.0e99


def gridded_metric_plot(yaml_data, model_data, sm, s, output_dir,
                        name='GriddedObs_Figure1', analysis='forward',
                        savefig=None, figsize=(10, 8), genfontsize=10,
                        axislabelfontsize=14, titlefontsize=14, boldlabels=True,
                        colormap=None, figure_dpi=100, min_value=0,
                        metricMarkerSize=12):
    """
    Makes a map-view figure showing the values of gridded metrics.
    """
    os.chdir(output_dir)

    GriddedMetric_yaml_input_dict = get_gridded_plot_yaml_input(
        yaml_data, name, plot_type='GriddedMetric')

    if analysis =='forward':
        sim_index = 0

    elif analysis in ['lhs', 'parstudy']:
        sim_index = GriddedMetric_yaml_input_dict['Realization']

        if sim_index is None:
            # For LHS or parstudy simulations, the gridded data from the
            # SealHorizon and FaultFlow components are saved with indices of 1
            # and higher. An index of 0 does not work.
            sim_index = 1

    components_name_list = GriddedMetric_yaml_input_dict['comp_name_list']

    metric_name = GriddedMetric_yaml_input_dict['metric_name']

    if colormap is None:
        colormap = COLORMAP_OPTIONS.get(metric_name, DEFAULT_COLORMAP)

    cmap = plt.cm.get_cmap(colormap)

    time_list = GriddedMetric_yaml_input_dict['TimeList']

    figure_dpi = GriddedMetric_yaml_input_dict['FigureDPI']

    plot_injection_sites = GriddedMetric_yaml_input_dict['plot_injection_sites']

    InjectionCoordx = GriddedMetric_yaml_input_dict['InjectionCoordx']
    InjectionCoordy = GriddedMetric_yaml_input_dict['InjectionCoordy']

    EnforceXandYLims = GriddedMetric_yaml_input_dict['EnforceXandYLims']
    xLims = GriddedMetric_yaml_input_dict['xLims']
    yLims = GriddedMetric_yaml_input_dict['yLims']

    equal_axes = GriddedMetric_yaml_input_dict['EqualAxes']

    plot_val_over_area = GriddedMetric_yaml_input_dict['PlotOverAreas']

    saveCSVFiles = GriddedMetric_yaml_input_dict['SaveCSVFiles']

    unit_str = UNIT_STR_OPTIONS.get(metric_name, '')

    # Currently, this plot type only handles mass metrics (e.g., mass_brine_aquifer)
    # or rate metrics (e.g., brine_aquifer), but I'm setting this up so other
    # types of metrics can be added in future updates.
    mass_metric = False
    rate_metric = False
    if metric_name in MASS_METRICS:
        mass_metric = True
    elif metric_name in RATE_METRICS:
        rate_metric = True

    time_array = sm.time_array

    num_times = len(time_array)

    if time_list == 'All':
        time_index_list = range(0, num_times)
    else:
        time_index_list = get_t_indices(time_list, time_array)

    components_to_use, component_types, component_cell_areas, \
        component_num_parts, component_xvals, component_yvals, \
            component_length, component_strike, component_dip, \
                res_comp_injX, res_comp_injY, x_range, y_range = \
                get_comps_and_limits(components_name_list, yaml_data, name, sm)

    if plot_injection_sites and InjectionCoordx is None:
        InjectionCoordx = res_comp_injX
        InjectionCoordy = res_comp_injY

    component_min_metric, component_max_metric = get_max_vals_over_time(
        components_to_use, num_times, metric_name, sim_index, min_value=min_value)

    # Check if the component types are accepted
    comp_type_check = True
    for comp in components_to_use:
        if not comp.class_type in GRIDDED_METRIC_COMPONENTS:
            comp_type_check = False

    if comp_type_check:
        if analysis == 'forward':
            file_name = '{}_{}_t{}{}'
            file_name_no_results = '{}_{}_All_t{}'

        elif analysis in ['lhs', 'parstudy']:
            file_name = '{}_{}_Sim{}_t{}{}'
            file_name_no_results = '{}_{}_Sim{}_All_t{}'

    # Not using this value right now, but leaving it for potential updates
    # min_value = min_value_updated

    # Check the file name for an extension like .tiff or .eps
    if '.' in name:
        name_extension = name[name.index('.'):None]
    else:
        name_extension = '.png'

    metric_vals_compiled = []
    x_vals_compiled = []
    y_vals_compiled = []
    t_vals_compiled = []

    plot_faults = False
    for compRef, compType in enumerate(component_types):
        comp_name = components_to_use[compRef].name

        if compType in FAULT_COMPONENTS:
            plot_faults = True
            plot_val_over_area = False

        elif compType in AREA_COMPONENTS:
            if plot_val_over_area:
                if GriddedMetric_yaml_input_dict['CellLengthX'] is not None and \
                        GriddedMetric_yaml_input_dict['CellLengthY'] is not None:
                    cell_length_x_km = GriddedMetric_yaml_input_dict[
                        'CellLengthX'] / 1000
                    cell_length_y_km = GriddedMetric_yaml_input_dict[
                        'CellLengthY'] / 1000
                else:
                    cell_length_x_km = (component_cell_areas[compRef] ** 0.5) / 1000
                    cell_length_y_km = (component_cell_areas[compRef] ** 0.5) / 1000

        no_results_check = False
        if np.max(component_max_metric[compRef]) <= min_value:
            no_results_check = True
            interval = (1 - min_value) / 100
            levels = np.arange(min_value, 1 + interval, interval)

        else:
            interval = (component_max_metric[compRef]
                        - component_min_metric[compRef]) / 100
            levels = np.arange(component_min_metric[compRef],
                               component_max_metric[compRef] + interval,
                               interval)

        if EnforceXandYLims:
            x_axis_limits = [xLims[0] / 1000, xLims[1] / 1000]
            y_axis_limits = [yLims[0] / 1000, yLims[1] / 1000]

        else:
            # These are used to adjust the x and y axis limits
            x_buffer = (x_range[1] - x_range[0]) / 20
            y_buffer = (y_range[1] - y_range[0]) / 20

            if x_buffer == 0:
                x_buffer = 10

            if y_buffer == 0:
                y_buffer = 10

            x_axis_limits = [(x_range[0] - x_buffer) / 1000,
                              (x_range[1] + x_buffer) / 1000]
            y_axis_limits = [(y_range[0] - y_buffer) / 1000,
                              (y_range[1] + y_buffer) / 1000]

        use_sci_formatting = False
        # If the distances are very large, use scientific notation to avoid
        # the ugly formatting used by default.
        if np.max(x_axis_limits) >= 1000 or np.max(y_axis_limits) >= 1000:
            use_sci_formatting = True

        if no_results_check:
            fig, ax = figure_setup(
                sim_index, 1, levels, comp_name, metric_name, x_axis_limits,
                y_axis_limits, genfontsize=genfontsize, axislabelfontsize=axislabelfontsize,
                titlefontsize=titlefontsize, boldlabels=boldlabels, figsize=figsize,
                colormap=colormap, equal_axes=equal_axes, min_value=min_value)

            if analysis == 'forward':
                title_str = 'Gridded Results, All Model Times'
            elif analysis in ['lhs', 'parstudy']:
                title_str = 'Gridded Results for Simulation {},'.format(
                    sim_index) + '\nAll Model Times'

            if min_value > 0:
                a, b = '{:.2e}'.format(min_value).split('e')
                b = int(b)
                min_val_str = r'$\leq {}\times10^{{{}}}$'.format(a, b)

                title_str += ' (Gray: Results {})'.format(
                    min_val_str)
            else:
                title_str += ' (Gray: Results = {})'.format(
                    min_value)

            plot_comp_inj_sites(
                ax, comp_name, compType, component_xvals[compRef], component_yvals[compRef],
                plot_injection_sites, InjectionCoordx, InjectionCoordy,
                genfontsize=genfontsize, include_labels=True, plot_faults=plot_faults,
                component_length=component_length[compRef],
                component_strike=component_strike[compRef],
                component_dip=component_dip[compRef])

            plt.title(title_str, fontsize = titlefontsize, fontweight = 'bold')

            if use_sci_formatting:
                ax.ticklabel_format(
                    style='sci', axis='both', useOffset=False,
                    scilimits=(0, 0), useMathText=True)

            if analysis == 'forward':
                if savefig:
                    plt.savefig(os.sep.join([output_dir, file_name_no_results.format(
                        comp_name, metric_name, name_extension)]),
                        dpi=figure_dpi)
                    plt.close()
                else:
                    fig.show()

            elif analysis in ['lhs', 'parstudy']:
                if savefig:
                    plt.savefig(os.sep.join([output_dir, file_name_no_results.format(
                        comp_name, metric_name, sim_index,
                        name_extension)]), dpi=figure_dpi)
                    plt.close()
                else:
                    fig.show()

        else:
            # case for no_results_check == False
            filename = FILENAME_OPTIONS.get(metric_name, DEFAULT_FILENAME)

            for t_index in time_index_list:
                if rate_metric or mass_metric:
                    # Reset the sum to 0 at each time step.
                    metricSum = 0

                fig, ax = figure_setup(
                    sim_index, t_index, levels, comp_name, metric_name,
                    x_axis_limits, y_axis_limits, genfontsize=genfontsize,
                    axislabelfontsize=axislabelfontsize, titlefontsize=titlefontsize,
                    boldlabels=boldlabels, figsize=figsize, colormap=colormap,
                    equal_axes=equal_axes, min_value=min_value)

                if analysis == 'forward':
                    title_str = 'Gridded Results for t = {} years,\n'.format(
                        time_array[t_index] / 365.25)
                elif analysis in ['lhs', 'parstudy']:
                    title_str = 'Simulation {}, Gridded Results for t = {} years,\n'.format(
                        sim_index, time_array[t_index] / 365.25)

                plot_comp_inj_sites(
                    ax, comp_name, compType, component_xvals[compRef], component_yvals[compRef],
                    plot_injection_sites, InjectionCoordx, InjectionCoordy,
                    genfontsize=genfontsize, include_labels=True, plot_faults=plot_faults,
                    component_length=component_length[compRef],
                    component_strike=component_strike[compRef],
                    component_dip=component_dip[compRef])

                if comp_type_check:
                    metric_file_data = np.load(filename.format(
                        comp_name=comp_name, metric_name=metric_name,
                        sim_index=sim_index, t_index=t_index))
                    metric = metric_file_data['data']
                    metric_file_data.close()
                else:
                    # Leaving this as a placeholder for obtaining the output for
                    # future components that produce gridded data.
                    pass

                min_val_actual = None
                max_val = None
                # The parts are cells for AREA_COMPONENTS or segments for FAULT_COMPONENTS
                for partRef, cellMetric in enumerate(metric):

                    x_center_km = component_xvals[compRef][partRef] / 1000
                    y_center_km = component_yvals[compRef][partRef] / 1000

                    if plot_val_over_area and compType in AREA_COMPONENTS:
                        node_x = [x_center_km - (cell_length_x_km / 2),
                                  x_center_km + (cell_length_x_km / 2),
                                  x_center_km + (cell_length_x_km / 2),
                                  x_center_km - (cell_length_x_km / 2)]
                        node_y = [y_center_km - (cell_length_y_km / 2),
                                  y_center_km - (cell_length_y_km / 2),
                                  y_center_km + (cell_length_y_km / 2),
                                  y_center_km + (cell_length_y_km / 2)]

                    if np.max(cellMetric) > min_value:
                        metric_vals = cellMetric[cellMetric > min_value]

                        for metricRef in range(0, len(metric_vals)):
                            metric_vals_compiled.append(metric_vals[metricRef])
                            x_vals_compiled.append(x_center_km * 1000)
                            y_vals_compiled.append(y_center_km * 1000)
                            t_vals_compiled.append(time_array[t_index] / 365.25)

                            if mass_metric or rate_metric:
                                metricSum += metric_vals[metricRef]

                            # Get the color by the value's proximity to the upper and
                            # lower limits used for levels (which set the colorbar limits).
                            rgba = cmap(
                                (metric_vals[metricRef] - np.min(levels))
                                / (np.max(levels) - np.min(levels)))

                            if plot_val_over_area:
                                plt.fill(node_x, node_y, color=rgba[0:3], zorder=1)

                            else:
                                plt.plot(x_center_km, y_center_km, color=rgba[0:3],
                                         marker='o', markerfacecolor=rgba[0:3],
                                         linewidth=1, linestyle='none',
                                         markersize=metricMarkerSize, zorder=1)

                        if min_val_actual is None:
                            min_val_actual = np.min(metric_vals)
                        elif min_val_actual > np.min(metric_vals):
                            min_val_actual = np.min(metric_vals)

                        if max_val is None:
                            max_val = np.max(metric_vals)
                        elif max_val < np.max(metric_vals):
                            max_val = np.max(metric_vals)

                metricSum_str = ''
                if mass_metric or rate_metric:
                    if metricSum == 0:
                        metricSum_str = '0' + unit_str
                    else:
                        a, b = '{:.2e}'.format(metricSum).split('e')
                        b = int(b)
                        metricSum_str = r'${}\times10^{{{}}}$'.format(a, b) + unit_str

                    metricSum_str = 'Total: {} '.format(metricSum_str)

                if min_value == 0:
                    min_val_str = '= 0' + unit_str
                else:
                    a, b = '{:.2e}'.format(min_value).split('e')
                    b = int(b)
                    min_val_str = r'$\leq {}\times10^{{{}}}$'.format(a, b) + unit_str

                if min_val_actual is not None and max_val is not None:
                    a, b = '{:.2e}'.format(min_val_actual).split('e')
                    b = int(b)
                    min_val_actual_str = r'${}\times10^{{{}}}$'.format(a, b) + unit_str

                    a, b = '{:.2e}'.format(max_val).split('e')
                    b = int(b)
                    max_val_str = r'${}\times10^{{{}}}$'.format(a, b) + unit_str

                    if min_val_actual == max_val:
                        title_str_edit = title_str + \
                            'Value: {},\n{}(Gray: Results {})'.format(
                                max_val_str, metricSum_str, min_val_str)
                    else:
                        title_str_edit = title_str + \
                            'Range: {} to {},\n{}(Gray: Results {})'.format(
                                min_val_actual_str, max_val_str,
                                metricSum_str, min_val_str)
                else:
                    if min_value == 0:
                        min_val_str = '= 0' + unit_str
                    else:
                        a, b = '{:.2e}'.format(min_value).split('e')
                        b = int(b)
                        min_val_str = r'$\leq {}\times10^{{{}}}$'.format(a, b) + unit_str

                    title_str_edit = title_str + '{}(Gray: Results {})'.format(
                        metricSum_str, min_val_str)

                plt.title(title_str_edit, fontsize = titlefontsize,
                          fontweight = 'bold')

                if use_sci_formatting:
                    ax.ticklabel_format(
                        style='sci', axis='both', useOffset=False,
                        scilimits=(0, 0), useMathText=True)

                if analysis == 'forward':
                    if savefig:
                        plt.savefig(os.sep.join([output_dir, file_name.format(
                            comp_name, metric_name, t_index,
                            name_extension)]), dpi=figure_dpi)
                        plt.close()
                    else:
                        fig.show()

                elif analysis in ['lhs', 'parstudy']:
                    if savefig:
                        plt.savefig(os.sep.join([output_dir, file_name.format(
                            comp_name, metric_name, sim_index,
                            t_index, name_extension)]), dpi=figure_dpi)
                        plt.close()
                    else:
                        fig.show()

        if saveCSVFiles:
            save_results_to_csv(
                output_dir, comp_name, analysis, sim_index, metric_name,
                x_vals_compiled, y_vals_compiled, t_vals_compiled,
                metric_vals_compiled, component_xvals[compRef],
                component_yvals[compRef], time_array)


def get_max_vals_over_time(components_to_use, num_times, metric_name, sim_index,
                           min_value=0):
    """
    Evaluates the distribution of gridded observations over time and returns
    the minimum and maximum values. These values are used to set the
    colorbar limits; otherwise, the colorbars for different model times would
    change, which impedes comparison.
    """
    component_min_metric = np.ones((len(components_to_use))) * DEFAULT_MIN
    component_max_metric = np.ones((len(components_to_use))) * DEFAULT_MAX

    filename = FILENAME_OPTIONS.get(metric_name, DEFAULT_FILENAME)

    for compRef, comp in enumerate(components_to_use):
        comp_name = comp.name

        for t_index in range(num_times):
            if comp.class_type in GRIDDED_METRIC_COMPONENTS:
                metric_file_data = np.load(filename.format(
                    comp_name=comp_name, metric_name=metric_name,
                    sim_index=sim_index, t_index=t_index))
                metric = metric_file_data['data']
                metric_file_data.close()
            else:
                # Leaving this as a placeholder for other component types
                pass

            metric_vals = metric[metric > min_value]

            if len(metric_vals) > 0:
                if component_min_metric[compRef] > np.min(metric_vals):
                    component_min_metric[compRef] = np.min(metric_vals)

                if component_max_metric[compRef] < np.max(metric_vals):
                    component_max_metric[compRef] = np.max(metric_vals)

    component_min_metric[component_min_metric == DEFAULT_MIN] = min_value
    component_max_metric[component_max_metric == DEFAULT_MAX] = min_value

    return component_min_metric, component_max_metric


def get_comps_and_limits(components_name_list, yaml_data, name, sm):
    """
    Function that produces lists of components to use and their characteristics
    (e.g., x and y values and the index for the corresponding aquifer, if
    applicable). By going through these components, it also obtains and provides
    the minimum x and y values used in the simulation(s).
    """
    # These are overwritten within the loop below
    min_x_val = DEFAULT_MIN
    max_x_val = DEFAULT_MAX
    min_y_val = DEFAULT_MIN
    max_y_val = DEFAULT_MAX

    components_to_use = []
    component_types = []

    # Areas for components in the AREA_COMPONENTS, None values for FAULT_COMPONENTS
    component_cell_areas = []

    # Number of segments for FAULT_COMPONENTS, number of cells for AREA_COMPONENTS
    component_num_parts = []

    # The x and y values at the center of AREA_COMPONENTS, the xStart and yStart
    # values for FAULT_COMPONENTS
    component_xvals = []
    component_yvals = []

    # Data needed for FAULT_COMPONENTS, None values for AREA_COMPONENTS
    component_strike = []
    component_dip = []
    component_length = []

    res_comp_injX = []
    res_comp_injY = []

    components = list(sm.component_models.values())
    for comp in components:
        if comp.class_type in GRIDDED_METRIC_COMPONENTS:
            for compRef, compName in enumerate(components_name_list):
                if compName in comp.name:
                    components_to_use.append(comp)
                    component_types.append(comp.class_type)

                    if comp.class_type in AREA_COMPONENTS:
                        cell_area = iamcommons.get_parameter_val(comp, 'area')
                        component_cell_areas.append(cell_area)

                        cell_loc = comp.cell_xy_centers

                        # The rows are for different cells, the two colums are for
                        # the x and y values of each cell.
                        num_cells = cell_loc.shape[0]

                        cell_coordx = cell_loc[:, 0]
                        cell_coordy = cell_loc[:, 1]

                        component_num_parts.append(num_cells)

                        component_xvals.append([])
                        component_yvals.append([])
                        for cellRef in range(num_cells):
                            component_xvals[-1].append(cell_coordx[cellRef])
                            component_yvals[-1].append(cell_coordy[cellRef])

                        component_length.append(None)
                        component_strike.append(None)
                        component_dip.append(None)

                        GriddedMetric_yaml_input_dict = get_gridded_plot_yaml_input(
                            yaml_data, name, plot_type='GriddedMetric')

                        plot_val_over_area = GriddedMetric_yaml_input_dict['PlotOverAreas']
                        # Update the min and max x and y to reflect the cell dimensions,
                        # if using plot_val_over_area
                        if plot_val_over_area:
                            if GriddedMetric_yaml_input_dict['CellLengthX'] is not None and \
                                    GriddedMetric_yaml_input_dict['CellLengthY'] is not None:
                                cell_length_x_m = GriddedMetric_yaml_input_dict[
                                    'CellLengthX']
                                cell_length_y_m = GriddedMetric_yaml_input_dict[
                                    'CellLengthY']
                            else:
                                cell_length_x_m = (component_cell_areas[compRef] ** 0.5) / 1000
                                cell_length_y_m = (component_cell_areas[compRef] ** 0.5) / 1000

                            dx = cell_length_x_m / 2
                            dy = cell_length_y_m / 2
                            x_vals_to_check = [np.min(component_xvals[-1]) - dx,
                                              np.max(component_xvals[-1]) + dx]
                            y_vals_to_check = [np.min(component_yvals[-1]) - dy,
                                              np.max(component_yvals[-1]) + dy]

                            min_x_val, max_x_val, min_y_val, max_y_val = update_min_max_xy(
                                    min_x_val, max_x_val, min_y_val, max_y_val,
                                    x_vals_to_check, y_vals_to_check)

                    elif comp.class_type in FAULT_COMPONENTS:
                        component_cell_areas.append(None)

                        nSegments = iamcommons.get_parameter_val(comp, 'nSegments')
                        component_num_parts.append(nSegments)

                        comp_data = yaml_data[compName]
                        loc_data = comp_data['Segments']['Locations']

                        component_xvals.append([])
                        component_yvals.append([])
                        for locRef in range(len(loc_data['coordx'])):
                            component_xvals[-1].append(loc_data['coordx'][locRef])
                            component_yvals[-1].append(loc_data['coordy'][locRef])

                        length = iamcommons.get_parameter_val(comp, 'length')
                        component_length.append(length)

                        strike = iamcommons.get_parameter_val(comp, 'strike')
                        component_strike.append(strike)

                        dip = iamcommons.get_parameter_val(comp, 'dip')
                        component_dip.append(dip)

                        # Adjust the min and ax x and y to account for the fault length
                        input_distance = length / 2
                        dx, dy = get_loc_from_strike(strike, input_distance)

                        # Here, dx and dy can be positive or negative (depends on strike)
                        x_vals_to_check = [np.min(component_xvals[-1]) + dx,
                                           np.min(component_xvals[-1]) - dx,
                                           np.max(component_xvals[-1]) + dx,
                                           np.max(component_xvals[-1]) - dx]
                        y_vals_to_check = [np.min(component_yvals[-1]) + dy,
                                           np.min(component_yvals[-1]) - dy,
                                           np.max(component_yvals[-1]) + dy,
                                           np.max(component_yvals[-1]) - dy]

                        min_x_val, max_x_val, min_y_val, max_y_val = update_min_max_xy(
                                min_x_val, max_x_val, min_y_val, max_y_val,
                                x_vals_to_check, y_vals_to_check)

                        input_distance = length / 2
                        dx, dy = get_loc_from_strike(strike, input_distance)

                    x_vals_to_check = [np.min(component_xvals[-1]),
                                       np.max(component_xvals[-1])]
                    y_vals_to_check = [np.min(component_yvals[-1]),
                                       np.max(component_yvals[-1])]

                    min_x_val, max_x_val, min_y_val, max_y_val = update_min_max_xy(
                            min_x_val, max_x_val, min_y_val, max_y_val,
                            x_vals_to_check, y_vals_to_check)

                else:
                    # Leaving this as a placeholder for other component types
                    pass

        if comp.class_type in RESERVOIR_COMPONENTS:
            # Get the injection sites
            if comp.class_type != 'LookupTableReservoir':
                res_comp_injX.append(comp.injX)
                res_comp_injY.append(comp.injY)

            else:
                res_comp_injX.append(None)
                res_comp_injY.append(None)

            GriddedMetric_yaml_input_dict = get_gridded_plot_yaml_input(
                yaml_data, name)

            if GriddedMetric_yaml_input_dict['plot_injection_sites']:
                if comp.class_type != 'LookupTableReservoir':
                    if comp.injX < min_x_val:
                        min_x_val = comp.injX

                    if comp.injX > max_x_val:
                        max_x_val = comp.injX

                    if comp.injY < min_y_val:
                        min_y_val = comp.injY

                    if comp.injY > max_y_val:
                        max_y_val = comp.injY
                else:
                    InjectionCoordx = GriddedMetric_yaml_input_dict['InjectionCoordx']
                    InjectionCoordy = GriddedMetric_yaml_input_dict['InjectionCoordy']

                    if np.min(InjectionCoordx) < min_x_val:
                        min_x_val = np.min(InjectionCoordx)

                    if np.max(InjectionCoordx) > max_x_val:
                        max_x_val = np.max(InjectionCoordx)

                    if np.min(InjectionCoordy) < min_y_val:
                        min_y_val = np.min(InjectionCoordy)

                    if np.max(InjectionCoordy) > max_y_val:
                        max_y_val = np.max(InjectionCoordy)

    x_range = [min_x_val, max_x_val]
    y_range = [min_y_val, max_y_val]

    if len(components_to_use) == 0:
        err_msg = ''.join([
            'In this version of NRAP-Open-IAM, the GriddedMetric ',
            'plot type can only handle output from the following ',
            'component types: ', str(GRIDDED_METRIC_COMPONENTS), '. ',
            'The ComponentNameList entry provided for the GriddedMetric ',
            'plot ', name, ' did not match any compatible components. ',
            'Check your input for the plot ', name, '.'])
        raise KeyError(err_msg)

    return components_to_use, component_types, component_cell_areas, \
        component_num_parts, component_xvals, component_yvals, \
            component_length, component_strike, component_dip, \
                res_comp_injX, res_comp_injY, x_range, y_range


def update_min_max_xy(min_x_val, max_x_val, min_y_val, max_y_val,
                      x_vals_to_check, y_vals_to_check):
    """
    Goes through lists of x and y values and updates the minimum and maximum
    x and y values.
    """
    min_x = min(min_x_val, min(x_vals_to_check))
    max_x = max(max_x_val, max(x_vals_to_check))

    min_y = min(min_y_val, min(y_vals_to_check))
    max_y = max(max_y_val, max(y_vals_to_check))

    return min_x, max_x, min_y, max_y


def get_gridded_plot_yaml_input(yaml_data, name, plot_type='GriddedMetric'):
    """
    Function that reads the gridded metric plot input provided in the .yaml file
    and returns a dictionary containing the input.

    :param yaml_data: Dictionary of input values from the entire .yaml file
    :type yaml_data: dict

    :param name: Name of the plot provided in the Plots section of the
        .yaml file
    :type name: str

    :returns: gridplot_yaml_input_dict
    """
    comp_name_list = yaml_data['Plots'][name][plot_type]['ComponentNameList']

    metric_name = yaml_data['Plots'][name][plot_type]['MetricName']

    defaultTimeList = 'All'
    TimeList = defaultTimeList

    defaultFigureDPI = 100
    FigureDPI = defaultFigureDPI

    defaultRealization = None
    Realization = defaultRealization

    InjectionCoord_debug_msg = ''.join([
        'InjectionCoord{} was provided for the GriddedMetric plot ', name,
        ', but not InjectionCoord{}. Check your input. Injection sites will not ',
        'be displayed.'])

    axis_lims_debug_msg = ''.join([
        'The {}-axis limits provided for the GriddedMetric plot ', name,
        ' (SpecifyXandYLims: {}Lims) are not of length 2 and will not be used. ',
        'Check your inputs in the .yaml file.'])

    bool_option_debug_msg = ''.join([
        'The input provided for {} in the GriddedMetric plot ', name,
        ' was not of type boolean. {} will be set to the default value of {}.'])

    defaultPlotInjectionSites = False
    plot_injection_sites = defaultPlotInjectionSites

    EnforceXandYLims = False
    xLims = None
    yLims = None

    default_equal_axes = True
    equal_axes = default_equal_axes

    default_plot_val_over_area = True
    plot_val_over_area = default_plot_val_over_area

    cell_length_x = None
    cell_length_y = None

    defaultSaveCSVFiles = True
    saveCSVFiles = defaultSaveCSVFiles

    if 'TimeList' in yaml_data['Plots'][name][plot_type]:
        TimeList = yaml_data['Plots'][name][plot_type]['TimeList']

    if 'FigureDPI' in yaml_data['Plots'][name][plot_type]:
        FigureDPI = yaml_data['Plots'][name][plot_type]['FigureDPI']

        if not isinstance(FigureDPI, (int, float)):
            debug_msg = ''.join([
                'The FigureDPI provided for the GriddedMetric plot ', name,
                ' was not of type int or float. The dpi (dots-per-inch) will be ',
                'set to the default value of {}']).format(defaultFigureDPI)
            logging.debug(debug_msg)
            FigureDPI = defaultFigureDPI

    if 'PlotInjectionSites' in yaml_data['Plots'][name][plot_type]:
        plot_injection_sites = yaml_data['Plots'][name][plot_type]['PlotInjectionSites']
        if not isinstance(plot_injection_sites, bool):
            debug_msg = bool_option_debug_msg.format(
                'PlotInjectionSites', 'PlotInjectionSites', str(defaultPlotInjectionSites))
            logging.debug(debug_msg)
            plot_injection_sites = defaultPlotInjectionSites

    InjectionCoordx = None
    InjectionCoordy = None
    if 'InjectionCoordx' in yaml_data['Plots'][name][plot_type]:
        InjectionCoordx = yaml_data['Plots'][name][plot_type]['InjectionCoordx']

    if 'InjectionCoordy' in yaml_data['Plots'][name][plot_type]:
        InjectionCoordy = yaml_data['Plots'][name][plot_type]['InjectionCoordy']
        if InjectionCoordx is None:
            debug_msg = InjectionCoord_debug_msg.format('y', 'x')
            logging.debug(debug_msg)
            plot_injection_sites = defaultPlotInjectionSites

        elif isinstance(InjectionCoordx, list) and isinstance(InjectionCoordy, list):

            if len(InjectionCoordx) != len(InjectionCoordy):
                debug_msg = ''.join(['The InjectionCoordy provided was not of ',
                                     'the same length as the InjectionCoordx ',
                                     'provided. Check your input. Injection ',
                                     'sites will not be displayed.'])
                logging.debug(debug_msg)
                plot_injection_sites = False

    if InjectionCoordx is not None and InjectionCoordy is None:
        debug_msg = InjectionCoord_debug_msg.format('x', 'y')
        logging.debug(debug_msg)
        plot_injection_sites = defaultPlotInjectionSites

    if plot_injection_sites and InjectionCoordx is not None \
            and InjectionCoordy is not None:
        try:
            InjectionCoordx = float(InjectionCoordx)
        except TypeError:
            InjectionCoordx = np.array(list(InjectionCoordx), dtype=float)

        try:
            InjectionCoordy = float(InjectionCoordy)
        except TypeError:
            InjectionCoordy = np.array(list(InjectionCoordy), dtype=float)

    elif not plot_injection_sites:
        InjectionCoordx = None
        InjectionCoordy = None

    if 'SpecifyXandYLims' in yaml_data['Plots'][name][plot_type]:
        EnforceXandYLims = True
        xLims = yaml_data['Plots'][name][plot_type]['SpecifyXandYLims']['xLims']
        yLims = yaml_data['Plots'][name][plot_type]['SpecifyXandYLims']['yLims']

        if len(xLims) != 2:
            debug_msg = axis_lims_debug_msg.format('x', 'x')
            logging.debug(debug_msg)
            EnforceXandYLims = False

        if len(yLims) != 2:
            debug_msg = axis_lims_debug_msg.format('y', 'y')
            logging.debug(debug_msg)
            EnforceXandYLims = False

    if 'Realization' in yaml_data['Plots'][name][plot_type]:
        Realization = int(yaml_data['Plots'][name][plot_type]['Realization'])

        # The standard practice in the control file interface is to refer to a
        # realization with the python indexing rules (0 for the 1st entry). The
        # files saved by the SealHorizon and FaultFlow components do not use the
        # python indexing rules when saving output .npz files, however, so the
        # Realization value is increased by one.
        Realization += 1

    if 'EqualAxes' in yaml_data['Plots'][name][plot_type]:
        equal_axes = yaml_data['Plots'][name][plot_type]['EqualAxes']

        if not isinstance(equal_axes, bool):
            debug_msg = bool_option_debug_msg.format(
                'EqualAxes', 'EqualAxes', str(default_equal_axes))
            logging.debug(debug_msg)
            equal_axes = default_equal_axes

    if 'PlotOverAreas' in yaml_data['Plots'][name][plot_type]:
        plot_val_over_area = yaml_data['Plots'][name][plot_type]['PlotOverAreas']

        if not isinstance(equal_axes, bool):
            debug_msg = bool_option_debug_msg.format(
                'PlotOverAreas', 'PlotOverAreas', str(default_plot_val_over_area))
            logging.debug(debug_msg)
            plot_val_over_area = default_plot_val_over_area

    if 'CellLengthX' in yaml_data['Plots'][name][plot_type]:
        cell_length_x = yaml_data['Plots'][name][plot_type]['CellLengthX']

    if 'CellLengthY' in yaml_data['Plots'][name][plot_type]:
        cell_length_x = yaml_data['Plots'][name][plot_type]['CellLengthY']

    if cell_length_x is None and cell_length_y is not None:
        cell_length_y = None

    elif cell_length_y is None and cell_length_x is not None:
        cell_length_x = None

    if 'SaveCSVFiles' in yaml_data['Plots'][name][plot_type]:
        saveCSVFiles = yaml_data['Plots'][name][plot_type]['SaveCSVFiles']

        if not isinstance(saveCSVFiles, bool):
            debug_msg = bool_option_debug_msg.format(
                'SaveCSVFiles', 'SaveCSVFiles', str(defaultSaveCSVFiles))
            logging.debug(debug_msg)
            saveCSVFiles = defaultSaveCSVFiles

    GriddedMetric_yaml_input_dict = {
        'TimeList': TimeList,
        'FigureDPI': FigureDPI,
        'Realization': Realization,
        'comp_name_list': comp_name_list,
        'metric_name': metric_name,
        'plot_injection_sites': plot_injection_sites,
        'InjectionCoordx': InjectionCoordx,
        'InjectionCoordy': InjectionCoordy,
        'EnforceXandYLims': EnforceXandYLims,
        'xLims': xLims,
        'yLims': yLims,
        'EqualAxes': equal_axes,
        'PlotOverAreas': plot_val_over_area,
        'CellLengthX': cell_length_x,
        'CellLengthY': cell_length_y,
        'SaveCSVFiles': saveCSVFiles
        }

    return GriddedMetric_yaml_input_dict


def get_t_indices(time_list, time_array):
    """
    Returns the time index corresponding to the time_array value closest to
    the times in time_list. Note that time_array is in days, while time_list is
    in years.
    """
    time_index_list = []

    corr_t_index = np.arange(0, len(time_array))

    if not isinstance(time_list, list):
        time_list = [time_list]

    for time in time_list:
        abs_diff = np.zeros(len(time_array))

        for t_ref, t_val in enumerate(time_array):
            abs_diff[t_ref] = np.abs(time - (t_val / 365.25))

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


def figure_setup(sim_index, fig_num, levels, comp_name, metric_name,
                 x_axis_limits, y_axis_limits, genfontsize=10,
                 axislabelfontsize=14, titlefontsize=14, boldlabels=True,
                 figsize=(10, 8), colormap='cividis', equal_axes=True, min_value=0):
    """
    Sets up the figure.
    """
    # Figures
    font = RC_FONT
    font['size'] = genfontsize
    plt.rc('font', **font)

    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${}\times10^{{{}}}$'.format(a, b)

    if boldlabels:
        selected_labelfontweight = 'bold'
    else:
        selected_labelfontweight = 'normal'

    fig = plt.figure(fig_num, figsize=figsize)
    ax = plt.gca()

    ax.set_facecolor(BACKGROUND_COLOR)

    cmap = plt.cm.get_cmap(colormap)
    cbar = fig.colorbar(cm.ScalarMappable(cmap = cmap), ax = ax,
                        values = levels, format=ticker.FuncFormatter(fmt))
    cbar.set_label(
        Y_LABEL_DICT.get(metric_name, metric_name),
        rotation = 90, fontsize = axislabelfontsize,
        fontweight = selected_labelfontweight)
    tick_locator = ticker.MaxNLocator(nbins = 5)
    cbar.locator = tick_locator
    cbar.update_ticks()

    # Remvove min_value from the colorbar. If the min_value is 0, then 0 is
    # excluded from the color range and it shouldn't be included in the labels.
    cbar_ticks = cbar.ax.get_yticks()
    cbar_ticks = cbar_ticks[cbar_ticks > min_value].tolist()
    cbar.ax.set_yticks(cbar_ticks)
    cbar.ax.set_ylim([np.min(levels), np.max(levels)])

    plt.xlabel('Easting (km)', fontsize = axislabelfontsize,
               fontweight = selected_labelfontweight)
    plt.ylabel('Northing (km)', fontsize = axislabelfontsize,
               fontweight = selected_labelfontweight)

    if equal_axes:
        x_range = x_axis_limits[1] - x_axis_limits[0]
        y_range = y_axis_limits[1] - y_axis_limits[0]

        if x_range > y_range:
            expand_ax_lim = (x_range - y_range) / 2

            y_axis_limits[1] += expand_ax_lim
            y_axis_limits[0] -= expand_ax_lim

        elif y_range > x_range:
            expand_ax_lim = (y_range - x_range) / 2

            x_axis_limits[1] += expand_ax_lim
            x_axis_limits[0] -= expand_ax_lim

    plt.xlim(x_axis_limits)
    plt.ylim(y_axis_limits)

    plt.subplots_adjust(left=0.125, bottom=0.125, right=0.875,
                        top=0.875, wspace=0.1, hspace=0.1)

    return fig, ax


def plot_comp_inj_sites(ax, comp_name, comp_type, cell_xvals, cell_yvals,
                        plot_injection_sites, InjectionCoordx, InjectionCoordy,
                        genfontsize=10, include_labels=True,
                        cellCenterMarkerSize=6, plot_faults=True,
                        component_length=None, component_strike=None,
                        component_dip=None):
    """
    Function that plots the wells and injection sites used in the simulation.
    """

    if plot_injection_sites:
        if isinstance(InjectionCoordx, float):
            if include_labels:
                plt.plot(InjectionCoordx / 1000, InjectionCoordy / 1000,
                         marker='s', color='k', linestyle='none',
                         markeredgewidth=2, markersize=8,
                         markerfacecolor='none', zorder=1e6,
                         label='Injection Site')
            else:
                plt.plot(InjectionCoordx / 1000, InjectionCoordy / 1000,
                         marker='s', color='k', linestyle='none',
                         markeredgewidth=2, markersize=8,
                         markerfacecolor='none', zorder=1e6)
        else:
            for injRef, (icoordX, icoordY) in enumerate(
                    zip(InjectionCoordx, InjectionCoordy)):
                if injRef == 0 and include_labels:
                    plt.plot(icoordX / 1000, icoordY / 1000,
                             marker='s', color='k', linestyle='none',
                             markeredgewidth=2, markersize=8,
                             markerfacecolor='none', zorder=1e6,
                             label='Injection Site')
                else:
                    plt.plot(icoordX / 1000, icoordY / 1000,
                             marker='s', color='k', linestyle='none',
                             markeredgewidth=2, markersize=8,
                             markerfacecolor='none', zorder=1e6)

    if include_labels:
        if '_' in comp_name:
            comp_name_edit = comp_name[0:comp_name.index('_')]
        else:
            comp_name_edit = comp_name

        loc_ind = time_series.is_location_in_name(comp_name)

        if loc_ind == -1:
            label = comp_name_edit
        else:
            label = comp_name_edit + ' at Location {:.0f},'.format(loc_ind)

        if comp_type in AREA_COMPONENTS:
            label += ' Cell Centers'
        elif comp_type in FAULT_COMPONENTS:
            label += ' Segment Centers'

        plt.plot(np.array(cell_xvals) / 1000, np.array(cell_yvals) / 1000,
                 linestyle='none', marker='o', color='k', markeredgewidth=1.5,
                 markersize=cellCenterMarkerSize, markerfacecolor='none',
                 zorder=1e6, label=label)
    else:
        plt.plot(np.array(cell_xvals) / 1000, np.array(cell_yvals) / 1000,
                 linestyle='none', marker='o', color='k', markeredgewidth=1.5,
                 markersize=cellCenterMarkerSize, markerfacecolor='none', zorder=1e6)

    if comp_type in FAULT_COMPONENTS and plot_faults:
        if component_length is not None and component_strike is not None \
                and component_dip is not None:
            for locRef, (center_x, center_y) in enumerate(zip(cell_xvals, cell_yvals)):

                # Using half the fault length but extending it in two directions,
                # making it appear with the full length.
                input_distance = component_length / 2
                dx, dy = get_loc_from_strike(component_strike, input_distance)

                # This is a point in the direction of strike
                x_along_strike = center_x + dx
                y_along_strike = center_y + dy

                # This is a point in the direction opposite of the reported
                # strike direction.
                x_opp_strike = center_x - dx
                y_opp_strike = center_y - dy

                fault_line_x = [x_along_strike, center_x, x_opp_strike]
                fault_line_y = [y_along_strike, center_y, y_opp_strike]

                if locRef == 0 and include_labels:
                    plt.plot(np.array(fault_line_x) / 1000,
                             np.array(fault_line_y) / 1000,
                             linestyle='-', color='k', linewidth=1.5,
                             label='Fault', zorder=1)
                else:
                    plt.plot(np.array(fault_line_x) / 1000,
                             np.array(fault_line_y) / 1000,
                             linestyle='-', color='k', linewidth=1.5, zorder=1)

    ax.legend(fancybox=False, fontsize=genfontsize,
              edgecolor=[0, 0, 0], loc='upper left',
              framealpha=0.5).set_zorder(2e6)


def get_loc_from_strike(strike, input_distance):
    """
    Function returns the x and y distances (dx, dy) between an input location
    and a location at a certain distance from that location in the direction of
    strike. Note that strike is taken as degrees clockwise from north in a
    map-view image. For example, 90 degrees would be east, 180 degrees would be
    south, and 270 degrees would be west.
    """
    if 0 <= strike <= 90:
        dx = input_distance * np.sin(np.radians(strike))
        dy = input_distance * np.cos(np.radians(strike))

    elif 90 < strike <= 180:
        dx = input_distance * np.cos(np.radians(strike - 90))
        dy = -input_distance * np.sin(np.radians(strike - 90))

    elif 180 < strike <= 270:
        dx = -input_distance * np.sin(np.radians(strike - 180))
        dy = -input_distance * np.cos(np.radians(strike - 180))

    elif 270 < strike <= 360:
        dx = -input_distance * np.cos(np.radians(strike - 270))
        dy = input_distance * np.sin(np.radians(strike - 270))

    return dx, dy


def save_results_to_csv(output_dir, comp_name, analysis, sim_index,
                        metric_name, x_vals_compiled, y_vals_compiled, t_vals_compiled,
                        metric_vals_compiled, all_x_vals, all_y_vals, time_array):
    """
    Function that saves the gridded radial metric data to a .csv file.
    """
    if not os.path.exists(os.path.join(output_dir, 'csv_files')):
        os.mkdir(os.path.join(output_dir, 'csv_files'))

    if analysis == 'forward':
        filename = os.path.join(output_dir, 'csv_files',
                                '{}_{}.csv'.format(
                                    metric_name, comp_name))
    elif analysis in ['lhs', 'parstudy']:
        filename = os.path.join(output_dir, 'csv_files',
                                '{}_{}_Sim{}.csv'.format(
                                    metric_name, comp_name, sim_index))

    num_rows = len(metric_vals_compiled)

    # This variable is None if the simulation doesn't have any values exceeding
    # MinValue
    if num_rows == 0:
        results_formatted = np.empty((2, 5), dtype=list)

        results_formatted[0, 0] = 'x (km)'
        results_formatted[0, 1] = 'y (km)'
        results_formatted[0, 2] = 't (years)'

        results_formatted[0, 3] = CSV_FILE_COLUMNS.get(metric_name, metric_name)

        results_formatted[1, 0] = 'All x values evaluated from {:.2e} km to {:.2e} km'.format(
            np.min(all_x_vals) / 1000, np.max(all_x_vals) / 1000)

        results_formatted[1, 1] = 'All y values evaluated from {:.2e} km to {:.2e} km'.format(
            np.min(all_y_vals) / 1000, np.max(all_y_vals) / 1000)

        results_formatted[1, 2] = 'All times evaluated from {} years to {} years'.format(
            np.min(time_array) / 365.25, np.max(time_array) / 365.25)

        results_formatted[1, 3] = 'No values'

        # Save the ouput for the simulation
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            for row_ref in range(2):
                writer.writerow(results_formatted[row_ref, :])
        f.close()

    else:
        results_formatted = np.empty(((num_rows + 1), 5), dtype=list)

        results_formatted[0, 0] = 'x (km)'
        results_formatted[0, 1] = 'y (km)'
        results_formatted[0, 2] = 't (years)'

        # Get the column name depending on the output
        # If name does not match any known expected observations use it as a column label
        results_formatted[0, 3] = CSV_FILE_COLUMNS.get(metric_name, metric_name)

        for row_ref in range(0, len(metric_vals_compiled)):
            results_formatted[row_ref + 1, 0] = str(x_vals_compiled[row_ref] / 1000)
            results_formatted[row_ref + 1, 1] = str(y_vals_compiled[row_ref] / 1000)
            results_formatted[row_ref + 1, 2] = str(t_vals_compiled[row_ref])
            results_formatted[row_ref + 1, 3] = str(metric_vals_compiled[row_ref])

        # Save the ouput for the simulation
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            for row_ref in range(num_rows + 1):
                writer.writerow(results_formatted[row_ref, :])
        f.close()
