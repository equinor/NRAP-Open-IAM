
import os
import logging
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import ticker

from openiam.visualize import time_series
import openiam.cfi.strata as iam_strata
from .label_setup import Y_LABEL_DICT


RC_FONT = {'family': 'Arial', 'weight': 'normal', 'size': None}

GRIDDED_METRIC_COMPONENTS = ['GenericAquifer']

RESERVOIR_COMPONENTS = ['LookupTableReservoir',
                        'SimpleReservoir',
                        'AnalyticalReservoir']

BACKGROUND_COLOR = [0.67, 0.67, 0.67]

COLORMAP_OPTIONS = {'Dissolved_CO2_mass_fraction': 'plasma',
                    'Dissolved_salt_mass_fraction': 'viridis'}

DEFAULT_FIGURE_DPI = 100
DEFAULT_SAVE_CSV_FILES = True
DEFAULT_REALIZATION = 0
DEFAULT_DEGREE_INTERVAL = 15
ACCEPTABLE_DEGREE_INTERVALS = [1, 5, 10, 15, 30, 45]
DEFAULT_PLOT_INJECTION_SITES = False
DEFAULT_MIN_VALUE = {'Dissolved_CO2_mass_fraction': 0.01,
                     'Dissolved_salt_mass_fraction': 0.002}
DEFAULT_TIME_LIST = 'All'
DEFAULT_MIN_XY = 9.0e99
DEFAULT_MAX_XY = -9.0e99
DEFAULT_MAX_METRIC = -9.0e99

CSV_FILE_COLUMNS = {
    'Dissolved_CO2_mass_fraction':
        'Dissolved CO2 mass fraction [-] exceeding MinValue of {:.2e}',
    'Dissolved_salt_mass_fraction':
        'Dissolved salt mass fraction [-] exceeding MinValue of {:.2e}'}

ZCOORD_FILE_NAME = '{comp_name}_z_coordinate_sim_{sim_index}_time_{t_index}.npz'
RCOORD_FILE_NAME = '{comp_name}_r_coordinate_sim_{sim_index}_time_{t_index}.npz'
METRIC_FILE_NAME = '{comp_name}_{metric_name}_sim_{sim_index}_time_{t_index}.npz'

DETERM_RUN_TITLE = 'Gridded Results for z = {} m and t = {} years,\n'
STOCH_RUN_TITLE = 'Simulation {}, Gridded Results for z = {} m and t = {} years,\n'
TITLE_END1 = 'Range: {} to {}\n(Gray: Results {})'
TITLE_END2 = 'Range: {} to {}\n(Gray: Results of 0)'


def gridded_radial_metric_plot(yaml_data, sm, output_dir,
                               name='RadialObs_Figure1', analysis='forward',
                               savefig=None, figsize=(10, 8), genfontsize=10,
                               axislabelfontsize=14, titlefontsize=14, boldlabels=True,
                               cmap='viridis', figure_dpi=100, min_value=0):
    """
    Produce gridded radial metric plot.
    """
    os.chdir(output_dir)

    grid_metric_yaml_input_dict = get_gridded_plot_yaml_input(
        yaml_data, name, plot_type='GriddedRadialMetric')

    if analysis =='forward':
        sim_index = 0

    elif analysis in ['lhs', 'parstudy']:
        sim_index = grid_metric_yaml_input_dict['Realization']

        if sim_index is None:
            # For LHS or parstudy simulations, the gridded radial data from the
            # GenericAquifer component are saved with indices of 1 and higher.
            # An index of 0 does not work.
            sim_index = 1

    components_name_list = grid_metric_yaml_input_dict['comp_name_list']

    metric_name = grid_metric_yaml_input_dict['metric_name']

    cmap = COLORMAP_OPTIONS.get(metric_name, 'viridis')

    min_value = grid_metric_yaml_input_dict['MinValue']

    z_list = grid_metric_yaml_input_dict['ZList']

    time_list = grid_metric_yaml_input_dict['TimeList']

    figure_dpi = grid_metric_yaml_input_dict['FigureDPI']

    plot_injection_sites = grid_metric_yaml_input_dict['plot_injection_sites']

    injection_coordx = grid_metric_yaml_input_dict['InjectionCoordx']
    injection_coordy = grid_metric_yaml_input_dict['InjectionCoordy']

    enforce_x_and_y_lims = grid_metric_yaml_input_dict['EnforceXandYLims']
    xLims = grid_metric_yaml_input_dict['xLims']
    yLims = grid_metric_yaml_input_dict['yLims']

    degree_interval = grid_metric_yaml_input_dict['DegreeInterval']

    equal_axes = grid_metric_yaml_input_dict['EqualAxes']

    saveCSVFiles = grid_metric_yaml_input_dict['saveCSVFiles']

    # Degrees clockwise from north in map view with north pointing up
    degrees = np.arange(0, 360 + degree_interval, degree_interval)

    time_array = sm.time_array

    num_times = len(time_array)

    if time_list == 'All':
        time_index_list = range(num_times)
    else:
        time_index_list = get_t_indices(time_list, time_array)

    components_to_use, components_types, aq_component_indices, \
        component_xvals, component_yvals, \
            res_comp_injX, res_comp_injY, x_range, y_range = \
                get_comps_and_limits(components_name_list, yaml_data, name, sm)

    # Check if string inputs were provided for the depths (e.g., 'aquifer2Depth').
    # Just in case spatially variable stratigraphy is used and unit depths vary
    # over space, a z_list is created for each component.
    z_list_orig = z_list.copy()
    comp_z_lists = get_comp_z_lists(sm, z_list_orig, name, components_to_use)

    if plot_injection_sites and injection_coordx is None:
        injection_coordx = res_comp_injX
        injection_coordy = res_comp_injY

    component_min_x, component_max_x, component_min_y, component_max_y, \
        component_max_metric, min_value_updated = \
            get_max_vals_over_time_radial(
                components_to_use, component_xvals, component_yvals, degrees,
                comp_z_lists, num_times, metric_name, sim_index,
                min_value=min_value)

    # Check if the component types are all GenericAquifer. Currently, GenericAquifer
    # is the only component type that provides a gridded radial output. This code
    # is being set up to more easily incorporate other component types.
    GenAq_comp_check = True
    for comp in components_to_use:
        if comp.class_type != 'GenericAquifer':
            GenAq_comp_check = False

    if GenAq_comp_check:
        if analysis == 'forward':
            file_name = '{}_{}_z{}_t{}{}'
            file_name_no_results = '{}_{}_z{}_All_t{}'

        elif analysis in ['lhs', 'parstudy']:
            file_name = '{}_{}_Sim{}_z{}_t{}{}'
            file_name_no_results = '{}_{}_Sim{}_z{}_All_t{}'

    # Not using this value right now, but leaving it for potential updates
    # min_value = min_value_updated

    # Check the file name for an extension like .tiff or .eps
    if '.' in name:
        name_extension = name[name.index('.'):None]
    else:
        name_extension = '.png'

    for comp_ref, comp in enumerate(components_to_use):
        z_list = comp_z_lists[comp_ref]

        for z_ref, z_val in enumerate(z_list):

            # Get the x and y location for the component
            xloc = component_xvals[comp_ref]
            yloc = component_yvals[comp_ref]

            # If the component is an aquifer component, get the aquifer index
            # so it can be displayed on figures.
            if GenAq_comp_check:
                aq_comp_index = aq_component_indices[comp_ref]
            else:
                aq_comp_index = None

            comp_name = comp.name

            # Initialize flag variable indicating whether all the results are
            # below the user defined minimum value
            no_results_check = False
            if np.max(component_max_metric[comp_ref, z_ref]) <= min_value:
                no_results_check = True
                interval = (1 - min_value) / 100
                levels = np.arange(min_value, 1 + interval, interval)

            else:
                interval = (component_max_metric[comp_ref, z_ref] - min_value) / 100
                levels = np.arange(min_value,
                                   component_max_metric[comp_ref, z_ref] + interval,
                                   interval)

            if enforce_x_and_y_lims:
                x_axis_limits = [xLims[0] / 1000, xLims[1] / 1000]
                y_axis_limits = [yLims[0] / 1000, yLims[1] / 1000]

            elif no_results_check:
                min_x = xloc
                max_x = xloc

                min_y = yloc
                max_y = yloc

                if plot_injection_sites:
                    if np.min(injection_coordx) < xloc:
                        min_x = np.min(injection_coordx)

                    if np.max(injection_coordx) > xloc:
                        max_x = np.max(injection_coordx)

                    if np.min(injection_coordy) < yloc:
                        min_y = np.min(injection_coordy)

                    if np.max(injection_coordy) > yloc:
                        max_y = np.max(injection_coordy)

                # These are used to adjust the x and y axis limits
                x_buffer = (max_x - min_x) / 20
                y_buffer = (max_y - min_y) / 20

                if x_buffer == 0:
                    x_buffer = 10

                if y_buffer == 0:
                    y_buffer = 10

                x_axis_limits = [(min_x - x_buffer) / 1000,
                                  (max_x + x_buffer) / 1000]
                y_axis_limits = [(min_y - y_buffer) / 1000,
                                  (max_y + y_buffer) / 1000]
            else:
                if plot_injection_sites:
                    if np.min(injection_coordx) < component_min_x[comp_ref, z_ref]:
                        component_min_x[comp_ref, z_ref] = np.min(injection_coordx)

                    if np.max(injection_coordx) > component_max_x[comp_ref, z_ref]:
                        component_max_x[comp_ref, z_ref] = np.max(injection_coordx)

                    if np.min(injection_coordy) < component_min_y[comp_ref, z_ref]:
                        component_min_y[comp_ref, z_ref] = np.min(injection_coordy)

                    if np.max(injection_coordy) > component_max_y[comp_ref, z_ref]:
                        component_max_y[comp_ref, z_ref] = np.max(injection_coordy)

                # These are used to adjust the x and y axis limits
                x_buffer = (component_max_x[comp_ref, z_ref]
                            - component_min_x[comp_ref, z_ref]) / 20
                y_buffer = (component_max_x[comp_ref, z_ref]
                            - component_min_x[comp_ref, z_ref]) / 20

                x_axis_limits = [(component_min_x[comp_ref, z_ref] - x_buffer) / 1000,
                                  (component_max_x[comp_ref, z_ref] + x_buffer) / 1000]
                y_axis_limits = [(component_min_y[comp_ref, z_ref] - y_buffer) / 1000,
                                  (component_max_y[comp_ref, z_ref] + y_buffer) / 1000]

            use_sci_formatting = False
            # If the distances are very large, use scientific notation to avoid
            # the ugly formatting used by default.
            if np.max(x_axis_limits) >= 1000 or np.max(y_axis_limits) >= 1000:
                use_sci_formatting = True

            if no_results_check:
                if GenAq_comp_check:
                    z_data = np.load(ZCOORD_FILE_NAME.format(
                        comp_name=comp_name, sim_index=sim_index, t_index=0))
                    z = z_data['data']
                    z_data.close()
                else:
                    # Leaving this as a placeholder for obtaining the z output
                    # from future components that produce gridded radial data.
                    pass

                z_index = get_z_index(z_val, z)

                z_value = z[0, z_index]
                z_value = int(z_value)

                fig, ax = figure_setup(
                    sim_index, 1, levels, comp_name, z_value, metric_name,
                    x_axis_limits, y_axis_limits, genfontsize=genfontsize,
                    axislabelfontsize=axislabelfontsize,
                    titlefontsize=titlefontsize, boldlabels=boldlabels,
                    figsize=figsize, cmap=cmap, equal_axes=equal_axes)

                if analysis == 'forward':
                    title_str = 'Gridded Results for z = {}'.format(
                        z_value) + ' m,\nAll Model Times'
                elif analysis in ['lhs', 'parstudy']:
                    title_str = 'Simulation {}, Gridded Results for z = {}'.format(
                        sim_index, z_value) + ' m,\nAll Model Times'

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
                    ax, comp_name, xloc, yloc, plot_injection_sites,
                    injection_coordx, injection_coordy, aq_comp_index,
                    genfontsize=genfontsize, include_labels=True)

                plt.title(title_str, fontsize=titlefontsize, fontweight='bold')

                if use_sci_formatting:
                    ax.ticklabel_format(
                        style='sci', axis='both', useOffset=False,
                        scilimits=(0, 0), useMathText=True)

                if analysis == 'forward':
                    if savefig:
                        plt.savefig(os.sep.join([output_dir, file_name_no_results.format(
                            comp_name, metric_name, z_index, name_extension)]),
                            dpi=figure_dpi)

                        plt.close()

                    else:
                        fig.show()
                elif analysis in ['lhs', 'parstudy']:
                    if savefig:
                        plt.savefig(os.sep.join([output_dir, file_name_no_results.format(
                            comp_name, metric_name, sim_index, z_index,
                            name_extension)]), dpi=figure_dpi)

                        plt.close()

                    else:
                        fig.show()

            else:
                # case for no_results_check == False
                for t_index in time_index_list:
                    if GenAq_comp_check:
                        r_data = np.load(RCOORD_FILE_NAME.format(
                            comp_name=comp_name, sim_index=sim_index, t_index=t_index))
                        r = r_data['data']
                        r_data.close()

                        z_data = np.load(ZCOORD_FILE_NAME.format(
                            comp_name=comp_name, sim_index=sim_index, t_index=t_index))
                        z = z_data['data']
                        z_data.close()
                    else:
                        # Leaving this as a placeholder for obtaining the z and
                        # r output from future components that produce gridded
                        # radial data.
                        pass

                    z_index = get_z_index(z_val, z)

                    z_value = z[0, z_index]
                    z_value = int(z_value)

                    fig, ax = figure_setup(
                        sim_index, t_index, levels, comp_name, z_value,
                        metric_name, x_axis_limits, y_axis_limits,
                        genfontsize=genfontsize, axislabelfontsize=axislabelfontsize,
                        titlefontsize=titlefontsize, boldlabels=boldlabels,
                        figsize=figsize, cmap=cmap)

                    if analysis == 'forward':
                        title_str = DETERM_RUN_TITLE.format(
                            z_value, time_array[t_index] / 365.25)
                    elif analysis in ['lhs', 'parstudy']:
                        title_str = STOCH_RUN_TITLE.format(
                            sim_index, z_value, time_array[t_index] / 365.25)

                    plot_comp_inj_sites(
                        ax, comp_name, xloc, yloc, plot_injection_sites,
                        injection_coordx, injection_coordy, aq_comp_index,
                        genfontsize=genfontsize, include_labels=True)

                    r_shape = r.shape

                    x_vals = np.zeros((len(degrees), r_shape[0]))
                    y_vals = np.zeros((len(degrees), r_shape[0]))
                    x_vals_list = []
                    y_vals_list = []

                    for r_ref in range(r_shape[0]):
                        for deg_ref, degree_val in enumerate(degrees):
                            x_vals[deg_ref, r_ref] = xloc + (r[r_ref, 0] * np.sin(
                                np.radians(degree_val)))
                            x_vals_list.append(x_vals[deg_ref, r_ref])

                            y_vals[deg_ref, r_ref] = yloc + (r[r_ref, 0] * np.cos(
                                np.radians(degree_val)))
                            y_vals_list.append(y_vals[deg_ref, r_ref])

                    x_vals_array = np.array(x_vals_list)
                    y_vals_array = np.array(y_vals_list)

                    if GenAq_comp_check:
                        radial_metric_data = np.load(METRIC_FILE_NAME.format(
                            comp_name=comp_name, metric_name=metric_name,
                            sim_index=sim_index, t_index=t_index))
                        radial_metric = radial_metric_data['data']
                        radial_metric_data.close()
                    else:
                        # Leaving this as a placeholder for obtaining the output for
                        # future components that produce gridded radial data.
                        pass

                    # This is here for the case where no values are > min_value
                    radial_metric_array = np.array([])

                    if np.max(radial_metric) > min_value:

                        radial_metric_vals = np.zeros((len(degrees), r_shape[0]))
                        radial_metric_list = []
                        for r_ref in range(r_shape[0]):
                            for deg_ref in range(len(degrees)):
                                radial_metric_vals[deg_ref, r_ref] = \
                                    radial_metric[r_ref, z_index]
                                radial_metric_list.append(
                                    radial_metric_vals[deg_ref, r_ref])

                        radial_metric_array = np.array(radial_metric_list)

                        plt.tricontourf(x_vals_array / 1000.0, y_vals_array / 1000.0,
                                        radial_metric_array, levels,
                                        locator=ticker.MaxNLocator(
                                            nbins=100, prune='lower'), cmap=cmap)

                        min_val_actual = np.min(radial_metric_array[
                            radial_metric_array > 0])
                        max_val = np.max(radial_metric_array)

                        a, b = '{:.2e}'.format(min_val_actual).split('e')
                        b = int(b)
                        min_val_actual_str = r'${}\times10^{{{}}}$'.format(a, b)

                        a, b = '{:.2e}'.format(max_val).split('e')
                        b = int(b)
                        max_val_str = r'${}\times10^{{{}}}$'.format(a, b)

                        if min_value > min_val_actual and min_value > 0:
                            a, b = '{:.2e}'.format(min_value).split('e')
                            b = int(b)
                            min_val_str = r'$\leq {}\times10^{{{}}}$'.format(a, b)

                            title_str_edit = title_str + TITLE_END1.format(
                                min_val_actual_str, max_val_str, min_val_str)

                        elif min_value > min_val_actual and min_value == 0:
                            title_str_edit = title_str + TITLE_END2.format(
                                min_val_actual_str, max_val_str)
                        else:
                            # This case shouldn't be needed, but it's included
                            # to be safe
                            a, b = '{:.2e}'.format(min_value).split('e')
                            b = int(b)
                            min_val_str = r'$\leq {}\times10^{{{}}}$'.format(a, b)

                            title_str_edit = title_str + TITLE_END1.format(
                                min_val_actual_str, max_val_str, min_val_str)

                    else:
                        if min_value > 0:
                            a, b = '{:.2e}'.format(min_value).split('e')
                            b = int(b)
                            min_val_str = r'$\leq {}\times10^{{{}}}$'.format(a, b)

                            title_str_edit = title_str + '(Gray: Results {})'.format(
                                min_val_str)
                        else:
                            title_str_edit = title_str + '(Gray: Results = {})'.format(
                                min_value)

                    plt.title(title_str_edit, fontsize=titlefontsize,
                              fontweight='bold')

                    if use_sci_formatting:
                        ax.ticklabel_format(
                            style='sci', axis='both', useOffset=False,
                            scilimits=(0, 0), useMathText=True)

                    if saveCSVFiles:
                        time = time_array[t_index] / 365.25

                        save_results_to_csv(
                            output_dir, comp_name, z_value, t_index, time,
                            analysis, sim_index, metric_name, x_vals_array,
                            y_vals_array, radial_metric_array, min_value)

                    if analysis == 'forward':
                        if savefig:
                            plt.savefig(os.sep.join([output_dir, file_name.format(
                                comp_name, metric_name, z_index, t_index,
                                name_extension)]), dpi=figure_dpi)

                            plt.close()

                        else:
                            fig.show()

                    elif analysis in ['lhs', 'parstudy']:
                        if savefig:
                            plt.savefig(os.sep.join([output_dir, file_name.format(
                                comp_name, metric_name, sim_index, z_index,
                                t_index, name_extension)]), dpi=figure_dpi)

                            plt.close()

                        else:
                            fig.show()


def get_max_vals_over_time_radial(components_to_use, component_xvals,
                                  component_yvals, degrees, comp_z_lists, num_times,
                                  metric_name, sim_index, min_value=0):
    """
    Evaluate the distribution of gridded radial observations over time and
    return the maximum value as well as the minimum and maximum x- and y-values
    exceeding min_value over the course of the simulation. This maximum radius
    is used to set the x- and y-axis limits; otherwise, the limits used might
    not allow the data to be seen clearly.
    """
    # If the component_max_metric value is less than min_value, there will be
    # an error when creating the variable levels. To avoid that error, the
    # min_value is updated to account for the component_max_metric values.
    min_value_updated = min_value

    # These store the minimum and maximum x and y values at which the metric
    # exceeds min_value in the simulation.
    num_comp_xvals = len(component_xvals)
    num_comp_yvals = len(component_yvals)
    component_min_x = np.ones(
        (num_comp_xvals, len(comp_z_lists[0]))) * DEFAULT_MIN_XY
    component_max_x = np.ones(
        (num_comp_xvals, len(comp_z_lists[0]))) * DEFAULT_MAX_XY

    component_min_y = np.ones(
        (num_comp_yvals, len(comp_z_lists[0]))) * DEFAULT_MIN_XY
    component_max_y = np.ones(
        (num_comp_yvals, len(comp_z_lists[0]))) * DEFAULT_MAX_XY

    # This stores the maximum output value for each component. These values are
    # used to set the colorbar limits - otherwise, the colorbars for different
    # model times would change, which impedes comparison.
    component_max_metric = np.ones(
        (num_comp_yvals, len(comp_z_lists[0]))) * DEFAULT_MAX_METRIC

    for comp_ref in range(num_comp_xvals):
        z_list = comp_z_lists[comp_ref]

        xloc = component_xvals[comp_ref]
        yloc = component_yvals[comp_ref]

        comp = components_to_use[comp_ref]
        comp_name = comp.name

        for z_ref, z_val in enumerate(z_list):
            for t_index in range(num_times):
                if comp.class_type == 'GenericAquifer':
                    r_data = np.load(RCOORD_FILE_NAME.format(
                        comp_name=comp_name, sim_index=sim_index, t_index=t_index))
                    r = r_data['data']
                    r_data.close()

                    z_data = np.load(ZCOORD_FILE_NAME.format(
                        comp_name=comp_name, sim_index=sim_index, t_index=t_index))
                    z = z_data['data']
                    r_data.close()
                else:
                    # Leaving this as a placeholder for obtaining the z and
                    # r output from future components that produce gridded
                    # radial data.
                    pass

                r_shape = r.shape

                z_index = get_z_index(z_val, z)

                x_vals = np.zeros((len(degrees), r_shape[0]))
                y_vals = np.zeros((len(degrees), r_shape[0]))
                x_vals_list = []
                y_vals_list = []

                for r_ref in range(r_shape[0]):
                    for deg_ref, deg_val in enumerate(degrees):
                        x_vals[deg_ref, r_ref] = xloc + (r[r_ref, 0] * np.sin(
                            np.radians(deg_val)))
                        x_vals_list.append(x_vals[deg_ref, r_ref])

                        y_vals[deg_ref, r_ref] = yloc + (r[r_ref, 0] * np.cos(
                            np.radians(deg_val)))
                        y_vals_list.append(y_vals[deg_ref, r_ref])

                x_vals_array = np.array(x_vals_list)
                y_vals_array = np.array(y_vals_list)

                if comp.class_type == 'GenericAquifer':
                    radial_metric_data = np.load(METRIC_FILE_NAME.format(
                        comp_name=comp_name, metric_name=metric_name,
                        sim_index=sim_index, t_index=t_index))
                    radial_metric = radial_metric_data['data']
                    radial_metric_data.close()
                else:
                    # Leaving this as a placeholder for obtaining the output for
                    # future components that produce gridded radial data.
                    pass

                radial_metric_vals = np.zeros((len(degrees), r_shape[0]))
                radial_metric_list = []
                for r_ref in range(r_shape[0]):
                    for deg_ref in range(len(degrees)):
                        radial_metric_vals[deg_ref, r_ref] = \
                            radial_metric[r_ref, z_index]
                        radial_metric_list.append(
                            radial_metric_vals[deg_ref, r_ref])

                radial_metric_array = np.array(radial_metric_list)

                for loc_ref, (x_val, y_val) in enumerate(zip(x_vals_array, y_vals_array)):
                    if radial_metric_array[loc_ref] > min_value:
                        if component_min_x[comp_ref, z_ref] > x_val:
                            component_min_x[comp_ref, z_ref] = x_val

                        if component_max_x[comp_ref, z_ref] < x_val:
                            component_max_x[comp_ref, z_ref] = x_val

                    if radial_metric_array[loc_ref] > min_value:
                        if component_min_y[comp_ref, z_ref] > y_val:
                            component_min_y[comp_ref, z_ref] = y_val

                        if component_max_y[comp_ref, z_ref] < y_val:
                            component_max_y[comp_ref, z_ref] = y_val

                    if component_max_metric[comp_ref, z_ref] < radial_metric_array[
                            loc_ref] and radial_metric_array[loc_ref] > 0:
                        component_max_metric[comp_ref, z_ref] = \
                            radial_metric_array[loc_ref]

                        if min_value_updated > component_max_metric[comp_ref, z_ref]:
                            min_value_updated = 0.9 * component_max_metric[comp_ref, z_ref]

    # If no results will be shown, adjust these so the graph doesn't have absurd limits
    if np.min(component_xvals) == np.max(component_xvals):
        if np.max(component_xvals) == 0:
            x_adjust = 1
        else:
            x_adjust = np.max(component_xvals) / 10
    else:
        x_adjust = (np.max(component_xvals) - np.min(component_xvals)) / 10

    component_min_x[component_min_x == DEFAULT_MIN_XY] = np.min(
        component_xvals) - x_adjust
    component_max_x[component_max_x == DEFAULT_MAX_XY] = np.max(
        component_xvals) + x_adjust

    if np.min(component_yvals) == np.max(component_yvals):
        if np.max(component_yvals) == 0:
            y_adjust = 1
        else:
            y_adjust = np.max(component_yvals) / 10
    else:
        y_adjust = (np.max(component_yvals) - np.min(component_yvals)) / 10

    component_min_y[component_min_x == DEFAULT_MIN_XY] = np.min(
        component_yvals) - y_adjust
    component_max_y[component_max_x == DEFAULT_MAX_XY] = np.max(
        component_yvals) + y_adjust

    component_max_metric[component_max_metric == DEFAULT_MAX_METRIC] = min_value

    return component_min_x, component_max_x, \
        component_min_y, component_max_y, component_max_metric, min_value_updated


def get_comps_and_limits(components_name_list, yaml_data, name, sm):
    """
    Function that produces lists of components to use and their characteristics
    (e.g., x- and y-values and the index for the corresponding aquifer, if
    applicable). By going through these components, it also obtains and provides
    the minimum x- and y-values used in the simulation(s).
    """

    # These are overwritten within the loop below
    # Start with default values
    min_x_val = DEFAULT_MIN_XY
    max_x_val = DEFAULT_MAX_XY
    min_y_val = DEFAULT_MIN_XY
    max_y_val = DEFAULT_MAX_XY

    components_to_use = []
    components_types = []
    aq_component_indices = []
    component_xvals = []
    component_yvals = []

    res_comp_injX = []
    res_comp_injY = []

    components = list(sm.component_models.values())
    for comp in components:
        if comp.class_type != 'Stratigraphy':
            comp_data = yaml_data[comp.name]
        else:
            continue

        if comp.class_type in GRIDDED_METRIC_COMPONENTS:
            for aq_name in components_name_list:

                if aq_name in comp.name:
                    if 'AquiferName' in comp_data:
                        aq_name_from_yaml = comp_data['AquiferName']
                        # See if the aq. number has one digit or two, then get it
                        try:
                            aq_number_from_yaml = int(aq_name_from_yaml[7:])
                            aq_index = aq_number_from_yaml
                        except ValueError:
                            err_msg = ''.join([
                                'The aquifer component {} has an entry of {} for ',
                                'the AquiferName field. This entry does not ',
                                'match the expected format of aquifer#, where ',
                                '# is an integer between 1 and 29. Check your ',
                                'input in the .yaml file.']).format(
                                    comp.name, comp_data['AquiferName'])
                            raise ValueError(err_msg) from None
                    else:
                        aq_index = None

                    components_to_use.append(comp)
                    components_types.append(comp.class_type)

                    try:
                        # If reservoir data to wellbore components come from
                        # reservoir component get locations from it
                        res_comp = sm.component_models[yaml_data[yaml_data[
                            comp.name]['Connection']]['Connection']]

                        component_xvals.append(res_comp.locX)
                        component_yvals.append(res_comp.locY)
                    except KeyError:
                        # Get location from the wellbore component itself
                        well_comp = yaml_data[comp.name]['Connection']
                        component_xvals.append(yaml_data[well_comp]['locX'])
                        component_yvals.append(yaml_data[well_comp]['locY'])

                    aq_component_indices.append(aq_index)

        if comp.class_type in RESERVOIR_COMPONENTS:
            # Get the injection sites
            if comp.class_type != 'LookupTableReservoir':
                res_comp_injX.append(comp.injX)
                res_comp_injY.append(comp.injY)

            else:
                res_comp_injX.append(None)
                res_comp_injY.append(None)

            grid_metric_yaml_input_dict = get_gridded_plot_yaml_input(yaml_data, name)

            if grid_metric_yaml_input_dict['plot_injection_sites']:
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
                    injection_coordx = grid_metric_yaml_input_dict['InjectionCoordx']
                    injection_coordy = grid_metric_yaml_input_dict['InjectionCoordy']

                    if np.min(injection_coordx) < min_x_val:
                        min_x_val = np.min(injection_coordx)

                    if np.max(injection_coordx) > max_x_val:
                        max_x_val = np.max(injection_coordx)

                    if np.min(injection_coordy) < min_y_val:
                        min_y_val = np.min(injection_coordy)

                    if np.max(injection_coordy) > max_y_val:
                        max_y_val = np.max(injection_coordy)

            # Get the well / aquifer component locations from the res. comp.
            x_vals = comp.locX
            y_vals = comp.locY

            if np.min(x_vals) < min_x_val:
                min_x_val = np.min(x_vals)

            if np.max(x_vals) > max_x_val:
                max_x_val = np.max(x_vals)

            if np.min(y_vals) < min_y_val:
                min_y_val = np.min(y_vals)

            if np.max(y_vals) > max_y_val:
                max_y_val = np.max(y_vals)

    x_range = [min_x_val, max_x_val]
    y_range = [min_y_val, max_y_val]

    if len(components_to_use) == 0:
        err_msg = ''.join([
            'In this version of NRAP-Open-IAM, the GriddedRadialMetric ',
            'plot type can only handle output from the following ',
            'component types: ', GRIDDED_METRIC_COMPONENTS, '. ',
            'The ComponentNameList entry provided for the GriddedRadialMetric ',
            'plot ', name, ' did not correspond with any compatible components. ',
            'Check your input for the plot ', name, '.'])
        raise KeyError(err_msg)

    return components_to_use, components_types, aq_component_indices, \
        component_xvals, component_yvals, \
            res_comp_injX, res_comp_injY, x_range, y_range


def get_gridded_plot_yaml_input(yaml_data, name, plot_type='GriddedRadialMetric'):
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

    time_list = DEFAULT_TIME_LIST

    min_value = DEFAULT_MIN_VALUE.get(metric_name, 0)

    figure_dpi = DEFAULT_FIGURE_DPI

    saveCSVFiles = DEFAULT_SAVE_CSV_FILES

    realization = DEFAULT_REALIZATION

    degree_interval = DEFAULT_DEGREE_INTERVAL

    InjectionCoord_debug_msg = ''.join([
        'InjectionCoord{} is provided for the GriddedRadialMetric plot ', name,
        ', but InjectionCoord{} is not. Check your input. Injection sites will not ',
        'be displayed.'])

    axis_lims_debug_msg = ''.join([
        'The {}-axis limits provided for the GriddedRadialMetric plot ', name,
        ' (SpecifyXandYLims: {}Lims) are not of length 2 and will not be used. ',
        'Check your inputs in the .yaml file.'])

    bool_option_debug_msg = ''.join([
        'The input provided for {} in the GriddedRadialMetric plot ', name,
        ' is not of type boolean. {} will be set to the default value of {}.'])

    plot_injection_sites = DEFAULT_PLOT_INJECTION_SITES

    enforce_x_and_y_lims = False
    xLims = None
    yLims = None

    default_equal_axes = True
    equal_axes = default_equal_axes

    if 'MinValue' in yaml_data['Plots'][name][plot_type]:
        min_value = float(yaml_data['Plots'][name][plot_type]['MinValue'])

        if min_value < 0:
            debug_msg = ''.join([
                'The MinValue provided for the {} plot '.format(plot_type),
                name, ' is < 0. MinValue must be >= 0. ',
                'The default MinValue of {} will be used.'.format(
                    DEFAULT_MIN_VALUE.get(metric_name, 0))])
            logging.debug(debug_msg)
            min_value = DEFAULT_MIN_VALUE.get(metric_name, 0)

    if 'ZList' in yaml_data['Plots'][name][plot_type]:
        z_list = yaml_data['Plots'][name][plot_type]['ZList']

        if not isinstance(z_list, list):
            if isinstance(z_list, (float, int)):
                z_list = [z_list]
            else:
                err_msg = ''.join([
                    'The ZList provided for the GriddedRadialMetric plot ', name,
                    ' is not a list (e.g., [-400, -500] or [aquifer2Depth, shale3Depth]) ',
                    'or a single value (e.g., -500). Since ZList is required for the ',
                    'GriddedRadialMetric plot type, the plot cannot be made. ',
                    'Check your input for the plot ', name, '.'])
                raise KeyError(err_msg)
    else:
        err_msg = ''.join([
            'The ZList entry is required for the GriddedRadialMetric plot type, ',
            'but it is not provided. The plot ', name, ' cannot be made.'])
        raise KeyError(err_msg)

    if 'TimeList' in yaml_data['Plots'][name][plot_type]:
        time_list = yaml_data['Plots'][name][plot_type]['TimeList']

    if 'DegreeInterval' in yaml_data['Plots'][name][plot_type]:
        degree_interval = int(yaml_data['Plots'][name][plot_type]['DegreeInterval'])

        if degree_interval not in ACCEPTABLE_DEGREE_INTERVALS:
            debug_msg = ''.join([
                'The DegreeInterval provided for the {} plot '.format(plot_type),
                name, ' is not one of the acceptable values (1, 5, 10, 15, 30, or 45). ',
                'The default DegreeInterval of 15 degrees will be used.'])
            logging.debug(debug_msg)
            degree_interval = DEFAULT_DEGREE_INTERVAL

    if 'FigureDPI' in yaml_data['Plots'][name][plot_type]:
        figure_dpi = yaml_data['Plots'][name][plot_type]['FigureDPI']

        if not isinstance(figure_dpi, (int, float)):
            debug_msg = ''.join([
                'The FigureDPI provided for the GriddedRadialMetric plot ', name,
                ' is not of integer or float type. The dpi (dots-per-inch) will be ',
                'set to the default value of {}']).format(DEFAULT_FIGURE_DPI)
            logging.debug(debug_msg)
            figure_dpi = DEFAULT_FIGURE_DPI

    if 'PlotInjectionSites' in yaml_data['Plots'][name][plot_type]:
        plot_injection_sites = yaml_data['Plots'][name][plot_type]['PlotInjectionSites']
        if not isinstance(plot_injection_sites, bool):
            debug_msg = ''.join(['The input provided for PlotInjectionSites ',
                                 'in the GriddedRadialMetric plot ', name,
                                 'is not of boolean type. PlotInjectionSites will ',
                                 'be set to the default value of False.'])
            logging.debug(debug_msg)
            plot_injection_sites = DEFAULT_PLOT_INJECTION_SITES

    injection_coordx = None
    injection_coordy = None
    if 'InjectionCoordx' in yaml_data['Plots'][name][plot_type]:
        injection_coordx = yaml_data['Plots'][name][plot_type]['InjectionCoordx']

    if 'InjectionCoordy' in yaml_data['Plots'][name][plot_type]:
        injection_coordy = yaml_data['Plots'][name][plot_type]['InjectionCoordy']
        if injection_coordx is None:
            debug_msg = InjectionCoord_debug_msg.format('y', 'x')
            logging.debug(debug_msg)
            plot_injection_sites = DEFAULT_PLOT_INJECTION_SITES

        elif isinstance(injection_coordx, list) and isinstance(injection_coordy, list):

            if len(injection_coordx) != len(injection_coordy):
                debug_msg = ''.join(['The InjectionCoordy provided is not of ',
                                     'the same length as the InjectionCoordx ',
                                     'provided. Check your input. Injection ',
                                     'sites will not be displayed.'])
                logging.debug(debug_msg)
                plot_injection_sites = False

    if injection_coordx is not None and injection_coordy is None:
        debug_msg = InjectionCoord_debug_msg.format('x', 'y')
        logging.debug(debug_msg)
        plot_injection_sites = DEFAULT_PLOT_INJECTION_SITES

    if plot_injection_sites and injection_coordx is not None \
            and injection_coordy is not None:
        try:
            injection_coordx = float(injection_coordx)
        except TypeError:
            injection_coordx = np.array(list(injection_coordx), dtype=float)

        try:
            injection_coordy = float(injection_coordy)
        except TypeError:
            injection_coordy = np.array(list(injection_coordy), dtype=float)

    elif not plot_injection_sites:
        injection_coordx = None
        injection_coordy = None

    if 'SpecifyXandYLims' in yaml_data['Plots'][name][plot_type]:
        enforce_x_and_y_lims = True
        xLims = yaml_data['Plots'][name][plot_type]['SpecifyXandYLims']['xLims']
        yLims = yaml_data['Plots'][name][plot_type]['SpecifyXandYLims']['yLims']

        if len(xLims) != 2:
            debug_msg = axis_lims_debug_msg.format('x', 'x')
            logging.debug(debug_msg)
            enforce_x_and_y_lims = False

        if len(yLims) != 2:
            debug_msg = axis_lims_debug_msg.format('y', 'y')
            logging.debug(debug_msg)
            enforce_x_and_y_lims = False

    if 'Realization' in yaml_data['Plots'][name][plot_type]:
        realization = int(yaml_data['Plots'][name][plot_type]['Realization'])

        # The standard practice in the control file interface is to refer to a
        # realization with the python indexing rules (0 for the 1st entry). The
        # files saved by the GenericAquifer component do not use the python
        # indexing rules when saving output .npz files, however, so the
        # realization value is increased by one.
        realization += 1

    if 'SaveCSVFiles' in yaml_data['Plots'][name][plot_type]:
        saveCSVFiles = yaml_data['Plots'][name][plot_type]['SaveCSVFiles']

        if not isinstance(saveCSVFiles, bool):
            debug_msg = bool_option_debug_msg.format(
                'SaveCSVFiles', 'SaveCSVFiles', str(DEFAULT_SAVE_CSV_FILES))
            logging.debug(debug_msg)
            saveCSVFiles = DEFAULT_SAVE_CSV_FILES

    if 'EqualAxes' in yaml_data['Plots'][name][plot_type]:
        equal_axes = yaml_data['Plots'][name][plot_type]['EqualAxes']

        if not isinstance(equal_axes, bool):
            debug_msg = bool_option_debug_msg.format(
                'EqualAxes', 'EqualAxes', str(default_equal_axes))
            logging.debug(debug_msg)
            equal_axes = default_equal_axes

    grid_metric_yaml_input_dict = {
        'ZList': z_list,
        'MinValue': min_value,
        'TimeList': time_list,
        'DegreeInterval': degree_interval,
        'FigureDPI': figure_dpi,
        'Realization': realization,
        'comp_name_list': comp_name_list,
        'metric_name': metric_name,
        'plot_injection_sites': plot_injection_sites,
        'InjectionCoordx': injection_coordx,
        'InjectionCoordy': injection_coordy,
        'EnforceXandYLims': enforce_x_and_y_lims,
        'xLims': xLims,
        'yLims': yLims,
        'EqualAxes': equal_axes,
        'saveCSVFiles': saveCSVFiles}

    return grid_metric_yaml_input_dict


def get_t_indices(time_list, time_array):
    """
    Returns the time index corresponding to the time_array value closest to
    the times in time_list. Note that time_array is in days, while time_list is
    in years.
    """
    time_index_list = []

    if not isinstance(time_list, list):
        time_list = [time_list]

    for time in time_list:
        abs_diff = np.abs(time - (time_array / 365.25))

        # np.where returns a tuple with the first element being an array
        # of indices where the match happens (the first [0]).
        # Then we take the first element in the array of indices (the second [0])
        closest_t_index = np.where(abs_diff == np.min(abs_diff))[0][0]

        time_index_list.append(closest_t_index)

    return time_index_list


def get_z_index(z_val, z):
    """
    Returns the z index corresponding with the z value closest to z_val.
    """
    abs_diff = np.zeros(z.shape[1])

    for z_ref in range(z.shape[1]):
        abs_diff[z_ref] = np.abs(z_val - z[0, z_ref])

    closest_z_index = np.where(abs_diff == np.min(abs_diff))[0][0]

    return closest_z_index


def get_comp_z_lists(sm, z_list, name, components_to_use):
    """
    Checks if the ZList input includes string inputs referring to unit depths
    (e.g., aquifer1Depth or shale2Depth). When spatially variable stratigraphy
    is used, unit depths need to be adjusted to correspond with different
    components. To account for that possibility, this function returns a list
    of length len(components_to_use) that in turn contains a list of length
    len(z_list). Note that the numeric values provided for depths through the
    ZList entry are taken as being negative when the depth is beneath the surface.
    This practice is used to conform with the .yaml inputs used for other plot
    types in control file interface (e.g., the TTFD plot type). This function
    accounts for the input depths being negative by making them positive.
    """
    comp_z_lists = []

    for comp_ref, comp in enumerate(components_to_use):
        comp_z_lists.append([])

        for z_ref, z_val in enumerate(z_list):
            if isinstance(z_list[z_ref], str) and 'Depth' in z_list[z_ref]:
                try:
                    strat_comp_temp = sm.component_models['strata']
                except KeyError:
                    # For simulations with spatially variable stratigraphy
                    strat_comp_temp = sm.component_models['strata' + comp.name]

                strata_info = iam_strata.get_strata_info_from_component(strat_comp_temp)

                numShaleLayers = strata_info['numberOfShaleLayers']

                if 'aquifer' in z_list[z_ref]:
                    unitType = 'aquifer'
                    # The first 7 characters are 'aquifer' and the last 5
                    # characters should be 'Depth'
                    unitNumber = int(z_list[z_ref][7:-5])
                elif 'shale' in z_list[z_ref]:
                    unitType = 'shale'
                    # The first 5 characters are 'shale' and the last 5
                    # characters should be 'Depth'
                    unitNumber = int(z_list[z_ref][5:-5])

                try:
                    depth_val = iam_strata.get_unit_depth_from_component(
                        numShaleLayers, strat_comp_temp, unitNumber=unitNumber,
                        unitType=unitType, top_or_bottom='bottom')
                except:
                    debug_msg = ''.join([
                        'For the GriddedRadialMetric plot ', name, ', the input ',
                        'provided for ZList (', z_list[z_ref], ') included a ',
                        'string input. The input does not correspond to any ',
                        'unit depths available (e.g., aquifer1Depth ',
                        'or shale2Depth). The ZList input will be set ',
                        'to a value of 500 m, and the depth closest to 500 m ',
                        'in the .npz z_coordinate files will be used.'])
                    logging.debug(debug_msg)
                    # To avoid an error, the depth is set to 500 m. The depth
                    # closest to 500 m in the saved depths (see z_data below)
                    # will be used.
                    depth_val = 500

                comp_z_lists[comp_ref].append(depth_val)
            else:
                comp_z_lists[comp_ref].append(-z_list[z_ref])

    return comp_z_lists


def figure_setup(sim_index, fig_num, levels, comp_name, z_value, metric_name,
                 x_axis_limits, y_axis_limits, genfontsize=10,
                 axislabelfontsize=14, titlefontsize=14, boldlabels=True,
                 figsize=(10, 8), cmap='viridis', equal_axes=True):
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

    cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                        values=levels, format=ticker.FuncFormatter(fmt))
    cbar.set_label(
        Y_LABEL_DICT.get(metric_name, metric_name),
        rotation=90, fontsize=axislabelfontsize,
        fontweight=selected_labelfontweight)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()

    plt.xlabel('Easting (km)', fontsize=axislabelfontsize,
               fontweight=selected_labelfontweight)
    plt.ylabel('Northing (km)', fontsize=axislabelfontsize,
               fontweight=selected_labelfontweight)

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


def plot_comp_inj_sites(ax, comp_name, xloc, yloc, plot_injection_sites,
                        injection_coordx, injection_coordy, aq_component_index=None,
                        genfontsize=10, include_labels=True):
    """
    Function that plots the wells and injection sites used in the simulation.
    """
    if plot_injection_sites:
        if isinstance(injection_coordx, float):
            if include_labels:
                plt.plot(injection_coordx / 1000, injection_coordy / 1000,
                         marker='s', color='k', linestyle='none',
                         markeredgewidth=2, markersize=8,
                         markerfacecolor='none', zorder=1e6,
                         label='Injection Site')
            else:
                plt.plot(injection_coordx / 1000, injection_coordy / 1000,
                         marker='s', color='k', linestyle='none',
                         markeredgewidth=2, markersize=8,
                         markerfacecolor='none', zorder=1e6)
        else:
            for injRef, (icoordX, icoordY) in enumerate(
                    zip(injection_coordx, injection_coordy)):
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
            label = comp_name_edit + ' at Location {}'.format(loc_ind)

        if not aq_component_index is None:
            label += ', Aquifer {}'.format(aq_component_index)

        plt.plot(xloc / 1000, yloc / 1000, linestyle='none', marker='o',
                 color='k', markeredgewidth=1.5, markersize=6,
                  markerfacecolor='none', zorder=1e6, label=label)
    else:
        plt.plot(xloc / 1000, yloc / 1000, linestyle='none', marker='o',
                 color='k', markeredgewidth=1.5, markersize=6,
                  markerfacecolor='none', zorder=1e6)

    ax.legend(fancybox=False, fontsize=genfontsize,
              edgecolor=[0, 0, 0], loc='upper left',
              framealpha=0.5).set_zorder(2e6)


def save_results_to_csv(output_dir, comp_name, z_value, t_index, time, analysis,
                        sim_index, metric_name, x_vals_array, y_vals_array,
                        radial_metric_array, min_value):
    """
    Function that saves the gridded radial metric data to a .csv file.
    """

    if not os.path.exists(os.path.join(output_dir, 'csv_files')):
        os.mkdir(os.path.join(output_dir, 'csv_files'))

    if analysis == 'forward':
        filename = os.path.join(output_dir, 'csv_files',
                                '{}_{}_z{}m_tIndex{}.csv'.format(
                                    metric_name, comp_name, z_value, t_index))
    elif analysis in ['lhs', 'parstudy']:
        filename = os.path.join(output_dir, 'csv_files',
                                '{}_{}_z{}m_tIndex{}_Sim{}.csv'.format(
                                    metric_name, comp_name, z_value, t_index,
                                    sim_index))

    num_rows = len(radial_metric_array[radial_metric_array > min_value])

    # This variable is None if the simulation didn't have any values exceeding
    # MinValue
    if num_rows == 0:
        results_formatted = np.empty((2, 5), dtype=list)

        results_formatted[0, 0] = 'x (km)'
        results_formatted[0, 1] = 'y (km)'
        results_formatted[0, 2] = 'z (m)'
        results_formatted[0, 3] = 't (years)'

        results_formatted[0, 4] = CSV_FILE_COLUMNS.get(
            metric_name, metric_name).format(min_value)

        results_formatted[1, 0] = 'All x values evaluated from {:.2e} km to {:.2e} km'.format(
            np.min(x_vals_array) / 1000, np.max(x_vals_array) / 1000)

        results_formatted[1, 1] = 'All y values evaluated from {:.2e} m to {:.2e} m'.format(
            np.min(y_vals_array), np.max(y_vals_array))

        results_formatted[1, 2] = -z_value

        results_formatted[1, 3] = time

        results_formatted[1, 4] = 'No values'

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
        results_formatted[0, 2] = 'z (m)'
        results_formatted[0, 3] = 't (years)'

        # Get the column name depending on the output
        # If name does not match known observations use it as a column label
        results_formatted[0, 4] = CSV_FILE_COLUMNS.get(
            metric_name, metric_name).format(min_value)

        row_ref2 = 1
        for row_ref, (x_val, y_val) in enumerate(zip(x_vals_array, y_vals_array)):
            if radial_metric_array[row_ref] > min_value:
                results_formatted[row_ref2, 0] = str(x_val / 1000)
                results_formatted[row_ref2, 1] = str(y_val / 1000)
                results_formatted[row_ref2, 4] = str(radial_metric_array[row_ref])
                row_ref2 += 1

        results_formatted[1:None, 2] = np.ones(num_rows) * -z_value
        results_formatted[1:None, 3] = np.ones(num_rows) * time

        # Save the ouput for the simulation
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            for row_ref in range(num_rows + 1):
                writer.writerow(results_formatted[row_ref, :])
        f.close()
