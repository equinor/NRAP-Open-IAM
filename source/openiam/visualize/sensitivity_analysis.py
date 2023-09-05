# -*- coding: utf-8 -*-
"""
Created: 02/13/2018
Last modified: 03/31/2023

@author: Seth King
AECOM supporting NETL

@author: Nate Mitchell
LRST supporting NETL

"""
import os
import sys
import random
import warnings
from collections.abc import Iterable

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

from matk.sampleset import corr
from openiam import SystemModel, SimpleReservoir, MultisegmentedWellbore

from openiam.visualize.time_series import get_colors

from openiam.visualize.sensitivity_labels import (PAR_NAME_DICT, OUTPUT_DICT)

warnings.simplefilter(action='ignore', category=FutureWarning)

RC_FONT = {'family': 'Arial', 'weight': 'normal', 'size': None}

GENFONTSIZE = 12
AXISLABELFONTSIZE = 16
TITLEFONTSIZE = 16
SELECTEDFONTWEIGHT = 'bold'
LEGEND_ITEM_THRESH1 = 5
LEGEND_ITEM_THRESH2 = 10
LEGEND_ITEM_THRESH3 = 15
LEGEND_ITEM_THRESH4 = 20

def array2string(array):
    """ Convert array to string representation."""
    return np.array2string(array, separator=',', max_line_width=999999).strip('[]')


def correlations_at_time(sampleset, time_array, capture_point=1, excludes=None,
                         ctype='pearson', plot=True, printout=False,
                         plotvals=True, figsize=(10, 10), title=None,
                         xrotation=90, savefig=None, outfile=None,
                         GUI_output=False, figure_dpi=100):
    """
    Calculate correlation coefficients of parameters and responses at
    capture_point time step

    :param sampleset: Sampleset from sampling openiam.SystemModel with
        observations (has been ran)
    :type sampleset: matk.SampleSet

    :param capture_point: Time point to calculate correlation coefficients on
    :type capture_point: int or list of integers

    :param excludes: Observation base names to be excluded from calculations
    :type excludes: list of strings or a string

    :param ctype: Type of correlation coefficient (pearson by default,
        spearman also available)
    :type ctype: str

    :param plot: If True, plot correlation matrix
    :type plot: bool

    :param printout: If True, print correlation matrix with row and column headings
    :type printout: bool

    :param plotvals: If True, print correlation coefficients on plot matrix
    :type plotvals: bool

    :param figsize: Width and height of figure in inches
    :type figsize: tuple(fl64,fl64)

    :param title: Title of plot
    :type title: str

    :param xrotation: Rotation for x axis tick labels (observations)
    :type xrotation: int, float, or str

    :param savefig: Filename for figure to be saved as; figure is not saved if None
    :type savefig: str

    :param outfile: Filename for csv file output of correlation coeff; not saved if None
    :type outfile: str

    :param GUI_output: If True, figures are not closed immediately after creation
    :type GUI_output: bool

    :param figure_dpi: dpi (dots-per-inch) for the figure
    :type figure_dpi: float or int

    :returns: ndarray(fl64) -- Correlation coefficients
    """
    # Figure formatting
    font = RC_FONT
    font['size'] = GENFONTSIZE
    plt.rc('font', **font)

    if not isinstance(capture_point, Iterable):
        capture_points = [capture_point]
    else:
        capture_points = capture_point

    if excludes is None:
        excludes = []
    elif isinstance(excludes, str):
        excludes = [excludes]

    for cp in capture_points:
        resp_names = []
        for name in sampleset.responses.names:
            if name.endswith('_{cp}'.format(cp=cp)):
                resp_names.append(name)

        exc_names = ['{ex}_{cp}'.format(ex=exc, cp=cp) for exc in excludes]

        resp_names = [nm for nm in resp_names if nm not in exc_names]

        rc2 = sampleset.responses.recarray[resp_names]

        corr_coeff = corr(sampleset.samples.recarray, rc2,
                          type=ctype, plot=plot, printout=printout,
                          plotvals=plotvals, figsize=figsize,
                          title=title.format(cp=cp, ct=time_array[cp]/365.25),
                          xrotation=xrotation,
                          filename=savefig.format(cp=cp, ct=time_array[cp]/365.25),
                          adjust_dict={'left': 0.275, 'right': 0.975,
                                       'bottom': 0.25, 'top': 0.95})

        if not GUI_output:
            if savefig:
                plt.close()

        if outfile:
            with open(outfile.format(cp=cp), 'w') as ofile:
                ofile.write(' , ' + ', '.join(resp_names) + '\n')
                for i, par in enumerate(sampleset.samples.names):
                    ofile.write(par + ', ' + array2string(corr_coeff[i]) + '\n')

    return corr_coeff


def top_sensitivities(sensitivities, num_sensitivities=5):
    """
    Return a list of arguments that can be used to sort parameters
    and sensitivities from largest to smallest impact.

    :param sensitivities: Dictionary of sensitivities with 'S1' entry
    :type sensitivities: dict

    :param num_sensitivities: Number of top sensitivities to return
    :type num_sensitivities: int

    :returns sorted_args: array of indices to sort sensitivities by
    :type sorted_args: np.array of size [num_sensitivities]
    """
    sorted_args = np.flipud(np.argsort(sensitivities['S1']))
    if num_sensitivities:
        sorted_args = sorted_args[:num_sensitivities]
    return sorted_args


def time_series_sensitivities(obs_base_name, sm, lhs_sample, time_array,
                              title=None, ylabel=None, num_sensitivities=5,
                              savefig=None, capture_point=None, outfile=None,
                              GUI_output=False, use_formatted_labels=False,
                              figure_dpi=100):
    """
    Plot a time_series of top sensitivity coefficients for the observation

    :param obs_base_name: Base observation name to append time index to
        for sensitivity analysis
    :type obs_base_name: str

    :param sm: System model
    :type sm: openiam.SystemModel object

    :param lhs_sample: Sample set for system model with observations (has been ran)
    :type lhs_sample: matk.SampleSet object

    :param time_array: Time values for system model (in days)
    :type time_array: np.array

    :param title: title for figure
    :type title: str

    :param ylabel: ylabel for figure
    :type ylabel: str

    :param num_sensitivities: Number of top sensitivities to plot
    :type num_sensitivities: int

    :param savefig: Filename for figure to be saved as; figure is not saved
        if savefig is None
    :type savefig: str

    :param capture_point: Time Index to use for top sensitivity calculation.
        If None, use end time point
    :type capture_point: int

    :param outfile: Filename for csv file output of top sensitivities;
        output is not saved if outfile is None
    :type outfile: str

    :param GUI_output: If True, figures are not closed immediately after creation
    :type GUI_output: bool

    :param use_formatted_labels: option to format parameter or output names
        with a more descriptive phrase (True) or to only present the default
        name itself (False).
    :type use_formatted_labels: bool

    :param figure_dpi: dpi (dots-per-inch) for the figure
    :type figure_dpi: float or int

    :returns: raw_sensitivities
    """
    # Figure formatting
    font = RC_FONT
    font['size'] = GENFONTSIZE
    plt.rc('font', **font)

    time_array = sm.time_array

    if capture_point is not None:
        obs_name = obs_base_name + '_{0}'.format(capture_point)
    else:
        num_time = len(time_array)
        # Find top sensitivities at end time
        obs_name = obs_base_name + '_{0}'.format(num_time-1)

    end_sensitivities = lhs_sample.rbd_fast(obsname=obs_name, print_to_console=False)

    array_args = top_sensitivities(end_sensitivities,
                                   num_sensitivities=num_sensitivities)
    raw_sensitivities = []
    for time_i in range(len(time_array[1:])):
        obs_name = obs_base_name + '_{0}'.format(time_i + 1)
        sensitivities = np.array(
            lhs_sample.rbd_fast(obsname=obs_name, print_to_console=False)['S1'])
        raw_sensitivities.append(sensitivities[array_args])

    raw_sensitivities = np.array(raw_sensitivities)

    last_top = np.zeros_like(raw_sensitivities[:, 0])
    last_bottom = np.zeros_like(raw_sensitivities[:, 0])

    color_index = 0
    # This isn't used here, it is used in time_series.py
    stats_index = 0
    # This is just used to force get_colors() to ignore the stats option
    subplots_data = {'type': 'real'}
    used_colors = []

    plt.figure(figsize=(10, 8))

    plt.xlim([time_array[1] / 365.25, time_array[-1] / 365.25])

    num_items = 0
    for pi, param_i in enumerate(array_args):
        tsens_pos = raw_sensitivities[:, pi].copy()
        tsens_pos[tsens_pos < 0] = 0
        tsens_pos = tsens_pos + last_top
        tsens_neg = raw_sensitivities[:, pi].copy()
        tsens_neg[tsens_neg > 0] = 0
        tsens_neg = tsens_neg + last_bottom

        par_name = sm.parnames[param_i]
        comp_names, labels, _, _ = get_labels(
            [par_name], PAR_NAME_DICT, use_formatted_labels=use_formatted_labels,
            remove_time_index=False, include_loc_in_comp_label=True,
            capitalize_loc=True)

        label = comp_names[0] + ', ' + labels[0]

        colorValReal, _, used_colors, _, _ = get_colors(
            color_index, stats_index, used_colors, subplots_data)

        plt.fill_between(time_array[1:]/365.25, last_top, tsens_pos,
                         label=label, color=colorValReal, alpha=0.5)

        plt.fill_between(time_array[1:]/365.25, tsens_neg, last_bottom,
                         color=colorValReal, alpha=0.5)

        color_index += 1
        num_items += 1
        last_top = tsens_pos
        last_bottom = tsens_neg

    plt.subplots_adjust(bottom=0.3, top=0.9)

    ncol, legend_fontsize = get_lgnd_setup(num_items)

    plt.legend(fontsize=legend_fontsize, ncol=ncol,
               loc='upper center', bbox_to_anchor=(0.5, -0.125))

    plt.grid(alpha=0.15)

    plt.xlabel('Time (yrs)', fontsize=AXISLABELFONTSIZE,
               fontweight=SELECTEDFONTWEIGHT)

    # The previous version of the code won't work for empty label since
    # empty string also evaluates to False
    if ylabel is not None:
        plt.ylabel(ylabel, fontsize=AXISLABELFONTSIZE,
                   fontweight=SELECTEDFONTWEIGHT)
    else:
        plt.ylabel('Stacked First Order Sensitivities', fontsize=AXISLABELFONTSIZE,
                   fontweight=SELECTEDFONTWEIGHT)

    comp_names, labels, _, _ = get_labels(
        [obs_base_name], PAR_NAME_DICT, use_formatted_labels=use_formatted_labels,
        remove_time_index=False, include_loc_in_comp_label=True,
        capitalize_loc=True, two_lines=False)

    if title is not None:  # see comment above for y-axis label
        plt.title(title.format(ob=labels[0], comp=comp_names[0]),
                  fontsize=TITLEFONTSIZE, fontweight=SELECTEDFONTWEIGHT)
    else:
        plt.title('Sensitivity Over Time for\nthe Output {ob} of Component {comp}'.format(
            ob=labels[0], comp=comp_names[0]), fontsize=TITLEFONTSIZE,
            fontweight=SELECTEDFONTWEIGHT)

    if savefig:
        plt.savefig(savefig, dpi=figure_dpi)

        if not GUI_output:
            plt.close()

    if outfile:
        with open(outfile, 'w') as ofile:
            ofile.write(' , ' + array2string(time_array[1:]/365.25) + '\n')
            for pi, param_i in enumerate(array_args):
                ofile.write(sm.parnames[param_i] + ',')
                ofile.write(array2string(raw_sensitivities[:, pi]) + '\n')

    return raw_sensitivities


def multi_sensitivities_barplot(obs_names, system_model, lhs_sample, title=None,
                                ylabel=None, savefig=None, outfile=None,
                                GUI_output=False, use_formatted_labels=False,
                                figure_dpi=100):
    """
    Calculate sensitivities and make side-by-side bar chart out of values

    :param obs_names: list of observation names to calculate sensitivities on
    :type obs_names: list of strings or a string

    :param system_model: System model for which sensitivity analysis was ran.
    :type system_model: openiam.SystemModel

    :param lhs_sample: Sample set for system model with observations (has been ran)
    :type lhs_sample: matk.SampleSet object

    :param title: title for figure
    :type title: str

    :param ylabel: ylabel for figure
    :type ylabel: str

    :param savefig: Filename for figure to be saved as; figure is not saved
        if savefig is None
    :type savefig: str

    :param outfile: Filename for csv file output of sensitivities;
        output is not saved if outfile is None
    :type outfile: str

    :param GUI_output: If True, figures are not closed immediately after creation
    :type GUI_output: bool

    :param use_formatted_labels: option to format parameter or output names
        with a more descriptive phrase (True) or to only present the default
        name itself (False).
    :type use_formatted_labels: bool

    :param figure_dpi: dpi (dots-per-inch) for the figure
    :type figure_dpi: float or int

    :returns: raw_sensitivities
    """
    # Figure formatting
    font = RC_FONT
    font['size'] = GENFONTSIZE
    plt.rc('font', **font)

    time_array = system_model.time_array

    if isinstance(obs_names, str):
        obs_names = [obs_names]

    total_obs = len(obs_names)
    total_params = len(system_model.parnames)

    width = 0.8

    fig = plt.figure(figsize=(10, 8))
    ax = plt.subplot(1, 1, 1)

    color_index = 0
    # This isn't used here, it is used in time_series.py
    stats_index = 0
    # This is just used to force get_colors() to ignore the stats option
    subplots_data = {'type': 'real'}
    used_colors = []

    num_items = 0

    output_time_indices = []
    raw_sensitivities = {}
    for obs_i, obs_name in enumerate(obs_names):
        sens = lhs_sample.rbd_fast(obsname=obs_name, print_to_console=False)['S1']
        x_values = [total_obs*x + width*obs_i for x in range(total_params)]

        comp_name, label, time_index, _ = get_labels(
            [obs_name], OUTPUT_DICT, use_formatted_labels=use_formatted_labels,
            include_comp_in_label=False, include_loc_in_comp_label=True,
            capitalize_loc=False, two_lines=False)

        output_time_indices.append(time_index[0])
        label = label[0] + ' from ' + comp_name[0]
        label = add_time_index_str(label, time_index, time_array,
                                   include_word_at=False, include_word_time=False,
                                   capitalize_time=False, include_comma=True)

        colorValReal, _, used_colors, _, _ = get_colors(
            color_index, stats_index, used_colors, subplots_data)

        plt.bar(x_values, sens, label=label, color=colorValReal)

        color_index += 1
        num_items += 1

        raw_sensitivities[obs_name] = np.array(sens)

    x_values = [total_obs*x + 0.8*(total_obs/2) for x in range(total_params)]

    ax.set_xticks(x_values)

    par_names = system_model.parnames
    _, labels, _, _ = get_labels(
        par_names, PAR_NAME_DICT, use_formatted_labels=use_formatted_labels,
        remove_time_index=False, include_comp_in_label=True,
        include_loc_in_comp_label=True, capitalize_loc=True, two_lines=True)

    ax.set_xticklabels(labels, rotation=90)

    if title is None:
        title = 'Output Sensitivity'
        title = add_time_index_str(title, output_time_indices, time_array,
                                   include_word_at=True, include_word_time=True,
                                   capitalize_time=True)

    plt.title(title, fontsize=TITLEFONTSIZE, fontweight=SELECTEDFONTWEIGHT)

    if ylabel is None:
        ylabel = 'First Order Sensitivity'

    plt.ylabel(ylabel, fontsize=AXISLABELFONTSIZE, fontweight=SELECTEDFONTWEIGHT)

    plt.subplots_adjust(bottom=0.35, top=0.8)

    ncol, legend_fontsize = get_lgnd_setup(num_items)

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=legend_fontsize, ncol=ncol)

    plt.grid(axis='y', linestyle='--', linewidth=1, alpha=0.3)

    if savefig:
        plt.savefig(savefig.format(), dpi=figure_dpi)

        if not GUI_output:
            plt.close()

    if outfile:
        with open(outfile, 'w') as ofile:
            ofile.write(' ,' + ', '.join(system_model.parnames) + '\n')
            for obs_name in obs_names:
                ofile.write(obs_name + ',' +
                            array2string(raw_sensitivities[obs_name]) + '\n')

    return raw_sensitivities


def simple_sensitivities_barplot(sensitivities, system_model, obs=None,
                                 title=None, ylabel=None, savefig=None,
                                 outfile=None, name=None, GUI_output=False,
                                 use_formatted_labels=False, figure_dpi=100):
    """
    Makes simple bar chart out of sensitivity values

    :param sensitivities: dictionary of sensitivity output from SAlib analysis
    :type sensitivities: dict

    :param system_model: System model sensitivity analysis were ran for.
    :type system_model: openiam.SystemModel

    :param title: title for figure
    :type title: str

    :param ylabel: ylabel for figure
    :type ylabel: str

    :param savefig: Filename for figure to be saved as; figure is not saved if None
    :type savefig: str

    :param outfile: Filename for text file output of sensitivies; output
        is not saved if None
    :type outfile: str

    :param GUI_output: If True, figures are not closed immediately after creation
    :type GUI_output: bool

    :param use_formatted_labels: option to format parameter or output names
        with a more descriptive phrase (True) or to only present the default
        name itself (False).
    :type use_formatted_labels: bool

    :param figure_dpi: dpi (dots-per-inch) for the figure
    :type figure_dpi: float or int

    :returns: sensitivities['S1']
    """
    # Figure formatting
    font = RC_FONT
    font['size'] = GENFONTSIZE
    plt.rc('font', **font)

    time_array = system_model.time_array

    if not name:
        name = ' '

    plt.figure(name, figsize=(10, 8))
    ax = plt.subplot(1, 1, 1)

    par_names = system_model.parnames
    x_vals = list(range(len(par_names)))

    plt.bar(x_vals, sensitivities['S1'])
    ax.set_xticks(x_vals)

    comp_names, labels, time_indices, _ = get_labels(
        par_names, PAR_NAME_DICT, use_formatted_labels=use_formatted_labels,
        remove_time_index=False, include_comp_in_label=True,
        include_loc_in_comp_label=True, capitalize_loc=True, two_lines=True)

    ax.set_xticklabels(labels, rotation=90)
    plt.subplots_adjust(bottom=0.4, top = 0.9)

    if title is None:
        if obs is not None:
            comp_names, labels, time_indices, _ = get_labels(
                [obs], OUTPUT_DICT, use_formatted_labels=use_formatted_labels,
                remove_time_index=True, include_comp_in_label=False,
                include_loc_in_comp_label=True, capitalize_loc=True,
                two_lines=False)

            title = labels[0] + ' from ' + comp_names[0] + ',\n'
            title = add_time_index_str(title, time_indices, time_array,
                                       include_word_at=False, include_word_time=True,
                                       capitalize_time=True)
        else:
            title = 'Sensitivity Coefficients'

    plt.title(title, fontsize=TITLEFONTSIZE, fontweight=SELECTEDFONTWEIGHT)

    if ylabel is None:
        ylabel = 'First Order Sensitivity\nof Output to Component Parameters'

    plt.ylabel(ylabel, fontsize=AXISLABELFONTSIZE, fontweight=SELECTEDFONTWEIGHT)

    plt.grid(axis='y', linestyle='--', linewidth=1, alpha=0.3)

    if savefig:
        plt.savefig(savefig, dpi=figure_dpi)

        if not GUI_output:
            plt.close()

    if outfile:
        with open(outfile, 'w') as ofile:
            if obs is not None:
                ofile.write('Parameter,Output_Sensitivity_to_Parameter\n')
            for i, par in enumerate(system_model.parnames):
                ofile.write(par + ',{0}\n'.format(sensitivities['S1'][i]))

    return sensitivities['S1']


def stacked_sensitivities_barplot(sensitivities, names, system_model,
                                  title=None, ylabel=None, savefig=None,
                                  outfile=None, GUI_output=False,
                                  use_formatted_labels=False, figure_dpi=100):
    """
    Makes simple bar chart out of sensitivity values

    :param sensitivities: list of dictionaries of sensitivity output
        from SAlib analysis
    :type sensitivities: lst

    :param names: List of names for sensitivity labels
    :type names: lst

    :param system_model: System model sensitivities were ran on.
    :type system_model: openiam.SystemModel

    :param title: title for figure
    :type title: str

    :param ylabel: ylabel for figure
    :type ylabel: str

    :param savefig: Filename for figure to be saved as; figure is not saved if None
    :type savefig: str

    :param outfile: Filename for csv file output of top sensitivities;
        output is not saved if outfile is None
    :type outfile: str

    :param GUI_output: If True, figures are not closed immediately after creation
    :type GUI_output: bool

    :param use_formatted_labels: option to format parameter or output names
        with a more descriptive phrase (True) or to only present the default
        name itself (False).
    :type use_formatted_labels: bool

    :param figure_dpi: dpi (dots-per-inch) for the figure
    :type figure_dpi: float or int

    :returns: None
    """
    # Figure formatting
    font = RC_FONT
    font['size'] = GENFONTSIZE
    plt.rc('font', **font)

    fig = plt.figure(figsize=(10, 8))
    ax = plt.subplot(1, 1, 1)

    x_vals = list(range(len(sensitivities)))

    par_names = system_model.parnames

    color_index = 0
    # This isn't used here, it is used in time_series.py
    stats_index = 0
    # This is just used to force get_colors() to ignore the stats option
    subplots_data = {'type': 'real'}
    used_colors = []

    bottom = np.array([0 for x in x_vals])
    bottom_pos = np.array([0 for x in x_vals])
    bottom_neg = np.array([0 for x in x_vals])
    for i, param in enumerate(par_names):
        colorValReal, _, used_colors, _, _ = get_colors(
            color_index, stats_index, used_colors, subplots_data)
        # print('param: ', param)
        # print('    colorValReal: ', colorValReal)
        # print('    used_colors: ', used_colors)
        _, labels, time_indices, _ = get_labels(
            [param], PAR_NAME_DICT, use_formatted_labels=use_formatted_labels,
            remove_time_index=False, include_comp_in_label=True,
            include_loc_in_comp_label=True, capitalize_loc=True, two_lines=False)

        y = np.array([x['S1'][i] for x in sensitivities])

        for yref, yval in enumerate(y):
            if yval >= 0:
                bottom[yref] = bottom_pos[yref]
                bottom_pos[yref] += yval

            else:
                bottom_neg[yref] += yval
                bottom[yref] = bottom_neg[yref]

        plt.bar(x_vals, y, bottom=bottom, label=labels[0], color=colorValReal)
        color_index +=1

    _, labels, time_indices, _ = get_labels(
        names, OUTPUT_DICT, use_formatted_labels=use_formatted_labels,
        remove_time_index=True, include_comp_in_label=True,
        include_loc_in_comp_label=True, capitalize_loc=True, two_lines=True)

    ax.set_xticks(x_vals)
    ax.set_xticklabels(labels, rotation=90)

    plt.subplots_adjust(bottom=0.35, top=0.8)

    num_items = len(par_names)
    ncol, legend_fontsize = get_lgnd_setup(num_items)

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=legend_fontsize, ncol=ncol)
               # loc='upper center', bbox_to_anchor=(0.5, 0.99))

    if title is None:
        time_array = system_model.time_array

        _, labels, time_indices, _ = get_labels(
            names, OUTPUT_DICT, use_formatted_labels=use_formatted_labels,
            remove_time_index=True, include_comp_in_label=False,
            include_loc_in_comp_label=True, capitalize_loc=True,
            two_lines=False)

        title = 'Output Sensitivity '
        title = add_time_index_str(title, time_indices, time_array,
                                   include_word_at=True, include_word_time=True,
                                   capitalize_time=True)

    plt.title(title, fontsize=TITLEFONTSIZE, fontweight=SELECTEDFONTWEIGHT)

    if ylabel is None:
        ylabel = 'Stacked First Order Sensitivities\nof Output to Component Parameters'

    plt.ylabel(ylabel, fontsize=AXISLABELFONTSIZE, fontweight=SELECTEDFONTWEIGHT)

    plt.grid(axis='y', linestyle='--', linewidth=1, alpha=0.3)

    if savefig:
        plt.savefig(savefig, dpi=figure_dpi)

        if not GUI_output:
            plt.close()

    if outfile:
        with open(outfile, 'w') as ofile:
            ofile.write('Output_Name,Parameter_Name,Output_Sensitivity_to_Parameter\n')
            for i, par in enumerate(par_names):
                for sens_ref, sens_val in enumerate(sensitivities):
                    ofile.write(names[sens_ref] + ',' + par + ',' \
                                + '{0}\n'.format(sens_val['S1'][i]))

def get_labels(names, label_dict, remove_comp_name=True, remove_loc_index=True,
               remove_time_index=True, use_formatted_labels=False,
               include_comp_in_label=False, include_loc_in_comp_label=False,
               capitalize_loc=False, two_lines=True):
    """
    Get label for a given parameter or output. Note that a parameter name is
    formatted like 'SimpleReservoir1_000.injRate, where SimpleReservoir1_000 is
    the component name, _000 is the location index (index 0), and injRate is a
    specific parameter of that component. Output names are formatted like
    'SimpleReservoir1_000.pressure_0,' where pressure is the output name and
    '_0' is the time index.

    :param par_names: list of parameter or output names
    :type par_name: list

    :param labels: dictionary containing labels corresponding to parameter
        or output names
    :type labels: dict()

    :param remove_comp_name: option to remove the component name from the
        names in the name list.
    :type remove_comp_name: bool

    :param remove_loc_index: option to remove the location index from the
        names in the name list.
    :type remove_loc_index: bool

    :param remove_time_index: option to remove the time indices from the
        names in the name list. This should be set to False when handling
        parameters and True when handlig outputs.
    :type remove_time_index: bool

    :param use_formatted_labels: option to format parameter or output names
        with a more descriptive phrase (True) or to only present the default
        name itself (False).
    :type use_formatted_labels: bool

    :param include_comp_in_label: option to include the component name in each
        label name. For example, 'SimleReservoir1,\ninjRate.'
    :type include_comp_in_label: bool

    Returns:
        string to use for parameter label
    """
    comp_names = []
    labels = []
    time_indices = []
    loc_indices = []

    for name in names:
        # example of a parameter name 'SimpleReservoir1_000.logResPerm'
        if '.' in name:
            comp_name = name[0:name.index('.')]
        else:
            comp_name = ''

        if '_' in comp_name:
            loc_index = comp_name[(comp_name.index('_') + 1):None]

            if remove_loc_index:
                comp_name = comp_name[0:comp_name.index('_')]
        else:
            loc_index = ''

        loc_indices.append(loc_index)

        # In the example above, SimpleReservoir1_000 is the component name
        if '.' in name and remove_comp_name:
            name = name[(name.index('.') + 1):None]

        if '_' in name:
            if name[-2] == '_':
                time_index = name[-1]
            elif name[-3] == '_':
                time_index = name[-2:None]
            elif name[-4] == '_':
                time_index = name[-3:None]
            else:
                time_index = ''

            time_indices.append(time_index)

        # Remove the time index (e.g., '_10')
        if remove_time_index:
            name = name[0:(-len(time_index) - 1)]

        if name in label_dict and use_formatted_labels:
            str_label = label_dict[name]
        elif name[0:-1] in label_dict and use_formatted_labels:
            # This is for metrics like 'CO2_aquifer1'
            str_label = label_dict[name[0:-1]] + ' ' + name[-1]
        elif name[0:-2] in label_dict and use_formatted_labels:
            # This is for metrics like 'brine_aquifer15'
            str_label = label_dict[name[0:-2]] + ' ' + name[-2:None]
        else:
            str_label = name

        if include_loc_in_comp_label:
            if not isinstance(loc_index, str):
                if capitalize_loc:
                    comp_name += ' at Location {:.0f}'.format(int(loc_index))
                else:
                    comp_name += ' at location {:.0f}'.format(int(loc_index))

        if include_comp_in_label:
            if two_lines:
                str_label = comp_name + ',\n' + str_label
            else:
                str_label = comp_name + ', ' + str_label

        comp_names.append(comp_name)
        labels.append(str_label)

    return comp_names, labels, time_indices, loc_indices


def get_comp_title_str(comp_names):
    """
    Get string to show the component names in a figure title. Note that component
    names are formatted like 'SimpleReservoir1_000', where 'SimpleReservoir1'
    is specified by the user and '_000' designates the location index as 0.

    :param comp_names: list of component names
    :type comp_name: list

    """
    if len(set(comp_names)) > 2:
        comp_title_str = ''
        used_comp_names = []
        for cm_name in comp_names:
            if cm_name not in used_comp_names:
                comp_title_str += cm_name + ', '
                used_comp_names.append(cm_name)
                last_added = cm_name

        # remove last comma and add ' and ' before the last comp_name
        comp_title_str = comp_title_str[
            0:(-2 - (len(last_added)))] + 'and ' + comp_title_str[
                (-2 - (len(last_added))):-2]

        comp_title_str += ' Sensitivities'

    elif len(set(comp_names)) == 2:
        comp_title_str = ''
        used_comp_names = []
        for cm_name in comp_names:
            if cm_name not in used_comp_names:
                if comp_title_str == '':
                    comp_title_str += cm_name + ' and '
                else:
                    comp_title_str += cm_name
                used_comp_names.append(cm_name)

        comp_title_str += ' Sensitivities'

    elif len(set(comp_names)) == 1:
        comp_title_str = comp_names[0] + ' Sensitivity'

    return comp_title_str


def add_time_index_str(initial_str, time_indices, time_array,
                       capitalize_time=True, include_word_at=False,
                       include_word_time=False, include_comma=False):
    """
    Function that adds the time index or indices to a given string. The index
    or indices come from an output name (e.g., '10' in
    SimpleReservoir1_000.pressure_10).

    :param initial_str: string to which the time index or indices will be added
    :type initial_str: str

    :param time_indices: list of time indices
    :type time_indices: list

    :param time_array: time array used for the simulation
    :type time_array: array
    """
    if include_word_at:
        time_indices_str_i = ' at '
    else:
        if include_comma:
            time_indices_str_i = ', '
        else:
            time_indices_str_i = ''

    time_indices_str = ''
    if len(time_indices) > 0:
        if len(set(time_indices)) > 2:
            time_indices_str = ''
            used_time_indices = []
            for time_ind in time_indices:
                if time_ind not in used_time_indices:
                    time_indices_str += '{}'.format(
                        time_array[int(time_ind)] / 365.25) + ', '
                    used_time_indices.append(time_ind)
                    last_added = time_ind

            # remove last ', ' and add ' and ' before the last time index
            time_indices_str = time_indices_str[
                0:(-2 - (len(last_added)))] + 'and ' + '{}'.format(
                    time_array[int(time_indices_str[(-2 - (len(
                        last_added))):-2])] / 365.25) + ' years'

            if include_word_time:
                if capitalize_time:
                    time_indices_str = time_indices_str_i + 'Time t = ' + time_indices_str
                else:
                    time_indices_str = time_indices_str_i + 'time t = ' + time_indices_str
            else:
                time_indices_str = time_indices_str_i + 't = ' + time_indices_str

        elif len(set(time_indices)) == 2:
            time_indices_str ='{}'.format(time_array[int(time_indices[
                0])] / 365.25) + ' and ' + '{}'.format(time_array[int(
                    time_indices[1])] / 365.25) + ' years'

            if include_word_time:
                if capitalize_time:
                    time_indices_str = time_indices_str_i + 'Time t = ' + time_indices_str
                else:
                    time_indices_str = time_indices_str_i + 'time t = ' + time_indices_str
            else:
                time_indices_str = time_indices_str_i + 't = ' + time_indices_str

        elif len(set(time_indices)) == 1:
            if time_indices[0]:
                if include_word_time:
                    if capitalize_time:
                        time_indices_str = time_indices_str_i + 'Time t = {}'.format(
                            time_array[int(time_indices[0])] / 365.25) + ' years'
                    else:
                        time_indices_str = time_indices_str_i + 'time t = {}'.format(
                            time_array[int(time_indices[0])] / 365.25) + ' years'
                else:

                    time_indices_str = time_indices_str_i + 't = {}'.format(
                        time_array[int(time_indices[0])] / 365.25) + ' years'

        initial_str += time_indices_str

    return initial_str


def get_lgnd_setup(num_items):
    """
    This function evaluates the number of items in the legend and returns the
    legend fontsize and number of columns used.
    """
    ncol = 1
    legend_fontsize = GENFONTSIZE

    if LEGEND_ITEM_THRESH1 < num_items <= LEGEND_ITEM_THRESH2:
        ncol = 2
        legend_fontsize = GENFONTSIZE - 2

    elif LEGEND_ITEM_THRESH2 < num_items <= LEGEND_ITEM_THRESH3:
        ncol = 2
        legend_fontsize = GENFONTSIZE - 4

    if LEGEND_ITEM_THRESH3 < num_items <= LEGEND_ITEM_THRESH4:
        ncol = 2
        legend_fontsize = GENFONTSIZE - 6

    elif num_items > LEGEND_ITEM_THRESH4:
        ncol = 3
        legend_fontsize = GENFONTSIZE - 8

    return ncol, legend_fontsize


if __name__ == "__main__":
    import os
    __spec__ = None
    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}  # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))

    # Add parameters of reservoir component model
    sres.add_par('numberOfShaleLayers', value=3, vary=False)
    sres.add_par('shale1Thickness', value=850.0, vary=False)
    sres.add_par('shale2Thickness', value=400.0, vary=False)
    sres.add_par('shale3Thickness', value=20.0, vary=False)
    sres.add_par('aquifer1Thickness', value=45.0, vary=False)
    sres.add_par('aquifer2Thickness', value=120, vary=False)
    sres.add_par('reservoirThickness', value=65.0, vary=False)
    sres.add_par('reservoirPorosity', value=0.3, min=0.2, max=0.4)
    sres.add_par('injRate', value=0.1, min=0.2, max=0.25)

    # Add observations of reservoir component model
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')
    sres.add_obs_to_be_linked('pressure')
    sres.add_obs_to_be_linked('CO2saturation')

    # Add open wellbore component
    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))
    # A lot of parameters of multisegmented wellbore component
    # are the same as for the reservoir component
    # Add parameters linked to the same parameters from reservoir model
    ms.add_par_linked_to_par(
        'numberOfShaleLayers', sres.deterministic_pars['numberOfShaleLayers'])
    ms.add_par_linked_to_par(
        'shale1Thickness', sres.deterministic_pars['shale1Thickness'])
    ms.add_par_linked_to_par(
        'shale2Thickness', sres.deterministic_pars['shale2Thickness'])
    ms.add_par_linked_to_par(
        'shale3Thickness', sres.deterministic_pars['shale3Thickness'])
    ms.add_par_linked_to_par(
        'aquifer1Thickness', sres.deterministic_pars['aquifer1Thickness'])
    ms.add_par_linked_to_par(
        'aquifer2Thickness', sres.deterministic_pars['aquifer2Thickness'])
    ms.add_par_linked_to_par(
        'reservoirThickness', sres.deterministic_pars['reservoirThickness'])
    ms.add_par_linked_to_par('datumPressure', sres.default_pars['datumPressure'])

    # Add parameters specific to the wellbore component
    ms.add_par('logWellPerm', min=-14., max=-11., value=-13.5)
    ms.add_par('logAquPerm', min=-13.8, max=-10.3, value=-12.0)

    # Add keyword arguments linked to reservoir data
    ms.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

    # Add observations of the wellbore component
    for ind in [1, 2]:
        ms.add_obs('CO2_aquifer{}'.format(ind))
        ms.add_obs('brine_aquifer{}'.format(ind))

    num_samples = 300
    ncpus = 5
    # Draw Latin hypercube samples of parameter values
    lhs_sample = sm.lhs(siz=num_samples, seed=random.randint(500, 1100))

    lhs_sample.run(cpus=ncpus, verbose=False)

    output_directory = 'test_sensitivity_analysis'
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # Run RBD_fast sensitivity analysis on leakage rates of CO2 into aquifer 1
    rates_sensitivity1 = lhs_sample.rbd_fast(obsname='ms.CO2_aquifer1_24',
                                             print_to_console=False)
    results1 = simple_sensitivities_barplot(
        rates_sensitivity1, sm, name=1,
        title='RBD-Fast Sensitivity\nfor Leakage Rates of $CO_2$ into the Aquifer 1',
        ylabel='First Order Sensitivity\n(Calculated using RBD-Fast Method)',
        savefig=os.path.join(output_directory,
                             'RBD-fast_co2_leakage_rates_plot.png'),
        GUI_output=True)

    # Add capture point so sensitivity is calculated at time other than final
    capture_point = 24  # in years
    cp_index = np.argmin(np.abs(time_array - capture_point*365.25))
    sensitivity_obs = 'ms.brine_aquifer1'
    sen_obs_name = sensitivity_obs + '_{0}'.format(cp_index)

    # Run sensitivity analysis on brine leakage rates to aquifer 1 at capture point
    rates_sensitivity2 = lhs_sample.rbd_fast(obsname=sen_obs_name,
                                             print_to_console=False)

    results2 = simple_sensitivities_barplot(
        rates_sensitivity2, sm, name=2,
        title='First Order Sensitivity\nfor Brine Leakage Rates to Aquifer 1',
        ylabel='First Order Sensitivity\n(Calculated using RBD-Fast Method)',
        savefig=os.path.join(output_directory,
                             'RBD-fast_brine_leakage_rates_plot.png'),
        outfile=os.path.join(output_directory,
                             'brine_aquifer1_sensitivities.csv'),
        GUI_output=True)

    results3 = time_series_sensitivities(
        'ms.brine_aquifer1', sm, lhs_sample, time_array,
        savefig=os.path.join(output_directory,
                             'RBD-fast_brine_leak_rates_time_series.png'),
        outfile=os.path.join(output_directory,
                             'brine_leak_rates_aquifer_time_sens.csv'),
        GUI_output=True)

    results4 = time_series_sensitivities(
        'ms.CO2_aquifer1', sm, lhs_sample, time_array,
        savefig=os.path.join(output_directory,
                             'RBD-fast_co2_leak_rates_time_series.png'),
        outfile=os.path.join(output_directory,
                             'co2_leak_rates_aquifer_time_sens.csv'),
        GUI_output=True)

    results5 = multi_sensitivities_barplot(
        ['sres.pressure_40', 'sres.CO2saturation_40',
         'ms.CO2_aquifer1_40', 'ms.brine_aquifer1_40',
         'ms.CO2_aquifer2_40', 'ms.brine_aquifer2_40'],
        sm, lhs_sample,
        title='Sensitivity Coefficients at 40 years',
        savefig=os.path.join(output_directory,
                             'multi_bar_sensitivities_time_index_40.png'),
        outfile=os.path.join(output_directory,
                             'multi_bar_sens_time_index_40.csv'),
        GUI_output=True)

    corr_coeffs = correlations_at_time(
        lhs_sample, time_array, capture_point=25,
        excludes=['sres.CO2saturation', 'sres.pressure', 'ms.brine_aquifer2'],
        plot=True, figsize=(15, 15), printout=False, xrotation=90,
        title='Pearson Correlation Coefficients at 25 years',
        savefig=os.path.join(output_directory,
                             'corr_coeff_at_time_index_25.png'),
        outfile=os.path.join(output_directory,
                             'corr_coeffs_at_time_index_25.csv'),
        GUI_output=True)
