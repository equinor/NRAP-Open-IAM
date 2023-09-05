# -*- coding: utf-8 -*-
"""
Code to create different time-series visualizations for NRAP-Open-IAM.

Last Modified: June, 2023

@author: Seth King
@author: Nate Mitchell (Nathaniel.Mitchell@NETL.DOE.GOV)
@author: Veronika Vasylkivska (Veronika.Vasylkivska@NETL.DOE.GOV)
LRST (Battelle/Leidos) supporting NETL
"""
import warnings
import logging
import math
import re

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cbook
import matplotlib.colors as clrs
from matk.sampleset import percentile, mean
from .label_setup import (LEGEND_DICT, Y_LABEL_DICT, Y_LABEL_SPLIT_DICT,
                          Y_LABEL_2ROWS_DICT, Y_LABEL_2ROWS_SPLIT_DICT,
                          TITLE_DICT, TITLE_SPLIT_DICT)

# Ignore futurewarning from numpy about record array subarray copies vs views.
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=matplotlib.cbook.mplDeprecation)

# Constants used to adjust figure formatting
DEFAULT_FIG_WIDTH = 13
DEFAULT_FIG_HEIGHT = 8
AXIS_LABEL_PAD_REF = 4
TITLE_PAD_REF = 3
SINGLE_PLOT_FONTSIZE_ADJUST = 1.5

AX_LINEWIDTH_REF = 1.5
LINEWIDTH_REF = 2.5
XTICK_SIZE_REF = 5
YTICK_SIZE_REF = 5
XTICK_WIDTH_REF = 1.5
YTICK_WIDTH_REF = 1.5

# These are used to check if the x and y labels are too long relative to the subplot.
# These values were found by trial and error.
MAX_YLABEL_HEIGHT_FRAC = 0.85
MAX_YLABEL_WIDTH_FRAC = 0.075

MAX_XLABEL_HEIGHT_FRAC = 0.4
MAX_XLABEL_WIDTH_FRAC = 0.5

MAX_TITLE_HEIGHT_FRAC = 0.1
MAX_TITLE_WIDTH_FRAC = 0.4

LEGEND_ITEM_THRESH1 = 5
LEGEND_ITEM_THRESH2 = 10
LEGEND_ITEM_THRESH3 = 15
LEGEND_ITEM_THRESH4 = 20

LEGEND_FRAME_ALPHA = 0.5
LEGEND_COLUMNS = 1

# I purposefully included more than 10 marker styles so that the rotation of
# marker styles (when plotting a large number of metrics) is desynchronized from
# the rotation of line colors. The default rotation of 10 line colors is used,
# then two rotations of darker and lighter versions of the same 10 colors are
# used. Having the rotations of marker styles and colors desynchronized may
# help clarify such distinctions. Most simulations wil use many years of data
# with a 1 year timestep, however, in which case markers may be too close together.
defaultMarkerList = ['o', '^', 's', 'd', 'X', 'P', '*', 'p', 'D', 'H',
                     '>', '<', 'v']
defaultLineStyleList = ['solid', 'dotted', 'dashed', 'dashdot']
hatchList = ['|', '-', '+', 'x', 'o', 'O', '.', '*', '/', '\\']

def time_series_plot(output_names, sm, s, plot_data, output_list, name='Figure 1',
                     analysis='forward', savefig=None, title=None, subplot=None,
                     plot_type=None, figsize=None, fontname='Arial',
                     gen_font_size=10, axis_label_font_size=12,
                     title_font_size=12, legend_font_size=10, bold_labels=True,
                     useMathTextOption=True, generate_title=True,
                     plot_grid_option=True, grid_alpha_val=0.15):
    """
    Makes a time series plots of statistics of data in output_names.
    Highlights percentiles, mean, and median values.

    :param output_names: List of observation names to match with component models and plot.
    Examples:
        output_names=['pressure']

        output_names=['pressure', 'mass_CO2_aquifer1']

        output_names=['pressure', 'CO2_aquifer1', 'CO2_aquifer2']
    :type output_names: list

    :param sm: OpenIAM System model for which plots are created
    :type sm: openiam.SystemModel object

    :param s: SampleSet with results of analysis performed on the system model.
        None for forward analysis.
    :type s: matk.SampleSet object

    :param plot_data: Dictionary of setup for a given plot
    :type plot_data: dict

    :param output_list: Dictionary mapping component models to observations
    Examples:
        This dictionary only includes pressure from a SimpleReservoir named sres.
        output_list={sres: 'pressure'}

        This dictionary includes pressure from a SimpleReservoir named sres as well
        as the mass of CO2 leaked to aquifer 1. The CO2 mass comes from a RateToMassAdapter
        named adapt.
        output_list={sres: 'pressure', adapt: 'mass_CO2_aquifer1'}

        This dictionary includes pressure from a SimpleReservior named sres as well
        as well as the CO2 leakage rates to aquifers 1 and 2. The CO2 leakage rates
        come from a MultisegmentedWellbore named ms.
        output_list={sres: 'pressure', ms: ['CO2_aquifer1', 'CO2_aquifer2']}
    :type output_list: dict

    :param name: Figure Name to be used/created.
    :type name: str

    :param analysis: Type of OpenIAM system analysis performed ('forward',
        'lhs', or 'parstudy')
    :type analysis: str

    :param savefig: Filename to save resulting figure to. No file saved if None.
    :type savefig: str

    :param title: Optional Title for figure
    :type title: str

    :param subplot: Dictionary for subplot controls, use=True will use
        subplotting (boolean default=False), ncols=n will use n columns
        of subplots (positive integer default 1); comp.var_name=sub_title
        will title the subplots of comp.var_name subtitle (string default=comp.var_name).
        Examples:
            This dictionary creates only one plot (i.e., no subplots). If you are
            plotting multiple types of metrics (e.g., pressure and CO2 leakage
            rates), all metrics will be displayed together and only the yaxis
            label for the metric plotted last will be shown. Do not use
            this setup in such a case: subplot={'use': False}

            This dictionary enables the creation of subplots and specifies
            the use of only one column: subplot={'use': True, 'ncols': 1}

            This dictionary includes figure titles for sres.pressure, ms.CO2_aquifer1, and
            ms.CO2_aquifer2.
            subplot={'use': True, 'ncols': 3, 'sres.pressure': 'Pressure at location',
                     'ms.CO2_aquifer1': 'CO$_2$ Leakage Rate to Aquifer 1',
                     'ms.CO2_aquifer2': 'CO$_2$ Leakage Rate to Aquifer 2'}
    :type subplot: dict

    :param plot_type: List of 'real' and/or 'stats' plot types to produce.
        'real' plots realization values
        'stats' plots quartiles, mean, and median values
    :type plot_type: list of 'real' and/or 'stats'

    :param figsize: width and height of the figure (width, height), in inches. Default is
        None, in which case the DEFAULT_FIG_WIDTH and DEFAULT_FIG_HEIGHT are used.
    :type figsize: tuple or NoneType

    :param fontname: name of the font type to be used
    :type fontname: str

    :param gen_font_size: fontsize for tick labels, etc.
    :type gen_font_size: float or int

    :param axis_label_font_size: fontsize for x and y axis labels. These font sizes are
        later updated depending on the figsize and number of subplots.
    :type axis_label_font_size: float or int

    :param title_font_size: fontsize for the title. This font size is later
        updated depending on the figsize and number of subplots.
    :type title_font_size: float or int

    :param legend_font_size: fontsize for the legend. This font size is later
        updated depending on the figsize.
    :type legend_font_size: float or int

    :param bold_labels: option to use bold x and y labels and bold titles. Set to
        True for bold labels, False for normal labels.
    :type bold_labels: bool

    :param useMathTextOption: option for the useMathText option for the y axis.
    :type useMathTextOption: bool

    :param generate_title: option to enable (True) or disable (False) figure titles.
        If no title is included in the subplot dictionary, a title is created. The
        created title will include location numbers if the variable name includes
        a location (e.g., '001' in 'msw1_001.C02_aquifer1').
    :type generate_title: bool

    :param plot_grid_option: option to display a grid (True) or not (False)
    :type plot_grid_option: bool

    :param grid_alpha_val: alpha value for the grid
    :type grid_alpha_val: float

    :return: None
    """
    # Dictionary with variables used to adjust figure formatting
    fig_setup = {'gen_font_size': gen_font_size,
             'axis_label_font_size': axis_label_font_size,
             'title_font_size': title_font_size,
             'legend_font_size': legend_font_size,
             'line_width': LINEWIDTH_REF,
             'ax_line_width': AX_LINEWIDTH_REF,
             'xtick_size': XTICK_SIZE_REF,
             'ytick_size': YTICK_SIZE_REF,
             'xtick_width': XTICK_WIDTH_REF,
             'ytick_width': YTICK_WIDTH_REF,
             'axis_label_pad': AXIS_LABEL_PAD_REF,
             'title_pad': TITLE_PAD_REF,
             'xaxis_font_size': axis_label_font_size,
             'yaxis_font_size': axis_label_font_size}

    selected_keys = ['gen_font_size', 'axis_label_font_size',
                     'title_font_size', 'legend_font_size', 'line_width',
                     'ax_line_width', 'xtick_size', 'ytick_size',
                     'xtick_width', 'ytick_width']

    if bold_labels:
        fig_setup['label_font_weight'] = 'bold'
    else:
        fig_setup['label_font_weight'] = 'normal'

    if figsize is None:
        figsize = (DEFAULT_FIG_WIDTH, DEFAULT_FIG_HEIGHT)
    else:
        # Scale font sizes to the specified figure size. Here, the updated
        # font size scales with the new height or width (relative to the default
        # height or width). This scaling uses the length scale with the largest change
        # relative to the default values (e.g., if height has a larger change, then
        # fontsize is scaled using the specified height value).
        dw_ratio = np.abs(figsize[0] - DEFAULT_FIG_WIDTH)/DEFAULT_FIG_WIDTH
        dh_ratio = np.abs(figsize[1] - DEFAULT_FIG_HEIGHT)/DEFAULT_FIG_HEIGHT

        if figsize != (DEFAULT_FIG_WIDTH, DEFAULT_FIG_HEIGHT):
            if dw_ratio > dh_ratio:
                L1 = DEFAULT_FIG_WIDTH
                L2 = figsize[0]
            elif dw_ratio < dh_ratio:
                L1 = DEFAULT_FIG_HEIGHT
                L2 = figsize[1]

            # Update the font sizes -  some of these are also updated later on, depending
            # on subplot sizes. The fontsize is adjusted depending on the input figsize -
            # sufficiently shrinking the fontsize can be important for small figures, but
            # increasing font sizes in the same way often leads to font that is way too large.
            if (L2 / L1) < 1:
                for key in selected_keys:
                    # The formula is simplified based on initial Nate's idea
                    fig_setup[key] = 0.5*(1 + L2/L1) * fig_setup[key]
            elif (L2 / L1) > 1:
                for key in selected_keys:
                    # The formula is simplified based on initial Nate's idea
                    fig_setup[key] = 0.25*(3 + L2/L1) * fig_setup[key]

    # These are updated separately depending on the number of rows and columns
    xaxis_font_size_ref = fig_setup['axis_label_font_size']
    yaxis_font_size_ref = fig_setup['axis_label_font_size']
    title_font_size_ref = fig_setup['title_font_size']
    legend_font_size_ref = fig_setup['legend_font_size']

    # Update matplotlib figure setup
    update_rc_setup(fig_setup, fontname)

    # Find number of subplots
    num_plots = get_number_of_subplots(output_names, output_list)

    # Process subplots and their type data
    # subplots_data has the following keys: 'use', 'nrows', 'ncols',
    # 'single_plot', 'type' + (possibly) keys representing titles of subplots
    # if defined by user
    subplots_data = setup_subplots_data(subplot, plot_type, num_plots)

    # Initialize indices
    subplot_ind = 1
    color_ind = 0
    reals_ind = 0
    used_colors = []
    markerRef = 0
    lineStyleRef = 0
    hatchRef = 0

    useMarkers, useLines, varyLineStyles, figdpi = get_plot_yaml_input(
        plot_data, name)

    if not useMarkers:
        markerList = ['None']
    else:
        markerList = defaultMarkerList

    if not useLines:
        lineStyleList = ['None']
        varyLineStyles = False
    else:
        lineStyleList = defaultLineStyleList

    # Transform time points from days to years
    time = sm.time_array / 365.25
    # Create figure
    fig = plt.figure(num=name, figsize=figsize, dpi=figdpi)
    # Loop over observation names in outputs
    for obs_to_plot in output_names:

        # Some of the sizes can be adjusted within the loop, so reset them each time
        fig_setup = reset_sizes(fig_setup, xaxis_font_size_ref, yaxis_font_size_ref,
                                title_font_size_ref, legend_font_size_ref)
        xaxis_fontsizes_used = []
        yaxis_fontsizes_used = []

        # If this figure has only one plot, slightly increase the font sizes
        # (no risk of overlap with other subplots)
        if subplots_data['single_plot']:
            fig_setup = update_single_plot_setup(fig_setup, fontname)

        # List of components providing given observation
        cmpnts_to_process = process_components(obs_to_plot, output_list)

        for cmpnt_obj in cmpnts_to_process:

            cmpnt_name = cmpnt_obj.name
            full_obs_nm = '.'.join([cmpnt_name, obs_to_plot])

            # Add subplot
            ax = plt.subplot(
                subplots_data['nrows'], subplots_data['ncols'], subplot_ind)

            if not subplots_data['single_plot']:
                subplot_ind += 1
                color_ind = 0
                reals_ind = 0
                used_colors = []
                markerRef = 0
                lineStyleRef = 0

            lgnd_label, label, loc_ind = generate_legend_setup(
                obs_to_plot, cmpnt_name, analysis)

            colorValReal, colorValStats, used_colors, color_ind, \
                hatch_check = get_colors(reals_ind, color_ind, used_colors,
                                         subplots_data)

            if useMarkers and not colorValReal is None:
                rgbReal = clrs.to_rgba(colorValReal[:])
                markerEdgeColor = np.array(list(rgbReal[:]))
                markerEdgeColor /= 2
                markerEdgeColor[-1] = 1
            else:
                markerEdgeColor = 'None'

            if analysis == 'forward':
                values = sm.collect_observations_as_time_series(
                    cmpnt_obj, obs_to_plot)
                ax.plot(time, values, '-', label=label, color=colorValReal,
                        marker=markerList[markerRef],
                        markeredgecolor = markerEdgeColor,
                        linestyle=lineStyleList[lineStyleRef], alpha=0.8,
                        linewidth=fig_setup['line_width'])
                reals_ind = reals_ind + 1

            elif analysis in ['lhs', 'parstudy']:
                ind_list = list(range(len(time)))
                obs_names = [full_obs_nm + '_{0}'.format(indd)
                             for indd in ind_list]
                obs_percentiles = percentile(s.recarray[obs_names],
                                             [0, 25, 50, 75, 100])
                obs_means = mean(s.recarray[obs_names])
                values = np.array(
                    [s.recarray[full_obs_nm + '_'
                                + str(indd)] for indd in ind_list])

                if 'real' in subplots_data['type']:
                    if 'stats' in subplots_data['type']:
                        ax.plot(time, values, color=colorValReal,
                                marker=markerList[markerRef],
                                markeredgecolor = markerEdgeColor,
                                linestyle=lineStyleList[lineStyleRef],
                                label=label, linewidth=fig_setup['line_width'],
                                alpha=0.33, zorder = 0)
                    else:
                        ax.plot(time, values, color=colorValReal,
                                marker=markerList[markerRef],
                                markeredgecolor = markerEdgeColor,
                                linestyle=lineStyleList[lineStyleRef],
                                label=label, linewidth=fig_setup['line_width'],
                                alpha=0.8)
                    reals_ind = reals_ind + 1

                if 'stats' in subplots_data['type']:
                    if not hatch_check:
                        ax.fill_between(time, obs_percentiles[3, :],
                                        obs_percentiles[4, :],
                                        label='Upper and lower quartiles' + lgnd_label,
                                        color=colorValStats, alpha=0.15)
                        ax.fill_between(time, obs_percentiles[1, :],
                                        obs_percentiles[3, :],
                                        label='Middle quartiles' + lgnd_label,
                                        color=colorValStats, alpha=0.3)
                        ax.fill_between(time, obs_percentiles[0, :],
                                        obs_percentiles[1, :],
                                        color=colorValStats, alpha=0.15)
                        ax.plot(time, obs_percentiles[2, :], color=colorValStats,
                                label='Median value' + lgnd_label,
                                linewidth = fig_setup['line_width'],
                                linestyle = ':', alpha=0.8)
                        ax.plot(time, obs_means, color=colorValStats,
                                label='Mean value' + lgnd_label,
                                linewidth = fig_setup['line_width'], alpha=0.8)

                    elif hatch_check:
                        ax.fill_between(time, obs_percentiles[3, :],
                                        obs_percentiles[4, :],
                                        label='Upper and lower quartiles' + lgnd_label,
                                        color=colorValStats, alpha=0.15,
                                        hatch = hatchList[hatchRef])
                        ax.fill_between(time, obs_percentiles[1, :],
                                        obs_percentiles[3, :],
                                        label='Middle quartiles' + lgnd_label,
                                        color=colorValStats, alpha=0.3,
                                        hatch = hatchList[hatchRef])
                        ax.fill_between(time, obs_percentiles[0, :],
                                        obs_percentiles[1, :],
                                        color=colorValStats, alpha=0.15,
                                        hatch = hatchList[hatchRef])
                        ax.plot(time, obs_percentiles[2, :], color=colorValStats,
                                label='Median value' + lgnd_label,
                                linewidth = fig_setup['line_width'],
                                linestyle = ':', alpha=0.8)
                        ax.plot(time, obs_means, color=colorValStats,
                                label='Mean value' + lgnd_label,
                                linewidth = fig_setup['line_width'], alpha=0.8)
                        hatch_check = False
                        hatchRef += 1
                        if hatchRef > (len(hatchList) - 1):
                            hatchRef = 0

                    color_ind += 1

            if useMarkers:
                markerRef += 1
                if markerRef > (len(defaultMarkerList) - 1):
                    markerRef = 0

            if varyLineStyles:
                lineStyleRef += 1
                if lineStyleRef > (len(defaultLineStyleList) - 1):
                    lineStyleRef = 0

            if plot_grid_option:
                ax.grid(alpha=grid_alpha_val)

            # X LABEL
            fig_setup = adjust_x_label(ax, fig, fig_setup, subplots_data)
            xaxis_fontsizes_used.append(fig_setup['xaxis_font_size'])

            # Y LABEL
            fig_setup = adjust_y_label(obs_to_plot, cmpnt_name, ax, fig,
                                        fig_setup, subplots_data, useMathTextOption)
            yaxis_fontsizes_used.append(fig_setup['yaxis_font_size'])

            # TITLE
            sub_title = get_title(obs_to_plot, cmpnt_name, subplots_data, loc_ind)

            if generate_title:
                fig_setup = adjust_title(sub_title, ax, fig, fig_setup, subplots_data)
            else:
                # No title
                pass

    min_fontsize = make_label_fontsizes_uniform(
        subplot_ind - 1, xaxis_fontsizes_used, yaxis_fontsizes_used,
        subplots_data)

    if fig_setup['gen_font_size'] > min_fontsize:
        fig_setup['gen_font_size'] = min_fontsize
        update_rc_setup(fig_setup, fontname)

    # Used to reset the legend font size within the loop through the axes.
    # Otherwise, the font size could shrink from one subplot to the next.
    legend_font_size_ref2 = fig_setup['legend_font_size']

    if analysis == 'forward':
        ax_list = fig.axes
        for ax in ax_list:
            handle_list = []
            label_list = []
            handles, labels = ax.get_legend_handles_labels()

            for handle, label in zip(handles, labels):
                if label not in label_list:
                    handle_list.append(handle)
                    label_list.append(label)

            fig_setup['legend_font_size'] = legend_font_size_ref2
            fig_setup = check_legend(handle_list, fig_setup, min_fontsize,
                                     subplots_data)

            ax.legend(handle_list, label_list, fancybox=False,
                      fontsize=fig_setup['legend_font_size'],
                      framealpha=fig_setup['legend_framealpha'],
                      ncol=fig_setup['legend_columns'])

    elif analysis in ['lhs', 'parstudy']:
        ax_list = fig.axes
        for ax in ax_list:
            handle_list = []
            label_list = []
            handles, labels = ax.get_legend_handles_labels()
            for handle, label in zip(handles, labels):
                if label not in label_list:
                    handle_list.append(handle)
                    label_list.append(label)

            fig_setup['legend_font_size'] = legend_font_size_ref2
            fig_setup = check_legend(handle_list, fig_setup, min_fontsize,
                                     subplots_data)

            ax.legend(handle_list, label_list, fancybox=False,
                      fontsize=fig_setup['legend_font_size'],
                      framealpha=fig_setup['legend_framealpha'],
                      ncol=fig_setup['legend_columns'])
    else:
        pass

    if title:
        fig.suptitle(title, fontweight=fig_setup['label_font_weight'],
                     fontsize=fig_setup['title_font_size'])

    # With 3 or more rows, the titles and x-axis labels often overlap.
    if subplots_data['single_plot']:
        fig.subplots_adjust(left=0.15, bottom=0.1, right=0.85, top=0.9,
                            wspace=0.1, hspace=0.1)
    elif subplots_data['nrows'] >= 3:
        fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9,
                            wspace=0.25, hspace=0.5)
    else:
        fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9,
                            wspace=0.25, hspace=0.33)

    if savefig:
        try:
            fig.savefig(savefig)
        except ValueError:
            # User has specified plot with a '.' in name but no extension.
            # Add .png as output format.
            savefig += '.png'
            fig.savefig(savefig)
    else:
        fig.show()

    plt.close()
    # Restore to default matplotlib settings
    # matplotlib.rcdefaults()


def reset_sizes(fig_setup, xaxis_font_size_ref, yaxis_font_size_ref,
                title_font_size_ref, legend_font_size_ref):
    """
    Reset some of the sizes to original values

    Parameters
    ----------
    fig_setup : dict
        Contains information about sizes of different figure elements.

    Returns updated dictionary sizes
    -------
    """
    fig_setup['xaxis_font_size'] = xaxis_font_size_ref
    fig_setup['yaxis_font_size'] = yaxis_font_size_ref
    fig_setup['title_font_size'] = title_font_size_ref
    fig_setup['legend_font_size'] = legend_font_size_ref
    fig_setup['legend_framealpha'] = LEGEND_FRAME_ALPHA
    fig_setup['legend_columns'] = LEGEND_COLUMNS
    fig_setup['axis_label_pad'] = AXIS_LABEL_PAD_REF
    fig_setup['title_pad'] = TITLE_PAD_REF

    return fig_setup


def update_single_plot_setup(fig_setup, fontname):
    """
    Update some of the relevant fig_setup entries if single plot (without subplots)
    is to be generated.

    Parameters
    ----------
    fig_setup : dict
        Contains information about sizes and properties of different figure elements.

    Returns updated dictionary fig_setup
    """
    fig_setup['gen_font_size'] *= SINGLE_PLOT_FONTSIZE_ADJUST
    fig_setup['xaxis_font_size'] *= SINGLE_PLOT_FONTSIZE_ADJUST
    fig_setup['yaxis_font_size'] *= SINGLE_PLOT_FONTSIZE_ADJUST
    fig_setup['title_font_size'] *= SINGLE_PLOT_FONTSIZE_ADJUST
    # The legend font sizes can easily become too big, so the initial
    # adjustment is scaled down
    fig_setup['legend_font_size'] = 0.5*fig_setup['legend_font_size']*(
        1 + SINGLE_PLOT_FONTSIZE_ADJUST)
    value = fig_setup['gen_font_size']
    font = {'family': fontname,
            'weight': 'normal',
            'size': value}
    plt.rc('font', **font)
    fig_setup['axis_label_pad'] *= SINGLE_PLOT_FONTSIZE_ADJUST
    fig_setup['title_pad'] *= SINGLE_PLOT_FONTSIZE_ADJUST

    return fig_setup


def update_rc_setup(fig_setup, fontname):
    """
    Update matplotlib.pyplot parameters using information in the fig_setup and fontname
    arguments.

    Parameters
    ----------
    fig_setup : dict
        Contains information about sizes and properties of different figure elements.

    fontname : str
        Name of font for figure elements

    Returns
    -------
    None.
    """
    font = {'family': fontname,
            'weight': 'normal',
            'size': fig_setup['gen_font_size']}
    plt.rc('font', **font)
    plt.rcParams['axes.linewidth'] = fig_setup['ax_line_width']
    plt.rcParams['xtick.major.size'] = fig_setup['xtick_size']
    plt.rcParams['ytick.major.size'] = fig_setup['ytick_size']
    plt.rcParams['xtick.major.width'] = fig_setup['xtick_width']
    plt.rcParams['ytick.major.width'] = fig_setup['ytick_width']


def get_number_of_subplots(output_names, cmpnt_to_output):
    """
    Calculate required number of subplots for the requested outputs.

    :param output_names: List of observation names to match with component models.
    Examples:
        output_names=['pressure']

        output_names=['pressure', 'mass_CO2_aquifer1']

        output_names=['pressure', 'CO2_aquifer1', 'CO2_aquifer2']
    :type output_names: list

    :param cmpnt_to_output: Dictionary mapping component models to observations
    Examples:
        This dictionary only includes pressure from a component saved in variable sres.
        cmpnt_to_output={sres: 'pressure'}

    Returns number of subplots required
    """
    nplots = 0
    for obs_nm in output_names:
        for cmpnt in cmpnt_to_output:
            if obs_nm in cmpnt_to_output[cmpnt]:
                nplots += 1

    return nplots


def setup_subplots_data(subplot, plot_type, num_plots):
    """
    Setup dictionary containing information about subplots:
        single subplot versus many, number of rows and columns,
        type of data (realizations and/or stats).

    Possible keys: 'use', 'nrows', 'ncols', 'single_plot', 'type' (plus data about
    titles corresponding to different cmpnt_name.obs_name subplots; these data
    might not be necessarily present). Note that this function is also set up to
    read capitalized versions of the in .yaml files (Use, and NumCols instead
    of use and ncols), as this approach matches the conventions of other plot
    types. The capitalized or non-capitaized inputs will have the same effects,
    however.
    """
    # Initialize subplots_data dictionary depending on the input arguments
    if subplot is None:
        subplots_data = {'use': False}
    else:
        subplots_data = subplot

    if 'Use' in subplots_data:
        subplots_data['use'] = subplots_data['Use']
        del subplots_data['Use']

    if 'NumCols' in subplots_data:
        subplots_data['ncols'] = subplots_data['NumCols']
        del subplots_data['NumCols']

    # Process plot type
    if plot_type is None:
        subplots_data['type'] = ['real']
    else:
        subplots_data['type'] = plot_type

    if not subplots_data['use']:  # if subplots are not to be used
        subplots_data['single_plot'] = True
        subplots_data['ncols'] = 1
        subplots_data['nrows'] = 1
    else:
        subplots_data['single_plot'] = False
        if 'ncols' not in subplots_data:
            if num_plots <= 3:
                subplots_data['ncols'] = 1
            elif num_plots > 3:
                subplots_data['ncols'] = 2

        subplots_data['nrows'] = math.ceil(float(num_plots) / subplots_data['ncols'])

        # If the user entered 'use': True in the subplot dictionary but there is
        # still only one row and one column in this plot, set single_plot to True.
        if subplots_data['ncols'] == 1 and subplots_data['nrows'] == 1:
            subplots_data['single_plot'] = True

    return subplots_data


def get_label(obs_name, labels, split_labels):
    """
    Get label for y-axis for a plot of a given observation

    :param obs_name: name of observation to be shown on a figure
    :type obs_name: str

    :param labels: dictionary containing y-labels corresponding to observation
    names
    :type labels: dict()

    :param split_labels: dictionary containing y-labels corresponding to particular
    parts of observation names
    :type split_labels: dict()

    Returns:
        string to use for y-label
    """
    out_flag = 1
    if obs_name in labels:
        str_label = labels[obs_name]
    else:
        # The observation name is split on numerical characters:
        # it can be 2 in CO2_aquifer or it can be 1 (or similar) in brine_aquifer1
        # This works on preexisting known observation names
        split_name = re.split('\d', obs_name)
        try:
            str_label = split_labels[split_name[0]]
        except KeyError:
            str_label = obs_name
            out_flag = 0

    return str_label, out_flag


def generate_legend_setup(obs_name, cmpnt_name, analysis):
    """
    Generate legend setup

    Returns lgnd_label, label and line_style
    """
    if '_' in cmpnt_name:
        cmpnt_name_edit = cmpnt_name[0:cmpnt_name.index('_')]
    else:
        cmpnt_name_edit = cmpnt_name

    # Determine whether location is provided in the component name and get index if it is
    loc_ind = is_location_in_name(cmpnt_name)

    # Get initial version of legend label and update it later if location index
    # is in the component name
    lgnd_label = get_legend_label(obs_name)

    # If location index is known
    if loc_ind != -1:
        if analysis == 'forward':
            if lgnd_label != '':
                lgnd_label =  '{} to {} at location {}'.format(
                    cmpnt_name_edit, lgnd_label, loc_ind)
            else:
                lgnd_label = '{} at location {}'.format(cmpnt_name_edit, loc_ind)
        else:
            # If displaying lhs results, the lines are displayed as 'Simulated values'
            # (when lots of lines are shown), 'Median value', or 'Mean value.' In the
            # latter two cases, the ' at location #" is added for clarification.
            if lgnd_label != '':
                # lgnd_label can be 'Aquifer 1' or 'Aquifer 2'
                lgnd_label = ', {} to {} at location {}'.format(
                    cmpnt_name_edit, lgnd_label, loc_ind)
            else:
                lgnd_label = ', {} at location {}'.format(cmpnt_name_edit, loc_ind)
    else:
        lgnd_label = cmpnt_name_edit

    if analysis == 'forward':
        label = lgnd_label
    elif analysis in ['lhs', 'parstudy']:
        if lgnd_label == cmpnt_name_edit:
            lgnd_label = ', ' + cmpnt_name_edit

        label = 'Simulated values'+lgnd_label

    return lgnd_label, label, loc_ind


def get_legend_label(obs_name):
    """
    Generate part of the legend based on the provided observation name
    """
    if obs_name in LEGEND_DICT:
        legend_label = LEGEND_DICT[obs_name]
    else:
        # Check if 'aquifer' in the name
        place_ind = obs_name.rfind('aquifer')
        if place_ind != -1:
            # Determine index of aquifer
            try:
                aquifer_ind = int(obs_name[place_ind+7:])
            except ValueError:
                # Possibly observations of Seal Horizon or Fault Flow components
                # TODO update with location being defined by cell or segment ind
                legend_label = ''
            else:
                legend_label = 'aquifer {}'.format(aquifer_ind)
        else:
            legend_label = ''

    return legend_label


def adjust_x_label(ax, fig, fig_setup, subplots_data):
    """
    Adjust font of x-label based on figure size.
    """
    h_xlabel = ax.set_xlabel(
        'Time, t [years]', fontsize=fig_setup['xaxis_font_size'],
        fontweight=fig_setup['label_font_weight'],
        labelpad=fig_setup['axis_label_pad'])

    continue_test = True
    while continue_test:
        height_frac, width_frac = width_and_depth_frac(fig, h_xlabel, subplots_data)

        # If the xlabel is too long, shrink the fontsize
        if (width_frac > MAX_XLABEL_WIDTH_FRAC) or (height_frac > MAX_XLABEL_HEIGHT_FRAC):
            fig_setup['xaxis_font_size'] = 0.95*fig_setup['xaxis_font_size']
            h_xlabel = ax.set_xlabel(
                'Time, t [years]', fontsize=fig_setup['xaxis_font_size'],
                fontweight=fig_setup['label_font_weight'],
                labelpad=fig_setup['axis_label_pad'])
        else:
            continue_test = False

    return fig_setup


def adjust_y_label(obs_name, cmpnt_name, ax, fig, fig_setup, subplots_data, useMathTextOption):
    """
    Adjust font of y-label based on figure size.
    """
    # Get y-label associated with given observation
    y_label, out_flag = get_label(obs_name, Y_LABEL_DICT, Y_LABEL_SPLIT_DICT)
    if out_flag != 1:
        y_label = '{}.{}'.format(cmpnt_name, obs_name)

    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(0, 0), useMathText=useMathTextOption)

    h_ylabel = ax.set_ylabel(y_label, fontsize=fig_setup['yaxis_font_size'],
                              fontweight=fig_setup['label_font_weight'],
                              labelpad=fig_setup['axis_label_pad'])

    if out_flag == 1:
        height_frac, _ = width_and_depth_frac(fig, h_ylabel, subplots_data)
        # If the ylabel is too long relative to the figure, use labels with two rows
        if height_frac > MAX_YLABEL_HEIGHT_FRAC:
            y_label, _ = get_label(
                obs_name, Y_LABEL_2ROWS_DICT, Y_LABEL_2ROWS_SPLIT_DICT)
            h_ylabel = ax.set_ylabel(y_label,
                                     fontsize=fig_setup['yaxis_font_size'],
                                     fontweight=fig_setup['label_font_weight'],
                                     labelpad=fig_setup['axis_label_pad'])

    # Check if the fontsize is too large
    continue_test = True
    while continue_test:
        height_frac, width_frac = width_and_depth_frac(fig, h_ylabel, subplots_data)
        # If the ylabels are still too long, shrink the fontsize
        if (height_frac > MAX_YLABEL_HEIGHT_FRAC) or (width_frac > MAX_YLABEL_WIDTH_FRAC):
            fig_setup['yaxis_font_size'] *= 0.95
            h_ylabel = ax.set_ylabel(y_label,
                                     fontsize=fig_setup['yaxis_font_size'],
                                     fontweight=fig_setup['label_font_weight'],
                                     labelpad=fig_setup['axis_label_pad'])
        else:
            continue_test = False

    return fig_setup


def get_title(obs_name, cmpnt_name, subplots_data, loc_ind):
    """
    Generate figure title based on the provided observation name
    """
    full_obs_name = '{}.{}'.format(cmpnt_name, obs_name)
    if full_obs_name in subplots_data:
        sub_title = subplots_data[full_obs_name]

    # This checks if the name includes a number like "_000," which indicates a location
    else:
        title_label, _ = get_label(obs_name, TITLE_DICT, TITLE_SPLIT_DICT)
        if loc_ind != -1:  # -1 means no location in component name
            # If it's a single plot, the results plotted could represent multiple locations.
            # In that scenario, you shouldn't have one location displayed in the title.
            if not subplots_data['single_plot']:
                sub_title = '{} at Location {}'.format(title_label, loc_ind)
            else:
                sub_title = title_label
        # If the location number is not in the name, just use the title dictionary
        else:
            sub_title = title_label

    return sub_title


def adjust_title(sub_title, ax, fig, fig_setup, subplots_data):
    """
    Adjust font of subplot title based on figure size.
    """
    h_title = ax.set_title(sub_title, fontsize=fig_setup['title_font_size'],
                           fontweight=fig_setup['label_font_weight'],
                           pad=fig_setup['title_pad'])

    height_frac, width_frac = width_and_depth_frac(fig, h_title, subplots_data)
    # If the title is too long, it can overlap the axis labels (e.g., 'x 10^6')
    # If there is a space in the title, find it and make a line break
    if (height_frac > MAX_TITLE_HEIGHT_FRAC) or (width_frac > MAX_TITLE_WIDTH_FRAC):
        if ' ' in sub_title:
            sub_title = split_at_space(sub_title)
            h_title = ax.set_title(sub_title, fontsize=fig_setup['title_font_size'],
                                   fontweight=fig_setup['label_font_weight'],
                                   pad=fig_setup['title_pad'])
        else:  # If there's no space in the title shrink the fontsize
            fig_setup['title_font_size'] *= 0.95
            h_title = ax.set_title(sub_title,
                                   fontsize=fig_setup['title_font_size'],
                                   fontweight=fig_setup['label_font_weight'],
                                   pad=fig_setup['title_pad'])

        continue_test = True
        while continue_test:
            height_frac, width_frac = width_and_depth_frac(fig, h_title, subplots_data)

            # If the title is still too long, shrink the fontsize
            if (height_frac > MAX_TITLE_HEIGHT_FRAC) or (width_frac > MAX_TITLE_WIDTH_FRAC):
                fig_setup['title_font_size'] *= 0.95
                h_title = ax.set_title(sub_title,
                                       fontsize=fig_setup['title_font_size'],
                                       fontweight=fig_setup['label_font_weight'],
                                       pad=fig_setup['title_pad'])
            else:
                continue_test = False

    return fig_setup


def is_location_in_name(cmpnt_name):
    """
    Determine whether name of component contains location: applicable only for
    control file interface created components with names in the format name_###

    Returns location index extracted from name if found and -1 if there is
    no location index in the component name.
    """
    # Determine index of last underscore in the name: if -1 is returned underscore
    # symbol is not present
    underscore_ind = cmpnt_name.rfind('_')

    # Default value of output that can be changed to location specified
    # the component name
    loc_ind = -1

    if underscore_ind != -1: # i.e., underscore is found
        # Check that all symbols after underscore are numeric
        if np.char.isnumeric(cmpnt_name[underscore_ind+1:]):
            # Transform what is after underscore symbol into location index
            loc_ind = int(cmpnt_name[underscore_ind+1:])

    return loc_ind


def split_at_space(label):
    """
    Split string at '_' (space symbol) trying for the two parts
    to be approximately the same.
    """
    # Find a space (' ') near the middle of the title
    center_of_label_index = int(np.ceil(len(label) / 2))

    if label[center_of_label_index] == ' ':
        final_label = '{}\n{}'.format(label[0:center_of_label_index],
                                      label[center_of_label_index + 1:])
    else:
        lower_index = center_of_label_index - 1
        upper_index = center_of_label_index + 1
        continue_test = True
        while continue_test:
            if label[lower_index] == ' ':
                final_label = '{}\n{}'.format(label[0:lower_index],
                                                   label[lower_index + 1:])
                continue_test = False
            elif label[upper_index] == ' ':
                final_label = '{}\n{}'.format(label[0:upper_index],
                                        label[upper_index + 1:])
                continue_test = False
            else:
                lower_index -= 1
                upper_index += 1

    return final_label


def width_and_depth_frac(fig, element, subplots_data):
    """
    Calculate width and height fractions for a given figure.
    """
    fig_renderer = fig.canvas.get_renderer()
    bb = element.get_window_extent(renderer=fig_renderer)
    ywidth = bb.width
    yheight = bb.height

    height_frac = yheight / (fig_renderer.height / subplots_data['nrows'])
    width_frac = ywidth / (fig_renderer.width / subplots_data['ncols'])

    return height_frac, width_frac


def process_components(obs_name, cmpnt_to_output):
    """
    Return list of components from the dictionary returning given observation.
    """
    comp_list = []
    for cmpnt in cmpnt_to_output:
        if obs_name in cmpnt_to_output[cmpnt]:
            comp_list.append(cmpnt)

    return comp_list


def make_label_fontsizes_uniform(num_subplots, xaxis_fontsizes_used,
                                 yaxis_fontsizes_used, subplots_data):
    """
    Function that ensures all x and y axis labels have the same fontsizes
    """

    min_fontsize = np.min(xaxis_fontsizes_used)

    if np.min(yaxis_fontsizes_used) < min_fontsize:
        min_fontsize = np.min(yaxis_fontsizes_used)

    if num_subplots == 0:
        ax = plt.gca()
        ax.xaxis.label.set_fontsize(min_fontsize)
        ax.yaxis.label.set_fontsize(min_fontsize)

    else:
        for subplotRef in range(0, num_subplots):
            ax = plt.subplot(
                subplots_data['nrows'], subplots_data['ncols'], subplotRef + 1)

            ax = plt.gca()
            ax.xaxis.label.set_fontsize(min_fontsize)
            ax.yaxis.label.set_fontsize(min_fontsize)

    return min_fontsize


def get_colors(reals_ind, color_ind, used_colors, subplots_data):
    """
    Function that checks the colors used previously and provides a color that
    has not been used yet.
    """
    hatch_check = False

    if 'real' in subplots_data['type']:
        colorRefReal = 'C' + str((reals_ind) % 10)
        colorValReal = colorRefReal[:]

        rgbReal = clrs.to_rgba('C' + str((reals_ind) % 10))
        darkColorRefReal = 'DC' + str((reals_ind) % 10)
        lightColorRefReal = 'LC' + str((reals_ind) % 10)

        if not colorRefReal in used_colors:
            used_colors.append(colorRefReal)

        else:
            if not darkColorRefReal in used_colors:
                dark_clr = np.array(list(rgbReal[:]))
                dark_clr *= (2 / 3)
                dark_clr[-1] = 1
                colorValReal = dark_clr

                used_colors.append(darkColorRefReal)

            elif not lightColorRefReal in used_colors:
                light_clr = np.array(list(rgbReal[:]))
                light_clr[0] = (light_clr[0] + 1) / 2
                light_clr[1] = (light_clr[1] + 1) / 2
                light_clr[2] = (light_clr[2] + 1) / 2
                light_clr[-1] = 1
                colorValReal = light_clr

                used_colors.append(lightColorRefReal)
    else:
        colorValReal = None

    if 'stats' in subplots_data['type']:
        colorRefStats = 'C' + str((color_ind) % 10)
        colorValStats = colorRefStats[:]

        rgbStats = clrs.to_rgba('C' + str((color_ind) % 10))
        darkColorRefStats = 'DC' + str((color_ind) % 10)
        lightColorRefStats = 'LC' + str((color_ind) % 10)

        if not colorRefStats in used_colors:
            used_colors.append(colorRefStats)

        else:
            color_ind += 1

            colorRefStats = 'C' + str((color_ind) % 10)
            colorValStats = colorRefStats[:]

            rgbStats = clrs.to_rgba('C' + str((color_ind) % 10))
            darkColorRefStats = 'DC' + str((color_ind) % 10)
            lightColorRefStats = 'LC' + str((color_ind) % 10)

            if not colorRefStats in used_colors:
                used_colors.append(colorRefStats)

            else:
                if not darkColorRefStats in used_colors:
                    dark_clr = np.array(list(rgbStats[:]))
                    dark_clr /= 2
                    dark_clr[-1] = 1
                    colorValStats = dark_clr

                    used_colors.append(darkColorRefStats)

                elif not lightColorRefStats in used_colors:
                    light_clr = np.array(list(rgbReal[:]))
                    light_clr[0] = (light_clr[0] + 1) / 2
                    light_clr[1] = (light_clr[1] + 1) / 2
                    light_clr[2] = (light_clr[2] + 1) / 2
                    light_clr[-1] = 1
                    colorValStats = light_clr

                    used_colors.append(lightColorRefStats)

                else:
                    # If there are no more colors to use, the hatches can help
                    # distinguish the filled areas.
                    hatch_check = True
    else:
        colorValStats = None

    return colorValReal, colorValStats, used_colors, color_ind, hatch_check


def check_legend(handle_list, fig_setup, min_fontsize, subplots_data):
    """
    This function checks the number of items in the legend (handle_list) and
    adjusts the legend fontsize and columns if there are too many items. The
    min_fontsize is the minimum fontsize used for x and y axis labels, which is
    adjusted based on the number of subplots. If the legend fontsize is larger
    than min_fontsize, it is set to min_fontsize.
    """

    if fig_setup['legend_font_size'] > min_fontsize:
        fig_setup['legend_font_size'] = min_fontsize

    if subplots_data['ncols'] == 2:
        fig_setup['legend_font_size'] *= 0.9
    elif subplots_data['ncols'] == 3:
        fig_setup['legend_font_size'] *= 0.75
    elif subplots_data['ncols'] >= 4:
        fig_setup['legend_font_size'] *= 0.67

    if subplots_data['single_plot']:
        if LEGEND_ITEM_THRESH1 <= len(handle_list) < LEGEND_ITEM_THRESH2:
            fig_setup['legend_font_size'] *= 0.9

        elif LEGEND_ITEM_THRESH2 <= len(handle_list) < LEGEND_ITEM_THRESH3:
            fig_setup['legend_font_size'] *= 0.75
            fig_setup['legend_framealpha'] *= 0.75

        elif LEGEND_ITEM_THRESH3 <= len(handle_list) < LEGEND_ITEM_THRESH4:
            fig_setup['legend_font_size'] *= 0.67
            fig_setup['legend_framealpha'] *= 0.67
            fig_setup['legend_columns'] = 2

        elif len(handle_list) >= LEGEND_ITEM_THRESH4:
            fig_setup['legend_font_size'] *= 0.5
            fig_setup['legend_framealpha'] *= 0.5
            fig_setup['legend_columns'] = 2

    else:
        if LEGEND_ITEM_THRESH1 <= len(handle_list) < LEGEND_ITEM_THRESH2:
            fig_setup['legend_font_size'] *= 0.75

        elif LEGEND_ITEM_THRESH2 <= len(handle_list) < LEGEND_ITEM_THRESH3:
            fig_setup['legend_font_size'] *= 0.67
            fig_setup['legend_framealpha'] *= 0.75

        elif LEGEND_ITEM_THRESH3 <= len(handle_list) < LEGEND_ITEM_THRESH4:
            fig_setup['legend_font_size'] *= 0.5
            fig_setup['legend_framealpha'] *= 0.5
            fig_setup['legend_columns'] = 2

        elif len(handle_list) >= LEGEND_ITEM_THRESH4:
            fig_setup['legend_font_size'] *= 0.33
            fig_setup['legend_framealpha'] *= 0.5
            fig_setup['legend_columns'] = 2

    return fig_setup

def get_plot_yaml_input(plot_data, name):
    """
    This function checks the plot's section within the .yaml file for any input
    related to markerstyles and linestyles.
    """
    # Each value under a particular key is a list: default value, type, type written as string
    default_values = {'UseMarkers': [False, bool, 'boolean'],
                      'UseLines': [True, bool, 'boolean'],
                      'VaryLineStyles': [False, bool, 'boolean'],
                      'FigureDPI': [100, (int, float), 'integer or float']}

    out_values = {key: val[0] for key, val in default_values.items()}

    warning_msg = ''.join([
        'In the .yaml file, the input provided for {} under the ', name, ' plot ',
        'was not of type {}. The default value of {} will be used.'])

    for key, val in default_values.items():
        if key in plot_data:
            out_values[key] = plot_data[key]
            if not isinstance(out_values[key], val[1]):
                msg = warning_msg.format(key, val[2], val[0])
                logging.warning(msg)
                out_values[key] = val[0]

    return out_values['UseMarkers'], out_values['UseLines'],\
        out_values['VaryLineStyles'], out_values['FigureDPI']
