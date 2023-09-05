# -*- coding: utf-8 -*-
"""
Code to create map-view figure of plume areas due to CO2 leaking to the atmosphere.

Examples illustrating applications or setup of atmosphere_plot method:
    atmosphere_plot.py (see at the end of this file)
    ControlFile_ex9a.yaml
    ControlFile_ex9b.yaml
    ControlFile_ex9c.yaml

Created on Thu Jan 25 12:16:47 2018
Last Modified: February 10th, 2023

@author: Seth King
AECOM supporting NETL
Seth.King@NETL.DOE.GOV

@author: Nate Mitchell
LRST (Battelle|Leidos) supporting NETL
Nathaniel.Mitchell@NETL.DOE.GOV

@contributor: Veronika Vasylkivska
LRST (Battelle|Leidos) supporting NETL
Veronika.Vasylkivska@NETL.DOE.GOV
"""
import os
import sys
import logging
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

source_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(source_dir)

from openiam import IAM_DIR

RC_FONT = {'family': 'Arial', 'weight': 'normal', 'size': None}
BACKGROUND_COLOR = [0.67, 0.67, 0.67]

ATM_MAP_PLOT_RESERVOIR_COMPONENTS = ['LookupTableReservoir', 'SimpleReservoir',
                                     'AnalyticalReservoir']


class Plume():
    """
    Object to package ellipse data together
    Does not include rotation, must be oriented
    with x and y directions.
    """
    def __init__(self, x, y, dx, dy):
        self.x = x
        self.y = y
        self.dx = dx
        self.dy = dy

    def plot(self, ax, alpha=0.2, label_option=True, extent=None):
        """ Plot plume location."""
        if label_option:
            ax.plot(self.x, self.y, color='r', label='Source', marker='o',
                    markerfacecolor='none', markersize=10, markeredgewidth=2,
                    linestyle='none')
            e1 = mpatches.Ellipse(
                (self.x, self.y), (self.dx * 2), (self.dy * 2), alpha=alpha,
                label='Critical Areas')
        else:
            ax.plot(self.x, self.y, color='r', marker='o',
                    markerfacecolor='none', markersize=10, markeredgewidth=2,
                    linestyle='none')
            e1 = mpatches.Ellipse((self.x, self.y), (self.dx * 2),
                                  (self.dy * 2), alpha=alpha)

        ax.add_artist(e1)

        if not extent is None:
            x_adjust_factor = (extent[0][1] - extent[0][0]) / 60
            y_adjust_factor = (extent[1][1] - extent[1][0]) / 60
        else:
            x_adjust_factor = 1
            y_adjust_factor = 1

        ax.annotate("{0:.1f} m".format(self.dx), xy=(self.x, self.y),
                    xytext=(self.x + x_adjust_factor, self.y + y_adjust_factor),
                    size=18)


def make_plume_plots(plumes, time=None, alpha=0.2, sample=0, extent=None,
                     analysis='lhs', figsize=(15, 10), genfontsize=18, 
                     axislabelfontsize=24, titlefontsize=20, boldlabels=True, 
                     savefig=None, title=None, receptors=None, yaml_input_dict=None, 
                     res_comp_injX=None, res_comp_injY=None):
    """
    Takes in list of plume objects and plot them on a figure.

    :param plumes: list of plume objects to plot.
    :type plumes: lst

    :param time: Time step to display on figure title.
    :type time: float, int, str

    :param alpha: alpha (opacity) value for plumes.
    :type alpha: float in [0, 1]

    :param sample: index of realization
    :type sample: int

    :param extent: list of lists of (min, max) extent of the form
        [[xmin, xmax], [ymin, ymax]].  If not given will estimate from first sample.
    :type extent: list

    :param analysis: type of analysis. One of 'forward', 'lhs', 'parstudy'
    :type analysis: str

    :param genfontsize: font size for x- and y-ticks
    :type genfontsize: int

    :param axislabelfontsize: font size for x- and y-labels
    :type axislabelfontsize: int

    :param titlefontsize: font size for a figure title
    :type titlefontsize: int

    :param boldlabels: Flag indicating whether label font is bold
    :type boldlabels: boolean

    :param savefig: Directory and Filename to save the figure to.
    :type savefig: str

    :param title: Title of figure
    :type title: str

    :param receptors: location of receptors to plot
    :type receptors: list

    :param yaml_input_dict: part of yaml data dictionary containing setup for
        the given plot
    :type yaml_input_dict: dict

    :param res_comp_injX: list with x-coordinate(s) of injection well(s)
    :type res_comp_injX: list

    :param res_comp_injY: list with y-coordinate(s) of injection well(s)
    :type res_comp_injY: list

    :returns: None
    """
    if analysis == 'lhs':
        analysis_abbrev = 'LHS'
    elif analysis == 'parstudy':
        analysis_abbrev = 'parstudy'
    elif analysis == 'forward':
        analysis_abbrev = False

    # Number of columns in the legend
    ncol_number = 1

    # Figure formatting
    font = RC_FONT
    font['size'] = genfontsize
    plt.rc('font', **font)

    if boldlabels:
        selected_labelfontweight = 'bold'
    else:
        selected_labelfontweight = 'normal'

    if not yaml_input_dict is None:
        FigureDPI = yaml_input_dict['FigureDPI']

        plot_injection_sites = yaml_input_dict['plot_injection_sites']

        InjectionCoordx = yaml_input_dict['InjectionCoordx']
        InjectionCoordy = yaml_input_dict['InjectionCoordy']

        ncol_number += 1
    else:
        plot_injection_sites = False

        InjectionCoordx = None
        InjectionCoordy = None

    fig = plt.figure(figsize=figsize)
    plt.clf()

    ax = fig.add_subplot(1, 1, 1)

    plt.grid()

    max_x = -9.0e+99
    min_x = 9.0e+99
    max_y = -9.0e+99
    min_y = 9.0e+99
    max_dx = 0

    if receptors is None:
        receptors = []

    first_time_check = True
    for r in receptors:
        if first_time_check:
            plt.plot(r[0], r[1], 'k.', label='Receptor', markersize=10)
            first_time_check = False
            ncol_number += 1
        else:
            plt.plot(r[0], r[1], 'k.', markersize=10)

        if r[0] > max_x:
            max_x = r[0]
        if r[0] < min_x:
            min_x = r[0]
        if r[1] > max_y:
            max_y = r[1]
        if r[1] < min_y:
            min_y = r[1]

    first_time_check = True
    for p in plumes:
        if first_time_check:
            p.plot(ax, alpha=alpha, label_option=True, extent=extent)
            first_time_check = False
            ncol_number += 1
        else:
            p.plot(ax, alpha=alpha, label_option=False, extent=extent)

        if p.x > max_x:
            max_x = p.x
        if p.x < min_x:
            min_x = p.x
        if p.y > max_y:
            max_y = p.y
        if p.y < min_y:
            min_y = p.y
        if p.dx > max_dx:
            max_dx = p.dx

    if plot_injection_sites:
        min_x, max_x, min_y, max_y = check_injection_sites_for_min_max_xy(
            min_x, max_x, min_y, max_y, InjectionCoordx, InjectionCoordy,
            res_comp_injX, res_comp_injY)

        plot_injection_sites_func(
            InjectionCoordx=InjectionCoordx, InjectionCoordy=InjectionCoordy,
            res_comp_injX=res_comp_injX, res_comp_injY=res_comp_injY)

    if extent is None:
        extent = [
            ((min_x - (1.5 * max_dx)), (max_x + (1.5 * max_dx))),
            ((min_y - (1.5 * max_dx)), (max_y + (1.5 * max_dx)))]

    # Dummy plot elements used to enforce the axis limits. When using
    # plt.axis('equal'), it can make somewhat unexpected changes.
    plt.plot([extent[0][0], extent[0][1]], [extent[1][0], extent[1][0]],
              color='none', linewidth=1, linestyle='none', zorder=0)
    plt.plot([extent[0][1], extent[0][1]], [extent[1][0], extent[1][1]],
              color='none', linewidth=1, linestyle='none', zorder=0)
    plt.plot([extent[0][1], extent[0][0]], [extent[1][1], extent[1][1]],
              color='none', linewidth=1, linestyle='none', zorder=0)
    plt.plot([extent[0][0], extent[0][0]], [extent[1][1], extent[1][0]],
              color='none', linewidth=1, linestyle='none', zorder=0)

    ax.set(xlim=extent[0], ylim=extent[1])

    if not yaml_input_dict['EnforceXandYLims']:
        plt.axis('equal')

    if title is None:
        if time is not None:
            if analysis_abbrev:
                plt.title(
                    'Map View of CO$_2$ Release\nat Time t = {} years, {} Sample {}'.format(
                        time, analysis_abbrev, sample), fontsize=titlefontsize,
                    fontweight=selected_labelfontweight)
            else:
                plt.title(
                    'Map View of CO$_2$ Release\nat Time t = {} years'.format(
                        time), fontsize=titlefontsize,
                    fontweight=selected_labelfontweight)
    else:
        plt.title(title, fontsize=titlefontsize, fontweight=selected_labelfontweight)

    plt.xlabel('x (m)', fontsize=axislabelfontsize,
               fontweight=selected_labelfontweight)
    plt.ylabel('y (m)', fontsize=axislabelfontsize,
               fontweight=selected_labelfontweight)

    plt.xticks(fontsize=genfontsize)
    plt.yticks(fontsize=genfontsize)

    # Shrink current axis's height by 5% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.05,
                     box.width, box.height * 0.95])

    ax.legend(fancybox=False, fontsize=genfontsize - 2, ncol=ncol_number,
              edgecolor=[0, 0, 0], loc='best', bbox_to_anchor=(0.75, -0.1),
              framealpha=0.67)

    if savefig:
        plt.savefig(savefig, dpi=FigureDPI)
        plt.close()

    return extent


def make_grid(extent, npoints=(100, 100)):
    """
    Makes a grid over the given extent.

    :param extent: list of lists of (min, max) extent of the form
        [[xmin, xmax], [ymin, ymax]]
    :type extent: lst

    :param npoints: tuple of number of points in x and y direction (nx, ny).
    :type npoints: tuple

    :returns: ndarray(fl64) of size [nx*ny, 2] with (x, y) grid points.
    """
    if hasattr(npoints, '__iter__'):
        npoints_x = npoints[0]
        npoints_y = npoints[-1]
    else:
        npoints_x = npoints
        npoints_y = npoints

    x = np.linspace(extent[0][0], extent[0][1], num=npoints_x)
    y = np.linspace(extent[1][0], extent[1][1], num=npoints_y)

    xx, yy = np.meshgrid(x, y)
    xy = np.column_stack((xx.flatten(), yy.flatten()))

    return xy


def get_plumes(sm, s, atm_comp, time_index=0, sample=0, analysis='lhs',
               dr_factor=1.01):
    """
    Reads atmosphere plumes from sample data set and returns list of plume objects.

    :param sm: system model of interest
    :type sm: instance of SystemModel class

    :param s: Sample data set from completed analysis run.
    :type s: matk.sampleset

    :param atm_comp: Atmosphere rom component instance used in system model
    :type atm_comp: openiam.AtmosphericRom

    :param time_index: Index of time array to produce plumes for.
    :type time_index: int

    :param sample: Realization number to compute plumes from.
    :type sample: int

    :param analysis: type of analysis. One of 'forward', 'lhs', 'parstudy'
    :type analysis: str

    :param dr_factor: Padding factor to multiply radius by when calculating extent.
    :type dr_factor: float

    :returns: ([plumes], extent) List of plume objects and recommended extent to plot them on.
    """
    atm_name = atm_comp.name

    if analysis in ['lhs', 'parstudy']:
        num_source_real = s.recarray['{atm}.num_sources_{time}'.format(
            time=time_index, atm=atm_comp.name)]
    elif analysis == 'forward':
        num_source_real = sm.collect_observations_as_time_series(
            atm_comp, 'num_sources')
        num_source_real = num_source_real[time_index]

    plumes = []

    max_x = -9e99
    min_x = 9e99
    max_y = -9e99
    min_y = 9e99

    if analysis == 'forward':
        for source in range(int(num_source_real)):
            x = sm.collect_observations_as_time_series(
                atm_comp, 'x_new_s{source:03}'.format(source=source))
            x = x[time_index]

            y = sm.collect_observations_as_time_series(
                atm_comp, 'y_new_s{source:03}'.format(source=source))
            y = y[time_index]

            cd = sm.collect_observations_as_time_series(
                atm_comp, 'critical_distance_s{source:03}'.format(source=source))
            cd = cd[time_index]

            plumes.append(Plume(x, y, cd, cd))

            if (x + (dr_factor * cd)) > max_x:
                max_x = x + (dr_factor * cd)

            if (x - (dr_factor * cd)) < min_x:
                min_x = x - (dr_factor * cd)

            if (y + (dr_factor * cd)) > max_y:
                max_y = y + (dr_factor * cd)

            if (y - (dr_factor * cd)) < min_y:
                min_y = y - (dr_factor * cd)
    else:
        for source in range(int(num_source_real[sample])):
            x = s.recarray['{atm}.x_new_s{source:03}_{time}'.format(
                source=source, time=time_index, atm=atm_name)][sample]
            y = s.recarray['{atm}.y_new_s{source:03}_{time}'.format(
                source=source, time=time_index, atm=atm_name)][sample]
            cd = s.recarray['{atm}.critical_distance_s{source:03}_{time}'.format(
                source=source, time=time_index, atm=atm_name)][sample]

            plumes.append(Plume(x, y, cd, cd))

            if (x + (dr_factor * cd)) > max_x:
                max_x = x + (dr_factor * cd)

            if (x - (dr_factor * cd)) < min_x:
                min_x = x - (dr_factor * cd)

            if (y + (dr_factor * cd)) > max_y:
                max_y = y + (dr_factor * cd)

            if (y - (dr_factor * cd)) < min_y:
                min_y = y - (dr_factor * cd)

    extent = [(min_x, max_x), (min_y, max_y)]

    return (plumes, extent)


def map_plume_plot_single(plot_data, name, sm, s, satm, time_array, output_dir, 
                          analysis='lhs', figsize=(15, 10), sample=0,
                          genfontsize=18, axislabelfontsize=24, titlefontsize=20,
                          boldlabels=True, savefig=None, extent=None, title=None):
    """
    Maps atmosphere plumes from sample data set for single realization sample.

    :param plot_data: Dictionary of setup for a given plot
    :type plot_data: dict

    :param name: Name of plot defined by user
    :type name: str

    :param sm: system model of interest
    :type sm: instance of SystemModel class

    :param s: Sample data set from completed analysis run.
    :type s: matk.sampleset

    :param satm: Atmosphere rom component instance used in system model
    :type satm: openiam.AtmosphericRom

    :param time_array: Time array (in days) the system model used to generate
        the sample set
    :type time_array: list or np.array

    :param output_dir: Path to the folder where plots will be saved
    :type output_dir: str

    :param analysis: type of analysis. One of 'forward', 'lhs', 'parstudy'
    :type analysis: str

    :param sample: Realization number to compute plumes from.
    :type sample: int

    :param genfontsize: font size for x- and y-ticks
    :type genfontsize: int

    :param axislabelfontsize: font size for x- and y-labels
    :type axislabelfontsize: int

    :param titlefontsize: font size for a figure title
    :type titlefontsize: int

    :param boldlabels: Flag indicating whether label font is bold
    :type boldlabels: boolean

    :param savefig: Filename to save figure to. Include keys {time_index} or
        {time} for time series plots.
    :type savefig: str

    :param extent: list of lists of (min, max) extent of the form
        [[xmin, xmax], [ymin, ymax]].  If not given will estimate from first sample.
    :type extent: lst

    :param title: Title of figure
    :type title: str

    :returns: None
    """
    yaml_input_dict = get_yaml_input(plot_data, name, 'AtmPlumeSingle')

    sample = yaml_input_dict['Realization']

    if '.' in name:
        base_name= name[0:name.index('.')]
        extension = name[name.index('.'):None]
        name = name[0:name.index('.')] + '_{time_index}' + extension
    else:
        base_name = name
        name = name + '_{time_index}.png'

    if savefig is None:
        savefig = os.path.join(IAM_DIR, output_dir, name)
    else:
        if 'time' not in savefig and 'time_index' not in savefig:
            savefig = savefig + '_{time_index}.png'
        if '.png' not in savefig:
            savefig = savefig + '.png'
        savefig = os.path.join(IAM_DIR, output_dir, savefig)

    if not yaml_input_dict['PlotReceptors']:
        receptors = None

    if yaml_input_dict['plot_injection_sites']:
        res_comp_injX, res_comp_injY = get_res_comp_inj_xy(sm)
    else:
        res_comp_injX = [None]
        res_comp_injY = [None]

    if yaml_input_dict['EnforceXandYLims']:
        xLims = yaml_input_dict['xLims']
        yLims = yaml_input_dict['yLims']
        extent = [(xLims[0], xLims[1]), (yLims[0], yLims[1])]
    else:
        # Get the maximum extent over time so that all graphs can use it. The limits
        # can change between the plots, which makes comparing them difficult.
        extent = get_max_extent(
            sm, s, satm, time_array, analysis=analysis, plot_type='single',
            sample=sample, numSamples=1, yaml_input_dict=yaml_input_dict,
            res_comp_injX=res_comp_injX, res_comp_injY=res_comp_injY)

    # Plot results over time
    for time_index, time in enumerate(time_array/365.25):
        plumes, _ = get_plumes(
            sm, s, satm, time_index=time_index, sample=sample, analysis=analysis)

        if yaml_input_dict['PlotReceptors']:
            receptors = list(zip(satm.model_kwargs['x_receptor'],
                                 satm.model_kwargs['y_receptor']))

        savefig_ti = None
        if savefig:
            savefig_ti = savefig.format(time_index=time_index,
                                        time='{:.3f}'.format(time))

        title_ti = 'Map View of CO$_2$ Release\nat Time t = {} years'.format(time)

        if analysis in ['lhs', 'parstudy']:
            title_ti += ', Realization {}'.format(sample)

        if title:
            title_ti = title.format(time_index=time_index,
                                    time='{:.3f}'.format(time))
        make_plume_plots(
            plumes, time=time, alpha=0.2, sample=sample, extent=extent,
            analysis=analysis, figsize=figsize, genfontsize=genfontsize,
            axislabelfontsize=axislabelfontsize,
            titlefontsize=titlefontsize, boldlabels=boldlabels,
            savefig=savefig_ti, receptors=receptors, title=title_ti,
            yaml_input_dict=yaml_input_dict, res_comp_injX=res_comp_injX,
            res_comp_injY=res_comp_injY)

        if yaml_input_dict['SaveCSVFiles'] and time > 0:
            output_dir = os.path.join(IAM_DIR, output_dir)
            save_single_results_to_csv(
                plumes, time, base_name, output_dir, analysis=analysis,
                realization=yaml_input_dict['Realization'])


def ellipse_to_grid(grid_xy, center_xy, dxdy):
    """
    Returns a boolean grid of the same size as grid_xy with 1's
    where the ellipse is and zeros everywhere else.

    :param grid_xy: Array of size [nx*ny, 2] with (x,y) pairs of coordinates.
    :type grid_xy: ndarray(fl64)

    :param center_xy: tuple of (x, y) coordinates for center of ellipse.
    :type center_xy: tuple

    :param dxdy: tuple of major and minor axis of ellipse (dx, dy).
    :typle dxdy: tuple

    :returns: ndarray(fl64) Boolean grid values of shape [nx*ny]
    """
    ellipse = np.zeros_like(grid_xy[:, 0])
    if dxdy[0] != 0.0 and dxdy[1] != 0.0:
        etf = ((grid_xy[:, 0]-center_xy[0])/dxdy[0])**2 + (
            (grid_xy[:, 1]-center_xy[1])/dxdy[1])**2 <= 1
        ellipse[etf] = 1

    return ellipse


def prob_plume_grid(sm, s, satm, time_index, analysis='lhs', extent=None, xy=None,
                    yaml_input_dict=None):
    """
    Reads plumes for all realizations at a time step and return grid of existence
    probabilities as percentages.

    :param s: Sample data set from completed analysis run.
    :type s: matk.sampleset

    :param satm: Atmosphere rom component instance used in system model
    :type satm: openiam.AtmosphericRom

    :param time_index: Index of time array to produce plumes for.
    :type time_index: int

    :param extent: list of lists of (min, max) extent of the form
        [[xmin, xmax], [ymin, ymax]].  If not given will estimate from first sample.
    :type extent: lst

    :param xy: Array of size [nx*ny, 2] with (x,y) pairs of coordinates.
        If not given, will be created from extent.
    :type xy: ndarray(fl64)

    :returns: (plumes_grid, xy) Probabilities of plumes on a grid, and the xy grid coordinates.
    """
    if not extent:
        plumes, extent = get_plumes(
            sm, s, satm, time_index=time_index, sample=0, analysis=analysis)

    if yaml_input_dict['xGridSpacing']:
        x_range = extent[0][1] - extent[0][0]
        numXPoints = int(np.ceil(x_range / yaml_input_dict['xGridSpacing']) + 1)
    else:
        numXPoints = 101

    if yaml_input_dict['yGridSpacing']:
        y_range = extent[1][1] - extent[1][0]
        numYPoints = int(np.ceil(y_range / yaml_input_dict['yGridSpacing']) + 1)
    else:
        numYPoints = 101

    if xy is None:
        xy = make_grid(extent, npoints=(numXPoints, numYPoints))

    output = np.zeros_like(xy[:, 0])
    num_samples = s.recarray.shape[0]

    for sample in range(num_samples):
        el = np.zeros_like(xy[:, 0])

        plumes, _ = get_plumes(
            sm, s, satm, time_index=time_index, sample=sample, analysis=analysis)

        for p in plumes:
            el = np.logical_or(el, ellipse_to_grid(xy, (p.x, p.y), (p.dx, p.dy)))

        output += el

    output = (output / num_samples) * 100

    return (output, xy, extent)


def plot_grid_plumes(xy, plumes_grid, extent, savefig=None, time=None,
                     figsize=(15, 10), title=None, receptors=None,
                     genfontsize=18, axislabelfontsize=24, titlefontsize=20,
                     boldlabels=True, yaml_input_dict=None,
                     res_comp_injX=None, res_comp_injY=None):
    """
    Plots grid of probabilities of plumes.

    :param xy: Array of size [nx*ny, 2] with (x,y) pairs of coordinates.
    :type xy: ndarray(fl64)

    :param plumes_grid: Array of size [nx*ny] with probabilities
        of plume existing at (x, y) coordinate.
    :type plumes_grid: ndarray(fl64)

    :param extent: list of lists of (min, max) extent of the form
        [[xmin, xmax], [ymin, ymax]].  If not given will estimate from first sample.
    :type extent: list

    :param savefig: Filename to save figure to.  If not given, no file will be saved.
    :type savefig: str

    :param time: Time step to print to figure title.
    :type time: float, int, str

    :param figsize: tuple of figure size.
    :type figsize: tuple

    :param title: Title of figure
    :type title: str

    :param receptors: location of receptors to plot
    :type receptors: list

    :param genfontsize: font size for x- and y-ticks
    :type genfontsize: int

    :param axislabelfontsize: font size for x- and y-labels
    :type axislabelfontsize: int

    :param titlefontsize: font size for a figure title
    :type titlefontsize: int

    :param boldlabels: Flag indicating whether label font is bold
    :type boldlabels: boolean

    :param yaml_input_dict: part of yaml data dictionary containing setup for
        the given plot
    :type yaml_input_dict: dict

    :returns: None
    """
    # Number of columns used in legend
    ncol_number = 1

    if boldlabels:
        selected_labelfontweight = 'bold'
    else:
        selected_labelfontweight = 'normal'

    if yaml_input_dict is not None:
        FigureDPI = yaml_input_dict['FigureDPI']

        plot_injection_sites = yaml_input_dict['plot_injection_sites']

        InjectionCoordx = yaml_input_dict['InjectionCoordx']
        InjectionCoordy = yaml_input_dict['InjectionCoordy']

    else:
        plot_injection_sites = False

        InjectionCoordx = None
        InjectionCoordy = None

    plt.figure(figsize=figsize)
    plt.clf()

    ax = plt.gca()
    ax.set_facecolor(BACKGROUND_COLOR)

    plt.tricontourf(xy[:, 0], xy[:, 1], plumes_grid,
                    levels=np.arange(0.5, 101, 0.5))

    plt.plot([extent[0][0], extent[0][1]], [extent[1][0], extent[1][0]],
             color='k', linewidth=3, label='Grid Boundaries',
             linestyle='--', zorder=1.0e6)
    plt.plot([extent[0][1], extent[0][1]], [extent[1][0], extent[1][1]],
             color='k', linewidth=3, linestyle='--', zorder=1.0e6)
    plt.plot([extent[0][1], extent[0][0]], [extent[1][1], extent[1][1]],
             color='k', linewidth=3, linestyle='--', zorder=1.0e6)
    plt.plot([extent[0][0], extent[0][0]], [extent[1][1], extent[1][0]],
             color='k', linewidth=3, linestyle='--', zorder=1.0e6)

    if yaml_input_dict['PlotReceptors'] and not receptors is None:
        ncol_number += 1
        first_time_check = True
        for r in receptors:
            if first_time_check:
                plt.plot(r[0], r[1], 'k.', label='Receptor', markersize = 10)
                first_time_check = False
            else:
                plt.plot(r[0], r[1], 'k.', markersize = 10)

    if yaml_input_dict['EnforceXandYLims']:
        xLims = yaml_input_dict['xLims']
        yLims = yaml_input_dict['yLims']
        extent = [(xLims[0], xLims[1]), (yLims[0], yLims[1])]
        ax.set(xlim=extent[0], ylim=extent[1])
    else:
        plt.axis('equal')

    ticks = np.arange(10, 110, 10).tolist()
    cb = plt.colorbar(format='%.0f', ticks=ticks)
    cb.set_label(
        'Probability (%)', rotation=90, fontsize=axislabelfontsize,
        fontweight=selected_labelfontweight, labelpad=4)
    cb.ax.tick_params(labelsize=genfontsize)

    if plot_injection_sites:
        ncol_number += 1
        plot_injection_sites_func(
            InjectionCoordx=InjectionCoordx, InjectionCoordy=InjectionCoordy,
            res_comp_injX=res_comp_injX, res_comp_injY=res_comp_injY)

    # Shrink current axis's height by 5% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.05,
                     box.width, box.height * 0.95])

    ax.legend(fancybox=False, fontsize=genfontsize - 2, ncol=ncol_number,
              edgecolor=[0, 0, 0], loc='best', bbox_to_anchor=(0.75, -0.1),
              framealpha=0.67)

    plt.xlabel('x (m)', fontsize=axislabelfontsize,
               fontweight=selected_labelfontweight)
    plt.ylabel('y (m)', fontsize=axislabelfontsize,
               fontweight=selected_labelfontweight)

    plt.xticks(fontsize=genfontsize)
    plt.yticks(fontsize=genfontsize)

    if title is None:
        if time is not None:
            title = ''.join([
                'Probability of Atmospheric CO$_2$ Concentrations\n',
                'Exceeding Threshold for Time t = {time} years'.format(
                    time='{:.3f}'.format(time))])
        else:
            title = ''.join(['Probability of Atmospheric CO$_2$\n',
                             'Concentrations Exceeding Threshold'])

        title += ' (Gray within Boundaries: 0%)'

    plt.title(title, fontsize=titlefontsize, fontweight=selected_labelfontweight)

    if savefig:
        plt.savefig(savefig, dpi=FigureDPI)
        plt.close()


def map_plume_plot_ensemble(plot_data, name, sm, s, satm,
                            time_array, output_dir,
                            analysis='lhs', extent=None,
                            xy=None, axislabelfontsize=24, titlefontsize=20,
                            boldlabels=True, savefig=None,
                            figsize=(15, 10), title=None):
    """
    Reads plumes for all realizations at all time steps and plots the probabilities of existence.

    :param plot_data: Dictionary of setup for a given plot
    :type plot_data: dict

    :param name: Name of plot defined by user
    :type name: str

    :param sm: system model of interest
    :type sm: instance of SystemModel class

    :param s: Sample data set from completed analysis run.
    :type s: matk.sampleset

    :param satm: Atmosphere ROM component instance used in system model
    :type satm: openiam.AtmosphericRom

    :param time_array: Time array (in days) the system model used to generate
        the sample set
    :type time_array: list or np.array

    :param output_dir: Path to the folder where plots will be saved
    :type output_dir: str

    :param analysis: type of analysis. One of 'forward', 'lhs', 'parstudy'
    :type analysis: str

    :param extent: list of lists of (min, max) extent of the form
        [[xmin, xmax], [ymin, ymax]].  If not given will estimate from first sample.
    :type extent: lst

    :param xy: Array of size [nx*ny, 2] with (x, y) pairs of coordinates.
        If not given, will be created from extent.
    :type xy: ndarray(fl64)

    :param axislabelfontsize: font size for x- and y-labels
    :type axislabelfontsize: int

    :param titlefontsize: font size for a figure title
    :type titlefontsize: int

    :param boldlabels: Flag indicating whether label font is bold
    :type boldlabels: boolean

    :param savefig: Directory and filename to save figure to. Include keys
        {time_index} or {time} for time series plots.
    :type savefig: str

    :param figsize: tuple of figure size.
    :type figsize: tuple

    :param title: Title of figure
    :type title: str

    :returns: (plumes_grid, xy) Probabilities of plumes on a grid,
        and the xy grid coordinates.
    """
    continue_check = True

    if analysis == 'forward':
        debug_msg = ''.join(['The AtmPlumeEnsemble plot type was setup for the ',
                             'plot ', name, ', but the analysis type was ',
                             'forward. The AtmPlumeEnsemble plot type can ',
                             'only be used when the analysis type is Latin ',
                             'Hypercube Sampling (lhs) or parstudy. The plot ',
                             name, ' will not be made.'])
        logging.debug(debug_msg)
        continue_check = False

    if continue_check:
        numSamples = len(s.indices)

        if analysis == 'lhs':
            analysis_abbrev = 'LHS'
        elif analysis == 'parstudy':
            analysis_abbrev = 'parstudy'

        yaml_input_dict = get_yaml_input(plot_data, name, 'AtmPlumeEnsemble')

        if '.' in name:
            base_name = name[0:name.index('.')]
            extension = name[name.index('.'):None]
            name = name[0:name.index('.')] + '_{time_index}' + extension
        else:
            base_name = name
            name = name + '_{time_index}.png'

        if savefig is None:
            savefig = os.path.join(IAM_DIR, output_dir, name)
        else:
            if 'time' not in savefig and 'time_index' not in savefig:
                savefig = savefig + '_{time_index}.png'
            if '.png' not in savefig:
                savefig = savefig + '.png'
            savefig = os.path.join(IAM_DIR, output_dir, savefig)

        if yaml_input_dict['plot_injection_sites']:
            res_comp_injX, res_comp_injY = get_res_comp_inj_xy(sm)
        else:
            res_comp_injX = [None]
            res_comp_injY = [None]

        if yaml_input_dict['PlotReceptors']:
            receptors = list(zip(satm.model_kwargs['x_receptor'],
                                 satm.model_kwargs['y_receptor']))
        else:
            receptors = None

        savefig_ti = None

        if yaml_input_dict['EnforceXandYGridLims']:
            gridXLims = yaml_input_dict['gridXLims']
            gridYLims = yaml_input_dict['gridYLims']
            extent = [(gridXLims[0], gridXLims[1]), (gridYLims[0], gridYLims[1])]
        else:
            extent = get_max_extent(
                sm, s, satm, time_array, analysis=analysis, plot_type='ensemble',
                numSamples=numSamples, yaml_input_dict=yaml_input_dict,
                res_comp_injX=res_comp_injX, res_comp_injY=res_comp_injY)

        for time_index, time in enumerate(time_array[1:]/365.25):
            plume_grid, xy, _ = prob_plume_grid(
                sm, s, satm, time_index+1, analysis=analysis, extent=extent, xy=xy,
                yaml_input_dict=yaml_input_dict)

            if savefig:
                savefig_ti = savefig.format(time_index=time_index,
                                            time='{:.3f}'.format(time))

            title_ti = ''.join([
                'Probability of Atmospheric CO$_2$ Concentrations ',
                'Exceeding Threshold\nfor Time t = {} years, {} {} Realizations ',
                '(Gray within Boundaries: 0%)']).format(time, numSamples,
                                                        analysis_abbrev)

            if title:
                title_ti = title.format(time_index=time_index,
                                        time='{:.3f}'.format(time))

            plot_grid_plumes(
                xy, plume_grid, extent, savefig=savefig_ti, time=time,
                figsize=figsize, title=title_ti, receptors=receptors,
                axislabelfontsize=axislabelfontsize, titlefontsize=titlefontsize,
                boldlabels=boldlabels, yaml_input_dict=yaml_input_dict,
                res_comp_injX=res_comp_injX, res_comp_injY=res_comp_injY)

            if yaml_input_dict['SaveCSVFiles']:
                output_dir_upd = os.path.join(IAM_DIR, output_dir)
                save_ensemble_results_to_csv(
                    xy, plume_grid, time, base_name, output_dir_upd)


def get_yaml_input(plot_data, name, plot_type):
    """
    Function that reads the plot's input in the .yaml file and returns a
    dictionary containing the input.
    """
    defaultRealization = 0
    Realization = defaultRealization

    defaultFigureDPI = 100
    FigureDPI = defaultFigureDPI

    defaultPlotReceptors = False
    PlotReceptors = defaultPlotReceptors

    defaultPlotInjectionSites = False
    plot_injection_sites = defaultPlotInjectionSites
    InjectionCoordx = None
    InjectionCoordy = None

    EnforceXandYLims = False
    xLims = None
    yLims = None

    EnforceXandYGridLims = False
    gridXLims = None
    gridYLims = None

    defaultXGridSpacing = None
    xGridSpacing = defaultXGridSpacing
    defaultYGridSpacing = None
    yGridSpacing = defaultYGridSpacing

    defaultSaveCSVFiles = True
    SaveCSVFiles = defaultSaveCSVFiles

    grid_limits_debug_msg = ''.join([
        'The limits provided for the {}_grid used for the ', plot_type,' plot ',
        name, ' is not a list of length 2 (SpecifyXandYGridLims: grid{}Lims). ',
        'The {}_grid will be created in the default manner. Check your inputs ',
        'in the .yaml file.'])

    InjectionCoord_debug_msg = ''.join([
        'InjectionCoord{} was provided for the ', plot_type, ' plot ', name,
        ', but not InjectionCoord{}. Check your input. Injection sites will ',
        'not be displayed.'])

    InjectionCoord_debug_msg = ''.join([
        'InjectionCoord{} was provided for the ', plot_type, ' plot ', name,
        ', but not InjectionCoord{}. Check your input. Injection sites will ',
        'not be displayed.'])

    axis_lims_debug_msg = ''.join([
        'The {}-axis limits provided for the ', plot_type, ' plot ', name,
        ' (SpecifyXandYLims: {}Lims) are not of length 2 and will not be used. ',
        'Check your inputs in the .yaml file.'])

    input_type_debug_msg = ''.join([
        'The {} provided for the ', plot_type,' plot ', name,
        ' was not of type {}. The {} will ',
        'be set to the default value of {}.'])

    if plot_data[plot_type] is not None:
        if 'FigureDPI' in plot_data[plot_type]:
            FigureDPI = plot_data[plot_type]['FigureDPI']

            if not isinstance(FigureDPI, (int, float)):
                debug_msg = input_type_debug_msg.format(
                    'FigureDPI', 'int or float', 'FigureDPI',
                    str(defaultFigureDPI))
                logging.debug(debug_msg)
                FigureDPI = defaultFigureDPI

        if 'Realization' in plot_data[plot_type]:
            Realization = plot_data[plot_type]['Realization']

            if not isinstance(Realization, (int, float)):
                debug_msg = input_type_debug_msg.format(
                    'Realization', 'int or float', 'Realization',
                    str(defaultRealization))
                logging.debug(debug_msg)
                Realization = defaultRealization

        if 'PlotReceptors' in plot_data[plot_type]:
            PlotReceptors = plot_data[plot_type]['PlotReceptors']

            if not isinstance(PlotReceptors, bool):
                debug_msg = input_type_debug_msg.format(
                    'PlotReceptors', 'boolean', 'PlotReceptors',
                    str(defaultPlotReceptors))
                logging.debug(debug_msg)
                PlotReceptors = defaultPlotReceptors

        if 'SpecifyXandYGridLims' in plot_data[plot_type]:
            EnforceXandYGridLims = True
            gridXLims = plot_data[plot_type]['SpecifyXandYGridLims'][
                'gridXLims']
            gridYLims = plot_data[plot_type]['SpecifyXandYGridLims'][
                'gridYLims']

            gridXLimsWarning = False
            if isinstance(gridXLims, list):
                if len(gridXLims) != 2:
                    gridXLimsWarning = True

            if gridXLimsWarning:
                debug_msg = grid_limits_debug_msg.format('x', 'X', 'x')
                logging.debug(debug_msg)
                EnforceXandYGridLims = False
                gridXLims = None

            gridYLimsWarning = False
            if isinstance(gridYLims, list):
                if len(gridYLims) != 2:
                    gridYLimsWarning = True

            if gridYLimsWarning:
                debug_msg = grid_limits_debug_msg.format('y', 'Y', 'y')
                logging.debug(debug_msg)
                EnforceXandYGridLims = False
                gridYLims = None

        if 'PlotInjectionSites' in plot_data[plot_type]:
            plot_injection_sites = plot_data[plot_type]['PlotInjectionSites']
            if not isinstance(plot_injection_sites, bool):
                debug_msg = input_type_debug_msg.format(
                    'PlotInjectionSites', 'boolean', 'PlotInjectionSites',
                    str(defaultPlotInjectionSites))
                logging.debug(debug_msg)
                plot_injection_sites = defaultPlotInjectionSites

        InjectionCoordx = None
        InjectionCoordy = None
        if 'InjectionCoordx' in plot_data[plot_type]:
            InjectionCoordx = plot_data[plot_type]['InjectionCoordx']

        if 'InjectionCoordy' in plot_data[plot_type]:
            InjectionCoordy = plot_data[plot_type]['InjectionCoordy']
            if InjectionCoordx is None:
                debug_msg = InjectionCoord_debug_msg.format('y', 'x')
                logging.debug(debug_msg)
                plot_injection_sites = defaultPlotInjectionSites

            elif isinstance(InjectionCoordx, list) and isinstance(InjectionCoordy, list):
                if len(InjectionCoordx) != len(InjectionCoordy):
                    debug_msg = ''.join([
                        'The InjectionCoordy provided for the ', plot_type, ' ',
                        name, 'was not of the same length as the InjectionCoordx ',
                        'provided. Check your input. Injection sites will not ',
                        'be displayed.'])
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

        if 'SpecifyXandYLims' in plot_data[plot_type]:
            EnforceXandYLims = True
            xLims = plot_data[plot_type]['SpecifyXandYLims']['xLims']
            yLims = plot_data[plot_type]['SpecifyXandYLims']['yLims']

            if len(xLims) != 2:
                debug_msg = axis_lims_debug_msg.format('x', 'x')
                logging.debug(debug_msg)
                EnforceXandYLims = False

            if len(yLims) != 2:
                debug_msg = axis_lims_debug_msg.format('y', 'Y')
                logging.debug(debug_msg)
                EnforceXandYLims = False

        if 'xGridSpacing' in plot_data[plot_type]:
            xGridSpacing = plot_data[plot_type]['xGridSpacing']

            if not isinstance(xGridSpacing, (int, float)):
                debug_msg = input_type_debug_msg.format(
                    'xGridSpacing', 'int or float', 'xGridSpacing',
                    str(defaultXGridSpacing))
                logging.debug(debug_msg)
                xGridSpacing = defaultXGridSpacing

        if 'yGridSpacing' in plot_data[plot_type]:
            yGridSpacing = plot_data[plot_type]['yGridSpacing']

            if not isinstance(FigureDPI, (int, float)):
                debug_msg = input_type_debug_msg.format(
                    'yGridSpacing', 'int or float', 'yGridSpacing',
                    str(defaultYGridSpacing))
                logging.debug(debug_msg)
                yGridSpacing = defaultYGridSpacing

        if 'SaveCSVFiles' in plot_data[plot_type]:
            SaveCSVFiles = plot_data[plot_type]['SaveCSVFiles']

            if not isinstance(SaveCSVFiles, bool):
                debug_msg = input_type_debug_msg.format(
                    'SaveCSVFiles', 'boolean', 'SaveCSVFiles',
                    str(defaultSaveCSVFiles))
                logging.debug(debug_msg)
                SaveCSVFiles = defaultSaveCSVFiles

    yaml_input_dict = dict()

    yaml_input_dict['FigureDPI'] = FigureDPI

    yaml_input_dict['Realization'] = Realization

    yaml_input_dict['PlotReceptors'] = PlotReceptors

    yaml_input_dict['plot_injection_sites'] = plot_injection_sites

    yaml_input_dict['InjectionCoordx'] = InjectionCoordx
    yaml_input_dict['InjectionCoordy'] = InjectionCoordy

    yaml_input_dict['EnforceXandYLims'] = EnforceXandYLims
    yaml_input_dict['xLims'] = xLims
    yaml_input_dict['yLims'] = yLims

    yaml_input_dict['EnforceXandYGridLims'] = EnforceXandYGridLims
    yaml_input_dict['gridXLims'] = gridXLims
    yaml_input_dict['gridYLims'] = gridYLims

    yaml_input_dict['xGridSpacing'] = xGridSpacing
    yaml_input_dict['yGridSpacing'] = yGridSpacing

    yaml_input_dict['SaveCSVFiles'] = SaveCSVFiles

    return yaml_input_dict


def get_res_comp_inj_xy(sm):
    """
    Function that returns the x and y location(s) for injection site(s). This
    approach does not work for LookupTableReservoirs. InjectionCoordx and
    InjectionCoordy have to be provided in the .yaml file for simulations using
    LookupTableReservoir components.

    :param sm: system model of interest
    :type sm: instance of SystemModel class
    """
    res_comp_injX = []
    res_comp_injY = []

    components = list(sm.component_models.values())

    for comp in components:
        if comp.class_type in ATM_MAP_PLOT_RESERVOIR_COMPONENTS:
            # Get the injection sites
            if comp.class_type != 'LookupTableReservoir':
                res_comp_injX.append(comp.injX)
                res_comp_injY.append(comp.injY)
            else:
                res_comp_injX.append(None)
                res_comp_injY.append(None)

    return res_comp_injX, res_comp_injY


def get_max_extent(sm, s, satm, time_array, analysis='lhs', plot_type='single',
                   sample=0, numSamples=1, yaml_input_dict=None,
                   res_comp_injX=None, res_comp_injY=None):
    """
    Function that evaluates all realizations over time and returns the maximum
    extent from all realizations. This maximum extent is used to set the x and
    y limits.
    """
    min_x_over_t = 9.0e99
    max_x_over_t = -9.0e99
    min_y_over_t = 9.0e99
    max_y_over_t = -9.0e99

    if analysis == 'forward':
        for time_index in range(len(time_array)):
            _, ex = get_plumes(
                sm, s, satm, time_index=time_index, sample=sample, analysis=analysis)

            if ex[0][0] < min_x_over_t:
                min_x_over_t = ex[0][0]
            if ex[0][1] > max_x_over_t:
                max_x_over_t = ex[0][1]

            if ex[1][0] < min_y_over_t:
                min_y_over_t = ex[1][0]
            if ex[1][1] > max_y_over_t:
                max_y_over_t = ex[1][1]
    else:
        if plot_type == 'single':
            for time_index in range(len(time_array)):
                _, ex = get_plumes(
                    sm, s, satm, time_index=time_index, sample=sample,
                    analysis=analysis)

                if ex[0][0] < min_x_over_t:
                    min_x_over_t = ex[0][0]
                if ex[0][1] > max_x_over_t:
                    max_x_over_t = ex[0][1]

                if ex[1][0] < min_y_over_t:
                    min_y_over_t = ex[1][0]
                if ex[1][1] > max_y_over_t:
                    max_y_over_t = ex[1][1]

        elif plot_type == 'ensemble':
            for sampleRef in range(0, numSamples):
                for time_index in range(len(time_array)):
                    _, ex = get_plumes(
                        sm, s, satm, time_index=time_index, sample=sampleRef,
                        analysis=analysis)

                    if ex[0][0] < min_x_over_t:
                        min_x_over_t = ex[0][0]
                    if ex[0][1] > max_x_over_t:
                        max_x_over_t = ex[0][1]

                    if ex[1][0] < min_y_over_t:
                        min_y_over_t = ex[1][0]
                    if ex[1][1] > max_y_over_t:
                        max_y_over_t = ex[1][1]

    if yaml_input_dict['PlotReceptors'] and plot_type == 'single':
        receptors = list(zip(satm.model_kwargs['x_receptor'],
                             satm.model_kwargs['y_receptor']))
        if receptors is None:
            receptors = []

        min_x_over_t_i = min_x_over_t
        max_x_over_t_i = max_x_over_t
        min_y_over_t_i = min_y_over_t
        max_y_over_t_i = max_y_over_t

        for r in receptors:
            # If necessary, change the min and max x and y values (with a buffer)
            if r[0] > max_x_over_t:
                max_x_over_t = r[0] + ((max_x_over_t_i - min_x_over_t_i) / 20)
            if r[0] < min_x_over_t:
                min_x_over_t = r[0] - ((max_x_over_t_i - min_x_over_t_i) / 20)
            if r[1] > max_y_over_t:
                max_y_over_t = r[1] + ((max_y_over_t_i- min_y_over_t_i) / 20)
            if r[1] < min_y_over_t:
                min_y_over_t = r[1] - ((max_y_over_t_i - min_y_over_t_i) / 20)

    if yaml_input_dict['plot_injection_sites']:
        min_x, max_x, min_y, max_y = check_injection_sites_for_min_max_xy(
            min_x_over_t, max_x_over_t, min_y_over_t, max_y_over_t,
            yaml_input_dict['InjectionCoordx'], yaml_input_dict['InjectionCoordy'],
            res_comp_injX, res_comp_injY)

        # If the injection location is outside the min and max x and y, reset
        # the min and max x and y (with a buffer).
        if min_x < min_x_over_t:
            min_x_over_t = min_x - ((max_x_over_t - min_x_over_t) / 20)

        if max_x > max_x_over_t:
            max_x_over_t = max_x + ((max_x_over_t - min_x_over_t) / 20)

        if min_y < min_y_over_t:
            min_y_over_t = min_y - ((max_y_over_t - min_y_over_t) / 20)

        if max_y > max_y_over_t:
            max_y_over_t = max_y + ((max_y_over_t - min_y_over_t) / 20)

    extent = [(min_x_over_t, max_x_over_t), (min_y_over_t, max_y_over_t)]

    return extent


def check_injection_sites_for_min_max_xy(min_x, max_x, min_y, max_y,
                                         InjectionCoordx, InjectionCoordy,
                                         res_comp_injX, res_comp_injY):
    """
    Function that checks if the minimum and maximum x and y values need to be
    updated when the injection locations are included.
    """
    if isinstance(InjectionCoordx, list) and isinstance(InjectionCoordy, list):
        if None not in InjectionCoordx and InjectionCoordy is not None:
            if np.min(InjectionCoordx) < min_x:
                min_x = np.min(InjectionCoordx)
            if np.max(InjectionCoordx) > max_x:
                max_x = np.max(InjectionCoordx)

            if np.min(InjectionCoordy) < min_y:
                min_y = np.min(InjectionCoordy)
            if np.max(InjectionCoordy) > max_y:
                max_y = np.max(InjectionCoordy)
    else:
        if InjectionCoordx is not None and InjectionCoordy is not None:
            if InjectionCoordx < min_x:
                min_x = InjectionCoordx
            if InjectionCoordx > max_x:
                max_x = InjectionCoordx

            if InjectionCoordy < min_y:
                min_y = InjectionCoordy
            if InjectionCoordy > max_y:
                max_y = InjectionCoordy

    if not None in res_comp_injX and not None in res_comp_injY:
        if np.min(res_comp_injX) < min_x:
            min_x = np.min(res_comp_injX)
        if np.max(res_comp_injX) > max_x:
            max_x = np.max(res_comp_injX)

        if np.min(res_comp_injY) < min_y:
            min_y = np.min(res_comp_injY)
        if np.max(res_comp_injY) > max_y:
            max_y = np.max(res_comp_injY)

    return min_x, max_x, min_y, max_y


def plot_injection_sites_func(InjectionCoordx=None, InjectionCoordy=None,
                              res_comp_injX=None, res_comp_injY=None):
    """
    Function that plots injection sites.
    """
    if InjectionCoordx is None:
        InjectionCoordx = res_comp_injX
        InjectionCoordy = res_comp_injY

    if isinstance(InjectionCoordx, list):
        for injRef, (icoordX, icoordY) in enumerate(
                zip(InjectionCoordx, InjectionCoordy)):
            if injRef == 0:
                plt.plot(icoordX, icoordY, marker='s',
                         color='k', linestyle='none', markeredgewidth=2.5,
                         markersize=10, markerfacecolor='none', label='Injection Sites')
            else:
                plt.plot(icoordX, icoordY, marker='s',
                         color='k', linestyle='none', markeredgewidth=2.5,
                         markersize=10, markerfacecolor='none')
    else:
        plt.plot(InjectionCoordx, InjectionCoordy, marker='s',
                 color='k', linestyle='none', markeredgewidth=2.5, markersize=10,
                 markerfacecolor='none', label='Injection Site')


def save_single_results_to_csv(plumes, time, name, output_dir, analysis='lhs',
                               realization=1, min_crit_distance=1e-6):
    """
    Function that saves the results for AtmPlumeSingle plots
    (map_plume_plot_single) to a .csv file.
    """

    if not os.path.exists(os.path.join(output_dir, 'csv_files')):
        os.mkdir(os.path.join(output_dir, 'csv_files'))

    if analysis in ['lhs', 'parstudy']:
        name_addition = '_Realization_' + str(realization)
    else:
        name_addition = ''

    filename = name + name_addition + '_t_{:.0f}'.format(time) + '_years.csv'
    filename = os.path.join(output_dir, 'csv_files', filename)

    first_row = ['Source_Easting_x_m', 'Source_Northing_y_m',
                 'Critical_Distance_m',
                 'Note that only critical_distance results above {:.2e}'.format(
                     min_crit_distance) + ' m are shown']

    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(first_row)

        for p in plumes:
            if p.dx > min_crit_distance:
                row_temp = [p.x, p.y, p.dx]
                writer.writerow(row_temp)


def save_ensemble_results_to_csv(xy, plume_grid, time, name, output_dir,
                                 min_probability=1.0e-6):
    """
    Function that saves the results for AtmPlumeEnsemble plots
    (map_plume_plot_ensemble) to a .csv file.
    """

    if not os.path.exists(os.path.join(output_dir, 'csv_files')):
        os.mkdir(os.path.join(output_dir, 'csv_files'))

    try:
        filename = name + '_x_y_grid.csv'
        filename_grid = os.path.join(output_dir, 'csv_files', filename)

        if not os.path.isfile(filename_grid):
            first_row = ['Grid_x_value_m', 'Grid_y_value_m']

            with open(filename_grid, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(first_row)

                xy_shape = xy.shape
                for rowRef in range(0, xy_shape[0]):
                    row_temp = [xy[rowRef, 0], xy[rowRef, 1]]
                    writer.writerow(row_temp)
    except:
        pass

    filename = name + '_t_{:.0f}'.format(time) + '_years.csv'
    filename = os.path.join(output_dir, 'csv_files', filename)

    first_row = ['Easting_x_m', 'Northing_y_m',
                 'Probability_of_Atmospheric_CO2_Plume_percent',
                 'Note that only grid locations wit probabilities above '
                 + '{:.2e}'.format(min_probability) + '% are shown. All x and '
                 + 'y grid values are saved in a separate .csv file.']

    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(first_row)

        xy_shape = xy.shape
        for rowRef in range(0, xy_shape[0]):
            if plume_grid[rowRef] > min_probability:
                row_temp = [xy[rowRef, 0], xy[rowRef, 1], plume_grid[rowRef]]
                writer.writerow(row_temp)


if __name__ == '__main__':
    from openiam import SystemModel, SimpleReservoir, OpenWellbore, AtmosphericROM
    from matk import pyDOE

    # Define keyword arguments of the system model
    num_years = 5
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}  # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Create randomly located Leaky well locations
    # within box defined by xmin,xmax,ymin,ymax
    xymins = np.array([200., 550.])
    xymaxs = np.array([600., 700.])
    well_xys = xymins + pyDOE.lhs(2, samples=5)*(xymaxs-xymins)

    sress = []
    mss = []
    ow = []
    co2_leakrates_collector = []
    for i, crds in enumerate(well_xys):
        # Add reservoir components
        sress.append(sm.add_component_model_object(SimpleReservoir(
            name='sres'+str(i), parent=sm, injX=0., injY=0.,
            locX=crds[0], locY=crds[1])))

        # Add parameters of reservoir component model
        sress[-1].add_par('numberOfShaleLayers', value=3, vary=False)
        sress[-1].add_par('injRate', value=0.05, vary=False)
        sress[-1].add_par('shale1Thickness', min=300.0, max=500., value=400.0)
        sress[-1].add_par('aquifer1Thickness', min=20.0, max=60., value=50.0)
        sress[-1].add_par('shale2Thickness', min=700.0, max=900., value=800.0)
        sress[-1].add_par('aquifer2Thickness', min=60.0, max=80., value=75.0)
        sress[-1].add_par('shale3Thickness', min=150.0, max=250., value=200.0)
        sress[-1].add_par('logResPerm', min=-12.5, max=-11.5, value=-12.)

        # Add observations of reservoir component model to be used by the next component
        sress[-1].add_obs_to_be_linked('pressure')
        sress[-1].add_obs_to_be_linked('CO2saturation')

        # Add observations of reservoir component model
        sress[-1].add_obs('pressure')
        sress[-1].add_obs('CO2saturation')

        # Add open wellbore component
        ow.append(sm.add_component_model_object(OpenWellbore(name='ow'+str(i), parent=sm)))

        # Add parameters of open wellbore component
        ow[-1].add_par('wellRadius', min=0.026, max=0.03, value=0.028)
        ow[-1].add_par('logReservoirTransmissivity', min=-11.0, max=-9.0, value=-10.0)
        ow[-1].add_par('logAquiferTransmissivity', min=-11.0, max=-9.0, value=-10.0)
        ow[-1].add_par('brineSalinity', value=0.1, vary=False)
        ow[-1].add_par('wellTop', value=0.0, vary=False)

        # Add keyword arguments of the open wellbore component model
        ow[-1].add_kwarg_linked_to_obs('pressure', sress[-1].linkobs['pressure'])
        ow[-1].add_kwarg_linked_to_obs('CO2saturation', sress[-1].linkobs['CO2saturation'])

        # Add composite parameter of open wellbore component
        ow[-1].add_composite_par('reservoirDepth', expr='sres0.shale1Thickness' +
                                 '+sres0.shale2Thickness + sres0.shale3Thickness' +
                                 '+sres0.aquifer1Thickness+ sres0.aquifer2Thickness')

        # Add observations of multisegmented wellbore component model
        ow[-1].add_obs_to_be_linked('CO2_atm')
        ow[-1].add_obs('CO2_aquifer')
        ow[-1].add_obs('brine_aquifer')
        ow[-1].add_obs('CO2_atm')
        ow[-1].add_obs('brine_atm')

        # Create Collector for leakrates
        co2_leakrates_collector.append(ow[-1].linkobs['CO2_atm'])

    # Add Atm ROM
    satm = sm.add_component_model_object(AtmosphericROM(name='satm', parent=sm))
    satm.add_par('T_amb', min=13.0, max=15, value=18.0, vary=False)
    satm.add_par('P_amb', min=0.99, max=1.02, value=1.01E+00, vary=False)
    satm.add_par('V_wind', min=3.0, max=8.0, value=5.0, vary=False)
    satm.add_par('C0_critical', value=0.01, vary=False)

    satm.model_kwargs['x_receptor'] = [200, 250, 400, 350]
    satm.model_kwargs['y_receptor'] = [580, 660, 525, 600]

    satm.model_kwargs['x_coor'] = well_xys[:, 0]
    satm.model_kwargs['y_coor'] = well_xys[:, 1]

    satm.add_kwarg_linked_to_collection('co2_leakrate', co2_leakrates_collector)

    # Add Observations for receptors
    for i, r in enumerate(satm.model_kwargs['x_receptor']):
        satm.add_obs('outflag_r{0:03}'.format(i))
    n_sources = len(satm.model_kwargs['x_coor'])
    for i in range(n_sources):
        satm.add_obs('x_new_s{0:03}'.format(i))
        satm.add_obs('y_new_s{0:03}'.format(i))
        satm.add_obs('critical_distance_s{0:03}'.format(i))

    satm.add_obs('num_sources')

    import random
    num_samples = 50
    ncpus = 1
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=random.randint(500, 1100))

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)
    # Map single result
    output_dir = os.path.join(IAM_DIR, 'output', 'test_atm_plot')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    plot_data = {'AtmPlumeSingle': {
        'Realization': 5,
        'FigureDPI': 100,
        'SaveCSVFiles': False}}
    map_plume_plot_single(plot_data, 'test1', sm, s, satm, time_array, output_dir,
        savefig='atm_plume_single_t{time_index}', extent=[[100, 700], [450, 800]])

    # Map lhs results
    plot_data = {'AtmPlumeEnsemble': {
        'FigureDPI': 100,
        'SaveCSVFiles': False}}
    map_plume_plot_ensemble(plot_data, 'test2', sm, s, satm, time_array,
                            output_dir, savefig='atm_prob_plume_t{time_index}',
                            extent=[[100, 700], [450, 800]])
