"""
Module contains several methods needed for creating tab (page) in GUI
for MultisegmentedWellbore component. Methods read and write dictionaries
needed for control file interface yaml files.
"""

import os
import sys
from re import split

import tkinter as tk
from tkinter import ttk
from tkinter import (StringVar, DoubleVar, IntVar, BooleanVar)

import Pmw

from dictionarydata import componentVars, componentChoices
from dictionarydata import connectionsDictionary, componentTypeDictionary
from dictionarydata import DISTRIBUTION_OPTIONS
from dictionarydata import connections

from dictionarydata import LABEL_FONT
from dictionarydata import (STRATA_PARAMETER_LABEL_WIDTH,
                            DISTRIBUTION_MENU_WIDTH,
                            DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            OUTPUT_LABEL_WIDTH1, PARAMETER_FRAME_PADX,
                            CB_PADX)

SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(SOURCE_DIR)


STRATA_PARAMETERS = ['numberOfShaleLayers', 'datumPressure', 'reservoirThickness']

# Set Multisegmented Wellbore parameter names and value, min, max, second value, mean, std
STRATA_PARAMETER_VALUES = {
    'shale': [100, 1, 100, 120, 200, 25, 1, 1600],
    'aquifer': [75, 1, 100, 100, 75, 10, 1, 1600],
    'reservoir': [50, 1, 100, 75, 50, 10, 1, 1600]}


def read_tab_vars(controller):
    """ Read values of tkinter variables associated with the component tab."""

    d = {}
    d['numberOfShaleLayers'] = {}
    d['reservoirThickness'] = {}
    num_shale_layers = componentVars['strata']['Params']['numberOfShaleLayers'].get()
    d['numberOfShaleLayers']['value'] = num_shale_layers
    d['numberOfShaleLayers']['vary'] = False
    d['datumPressure'] = componentVars['strata']['Params']['datumPressure'].get()

    par_names = ['reservoirThickness'] + [
        'shale{}Thickness'.format(ind) for ind in range(1, num_shale_layers+1)] + [
            'aquifer{}Thickness'.format(ind) for ind in range(1, num_shale_layers)]

    for par_nm in par_names:
        # Get the type of the distribution of the parameter
        distr_type = componentVars['strata']['Params'][par_nm]['distribution'].get()
        d[par_nm] = {}
        if distr_type == 'Fixed Value':
            d[par_nm]['value'] = componentVars['strata']['Params'][par_nm]['value'].get()
            d[par_nm]['vary'] = False

        if distr_type == 'Uniform':
            d[par_nm]['dist'] = 'uniform'
            d[par_nm]['min'] = componentVars['strata']['Params'][par_nm]['min'].get()
            d[par_nm]['max'] = componentVars['strata']['Params'][par_nm]['max'].get()

        if distr_type == 'Normal':
            d[par_nm]['dist'] = 'norm'
            d[par_nm]['mean'] = componentVars['strata']['Params'][par_nm]['mean'].get()
            d[par_nm]['std'] = componentVars['strata']['Params'][par_nm]['std'].get()

        if distr_type == 'Truncated':
            d[par_nm]['dist'] = 'truncnorm'
            kwargs = {
                'min': componentVars['strata']['Params'][par_nm]['min'].get(),
                'max': componentVars['strata']['Params'][par_nm]['max'].get(),
                'mean': componentVars['strata']['Params'][par_nm]['mean'].get(),
                'std': componentVars['strata']['Params'][par_nm]['std'].get()}
            dist_pars = controller.reparametrize_truncated_distribution(kwargs)
            d[par_nm]['dist_pars'] = [dist_pars['a'], dist_pars['b'],
                                      dist_pars['loc'], dist_pars['scale']]

        if distr_type == 'Lognormal':
            d[par_nm]['dist'] = 'lognorm'
            kwargs = {'mean': componentVars['strata']['Params'][par_nm]['mean'].get(),
                      'std': componentVars['strata']['Params'][par_nm]['std'].get()}
            dist_pars = controller.reparametrize_lognorm_distribution(kwargs)
            d[par_nm]['dist_pars'] = [
                dist_pars['s'], dist_pars['loc'], dist_pars['scale']]

        if distr_type == 'Triangular':
            d[par_nm]['dist'] = 'triang'
            kwargs = {'min': componentVars['strata']['Params'][par_nm]['min'].get(),
                      'max': componentVars['strata']['Params'][par_nm]['max'].get(),
                      'mode': componentVars['strata']['Params'][par_nm]['mode'].get()}
            dist_pars = controller.reparametrize_triang_distribution(kwargs)
            d[par_nm]['dist_pars'] = [
                dist_pars['c'], dist_pars['loc'], dist_pars['scale']]

        if distr_type == 'Discrete':
            d[par_nm]['discrete_vals'] = []
            discrete_values = []
            discrete_weights = []

            for value in componentVars['strata']['Params'][par_nm]['values'].get().split(','):
                discrete_values.append(float(value.strip()))

            for weight in componentVars['strata']['Params'][par_nm]['weights'].get().split(','):
                discrete_weights.append(float(weight.strip()))

            d[par_nm]['discrete_vals'].append(
                discrete_values)
            d[par_nm]['discrete_vals'].append(
                discrete_weights)

    return d

def add_widgets(controller, tab, toolTip):
    """ Add widgets to the Stratigraphy tab."""

    # Setup stratigraphy related global variables
    componentVars['strata']['Params'] = {}
    componentVars['strata']['Params']['numberOfShaleLayers'] = IntVar()
    counts = list(range(3, 31))
    componentVars['strata']['Params']['numberOfShaleLayers'].set(counts[0])

    componentVars['strata']['Params']['datumPressure'] = DoubleVar()
    componentVars['strata']['Params']['datumPressure'].set(101325)

    layers_par_names = [
        'shale3Thickness', 'aquifer2Thickness', 'shale2Thickness',
        'aquifer1Thickness', 'shale1Thickness', 'reservoirThickness']

    for par_name in layers_par_names:
        short_par_name = get_par_short_name(par_name)

        componentVars['strata']['Params'][par_name] = {}
        for key in ['distribution', 'values', 'weights']:
            componentVars['strata']['Params'][par_name][key] = StringVar()

        componentVars['strata']['Params'][par_name]['distribution'].set(
            DISTRIBUTION_OPTIONS[0])
        componentVars['strata']['Params'][par_name]['weights'].set('0.5, 0.5')
        componentVars['strata']['Params'][par_name]['values'].set(
            '{}, {}'.format(STRATA_PARAMETER_VALUES[short_par_name][0],
                            STRATA_PARAMETER_VALUES[short_par_name][3]))

        for key in ['value', 'min', 'max', 'mode', 'mean', 'std']:
            componentVars['strata']['Params'][par_name][key] = DoubleVar()

        componentVars['strata']['Params'][par_name]['value'].set(
            STRATA_PARAMETER_VALUES[short_par_name][0])
        componentVars['strata']['Params'][par_name]['min'].set(
            STRATA_PARAMETER_VALUES[short_par_name][1])
        componentVars['strata']['Params'][par_name]['max'].set(
            STRATA_PARAMETER_VALUES[short_par_name][2])
        componentVars['strata']['Params'][par_name]['mode'].set(
            STRATA_PARAMETER_VALUES[short_par_name][0])
        componentVars['strata']['Params'][par_name]['mean'].set(
            STRATA_PARAMETER_VALUES[short_par_name][4])
        componentVars['strata']['Params'][par_name]['std'].set(
            STRATA_PARAMETER_VALUES[short_par_name][5])

    # Stratigraphy tab label
    strata_label = ttk.Label(
        tab, text="Stratigraphy", font=LABEL_FONT)
    strata_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    par_frames['numberOfShaleLayers'] = tk.Frame(tab)
    par_frames['numberOfShaleLayers'].grid(row=1, column=0, sticky='w',
                                           padx=PARAMETER_FRAME_PADX)

    par_label = ttk.Label(par_frames['numberOfShaleLayers'],
                          width=STRATA_PARAMETER_LABEL_WIDTH,
                          text="Number of shale layers:")
    val_menu = tk.OptionMenu(
        par_frames['numberOfShaleLayers'],
        componentVars['strata']['Params']['numberOfShaleLayers'], *counts,
        command=lambda n: add_stratigraphy_layers(n, controller))
    val_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    par_label.grid(row=0, column=0, sticky='w', padx=5)
    val_menu.grid(row=0, column=1, padx=5)
    toolTip.bind(val_menu,
                 'Select the total number of shale layers for simulation.')

    hint_frame = tk.Frame(tab)
    hint_frame.grid(row=2, column=0, sticky='w', pady=10)
    controller.strata_layers_image = tk.PhotoImage(
        file=os.path.join(SOURCE_DIR, 'GUI', 'images', 'ShaleLayers.gif'))
    options = [controller.strata_layers_image]
    hint_menu_var = StringVar()
    hint_menu_var.set('Stratigraphy layers')
    hint_menu = tk.OptionMenu(hint_frame, hint_menu_var, *options)
    hint_menu.grid(row=0, column=0, sticky='w', padx=5)
    menu = hint_menu['menu']
    menu.delete("0", tk.END)
    menu.add_command(image=controller.strata_layers_image)
    toolTip.bind(hint_menu,
                 'Click to check the stratigraphy layers enumeration.')

    # Surface pressure frame
    par_frames['datumPressure'] = tk.Frame(tab)
    par_frames['datumPressure'].grid(row=3, column=0, sticky='w', pady=5,
                                     padx=PARAMETER_FRAME_PADX)

    par_label = ttk.Label(par_frames['datumPressure'],
                          width=STRATA_PARAMETER_LABEL_WIDTH,
                          text="Land surface pressure [Pa]:")
    val_field = tk.Entry(
        par_frames['datumPressure'],
        textvariable=componentVars['strata']['Params']['datumPressure'])
    val_field.config(width=DISTRIBUTION_ARG_TEXTFIELD_WIDTH)
    par_label.grid(row=0, column=0, sticky='w', padx=5)
    val_field.grid(row=0, column=1, padx=5)
    toolTip.bind(val_field,
                 'Enter pressure at the top of the uppermost shale layer.')

    # Create frames for other parameters
    for ind, par_name in enumerate(layers_par_names):
        short_par_name = get_par_short_name(par_name)

        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+4, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': STRATA_PARAMETER_VALUES[short_par_name][6],
             'upper_bound': STRATA_PARAMETER_VALUES[short_par_name][7]},
            get_par_label_text(par_name),
            get_tooltip_part(par_name),
            STRATA_PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars['strata']['Params'][par_name],
            'strata', toolTip)

    controller.strata_par_frames = par_frames
    controller.strata_tab = tab

def add_stratigraphy_layers(num_shale_layers, controller):
    """
    Add the number of shale and aquifer layers based on user selection.
    """
    tool_tip = Pmw.Balloon(controller)
    if num_shale_layers >= 7:
        controller.strata_scanv.config(
            scrollregion=(0, 0, 0, 550+(num_shale_layers-6)*63))
    else:
        controller.strata_scanv.config(
            scrollregion=(0, 0, 0, 0))

    par_list = list(controller.strata_par_frames.keys())
    for par_nm in par_list:
        # If parameter describes thickness of shale or aquifer layer
        if 'Thickness' in par_nm:
            for widget in controller.strata_par_frames[par_nm].winfo_children():
                # Delete all widgets from parameter (aquifer or shale thickness) frame
                tool_tip.unbind(widget)
                widget.destroy()
            frame = controller.strata_par_frames.pop(par_nm)
            frame.destroy()

    new_layers_par_names = [
        'shale{}Thickness'.format(ind) for ind in range(4, num_shale_layers+1)] + [
            'aquifer{}Thickness'.format(ind) for ind in range(3, num_shale_layers)]

    # Create variables for (possibly) new parameters
    for par_name in new_layers_par_names:
        short_par_name = get_par_short_name(par_name)

        componentVars['strata']['Params'][par_name] = {}
        for key in ['distribution', 'values', 'weights']:
            componentVars['strata']['Params'][par_name][key] = StringVar()

        componentVars['strata']['Params'][par_name]['distribution'].set(
            DISTRIBUTION_OPTIONS[0])
        componentVars['strata']['Params'][par_name]['weights'].set('0.5, 0.5')
        componentVars['strata']['Params'][par_name]['values'].set(
            '{}, {}'.format(STRATA_PARAMETER_VALUES[short_par_name][0],
                            STRATA_PARAMETER_VALUES[short_par_name][3]))

        for key in ['value', 'min', 'max', 'mode', 'mean', 'std']:
            componentVars['strata']['Params'][par_name][key] = DoubleVar()

        componentVars['strata']['Params'][par_name]['value'].set(
            STRATA_PARAMETER_VALUES[short_par_name][0])
        componentVars['strata']['Params'][par_name]['min'].set(
            STRATA_PARAMETER_VALUES[short_par_name][1])
        componentVars['strata']['Params'][par_name]['max'].set(
            STRATA_PARAMETER_VALUES[short_par_name][2])
        componentVars['strata']['Params'][par_name]['mode'].set(
            STRATA_PARAMETER_VALUES[short_par_name][0])
        componentVars['strata']['Params'][par_name]['mean'].set(
            STRATA_PARAMETER_VALUES[short_par_name][4])
        componentVars['strata']['Params'][par_name]['std'].set(
            STRATA_PARAMETER_VALUES[short_par_name][5])

    # Define names of layers parameters
    layers_par_names = ['shale{}Thickness'.format(num_shale_layers)]
    for ind in range(1, num_shale_layers):
        layers_par_names = layers_par_names + [
            'aquifer{}Thickness'.format(num_shale_layers-ind),
            'shale{}Thickness'.format(num_shale_layers-ind)]
    layers_par_names = layers_par_names + ['reservoirThickness']

    # Create frames for other parameters
    for ind, par_name in enumerate(layers_par_names):
        short_par_name = get_par_short_name(par_name)
        controller.strata_par_frames[par_name] = tk.Frame(controller.strata_tab)
        controller.strata_par_frames[par_name].grid(
            row=ind+4, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            controller.strata_par_frames[par_name], par_name,
            {'lower_bound': STRATA_PARAMETER_VALUES[short_par_name][6],
             'upper_bound': STRATA_PARAMETER_VALUES[short_par_name][7]},
            get_par_label_text(par_name),
            get_tooltip_part(par_name),
            STRATA_PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars['strata']['Params'][par_name],
            'strata', tool_tip)

def get_par_label_text(name):
    if 'aquifer' in name:
        ind = split('Thickness', name)[0][7:]
        return 'Aquifer {} thickness [m]:'.format(ind)
    elif 'shale' in name:
        ind = split('Thickness', name)[0][5:]
        return 'Shale {} thickness [m]:'.format(ind)
    else:
        return 'Reservoir thickness [m]:'

def get_tooltip_part(name):
    if 'aquifer' in name:
        return 'aquifer thickness'
    elif 'shale' in name:
        return 'shale thickness'
    else:
        return 'reservoir thickness'

def get_par_short_name(name):
    if 'aquifer' in name:
        return 'aquifer'
    elif 'shale' in name:
        return 'shale'
    else:
        return 'reservoir'