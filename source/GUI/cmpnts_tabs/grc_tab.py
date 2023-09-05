"""
Module contains several methods needed for creating tab (page) in GUI
for GenericReservoir component. Methods read and write dictionaries
needed for control file interface yaml files.
"""

import tkinter as tk
from tkinter import ttk
from tkinter import StringVar
from tkinter import BooleanVar

from dictionarydata import componentVars, componentChoices
from dictionarydata import connectionsDictionary, componentTypeDictionary
from dictionarydata import DISTRIBUTION_OPTIONS

from dictionarydata import LABEL_FONT
from dictionarydata import (PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                            DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            OUTPUT_LABEL_WIDTH1, PARAMETER_FRAME_PADX, CB_PADX)

from cmpnts_tabs.locations import (add_inj_well_frame_widgets,
                                   add_obs_locs_frame_widgets,
                                   read_obs_locations_data)
from cmpnts_tabs.commons import commons_read_tab_vars


GRC_PARAMETERS = [
    'logResPerm', 'reservoirPorosity', 'resTempGradient',
    'initialSalinity', 'injRate']

GRC_PARAMETERS_SETUP = {
    'logResPerm': ["Reservoir permeability [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B2'),
                   'reservoir permeability'],
    'reservoirPorosity': ["Reservoir porosity [-]:",
                          'reservoir porosity'],
    'resTempGradient': ["Temperature gradient [{}C/km]:".format(u'\u00B0'),
                        'reservoir temperature gradient'],
    'initialSalinity': ["Initial salinity [-]:",
                        'reservoir initial salinity'],
    'injRate': ["CO{} injection rate [kg/s]:".format(u'\u2082'),
                "CO{} injection rate".format(u'\u2082')]}

# Set Simple Reservoir parameters names and value, min, max, second value, mean, std, bounds
GRC_PARAMETER_VALUES = {
    'logResPerm': [-14, -14.5, -13, -13, -13, 1, -15, -12],
    'reservoirPorosity': [0.15, 0.1, 0.3, 0.2, 0.15, 0.05, 0.08, 0.4],
    'resTempGradient': [22, 20, 24, 23, 25, 1, 18, 32],
    'initialSalinity': [0.003, 0.002, 0.005, 0.004, 0.01, 0.001, 0.001, 0.05],
    'injRate': [75, 40, 90, 50, 50, 5, 29, 179]}

GRC_OBSERVATIONS = ['pressure', 'CO2saturation']

GRC_OBSERVATIONS_SETUP = {
    'CO2saturation': ["CO{} saturation [-]".format(u'\u2082'),
                      'Enable CO{} saturation as output.'.format(u'\u2082')],
    'pressure': ["Pressure [Pa]", 'Enable pressure as output.']}


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""

    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'GenericReservoir',
        parameter_names=GRC_PARAMETERS, observation_names=GRC_OBSERVATIONS)

    cmpnt_data = read_obs_locations_data(cmpnt_data, cmpnt_nm)

    return cmpnt_data


def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip, *args):
    """ Add widgets to the component tab."""

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('GenericReservoir')
    connectionsDictionary.append('none')

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set('none')

    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    # Populate dictionary
    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(GRC_PARAMETER_VALUES)

    for output_key in GRC_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    comp_type_label = ttk.Label(
        tab, text="Generic Reservoir Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(GRC_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': GRC_PARAMETER_VALUES[par_name][6],
             'upper_bound': GRC_PARAMETER_VALUES[par_name][7]},
            GRC_PARAMETERS_SETUP[par_name][0],
            GRC_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, tool_tip)

    # Injection well location
    inj_well_frame = tk.Frame(tab)
    inj_well_frame.grid(row=6, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)
    add_inj_well_frame_widgets(controller, cmpnt_nm, 'GenericReservoir',
                               inj_well_frame, tool_tip)

    # Observation locations frame
    obs_locs_frame = tk.Frame(tab)
    obs_locs_frame.grid(row=7, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)
    add_obs_locs_frame_widgets(controller, cmpnt_nm, obs_locs_frame, tool_tip)

    outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    outputs_label.grid(row=8, column=0, sticky='w', pady=(5, 10))

    # Outputs
    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(row=9, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    # Create and place outputs widgets
    output_nms_labels = []
    output_nms_checkboxes = []
    for ind in range(2):
        obs_nm = GRC_OBSERVATIONS[ind]
        # Create output checkbox
        output_nms_checkboxes.append(tk.Checkbutton(
            outputs_frame, variable=componentVars[cmpnt_nm][obs_nm]))
        # Place checkbox
        output_nms_checkboxes[-1].grid(
            row=1, column=2*ind, pady=5, padx=CB_PADX, sticky='w')
        # Create output label
        output_nms_labels.append(ttk.Label(
            outputs_frame, text=GRC_OBSERVATIONS_SETUP[obs_nm][0],
            width=OUTPUT_LABEL_WIDTH1, anchor='w'))
        # Place label
        output_nms_labels[-1].grid(row=1, column=2*ind+1, pady=5, sticky='w')
        # Bind checkbox to the tip
        tool_tip.bind(output_nms_checkboxes[-1],
                      GRC_OBSERVATIONS_SETUP[obs_nm][1])
