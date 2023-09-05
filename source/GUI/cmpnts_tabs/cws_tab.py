"""
Module contains several methods needed for creating tab (page) in GUI
for ChemicalWellSealing component. Methods read and write dictionaries
needed for control file interface yaml files.
"""
import os
import sys

import tkinter as tk
from tkinter import ttk
from tkinter import StringVar, BooleanVar

# Save location of GUI folder
GUI_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(GUI_DIR)

from dictionarydata import componentVars, componentChoices
from dictionarydata import connectionsDictionary, componentTypeDictionary
from dictionarydata import DISTRIBUTION_OPTIONS

from dictionarydata import LABEL_FONT
from dictionarydata import (PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                            DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            OUTPUT_LABEL_WIDTH1, PARAMETER_FRAME_PADX, CB_PADX)
from cmpnts_tabs.commons import commons_read_tab_vars

CWS_PARAMETERS = ['fractureAperture', 'fractureLength', 'maxOverpressure']

CWS_PARAMETERS_SETUP = {
    'fractureAperture': ["Fracture aperture [m]:",
                         'aperture of the fractured leakage path'],
    'fractureLength': ["Fracture length [m]:",
                       'length of the fractured leakage path'],
    'maxOverpressure': ["Maximum overpressure [Pa]:",
                        'maximum overpressure']}

# Set ChemicalWellSealing parameter names and value, min, max, second value,
# mean, std, bounds
CWS_PARAMETER_VALUES = {
    'fractureAperture': [2.0e-5, 1.0e-5, 1.2e-4, 8.0e-5, 1.0e-4, 1.0e-5, 1.0e-5, 2.0e-3],
    'fractureLength': [20, 10, 100, 60, 100, 20, 10, 400],
    'maxOverpressure': [5.0e+6, 1.0e+6, 8.0e+6, 2.0e+6, 2.5e+6, 1.0e+5, 1.0e+6, 1.5e+7]}

CWS_OBSERVATIONS = ['seal_flag', 'seal_time']

CWS_OBSERVATIONS_SETUP = {
    'seal_flag': ["Seal flag", 'Enable flag for fracture sealing as output.'],
    'seal_time': ["Sealing time",
                  'Enable predicted time for fracture sealing as output.']}


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""

    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'ChemicalWellSealing', parameter_names=CWS_PARAMETERS,
        observation_names=CWS_OBSERVATIONS)

    return cmpnt_data

def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip, *args):
    """ Add widgets to the component tab."""

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('ChemicalWellSealing')
    connectionsDictionary.append('none')

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set('none')

    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(
        CWS_PARAMETER_VALUES)

    for obs_nm in CWS_OBSERVATIONS:
        componentVars[cmpnt_nm][obs_nm] = BooleanVar()
        componentVars[cmpnt_nm][obs_nm].set(0)

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Chemical Well Sealing Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(CWS_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': CWS_PARAMETER_VALUES[par_name][6],
             'upper_bound': CWS_PARAMETER_VALUES[par_name][7]},
            CWS_PARAMETERS_SETUP[par_name][0],
            CWS_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, tool_tip)

    # Outputs
    outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    outputs_label.grid(row=12, column=0, sticky='w', pady=(5, 10))

    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(row=13, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    output_nms_labels = []
    output_nms_checkboxes = []
    for ind in range(2):
        obs_nm = CWS_OBSERVATIONS[ind]
        # Create output checkbox
        output_nms_checkboxes.append(tk.Checkbutton(
            outputs_frame, variable=componentVars[cmpnt_nm][obs_nm]))
        # Place checkbox
        output_nms_checkboxes[-1].grid(
            row=1, column=2*ind, pady=5, padx=CB_PADX, sticky='w')
        # Create output label
        output_nms_labels.append(ttk.Label(
            outputs_frame, text=CWS_OBSERVATIONS_SETUP[obs_nm][0],
            width=OUTPUT_LABEL_WIDTH1, anchor='w'))
        # Place label
        output_nms_labels[-1].grid(row=1, column=2*ind+1, pady=5, sticky='w')
        # Bind checkbox to the tip
        tool_tip.bind(output_nms_checkboxes[-1],
                      CWS_OBSERVATIONS_SETUP[obs_nm][1])
