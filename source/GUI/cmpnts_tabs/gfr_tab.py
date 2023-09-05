"""
Module contains several methods needed for creating tab (page) in GUI
for GeneralizedFlowRate component. Methods read and write dictionaries
needed for control file interface yaml files.
"""
import os
import sys

import numpy as np

import tkinter as tk
from tkinter import ttk
from tkinter import StringVar, IntVar, BooleanVar

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
from dictionarydata import GFR_PARAMETER_LABEL_WIDTH

from cmpnts_tabs.locations import read_locations_data, add_wellbore_frame_widgets
from cmpnts_tabs.commons import commons_read_tab_vars


GFR_PARAMETERS = [
    'logPeakCO2Rate', 'timePeakCO2Rate', 'durationPeakCO2Rate',
    'durationPeakZeroCO2Rate', 'logInitBrineRate', 'logFinalBrineRate',
    'durationInitBrineRate', 'durationInitFinalBrineRate', 'mitigationTime']

GFR_PARAMETERS_SETUP = {
    'logPeakCO2Rate': [
        "Logarithm of the max CO{} rate [log{} kg/s]:".format(
            u'\u2082', u'\u2081'u'\u2080'),
        'logarithm of the maximum CO{} rate'.format(u'\u2082')],
    'timePeakCO2Rate': ["Time to the max CO{} rate [years]:".format(u'\u2082'),
                        'time to the maximum CO{} rate'.format(u'\u2082')],
    'durationPeakCO2Rate': ["Duration of the max CO{} rate [years]:".format(u'\u2082'),
                            'duration of the maximum CO{} rate'.format(u'\u2082')],
    'durationPeakZeroCO2Rate': [
        "Duration of CO{} rate decline [years]:".format(u'\u2082'),
        'duration of CO{} rate decline'.format(u'\u2082')],
    'logInitBrineRate': [
        "Logarithm of initial brine rate [log{} kg/s]:".format(u'\u2081'u'\u2080'),
        'logarithm of initial brine rate'],
    'logFinalBrineRate': [
        "Logarithm of final brine rate [log{} kg/s]:".format(u'\u2081'u'\u2080'),
        'logarithm of final brine rate'],
    'durationInitBrineRate': ["Duration of initial brine rate [years]:",
                              'duration of initial brine rate'],
    'durationInitFinalBrineRate': [
        "Duration of transient brine rate [years]:",
        'duration of transition from initial to final brine rate'],
    'mitigationTime': ["Mitigation time [years]:",
                       'mitigation time']}

# Set Generalized Flow rate parameter names and value, min, max, second value,
# mean, std, bounds
GFR_PARAMETER_VALUES = {
    'logPeakCO2Rate': [-5, -10, 5, -6, -5.5, -1, -np.inf, 5],
    'timePeakCO2Rate': [5, 1, 100, 8, 6, 3, 0, 1000],
    'durationPeakCO2Rate': [10, 1, 100, 15, 10, 1, 0, 1000],
    'durationPeakZeroCO2Rate': [10, 1, 100, 30, 15, 5, 0, 1000],
    'logInitBrineRate': [-10, -15, 5, -5, -5, 2, -np.inf, 5],
    'logFinalBrineRate': [-11, -15, 5, -5, -5, 5, -np.inf, 5],
    'durationInitBrineRate': [10, 1, 100, 5, 8, 5, 0, 1000],
    'durationInitFinalBrineRate': [10, 1, 100, 20, 10, 5, 0, 1000],
    'mitigationTime': [100, 5, 50, 200, 50, 10, 0, np.inf]}

GFR_OBSERVATIONS = ['CO2_aquifer', 'brine_aquifer']

GFR_OBSERVATIONS_SETUP = {
    'CO2_aquifer': [
        "CO{} aquifer [kg/s]".format(u'\u2082'),
        'Enable CO{} leakage rates to aquifer as output.'.format(u'\u2082')],
    'brine_aquifer': [
        "Brine aquifer [kg/s]",
        'Enable brine leakage rates to aquifer as output.']}

def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'GeneralizedFlowRate', parameter_names=GFR_PARAMETERS,
        add_number=True, add_leak_to=True)

    cmpnt_data['numberOfShaleLayers'] = componentVars['strata']['Params'][
        'numberOfShaleLayers'].get()

    # Read information about locations associated with component
    read_locations_data(cmpnt_data, cmpnt_nm)

    model_outputs = []
    for key in componentVars[cmpnt_nm]['outputs']:
        if componentVars[cmpnt_nm]['outputs'][key].get():
            model_outputs.append(key+cmpnt_data['LeakTo'][7:])

    cmpnt_data['Outputs'] = model_outputs
    return cmpnt_data

def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip,
                cnctn_nm, dyn_data, aquifer_nm, *args):
    """ Add widgets to the component tab."""

    aquifers = ['aquifer{}'.format(ind) for ind in range(
        1, componentVars['strata']['Params']['numberOfShaleLayers'].get())]

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('GeneralizedFlowRate')
    connectionsDictionary.append('none')

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set('none')

    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    componentVars[cmpnt_nm]['number'] = IntVar()
    componentVars[cmpnt_nm]['number'].set(1)
    componentVars[cmpnt_nm]['LeakTo'] = StringVar()
    componentVars[cmpnt_nm]['LeakTo'].set(aquifer_nm)

    componentVars[cmpnt_nm]['useRandomLocDomain'] = BooleanVar()
    componentVars[cmpnt_nm]['useRandomLocDomain'].set(0)

    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(
        GFR_PARAMETER_VALUES)

    componentVars[cmpnt_nm]['outputs'] = {}
    for obs_nm in GFR_OBSERVATIONS:
        componentVars[cmpnt_nm]['outputs'][obs_nm] = BooleanVar()
        componentVars[cmpnt_nm]['outputs'][obs_nm].set(0)

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Generalized Flow Rate Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(GFR_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': GFR_PARAMETER_VALUES[par_name][6],
             'upper_bound': GFR_PARAMETER_VALUES[par_name][7]},
            GFR_PARAMETERS_SETUP[par_name][0],
            GFR_PARAMETERS_SETUP[par_name][1],
            GFR_PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, tool_tip)

    # Leak to frame
    leak_to_frame = tk.Frame(tab)
    leak_to_frame.grid(row=10, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)
    # Leaked to aquifer name
    leak_to_label = ttk.Label(leak_to_frame, text='Leak to:',
                              width=GFR_PARAMETER_LABEL_WIDTH)
    leak_to_menu = tk.OptionMenu(
        leak_to_frame, componentVars[cmpnt_nm]['LeakTo'], *aquifers)
    leak_to_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    tool_tip.bind(
        leak_to_menu,
        'Select which of the aquifer layers this component model will represent.')
    leak_to_label.grid(row=0, column=0, padx=5, sticky='w')
    leak_to_menu.grid(row=0, column=1, padx=5)

    sources_number_label = ttk.Label(
        leak_to_frame, text='Number of sources:', width=PARAMETER_LABEL_WIDTH)
    number_spinbox = tk.Spinbox(leak_to_frame, from_=1, to=1000,
                                textvariable=componentVars[cmpnt_nm]['number'])
    sources_number_label.grid(row=1, column=0, sticky='w', padx=5)
    number_spinbox.grid(row=1, column=1, padx=5, pady=2, sticky='ew')
    tool_tip.bind(number_spinbox,
                  ''.join(['Set the total number of leakage sources for this ',
                           'component including random locations.']))

    # Wellbore locations frame
    well_locs_frame = tk.Frame(tab)
    well_locs_frame.grid(row=11, column=0, sticky='w',
                         padx=PARAMETER_FRAME_PADX, pady=(5, 10))
    add_wellbore_frame_widgets(
        controller, cmpnt_nm, well_locs_frame, tool_tip, gfr_cmpnt=True)

    # Outputs
    outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    outputs_label.grid(row=12, column=0, sticky='w', pady=(5, 10))

    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(row=13, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    output_nms_labels = []
    output_nms_checkboxes = []
    for ind in range(2):
        obs_nm = GFR_OBSERVATIONS[ind]
        # Create output checkbox
        output_nms_checkboxes.append(tk.Checkbutton(
            outputs_frame, variable=componentVars[cmpnt_nm]['outputs'][obs_nm]))
        # Place checkbox
        output_nms_checkboxes[-1].grid(
            row=1, column=2*ind, pady=5, padx=CB_PADX, sticky='w')
        # Create output label
        output_nms_labels.append(ttk.Label(
            outputs_frame, text=GFR_OBSERVATIONS_SETUP[obs_nm][0],
            width=OUTPUT_LABEL_WIDTH1, anchor='w'))
        # Place label
        output_nms_labels[-1].grid(row=1, column=2*ind+1, pady=5, sticky='w')
        # Bind checkbox to the tip
        tool_tip.bind(output_nms_checkboxes[-1],
                      GFR_OBSERVATIONS_SETUP[obs_nm][1])
