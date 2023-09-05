"""
Module contains several methods needed for creating tab (page) in GUI
for FaultLeakage component. Methods read and write dictionaries
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
from dictionarydata import (FL_PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                            DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            OUTPUT_LABEL_WIDTH1, PARAMETER_FRAME_PADX, CB_PADX)
from cmpnts_tabs.commons import commons_read_tab_vars

FL_PARAMETERS = ['damage_zone_perm', 'damage_zone_por', 'dip_angle',
                 'shallow_aquifer_perm', 'shallow_aquifer_por',
                 'deep_aquifer_perm', 'deep_aquifer_por',
                 'well_index', 'well_rate', 'injection_time',
                 'geothermal_gradient']

FL_PARAMETERS_SETUP = {
    'damage_zone_perm': ["Fault permeability [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B2'), 'permeability of the fault'],
    'damage_zone_por': ["Fault porosity [-]:", 'porosity of the fault'],
    'dip_angle': ["Dip angle [deg]:", 'dip angle of the fault'],
    'shallow_aquifer_perm': ["Shallow aquifer permeability [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B2'), 'permeability of the shallow aquifer'],
    'shallow_aquifer_por': ["Shallow aquifer porosity [-]:",
                            'porosity of the shallow aquifer'],
    'deep_aquifer_perm': ["Deep aquifer permeability [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B2'),
                          'permeability of the deep aquifer (storage/reservoir)'],
    'deep_aquifer_por': ["Deep aquifer porosity [-]:",
                         'porosity of the deep aquifer (storage/reservoir)'],
    'well_index': ['Well index [-]:',
                   'well index (proxy for the horizontal distance of the well from the fault)'],
    'well_rate': ['Injection rate [kg/s]:', 'injection rate of the well'],
    'injection_time': ['Injection time [years]:', 'duration of injection'],
    'geothermal_gradient': ['Geothermal gradient [{}C/km]:'.format(u'\u00B0'),
                            'geothermal gradient in the formation']}

# Set FaultLeakage parameter names and value, min, max, second value,
# mean, std, and lower and upper bounds
FL_PARAMETER_VALUES1 = {
    'damage_zone_perm': [-13.5, -14.0, -13.0, -13.1, -13.0, 1.0, -15.0, -12.0],
    'damage_zone_por': [0.01, 0.005, 0.05, 0.02, 0.05, 0.002, 0.001, 0.1],
    'shallow_aquifer_perm': [-13.0, -13.5, -12.5, -12.8, -13, 0.5, -14.0, -12.0],
    'shallow_aquifer_por': [0.25, 0.1, 0.5, 0.3, 0.25, 0.01, 0.05, 0.5],
    'deep_aquifer_perm': [-13.0, -13.5, -12.5, -12.8, -13, 0.5, -14.0, -12.0],
    'deep_aquifer_por': [0.2, 0.1, 0.3, 0.25, 0.2, 0.025, 0.05, 0.35],
    'well_rate': [15.0, 5.0, 20.0, 10.0, 15.0, 2.5, 0.5, 25.0],
    'injection_time': [30.0, 20.0, 40.0, 25.0, 20.0, 5.0, 10.0, 50.0],
    'geothermal_gradient': [30.0, 10.0, 40.0, 25.0, 28.0, 2.0, 8.0, 44.0]}

FL_PARAMETER_VALUES2 = {
    'dip_angle': [60, 40, 140, 40, 80, 10, 40, 140],
    'well_index': [0, 0, 2, 1, 1, 0, 0, 2]}

FL_DISCRETE_PARAMETER_BOUNDS = {
    'dip_angle': [40, 60, 80, 100, 120, 140],
    'well_index': [0, 1, 2]}

FL_OBSERVATIONS = ['CO2_aquifer', 'brine_aquifer',
                   'mass_CO2_aquifer','mass_brine_aquifer']

FL_OBSERVATIONS_SETUP = {
    'CO2_aquifer': [
        "CO{} aquifer [kg/s]".format(u'\u2082'),
        'Enable CO{} leakage rates to aquifer as output.'.format(u'\u2082')],
    'brine_aquifer': [
        "Brine aquifer [kg/s]",
        'Enable brine leakage rates to aquifer as output.'],
    'mass_CO2_aquifer': [
        'Mass CO{} aquifer [kg]'.format(u'\u2082'),
        'Enable CO{} mass leaked to aquifer as output.'.format(u'\u2082')],
    'mass_brine_aquifer': [
        'Mass brine aquifer [kg]',
        'Enable brine mass leaked to aquifer as output.']}


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""

    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'FaultLeakage', parameter_names=FL_PARAMETERS,
        observation_names=FL_OBSERVATIONS)

    return cmpnt_data

def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip, *args):
    """ Add widgets to the component tab."""

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('FaultLeakage')
    connectionsDictionary.append('none')

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set('none')

    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    # For parameters that can assume any of the distributions
    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(
        FL_PARAMETER_VALUES1)
    # For parameters that can assume only discrete distribution
    componentVars[cmpnt_nm]['Params'].update(controller.populate_params_dict(
        FL_PARAMETER_VALUES2, distr_options=['Fixed Value', 'Discrete'],
        options = {1: ['distribution', 'values', 'weights'],
                   2: ['value']}))

    for obs_nm in FL_OBSERVATIONS:
        componentVars[cmpnt_nm][obs_nm] = BooleanVar()
        componentVars[cmpnt_nm][obs_nm].set(0)

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Fault Leakage Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(FL_PARAMETERS):

        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)
        if ind not in [2, 7]:
            controller.setup_parameter_frame(
                par_frames[par_name], par_name,
                {'lower_bound': FL_PARAMETER_VALUES1[par_name][6],
                 'upper_bound': FL_PARAMETER_VALUES1[par_name][7]},
                FL_PARAMETERS_SETUP[par_name][0],
                FL_PARAMETERS_SETUP[par_name][1],
                FL_PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
                cmpnt_nm, tool_tip)
        else:  # for parameters dip_angle and well_index
            controller.setup_parameter_frame(
                par_frames[par_name], par_name,
                {'lower_bound': FL_PARAMETER_VALUES2[par_name][-2],
                 'upper_bound': FL_PARAMETER_VALUES2[par_name][-1],
                 'discrete_bounds': FL_DISCRETE_PARAMETER_BOUNDS[par_name]},
                FL_PARAMETERS_SETUP[par_name][0],
                FL_PARAMETERS_SETUP[par_name][1],
                FL_PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                ['Fixed Value', 'Discrete'], componentVars[cmpnt_nm]['Params'][par_name],
                cmpnt_nm, tool_tip)

    # Outputs
    outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    outputs_label.grid(row=12, column=0, sticky='w', pady=(5, 10))

    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(row=13, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    output_nms_labels = []
    output_nms_checkboxes = []
    for ind in range(4):
        obs_nm = FL_OBSERVATIONS[ind]
        # Create output checkbox
        output_nms_checkboxes.append(tk.Checkbutton(
            outputs_frame, variable=componentVars[cmpnt_nm][obs_nm]))
        # Place checkbox
        output_nms_checkboxes[-1].grid(
            row=ind//2+1, column=2*(ind%2), pady=5, padx=CB_PADX, sticky='w')
        # Create output label
        output_nms_labels.append(ttk.Label(
            outputs_frame, text=FL_OBSERVATIONS_SETUP[obs_nm][0],
            width=OUTPUT_LABEL_WIDTH1, anchor='w'))
        # Place label
        output_nms_labels[-1].grid(row=ind//2+1, column=2*(ind%2)+1, pady=5, sticky='w')
        # Bind checkbox to the tip
        tool_tip.bind(output_nms_checkboxes[-1],
                      FL_OBSERVATIONS_SETUP[obs_nm][1])
