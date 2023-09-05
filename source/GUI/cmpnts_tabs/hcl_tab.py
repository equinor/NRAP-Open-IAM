"""
Module contains several methods needed for creating tab (page) in GUI
for HydrocarbonLeakage component. Methods read and write dictionaries
needed for control file interface yaml files.
"""
import os
import sys

import tkinter as tk
from tkinter import ttk
from tkinter import StringVar, BooleanVar
import numpy as np

# Save location of GUI folder
GUI_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(GUI_DIR)

from dictionarydata import componentVars, componentChoices
from dictionarydata import connectionsDictionary, componentTypeDictionary
from dictionarydata import DISTRIBUTION_OPTIONS

from dictionarydata import LABEL_FONT
from dictionarydata import (HCL_PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                            DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            OUTPUT_LABEL_WIDTH1, PARAMETER_FRAME_PADX, CB_PADX)
from cmpnts_tabs.commons import commons_read_tab_vars

HCL_PARAMETERS = ['NTG', 'logResPerm',   'reservoirPressureMult', 'logWellPerm',
                  'avgWaterSaturation', 'FCO2', 'FC1', 'FC4', 'FC7Plus']

HCL_PARAMETERS_SETUP = {
    'NTG': ["Net-to-gross ratio:", 'fraction of the reservoir contributing to flow'],
    'logResPerm': ["Reservoir permeability [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B2'), 'reservoir permeability'],
    'reservoirPressureMult': ["Reservoir pressure multiplier [-]:",
                              ''.join(['factor representing state of reservoir ',
                                       'pressurization post-injection'])],
    'logWellPerm': ["Well permeability [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B2'), 'wellbore permeability'],
    'avgWaterSaturation': ["Average water saturation [-]:",
                           'average water saturation'],
    'FCO2': ["Mole fraction of CO{} [-]:".format(u'\u2082'),
             ''.join(['mole fraction of CO{} in the reservoir post-injection.\n',
                      'The four mole fraction parameters must sum up to 1']).format(u'\u2082')],
    'FC1': ["Mole fraction of methane [-]:",
            ''.join(['mole fraction of methane in the reservoir post-injection.\n',
                     'The four mole fraction parameters must sum up to 1'])],
    'FC4': ["Intermediate hydrocarbons [-]:",
            ''.join(['mole fraction of intermediate hydrocarbons in the reservoir ',
                     'post-injection.\nThe four fraction parameters must sum up to 1'])],
    'FC7Plus': ["Heavy hydrocarbons [-]:",
                ''.join(['mole fraction of heavy hydrocarbons in the reservoir ',
                         'post-injection.\nThe four mole fraction parameters ',
                         'must sum up to 1'])]}

# Set HydrocarbonLeakage parameter names and value, min, max, second value,
# mean, std, and lower and upper bounds
HCL_PARAMETER_VALUES = {
    'NTG': [0.6, 0.45, 0.95, 0.75, 0.7, 1, 0.4, 1.0],
    'logResPerm': [-13.5, -13.9, -13.1, -13.3, -13.5057, 0.2, -14.0057, -13.0057],
    'reservoirPressureMult': [1.1, 1.005, 1.195, 1.05, 1.1, 0.05, 1.0, 1.2],
    'logWellPerm': [-13.0, -16.9, -12.2, -15, -14.5057, 0.5, -17.0057, -12.0057],
    'avgWaterSaturation': [0.5, 0.48, 0.78, 0.6, 0.6315, 0.05, 0.471, 0.792],
    'FCO2': [0.55, 0.445, 0.68, 0.49, 0.5625, 0.05, 0.432, 0.693],
    'FC1': [0.05, 0.02, 0.1, 0.08, 0.0615, 0.01, 0.010, 0.113],
    'FC4': [0.05,0.012, 0.1, 0.08, 0.0605, 0.01, 0.010, 0.111],
    'FC7Plus': [0.35, 0.13, 0.48, 0.25, 0.3115, 0.05, 0.123, 0.500]}

HCL_OBSERVATIONS = ['mass_CO2_aquifer', 'mass_CO2_gas_aquifer',
                    'mass_methane_oil_aquifer', 'mass_methane_gas_aquifer',
                    'mass_oil_aquifer', 'mass_gas_aquifer']

HCL_OBSERVATIONS_SETUP = {
    'mass_CO2_aquifer': [
        'Liquid CO{} mass in aquifer [kg]'.format(u'\u2082'),
        'Enable liquid CO{} mass leaked to aquifer as output.'.format(u'\u2082')],
    'mass_CO2_gas_aquifer': [
        'Gas CO{} mass in aquifer [kg]'.format(u'\u2082'),
        'Enable gas CO{} mass leaked to aquifer as output.'.format(u'\u2082')],
    'mass_methane_oil_aquifer': [
        'Liquid methane mass in aquifer [kg]',
        'Enable liquid methane mass leaked to aquifer as output.'],
    'mass_methane_gas_aquifer': [
        'Gas methane mass in aquifer [kg]',
        'Enable gas methane mass leaked to aquifer as output.'],
    'mass_oil_aquifer': [
        'Oil mass in aquifer [kg]',
        'Enable total oil mass leaked to aquifer as output.'],
    'mass_gas_aquifer': [
        'Gas mass in aquifer [kg]',
        'Enable total gas mass leaked to aquifer as output.'],}


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""

    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'HydrocarbonLeakage', parameter_names=HCL_PARAMETERS,
        observation_names=HCL_OBSERVATIONS)

    return cmpnt_data

def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip, *args):
    """ Add widgets to the component tab."""

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('HydrocarbonLeakage')
    connectionsDictionary.append('none')

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set('none')

    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    # For parameters that can assume any of the distributions
    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(
        HCL_PARAMETER_VALUES)

    for obs_nm in HCL_OBSERVATIONS:
        componentVars[cmpnt_nm][obs_nm] = BooleanVar()
        componentVars[cmpnt_nm][obs_nm].set(0)

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Hydrocarbon Leakage Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(HCL_PARAMETERS):

        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': HCL_PARAMETER_VALUES[par_name][6],
             'upper_bound': HCL_PARAMETER_VALUES[par_name][7]},
            HCL_PARAMETERS_SETUP[par_name][0],
            HCL_PARAMETERS_SETUP[par_name][1],
            HCL_PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
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
    for ind in range(6):
        obs_nm = HCL_OBSERVATIONS[ind]
        # Create output checkbox
        cbox = tk.Checkbutton(
            outputs_frame, variable=componentVars[cmpnt_nm][obs_nm])

        output_nms_checkboxes.append(cbox)
        # Place checkbox
        output_nms_checkboxes[-1].grid(
            row=ind//2+1, column=2*(ind%2), pady=5, padx=CB_PADX, sticky='w')
        # Create output label
        output_nms_labels.append(ttk.Label(
            outputs_frame, text=HCL_OBSERVATIONS_SETUP[obs_nm][0],
            width=OUTPUT_LABEL_WIDTH1+20, anchor='w'))
        # Place label
        output_nms_labels[-1].grid(row=ind//2+1, column=2*(ind%2)+1, pady=5, sticky='w')
        # Bind checkbox to the tip
        tool_tip.bind(output_nms_checkboxes[-1],
                      HCL_OBSERVATIONS_SETUP[obs_nm][1])
