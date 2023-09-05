"""
Module contains several methods needed for creating tab (page) in GUI
for AtmosphericROM component. Methods read and write dictionaries
needed for control file interface yaml files.
"""

import tkinter as tk
from tkinter import ttk
from tkinter import StringVar
from tkinter import BooleanVar

from dictionarydata import componentVars, componentChoices
from dictionarydata import connectionsDictionary, componentTypeDictionary
from dictionarydata import DISTRIBUTION_OPTIONS, connections

from dictionarydata import LABEL_FONT
from dictionarydata import (PARAMETER_LABEL_WIDTH,
                            DISTRIBUTION_MENU_WIDTH,
                            DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            OUTPUT_LABEL_WIDTH1, PARAMETER_FRAME_PADX,
                            CB_PADX)

from cmpnts_tabs.commons import commons_read_tab_vars


ATM_PARAMETERS = ['T_amb', 'P_amb', 'V_wind', 'C0_critical', 'T_source']

ATM_PARAMETERS_SETUP = {
    'T_amb': ["Ambient temperature [{}C]:".format(u'\u00B0'),
              'ambient temperature'],
    'P_amb': ["Ambient pressure [atm]:",
              'ambient pressure'],
    'V_wind': ["Wind velocity [m/s]:",
               'wind velocity'],
    'C0_critical': ["Critical concentration [-]:",
                    'critical concentration'],
    'T_source': ["Released CO{} temperature [{}C]:".format(u'\u2082', u'\u00B0'),
                 "released CO{} temperature".format(u'\u2082')]}

ATM_OBSERVATIONS = [
    'outflag', 'num_sources', 'x_new', 'y_new', 'critical_distance']

# Set Atmospheric ROM parameters names and value, min, max, second value, mean, std
ATM_PARAMETER_VALUES = {'T_amb': [15, 5, 40, 25, 25, 5, 5, 40],
                        'P_amb': [1.0, 0.7, 1.08, 0.8, 0.8, 0.05, 0.7, 1.08],
                        'V_wind': [5, 1.0e-10, 20, 10, 10, 3, 1.0e-10, 20],
                        'C0_critical': [0.01, 0.002, 0.1, 0.05, 0.05, 0.025, 0.002, 0.1],
                        'T_source': [15, 5, 50, 25, 30, 10, 5, 50]}

ATM_KWARGS = ['x_receptor', 'y_receptor']

ATM_KWARG_VALUES = {'x_receptor': 10.0, 'y_receptor': 10.0}

ATM_DYNAMIC_KWARGS = ['co2_leakrate']


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'AtmosphericROM', parameter_names=ATM_PARAMETERS,
        dynamic_kwarg_names=ATM_DYNAMIC_KWARGS,
        observation_names=ATM_OBSERVATIONS)

    # Read values of entered keyword arguments
    for key in ATM_KWARGS:
        cmpnt_data[key] = []
        data = componentVars[cmpnt_nm][key].get().split(',')
        for val in data:
            cmpnt_data[key].append(float(val.strip()))

    return cmpnt_data

def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip,
                cnctn_nm, dyn_data, *args):
    """ Add widgets to the component tab."""

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('AtmosphericROM')
    connectionsDictionary.append(cnctn_nm)

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set(cnctn_nm)

    # Populate dictionary with parameters
    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(
        ATM_PARAMETER_VALUES)

    # Setup default kwargs values
    for key_arg in ATM_KWARGS:
        componentVars[cmpnt_nm][key_arg] = StringVar()
        componentVars[cmpnt_nm][key_arg].set(ATM_KWARG_VALUES[key_arg])

    # Setup default observation selection
    for output_key in ATM_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    if 'Dynamic' in cnctn_nm:  # dynamic kwargs are provided instead of connected component
        for ind, key in enumerate(ATM_DYNAMIC_KWARGS):
            componentVars[cmpnt_nm][key] = dyn_data[ind]

    componentVars[cmpnt_nm]['componentName'] = StringVar()
    componentVars[cmpnt_nm]['componentType'] = StringVar()
    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Atmospheric ROM Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(ATM_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': ATM_PARAMETER_VALUES[par_name][6],
             'upper_bound': ATM_PARAMETER_VALUES[par_name][7]},
            ATM_PARAMETERS_SETUP[par_name][0],
            ATM_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, tool_tip)

    x_receptor = tk.Frame(tab)
    x_receptor.grid(row=6, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    x_receptor_label = ttk.Label(
        x_receptor, text="x-coordinate of receptor [m]:",
        width=PARAMETER_LABEL_WIDTH)
    x_receptor_txtField = tk.Entry(
        x_receptor, textvariable=componentVars[cmpnt_nm]['x_receptor'])
    x_receptor_txtField.config(width=DISTRIBUTION_MENU_WIDTH+6)
    x_receptor.labelText = "x-coordinate of receptor [m]:"
    x_receptor.component = componentVars[cmpnt_nm]['componentName']
    x_receptor.toolTipText = "x-coordinate of receptor"
    x_receptor.text = 'x_receptor'

    x_receptor_label.grid(row=0, column=0, sticky='w', padx=5)
    x_receptor_txtField.grid(row=0, column=1, padx=5, pady=8)
    tool_tip.bind(x_receptor_txtField,
                 'Enter a list of x-coordinates for receptors.')
    controller.setvar(name='x_receptor', value=x_receptor)

    y_receptor = tk.Frame(tab)
    y_receptor.grid(row=7, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    y_receptor_label = ttk.Label(
        y_receptor,
        text="y-coordinate of receptor [m]:", width=PARAMETER_LABEL_WIDTH)
    y_receptor_txtField = tk.Entry(
        y_receptor, textvariable=componentVars[cmpnt_nm]['y_receptor'])
    y_receptor_txtField.config(width=DISTRIBUTION_MENU_WIDTH+6)
    y_receptor.labelText = "y-coordinate of receptor [m]:"
    y_receptor.component = componentVars[cmpnt_nm]['componentName']
    y_receptor.toolTipText = "y-coordinate of receptor"
    y_receptor.text = 'y_receptor'

    y_receptor_label.grid(row=0, column=0, sticky='w', padx=5)
    y_receptor_txtField.grid(row=0, column=1, padx=5, pady=8)
    tool_tip.bind(y_receptor_txtField,
                 'Enter a list of y-coordinates for receptors.')
    controller.setvar(name='y_receptor', value=y_receptor)

    con_frame = tk.Frame(tab)
    con_frame.grid(row=8, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    connection_label = ttk.Label(
        con_frame, text="Connection:", width=PARAMETER_LABEL_WIDTH)
    connection_menu = tk.OptionMenu(
        con_frame, componentVars[cmpnt_nm]['connection'], *connections)
    connection_label.grid(row=0, column=0, sticky='w', padx=5)
    connection_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    connection_menu.grid(row=0, column=1, padx=5)
    tool_tip.bind(connection_menu,
                 'Set connection for this component.')

    # Outputs
    outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    outputs_label.grid(row=9, column=0, sticky='w', pady=(5, 10))

    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(row=10, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    num_sources_label = ttk.Label(
        outputs_frame, text="Number of sources [-]",
        width=OUTPUT_LABEL_WIDTH1, anchor='w')
    num_sources_checkbox = tk.Checkbutton(
        outputs_frame, variable=componentVars[cmpnt_nm]['num_sources'])
    num_sources_label.grid(row=0, column=1, pady=5, sticky='w')
    num_sources_checkbox.grid(
        row=0, column=0, pady=5, padx=CB_PADX, sticky='w')
    tool_tip.bind(num_sources_checkbox,
                 'Enable number of sources as output.')

    x_new_label = ttk.Label(
        outputs_frame, text="Source x-coord [m]",
        width=OUTPUT_LABEL_WIDTH1, anchor='w')
    x_new_checkbox = tk.Checkbutton(
        outputs_frame, variable=componentVars[cmpnt_nm]['x_new'])
    x_new_label.grid(row=0, column=3, pady=5, sticky='w')
    x_new_checkbox.grid(
        row=0, column=2, pady=5, padx=CB_PADX, sticky='w')
    tool_tip.bind(x_new_checkbox,
                 'Enable x-coordinates of leakage sources as output.')

    y_new_label = ttk.Label(
        outputs_frame, text="Source y-coord [m]",
        width=OUTPUT_LABEL_WIDTH1, anchor='w')
    y_new_checkbox = tk.Checkbutton(
        outputs_frame, variable=componentVars[cmpnt_nm]['y_new'])
    y_new_label.grid(row=0, column=5, pady=5, sticky='w')
    y_new_checkbox.grid(
        row=0, column=4, pady=5, padx=CB_PADX, sticky='w')
    tool_tip.bind(y_new_checkbox,
                 'Enable y-coordinates of leakage sources as output.')

    outflag_label = ttk.Label(
        outputs_frame, text="Outflag [-]", width=OUTPUT_LABEL_WIDTH1, anchor='w')
    outflag_checkbox = tk.Checkbutton(
        outputs_frame, variable=componentVars[cmpnt_nm]['outflag'])
    outflag_label.grid(row=1, column=1, pady=5, sticky='w')
    outflag_checkbox.grid(
        row=1, column=0, pady=5, padx=CB_PADX, sticky='w')
    tool_tip.bind(outflag_checkbox, 'Enable outflag as output.')

    critical_distance_label = ttk.Label(
        outputs_frame, text="Critical distance [m]",
        width=OUTPUT_LABEL_WIDTH1, anchor='w')
    critical_distance_checkbox = tk.Checkbutton(
        outputs_frame, variable=componentVars[cmpnt_nm]['critical_distance'])
    critical_distance_label.grid(
        row=1, column=3, pady=5, sticky='w')
    critical_distance_checkbox.grid(
        row=1, column=2, pady=5, padx=CB_PADX, sticky='w')
    tool_tip.bind(critical_distance_checkbox,
                 'Enable critical distance as output.')
