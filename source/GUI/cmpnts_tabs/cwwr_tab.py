"""
Module contains several methods needed to create a tab (page) in GUI
for CementedWellboreWR component. Methods read and write dictionaries
needed for control file interface yaml files.
"""
import os
import sys

import tkinter as tk
from tkinter import ttk
from tkinter import StringVar, IntVar, BooleanVar

# Save location of GUI folder
GUI_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(GUI_DIR)

from dictionarydata import componentVars, componentChoices
from dictionarydata import connectionsDictionary, componentTypeDictionary
from dictionarydata import DISTRIBUTION_OPTIONS, connections

from dictionarydata import LABEL_FONT
from dictionarydata import (PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                            DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            OUTPUT_LABEL_WIDTH1, PARAMETER_FRAME_PADX, CB_PADX)

from cmpnts_tabs.locations import read_locations_data, add_wellbore_frame_widgets
from cmpnts_tabs.commons import commons_read_tab_vars


CW_WR_PARAMETERS = ['logWellPerm', 'logThiefPerm', 'wellRadius']

CW_WR_PARAMETERS_SETUP = {
    'logWellPerm': ["Well permeability [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B2'), 'well permeability'],
    'logThiefPerm': ["Thief zone permeability [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B2'), 'thief zone permeability'],
    'wellRadius': ["Well radius [m]:",
                   'well radius']}

# Set Cemented Wellbore WR parameters names and value, min, max, second value,
# mean, std, bounds
CW_WR_PARAMETER_VALUES = {
    'logWellPerm': [-13, -13.5, -10.5, -12, -12, 1, -13.95, -10.1],
    'logThiefPerm': [-12.5, -13.9, -12.1, -13.5, -13, 0.5, -13.986, -12.023],
    'wellRadius': [0.05, 0.025, 0.25, 0.1, 0.1, 0.04, 0.025, 0.25]}

CW_WR_DYNAMIC_KWARGS = ['pressure', 'CO2saturation']

CW_WR_OBSERVATIONS = ['CO2_aquifer1', 'brine_aquifer1', 'mass_CO2_aquifer1',
                      'CO2_aquifer2', 'brine_aquifer2', 'mass_CO2_aquifer2',
                      'CO2_atm', 'brine_atm']

CW_WR_OBSERVATIONS_SETUP = {
    'CO2_aquifer1': [
        "CO{} aquifer 1 [kg/s]".format(u'\u2082'),
        'Enable CO{} leakage rates to aquifer 1 as output.'.format(u'\u2082')],
    'brine_aquifer1': ["Brine aquifer 1 [kg/s]",
                       'Enable brine leakage rates to aquifer 1 as output.'],
    'mass_CO2_aquifer1': [
        "Mass CO{} aquifer 1 [kg]".format(u'\u2082'),
        'Enable CO{} mass leaked to aquifer 1 as output.'.format(u'\u2082')],
    'CO2_aquifer2': [
        "CO{} aquifer 2 [kg/s]".format(u'\u2082'),
        'Enable CO{} leakage rates to aquifer 2 as output.'.format(u'\u2082')],
    'brine_aquifer2': ["Brine aquifer 2 [kg/s]",
                       'Enable brine leakage rates to aquifer 2 as output.'],
    'mass_CO2_aquifer2': [
        "Mass CO{} aquifer 2 [kg]".format(u'\u2082'),
        'Enable CO{} mass leaked to aquifer 2 as output.'.format(u'\u2082')],
    'CO2_atm': [
        "CO{} atm [kg/s]".format(u'\u2082'),
        'Enable CO{} leakage rates to atmosphere as output.'.format(u'\u2082')],
    'brine_atm': ["Brine atm [kg/s]",
                  'Enable brine leakage rates to atmosphere as output.']}

def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""

    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'CementedWellboreWR', parameter_names=CW_WR_PARAMETERS,
        dynamic_kwarg_names=CW_WR_DYNAMIC_KWARGS,
        observation_names=CW_WR_OBSERVATIONS,
        add_number=True)

    # Read information about locations associated with component
    read_locations_data(cmpnt_data, cmpnt_nm)

    return cmpnt_data

def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip,
                cnctn_nm, dyn_data, *args):
    """ Add widgets to the component tab."""

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('CementedWellboreWR')
    connectionsDictionary.append(cnctn_nm)

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set(cnctn_nm)

    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    componentVars[cmpnt_nm]['number'] = IntVar()
    componentVars[cmpnt_nm]['number'].set(1)

    componentVars[cmpnt_nm]['useRandomLocDomain'] = BooleanVar()
    componentVars[cmpnt_nm]['useRandomLocDomain'].set(0)

    if 'Dynamic' in cnctn_nm:  # dynamic kwargs are provided instead of connected component
        for ind, key in enumerate(CW_WR_DYNAMIC_KWARGS):
            componentVars[cmpnt_nm][key] = dyn_data[ind]

    # Populate dictionary
    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(CW_WR_PARAMETER_VALUES)

    for output_key in CW_WR_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Cemented Wellbore (WR) Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(CW_WR_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': CW_WR_PARAMETER_VALUES[par_name][6],
             'upper_bound': CW_WR_PARAMETER_VALUES[par_name][7]},
            CW_WR_PARAMETERS_SETUP[par_name][0],
            CW_WR_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, tool_tip)

    # Number of wellbores defined by the same parameters
    number_frame = tk.Frame(tab)
    number_frame.grid(row=7, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    connection_label = ttk.Label(
        number_frame, text="Connection:", width=PARAMETER_LABEL_WIDTH)
    connection_menu = tk.OptionMenu(
        number_frame, componentVars[cmpnt_nm]['connection'], *connections)
    connection_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    connection_label.grid(row=0, column=0, sticky='w', padx=5)
    connection_menu.grid(row=0, column=1, padx=5)
    tool_tip.bind(connection_menu,
                  'Set connection for this component.')

    number_label = ttk.Label(
        number_frame, text='Number of wellbores:', width=PARAMETER_LABEL_WIDTH)
    number_spinbox = tk.Spinbox(
        number_frame, from_=1, to=1000, textvariable=componentVars[cmpnt_nm]['number'])
    number_label.grid(row=1, column=0, padx=5, sticky='w')
    number_spinbox.grid(row=1, column=1, padx=5, pady=2, sticky='ew')
    tool_tip.bind(number_spinbox,
                  ''.join(['Set the total number of wells for this wellbore ',
                           'component including random locations.']))

    # Wellbore locations frame
    well_locs_frame = tk.Frame(tab)
    well_locs_frame.grid(row=8, column=0, sticky='w',
                         padx=PARAMETER_FRAME_PADX, pady=(5, 10))
    add_wellbore_frame_widgets(controller, cmpnt_nm, well_locs_frame, tool_tip)

    # Outputs
    outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    outputs_label.grid(row=9, column=0, sticky='w', pady=(5, 10))

    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(row=10, column=0, columnspan=10, sticky='w',
                       padx=PARAMETER_FRAME_PADX)

    # Create and place outputs widgets
    output_nms_labels = []
    output_nms_checkboxes = []

    ind2_range = [3, 3, 2]
    for ind1 in range(3):
        for ind2 in range(ind2_range[ind1]):
            ind = ind1*3 + ind2
            obs_nm = CW_WR_OBSERVATIONS[ind]
            # Create output checkbox
            output_nms_checkboxes.append(
                tk.Checkbutton(outputs_frame, variable=componentVars[cmpnt_nm][obs_nm]))
            # Place checkbox
            output_nms_checkboxes[-1].grid(
                row=ind1+1, column=2*ind2, pady=5, padx=CB_PADX, sticky='w')
            # Create output label
            output_nms_labels.append(
                ttk.Label(outputs_frame, text=CW_WR_OBSERVATIONS_SETUP[obs_nm][0],
                          width=OUTPUT_LABEL_WIDTH1, anchor='w'))
            # Place label
            output_nms_labels[-1].grid(row=ind1+1, column=2*ind2+1, pady=5, sticky='w')
            # Bind checkbox to the tip
            tool_tip.bind(output_nms_checkboxes[-1],
                          CW_WR_OBSERVATIONS_SETUP[obs_nm][1])
