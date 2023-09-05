"""
Module contains several methods needed for creating tab (page) in GUI
for DeepAlluviumAquiferML component. Methods read and write dictionaries
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
from dictionarydata import (PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                            DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            OUTPUT_LABEL_WIDTH2, PARAMETER_FRAME_PADX, CB_PADX)

from cmpnts_tabs.commons import commons_read_tab_vars


DAAML_PARAMETERS = [
    'logK_sand1', 'logK_sand2', 'logK_sand3', 'logK_caprock',
    'correlationLengthX', 'correlationLengthZ',
    'sandFraction', 'groundwater_gradient', 'leak_depth']

DAAML_PARAMETERS_SETUP = {
    'logK_sand1': [
        "Permeability of layer 1 [log{} m{}]:".format(
            u'\u2081'u'\u2080', u'\u00B2'),
        'permeability of layer 1'],
    'logK_sand2': [
        "Permeability of layer 2 [log{} m{}]:".format(
            u'\u2081'u'\u2080', u'\u00B2'),
        'permeability of layer 2'],
    'logK_sand3': [
        "Permeability of layer 3 [log{} m{}]:".format(
            u'\u2081'u'\u2080', u'\u00B2'),
        'permeability of layer 3'],
    'logK_caprock': [
        "Caprock permeability [log{} m{}]:".format(
            u'\u2081'u'\u2080', u'\u00B2'),
        'caprock permeability'],
    'correlationLengthX': ["Correlation length x [m]:",
                           'correlation length in x-direction'],
    'correlationLengthZ': ["Correlation length z [m]:",
                           'correlation length in z-direction'],
    'sandFraction': ["Sand volume fraction [-]:",
                     'sand volume fraction'],
    'groundwater_gradient': ["Groundwater gradient [-]:",
                             'groundwater gradient'],
    'leak_depth': ["Depth of leakage [m]:",
                   'depth of leakage']}

# Set Deep Alluvium Aquifer ML parameters names and value, min, max, second value,
# mean, stddev, low_bound, upper_bound
DAAML_PARAMETER_VALUES = {
    'logK_sand1': [-11.91, -12.92, -10.92, -12.1, -11.5, 1, -12.92, -10.92],
    'logK_sand2': [-11.71, -12.72, -10.72, -12.1, -11.5, 1, -12.72, -10.72],
    'logK_sand3': [-11.69, -12.7, -10.7, -12.1, -11.5, 1, -12.7, -10.7],
    'logK_caprock': [-15.70, -16.7, -14.7, -15, -15, 1, -16.7, -14.7],
    'correlationLengthX': [1098.235, 200, 2000, 750, 1000, 450, 200, 2000],
    'correlationLengthZ': [79.827, 10, 150, 50, 75, 25, 10, 150],
    'sandFraction': [0.8, 0.7, 0.9, 0.75, 0.8, 0.05, 0.7, 0.9],
    'groundwater_gradient': [0.0013, 0.001, 0.0017, 0.0015, 0.0014, 0.0001, 0.001, 0.0017],
    'leak_depth': [883.3, 424.4, 1341.5, 550, 900, 150, 424.4, 1341.5]}

DAAML_DYNAMIC_KWARGS = ['brine_rate', 'co2_rate', 'brine_mass', 'co2_mass']

DAAML_OBSERVATIONS = [
    'TDS_volume', 'Pressure_volume', 'pH_volume',
    'TDS_dx', 'Pressure_dx', 'pH_dx',
    'TDS_dy', 'Pressure_dy', 'pH_dy',
    'TDS_dz', 'Pressure_dz', 'pH_dz']

DAAML_OBSERVATIONS_SETUP = {
    'TDS_volume': [
        "Volume above baseline TDS [m{}]".format(u'\u00B3'),
        'Enable aquifer volume above baseline TDS as output.'],
    'Pressure_volume': [
         "Volume above baseline {}P [m{}]".format(u'\u0394', u'\u00B3'),
         'Enable aquifer volume above baseline {}P as output.'.format(u'\u0394')],
    'pH_volume': [
         "Volume below pH threshold [m{}]".format(u'\u00B3'),
         'Enable aquifer volume below pH threshold as output.'],
    'TDS_dx': [
         "Length x of volume above baseline TDS [m]",
         'Enable length x of aquifer volume above baseline TDS as output.'],
    'Pressure_dx': [
         "Length x of volume above baseline {}P [m]".format(u'\u0394'),
         'Enable length x of aquifer volume above baseline {}P as output.'.format(u'\u0394')],
    'pH_dx': [
         "Length x of volume below pH threshold [m]",
         'Enable length x of aquifer volume below pH threshold as output.'],
    'TDS_dy': [
         "Width y of volume above baseline TDS [m]",
         'Enable width y of aquifer volume above baseline TDS as output.'],
    'Pressure_dy': [
         "Width y of volume above baseline {}P [m]".format(u'\u0394'),
         'Enable width y of aquifer volume above baseline {}P as output.'.format(u'\u0394')],
    'pH_dy': [
         "Width y of volume below pH threshold [m]",
         'Enable width y of aquifer volume below pH threshold as output.'],
    'TDS_dz': [
         "Height z of volume above baseline TDS [m]",
         'Enable height z of aquifer volume above baseline TDS as output.'],
    'Pressure_dz': [
         "Height z of volume above baseline {}P [m]".format(u'\u0394'),
         'Enable height z of aquifer volume above baseline {}P as output.'.format(u'\u0394')],
    'pH_dz': [
         "Height z of volume below pH threshold [m]",
         'Enable height z of aquifer volume below pH threshold as output.']}


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'DeepAlluviumAquiferML', parameter_names=DAAML_PARAMETERS,
        dynamic_kwarg_names=DAAML_PARAMETERS,
        observation_names=DAAML_OBSERVATIONS, add_aquifer=True)

    return cmpnt_data

def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, toolTip,
                cnctn_nm, dyn_data, aquifer_nm, *args):
    """ Add widgets to the component tab."""

    aquifers = ['aquifer{}'.format(ind) for ind in range(
        1, componentVars['strata']['Params']['numberOfShaleLayers'].get())]

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('DeepAlluviumAquiferML')
    connectionsDictionary.append(cnctn_nm)

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set(cnctn_nm)

    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    componentVars[cmpnt_nm]['aquiferName'] = StringVar()
    componentVars[cmpnt_nm]['aquiferName'].set(aquifer_nm)

    # Populate dictionary
    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(
        DAAML_PARAMETER_VALUES)

    for output_key in DAAML_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    if 'Dynamic' in cnctn_nm:  # dynamic kwargs are provided instead of connected component
        for ind, key in enumerate(DAAML_DYNAMIC_KWARGS):
            componentVars[cmpnt_nm][key] = dyn_data[ind]

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Deep Alluvium Aquifer ML Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(DAAML_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': DAAML_PARAMETER_VALUES[par_name][6],
             'upper_bound': DAAML_PARAMETER_VALUES[par_name][7]},
            DAAML_PARAMETERS_SETUP[par_name][0],
            DAAML_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, toolTip)

    con_frame = tk.Frame(tab)
    con_frame.grid(row=10, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    connection_label = ttk.Label(
        con_frame, text="Connection:", width=PARAMETER_LABEL_WIDTH)
    connection_menu = tk.OptionMenu(
        con_frame, componentVars[cmpnt_nm]['connection'], *connections)
    connection_label.grid(
        row=0, column=0, sticky='e', padx=5)
    connection_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    connection_menu.grid(row=0, column=1, padx=5)
    toolTip.bind(connection_menu,
                 'Set connection for this component.')

    aq_name_label = ttk.Label(
        con_frame, width=PARAMETER_LABEL_WIDTH, text="Aquifer name:")
    aq_name_menu = tk.OptionMenu(
        con_frame, componentVars[cmpnt_nm]['aquiferName'], *aquifers)
    aq_name_label.grid(row=1, column=0, sticky='e', padx=5)
    aq_name_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    aq_name_menu.grid(row=1, column=1, padx=5)
    toolTip.bind(aq_name_menu,
                 'Set the aquifer name for this component.')

    # Outputs
    outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    outputs_label.grid(row=18, column=0, sticky='w', pady=(5, 10))

    outputFrame = ttk.Frame(tab)
    outputFrame.grid(row=19, column=0, sticky='w',
                     padx=PARAMETER_FRAME_PADX, columnspan=20)

    obs_labels = []
    obs_checkboxes = []
    for ind, key in enumerate(DAAML_OBSERVATIONS):
        # Create checkbox for a given observation
        obs_checkboxes.append(tk.Checkbutton(
            outputFrame, variable=componentVars[cmpnt_nm][key]))
        # Create observation label
        obs_labels.append(ttk.Label(
            outputFrame, text=DAAML_OBSERVATIONS_SETUP[key][0],
            width=OUTPUT_LABEL_WIDTH2, anchor='w'))
        # Determine location on the grid and place widgets on the frame
        row_ind = ind//3 + 1
        col_ind = 2*(ind%3)
        obs_checkboxes[-1].grid(
            row=row_ind, column=col_ind, pady=5, padx=CB_PADX, sticky='w')
        obs_labels[-1].grid(
            row=row_ind, column=col_ind+1, pady=5, sticky='w')
        # Create tip for the checkbox
        toolTip.bind(obs_checkboxes[-1], DAAML_OBSERVATIONS_SETUP[key][1])
