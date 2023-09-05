"""
Module contains several methods needed for creating tab (page) in GUI
for AlluviumAquiferLF component. Methods read and write dictionaries
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


AALF_PARAMETERS = [
    'logK_sand', 'logK_clay', 'correlationLengthX', 'correlationLengthZ',
    'sandFraction', 'NaMolality', 'PbMolality', 'benzeneMolality', 'tMitigation']

AALF_PARAMETERS_SETUP = {
    'logK_sand': [
        "Permeability of sand layer [log{} m{}]:".format(
            u'\u2081'u'\u2080', u'\u00B2'),
        'permeability of sand layer'],
    'logK_clay': [
        "Permeability of clay layer [log{} m{}]:".format(
            u'\u2081'u'\u2080', u'\u00B2'),
        'permeability of clay layer'],
    'correlationLengthX': ["Correlation length x [m]:",
                           'correlation length in x-direction'],
    'correlationLengthZ': ["Correlation length z [m]:",
                           'correlation length in z-direction'],
    'sandFraction': ["Sand volume fraction [-]:",
                     'sand volume fraction'],
    'NaMolality': ["Sodium molality [log{} mol/kg]".format(u'\u2081'u'\u2080'),
                   'sodium molality'],
    'PbMolality': ["Lead molality [log{} mol/kg]".format(u'\u2081'u'\u2080'),
                   'lead molality'],
    'benzeneMolality': ["Benzene molality [log{} mol/kg]".format(u'\u2081'u'\u2080'),
                        'benzene molality'],
    'tMitigation': ["Mitigation time [years]:",
                    'mitigation time']}

# Set Alluvium Aquifer LF parameters names and value, min, max, second value,
# mean, stddev, low_bound, upper_bound
AALF_PARAMETER_VALUES = {
    'logK_sand': [-13, -14, -10, -12, -12, 1, -14, -10],
    'logK_clay': [-16, -18, -15, -17, -16.5, 1.5, -18, -10],
    'correlationLengthX': [1100, 200, 2500, 800, 1200, 300, 200, 2500],
    'correlationLengthZ': [15, 0.5, 25, 10, 12, 5, 0.5, 25],
    'sandFraction': [0.8, 0.6, 0.9, 0.7, 0.7, 0.05, 0.6, 0.9],
    'NaMolality': [0.5, -3, 1, -1, 0.5, 0.05, -3, 1],
    'PbMolality': [-6.5, -9, -6, -7, -7, 1, -9.1836, -5.6836],
    'benzeneMolality': [-7, -10, -6, -8, -8, 1, -10, -6],
    'tMitigation': [10, 1, 100, 15, 12, 5, 1, 100]}

AALF_DYNAMIC_KWARGS = ['brine_rate', 'co2_rate', 'brine_mass', 'co2_mass']

AALF_OBSERVATIONS = [
    'TDS_volume', 'pH_volume', 'As_volume', 'Ba_volume', 'Cd_volume',
    'Pb_volume', 'Benzene_volume', 'Naphthalene_volume', 'Phenol_volume']

AALF_OBSERVATIONS_SETUP = {
    'TDS_volume': [
        "Volume above baseline TDS [m{}]".format(u'\u00B3'),
        'Enable aquifer volume above baseline TDS as output.'],
    'pH_volume': [
         "Volume below pH threshold [m{}]".format(u'\u00B3'),
         'Enable aquifer volume below pH threshold as output.'],
    'As_volume': [
         "Volume above arsenic threshold [m{}]".format(u'\u00B3'),
         'Enable aquifer volume above arsenic threshold as output.'],
    'Ba_volume': [
         "Volume above barium threshold [m{}]".format(u'\u00B3'),
         'Enable aquifer volume above barium threshold as output.'],
    'Cd_volume': [
         "Volume above cadmium threshold [m{}]".format(u'\u00B3'),
         'Enable aquifer volume above cadmium threshold as output.'],
    'Pb_volume': [
         "Volume above lead threshold [m{}]".format(u'\u00B3'),
         'Enable aquifer volume above lead threshold as output.'],
    'Benzene_volume': [
         "Volume above benzene threshold [m{}]".format(u'\u00B3'),
         'Enable aquifer volume above benzene threshold as output.'],
    'Naphthalene_volume': [
         "Volume above naphthalene threshold [m{}]".format(u'\u00B3'),
         'Enable aquifer volume above naphthalene threshold as output.'],
    'Phenol_volume': [
         "Volume above phenol threshold [m{}]".format(u'\u00B3'),
         'Enable aquifer volume above phenol threshold as output.']}


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'AlluviumAquiferLF', parameter_names=AALF_PARAMETERS,
        dynamic_kwarg_names=AALF_DYNAMIC_KWARGS,
        observation_names=AALF_OBSERVATIONS, add_aquifer=True)

    return cmpnt_data

def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, toolTip,
                cnctn_nm, dyn_data, aquifer_nm, *args):
    """ Add widgets to the component tab."""

    aquifers = ['aquifer{}'.format(ind) for ind in range(
        1, componentVars['strata']['Params']['numberOfShaleLayers'].get())]

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('AlluviumAquiferLF')
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
        AALF_PARAMETER_VALUES)

    for output_key in AALF_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    if 'Dynamic' in cnctn_nm:  # dynamic kwargs are provided instead of connected component
        for ind, key in enumerate(AALF_DYNAMIC_KWARGS):
            componentVars[cmpnt_nm][key] = dyn_data[ind]

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Alluvium Aquifer LF Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(AALF_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': AALF_PARAMETER_VALUES[par_name][6],
             'upper_bound': AALF_PARAMETER_VALUES[par_name][7]},
            AALF_PARAMETERS_SETUP[par_name][0],
            AALF_PARAMETERS_SETUP[par_name][1],
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
    connection_label.grid(row=0, column=0, sticky='e', padx=5)
    connection_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    connection_menu.grid(row=0, column=1, padx=5)
    toolTip.bind(connection_menu, 'Set connection for this component.')

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
    for ind, key in enumerate(AALF_OBSERVATIONS):
        # Create checkbox for a given observation
        obs_checkboxes.append(tk.Checkbutton(
            outputFrame, variable=componentVars[cmpnt_nm][key]))
        # Create observation label
        obs_labels.append(ttk.Label(
            outputFrame, text=AALF_OBSERVATIONS_SETUP[key][0],
            width=OUTPUT_LABEL_WIDTH2, anchor='w'))
        # Determine location on the grid and place widgets on the frame
        row_ind = ind//3 + 1
        col_ind = 2*(ind%3)
        obs_checkboxes[-1].grid(
            row=row_ind, column=col_ind, pady=5, padx=CB_PADX, sticky='w')
        obs_labels[-1].grid(
            row=row_ind, column=col_ind+1, pady=5, sticky='w')
        # Create tip for the checkbox
        toolTip.bind(obs_checkboxes[-1], AALF_OBSERVATIONS_SETUP[key][1])
