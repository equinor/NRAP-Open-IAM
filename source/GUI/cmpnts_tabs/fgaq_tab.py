"""
Module contains several methods needed for creating tab (page) in GUI
for FutureGen2Aquifer component. Methods read and write dictionaries
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
                            OUTPUT_LABEL_WIDTH2, PARAMETER_FRAME_PADX,
                            CB_PADX)

from cmpnts_tabs.commons import commons_read_tab_vars


FGAQ_PARAMETERS = ['por', 'log_permh', 'log_aniso', 'rel_vol_frac_calcite']

FGAQ_PARAMETERS_SETUP = {
    'por': ["Aquifer porosity [-]:",
            'aquifer porosity'],
    'log_permh': [
        "Horizontal permeability [log{} m{}]:".format(
            u'\u2081'u'\u2080', u'\u00B2'),
        'horizontal permeability'],
    'log_aniso': ["Anisotropy ratio [log{}]:".format(u'\u2081'u'\u2080'),
                  'anisotropy ratio'],
    'rel_vol_frac_calcite': ["Volume fraction of calcite [-]:",
                             'relative volume fraction of calcite']}

# Set FutureGen 2 Aquifer parameters names and value, min, max, second value,
# mean, std, low_bound, upper_bound
FGAQ_PARAMETER_VALUES = {
    'por': [0.15, 0.02, 0.2, 0.18, 0.1, 0.005, 0.02, 0.2],
    'log_permh': [-13, -14, -11, -12, -12, 0.2, -14, -11],
    'log_aniso': [0.3, 0, 3, 0.25, 0.4, 0.05, 0, 3],
    'rel_vol_frac_calcite': [0.3, 0, 1, 0.5, 0.5, 0.05, 0, 1]}

FGAQ_DYNAMIC_KWARGS = ['brine_rate', 'co2_rate', 'brine_mass', 'co2_mass']

FGAQ_OBSERVATIONS = [
    'TDS_volume', 'TDS_dx', 'TDS_dy', 'TDS_dz',
    'pH_volume', 'pH_dx', 'pH_dy', 'pH_dz',
    'Pressure_volume', 'Pressure_dx', 'Pressure_dy', 'Pressure_dz',
    'Dissolved_CO2_volume', 'Dissolved_CO2_dx', 'Dissolved_CO2_dy', 'Dissolved_CO2_dz']

FGAQ_OBSERVATIONS_SETUP = {
    'TDS_volume': ["Volume above baseline TDS [m{}]".format(u'\u00B3'),
                   'Enable aquifer volume above baseline TDS as output.'],
    'Pressure_volume': [
        "Volume above baseline {}P [m{}]".format(u'\u0394', u'\u00B3'),
        'Enable aquifer volume above baseline {}P as output.'.format(u'\u0394')],
    'pH_volume': ["Volume below pH threshold [m{}]".format(u'\u00B3'),
                  'Enable aquifer volume below pH threshold as output.'],
    'TDS_dx': ["Length x of volume above baseline TDS [m]",
               'Enable length x of aquifer volume above baseline TDS as output.'],
    'Pressure_dx': [
        "Length x of volume above baseline {}P [m]".format(u'\u0394'),
        'Enable length x of aquifer volume above baseline {}P as output.'.format(u'\u0394')],
    'pH_dx': ["Length x of volume below pH threshold [m]",
              'Enable length x of aquifer volume below pH threshold as output.'],
    'TDS_dy': ["Width y of volume above baseline TDS [m]",
               'Enable width y of aquifer volume above baseline TDS as output.'],
    'Pressure_dy': [
        "Width y of volume above baseline {}P [m]".format(u'\u0394'),
        'Enable width y of aquifer volume above baseline {}P as output.'.format(u'\u0394')],
    'pH_dy': ["Width y of volume below pH threshold [m]",
              'Enable width y of aquifer volume below pH threshold as output.'],
    'TDS_dz': ["Height z of volume above baseline TDS [m]",
               'Enable height z of aquifer volume above baseline TDS as output.'],
    'Pressure_dz': [
        "Height z of volume above baseline {}P [m]".format(u'\u0394'),
        'Enable height z of aquifer volume above baseline {}P as output.'.format(u'\u0394')],
    'pH_dz': ["Height z of volume below pH threshold [m]",
              'Enable height z of aquifer volume below pH threshold as output.'],
    'Dissolved_CO2_volume': [
        "Volume above baseline CO{} [m{}]".format(u'\u2082', u'\u00B3'),
        ''.join(['Enable aquifer volume above baseline ',
                 'dissolved CO{} as output.']).format(u'\u2082')],
    'Dissolved_CO2_dx': [
        "Length x of volume above baseline CO{} [m]".format(u'\u2082'),
        ''.join(['Enable length x of aquifer volume above baseline ',
                 'dissolved CO{} as output.']).format(u'\u2082')],
    'Dissolved_CO2_dy': [
        "Width y of volume above baseline CO{} [m]".format(u'\u2082'),
        ''.join(['Enable width y of aquifer volume above baseline ',
                 'dissolved CO{} as output.']).format(u'\u2082')],
    'Dissolved_CO2_dz': [
        "Height z of volume above baseline CO{} [m]".format(u'\u2082'),
        ''.join(['Enable height z of aquifer volume above baseline ',
                 'dissolved CO{} as output.']).format(u'\u2082')]}


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'FutureGen2Aquifer', parameter_names=FGAQ_PARAMETERS,
        dynamic_kwarg_names=FGAQ_DYNAMIC_KWARGS,
        observation_names=FGAQ_OBSERVATIONS, add_aquifer=True)

    return cmpnt_data

def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip,
                cnctn_nm, dyn_data, aquifer_nm, *args):
    """ Add widgets to the component tab."""

    aquifers = ['aquifer{}'.format(ind) for ind in range(
        1, componentVars['strata']['Params']['numberOfShaleLayers'].get())]

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('FutureGen2Aquifer')
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
        FGAQ_PARAMETER_VALUES)

    for output_key in FGAQ_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    if 'Dynamic' in cnctn_nm:  # dynamic kwargs are provided instead of connected component
        for ind, key in enumerate(FGAQ_DYNAMIC_KWARGS):
            componentVars[cmpnt_nm][key] = dyn_data[ind]

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="FutureGen 2 Aquifer Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(FGAQ_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': FGAQ_PARAMETER_VALUES[par_name][6],
             'upper_bound': FGAQ_PARAMETER_VALUES[par_name][7]},
            FGAQ_PARAMETERS_SETUP[par_name][0],
            FGAQ_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, tool_tip)

    # Connection frame
    con_frame = tk.Frame(tab)
    con_frame.grid(row=14, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    connection_label = ttk.Label(
        con_frame, text="Connection:", width=PARAMETER_LABEL_WIDTH)
    connection_menu = tk.OptionMenu(
        con_frame, componentVars[cmpnt_nm]['connection'], *connections)
    connection_label.grid(row=0, column=0, sticky='e', padx=5)
    connection_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    connection_menu.grid(row=0, column=1, padx=5)
    tool_tip.bind(connection_menu,
                  'Set connection for this component.')

    aq_name_label = ttk.Label(
        con_frame, width=PARAMETER_LABEL_WIDTH, text="Aquifer name:")
    aq_name_menu = tk.OptionMenu(
        con_frame, componentVars[cmpnt_nm]['aquiferName'], *aquifers)
    aq_name_label.grid(row=1, column=0, sticky='e', padx=5)
    aq_name_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    aq_name_menu.grid(row=1, column=1, padx=5)
    tool_tip.bind(aq_name_menu,
                  'Set the aquifer name for this component.')

    # Outputs
    outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    outputs_label.grid(row=18, column=0, sticky='w', pady=(5, 10))

    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(
        row=19, column=0, sticky='w', padx=PARAMETER_FRAME_PADX, columnspan=20)

    # Create and place outputs widgets
    output_nms_labels = []
    output_nms_checkboxes = []

    for ind1 in range(2):
        for ind2 in range(8):
            ind = ind1*8 + ind2
            obs_nm = FGAQ_OBSERVATIONS[ind]
            # Create output checkbox
            output_nms_checkboxes.append(
                tk.Checkbutton(outputs_frame, variable=componentVars[cmpnt_nm][obs_nm]))
            # Place checkbox
            output_nms_checkboxes[-1].grid(
                row=ind2+1, column=2*ind1, pady=5, padx=CB_PADX, sticky='w')
            # Create output label
            output_nms_labels.append(
                ttk.Label(outputs_frame, text=FGAQ_OBSERVATIONS_SETUP[obs_nm][0],
                          width=OUTPUT_LABEL_WIDTH2, anchor='w'))
            # Place label
            output_nms_labels[-1].grid(row=ind2+1, column=2*ind1+1, pady=5, sticky='w')
            # Bind checkbox to the tip
            tool_tip.bind(output_nms_checkboxes[-1],
                          FGAQ_OBSERVATIONS_SETUP[obs_nm][1])
