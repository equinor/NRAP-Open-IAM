"""
Module contains several methods needed for creating tab (page) in GUI
for GenericAquifer component. Methods read and write dictionaries
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
                            OUTPUT_LABEL_WIDTH3, PARAMETER_FRAME_PADX,
                            CB_PADX)

from cmpnts_tabs.commons import commons_read_tab_vars


GA_PARAMETERS = ['por', 'log_permh', 'log_aniso',
                 'aquifer_salinity', 'reservoir_salinity',
                 'dissolved_salt_threshold', 'dissolved_co2_threshold']

GA_PARAMETERS_SETUP = {
    'por': ["Aquifer porosity [-]:",
            'aquifer porosity'],
    'log_permh': [
        "Horizontal permeability [log{} m{}]:".format(
            u'\u2081'u'\u2080', u'\u00B2'),
        'horizontal permeability'],
    'log_aniso': ["Anisotropy ratio [log{}]:".format(u'\u2081'u'\u2080'),
                  'anisotropy ratio'],
    'aquifer_salinity': ['Aquifer salt mass fraction [-]',
                         'salt mass fraction in aquifer water '],
    'reservoir_salinity': ['Reservoir salt mass fraction [-]',
                           'salt mass fraction in reservoir water'],
    'dissolved_salt_threshold': ['Salt mass fraction threshold [-]',
                                 'threshold for salt mass fraction'],
    'dissolved_co2_threshold': ['CO{} mass fraction threshold'.format(u'\u2082'),
                                'threshold for CO{} mass fraction'.format(u'\u2082')]}

# Set Generic Aquifer parameters names and value, min, max, second value, mean,
# std, low_bound, upper_bound
GA_PARAMETER_VALUES = {
    'por': [0.18, 0.05, 0.2, 0.1, 0.1, 0.005, 0.02, 0.25],
    'log_permh': [-13, -12, -11, -12, -12, 0.1, -14, -10],
    'log_aniso': [0.3, 0.2, 2.2, 0.25, 0.4, 0.05, 0, 3],
    'aquifer_salinity': [0.005, 0.005, 0.01, 0.007, 0.01, 0.0005, 0, 0.015],
    'reservoir_salinity': [0.03, 0.02, 0.04, 0.02, 0.03, 0.005, 0.015, 0.05],
    'dissolved_salt_threshold': [0.02, 0.1, 0.5, 0.25, 0.5, 0.01, 0, 1],
    'dissolved_co2_threshold': [0.01, 0.1, 0.75, 0.1, 0.5, 0.02, 0, 1]}

GA_DYNAMIC_KWARGS = ['brine_mass', 'co2_mass']

GA_OBSERVATIONS = [
    'Dissolved_CO2_volume', 'Dissolved_salt_volume',
    'Dissolved_CO2_dr', 'Dissolved_salt_dr',
    'Dissolved_CO2_dz', 'Dissolved_salt_dz',
    'Dissolved_CO2_mass_fraction', 'Dissolved_salt_mass_fraction']

GA_OBSERVATIONS_SETUP = {
    'Dissolved_CO2_volume': [
        "Volume above CO{} threshold [m{}]".format(u'\u2082', u'\u00B3'),
        ''.join(['Enable aquifer volume above dissolved CO{} threshold ',
                 'as output.']).format(u'\u2082')],
    'Dissolved_CO2_dr': [
        "Radius r of volume above CO{} threshold [m]".format(u'\u2082'),
        ''.join(['Enable radius r of aquifer volume above dissolved CO{} ',
                 'threshold as output.']).format(u'\u2082')],
    'Dissolved_CO2_dz': [
        "Height z of volume above CO{} threshold [m]".format(u'\u2082'),
        ''.join(['Enable height z of aquifer volume above dissolved CO{} ',
                 'threshold as output.']).format(u'\u2082')],
    'Dissolved_salt_volume': [
        "Volume above salt mass threshold [m{}]".format(u'\u00B3'),
        'Enable aquifer volume above salt mass threshold as output.'],
    'Dissolved_salt_dr': [
        "Radius r of volume above salt mass threshold [m]",
        'Enable radius r of aquifer volume above salt mass threshold as output.'],
    'Dissolved_salt_dz': [
        "Height z of volume above salt mass threshold [m]",
        'Enable height z of aquifer volume above salt mass threshold as output.'],
    'Dissolved_CO2_mass_fraction': [
        'Mass fraction of CO{} in aquifer [-]'.format(u'\u2082'),
        'Enable mass fraction of CO{} in aquifer pore water as output'.format(u'\u2082')],
    'Dissolved_salt_mass_fraction': [
        'Mass fraction of salt in aquifer [-]',
        'Enable mass fraction of salt in aquifer pore water as output']}


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'GenericAquifer', parameter_names=GA_PARAMETERS,
        dynamic_kwarg_names=GA_DYNAMIC_KWARGS,
        observation_names=GA_OBSERVATIONS, add_aquifer=True)

    return cmpnt_data

def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip,
                cnctn_nm, dyn_data, aquifer_nm, *args):
    """ Add widgets to the component tab."""

    aquifers = ['aquifer{}'.format(ind) for ind in range(
        1, componentVars['strata']['Params']['numberOfShaleLayers'].get())]

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('GenericAquifer')
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
        GA_PARAMETER_VALUES)

    for output_key in GA_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    if 'Dynamic' in cnctn_nm:  # dynamic kwargs are provided instead of connected component
        for ind, key in enumerate(GA_DYNAMIC_KWARGS):
            componentVars[cmpnt_nm][key] = dyn_data[ind]

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Generic Aquifer Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(GA_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': GA_PARAMETER_VALUES[par_name][6],
             'upper_bound': GA_PARAMETER_VALUES[par_name][7]},
            GA_PARAMETERS_SETUP[par_name][0],
            GA_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, tool_tip)

    # Connection frame
    con_frame = tk.Frame(tab)
    con_frame.grid(row=8, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

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
    outputs_label.grid(row=9, column=0, sticky='w', pady=(5, 10))

    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(
        row=10, column=0, sticky='w', padx=PARAMETER_FRAME_PADX, columnspan=20)

    # Create and place outputs widgets
    output_nms_labels = []
    output_nms_checkboxes = []

    for ind in range(8):
        ind1 = ind // 2  # row index
        ind2 = ind % 2   # column index

        obs_nm = GA_OBSERVATIONS[ind]
        # Create output checkbox
        output_nms_checkboxes.append(
            tk.Checkbutton(outputs_frame, variable=componentVars[cmpnt_nm][obs_nm]))
        # Place checkbox
        output_nms_checkboxes[-1].grid(
            row=ind1, column=2*ind2, pady=5, padx=CB_PADX, sticky='w')
        # Create output label
        output_nms_labels.append(
            ttk.Label(outputs_frame, text=GA_OBSERVATIONS_SETUP[obs_nm][0],
                      width=OUTPUT_LABEL_WIDTH3, anchor='w'))
        # Place label
        output_nms_labels[-1].grid(row=ind1, column=2*ind2+1, pady=5, sticky='w')
        # Bind checkbox to the tip
        tool_tip.bind(output_nms_checkboxes[-1],
                      GA_OBSERVATIONS_SETUP[obs_nm][1])
