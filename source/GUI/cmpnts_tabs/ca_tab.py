"""
Module contains several methods needed for creating tab (page) in GUI
for CarbonateAquifer component. Methods read and write dictionaries
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
                            OUTPUT_LABEL_WIDTH2,
                            PARAMETER_FRAME_PADX, CB_PADX)

from cmpnts_tabs.commons import commons_read_tab_vars


CA_PARAMETERS = [
    'rmin', 'mean_perm', 'perm_var', 'corr_len', 'aniso', 'hyd_grad',
    'calcite_ssa', 'organic_carbon', 'benzene_kd', 'benzene_decay',
    'nap_kd', 'nap_decay', 'phenol_kd', 'phenol_decay', 'cl']

CA_PARAMETERS_SETUP = {
    'rmin': ["Distance between leaks [m]:",
             'distance between leaks'],
    'mean_perm': [
        "Permeability mean [log{} m{}]:".format(u'\u2081'u'\u2080', u'\u00B2'),
        'permeability mean'],
    'perm_var': [
        "Permeability variance [log{} m{}]:".format(u'\u2081'u'\u2080', u'\u2074'),
        'permeability variance'],
    'corr_len': ["Correlation length [m]:",
                 'correlation length'],
    'aniso': ["Anisotropy factor [-]:",
              'anisotropy factor'],
    'hyd_grad': ["Hydraulic gradient [-]:",
                 'hydraulic gradient'],
    'calcite_ssa': ["Calcite surface [m{}/g]:".format(u'\u00B2'),
                    'calcite surface'],
    'organic_carbon': ["Organic carbon [-]:",
                       'organic carbon'],
    'benzene_kd': [
        "Benzene distribution [log{} K_oc]:".format(u'\u2081'u'\u2080'),
        'benzene distribution'],
    'benzene_decay': ["Benzene decay [log{} day]:".format(u'\u2081'u'\u2080'),
                      'benzene decay constant'],
    'nap_kd': [
        "Naphthalene distribution [log{} K_oc]:".format(u'\u2081'u'\u2080'),
        'naphthalene distribution'],
    'nap_decay': [
        "Naphthalene decay [log{} day]:".format(u'\u2081'u'\u2080'),
        'naphthalene decay constant'],
    'phenol_kd': [
        "Phenol distribution [log{} K_oc]:".format(u'\u2081'u'\u2080'),
        'phenol distribution'],
    'phenol_decay': ["Phenol decay [log{} day]:".format(u'\u2081'u'\u2080'),
                     'phenol decay constant'],
    'cl': ["Brine salinity [log{} molality]:".format(u'\u2081'u'\u2080'),
           'brine salinity']}

# Set Carbonate Aquifer parameter names and value, min, max, second value, mean,
# stddev, lower_bound, upper_bound
CA_PARAMETER_VALUES = {
    'rmin': [15, 0, 100, 10, 20, 10, 0, 100],
    'perm_var': [0.9535, 0.17, 1.89, 0.8, 1.4, 0.3, 1.7e-2, 1.89],
    'corr_len': [2.475, 1, 3.95, 2, 2.5, 1, 1, 3.95],
    'aniso': [25.1, 1.1, 49.1, 20, 28, 6, 1.1, 49.1],
    'mean_perm': [-12.5, -13.8, -10.3, -12, -12, 1, -13.8, -10.3],
    'hyd_grad': [9.59e-03, 2.88e-4, 1.89e-2, 5.0e-3, 5.0e-3, 4.0e-3, 2.88e-4, 1.89e-2],
    'calcite_ssa': [5.5e-03, 0, 1.0e-2, 2.0e-3, 5.0e-3, 2.0e-3, 0, 1.0e-2, 0, 1.0e-2],
    'organic_carbon': [5.5e-03, 0, 1.0e-2, 2.0e-3, 5.0e-3, 2.0e-3, 0, 1.0e-2],
    'benzene_kd': [1.61, 1.49, 1.73, 1.55, 1.6, 0.07, 1.49, 1.73],
    'benzene_decay': [0.595, 0.15, 2.84, 0.37, 1.3, 0.5, 0.15, 2.84],
    'nap_kd': [2.98, 2.78, 3.18, 2.85, 2.9, 0.05, 2.78, 3.18],
    'nap_decay': [0.595, -0.85, 2.04, 0.76, 0.52, 0.1, -0.85, 2.04],
    'phenol_kd': [1.35, 1.21, 1.48, 1.27, 1.35, 0.02, 1.21, 1.48],
    'phenol_decay': [0.42, -1.22, 2.06, 0.78, 0.78, 0.5, -1.22, 2.06],
    'cl': [0.776, 0.1, 6.025, 2.3, 0.75, 0.08, 0.1, 6.025]}

CA_DYNAMIC_KWARGS = ['brine_rate', 'co2_rate', 'brine_mass', 'co2_mass']

CA_OBSERVATIONS = ['dx', 'dy', 'Flux', 'pH_volume', 'TDS_volume', 'As_volume',
                   'Pb_volume', 'Cd_volume', 'Ba_volume', 'Benzene_volume',
                   'Naphthalene_volume', 'Phenol_volume']

CA_OBSERVATIONS_SETUP = {
    'dx': ["Length x of impacted volume [m]",
           'Enable length x of impacted aquifer volume as output.'],
    'dy': ["Width y of impacted volume [m]",
           'Enable width y of impacted aquifer volume as output.'],
    'Flux': ["CO{} flux to atmosphere [kg/s]".format(u'\u2082'),
             'Enable CO{} flux to atmosphere as output.'.format(u'\u2082')],
    'pH_volume': ["Volume below pH threshold [m{}]".format(u'\u00B3'),
                  'Enable aquifer volume below pH threshold as output.'],
    'TDS_volume': ["Volume above baseline TDS [m{}]".format(u'\u00B3'),
                   'Enable aquifer volume above baseline TDS as output.'],
    'As_volume': ["Volume above arsenic threshold [m{}]".format(u'\u00B3'),
                  'Enable aquifer volume above arsenic threshold as output.'],
    'Pb_volume': ["Volume above lead threshold [m{}]".format(u'\u00B3'),
                  'Enable aquifer volume above lead threshold as output.'],
    'Cd_volume': ["Volume above cadmium threshold [m{}]".format(u'\u00B3'),
                  'Enable aquifer volume above cadmium threshold as output'],
    'Ba_volume': ["Volume above barium threshold [m{}]".format(u'\u00B3'),
                  'Enable aquifer volume above barium threshold as output.'],
    'Benzene_volume': ["Volume above benzene threshold [m{}]".format(u'\u00B3'),
                       'Enable aquifer volume above benzene threshold as output.'],
    'Naphthalene_volume': [
        "Volume above naphthalene threshold [m{}]".format(u'\u00B3'),
        'Enable aquifer volume above naphthalene threshold as output.'],
    'Phenol_volume': ["Volume above phenol threshold [m{}]".format(u'\u00B3'),
                      'Enable aquifer volume above phenol threshold as output.']}


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'CarbonateAquifer', parameter_names=CA_PARAMETERS,
        dynamic_kwarg_names=CA_DYNAMIC_KWARGS,
        observation_names=CA_OBSERVATIONS, add_aquifer=True)

    if componentVars[cmpnt_nm]['ithresh'].get().find('MCL'):
        cmpnt_data['Parameters']['ithresh'] = 2
    if componentVars[cmpnt_nm]['ithresh'].get().find('No'):
        cmpnt_data['Parameters']['ithresh'] = 1

    if componentVars[cmpnt_nm]['logf'].get() == 'Linear':
        cmpnt_data['Parameters']['logf'] = 0
    if componentVars[cmpnt_nm]['logf'].get() == 'Log':
        cmpnt_data['Parameters']['logf'] = 1

    return cmpnt_data

def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, toolTip,
                cnctn_nm, dyn_data, aquifer_nm, *args):
    """ Add widgets to the component tab."""

    aquifers = ['aquifer{}'.format(ind) for ind in range(
        1, componentVars['strata']['Params']['numberOfShaleLayers'].get())]

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('CarbonateAquifer')
    connectionsDictionary.append(cnctn_nm)

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set(cnctn_nm)

    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    componentVars[cmpnt_nm]['ithresh'] = StringVar()
    componentVars[cmpnt_nm]['ithresh'].set('No-Impact')
    componentVars[cmpnt_nm]['logf'] = StringVar()
    componentVars[cmpnt_nm]['logf'].set('Linear')
    componentVars[cmpnt_nm]['aquiferName'] = StringVar()
    componentVars[cmpnt_nm]['aquiferName'].set(aquifer_nm)

    if 'Dynamic' in cnctn_nm:  # dynamic kwargs are provided instead of connected component
        for ind, key in enumerate(CA_DYNAMIC_KWARGS):
            componentVars[cmpnt_nm][key] = dyn_data[ind]

    componentVars[cmpnt_nm]['Params'] = (
        controller.populate_params_dict(CA_PARAMETER_VALUES))

    for output_key in CA_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Carbonate Aquifer Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    ithresh = tk.Frame(tab)
    ithresh.grid(row=1, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    ithresh_label = ttk.Label(
        ithresh, text="Threshold [-]:", width=PARAMETER_LABEL_WIDTH)
    ithresh_menu = tk.OptionMenu(
        ithresh, componentVars[cmpnt_nm]['ithresh'], *['MCL', 'No-Impact'])
    ithresh_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    ithresh_label.grid(row=0, column=0, sticky='e', padx=5)
    ithresh_menu.grid(row=0, column=1, padx=5)
    toolTip.bind(ithresh_menu,
                 'Select which threshold you wish to utilize.')

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(CA_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+2, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': CA_PARAMETER_VALUES[par_name][6],
             'upper_bound': CA_PARAMETER_VALUES[par_name][7]},
            CA_PARAMETERS_SETUP[par_name][0],
            CA_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, toolTip)


    logf = tk.Frame(tab)
    logf.grid(row=17, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    logf_label = ttk.Label(
        logf, text="Log transform:", width=PARAMETER_LABEL_WIDTH)
    logf_menu = tk.OptionMenu(
        logf, componentVars[cmpnt_nm]['logf'], *['Linear', 'Log'])

    logf_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    logf_label.grid(row=0, column=0, sticky='w', padx=5)
    logf_menu.grid(row=0, column=1, padx=5)
    toolTip.bind(logf_menu,
                 ''.join(['Select the type of transform you wish ',
                          'to use for this component model.']))

    connection_label = ttk.Label(
        logf, text="Connection:", width=PARAMETER_LABEL_WIDTH)
    connection_menu = tk.OptionMenu(
        logf, componentVars[cmpnt_nm]['connection'], *connections)
    connection_label.grid(row=1, column=0, sticky='w', padx=5)
    connection_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    connection_menu.grid(row=1, column=1, padx=5)
    toolTip.bind(connection_menu,
                 'Set connection for this component.')

    aquifer_name_label = ttk.Label(
        logf, width=PARAMETER_LABEL_WIDTH, text="Aquifer name:")
    aquifer_name_menu = tk.OptionMenu(
        logf, componentVars[cmpnt_nm]['aquiferName'], *aquifers)
    aquifer_name_label.grid(row=2, column=0, sticky='w', padx=5)
    aquifer_name_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    aquifer_name_menu.grid(row=2, column=1, padx=5)
    toolTip.bind(aquifer_name_menu,
                 'Set the aquifer name for this component.')

    outputs_label = ttk.Label(
        tab, text="Outputs", font=LABEL_FONT)
    outputs_label.grid(row=18, column=0, sticky='w', pady=(5, 10))

    # Outputs
    outputFrame = ttk.Frame(tab)
    outputFrame.grid(row=19, column=0, sticky='w',
                     padx=PARAMETER_FRAME_PADX, columnspan=20)

    # Create and place outputs widgets
    output_nms_labels = []
    output_nms_checkboxes = []

    for ind1 in range(4):
        for ind2 in range(3):
            ind = ind1*3 + ind2
            obs_nm = CA_OBSERVATIONS[ind]
            # Create output checkbox
            output_nms_checkboxes.append(
                tk.Checkbutton(outputFrame, variable=componentVars[cmpnt_nm][obs_nm]))
            # Place checkbox
            output_nms_checkboxes[-1].grid(
                row=ind1+1, column=2*ind2, pady=5, padx=CB_PADX, sticky='w')
            # Create output label
            output_nms_labels.append(
                ttk.Label(outputFrame, text=CA_OBSERVATIONS_SETUP[obs_nm][0],
                          width=OUTPUT_LABEL_WIDTH2, anchor='w'))
            # Place label
            output_nms_labels[-1].grid(row=ind1+1, column=2*ind2+1, pady=5, sticky='w')
            # Bind checkbox to the tip
            toolTip.bind(output_nms_checkboxes[-1],
                         CA_OBSERVATIONS_SETUP[obs_nm][1])
