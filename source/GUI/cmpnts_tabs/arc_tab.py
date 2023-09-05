"""
Module contains several methods needed for creating tab (page) in GUI
for AnalyticalReservoir component. Methods read and write dictionaries
needed for control file interface yaml files.
"""

import tkinter as tk
from tkinter import ttk
from tkinter import StringVar
from tkinter import BooleanVar

from dictionarydata import componentVars, componentChoices
from dictionarydata import connectionsDictionary, componentTypeDictionary
from dictionarydata import DISTRIBUTION_OPTIONS

from dictionarydata import LABEL_FONT
from dictionarydata import (PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                            DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            OUTPUT_LABEL_WIDTH1, PARAMETER_FRAME_PADX, CB_PADX)

from cmpnts_tabs.locations import (add_inj_well_frame_widgets,
                                   add_obs_locs_frame_widgets,
                                   read_obs_locations_data)
from cmpnts_tabs.commons import commons_read_tab_vars


AR_PARAMETERS = [
    'logResPerm', 'reservoirPorosity', 'reservoirRadius',
    'brineDensity', 'CO2Density',
    'brineViscosity', 'CO2Viscosity', 'brineResSaturation',
    'brineCompressibility', 'injRate']

AR_PARAMETERS_SETUP = {
    'logResPerm': ["Reservoir permeability [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B2'),
                   'reservoir permeability'],
    'reservoirPorosity': ["Reservoir porosity [-]:",
                          'reservoir porosity'],
    'reservoirRadius': ["Reservoir radius [m]:",
                        'reservoir radius'],
    'brineDensity': ["Brine density [kg/m{}]:".format(u'\u00B3'),
                     'brine density'],
    'CO2Density': ["CO{} density [kg/m{}]:".format(u'\u2082', u'\u00B3'),
                   "CO{} density".format(u'\u2082')],
    'brineViscosity': ["Brine viscosity [Pa{}s]:".format(u'\u22C5'),
                       'brine viscosity'],
    'CO2Viscosity': ["CO{} viscosity [Pa{}s]:".format(u'\u2082', u'\u22C5'),
                     "CO{} viscosity".format(u'\u2082')],
    'brineResSaturation': ["Brine saturation [-]:",
                           'brine residual saturation'],
    'brineCompressibility': ["Brine compressibility [Pa{}{}]:".format(
        u'\u207B', u'\u00B9'), 'brine compressibility'],
    'injRate': ["CO{} injection rate [m{}/s]:".format(
        u'\u2082', u'\u00B3'),
                "CO{} injection rate".format(u'\u2082')]}

# Set Simple Reservoir parameters names and value, min, max, second value, mean, std, bounds
AR_PARAMETER_VALUES = {
    'logResPerm': [-13, -14, -12, -14, -14, 1, -15.3, -12],
    'reservoirPorosity': [0.15, 0.1, 0.2, 0.2, 0.2, 0.09, 0.1, 0.3],
    'reservoirRadius': [600, 500, 600, 550, 500, 50, 500, 100000],
    'brineDensity': [1000, 965, 1100, 1100, 1000, 25, 965, 1195],
    'brineViscosity': [2.535e-4, 3.0e-4, 1.0e-3, 3.0e-4,
                       1.0e-3, 1.0e-4, 2.3e-4, 1.59e-3],
    'CO2Density': [479, 500, 900, 550, 600, 50, 450, 976],
    'CO2Viscosity': [3.95e-5, 1.0e-5, 1.0e-4, 6.0e-5,
                     5.0e-5, 1.5e-6, 0.455e-6, 1.043e-4],
    'brineResSaturation': [0, 0, 0.2, 0.1, 0.2, 0.02, 0, 0.25],
    'brineCompressibility': [4.5e-12, 3.9e-12, 2.0e-11, 2.0e-11,
                             5.0e-12, 1.0e-13, 3.63e-12, 2.31e-11],
    'injRate': [0.0185, 1.0e-2, 2, 0.05, 1, 0.01, 0.0024, 3.77]}

AR_OBSERVATIONS = ['pressure', 'CO2saturation', 'pressureAve', 'mass_CO2_reservoir']

AR_OBSERVATIONS_SETUP = {
    'CO2saturation': ["CO{} saturation [-]".format(u'\u2082'),
                      'Enable vertically averaged CO{} saturation as output.'.format(u'\u2082')],
    'pressure': ["Pressure [Pa]", 'Enable pressure as output.'],
    'pressureAve': ['Averaged pressure [Pa]', 'Enable vertically averaged pressure as output.'],
    'mass_CO2_reservoir': ["CO{} mass [kg]".format(u'\u2082'),
                           'Enable mass of CO{} as output.'.format(u'\u2082')]}


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""

    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'AnalyticalReservoir',
        parameter_names=AR_PARAMETERS, observation_names=AR_OBSERVATIONS)

    cmpnt_data = read_obs_locations_data(cmpnt_data, cmpnt_nm)

    return cmpnt_data


def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip, *args):
    """ Add widgets to the component tab."""

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('AnalyticalReservoir')
    connectionsDictionary.append('none')

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set('none')

    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    # Populate dictionary
    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(AR_PARAMETER_VALUES)

    for output_key in AR_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    comp_type_label = ttk.Label(
        tab, text="Analytical Reservoir Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(AR_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': AR_PARAMETER_VALUES[par_name][6],
             'upper_bound': AR_PARAMETER_VALUES[par_name][7]},
            AR_PARAMETERS_SETUP[par_name][0],
            AR_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, tool_tip)

    # Injection well location
    inj_well_frame = tk.Frame(tab)
    inj_well_frame.grid(row=11, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)
    add_inj_well_frame_widgets(controller, cmpnt_nm, 'AnalyticalReservoir',
                               inj_well_frame, tool_tip)

    # Observation locations frame
    obs_locs_frame = tk.Frame(tab)
    obs_locs_frame.grid(row=12, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)
    add_obs_locs_frame_widgets(controller, cmpnt_nm, obs_locs_frame, tool_tip)

    outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    outputs_label.grid(row=13, column=0, sticky='w', pady=(5, 10))

    # Outputs
    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(row=14, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    # Create and place outputs widgets
    output_nms_labels = []
    output_nms_checkboxes = []
    for ind1 in range(2):
        for ind2 in range(2):
            ind = ind1*2 + ind2
            obs_nm = AR_OBSERVATIONS[ind]
            # Create output checkbox
            output_nms_checkboxes.append(tk.Checkbutton(
                outputs_frame, variable=componentVars[cmpnt_nm][obs_nm]))
            # Place checkbox
            output_nms_checkboxes[-1].grid(
                row=ind1+1, column=2*ind2, pady=5, padx=CB_PADX, sticky='w')
            # Create output label
            output_nms_labels.append(ttk.Label(
                outputs_frame, text=AR_OBSERVATIONS_SETUP[obs_nm][0],
                width=OUTPUT_LABEL_WIDTH1, anchor='w'))
            # Place label
            output_nms_labels[-1].grid(row=ind1+1, column=2*ind2+1, pady=5, sticky='w')
            # Bind checkbox to the tip
            tool_tip.bind(output_nms_checkboxes[-1],
                          AR_OBSERVATIONS_SETUP[obs_nm][1])
