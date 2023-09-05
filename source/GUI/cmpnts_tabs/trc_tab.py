"""
Module contains several methods needed for creating tab (page) in GUI
for TheisReservoir component. Methods read and write dictionaries
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
                                   add_theis_inj_times_rates_frame_widgets,
                                   add_obs_locs_frame_widgets,
                                   read_obs_locations_data)
from cmpnts_tabs.commons import commons_read_tab_vars


TR_PARAMETERS = [
    'logResPerm', 'reservoirPorosity', 'compressibility',
    'brineDensity', 'CO2Density', 'brineViscosity', 'initialPressure']

TR_PARAMETERS_SETUP = {
    'initialPressure': ["Initial reservoir pressure [Pa]:",
                          'initial pressure'],
    'logResPerm': ["Reservoir permeability [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B2'),
                   'reservoir permeability'],
    'reservoirPorosity': ["Reservoir porosity [-]:",
                          'reservoir porosity'],
    'compressibility': ["Brine compressibility [Pa{}{}]:".format(
        u'\u207B', u'\u00B9'), 'brine compressibility'],
    'CO2Density': ["CO{} density [kg/m{}]:".format(u'\u2082', u'\u00B3'),
                   "CO{} density".format(u'\u2082')],
    'brineDensity': ["Brine density [kg/m{}]:".format(u'\u00B3'),
                     'brine density'],
    'brineViscosity': ["Brine viscosity [Pa{}s]:".format(u'\u22C5'),
                       'brine viscosity']}

# Set Theis Reservoir parameters names and value, min, max, second value, mean, std, bounds
TR_PARAMETER_VALUES = {
    'initialPressure': [1.0e+6, 1.0e+5, 5.0e+6, 2.5e+06, 2.5e+06, 1.0e+05,
                        8.0e+4, 1.0e+7],
    'logResPerm': [-12, -13, -11, -11.5, -12, 1, -14, -9],
    'reservoirPorosity': [0.3, 0.1, 0.33, 0.25, 0.2, 0.05, 0.01, 1.0],
    'compressibility': [1.0e-10, 5.0e-11, 5.0e-9, 1.0e-9, 2.5e-10, 1.0e-11,
                        1.0e-11, 1.0e-6],
    'CO2Density': [479, 400, 600, 500, 450, 25, 100, 1500],
    'brineDensity': [1000, 950, 1200, 1100, 1000, 25, 900, 1500],
    'brineViscosity': [2.535e-3, 3.0e-4, 3.0e-3, 1.0e-3,
                       1.0e-4, 5.0e-5, 1.0e-4, 5.0e-3]}

TR_OBSERVATIONS = ['pressure', 'CO2saturation']

TR_OBSERVATIONS_SETUP = {
    'pressure': ["Pressure [Pa]", 'Enable pressure as output.'],
    'CO2saturation': ["CO{} saturation [-]".format(u'\u2082'),
                      'Enable vertically averaged CO{} saturation as output.'.format(
                          u'\u2082')]}


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""

    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'TheisReservoir',
        parameter_names=TR_PARAMETERS, observation_names=TR_OBSERVATIONS)

    cmpnt_data = read_obs_locations_data(cmpnt_data, cmpnt_nm)

    return cmpnt_data


def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip, *args):
    """ Add widgets to the component tab."""

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('TheisReservoir')
    connectionsDictionary.append('none')

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set('none')

    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    # Populate dictionary
    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(TR_PARAMETER_VALUES)

    for output_key in TR_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    comp_type_label = ttk.Label(
        tab, text="Theis Reservoir Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(TR_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': TR_PARAMETER_VALUES[par_name][6],
             'upper_bound': TR_PARAMETER_VALUES[par_name][7]},
            TR_PARAMETERS_SETUP[par_name][0],
            TR_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, tool_tip)

    # Injection well location
    inj_well_frame = tk.Frame(tab)
    inj_well_frame.grid(row=11, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)
    add_inj_well_frame_widgets(controller, cmpnt_nm, 'TheisReservoir', inj_well_frame, tool_tip)

    # Injection times and rates
    inj_times_rates_frame = tk.Frame(tab)
    inj_times_rates_frame.grid(row=12, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)
    add_theis_inj_times_rates_frame_widgets(controller, cmpnt_nm,
                                            inj_times_rates_frame, tool_tip)

    # Observation locations frame
    obs_locs_frame = tk.Frame(tab)
    obs_locs_frame.grid(row=13, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)
    add_obs_locs_frame_widgets(controller, cmpnt_nm, obs_locs_frame, tool_tip)

    outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    outputs_label.grid(row=14, column=0, sticky='w', pady=(5, 10))

    # Outputs
    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(row=15, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    # Create and place outputs widgets
    output_nms_labels = []
    output_nms_checkboxes = []
    for ind in range(len(TR_OBSERVATIONS)):
        obs_nm = TR_OBSERVATIONS[ind]
        # Create output checkbox
        output_nms_checkboxes.append(tk.Checkbutton(
            outputs_frame, variable=componentVars[cmpnt_nm][obs_nm]))
        # Place checkbox
        output_nms_checkboxes[-1].grid(
            row=1, column=2*ind, pady=5, padx=CB_PADX, sticky='w')
        # Create output label
        output_nms_labels.append(ttk.Label(
            outputs_frame, text=TR_OBSERVATIONS_SETUP[obs_nm][0],
            width=OUTPUT_LABEL_WIDTH1, anchor='w'))
        # Place label
        output_nms_labels[-1].grid(row=1, column=2*ind+1, pady=5, sticky='w')
        # Bind checkbox to the tip
        tool_tip.bind(output_nms_checkboxes[-1],
                      TR_OBSERVATIONS_SETUP[obs_nm][1])
