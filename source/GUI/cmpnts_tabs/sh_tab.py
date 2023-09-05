"""
Module contains several methods needed for creating tab (page) in GUI
for SealHorizon component. Methods read and write dictionaries
needed for control file interface yaml files.
"""

import os

import tkinter as tk
from tkinter import ttk
from tkinter import StringVar, IntVar, BooleanVar

from dictionarydata import componentVars, componentChoices
from dictionarydata import connectionsDictionary, componentTypeDictionary
from dictionarydata import DISTRIBUTION_OPTIONS, SPEC_DISTRIBUTION_OPTIONS

from dictionarydata import LABEL_FONT
from dictionarydata import (PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                            DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            OUTPUT_LABEL_WIDTH1,
                            PARAMETER_FRAME_PADX, CB_PADX)
from cmpnts_tabs.locations import add_cell_locs_frame_widgets, disable_cell_locs_frame_widgets
from cmpnts_tabs.commons import commons_read_tab_vars


SH_CELL_PARAMETERS = ['area', 'thickness', 'permeability', 'baseDepth',
                      'entryPressure', 'influence']

SH_PARAMETERS = {
    1: ['influenceModel', 'aveBaseDepth', 'aveBasePressure', 'aveTemperature', 'salinity',
        'staticDepth', 'staticPressure', 'brineDensity', 'CO2Density', 'brineViscosity',
        'CO2Viscosity', 'CO2Solubility', 'brineResSaturation', 'CO2ResSaturation',
        'permRatio'],  # relative model setup follows
    3: ['totalEffect', 'rateEffect'],
    4: ['reactivity', 'carbonateContent', 'clayContent']}    # clay type setup follows

CLAY_TYPES = ['smectite', 'illite', 'chlorite']

RELATIVE_MODEL_TYPES = ['LET', 'BC']

BC_MODEL_PARAMETERS = ['lambda']

LET_MODEL_PARAMETERS = ['wetting1', 'wetting2', 'wetting3',
                        'nonwet1', 'nonwet2', 'nonwet3',
                        'capillary1', 'capillary2', 'capillary3', 'maxCapillary']

ALL_SH_PARAMETERS = [val for ind in [1, 3, 4] for val in SH_PARAMETERS[ind]]+\
    LET_MODEL_PARAMETERS+BC_MODEL_PARAMETERS

SH_SPECIAL_PARAMETERS = SH_PARAMETERS[3]+SH_PARAMETERS[4]+BC_MODEL_PARAMETERS+LET_MODEL_PARAMETERS

INFLUENCE_MODEL_TYPES = [0, 1, 2]

LOC_ARG_NAMES = ['coordx', 'coordy']

# Set Seal Horizon parameter names and value, min, max, second value,
# mean, std, min and max bounds
SH_CELL_PARAMETERS_VALUES = {
    'area': [100, 100, 1000, 120, 100, 5, 1.0, 260000],
    'thickness': [100, 50, 250, 75, 100, 25, 10, 1000],
    'permeability': [1.0e-20, 1.0e-21, 1.0e-18, 1.0e-19,
                     1.0e-20, 1.0e-21, 1.0e-22, 1.0e-16],
    'baseDepth': [1007, 900, 1100, 950, 1200, 200, 800, 9500],
    'influence': [1, 0, 1, 0.9, 0.8, 0.05, 0, 1],
    'entryPressure': [5000, 4500, 5500, 5300, 5500, 500, 100, 2000000]}

SH_PARAMETER_VALUES = {
    'aveBaseDepth': [1130, 900, 1200, 1100, 1200, 200, 800, 9500],
    'aveBasePressure': [3.2e+7, 3.0e+7, 4.0e+7, 3.5e+7,
                        3.0e+7, 1.0e+5, 1.0e+6, 6.0e+7],
    'aveTemperature': [50, 45, 80, 60, 50, 10, 31, 180],
    'salinity': [2.0e+4, 1.0e+4, 3.0e+4, 1.5e+4,
                 3.0e+4, 1.0e+3, 0, 8.0e+4],
    'staticDepth': [1000, 500, 1200, 800, 1000, 100, 80, 9500],
    'staticPressure': [1.0e+7, 9.0e+6, 1.3e+7, 1.1e+7,
                       2.0e+7, 1.0e+6, 1.0e+6, 6.0e+7],
    'brineDensity': [1007, 900, 1000, 950, 950, 10, 880, 1080],
    'CO2Density': [583, 450, 700, 600, 550, 25, 93, 1050],
    'brineViscosity': [2.535e-4, 1.9e-4, 3.0e-4, 2.0e-4,
                       2.0e-4, 1.0e-5, 1.5e-4, 1.6e-3],
    'CO2Viscosity': [4.387e-5, 2.0e-5, 8.0e-5, 5.0e-5,
                     3.0e-5, 1.0e-6, 1.8e-5, 1.4e-4],
    'CO2Solubility': [1, 0.2, 1.4, 0.8, 1, 0.1, 0, 2.04],
    'brineResSaturation': [0.2, 0.1, 0.3, 0.25, 0.1, 0.02, 0.01, 0.35],
    'CO2ResSaturation': [0.01, 0.01, 0.3, 0.02, 0.1, 0.01, 0.01, 0.35],
    'permRatio': [1, 0.9, 1.2, 1.1, 0.9, 0.02, 0.1, 1.5],
    'totalEffect': [0.1, 0.05, 1.5, 0.5, 1, 0.01, 0.01, 200],
    'rateEffect': [0.1, 0.05, 0.6, 0.5, 0.6, 0.01, 0.01, 0.65],
    'reactivity': [0, 0, 2, 1, 2, 0.25, 0, 10],
    'carbonateContent': [8, 5, 20, 10, 5, 1, 0, 100],
    'clayContent': [60, 10, 75, 50, 60, 10, 0, 100],
    'lambda': [2.5, 1, 4, 3, 2.5, 0.1, 0, 5],
    'wetting1': [1, 0.5, 3, 2, 2, 0.25, 0.5, 5],
    'wetting2': [10, 2, 15, 5, 10, 1, 1, 30],
    'wetting3': [1.25, 1, 2.5, 2, 1.5, 0.1, 0.2, 3],
    'nonwet1': [1.05, 0.5, 3, 2, 2, 0.25, 0.5, 5],
    'nonwet2': [10, 2, 15, 5, 10, 1, 1, 30],
    'nonwet3': [1.25, 1, 2.5, 2, 1.5, 0.1, 0, 3],
    'capillary1': [0.2, 0.1, 0.5, 0.3, 0.5, 0.05, 0.01, 5],
    'capillary2': [2.8, 0.5, 3.5, 2, 2.5, 0.1, 0.01, 30],
    'capillary3': [0.43, 0.1, 1, 0.5, 1, 0.01, 0.01, 3],
    'maxCapillary': [1.0e+7, 9.8e+6, 5.0e+7, 1.2e+7, 1.1e+7, 1.0e+5, 100, 2.0e+8]}

SH_CELL_PARAMETERS_SETUP = {
    'area': ['Cell area [m{}]:'.format(u'\u00B2'), 'cell area'],
    'status': ['Cell status [-]:', 'cell status'],
    'thickness': ['Cell thickness [m]:', 'cell thickness'],
    'permeability': ['Cell permeability [m{}]:'.format(u'\u00B2'),
                     'cell permeability'],
    'baseDepth': ['Seal base depth [m]:',
                  'depth to the base of the seal'],
    'influence': ['Influence factor:',
                  'permeability influence factor'],
    'entryPressure': ['Entry pressure [Pa]:', 'entry pressure']}

SH_PARAMETERS_SETUP = {
    'aveBaseDepth': ['Average seal base depth [m]:', 'average seal base depth'],
    'aveBasePressure': ['Average seal base pressure [Pa]:',
                        'average seal base pressure'],
    'aveTemperature': ['Average temperature [{}C]:'.format(u'\u00B0'),
                       'average temperature'],
    'salinity': ['Salinity [ppm]:', 'salinity'],
    'staticDepth': ['Reference depth [m]:', 'reference depth'],
    'staticPressure': ['Reference pressure [Pa]:', 'reference pressure'],
    'CO2Density': ["CO{} density [kg/m{}]:".format(u'\u2082', u'\u00B3'),
                   "CO{} density".format(u'\u2082')],
    'CO2Viscosity': ["CO{} viscosity [Pa{}s]:".format(u'\u2082', u'\u22C5'),
                     "CO{} viscosity".format(u'\u2082')],
    'brineDensity': ["Brine density [kg/m{}]:".format(u'\u00B3'),
                     'brine density'],
    'brineViscosity': ["Brine viscosity [Pa{}s]:".format(u'\u22C5'),
                       'brine viscosity'],
    'CO2Solubility': ["CO{} solubility [kg/kg]:".format(u'\u2082'),
                      "CO{} solubility".format(u'\u2082')],
    'brineResSaturation': ["Brine saturation [-]:",
                           'brine residual saturation'],
    'CO2ResSaturation': ["CO{} saturation [-]:".format(u'\u2082'),
                         'CO{} residual saturation'.format(u'\u2082')],
    'permRatio': ['Permeability ratio [-]:', 'permeability ratio'],
    'totalEffect': ['Total change [-]:', 'total change in permeability factor'],
    'rateEffect': ['Rate factor [-]:', 'rate factor for permeability change'],
    'reactivity': ['Reactivity [-]:', 'reactivity of time model'],
    'carbonateContent': ['Carbonate content [%]:', 'carbonate content of rock matrix'],
    'clayContent': ['Clay content [%]:', 'clay mineral content of rock'],
    'clayType': ['Clay type:', 'Select clay type.'],
    'relativeModel': ['Relative model:', 'Select relative permeability model.'],
    'lambda': ['Lambda [-]:', 'lambda parameter in Brooks-Corey model'],
    'wetting1': ['Parameter L (WP) [-]:', 'wetting phase parameter L'],
    'wetting2': ['Parameter E (WP) [-]:', 'wetting phase parameter E'],
    'wetting3': ['Parameter T (WP) [-]:', 'wetting phase parameter T'],
    'nonwet1': ['Parameter L (NWP) [-]:', 'nonwetting phase parameter L'],
    'nonwet2': ['Parameter E (NWP) [-]:', 'nonwetting phase parameter E'],
    'nonwet3': ['Parameter T (NWP) [-]:', 'nonwetting phase parameter T'],
    'capillary1': ['Parameter L (CP) [-]:',
                   'LET-model parameter L for capillary pressure'],
    'capillary2': ['Parameter E (CP) [-]:',
                   'LET-model parameter E for capillary pressure'],
    'capillary3': ['Parameter E (CP) [-]:',
                   'LET-model parameter E for capillary pressure'],
    'maxCapillary': ['Maximum capillary pressure [Pa]:',
                     'maximum capillary pressure']}

SH_OBSERVATIONS = ['CO2_aquifer', 'mass_CO2_aquifer',
                   'CO2_aquifer_total', 'mass_CO2_aquifer_total',
                   'brine_aquifer', 'mass_brine_aquifer',
                   'brine_aquifer_total', 'mass_brine_aquifer_total']

SH_OBSERVATIONS_SETUP = {
    'CO2_aquifer': [
        "CO{} aquifer [kg/s]".format(u'\u2082'),
        'Enable individual CO{} leakage rates to aquifer as output.'.format
            (u'\u2082')],
    'brine_aquifer': [
        "Brine aquifer [kg/s]",
        'Enable individual brine leakage rates to aquifer as output.'],
    'mass_CO2_aquifer': [
        'Mass CO{} aquifer [kg]'.format(u'\u2082'),
        'Enable local CO{} mass leaked to aquifer as output.'.format(u'\u2082')],
    'mass_brine_aquifer': [
        'Mass brine aquifer [kg]',
        'Enable local brine mass leaked to aquifer as output.'],
    'CO2_aquifer_total': [
        "CO{} aquifer total [kg/s]".format(u'\u2082'),
        'Enable total CO{} leakage rates to aquifer as output.'.format
            (u'\u2082')],
    'brine_aquifer_total': [
        "Brine aquifer total [kg/s]",
        'Enable total brine leakage rates to aquifer as output.'],
    'mass_CO2_aquifer_total': [
        'Mass CO{} aquifer total [kg]'.format(u'\u2082'),
        'Enable total CO{} mass leaked to aquifer as output.'.format(u'\u2082')],
    'mass_brine_aquifer_total': [
        'Mass brine aquifer total [kg]',
        'Enable total brine mass leaked to aquifer as output.']}

SH_KWARGS = ['clayType', 'relativeModel']

SH_DYNAMIC_KWARGS = ['pressure', 'CO2saturation']


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'SealHorizon', parameter_names=ALL_SH_PARAMETERS,
        dynamic_kwarg_names=SH_DYNAMIC_KWARGS, observation_names=SH_OBSERVATIONS)

    # Get cell parameters
    cmpnt_data['Cells'] = {}
    cmpnt_data['Cells']['Number'] = componentVars[cmpnt_nm]['Cells']['Number'].get()
    cmpnt_data['Cells']['Locations'] = {}
    if componentVars[cmpnt_nm]['Cells']['Locations']['fileinput'].get():  # == 1
        cmpnt_data['Cells']['Locations']['file'] = componentVars[cmpnt_nm]['Cells'][
            'Locations']['filename'].get()
    else:
        for nm in LOC_ARG_NAMES:
            loc_data = componentVars[cmpnt_nm]['Cells']['Locations'][nm].get().split(',')
            cmpnt_data['Cells']['Locations'][nm] = []
            for val in loc_data:
                cmpnt_data['Cells']['Locations'][nm].append(float(val.strip()))

    for cell_par_nm in SH_CELL_PARAMETERS:
        cmpnt_data['Cells'][cell_par_nm] = {}

    # Additional non-scalar parameters/keyword arguments
    for kw_nm in ['relativeModel', 'clayType']:
        cmpnt_data['Parameters'][kw_nm] = componentVars[cmpnt_nm][kw_nm].get()

    # # TODO DO NOT DELETE. Might be needed later to connect seal horizon component
    # # to aquifer components
    # cmpnt_data['AquiferName'] = componentVars[cmpnt_nm]['aquiferName'].get()

    return cmpnt_data

def add_widgets(app, tab, cmpnt_nm, cmpnt_type, tool_tip,
                cnctn_nm, dyn_data, aquifer_nm, *args):
    """ Add widgets to the component tab."""

    # # TODO DO NOT DELETE. Might be needed later to connect seal horizon component
    # # to aquifer components
    # aquifers = ['aquifer{}'.format(ind) for ind in range(
    #     1, componentVars['strata']['Params']['numberOfShaleLayers'].get())]

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('SealHorizon')
    connectionsDictionary.append(cnctn_nm)

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set(cnctn_nm)

    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    # # TODO DO NOT DELETE. Might be needed later to connect seal horizon component
    # # to aquifer components
    # componentVars[cmpnt_nm]['aquiferName'] = StringVar()
    # componentVars[cmpnt_nm]['aquiferName'].set(aquifer_nm)

    # Relative model and clay type variables
    componentVars[cmpnt_nm]['relativeModel'] = StringVar()
    componentVars[cmpnt_nm]['relativeModel'].set('LET')
    componentVars[cmpnt_nm]['clayType'] = StringVar()
    componentVars[cmpnt_nm]['clayType'].set(CLAY_TYPES[0]) # 'smectite' is the first clay type

    # Cells setup variables
    componentVars[cmpnt_nm]['Cells'] = app.populate_params_dict(
        SH_CELL_PARAMETERS_VALUES,
        distr_options=SPEC_DISTRIBUTION_OPTIONS,
        options={1: ['distribution', 'ordvalues', 'filename'], 2: ['value']})

    componentVars[cmpnt_nm]['Cells']['Number'] = IntVar()
    componentVars[cmpnt_nm]['Cells']['Number'].set(1)
    componentVars[cmpnt_nm]['Cells']['Locations'] = {}

    for ind in range(2):
        # Create variable to keep value of the coordinates
        componentVars[cmpnt_nm]['Cells']['Locations'][LOC_ARG_NAMES[ind]] = StringVar()
        componentVars[cmpnt_nm]['Cells']['Locations'][LOC_ARG_NAMES[ind]].set("")

    componentVars[cmpnt_nm]['Cells']['Locations']['fileinput'] = BooleanVar()
    componentVars[cmpnt_nm]['Cells']['Locations']['fileinput'].set(0)
    componentVars[cmpnt_nm]['Cells']['Locations']['filename'] = StringVar()
    componentVars[cmpnt_nm]['Cells']['Locations']['filename'].set("")

    # Regular parameters
    componentVars[cmpnt_nm]['Params'] = (
        app.populate_params_dict(SH_PARAMETER_VALUES))

    # Setup influenceModel parameter part of the dictionary componentVars[cmpnt_nm]['Params']
    componentVars[cmpnt_nm]['Params']['influenceModel'] = {}
    componentVars[cmpnt_nm]['Params']['influenceModel']['distribution'] = StringVar()
    componentVars[cmpnt_nm]['Params']['influenceModel']['distribution'].set('Fixed Value')
    componentVars[cmpnt_nm]['Params']['influenceModel']['value'] = IntVar()
    componentVars[cmpnt_nm]['Params']['influenceModel']['value'].set(0)

    # Dynamic keyword arguments
    if 'Dynamic' in cnctn_nm:  # dynamic kwargs are provided instead of connected component
        for ind, key in enumerate(SH_DYNAMIC_KWARGS):
            componentVars[cmpnt_nm][key] = dyn_data[ind]

    # Outputs of the component
    for output_key in SH_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Seal Horizon Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Cells setup frame
    app.cells_setup_frame = tk.Frame(tab)
    app.cells_setup_frame.grid(row=1, column=0, sticky='w', pady=0)
    add_cells_setup_frame_widgets(app, cmpnt_nm, tool_tip)

    # Influence model frame and its widgets
    infl_model_frame = tk.Frame(tab)
    infl_model_frame.grid(row=2, column=0, sticky='w',
                          padx=PARAMETER_FRAME_PADX, pady=0)
    infl_model_label = ttk.Label(
       infl_model_frame, text='Influence model:', width=PARAMETER_LABEL_WIDTH)
    infl_model_label.grid(row=0, column=0, pady=5, padx=5, sticky='w')
    infl_model_menu = tk.OptionMenu(
        infl_model_frame, componentVars[cmpnt_nm]['Params']['influenceModel']['value'],
        *INFLUENCE_MODEL_TYPES,
        command=lambda option: add_time_model_setup_frame_widgets(
            option, app, cmpnt_nm))
    tool_tip.bind(infl_model_menu,
                  'Set type of influence model.')
    infl_model_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    infl_model_menu.grid(row=0, column=1, padx=5)

    app.time_model_pars_frame = tk.Frame(tab)
    app.time_model_pars_frame.grid(row=3, column=0, sticky='w',
                                   padx=PARAMETER_FRAME_PADX)
    app.time_model_pars_frame.tool_tip = tool_tip
    add_time_model_setup_frame_widgets(0, app, cmpnt_nm)

    # Parameters frames
    par_frames = {}
    row_ind = 4
    for ind, par_name in enumerate(SH_PARAMETERS[1][1:]): # excluding influenceModel
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=row_ind, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        app.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': SH_PARAMETER_VALUES[par_name][6],
             'upper_bound': SH_PARAMETER_VALUES[par_name][7]},
            SH_PARAMETERS_SETUP[par_name][0],
            SH_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, tool_tip)
        row_ind = row_ind+1

    # Setup of relative model part
    # Add option menu
    rel_model_frame = tk.Frame(tab)
    rel_model_frame.grid(row=row_ind, column=0, sticky='w',
                         padx=PARAMETER_FRAME_PADX)

    rm_label = ttk.Label(
        rel_model_frame, text=SH_PARAMETERS_SETUP['relativeModel'][0],
        width=PARAMETER_LABEL_WIDTH)
    rm_menu = tk.OptionMenu(
        rel_model_frame, componentVars[cmpnt_nm]['relativeModel'],
        *RELATIVE_MODEL_TYPES,
        command=lambda option: add_rel_model_pars_frame_widgets(
            option, app, cmpnt_nm))
    rm_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    rm_label.grid(row=0, column=0, sticky='e', padx=5)
    rm_menu.grid(row=0, column=1, padx=5)
    tool_tip.bind(rm_menu, SH_PARAMETERS_SETUP['relativeModel'][1])

    row_ind = row_ind + 1

    # Add relative model parameters frame and its widgets
    app.rel_model_pars_frame = tk.Frame(tab)
    app.rel_model_pars_frame.grid(row=row_ind, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)
    app.rel_model_pars_frame.tool_tip = tool_tip
    add_rel_model_pars_frame_widgets('LET', app, cmpnt_nm)
    row_ind = row_ind + 1

    # Outputs
    cmpnt_outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    cmpnt_outputs_label.grid(row=row_ind, column=0, sticky='w', pady=(5, 10))

    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(row=row_ind+2, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    # Create and place outputs widgets
    output_nms_labels = []
    output_nms_checkboxes = []

    for ind1 in range(2):
        for ind2 in range(4):
            ind = ind1*4 + ind2
            obs_nm = SH_OBSERVATIONS[ind]
            # Create output checkbox
            output_nms_checkboxes.append(
                tk.Checkbutton(outputs_frame, variable=componentVars[cmpnt_nm][obs_nm]))
            # Place checkbox
            output_nms_checkboxes[-1].grid(
                row=ind2+1, column=2*ind1, pady=5, padx=CB_PADX, sticky='w')
            # Create output label
            output_nms_labels.append(
                ttk.Label(outputs_frame, text=SH_OBSERVATIONS_SETUP[obs_nm][0],
                          width=OUTPUT_LABEL_WIDTH1, anchor='w'))
            # Place label
            output_nms_labels[-1].grid(row=ind2+1, column=2*ind1+1, pady=5, sticky='w')
            # Bind checkbox to the tip
            tool_tip.bind(output_nms_checkboxes[-1],
                          SH_OBSERVATIONS_SETUP[obs_nm][1])


def add_rel_model_pars_frame_widgets(option, app, cmpnt_nm, change_distribution=False):
    for widget in app.rel_model_pars_frame.winfo_children():
        # Delete all widgets from parameters frame
        app.rel_model_pars_frame.tool_tip.unbind(widget)
        widget.destroy()

    par_dict = {'LET': LET_MODEL_PARAMETERS, 'BC': BC_MODEL_PARAMETERS}
    par_frames = {}

    for ind, par_name in enumerate(par_dict[option]):
        par_frames[par_name] = tk.Frame(app.rel_model_pars_frame)
        par_frames[par_name].grid(row=ind, column=0, sticky='w',
                                  padx=20)

        app.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': SH_PARAMETER_VALUES[par_name][6],
             'upper_bound': SH_PARAMETER_VALUES[par_name][7]},
            SH_PARAMETERS_SETUP[par_name][0],
            SH_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, app.rel_model_pars_frame.tool_tip)

        if change_distribution:
            app.change_distribution(par_frames[par_name])


def add_time_model_setup_frame_widgets(option, app, cmpnt_nm, change_distribution=False):
    for widget in app.time_model_pars_frame.winfo_children():
        # Delete all widgets from parameters frame
        app.time_model_pars_frame.tool_tip.unbind(widget)
        widget.destroy()

    if option == 0:
        # Hide frame containing widgets
        app.time_model_pars_frame.grid_remove()
    else:
        # Restore frame containing widgets to the previously placed location
        app.time_model_pars_frame.grid()

    par_frames = {}
    row_ind = 0
    # Setup frames of parameters common for both types of time model
    if option in [1, 2]:
        for ind, par_name in enumerate(SH_PARAMETERS[3]):
            par_frames[par_name] = tk.Frame(app.time_model_pars_frame)
            par_frames[par_name].grid(row=row_ind, column=0, sticky='w',
                                      padx=PARAMETER_FRAME_PADX)
            # Reset corresponding variables for distribution types
            componentVars[cmpnt_nm]['Params'][par_name]['distribution'].set(
                DISTRIBUTION_OPTIONS[0])
            # Add parameter frame
            app.setup_parameter_frame(
                par_frames[par_name], par_name,
                {'lower_bound': SH_PARAMETER_VALUES[par_name][6],
                 'upper_bound': SH_PARAMETER_VALUES[par_name][7]},
                SH_PARAMETERS_SETUP[par_name][0],
                SH_PARAMETERS_SETUP[par_name][1],
                PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
                cmpnt_nm, app.time_model_pars_frame.tool_tip)
            row_ind = row_ind  + 1

            # Update widgets if distribution is different
            if change_distribution:
                app.change_distribution(par_frames[par_name])

    # Setup frames of parameters typical for time model 2
    if option == 2:
        for ind, par_name in enumerate(SH_PARAMETERS[4]):
            par_frames[par_name] = tk.Frame(app.time_model_pars_frame)
            par_frames[par_name].grid(row=row_ind, column=0, sticky='w',
                                      padx=PARAMETER_FRAME_PADX)

            app.setup_parameter_frame(
                par_frames[par_name], par_name,
                {'lower_bound': SH_PARAMETER_VALUES[par_name][6],
                 'upper_bound': SH_PARAMETER_VALUES[par_name][7]},
                SH_PARAMETERS_SETUP[par_name][0],
                SH_PARAMETERS_SETUP[par_name][1],
                PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
                cmpnt_nm, app.time_model_pars_frame.tool_tip)

            # Update widgets if distribution is different
            if change_distribution:
                app.change_distribution(par_frames[par_name])

            row_ind = row_ind  + 1

        # Add remaining parameter
        clay_type_frame = tk.Frame(app.time_model_pars_frame)
        clay_type_frame.grid(row=row_ind, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

        ct_label = ttk.Label(clay_type_frame, text=SH_PARAMETERS_SETUP['clayType'][0],
                             width=PARAMETER_LABEL_WIDTH)
        ct_menu = tk.OptionMenu(
            clay_type_frame, componentVars[cmpnt_nm]['clayType'], *CLAY_TYPES)
        ct_menu.config(width=DISTRIBUTION_MENU_WIDTH)
        ct_label.grid(row=0, column=0, sticky='e', padx=5)
        ct_menu.grid(row=0, column=1, padx=5)
        app.time_model_pars_frame.tool_tip.bind(
            ct_menu, SH_PARAMETERS_SETUP['clayType'][1])

def add_cells_setup_frame_widgets(app, cmpnt_nm, tool_tip):
    LABEL_WIDTH = PARAMETER_LABEL_WIDTH

    number_frame = tk.Frame(app.cells_setup_frame)
    number_frame.grid(row=0, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)
    number_label = ttk.Label(
       number_frame, text='Number of cells:', width=LABEL_WIDTH)
    number_label.grid(row=0, column=0, pady=5, padx=5, sticky='w')
    number_field = tk.Entry(
        number_frame, width=2*DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
        textvariable=componentVars[cmpnt_nm]['Cells']['Number'])
    number_field.grid(row=0, column=1, pady=5, padx=10, sticky='w')
    tool_tip.bind(number_field, 'Set number of cells.')

    app.cells_setup_frame.cell_locs_frame = tk.Frame(app.cells_setup_frame)
    app.cells_setup_frame.cell_locs_frame.grid(
        row=1, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)
    add_cell_locs_frame_widgets(
        app, cmpnt_nm, app.cells_setup_frame.cell_locs_frame,
        componentVars[cmpnt_nm]['Cells']['Locations'], tool_tip)

    pars_frame = ttk.Frame(app.cells_setup_frame)
    pars_frame.grid(row=2, column=0, sticky='w', padx=0, pady=0)

    par_frames = {}
    # Create and place label and entry widgets
    for ind, par_name  in enumerate(SH_CELL_PARAMETERS):
        par_frames[par_name] = tk.Frame(pars_frame)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        app.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': SH_CELL_PARAMETERS_VALUES[par_name][6],
             'upper_bound': SH_CELL_PARAMETERS_VALUES[par_name][7]},
            SH_CELL_PARAMETERS_SETUP[par_name][0],
            SH_CELL_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            SPEC_DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Cells'][par_name],
            cmpnt_nm, tool_tip)

def load_additional_parameters(app, cmpnt_data, cmpnt_nm):
    """
    Load cell parameters data from dictionary into variables.
    """
    componentVars[cmpnt_nm]['Cells']['Number'].set(cmpnt_data['Cells']['Number'])

    if 'file' in cmpnt_data['Cells']['Locations']:
        componentVars[cmpnt_nm]['Cells']['Locations']['fileinput'].set(1)
        componentVars[cmpnt_nm]['Cells']['Locations']['filename'].set(
            cmpnt_data['Cells']['Locations']['file'])
    else:
        componentVars[cmpnt_nm]['Cells']['Locations']['fileinput'].set(0)
        for nm in LOC_ARG_NAMES:
            loc_data = ", ".join(
                str(item) for item in cmpnt_data['Cells']['Locations'][nm])
            componentVars[cmpnt_nm]['Cells']['Locations'][nm].set(loc_data)

    # Disable corresponding widgets
    disable_cell_locs_frame_widgets(app.cells_setup_frame.cell_locs_frame.coords_frame)

    # Populate cell parameters entry fields
    for par_nm in SH_CELL_PARAMETERS:
        # Get provided value
        value = cmpnt_data['Cells'][par_nm]
        if isinstance(value, str): # filename
            componentVars[cmpnt_nm]['Cells'][par_nm]['distribution'].set('File Input')
            componentVars[cmpnt_nm]['Cells'][par_nm]['filename'].set(value)
        elif isinstance(value, list): # array
            componentVars[cmpnt_nm]['Cells'][par_nm]['distribution'].set('List')
            # Create string representation of list
            str_values = ", ".join(str(item) for item in value)
            componentVars[cmpnt_nm]['Cells'][par_nm]['ordvalues'].set(str_values)
        else: # just one value
            componentVars[cmpnt_nm]['Cells'][par_nm]['distribution'].set('Fixed Value')
            componentVars[cmpnt_nm]['Cells'][par_nm]['value'].set(value)

        # Change distribution
        par_frame_name = '.'.join([cmpnt_nm, par_nm, 'frame'])
        app.change_distribution(app.nametowidget(app.getvar(par_frame_name)))

    # Setup widgets related to different values of relative model
    rel_model_value = cmpnt_data['Parameters']['relativeModel']
    add_rel_model_pars_frame_widgets(
        rel_model_value, app, cmpnt_nm, change_distribution=True)

    model_value = cmpnt_data['Parameters']['influenceModel']['value']
    add_time_model_setup_frame_widgets(
        model_value, app, cmpnt_nm, change_distribution=True)

def clean_data(d, cmpnt_nm):
    """
    Remove data about unnecessary parameters.
    """
    if d[cmpnt_nm]['Parameters']['relativeModel'] == 'LET':
        del d[cmpnt_nm]['Parameters']['lambda']
    elif d[cmpnt_nm]['Parameters']['relativeModel'] == 'BC':
        for par_nm in LET_MODEL_PARAMETERS:
            del d[cmpnt_nm]['Parameters'][par_nm]

    if d[cmpnt_nm]['Parameters']['influenceModel']['value'] in [0, 1]:
        for nm in SH_PARAMETERS[4]:
            del d[cmpnt_nm]['Parameters'][nm]
        del d[cmpnt_nm]['Parameters']['clayType']

    if d[cmpnt_nm]['Parameters']['influenceModel']['value'] == 0:
        for nm in SH_PARAMETERS[3]:
            del d[cmpnt_nm]['Parameters'][nm]
