"""
Module contains several methods needed for creating tab (page) in GUI
for FaultFlow component. Methods read and write dictionaries
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
from dictionarydata import DISTRIBUTION_OPTIONS
from dictionarydata import connections

from dictionarydata import LABEL_FONT
from dictionarydata import (PARAMETER_LABEL_WIDTH,
                            DISTRIBUTION_MENU_WIDTH,
                            DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            OUTPUT_LABEL_WIDTH1, PARAMETER_FRAME_PADX,
                            CB_PADX)
from cmpnts_tabs.parameter_entry import ParameterEntry
from cmpnts_tabs.commons import commons_read_tab_vars


FF_PARAMETERS = {
    1: ['xStart', 'yStart', 'length', 'nSegments', 'strike', 'dip'],
    2: ['SGR', 'stateVariable', 'aperture'],
    3: ['aquiferDepth', 'injectDepth', 'aquiferPressure',
        'fieldPressure', 'injectPressure', 'finalPressure', 'aquiferTemperature',
        'injectTemperature', 'injectX', 'injectY'],
    4: ['brineDensity', 'CO2Density', 'brineViscosity', 'CO2Viscosity',
        'salinity', 'CO2Solubility'],
    5: ['brineResSaturation', 'CO2ResSaturation'],
    7: ['permRatio', 'entryPressure'],
    8: ['maxHorizontal', 'minHorizontal', 'maxTrend']}

ALL_FF_PARAMETERS = [val for ind in [1, 2, 3, 4, 5, 7, 8] for val in FF_PARAMETERS[ind]]

RELATIVE_MODEL_TYPES = ['LET', 'BC']

BC_MODEL_PARAMETERS = ['lambda']

LET_MODEL_PARAMETERS = ['wetting1', 'wetting2', 'wetting3',
                        'nonwet1', 'nonwet2', 'nonwet3',
                        'capillary1', 'capillary2', 'capillary3', 'maxCapillary']

# Parameter labels names and parts of tooltip explanation
FF_PARAMETERS_SETUP = {
    'strike': ["Strike [deg]:",
               "fault strike"],
    'dip': ["Dip [deg]:",
            "fault dip"],
    'length': ["Trace length [m]:",
               "fault trace length"],
    'xStart': ["Start point x-coordinate [m]:",
               "x-coordinate"],
    'yStart': ["Start point y-coordinate [m]:",
               "y-coordinate"],
    'nSegments': ['Number of segments:',
                  "number of segments"],
    'SGR': ["Shale gauge ratio [%]:",
            "shale gauge ratio"],
    'stateVariable': ["Correction factor [-]",
                      "correction factor"],
    'aperture': ["Effective aperture [m]:",
                 "effective aperture"],
    # TODO we need to make aquiferDepth composite
    'aquiferDepth': ["Depth to aquifer above fault [m]:",
                     "depth to base of deepest aquifer along/above fault"],
    'injectDepth': ["Injection depth [m]:", "injection depth"],
    'aquiferPressure': ["Aquifer pressure [Pa]:", "aquifer pressure"],
    'fieldPressure': ["Initial pressure [Pa]:", "initial pressure at injection depth"],
    'injectPressure': ["Injection pressure [Pa]:",
                       "pressure at base during injection period"],
    'finalPressure': ["Final pressure [Pa]:",
                      "final pressure at injection depth"],
    'aquiferTemperature': ["Aquifer temperature [{}C]:".format(u'\u00B0'),
                           "temperature of brine in the aquifer"],
    'injectTemperature': ["Reservoir temperature [{}C]:".format(u'\u00B0'),
                          "temperature of brine at injection depth in reservoir"],
    'injectX': ["Injection well x-coordinate [m]:",
                "injection well x-coordinate"],
    'injectY': ["Injection well y-coordinate [m]:",
                "Injection well y-coordinate"],
    'salinity': ["Brine salinity [ppm]:",
                 "brine salinity"],
    'CO2Viscosity': ["CO{} viscosity [Pa{}s]:".format(u'\u2082', u'\u22C5'),
                     "CO{} viscosity".format(u'\u2082')],
    'CO2Density': ["CO{} density [kg/m{}]:".format(u'\u2082', u'\u00B3'),
                   "CO{} density".format(u'\u2082')],
    'brineDensity': ["Brine density [kg/m{}]:".format(u'\u00B3'),
                     "brine density"],
    'brineViscosity': ["Brine viscosity [Pa{}s]:".format(u'\u22C5'),
                       "brine viscosity"],
    'CO2Solubility': ["CO{} solubility [mol/kg]:".format(u'\u2082'),
                      "CO{} solubility".format(u'\u2082')],
    'brineResSaturation': ["Brine residual saturation [-]:", "brine residual saturation"],
    'CO2ResSaturation': ["CO{} residual saturation [-]:".format(u'\u2082'),
                         "CO{} residual saturation".format(u'\u2082')],
    'relativeModel': ["Relative model:", "relative model"],
    'permRatio': ["Permeability ratio [-]:",
                  "ratio of nonwetting permeability to the wetting permeability"],
    'entryPressure': ["Entry pressure [Pa]:", "entry pressure"],
    'lambda': ["Lambda [-]:", "lambda"],
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
                     'maximum capillary pressure'],
    'maxHorizontal': ["Maximum horizontal stress [Pa]:",
                      "secondary maximum horizontal principal stress"],
    'minHorizontal': ["Minimum horizontal stress [Pa]:",
                      "secondary minimum horizontal principal stress"],
    'maxTrend': ["Stress strike [deg]:",
                 "strike of secondary maximum horizontal stress"]}

# Set Fault Flow parameter names and value, min, max, second value, mean, std,
# low_bound, upper_bound
FF_PARAMETER_VALUES = {
    'strike': [30, 25, 50, 45, 30, 10, 0, 360],
    'dip': [70, 45, 90, 50, 45, 15, 0, 180],
    'length': [1000, 500, 1500, 750, 1000, 250, 0, 10000],
    'xStart': [500, 100, 600, 400, 300, 100, -5.0e+07, 5.0e+07],
    'yStart': [300, 100, 600, 400, 300, 100, -5.0e+07, 5.0e+07],
    'nSegments': [1, 4, 8, 5, 4, 2, 1, 100],
    'SGR': [0, 15, 30, 5, 10, 5, 0, 100],
    'stateVariable': [1, 0.6, 1, 0.9, 0.7, 0.1, 0, 1],
    'aperture': [2.5e-6, 1.0e-6, 3.0e-6, 1.5e-6, 2.0e-6, 1.0e-7, 0, 0.025],
    'aquiferDepth': [240, 200, 270, 150, 150, 25, 0, 2000],
    'injectDepth': [1880, 1800, 1900, 1700, 1500, 100, 800, 10000],
    'aquiferPressure': [2.455e+6, 2.4e+6, 2.5e+6, 2.5e+6, 2.4e+6, 1.0e+4, 0, 1.0e+7],
    'fieldPressure': [1.872e+7, 8.0e+6, 1.0e+7, 2.0e+7, 1.9e+7, 1.0e+5, 1.0e+5, 6.0e+7],
    'injectPressure': [2.89e+7, 1.0e+7, 2.0e+7, 2.0e+7, 3.0e+7, 1.0e+6, 7.0e+6, 6.0e+8],
    'finalPressure': [3.9e+7, 1.0e+7, 4.0e+7, 2.0e+7, 3.0e+7, 1.0e+6, 1.0e+5, 6.0e+7],
    'aquiferTemperature': [30, 10, 100, 50, 50, 15, 5, 180],
    'injectTemperature': [95, 70, 100, 80, 75, 15, 31, 180],
    'injectX': [0, 500, 1000, 750, 0, 250, -5.0e+07, 5.0e+07],
    'injectY': [0, 500, 1000, 750, 0, 250, -5.0e+07, 5.0e+07],
    'salinity': [15000, 1000, 10000, 9000, 10000, 1000, 0, 80000],
    'CO2Density': [430, 400, 650, 500, 450, 50, 93, 1050],
    'CO2Viscosity': [3.75e-5, 2.0e-5, 1.0e-4, 3.0e-5, 1.0e-4, 1.0e-6, 1.8e-05, 1.4e-04],
    'brineDensity': [988, 900, 1000, 950, 950, 25, 880, 1080],
    'brineViscosity': [4.36e-4, 2.0e-4, 8.0e-4, 5.0e-4,
                       8.0e-4, 1.0e-5, 1.5e-04, 1.6e-03],
    'CO2Solubility': [1.1, 0.5, 1, 1.5, 1, 0.2, 0, 2],
    'brineResSaturation': [0.15, 0.1, 0.3, 0.2, 0.2, 0.025, 0.01, 0.35],
    'CO2ResSaturation': [0.01, 0.1, 0.3, 0.2, 0.2, 0.025, 0.01, 0.35],
    # 'relativeModel': ['BC'],
    'permRatio': [0.6, 0.1, 1, 0.5, 0.7, 0.02, 0.1, 1.5],
    'entryPressure': [5000, 4500, 8000, 6000, 6000, 250, 100, 2.0e+6],
    'lambda': [2.5, 2, 4, 3, 2.5, 0.025, 0, 5],
    'wetting1': [1, 0.5, 3, 2, 2, 0.25, 0.5, 5],
    'wetting2': [10, 2, 15, 5, 10, 1, 1, 30],
    'wetting3': [1.25, 1, 2.5, 2, 1.5, 0.1, 0.2, 3],
    'nonwet1': [1.05, 0.5, 3, 2, 2, 0.25, 0.5, 5],
    'nonwet2': [10, 2, 15, 5, 10, 1, 1, 30],
    'nonwet3': [1.25, 1, 2.5, 2, 1.5, 0.1, 0, 3],
    'capillary1': [0.2, 0.1, 0.5, 0.3, 0.5, 0.05, 0.01, 5],
    'capillary2': [2.8, 0.5, 3.5, 2, 2.5, 0.1, 0.01, 30],
    'capillary3': [0.43, 0.1, 1, 0.5, 1, 0.01, 0.01, 3],
    'maxCapillary': [1.0e+7, 9.8e+6, 5.0e+7, 1.2e+7, 1.1e+7, 1.0e+5, 100, 2.0e+8],
    'maxHorizontal': [3.0e+7, 1.0e+7, 4.0e+7, 2.0e+7, 2.5e+7, 1.0e+6, 0, 5.0e+7],
    'minHorizontal': [2.0e+7, 1.0e+7, 4.0e+7, 1.5e+7, 2.5e+7, 1.0e+6, 0, 5.0e+7],
    'maxTrend': [55, 35, 90, 45, 55, 10, 0, 180]}

FF_DYNAMIC_KWARGS = ['pressure', 'CO2saturation']

FF_OBSERVATIONS = ['CO2_aquifer', 'mass_CO2_aquifer',
                   'CO2_aquifer_total', 'mass_CO2_aquifer_total',
                   'brine_aquifer', 'mass_brine_aquifer',
                   'brine_aquifer_total', 'mass_brine_aquifer_total']

FF_OBSERVATIONS_SETUP = {
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


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'FaultFlow', parameter_names=ALL_FF_PARAMETERS,
        dynamic_kwarg_names=FF_DYNAMIC_KWARGS,
        observation_names=FF_OBSERVATIONS)

    # # TODO DO NOT DELETE. Possibly needed to connect fault flow component
    # # to aquifer components
    # cmpnt_data['AquiferName'] = componentVars[cmpnt_nm]['aquiferName'].get()

    # Additional non-scalar parameters/keyword arguments
    cmpnt_data['Parameters']['relativeModel'] = componentVars[cmpnt_nm]['relativeModel'].get()
    for par_nm in LET_MODEL_PARAMETERS+BC_MODEL_PARAMETERS:
        cmpnt_data['Parameters'][par_nm] = {}

    # Get Segments parameters
    cmpnt_data['Segments'] = {}
    cmpnt_data['Segments']['Number'] = (
        componentVars[cmpnt_nm]['Params']['nSegments']['value'].get())

    # TODO Need to think whether we need something that reads location
    # or we're going to extract this from length, x-start, y-start variables
    # Read information about locations associated with component
    # read_locations_data(d, cmpnt_nm)

    return cmpnt_data


def add_widgets(app, tab, cmpnt_nm, cmpnt_type, tool_tip,
                cnctn_nm, dyn_data, aquifer_nm, *args):
    """ Add widgets to the component tab."""
    # # TODO DO NOT DELETE. Possibly needed to connect fault flow component
    # # to aquifer components
    # aquifers = ['aquifer{}'.format(ind) for ind in range(
    #     1, componentVars['strata']['Params']['numberOfShaleLayers'].get())]

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('FaultFlow')
    connectionsDictionary.append(cnctn_nm)

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set(cnctn_nm)

    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    # # TODO DO NOT DELETE. Possibly needed to connect fault flow component
    # # to aquifer components
    # componentVars[cmpnt_nm]['aquiferName'] = StringVar()
    # componentVars[cmpnt_nm]['aquiferName'].set(aquifer_nm)

    # Relative model variable
    componentVars[cmpnt_nm]['relativeModel'] = StringVar()
    componentVars[cmpnt_nm]['relativeModel'].set('LET')

    if 'Dynamic' in cnctn_nm:  # dynamic kwargs are provided instead of connected component
        for ind, key in enumerate(FF_DYNAMIC_KWARGS):
            componentVars[cmpnt_nm][key] = dyn_data[ind]

    componentVars[cmpnt_nm]['Params'] = app.populate_params_dict(
        FF_PARAMETER_VALUES)
    # Update nSegments variable
    componentVars[cmpnt_nm]['Params']['nSegments']['value'] = IntVar()
    componentVars[cmpnt_nm]['Params']['nSegments']['value'].set(1)

    for output_key in FF_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Fault Flow Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Fault setup frame
    fault_setup_frame = tk.Frame(tab)
    fault_setup_frame.grid(row=1, column=0, sticky='w', pady=0)
    add_fault_setup_frame_widgets(app, cmpnt_nm, fault_setup_frame, tool_tip)

    # Parameters frames
    par_frames = {}
    row_ind = 2
    for group_id in range(3, 9):
        if group_id == 6: # corresponds to "relative model" parameter
            # Add option menu
            rel_model_frame = tk.Frame(tab)
            rel_model_frame.grid(row=row_ind, column=0, sticky='w',
                                 padx=PARAMETER_FRAME_PADX)

            rm_label = ttk.Label(
                rel_model_frame, text=FF_PARAMETERS_SETUP['relativeModel'][0],
                width=PARAMETER_LABEL_WIDTH)
            rm_menu = tk.OptionMenu(
                rel_model_frame, componentVars[cmpnt_nm]['relativeModel'],
                *RELATIVE_MODEL_TYPES,
                command=lambda option: add_rel_model_pars_frame_widgets(
                    option, app, cmpnt_nm))
            rm_menu.config(width=DISTRIBUTION_MENU_WIDTH)
            rm_label.grid(row=0, column=0, sticky='e', padx=5)
            rm_menu.grid(row=0, column=1, padx=5)
            tool_tip.bind(rm_menu, FF_PARAMETERS_SETUP['relativeModel'][1])

            # Add relative model parameters frame and its widgets
            app.rel_model_pars_frame = tk.Frame(tab)
            app.rel_model_pars_frame.grid(row=row_ind+1, column=0, sticky='w',
                                            padx=PARAMETER_FRAME_PADX)
            app.rel_model_pars_frame.tool_tip = tool_tip
            add_rel_model_pars_frame_widgets('LET', app, cmpnt_nm)
            row_ind = row_ind+2
        else:
            for ind, par_name in enumerate(FF_PARAMETERS[group_id]):
                par_frames[par_name] = tk.Frame(tab)
                par_frames[par_name].grid(row=row_ind, column=0, sticky='w',
                                          padx=PARAMETER_FRAME_PADX)

                app.setup_parameter_frame(
                    par_frames[par_name], par_name,
                    {'lower_bound': FF_PARAMETER_VALUES[par_name][6],
                     'upper_bound': FF_PARAMETER_VALUES[par_name][7]},
                    FF_PARAMETERS_SETUP[par_name][0],
                    FF_PARAMETERS_SETUP[par_name][1],
                    PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                    DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                    DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
                    cmpnt_nm, tool_tip)
                row_ind = row_ind+1

    # Connection frame
    con_frame = tk.Frame(tab)
    con_frame.grid(row=34, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    connection_label = ttk.Label(
        con_frame, text="Connection:", width=PARAMETER_LABEL_WIDTH)
    connection_menu = tk.OptionMenu(
        con_frame, componentVars[cmpnt_nm]['connection'], *connections)
    connection_label.grid(row=0, column=0, sticky='e', padx=5)
    connection_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    connection_menu.grid(row=0, column=1, padx=5)
    tool_tip.bind(connection_menu,
                  'Set connection for this component.')

    # # TODO DO NOT DELETE. Possibly needed to connect fault flow component
    # # to aquifer components
    # aq_name_label = ttk.Label(
    #     con_frame, width=PARAMETER_LABEL_WIDTH, text="Aquifer name:")
    # aq_name_menu = tk.OptionMenu(
    #     con_frame, componentVars[cmpnt_nm]['aquiferName'], *aquifers)
    # aq_name_label.grid(row=1, column=0, sticky='e', padx=5)
    # aq_name_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    # aq_name_menu.grid(row=1, column=1, padx=5)
    # tool_tip.bind(aq_name_menu,
    #               'Set the aquifer name for this component.')

    # Outputs
    cmpnt_outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    cmpnt_outputs_label.grid(row=36, column=0, sticky='w', pady=(5, 10))

    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(row=37, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    # Create and place outputs widgets
    output_nms_labels = []
    output_nms_checkboxes = []

    for ind1 in range(2):
        for ind2 in range(4):
            ind = ind1*4 + ind2
            obs_nm = FF_OBSERVATIONS[ind]
            # Create output checkbox
            output_nms_checkboxes.append(
                tk.Checkbutton(outputs_frame, variable=componentVars[cmpnt_nm][obs_nm]))
            # Place checkbox
            output_nms_checkboxes[-1].grid(
                row=ind2+1, column=2*ind1, pady=5, padx=CB_PADX, sticky='w')
            # Create output label
            output_nms_labels.append(
                ttk.Label(outputs_frame, text=FF_OBSERVATIONS_SETUP[obs_nm][0],
                          width=OUTPUT_LABEL_WIDTH1, anchor='w'))
            # Place label
            output_nms_labels[-1].grid(row=ind2+1, column=2*ind1+1, pady=5, sticky='w')
            # Bind checkbox to the tip
            tool_tip.bind(output_nms_checkboxes[-1],
                          FF_OBSERVATIONS_SETUP[obs_nm][1])


def add_rel_model_pars_frame_widgets(option, app, cmpnt_nm):
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
            {'lower_bound': FF_PARAMETER_VALUES[par_name][6],
             'upper_bound': FF_PARAMETER_VALUES[par_name][7]},
            FF_PARAMETERS_SETUP[par_name][0],
            FF_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, app.rel_model_pars_frame.tool_tip)


def add_fault_setup_frame_widgets(app, cmpnt_nm, frame, tool_tip):
    LABEL_WIDTH = PARAMETER_LABEL_WIDTH
    # setup_label = ttk.Label(frame, text="Fault setup:")
    # setup_label.grid(row=0, column=0, columnspan=2, sticky='w', padx=0)

    pars_frame = ttk.Frame(frame)
    pars_frame.grid(row=1, column=0, sticky='w', padx=0)

    par_labels = []
    par_fields = []
    par_frames = {}
    # Create and place label and entry widgets
    for ind, par_name  in enumerate(FF_PARAMETERS[1]):
        par_frames[par_name] = tk.Frame(pars_frame)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)
        par_labels.append(ttk.Label(
            par_frames[par_name], text=FF_PARAMETERS_SETUP[par_name][0], width=LABEL_WIDTH))

        if par_name == 'nSegments':
            standard_tooltip_text = ''.join([
                'Set value of {}.\nPossible values are integers between ',
                '{} and {}.']).format(FF_PARAMETERS_SETUP[par_name][1],
                                      FF_PARAMETER_VALUES[par_name][6],
                                      FF_PARAMETER_VALUES[par_name][7])
            par_fields.append(ParameterEntry(
                par_frames[par_name], par_name, width=2*DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                text_variable=componentVars[cmpnt_nm]['Params'][par_name]['value'],
                tool_tip_link=tool_tip,
                standard_tooltip_text=standard_tooltip_text,
                discrete_bounds=list(range(1, 101)),
                discrete_bounds_msg=['', '', '', ''.join([
                    'The parameter value must belong to the list of integers ',
                    'between {} and {}.']).format(FF_PARAMETER_VALUES[par_name][6],
                                                  FF_PARAMETER_VALUES[par_name][7])]))
        else:
            standard_tooltip_text = ''.join([
                'Set value of {}.\nPossible values are between ',
                '{} and {}.']).format(FF_PARAMETERS_SETUP[par_name][1],
                                     FF_PARAMETER_VALUES[par_name][6],
                                     FF_PARAMETER_VALUES[par_name][7])
            par_fields.append(ParameterEntry(
                par_frames[par_name], par_name, width=2*DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                text_variable=componentVars[cmpnt_nm]['Params'][par_name]['value'],
                tool_tip_link=tool_tip,
                standard_tooltip_text=standard_tooltip_text,
                lower_bound=FF_PARAMETER_VALUES[par_name][6],
                upper_bound=FF_PARAMETER_VALUES[par_name][7]))

        par_labels[-1].grid(row=0, column=0, pady=5, padx=5, sticky='w')
        par_fields[-1].grid(row=0, column=1, pady=5, padx=10, sticky='w')
        tool_tip.bind(
            par_fields[-1],
            'Set value of {}.\nPossible values are between {} and {}.'.format(
                FF_PARAMETERS_SETUP[par_name][1],
                FF_PARAMETER_VALUES[par_name][6],
                FF_PARAMETER_VALUES[par_name][7]))

    row_ind = ind+2

    for ind, par_name in enumerate(FF_PARAMETERS[2]):
        par_frames[par_name] = tk.Frame(pars_frame)
        par_frames[par_name].grid(row=row_ind, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX, pady=0)

        app.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': FF_PARAMETER_VALUES[par_name][6],
             'upper_bound': FF_PARAMETER_VALUES[par_name][7]},
            FF_PARAMETERS_SETUP[par_name][0],
            FF_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, tool_tip)
        row_ind = row_ind+1


def clean_data(d, cmpnt_nm):
    """
    Remove data about unnecessary parameters.
    """
    if d[cmpnt_nm]['Parameters']['relativeModel'] == 'LET':
        del d[cmpnt_nm]['Parameters']['lambda']
    elif d[cmpnt_nm]['Parameters']['relativeModel'] == 'BC':
        for par_nm in LET_MODEL_PARAMETERS:
            del d[cmpnt_nm]['Parameters'][par_nm]
