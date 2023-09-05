"""
Module contains several methods needed for creating tab (page) in GUI
for MultisegmentedWellbore component. Methods read and write dictionaries
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

from cmpnts_tabs.locations import read_locations_data, add_wellbore_frame_widgets
from cmpnts_tabs.commons import commons_read_tab_vars


MSW_PARAMETERS = [
    'logWellPerm', 'logAquaPerm', 'brineDensity',
    'CO2Density', 'brineViscosity', 'CO2Viscosity',
    'aquBrineResSaturation', 'compressibility', 'wellRadius']

# Parameter labels names and parts of tooltip explanation
MSW_PARAMETERS_SETUP = {
    'logWellPerm': ["Well permeability [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B2'),
                    'well permeability'],
    'logAquaPerm': ["Aquifer permeability [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B2'),
                    'aquifer permeability'],
    'brineDensity': ["Brine density [kg/m{}]:".format(u'\u00B3'),
                     'brine density'],
    'CO2Density': ["CO{} density [kg/m{}]:".format(u'\u2082', u'\u00B3'),
                   "CO{} density".format(u'\u2082')],
    'brineViscosity': ["Brine viscosity [Pa{}s]:".format(u'\u22C5'),
                       'brine viscosity'],
    'CO2Viscosity': ["CO{} viscosity [Pa{}s]:".format(u'\u2082', u'\u22C5'),
                     "CO{} viscosity".format(u'\u2082')],
    'aquBrineResSaturation': ["Brine saturation [-]:",
                              'brine residual saturation in aquifer'],
    'compressibility': ["Compressibility [Pa{}{}]:".format(u'\u207B', u'\u00B9'),
                        'compressibility'],
    'wellRadius': ["Well radius [m]:",
                   'well radius']}

# Set Multisegmented Wellbore parameter names and value, min, max, second value,
# mean, std, low_bound, upper_bound
MSW_PARAMETER_VALUES = {
    'logWellPerm': [-13, -17, -9, -11, -13, 2, -17, -9],
    'logAquaPerm': [-11, -14, -9, -12, -11, 1, -14, -9],
    'wellRadius': [0.01, 0.01, 0.5, 0.05, 0.1, 0.05, 0.01, 0.5],
    'brineDensity': [1000, 900, 1500, 1200, 1200, 150, 900, 1500],
    'CO2Density': [479, 100, 1500, 750, 800, 250, 100, 1500],
    'brineViscosity': [2.535e-3, 1.0e-4, 5.0e-3, 3.0e-4,
                       2.5e-3, 1.0e-3, 1.0e-4, 5.0e-3],
    'CO2Viscosity': [3.95e-5, 1.0e-5, 1.0e-4, 6.0e-5,
                     5.5e-5, 1.5e-5, 1.0e-6, 1.0e-4],
    'aquBrineResSaturation': [0.1, 0, 0.7, 0.25, 0.5, 0.2, 0, 0.7],
    'compressibility': [5.1e-11, 5.0e-11, 1.0e-9, 1.0e-10,
                        5.25e-10, 3.5e-10, 5.0e-11, 1.0e-9]}

MSW_DYNAMIC_KWARGS = ['pressure', 'CO2saturation']


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'MultisegmentedWellbore', parameter_names=MSW_PARAMETERS,
        dynamic_kwarg_names=MSW_DYNAMIC_KWARGS, add_number=True)

    # Read information about locations associated with component
    read_locations_data(cmpnt_data, cmpnt_nm)

    modelOutputs = []
    for key in componentVars[cmpnt_nm]['outputs']:
        if componentVars[cmpnt_nm]['outputs'][key].get():
            modelOutputs.append(key)

    cmpnt_data['Outputs'] = modelOutputs
    return cmpnt_data

def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, toolTip,
                cnctn_nm, dyn_data, *args):
    """ Add widgets to the component tab."""

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('MultisegmentedWellbore')
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
        for ind, key in enumerate(MSW_DYNAMIC_KWARGS):
            componentVars[cmpnt_nm][key] = dyn_data[ind]

    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(
        MSW_PARAMETER_VALUES)

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Multisegmented Wellbore Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(MSW_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': MSW_PARAMETER_VALUES[par_name][6],
             'upper_bound': MSW_PARAMETER_VALUES[par_name][7]},
            MSW_PARAMETERS_SETUP[par_name][0],
            MSW_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, toolTip)

    # Number of wellbores defined by the same parameters
    number = tk.Frame(tab)
    number.grid(row=10, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    # Connection
    connection_label = ttk.Label(
        number, text="Connection:", width=PARAMETER_LABEL_WIDTH)
    connection_menu = tk.OptionMenu(
        number, componentVars[cmpnt_nm]['connection'], *connections)
    connection_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    connection_label.grid(row=0, column=0, sticky='w', padx=5)
    connection_menu.grid(row=0, column=1, padx=5)
    toolTip.bind(connection_menu,
                 'Set connection for this component.')

    number_label = ttk.Label(number, width=PARAMETER_LABEL_WIDTH,
                             text='Number of wellbores:')
    number_spinbox = tk.Spinbox(
        number, from_=1, to=1000,
        textvariable=componentVars[cmpnt_nm]['number'])
    # number_spinbox.config(width=DISTRIBUTION_MENU_WIDTH)
    number_label.grid(row=1, column=0, sticky='w', padx=5)
    number_spinbox.grid(row=1, column=1, padx=5, pady=2, sticky='ew')
    toolTip.bind(number_spinbox,
                 ''.join(['Set the total number of wells for this wellbore ',
                          'component including random locations.']))

    # Wellbore locations frame
    well_locs_frame = tk.Frame(tab)
    well_locs_frame.grid(row=11, column=0, sticky='w',
                         padx=PARAMETER_FRAME_PADX, pady=(5, 10))
    add_wellbore_frame_widgets(controller, cmpnt_nm, well_locs_frame, toolTip)

    # Outputs
    outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    outputs_label.grid(row=12, column=0, sticky='w', pady=(5, 10))

    componentVars[cmpnt_nm]['outputs'] = {}
    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(row=13, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    row = 0
    col = 0
    for ind in range(1, componentVars['strata']['Params']['numberOfShaleLayers'].get()):
        tr_key = 'aquifer{}'.format(ind)

        componentVars[cmpnt_nm]['outputs']['CO2_{}'.format(tr_key)] = BooleanVar()
        componentVars[cmpnt_nm]['outputs']['brine_{}'.format(tr_key)] = BooleanVar()
        componentVars[cmpnt_nm]['outputs']['mass_CO2_{}'.format(tr_key)] = BooleanVar()
        componentVars[cmpnt_nm]['outputs']['CO2_{}'.format(tr_key)].set(0)
        componentVars[cmpnt_nm]['outputs']['brine_{}'.format(tr_key)].set(0)
        componentVars[cmpnt_nm]['outputs']['mass_CO2_{}'.format(tr_key)].set(0)

        label = ttk.Label(outputs_frame,
                          width=OUTPUT_LABEL_WIDTH1, anchor='w',
                          text='CO{} aquifer {} [kg/s]'.format(u'\u2082', ind))
        checkbox = tk.Checkbutton(
            outputs_frame, variable=componentVars[
                cmpnt_nm]['outputs']['CO2_{}'.format(tr_key)])
        brine_label = ttk.Label(outputs_frame,
                                width=OUTPUT_LABEL_WIDTH1, anchor='w',
                                text='Brine aquifer {} [kg/s]'.format(ind))
        brine_checkbox = tk.Checkbutton(
            outputs_frame, variable=componentVars[
                cmpnt_nm]['outputs']['brine_{}'.format(tr_key)])
        mass_CO2_label = ttk.Label(
            outputs_frame, width=OUTPUT_LABEL_WIDTH1, anchor='w',
            text='Mass CO{} aquifer {} [kg]'.format(u'\u2082', ind))
        mass_CO2_checkbox = tk.Checkbutton(
            outputs_frame, variable=componentVars[
                cmpnt_nm]['outputs']['mass_CO2_{}'.format(tr_key)])

        label.grid(row=row, column=col+1, pady=5, sticky='w')
        checkbox.grid(row=row, column=col, pady=5, padx=CB_PADX, sticky='w')
        toolTip.bind(checkbox,
                     'Enable CO{} leakage rates to aquifer {} as output.'.format(
                         u'\u2082', ind))

        brine_label.grid(row=row, column=col+3, pady=5, sticky='w')
        brine_checkbox.grid(row=row, column=col+2, pady=5, padx=CB_PADX, sticky='w')
        toolTip.bind(brine_checkbox,
                     'Enable brine leakage rates to aquifer {} as output.'.format(ind))

        mass_CO2_label.grid(row=row, column=col+5, pady=5, sticky='w')
        mass_CO2_checkbox.grid(row=row, column=col+4, pady=5, padx=CB_PADX, sticky='w')
        toolTip.bind(mass_CO2_checkbox,
                     'Enable CO{} mass leaked to aquifer {} as output.'.format(
                         u'\u2082', ind))

        row = row + 1

    componentVars[cmpnt_nm]['outputs']['CO2_atm'] = BooleanVar()
    componentVars[cmpnt_nm]['outputs']['brine_atm'] = BooleanVar()
    componentVars[cmpnt_nm]['outputs']['CO2_atm'].set(0)
    componentVars[cmpnt_nm]['outputs']['brine_atm'].set(0)

    CO2_atm_label = ttk.Label(outputs_frame, width=OUTPUT_LABEL_WIDTH1, anchor='w',
                              text='CO{} atm [kg/s]'.format(u'\u2082'))
    CO2_atm_checkbox = tk.Checkbutton(
        outputs_frame, variable=componentVars[cmpnt_nm]['outputs']['CO2_atm'])
    brine_atm_label = ttk.Label(
        outputs_frame, width=OUTPUT_LABEL_WIDTH1, anchor='w', text='Brine atm [kg/s]')
    brine_atm_checkbox = tk.Checkbutton(
        outputs_frame, variable=componentVars[cmpnt_nm]['outputs']['brine_atm'])

    CO2_atm_label.grid(row=row, column=1, pady=5, sticky='w')
    CO2_atm_checkbox.grid(row=row, column=0, pady=5, padx=CB_PADX, sticky='w')
    toolTip.bind(CO2_atm_checkbox,
                 'Enable CO{} leakage rates to atmosphere as output.'.format(u'\u2082'))
    brine_atm_label.grid(row=row, column=3, pady=5, sticky='w')
    brine_atm_checkbox.grid(row=row, column=2, pady=5, padx=CB_PADX, sticky='w')
    toolTip.bind(brine_atm_checkbox,
                 'Enable brine leakage rates to atmosphere as output.')
