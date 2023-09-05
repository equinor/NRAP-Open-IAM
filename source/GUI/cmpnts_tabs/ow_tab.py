"""
Module contains several methods needed for creating tab (page) in GUI
for OpenWellbore component. Methods read and write dictionaries
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
from dictionarydata import DISTRIBUTION_OPTIONS, connections

from dictionarydata import LABEL_FONT
from dictionarydata import (PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
                            DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            OUTPUT_LABEL_WIDTH1, PARAMETER_FRAME_PADX, CB_PADX)

from cmpnts_tabs.locations import read_locations_data, add_wellbore_frame_widgets
from cmpnts_tabs.commons import commons_read_tab_vars


OW_PARAMETERS = ['logReservoirTransmissivity', 'logAquiferTransmissivity',
                 'brineSalinity', 'brineDensity', 'wellRadius', 'critPressure']

OW_PARAMETERS_SETUP = {
    'logReservoirTransmissivity': ["Reservoir transmissivity [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B3'), 'reservoir transmissivity'],
    'logAquiferTransmissivity': ["Aquifer transmissivity [log{} m{}]:".format(
        u'\u2081'u'\u2080', u'\u00B3'), 'aquifer transmissivity'],
    'brineSalinity': ["Brine salinity [-]:", 'brine salinity'],
    'brineDensity': ["Brine density [kg/m{}]".format(u'\u00B3'), 'brine density'],
    'wellRadius': ["Well radius [m]:", 'well radius'],
    'critPressure': ["Critical pressure [Pa]:", 'critical pressure']}

# Set Open Wellbore parameters names and value, min, max, second value, mean, std, bounds
OW_PARAMETER_VALUES = {
    'logReservoirTransmissivity': [-10, -11, -9, -10.5, -10, 0.2, -11.27, -8.4],
    'logAquiferTransmissivity': [-10, -11, -9, -10.5, -10, 0.2, -11.27, -8.4],
    'brineSalinity': [0.1, 0, 0.2, 0.15, 0.1, 0.05, 0, 0.2],
    'brineDensity': [1012, 1000, 1100, 1050, 1025, 50, 900, 1200],
    'wellRadius': [0.05, 0.025, 0.25, 0.1, 0.1, 0.05, 0.025, 0.25],
    'critPressure': [5.0e+6, 1.0e+6, 2.0e+7, 1.0e+7, 5.0e+6, 1.0e+4, 1.0e+5, 9.0e+7]}

OW_DYNAMIC_KWARGS = ['pressure', 'CO2saturation']

# Observations
OW_OBSERVATIONS = ['CO2_aquifer', 'brine_aquifer', 'CO2_atm', 'brine_atm']

OW_OBSERVATIONS_SETUP = {
    'CO2_aquifer': [
        "CO{} aquifer [kg/s]".format(u'\u2082'),
        'Enable CO{} leakage rates to aquifer as output.'.format(u'\u2082')],
    'CO2_atm': [
        "CO{} atm [kg/s]".format(u'\u2082'),
        'Enable CO{} leakage rates to atmosphere as output.'.format(u'\u2082')],
    'brine_aquifer': ["Brine aquifer [kg/s]",
                      'Enable brine leakage rates to aquifer as output.'],
    'brine_atm': ["Brine atm [kg/s]",
                  'Enable brine leakage rates to atmosphere as output.']}


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    cmpnt_data = commons_read_tab_vars(
        cmpnt_nm, 'OpenWellbore', parameter_names=OW_PARAMETERS,
        dynamic_kwarg_names=OW_DYNAMIC_KWARGS,
        observation_names=OW_OBSERVATIONS, add_number=True, add_leak_to=True)

    # Read information about locations associated with component
    read_locations_data(cmpnt_data, cmpnt_nm)

    return cmpnt_data

def add_widgets(controller, tab, cmpnt_nm, cmpnt_type, tool_tip,
                cnctn_nm, dyn_data, aquiferName, controls, *args):
    """ Add widgets to the component tab. Note that aquiferName will be 'none'
    and is unused."""

    aquifers = ['aquifer{}'.format(ind) for ind in range(
        1, componentVars['strata']['Params']['numberOfShaleLayers'].get())]

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('OpenWellbore')
    connectionsDictionary.append(cnctn_nm)

    componentVars[cmpnt_nm] = {}
    componentVars[cmpnt_nm]['connection'] = StringVar()
    componentVars[cmpnt_nm]['connection'].set(cnctn_nm)
    componentVars[cmpnt_nm]['LeakTo'] = StringVar()
    componentVars[cmpnt_nm]['LeakTo'].set(aquifers[0])

    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    componentVars[cmpnt_nm]['number'] = IntVar()
    componentVars[cmpnt_nm]['number'].set(1)

    componentVars[cmpnt_nm]['useRandomLocDomain'] = BooleanVar()
    componentVars[cmpnt_nm]['useRandomLocDomain'].set(0)

    # Populate dictionary
    componentVars[cmpnt_nm]['Params'] = controller.populate_params_dict(
        OW_PARAMETER_VALUES)

    if 'Dynamic' in cnctn_nm:  # dynamic kwargs are provided instead of connected component
        for ind, key in enumerate(OW_DYNAMIC_KWARGS):
            componentVars[cmpnt_nm][key] = dyn_data[ind]

    # Set outputs of Open Wellbore
    for output_key in OW_OBSERVATIONS:
        componentVars[cmpnt_nm][output_key] = BooleanVar()
        componentVars[cmpnt_nm][output_key].set(0)

    # Tab title label
    comp_type_label = ttk.Label(
        tab, text="Open Wellbore Component", font=LABEL_FONT)
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    # Parameters frames
    par_frames = {}
    for ind, par_name in enumerate(OW_PARAMETERS):
        par_frames[par_name] = tk.Frame(tab)
        par_frames[par_name].grid(row=ind+1, column=0, sticky='w',
                                  padx=PARAMETER_FRAME_PADX)

        controller.setup_parameter_frame(
            par_frames[par_name], par_name,
            {'lower_bound': OW_PARAMETER_VALUES[par_name][6],
             'upper_bound': OW_PARAMETER_VALUES[par_name][7]},
            OW_PARAMETERS_SETUP[par_name][0],
            OW_PARAMETERS_SETUP[par_name][1],
            PARAMETER_LABEL_WIDTH, DISTRIBUTION_MENU_WIDTH,
            DISTRIBUTION_ARG_LABEL_WIDTH, DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            DISTRIBUTION_OPTIONS, componentVars[cmpnt_nm]['Params'][par_name],
            cmpnt_nm, tool_tip)

    # Critical pressure input in the Controls section
    componentVars[cmpnt_nm]['Controls'] = {}

    componentVars[cmpnt_nm]['Controls']['critPressureApproach'] = BooleanVar()
    componentVars[cmpnt_nm]['Controls']['enforceCritPressure'] = BooleanVar()

    if len(controls) > 0:
        if 'critPressureApproach' in controls:
            componentVars[cmpnt_nm]['Controls']['critPressureApproach'].set(
                controls['critPressureApproach'])
        else:
            componentVars[cmpnt_nm]['Controls']['critPressureApproach'].set(0)

        if 'enforceCritPressure' in controls:
            componentVars[cmpnt_nm]['Controls']['enforceCritPressure'].set(
                controls['enforceCritPressure'])
        else:
            componentVars[cmpnt_nm]['Controls']['enforceCritPressure'].set(0)
    else:
        componentVars[cmpnt_nm]['Controls']['critPressureApproach'].set(0)
        componentVars[cmpnt_nm]['Controls']['enforceCritPressure'].set(0)

    use_critP_frame = tk.Frame(tab)
    use_critP_frame.grid(row=7, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    use_critP_label = ttk.Label(
        use_critP_frame,
        text="Use critical pressure approach:",
        width=PARAMETER_LABEL_WIDTH)
    use_critP_label.grid(row=0, column=0, sticky='w', padx=5)

    specify_critP_frame = tk.Frame(tab)
    specify_critP_frame.grid(row=8, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    # critP_frame = tk.Frame(tab)
    # critP_frame.grid(row=8, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    use_critP_checkbox = tk.Checkbutton(
        use_critP_frame,
        variable=componentVars[cmpnt_nm]['Controls']['critPressureApproach'],
        command=lambda: use_crit_pressure_approach(
            componentVars[cmpnt_nm]['Controls'], specify_critP_frame,
            par_frames['brineDensity'], par_frames['critPressure']))
    use_critP_checkbox.grid(row=0, column=1, sticky='w', padx=5)
    tool_tip.bind(use_critP_checkbox,
                  ''.join(['Check to enable the use of a critical pressure in\n',
                           'leakage calculations.']))

    specify_critP_label = ttk.Label(
        specify_critP_frame,
        text="Use the critical pressure parameter (above), otherwise\n calculate it using the default approach:",
        width = 1.5 * PARAMETER_LABEL_WIDTH)
    specify_critP_label.grid(row=0, column=0, sticky='w', padx=40)

    specify_critP_checkbox = tk.Checkbutton(
        specify_critP_frame,
        variable=componentVars[cmpnt_nm]['Controls']['enforceCritPressure'],
        command=lambda: enforce_crit_pressure_approach(
            componentVars[cmpnt_nm]['Controls'],
            par_frames['brineDensity'], par_frames['critPressure']))
    specify_critP_checkbox.grid(row=0, column=1, sticky='w', padx=40)
    tool_tip.bind(specify_critP_checkbox,
                  ''.join(['Check to enable the use of a specific critical pressure in\n',
                           'leakage calculations (the critical pressure parameter\n',
                           'shown above). Otherwise, critical pressure is calculated\n',
                           'in the default manner (see the User Guide).']))

    use_crit_pressure_approach(
        componentVars[cmpnt_nm]['Controls'], specify_critP_frame,
        par_frames['brineDensity'], par_frames['critPressure'])

    enforce_crit_pressure_approach(
        componentVars[cmpnt_nm]['Controls'],
        par_frames['brineDensity'], par_frames['critPressure'])

    # Number of wellbores defined by the same parameters
    number_frame = tk.Frame(tab)
    number_frame.grid(row=9, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    # Leaked to aquifer name
    leak_to_label = ttk.Label(number_frame, text='Leak to:', width=PARAMETER_LABEL_WIDTH)
    leak_to_menu = tk.OptionMenu(number_frame, componentVars[cmpnt_nm]['LeakTo'],
                                 *['Atmosphere', *aquifers])
    leak_to_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    tool_tip.bind(leak_to_menu,
                  'Select an aquifer that this wellbore will leak to.')
    leak_to_label.grid(row=0, column=0, padx=5, sticky='w')
    leak_to_menu.grid(row=0, column=1, padx=5)

    # Connection
    connection_label = ttk.Label(
        number_frame, text="Connection:", width=PARAMETER_LABEL_WIDTH)
    connection_menu = tk.OptionMenu(
        number_frame, componentVars[cmpnt_nm]['connection'], *connections)
    connection_label.grid(row=1, column=0, sticky='w', padx=5)
    connection_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    connection_menu.grid(row=1, column=1, padx=5)
    tool_tip.bind(connection_menu, 'Set connection for this component.')

    number_label = ttk.Label(
        number_frame, text='Number of wellbores:', width=PARAMETER_LABEL_WIDTH)
    number_spinbox = tk.Spinbox(
        number_frame, from_=1, to=1000, textvariable=componentVars[cmpnt_nm]['number'])
    number_label.grid(row=2, column=0, sticky='w', padx=5)
    number_spinbox.grid(row=2, column=1, padx=5, pady=2, sticky='ew')
    tool_tip.bind(number_spinbox,
                  ''.join(['Set the total number of wells for this wellbore ',
                           'component including random locations.']))

    # Wellbore locations frame
    well_locs_frame = tk.Frame(tab)
    well_locs_frame.grid(row=10, column=0, sticky='w',
                         padx=PARAMETER_FRAME_PADX, pady=(5, 10))
    add_wellbore_frame_widgets(controller, cmpnt_nm, well_locs_frame, tool_tip)

    # Outputs
    outputs_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    outputs_label.grid(row=11, column=0, sticky='w', pady=(5, 10))

    outputs_frame = ttk.Frame(tab)
    outputs_frame.grid(row=12, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    # Create and place outputs widgets
    output_nms_labels = []
    output_nms_checkboxes = []

    for ind1 in range(2):
        for ind2 in range(2):
            ind = ind1*2 + ind2
            obs_nm = OW_OBSERVATIONS[ind]
            # Create output checkbox
            output_nms_checkboxes.append(
                tk.Checkbutton(outputs_frame, variable=componentVars[cmpnt_nm][obs_nm]))
            # Place checkbox
            output_nms_checkboxes[-1].grid(
                row=ind1+1, column=2*ind2, pady=5, padx=CB_PADX, sticky='w')
            # Create output label
            output_nms_labels.append(
                ttk.Label(outputs_frame, text=OW_OBSERVATIONS_SETUP[obs_nm][0],
                          width=OUTPUT_LABEL_WIDTH1, anchor='w'))
            # Place label
            output_nms_labels[-1].grid(row=ind1+1, column=2*ind2+1, pady=5, sticky='w')
            # Bind checkbox to the tip
            tool_tip.bind(output_nms_checkboxes[-1],
                          OW_OBSERVATIONS_SETUP[obs_nm][1])

def use_crit_pressure_approach(variables, check_frame, brine_density_par_frame,
                               crit_pressure_par_frame):
    """ Disable/enable widgets for the relevant frames."""

    check_button_state = variables['critPressureApproach'].get()

    if not check_button_state:
        variables['enforceCritPressure'].set(0)
        for frame in [check_frame, brine_density_par_frame,
                      crit_pressure_par_frame]:
            for widget in frame.winfo_children():
                widget.configure(state='disabled')
    else:
        for frame in [check_frame, brine_density_par_frame]:
            for widget in frame.winfo_children():
                widget.configure(state='normal')

def enforce_crit_pressure_approach(variables, brine_density_par_frame,
                                   crit_pressure_par_frame):
    """ Disable/enable widgets for the relevant frames."""
    par_frames_state = {0: 'normal', 1: 'disabled'}
    check_button_state = variables['enforceCritPressure'].get()

    if variables['critPressureApproach'].get():
        for widget in brine_density_par_frame.winfo_children():
            widget.configure(state=par_frames_state[check_button_state])

        for widget in crit_pressure_par_frame.winfo_children():
            widget.configure(state=par_frames_state[1-check_button_state])

def process_crit_pressure_approach_pars(app, cmpnt_data, cmpnt_nm):
    if 'Controls' in cmpnt_data:
        critPressureApproach = cmpnt_data['Controls'].get(
            'critPressureApproach', False)
        enforceCritPressure = cmpnt_data['Controls'].get(
            'enforceCritPressure', False)
    else:
        critPressureApproach = False
        enforceCritPressure = False

    brine_density_par_frame_name = '.'.join([cmpnt_nm, 'brineDensity', 'frame'])
    brine_density_par_frame = app.nametowidget(
        app.getvar(brine_density_par_frame_name))
    crit_pressure_par_frame_name = '.'.join([cmpnt_nm, 'critPressure', 'frame'])
    crit_pressure_par_frame = app.nametowidget(
        app.getvar(crit_pressure_par_frame_name))

    for frame in [brine_density_par_frame,
                  crit_pressure_par_frame]:
        for widget in frame.winfo_children():
            widget.configure(state='disabled')

    if critPressureApproach:
        if enforceCritPressure:
            for widget in crit_pressure_par_frame.winfo_children():
                widget.configure(state='normal')
        else:
            for widget in brine_density_par_frame.winfo_children():
                widget.configure(state='normal')
