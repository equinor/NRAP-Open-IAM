"""
Module contains several methods needed for creating tab (page) in GUI
for PlumeStability component. Methods read and write dictionaries
needed for control file interface yaml files.
"""
import os
import sys
from re import split

import tkinter as tk
from tkinter import ttk
from tkinter import (StringVar, DoubleVar, BooleanVar, IntVar)
from tkinter.filedialog import askdirectory
from tkinter.filedialog import askopenfilename
from tkinter import messagebox

import numpy as np

import Pmw

SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(SOURCE_DIR)

from dictionarydata import componentVars, componentChoices
from dictionarydata import connectionsDictionary, componentTypeDictionary
from dictionarydata import LABEL_FONT
from dictionarydata import (DISTRIBUTION_MENU_WIDTH, DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            PARAMETER_LABEL_WIDTH, BUTTON_WIDTH,
                            OUTPUT_LABEL_WIDTH2, CB_PADX, FILE_ENTRY_WIDTH)
from cmpnts_tabs.parameter_entry import ParameterEntry
from openiam import IAM_DIR


PSA_OBSERVATIONS = ['areas', 'areas_dt', 'mobility', 'mobility_angles',
                    'spreading', 'spreading_angles']

PSA_OBS_UNITS = ['[m{}]'.format(u'\u00B2'), '[m{}/year]'.format(u'\u00B2'),
                 '[m/year]', '[-]',
                 '[m{}/year]'.format(u'\u00B2'), '[-]']

PSA_OBS_DESCRIPTIONS = {
    'pressure_areas': 'pressure plume area',
    'pressure_areas_dt': 'change in pressure plume area',
    'pressure_mobility': 'velocity of pressure plume centroid',
    'pressure_mobility_angles': 'direction of pressure plume centroid',
    'pressure_spreading': 'dispersion of pressure plume',
    'pressure_spreading_angles': 'direction of pressure plume dispersion',
    'CO2saturation_areas': 'CO{} saturation plume area'.format(u'\u2082'),
    'CO2saturation_areas_dt':
        'change in CO{} saturation plume area'.format(u'\u2082'),
    'CO2saturation_mobility':
        'velocity of CO{} saturation plume centroid'.format(u'\u2082'),
    'CO2saturation_mobility_angles':
        'direction of CO{} saturation plume centroid'.format(u'\u2082'),
    'CO2saturation_spreading':
        'dispersion of CO{} saturation plume'.format(u'\u2082'),
    'CO2saturation_spreading_angles':
        'direction of CO{} saturation plume dispersion'.format(u'\u2082')}


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    d = {}
    d['type'] = 'PlumeStability'
    d['connection'] = componentVars[cmpnt_nm]['connection'].get()

    d['FileDirectory'] = componentVars[cmpnt_nm]['file_dir'].get()
    d['TimeFile'] = componentVars[cmpnt_nm]['time_file'].get()
    d['ParameterFilename'] = componentVars[cmpnt_nm]['param_file'].get()

    d['Variables'] = []
    for var_nm in componentVars[cmpnt_nm]['variables']:
        if componentVars[cmpnt_nm][var_nm+'_cb_var'].get():
            d['Variables'].append(var_nm)

    d['Thresholds'] = {}
    for var_nm in d['Variables']:
        d['Thresholds'][var_nm] = componentVars[cmpnt_nm][var_nm+'_threshold'].get()

    pars_keys = ['index']

    d['Parameters'] = {}
    for key in pars_keys:
        d['Parameters'][key] = {}
        # Determine distribution of the parameter
        distr_type = componentVars[cmpnt_nm][key]['distribution'].get()
        if distr_type == 'Fixed Value':
            d['Parameters'][key]['value'] = (
                componentVars[cmpnt_nm][key]['value'].get())
            d['Parameters'][key]['vary'] = False

        if distr_type == 'Discrete':
            discrete_values = []
            discrete_weights = []

            for value in componentVars[cmpnt_nm][key]['values'].get().split(','):
                # We transform values to integer values since PSA has only one parameter
                discrete_values.append(int(value.strip()))

            for weight in componentVars[cmpnt_nm][key]['weights'].get().split(','):
                discrete_weights.append(float(weight.strip()))

            d['Parameters'][key]['Values'] = discrete_values
            d['Parameters'][key]['Weights'] = discrete_weights

    # Save observations
    model_outputs = []
    for obs_nm in componentVars[cmpnt_nm]['outputs']:
        if componentVars[cmpnt_nm]['outputs'][obs_nm].get():
            model_outputs.append(obs_nm)
    d['Outputs'] = model_outputs

    return d


def add_widgets(app, tab, cmpnt_nm, cmpnt_type, tool_tip, *args):
    """ Add widgets to the component tab."""

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('PlumeStability')
    connectionsDictionary.append(cmpnt_nm)

    componentVars[cmpnt_nm] = {}
    for key in ['connection', 'componentName', 'componentType']:
        componentVars[cmpnt_nm][key] = StringVar()

    componentVars[cmpnt_nm]['connection'].set('none')
    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type
    componentVars[cmpnt_nm]['variables'] = []
    componentVars[cmpnt_nm]['thresholds'] = {}

    for key in ['time_file', 'param_file', 'file_dir']:
        componentVars[cmpnt_nm][key] = StringVar()
        componentVars[cmpnt_nm][key].set('')

    app.psa_setup_frame = ttk.Frame(tab)
    app.psa_setup_frame.grid(row=4, column=0, columnspan=1, sticky='w', padx=0)

    # Tab title label
    comp_type_label = ttk.Label(tab, font=LABEL_FONT,
                                text="Plume Stability Analysis Component")
    comp_type_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    files_frame = ttk.Frame(tab)
    files_frame.grid(row=1, column=0, columnspan=6, sticky='w', padx=15)

    file_dir_label = ttk.Label(
        files_frame, text="Directory of input files:", width=PARAMETER_LABEL_WIDTH)
    file_dir_textfile = tk.Entry(
        files_frame, width=FILE_ENTRY_WIDTH,
        textvariable=componentVars[cmpnt_nm]['file_dir'])
    file_dir_button = ttk.Button(
        files_frame, text="Browse", width=BUTTON_WIDTH,
        command=lambda: choose_input_dir(componentVars[cmpnt_nm]))
    file_dir_label.grid(row=1, column=0, sticky='w', pady=5, padx=5)
    file_dir_textfile.grid(row=1, column=1, sticky='w', pady=5, padx=5)
    tool_tip.bind(file_dir_textfile,
                  'Enter directory that contains all data files.')
    file_dir_button.grid(row=1, column=2, padx=5)
    tool_tip.bind(file_dir_button,
                  'Select a directory containing all data files.')

    time_file_label = ttk.Label(
        files_frame, text="Input time points file:", width=PARAMETER_LABEL_WIDTH)
    time_file_textfile = tk.Entry(
        files_frame, width=FILE_ENTRY_WIDTH,
        textvariable=componentVars[cmpnt_nm]['time_file'])
    time_file_button = ttk.Button(
        files_frame, text="Browse", width=BUTTON_WIDTH,
        command=lambda: choose_time_points_file(componentVars[cmpnt_nm]))
    time_file_label.grid(row=2, column=0, sticky='w', pady=5, padx=5)
    time_file_textfile.grid(row=2, column=1, sticky='w', pady=5, padx=5)
    tool_tip.bind(time_file_textfile,
                  'Enter path to input time points file.')
    time_file_button.grid(row=2, column=2, padx=5)
    tool_tip.bind(time_file_button,
                  'Select input time points file.')

    param_file_label = ttk.Label(
        files_frame, text="Input parameters file:", width=PARAMETER_LABEL_WIDTH)
    param_file_textfile = tk.Entry(
        files_frame, width=FILE_ENTRY_WIDTH,
        textvariable=componentVars[cmpnt_nm]['param_file'])
    param_file_button = ttk.Button(
        files_frame, text="Browse", width=BUTTON_WIDTH,
        command=lambda: choose_parameters_filenames_file(
            app, componentVars[cmpnt_nm], app.psa_setup_frame))
    param_file_label.grid(row=3, column=0, sticky='w', pady=5, padx=5)
    param_file_textfile.grid(row=3, column=1, sticky='w', pady=5, padx=5)
    tool_tip.bind(param_file_textfile,
                  'Enter path to input parameters file.')
    param_file_button.grid(row=3, column=2, padx=5)
    tool_tip.bind(param_file_button,
                  'Select input parameters file.')


def get_possible_variables(file_directory, parameter_filename):
    """ Read data provided in input files and get names of provided observations
        for which variables can be created and for which metrics can be calculated.
    """
    # Read file with parameters and names of files with the corresponding data
    par_file_data = np.genfromtxt(os.path.join(file_directory, parameter_filename),
                                  delimiter=",", dtype='str')
    file_names = par_file_data[1:, -1] # last column has file names

    # Check whether the first column in the file parameter_filename has indices
    # of the lookup tables
    if par_file_data[0, 0] != 'index':
        raise NameError()

    # Get indices of linked lookup tables
    indices = [int(val) for val in par_file_data[1:, 0]]

    # Assume that all data files have identical variables
    with open(os.sep.join([file_directory, file_names[0]])) as fh:
        names = fh.readline().strip().split(',')
        names = names[2:] # exclude x and y

        # Remove non-time varying variables assuming that
        # they will not end in '_[0-9]+' (underscore followed by an integer)
        var_names = []
        for nm in names:
            split_result = split("_[0-9]+$", nm)
            if (len(split_result) == 2) and (split_result[0] not in var_names):
                var_names.append(split_result[0])

    return var_names, indices


def finish_load_setup(app, comp_data, cmpnt_nm):
    """ Add additional widgets after loading component data."""

    componentVars[cmpnt_nm]['time_file'].set(comp_data['TimeFile'])
    componentVars[cmpnt_nm]['param_file'].set(comp_data['ParameterFilename'])
    componentVars[cmpnt_nm]['file_dir'].set(comp_data['FileDirectory'])

    # Read user provided file with names of data files
    # Obtain names of variables for which the metrics can be calculated
    # and indices of linked lookup tables
    try:
        psa_variables, indices = get_possible_variables(
            os.path.join(IAM_DIR, comp_data['FileDirectory']),
            comp_data['ParameterFilename'])
    except NameError:
        full_par_file_path = os.path.join(
            IAM_DIR, comp_data['FileDirectory'], comp_data['ParameterFilename'])
        err_msg = "".join([
            "Please check setup of the component {}.",
            "\nThe first column of the file {} with path \n{} \nprovided in the field ",
            "'Input parameters file:' should contain the indices ",
            "of the lookup tables specified in the last column ('filename') and ",
            "be named 'index'.\nPlease update the file by adding ",
            "the corresponding column and try loading the simulation again. ",
            "Alternatively, one can download the updated data set from its original ",
            "location."]).format(
                cmpnt_nm, comp_data['ParameterFilename'], full_par_file_path)
        return 0, err_msg

    componentVars[cmpnt_nm]['variables'] = psa_variables.copy()

    psa_setup_subframes = add_psa_setup_frame_widgets(
        app, componentVars[cmpnt_nm], app.psa_setup_frame,
        psa_variables, indices)

    # Update values of global variables linked to the values from the tab textfields
    # Plume Stability component has only one parameter
    par_key = 'index'
    if 'value' in comp_data['Parameters'][par_key]:
        componentVars[cmpnt_nm][par_key]['distribution'].set(
            'Fixed Value')
        componentVars[cmpnt_nm][par_key]['value'].set(
            comp_data['Parameters'][par_key]['value'])

    elif 'Values' in comp_data['Parameters'][par_key]:
        componentVars[cmpnt_nm][par_key]['distribution'].set(
            'Discrete')
        values_list = comp_data['Parameters'][par_key][
            'Values']
        weights_list = comp_data['Parameters'][par_key][
            'Weights']
        componentVars[cmpnt_nm][par_key]['values'].set(
            ', '.join([str(val) for val in values_list]))
        componentVars[cmpnt_nm][par_key]['weights'].set(
            ', '.join([str(val) for val in weights_list]))
    change_psa_pars_distribution(
        app, par_key, componentVars[cmpnt_nm], psa_setup_subframes[par_key])

    # Setup global variables associated with PSA variables and their thresholds
    if 'Variables' in comp_data:
        for var_nm in psa_variables:
            if var_nm in comp_data['Variables']:
                componentVars[cmpnt_nm][var_nm+'_cb_var'].set(1)
            else:
                componentVars[cmpnt_nm][var_nm+'_cb_var'].set(0)
                enable_variable_metrics_calculation(
                    var_nm, psa_setup_subframes, componentVars[cmpnt_nm])

    if 'Thresholds' in comp_data:
        for var_nm in comp_data['Thresholds']:
            componentVars[cmpnt_nm][var_nm+'_threshold'].set(
                comp_data['Thresholds'][var_nm])

    # Setup outputs related widgets and related variables
    for obs_nm in comp_data['Outputs']:
        componentVars[cmpnt_nm]['outputs'][obs_nm].set(1)

    return 1, ''


def sample_all_psa_datafile_indices(app, frames, variables, indices):
    """ Check whether datafileIndexSample_checkbox is checked.

    If the button is checked, sampling of index will be done,
    and widgets corresponding to the parameters frames will be disabled.
    """
    # If check box is checked
    if variables['datafileSampleAll'].get():
        variables['index']['distribution'].set('Discrete')
        str_values = ', '.join([str(val) for val in indices])
        variables['index']['values'].set(str_values)
        str_weights = ', '.join(len(indices)*['1'])
        variables['index']['weights'].set(str_weights)
        change_psa_pars_distribution(app, 'index', variables,
                                     frames['index'])


def change_psa_pars_distribution(app, par_key, variables, frame):
    """ Change distribution for the parameters of the plume stability analysis component."""
    tool_tip = Pmw.Balloon(app)

    for widget in frame.winfo_children():
        widget.destroy()

    # Create and configure common widgets
    label = ttk.Label(frame, width=PARAMETER_LABEL_WIDTH, text=frame.labelText)
    value_menu = tk.OptionMenu(
        frame, variables[par_key]['distribution'], *['Fixed Value', 'Discrete'],
        command=lambda _: change_psa_pars_distribution(
            app, par_key, variables, frame))
    label.grid(row=0, column=0, sticky='w', padx=5)
    value_menu.grid(row=0, column=1, sticky='ew', padx=5)
    value_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    tool_tip.bind(value_menu, 'Select distribution for {}.'.format(
        frame.distToolTipText))

    # If distribution is chosen to be deterministic
    # we don't want user to be confused with sample all checkbox being selected
    if variables[par_key]['distribution'].get() == 'Fixed Value':
        variables['datafileSampleAll'].set(0)

    app.add_remaining_widgets(
        frame.toolTipText, variables[par_key]['distribution'].get(), frame,
        tool_tip, variables[par_key])


def choose_input_dir(variables):
    """
    Set the input directory for the plume stability analysis component.
    """
    file_dialog = tk.Tk()

    file_dialog.withdraw()

    # Determine initial directory for file dialog
    initial_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), 'components', 'reservoir', 'lookuptables')

    if variables['file_dir'].get():
        initial_dir = os.path.join(IAM_DIR, variables['file_dir'].get())

    try:
        directory = askdirectory(
            initialdir=initial_dir,
            title="Choose input directory")
    except:
        file_dialog.destroy()
    else:
        if directory != '':
            variables['file_dir'].set(directory)

        file_dialog.destroy()


def choose_time_points_file(variables):
    """
    Set the file containing time points at which data for
    plume stability analysis component is available.
    """
    file_dialog = tk.Tk()

    file_dialog.withdraw()

    # Determine initial directory for file dialog
    initial_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), 'components', 'reservoir', 'lookuptables')

    if variables['file_dir'].get():
        initial_dir = os.path.join(IAM_DIR, variables['file_dir'].get())

    try:
        time_file_name = askopenfilename(
            initialdir=initial_dir,
            title="Choose file with time points")
    except:
        file_dialog.destroy()
    else:
        if time_file_name != '':
            variables['time_file'].set(time_file_name)

        file_dialog.destroy()


def choose_parameters_filenames_file(app, variables, psa_setup_frame):
    """
    Read file containing parameters values and corresponding files with data.

    After reading the file, the weights are calculated for each input
    parameter value that exists for each parameter.
    """
    file_dialog = tk.Tk()

    file_dialog.withdraw()

    # Determine initial directory for file dialog
    initial_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), 'components', 'reservoir', 'lookuptables')

    if variables['file_dir'].get():
        initial_dir = os.path.join(IAM_DIR, variables['file_dir'].get())

    try:
        parameters_file_name = askopenfilename(
            initialdir=initial_dir,
            title="Choose parameters file")
    except:
        file_dialog.destroy()
    else:
        try:
            # Determine variables and indices of the lookup tables
            psa_variables, indices = get_possible_variables(
                os.path.join(IAM_DIR, variables['file_dir'].get()),
                parameters_file_name)
        except NameError:
            err_msg = "".join([
                "The first column of the file {} should contain the indices ",
                "of the lookup tables specified in the last column (filename) and ",
                "be named 'index'. Please update the file by adding ",
                "the corresponding column and try again. Alternatively, ",
                "one can download the updated data set from its original ",
                "location."]).format(parameters_file_name)
            messagebox.showerror("Error", err_msg)
            return

        variables['variables'] = psa_variables.copy()

        # Add parameters widgets
        psa_setup_subframes = add_psa_setup_frame_widgets(
            app, variables, psa_setup_frame, psa_variables, indices)

        for var_nm in psa_variables:
            # Check all boxes corresponding to the available variables
            variables[var_nm+'_cb_var'].set(1)

        # Save the name of parameter file in the global variables
        variables['param_file'].set(parameters_file_name)
        file_dialog.destroy()


def enable_variable_metrics_calculation(var_nm, psa_setup_subframes, variables):
    """ Enable widgets for outputs corresponding to a particular variable. """
    # If checkbox corresponding to the variable is checked

    states = {0: 'normal', 1: 'disabled'}

    check_button_state = variables[var_nm+'_cb_var'].get()

    # Enable variable frame widgets
    psa_setup_subframes[var_nm].var_label.configure(
        state=states[1-check_button_state])
    psa_setup_subframes[var_nm].threshold_label.configure(
        state=states[1-check_button_state])
    psa_setup_subframes[var_nm].threshold_field.configure(
        state=states[1-check_button_state])

    # Uncheck observation checkboxes corresponding to the unchecked variables
    if not check_button_state:
        for base_obs_nm in PSA_OBSERVATIONS:
            obs_nm = '{}_{}'.format(var_nm, base_obs_nm)
            variables['outputs'][obs_nm].set(0)

    # Enable/disable widgets on the output frame corresponding to the variable
    for widget in psa_setup_subframes[var_nm+'_outputs'].winfo_children():
        widget.configure(state=states[1-check_button_state])


def get_var_unit(var_nm):
    var_unit = ' [-]'

    if var_nm == 'pressure':
        var_unit = ' [Pa]'

    return var_unit


def add_psa_setup_frame_widgets(app, variables, psa_setup_frame,
                                psa_variables, indices):
    """
    Add widgets after loading file containing information about PSA parameters
    and related data files.
    """
    tool_tip = Pmw.Balloon(app)

    for widget in psa_setup_frame.winfo_children():
        widget.destroy()

    psa_setup_subframes = {}

    # Define possible values of index parameter
    values = {'index': indices}
    variables['Params'] = ['index']

    # Set variables responsible for sampling all linked lookup tables
    variables['datafileSampleAll'] = BooleanVar()
    variables['datafileSampleAll'].set(0)

    # Dictionary with default values of index
    if len(values['index']) > 1:
        setup_dict = {
            'distribution': 'Fixed Value',
            'values': '{}, {}'.format(
                values['index'][0], values['index'][1]),
            'weights': '0.5, 0.5'}
    else:
        setup_dict = {
            'distribution': 'Fixed Value',
            'values': '{}'.format(values['index'][0]),
            'weights': '1'}

    # Define and setup variables
    variables['index'] = {}
    for item in setup_dict:
        variables['index'][item] = StringVar()
        variables['index'][item].set(setup_dict[item])

    # Set variable for deterministic scenario
    variables['index']['value'] = IntVar()
    variables['index']['value'].set(int(values['index'][0]))

    # Create parameter frame
    psa_setup_subframes['index'] = tk.Frame(psa_setup_frame)
    label_text1 = "Data file index:"
    label_text2 = 'data file index'
    label_text3 = "data file index"

    # Place frame on grid
    psa_setup_subframes['index'].grid(
        row=0, column=0, columnspan=6, sticky='w', padx=15)

    # Configure frame attributes
    psa_setup_subframes['index'].labelText = label_text1
    psa_setup_subframes['index'].toolTipText = label_text2
    psa_setup_subframes['index'].distToolTipText = label_text3
    psa_setup_subframes['index'].distType = (
        variables['index']['distribution'])

    psa_setup_subframes['index'].par_bounds = {
        'lower_bound': None, 'upper_bound': None,
        'discrete_bounds': indices}

    # Define function for option menu
    def get_psa_distr_fun(key):
        return lambda _: change_psa_pars_distribution(
            app, key, variables, psa_setup_subframes[key])

    # Setup widgets for datafile index
    par_label = ttk.Label(psa_setup_subframes['index'],
                          text=label_text1, width=PARAMETER_LABEL_WIDTH)

    par_menu = tk.OptionMenu(psa_setup_subframes['index'],
                             variables['index']['distribution'],
                             *['Fixed Value', 'Discrete'],
                             command=get_psa_distr_fun('index'))

    value_label = ttk.Label(psa_setup_subframes['index'], text='Value:',
                            width=DISTRIBUTION_ARG_LABEL_WIDTH)

    value_field = ParameterEntry(
        psa_setup_subframes['index'], 'index',
        variables['index']['value'],
        DISTRIBUTION_ARG_TEXTFIELD_WIDTH, tool_tip,
        standard_tooltip_text=''.join([
            'Set value of data file index.\nPossible values are \n{}.'.format(
                app.reformat_list_presentation(indices))]),
        **psa_setup_subframes['index'].par_bounds)

    # Place widgets on the frame
    par_label.grid(row=0, column=0, sticky='w', padx=5)
    par_menu.grid(row=0, column=1, sticky='ew', padx=5)
    value_label.grid(row=0, column=2, sticky='e', padx=5, pady=0)
    value_field.grid(row=0, column=3, sticky='e', padx=5, pady=0)

    # Configure the widths of menu widgets and entry fields
    par_menu.config(width=DISTRIBUTION_MENU_WIDTH)

    # Assign tooltips for the menu widgets
    tool_tip.bind(par_menu, 'Select distribution for {}.'.format(
        psa_setup_subframes['index'].distToolTipText))

    # Sampling all data files frame
    psa_setup_subframes['sample_all'] = tk.Frame(psa_setup_frame)
    psa_setup_subframes['sample_all'].grid(
        row=1, column=0, columnspan=6, sticky='w', padx=15)
    datafile_sample_all_label = ttk.Label(
        psa_setup_subframes['sample_all'], text="Sample all data file indices:",
        width=PARAMETER_LABEL_WIDTH)
    datafile_sample_all_checkbox = tk.Checkbutton(
        psa_setup_subframes['sample_all'],
        variable=variables['datafileSampleAll'],
        command=lambda: sample_all_psa_datafile_indices(
            app, psa_setup_subframes, variables, indices))
    datafile_sample_all_label.grid(
        row=0, column=0, padx=5, sticky='w')
    datafile_sample_all_checkbox.grid(
        row=0, column=1, padx=5, sticky='e')
    tool_tip.bind(datafile_sample_all_checkbox,
                  'Check to sample indices of all linked lookup tables.')

    # Variables label frame
    psa_setup_subframes['variables'] = tk.Frame(psa_setup_frame)
    psa_setup_subframes['variables'].grid(
        row=2, column=0, columnspan=6, sticky='w', padx=15)
    variables_label = ttk.Label(psa_setup_subframes['variables'], text='Variables:',
                                width=16)
    variables_label.grid(row=0, column=0, sticky='w', padx=5, pady=5)

    # Define function for checkboxes associated with PSA variables
    def get_variable_obs_fun(nm):
        return lambda: enable_variable_metrics_calculation(
            nm, psa_setup_subframes, variables)

    # Frames for variables available in the linked lookup tables
    for ind, var_nm in enumerate(psa_variables):
        psa_setup_subframes[var_nm] = tk.Frame(psa_setup_frame)
        psa_setup_subframes[var_nm].grid(
            row=3+ind, column=0, columnspan=3, sticky='w', padx=15)

        # tkinter variables to keep track of what PSA variables are loaded
        variables[var_nm+'_cb_var'] = BooleanVar()
        variables[var_nm+'_cb_var'].set(0)
        variables[var_nm+'_threshold'] = DoubleVar()
        variables[var_nm+'_threshold'].set(0.0)

        var_checkbox = tk.Checkbutton(
            psa_setup_subframes[var_nm],
            variable=variables[var_nm+'_cb_var'],
            command=get_variable_obs_fun(var_nm))

        var_label = ttk.Label(psa_setup_subframes[var_nm], text=var_nm,
                              width=20)

        var_threshold_label = ttk.Label(
            psa_setup_subframes[var_nm], text="Threshold{}:".format(
                get_var_unit(var_nm)),
            width=15)

        var_threshold_field = tk.Entry(psa_setup_subframes[var_nm],
                                       textvariable=variables[var_nm+'_threshold'])
        var_threshold_field.config(width=DISTRIBUTION_ARG_TEXTFIELD_WIDTH)

        # Save widgets
        psa_setup_subframes[var_nm].var_label = var_label
        psa_setup_subframes[var_nm].threshold_label = var_threshold_label
        psa_setup_subframes[var_nm].threshold_field = var_threshold_field

        # Place widgets on psa_setup_subframes[var_nm] frame
        var_checkbox.grid(row=0, column=0, sticky='w', padx=(15, 5))
        var_label.grid(row=0, column=1, sticky='w', padx=(0, 5))
        var_threshold_label.grid(row=0, column=2, sticky='e', padx=5)
        var_threshold_field.grid(row=0, column=3, sticky='w', padx=5)
        if 'CO2' in var_nm:
            tool_tip_var_nm = var_nm.replace('CO2', 'CO{} '.format(u'\u2082'))
        else:
            tool_tip_var_nm = var_nm
        tool_tip.bind(var_checkbox,
                      'Check to enable calculation of metrics for {}.'.format(
                          tool_tip_var_nm))
        tool_tip.bind(var_threshold_field,
                      'Select value of threshold for {}.'.format(
                          tool_tip_var_nm))

    # Outputs label frame
    frame_ind = len(psa_variables) + 3
    psa_setup_subframes['outputs'] = tk.Frame(psa_setup_frame)
    psa_setup_subframes['outputs'].grid(
        row=frame_ind, column=0, columnspan=6, sticky='w', padx=0)
    outputs_label = ttk.Label(
        psa_setup_subframes['outputs'], text='Outputs', font=LABEL_FONT)
    outputs_label.grid(row=0, column=0, sticky='w', pady=(5, 10), padx=0)

    variables['outputs'] = {}
    # Variables associated observation frames
    for ind1, var_nm in enumerate(psa_variables):
        psa_setup_subframes[var_nm+'_outputs'] = tk.Frame(psa_setup_frame)
        psa_setup_subframes[var_nm+'_outputs'].grid(
            row=frame_ind+ind1+1, column=0, columnspan=6, sticky='w', padx=15)

        for ind2, base_obs_nm in enumerate(PSA_OBSERVATIONS):
            obs_nm = '{}_{}'.format(var_nm, base_obs_nm)
            variables['outputs'][obs_nm] = BooleanVar()
            variables['outputs'][obs_nm].set(0)
            obs_checkbox = tk.Checkbutton(
                psa_setup_subframes[var_nm+'_outputs'],
                variable=variables['outputs'][obs_nm])
            label_text = '{} {} {}'.format(var_nm, base_obs_nm, PSA_OBS_UNITS[ind2])
            if 'CO2' in var_nm:
                label_text = label_text.replace('CO2', 'CO{} '.format(u'\u2082'))
            elif 'pressure' in var_nm:
                label_text = label_text.replace('pressure', 'Pressure')
            elif 'salinity' in var_nm:
                label_text = label_text.replace('salinity', 'Salinity')
            obs_label = ttk.Label(psa_setup_subframes[var_nm+'_outputs'],
                                  text=label_text,
                                  width=OUTPUT_LABEL_WIDTH2, anchor='w')
            obs_checkbox.grid(row=ind2//3, column=2*(ind2%3),
                              pady=5, padx=CB_PADX, sticky='w')
            obs_label.grid(row=ind2//3, column=2*(ind2%3)+1,
                           pady=5, sticky='w')
            if obs_nm in PSA_OBS_DESCRIPTIONS:
                tool_tip.bind(obs_checkbox,
                              'Enable {} as output.'.format(PSA_OBS_DESCRIPTIONS[obs_nm]))
            else:
                tool_tip.bind(obs_checkbox,
                              'Enable {} {} as output.'.format(var_nm, base_obs_nm))

    return psa_setup_subframes


if __name__ == "__main__":
    file_directory = os.sep.join(['..', '..', 'components', 'reservoir',
                                  'lookuptables', 'Kimb_54_sims'])
    parameter_filename = 'parameters_and_filenames.csv'
    var_names, indices = get_possible_variables(file_directory, parameter_filename)
    print(var_names)
    print(indices)
