"""
Module contains several methods needed for creating tab (page) in GUI
for LookupTableReservoir component. Methods read and write dictionaries
needed for control file interface yaml files.
"""
import os
import sys
import csv

import tkinter as tk
from tkinter import ttk
from tkinter import (StringVar, DoubleVar, IntVar, BooleanVar)
from tkinter.filedialog import askdirectory
from tkinter.filedialog import askopenfilename
from tkinter import messagebox

import Pmw

SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(SOURCE_DIR)

from dictionarydata import componentVars, componentChoices
from dictionarydata import connectionsDictionary, componentTypeDictionary

from dictionarydata import LABEL_FONT
from dictionarydata import (DISTRIBUTION_MENU_WIDTH, DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            PARAMETER_LABEL_WIDTH, BUTTON_WIDTH,
                            OUTPUT_LABEL_WIDTH1, PARAMETER_FRAME_PADX,
                            CB_PADX, FILE_ENTRY_WIDTH)

from cmpnts_tabs.locations import (add_obs_locs_frame_widgets,
                                   read_obs_locations_data)

from cmpnts_tabs.parameter_entry import ParameterEntry

from openiam import IAM_DIR


def read_tab_vars(cmpnt_nm):
    """ Read values of tkinter variables associated with the component tab."""
    d = {}
    d['type'] = 'LookupTableReservoir'
    d['connection'] = componentVars[cmpnt_nm]['connection'].get()

    d['FileDirectory'] = componentVars[cmpnt_nm]['file_dir'].get()
    d['TimeFile'] = componentVars[cmpnt_nm]['time_file'].get()
    d['ParameterFilename'] = componentVars[cmpnt_nm]['param_file'].get()

    # Check what parameters will be added to the parameters dictionary
    if componentVars[cmpnt_nm]['parSample'].get() == 1:
        # If sampling of lookup data files index is turned off
        if 'Params' in componentVars[cmpnt_nm]:
            pars_keys = componentVars[cmpnt_nm]['Params'][1:]
        else:
            pars_keys = {}
    else:
        # If sampling of data file index is needed
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

            if key == 'index':
                for value in componentVars[cmpnt_nm][key]['values'].get().split(','):
                    discrete_values.append(int(value.strip()))
            else:
                for value in componentVars[cmpnt_nm][key]['values'].get().split(','):
                    discrete_values.append(float(value.strip()))

            for weight in componentVars[cmpnt_nm][key]['weights'].get().split(','):
                discrete_weights.append(float(weight.strip()))

            d['Parameters'][key]['Values'] = discrete_values
            d['Parameters'][key]['Weights'] = discrete_weights

    read_obs_locations_data(d, cmpnt_nm)

    model_outputs = []
    if componentVars[cmpnt_nm]['pressure'].get():
        model_outputs.append('pressure')
    if componentVars[cmpnt_nm]['CO2saturation'].get():
        model_outputs.append('CO2saturation')

    d['Outputs'] = model_outputs

    return d


def add_widgets(app, tab, cmpnt_nm, cmpnt_type, tool_tip, *args):
    """ Add widgets to the component tab."""

    componentChoices.append(cmpnt_nm)
    componentTypeDictionary.append('LookupTableReservoir')
    connectionsDictionary.append(cmpnt_nm)

    componentVars[cmpnt_nm] = {}
    for key in ['connection', 'componentName', 'componentType']:
        componentVars[cmpnt_nm][key] = StringVar()

    componentVars[cmpnt_nm]['connection'].set('none')
    componentVars[cmpnt_nm]['componentName'] = cmpnt_nm
    componentVars[cmpnt_nm]['componentType'] = cmpnt_type

    for key in ['time_file', 'param_file', 'file_dir']:
        componentVars[cmpnt_nm][key] = StringVar()
        componentVars[cmpnt_nm][key].set('')

    for key in ['pressure', 'CO2saturation']:
        componentVars[cmpnt_nm][key] = BooleanVar()
        componentVars[cmpnt_nm][key].set(0)

    app.lutr_setup_frame = ttk.Frame(tab)
    app.lutr_setup_frame.grid(row=5, column=0, columnspan=1, sticky='w',
                              padx=PARAMETER_FRAME_PADX)

    # Observation locations frame
    app.obs_locs_frame = tk.Frame(tab)
    app.obs_locs_frame.grid(row=4, column=0, sticky='w',
                            padx=PARAMETER_FRAME_PADX)
    add_obs_locs_frame_widgets(app, cmpnt_nm, app.obs_locs_frame, tool_tip)
    # TODO Disable locations for now

    # Tab title label
    lookup_table_tab_label = ttk.Label(tab, font=LABEL_FONT,
                                       text="Lookup Table Reservoir Component")
    lookup_table_tab_label.grid(row=0, column=0, sticky='w', pady=(5, 10))

    files_frame = ttk.Frame(tab)
    files_frame.grid(row=1, column=0, columnspan=6, sticky='w',
                     padx=PARAMETER_FRAME_PADX)

    lookup_table_file_dir_label = ttk.Label(files_frame,
                                            text="Directory of input files:",
                                            width=PARAMETER_LABEL_WIDTH)
    lookup_table_file_dir_textfile = tk.Entry(
        files_frame, width=FILE_ENTRY_WIDTH,
        textvariable=componentVars[cmpnt_nm]['file_dir'])
    lookup_table_file_dir_button = ttk.Button(
        files_frame, text="Browse", width=BUTTON_WIDTH,
        command=lambda: choose_input_dir(componentVars[cmpnt_nm]))
    lookup_table_file_dir_label.grid(row=1, column=0, sticky='w', pady=5, padx=5)
    lookup_table_file_dir_textfile.grid(row=1, column=1, sticky='w', pady=5, padx=5)
    tool_tip.bind(lookup_table_file_dir_textfile,
                  'Enter directory that contains all lookup table files.')
    lookup_table_file_dir_button.grid(row=1, column=2, padx=5)
    tool_tip.bind(lookup_table_file_dir_button,
                  'Select a directory containing all lookup table files.')

    lookup_table_time_file_label = ttk.Label(
        files_frame, text="Input time points file:", width=PARAMETER_LABEL_WIDTH)
    lookup_table_time_file_textfile = tk.Entry(
        files_frame, width=FILE_ENTRY_WIDTH,
        textvariable=componentVars[cmpnt_nm]['time_file'])
    lookup_table_time_file_button = ttk.Button(
        files_frame, text="Browse", width=BUTTON_WIDTH,
        command=lambda: choose_time_points_file(componentVars[cmpnt_nm]))
    lookup_table_time_file_label.grid(row=2, column=0, sticky='w', pady=5, padx=5)
    lookup_table_time_file_textfile.grid(row=2, column=1, sticky='w', pady=5, padx=5)
    tool_tip.bind(lookup_table_time_file_textfile,
                  'Enter path to input time points file.')
    lookup_table_time_file_button.grid(row=2, column=2, padx=5)
    tool_tip.bind(lookup_table_time_file_button,
                  'Select input time points file.')

    lookup_table_param_file_label = ttk.Label(
        files_frame, text="Input parameters file:", width=PARAMETER_LABEL_WIDTH)
    lookup_table_param_file_textfile = tk.Entry(
        files_frame, width=FILE_ENTRY_WIDTH,
        textvariable=componentVars[cmpnt_nm]['param_file'])
    lookup_table_param_file_button = ttk.Button(
        files_frame, text="Browse", width=BUTTON_WIDTH,
        command=lambda: choose_parameters_filenames_file(
            app, componentVars[cmpnt_nm], app.lutr_setup_frame))
    lookup_table_param_file_label.grid(row=3, column=0, sticky='w', pady=5, padx=5)
    lookup_table_param_file_textfile.grid(row=3, column=1, sticky='w', pady=5, padx=5)
    tool_tip.bind(lookup_table_param_file_textfile,
                  'Enter path to input parameters file.')
    lookup_table_param_file_button.grid(row=3, column=2, padx=5)
    tool_tip.bind(lookup_table_param_file_button,
                  'Select input parameters file.')

    # Outputs
    output_label = ttk.Label(tab, text="Outputs", font=LABEL_FONT)
    output_label.grid(row=6, column=0, sticky='w', pady=(5, 10))

    output_frame = ttk.Frame(tab)
    output_frame.grid(row=7, column=0, sticky='w', padx=PARAMETER_FRAME_PADX)

    lookup_table_pressure_label = ttk.Label(
        output_frame, text="Pressure [Pa]", width=OUTPUT_LABEL_WIDTH1, anchor='w')
    lookup_table_pressure_checkbox = tk.Checkbutton(
        output_frame, variable=componentVars[cmpnt_nm]['pressure'])
    lookup_table_pressure_label.grid(row=1, column=1, pady=5, sticky='w')
    lookup_table_pressure_checkbox.grid(
        row=1, column=0, pady=5, padx=CB_PADX, sticky='w')
    tool_tip.bind(lookup_table_pressure_checkbox, 'Enable pressure as output.')

    lookup_table_co2sat_label = ttk.Label(
        output_frame, text="CO{} saturation [-]".format(u'\u2082'),
        width=OUTPUT_LABEL_WIDTH1, anchor='w')
    lookup_table_co2sat_checkbox = tk.Checkbutton(
        output_frame, variable=componentVars[cmpnt_nm]['CO2saturation'])
    lookup_table_co2sat_label.grid(row=1, column=3, pady=5, sticky='w')
    lookup_table_co2sat_checkbox.grid(
        row=1, column=2, pady=5, padx=CB_PADX, sticky='w')
    tool_tip.bind(lookup_table_co2sat_checkbox,
                  'Enable CO{} saturation as output.'.format(u'\u2082'))

def finish_load_setup(app, comp_data, cmpnt_nm):
    """ Add additional widgets after loading component data."""

    componentVars[cmpnt_nm]['time_file'].set(comp_data['TimeFile'])
    componentVars[cmpnt_nm]['param_file'].set(comp_data['ParameterFilename'])
    componentVars[cmpnt_nm]['file_dir'].set(comp_data['FileDirectory'])

    for widget in app.lutr_setup_frame.winfo_children():
        widget.destroy()

    full_par_file_path = os.path.join(
         IAM_DIR, comp_data['FileDirectory'], comp_data['ParameterFilename'])
    dictionary = csv.DictReader(open(full_par_file_path))

    try:
        lut_pars_frames = add_lut_pars_widgets(
            app, componentVars[cmpnt_nm], app.lutr_setup_frame, dictionary)
    except NameError:
        err_msg = "".join([
            "Please check setup of the component {}. ",
            "\nThe first column of the file {} with path \n{} \nprovided in the field ",
            "'Input parameters file:' should contain the indices ",
            "of the lookup tables specified in the last column ('filename') and ",
            "be named 'index'.\nPlease update the file by adding ",
            "the corresponding column and try loading the simulation again. ",
            "Alternatively, one can download the updated data set from its original ",
            "location."]).format(
                cmpnt_nm, comp_data['ParameterFilename'], full_par_file_path)

        return 0, err_msg

    for par_key in comp_data['Parameters']:
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
        change_lut_pars_distribution(
            app, par_key, componentVars[cmpnt_nm], lut_pars_frames[par_key])

    if 'index' in comp_data['Parameters']:
        componentVars[cmpnt_nm]['parSample'].set(0)
        for par_key in lut_pars_frames:
            if par_key != 'index':
                for widget in lut_pars_frames[par_key].winfo_children():
                    widget.configure(state='disabled')
    else:
        for widget in lut_pars_frames['index'].winfo_children():
            widget.configure(state='disabled')
        componentVars[cmpnt_nm]['parSample'].set(1)

    if 'pressure' in comp_data['Outputs']:
        componentVars[cmpnt_nm]['pressure'].set(1)
    if 'CO2saturation' in comp_data['Outputs']:
        componentVars[cmpnt_nm]['CO2saturation'].set(1)

    return 1, ''

def sample_lut_parameters(variables, frames):
    """ Check whether par_sample_checkbox is checked.

    If the button is checked, sampling of parameters will be done,
    and widgets corresponding to the index frame will be disabled.
    """
    par_frames_state = {0: 'normal', 1: 'disabled'}

    check_button_state = variables['parSample'].get()
    if check_button_state:  # means sampling of parameters is needed
        variables['indexSampleAll'].set(0)  # sampling of index is not needed

    for key in frames:
        if key != 'index':
            for widget in frames[key].winfo_children():
                widget.configure(state=par_frames_state[1-check_button_state])
        else:
            for widget in frames[key].winfo_children():
                widget.configure(state=par_frames_state[check_button_state])


def sample_all_lut_indices(app, frames, variables, index_values):
    """ Check whether sample_all_checkbox is checked.

    If the button is checked, sampling of index will be done,
    and widgets corresponding to the parameters frames will be disabled.
    """
    # If check box is checked
    if variables['indexSampleAll'].get():
        variables['parSample'].set(0)
        sample_lut_parameters(variables, frames)
        variables['index']['distribution'].set('Discrete')
        str_values = ', '.join([str(val) for val in index_values])
        variables['index']['values'].set(str_values)
        str_weights = ', '.join(len(index_values)*['1'])
        variables['index']['weights'].set(str_weights)
        change_lut_pars_distribution(app, 'index', variables,
                                     frames['index'])


def change_lut_pars_distribution(app, par_key, variables, frame):
    """ Change distribution for the parameters of the LUT reservoir component."""
    tool_tip = Pmw.Balloon(app)

    for widget in frame.winfo_children():
        widget.destroy()

    # Create and configure common widgets
    label = ttk.Label(frame, width=PARAMETER_LABEL_WIDTH, text=frame.labelText)
    value_menu = tk.OptionMenu(
        frame, variables[par_key]['distribution'], *['Fixed Value', 'Discrete'],
        command=lambda _: change_lut_pars_distribution(app, par_key, variables, frame))
    label.grid(row=0, column=0, sticky='w', padx=5)
    value_menu.grid(row=0, column=1, sticky='ew', padx=5)
    value_menu.config(width=DISTRIBUTION_MENU_WIDTH)
    tool_tip.bind(value_menu, 'Select distribution for {}.'.format(
        frame.distToolTipText))

    # If distribution is chosen to be deterministic
    # we don't want user to be confused with "sample all" checkbox being selected
    if (par_key == 'index') and (
            variables[par_key]['distribution'].get() == 'Fixed Value'):
        variables['indexSampleAll'].set(0)

    app.add_remaining_widgets(
        frame.toolTipText, variables[par_key]['distribution'].get(), frame,
        tool_tip, variables[par_key])


def choose_input_dir(variables):
    """
    Set the input directory for the LUT reservoir component.
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
    Set the file containing time points at which the data for LUT component
    is available.
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
            variables['time_file'].set(os.path.abspath(time_file_name))

            if not variables['file_dir'].get():
                variables['file_dir'].set(
                    os.path.abspath(os.path.dirname(time_file_name)))

        file_dialog.destroy()


def choose_parameters_filenames_file(app, variables, lutr_setup_frame):
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
        if parameters_file_name != '':
            # Read the content of the file
            dictionary = csv.DictReader(open(parameters_file_name))

            # Add parameters widgets
            try:
                lut_pars_frames = add_lut_pars_widgets(
                    app, variables, lutr_setup_frame, dictionary)
            except NameError:
                err_msg = "".join([
                    "The first column of the file {} should contain the indices ",
                    "of the lookup tables specified in the last column (filename) and ",
                    "be named 'index'. Please update the file by adding ",
                    "the corresponding column and try again. Alternatively, ",
                    "one can download the updated data set from its original ",
                    "location."]).format(parameters_file_name)
                messagebox.showerror("Error: Update Required", err_msg)
                return

            # Uncheck checkbox corresponding to parameters sampling
            variables['parSample'].set(0)

            # Disable widgets corresponding to the lookup table data parameters sampling
            for key in lut_pars_frames:
                if key != 'index':
                    for widget in lut_pars_frames[key].winfo_children():
                        widget.configure(state='disabled')
                else:
                    for widget in lut_pars_frames['index'].winfo_children():
                        widget.configure(state='normal')

            variables['param_file'].set(os.path.abspath(parameters_file_name))

            if not variables['file_dir'].get():
                variables['file_dir'].set(
                    os.path.abspath(os.path.dirname(parameters_file_name)))

        file_dialog.destroy()


def add_lut_pars_widgets(app, variables, lutr_setup_frame, dictionary):
    """
    Add widgets after loading file containing information about LUT parameters
    and related LUTs.
    """
    tool_tip = Pmw.Balloon(app)

    for widget in lutr_setup_frame.winfo_children():
        widget.destroy()

    # Determine keys and check whether the first key is index
    # It is possible index is the only one
    keys = dictionary.fieldnames[0:-1]
    if keys[0] != 'index':
        # Raise error that will be dealt accordingly to the scenario.
        raise NameError()

    values = {}

    for key in keys:
        values[key] = []

    # Determine unique values for each parameter including index
    for row in dictionary:
        for key in keys:
            try:
                # Check if the value is already in the list of unique values
                values[key].index(row[key])
            except ValueError:
                # If we find a new value add it to the list
                values[key].append(row[key])
            else:
                pass

    lut_pars_frames = {}

    def get_lut_distr_fun(key):
        return lambda _: change_lut_pars_distribution(
            app, key, variables, lut_pars_frames[key])

    variables['Params'] = list(keys)

    # Cycle over all parameters starting with the index parameter
    for i, key in enumerate(keys):
        # Define variables associated with widgets
        variables[key] = {}

        if len(values[key]) > 1:
            setup_dict = {
                'distribution': 'Fixed Value',
                'values': '{}, {}'.format(
                    values[key][0], values[key][1]),
                'weights': '0.5, 0.5'}
        else:
            setup_dict = {
                'distribution': 'Fixed Value',
                'values': '{}'.format(values[key][0]),
                'weights': '1'}

        # Define and setup variables
        for item in setup_dict:
            variables[key][item] = StringVar()
            variables[key][item].set(setup_dict[item])

        # Set variable for deterministic scenario
        if key == 'index':
            variables[key]['value'] = IntVar()
        else:
            variables[key]['value'] = DoubleVar()
        variables[key]['value'].set(values[key][0])

        # Create parameter frame
        lut_pars_frames[key] = tk.Frame(lutr_setup_frame)
        if i > 0:
            row_ind = 2+i
            label_text1 = "Parameter {} ({}):".format(i, key.strip())
            float_vals = [float(val) for val in values[key]]
            label_text2 = 'parameter {}'.format(key.strip())
            label_text3 = "parameter {}".format(key.strip())
            par_bounds_kwargs = {'lower_bound': None, 'upper_bound': None,
                                 'discrete_bounds': float_vals}
            tool_tip_end = '\nPossible values are \n{}.'.format(
                app.reformat_list_presentation(float_vals))

        else:
            row_ind = 1
            int_vals = [int(val) for val in values[key]]
            label_text1 = "Data file index:"
            label_text2 = 'data file index'
            label_text3 = "data file index"
            par_bounds_kwargs = {'lower_bound': None, 'upper_bound': None,
                                 'discrete_bounds': int_vals}
            tool_tip_end = '\nPossible values are \n{}.'.format(
                app.reformat_list_presentation(int_vals))

        lut_pars_frames[key].grid(
            row=row_ind, column=0, columnspan=6, sticky='w')
        lut_pars_frames[key].labelText = label_text1
        lut_pars_frames[key].toolTipText = label_text2
        lut_pars_frames[key].distToolTipText = label_text3
        lut_pars_frames[key].distType = (
            variables[key]['distribution'])
        lut_pars_frames[key].parIndex = i+1  # save frame index
        lut_pars_frames[key].par_bounds = par_bounds_kwargs

        # Setup widgets
        par_label = ttk.Label(
            lut_pars_frames[key], text=label_text1, width=PARAMETER_LABEL_WIDTH)
        par_menu = tk.OptionMenu(
            lut_pars_frames[key], variables[key]['distribution'],
            *['Fixed Value', 'Discrete'], command=get_lut_distr_fun(key))

        value_label = ttk.Label(lut_pars_frames[key], text='Value:',
                                width=DISTRIBUTION_ARG_LABEL_WIDTH)
        value_field = ParameterEntry(
            lut_pars_frames[key], key, variables[key]['value'],
            DISTRIBUTION_ARG_TEXTFIELD_WIDTH, tool_tip,
            standard_tooltip_text='Select value of {}.'.format(
                lut_pars_frames[key].toolTipText)+tool_tip_end,
            **par_bounds_kwargs)

        # Place widgets on the frame
        par_label.grid(row=0, column=0, sticky='w', padx=5)
        par_menu.grid(row=0, column=1, sticky='ew', padx=5)
        value_label.grid(row=0, column=2, sticky='e', padx=5, pady=0)
        value_field.grid(row=0, column=3, sticky='e', padx=5, pady=0)

        # Configure the widths of menu widgets and entry fields
        par_menu.config(width=DISTRIBUTION_MENU_WIDTH)

        # Assign tooltips for the menu widgets
        tool_tip.bind(par_menu, 'Select distribution for {}.'.format(
            lut_pars_frames[key].distToolTipText))

    for key in ['parSample', 'indexSampleAll']:
        variables[key] = BooleanVar()
        variables[key].set(0)

    # Sampling all data files frame
    sample_all_cb_frame = tk.Frame(lutr_setup_frame)
    sample_all_cb_frame.grid(
        row=0, column=0, columnspan=6, sticky='w')
    index_sample_all_label = ttk.Label(
        sample_all_cb_frame, text="Sample all data file indices:",
        width=PARAMETER_LABEL_WIDTH)
    index_values = [int(val) for val in values['index']]
    index_sample_all_checkbox = tk.Checkbutton(
        sample_all_cb_frame,
        variable=variables['indexSampleAll'],
        command=lambda: sample_all_lut_indices(
            app, lut_pars_frames, variables, index_values))
    index_sample_all_label.grid(
        row=0, column=0, padx=5, sticky='w')
    index_sample_all_checkbox.grid(
        row=0, column=1, padx=5, sticky='e')
    tool_tip.bind(index_sample_all_checkbox,
                  'Check to sample indices of all linked lookup tables.')

    # The checkbox is needed only if parameters other than index parameter are also present
    if len(keys) > 1:
        # Sampling parameter as parameter frame
        sampling_cb_frame = tk.Frame(lutr_setup_frame)
        sampling_cb_frame.grid(
            row=2, column=0, columnspan=6, sticky='w')
        par_sample_label = ttk.Label(
            sampling_cb_frame, text="Sample parameters:",
            width=PARAMETER_LABEL_WIDTH)
        par_sample_checkbox = tk.Checkbutton(
            sampling_cb_frame,
            variable=variables['parSample'],
            command=lambda: sample_lut_parameters(
                variables, lut_pars_frames))
        par_sample_label.grid(
            row=0, column=0, padx=5, sticky='w')
        par_sample_checkbox.grid(
            row=0, column=1, padx=5, sticky='e')
        tool_tip.bind(par_sample_checkbox,
                      ''.join(['Check to sample parameters of lookup table ',
                               'data files instead of indices.']))


    return lut_pars_frames
