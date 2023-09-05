# -*- coding: utf-8 -*-
"""
Module contains methods needed to read data from wellbore type component tabs.
"""
import os
import sys
import ast
import numpy as np

import tkinter as tk
from tkinter import ttk
from tkinter import StringVar, IntVar, DoubleVar
from tkinter import messagebox

# Save location of GUI folder
GUI_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(GUI_DIR)

from dictionarydata import componentVars
from dictionarydata import (PARAMETER_LABEL_WIDTH, MODEL_TAB_LABEL_WIDTH1,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
                            SETUP_ENTRY_WIDTH, BUTTON_WIDTH,
                            FILE_ENTRY_WIDTH)

LABEL_WIDTH = 17
LABEL_WIDTH_THEIS = 37
THEIS_DEFAULTS = {'injTimes': 0, 'injRates': 0.1}


def read_locations_data(data, cmpnt_nm):
    """ Read data located in the well/source location frames."""
    # Get total number of locations
    num_locs = componentVars[cmpnt_nm]['number'].get()
    data['number'] = num_locs

    # Check if the known well locations are provided
    if componentVars[cmpnt_nm]['xCoordinates'].get() and (
            componentVars[cmpnt_nm]['yCoordinates'].get()):
        data['Locations'] = {}
        data['Locations']['coordx'] = []
        data['Locations']['coordy'] = []
        for numx in componentVars[cmpnt_nm]['xCoordinates'].get().split(','):
            data['Locations']['coordx'].append(float(numx.strip()))

        for numy in componentVars[cmpnt_nm]['yCoordinates'].get().split(','):
            data['Locations']['coordy'].append(float(numy.strip()))

    elif componentVars[cmpnt_nm]['xCoordinates'].get():
        messagebox.showerror('Error', ''.join([
            'The y-coordinates of the locations for component {} ',
            'are not provided.']).format(cmpnt_nm))
        return

    elif componentVars[cmpnt_nm]['yCoordinates'].get():
        messagebox.showerror('Error', ''.join([
            'The x-coordinates of the locations for component {} ',
            'are not provided.']).format(cmpnt_nm))
        return

    if len(data['Locations']['coordx']) != len(data['Locations']['coordy']):
        messagebox.showerror('Error', ''.join([
            'Lengths of the provided lists of x-coordinates and y-coordinates ',
            'are not the same.']))
        return

    if num_locs < len(data['Locations']['coordx']):
        messagebox.showwarning('Warning', ''.join([
            'Number of provided locations (through x- and y-coordinates) exceeds ',
            'specified value in parameter Number. Only the first {} locations ',
            'will be used for simulation.']).format(num_locs))
        data['Locations']['coordx'] = data['Locations']['coordx'][0:num_locs]
        data['Locations']['coordy'] = data['Locations']['coordy'][0:num_locs]
        return

    if num_locs > len(data['Locations']['coordx']):
        # Check whether random locations domain is involved
        if componentVars[cmpnt_nm]['useRandomLocDomain'].get():
            data['RandomLocDomain'] = {}
            for key in ['xmin', 'ymin', 'xmax', 'ymax', 'seed']:
                data['RandomLocDomain'][key] = (
                    componentVars[cmpnt_nm]['RandomLocDomain'][key].get())
        else:
            messagebox.showwarning('Warning', ''.join([
                'Number of provided locations (through x- and y-coordinates) ',
                'is less than the specified value in parameter Number. Only provided ',
                'locations will be used. Consider using random locations domain ',
                'to generate additional locations.']))
            data['number'] = len(data['Locations']['coordx'])


def load_locations_data(controller, comp_data, cmpnt_nm):
    componentVars[cmpnt_nm]['number'].set(comp_data['number'])
    # if known locations were provided
    if 'Locations' in comp_data:
        # TODO update values with what is inside
        xcoords = ", ".join(
            str(item) for item in comp_data['Locations']['coordx'])
        ycoords = ", ".join(
            str(item) for item in comp_data['Locations']['coordy'])
        componentVars[cmpnt_nm]['xCoordinates'].set(xcoords)
        componentVars[cmpnt_nm]['yCoordinates'].set(ycoords)
    else:
        componentVars[cmpnt_nm]['xCoordinates'].set('')
        componentVars[cmpnt_nm]['yCoordinates'].set('')

    if 'RandomLocDomain' in comp_data:
        for key_arg in ['xmin', 'ymin', 'xmax', 'ymax', 'seed']:
            componentVars[cmpnt_nm]['RandomLocDomain'][key_arg].set(
                str(comp_data['RandomLocDomain'][key_arg]))
            componentVars[cmpnt_nm]['useRandomLocDomain'].set(1)
    else:
        componentVars[cmpnt_nm]['useRandomLocDomain'].set(0)

    # Get name of frame containing random locations
    parent_frame = controller.nametowidget(controller.getvar(cmpnt_nm+'_locs'))
    randloc_frame = parent_frame.rand_locs_subframe
    disable_rand_locs_frame_widgets(componentVars[cmpnt_nm], randloc_frame)

    if comp_data['type'] in ['OpenWellbore', 'GeneralizedFlowRate']:
        componentVars[cmpnt_nm]['LeakTo'].set(comp_data['LeakTo'])


def add_wellbore_frame_widgets(controller, cmpnt_nm, frame,
                               tool_tip, gfr_cmpnt=False):
    """ Add frame allowing to setup locations for wellbore-type of components."""
    if not gfr_cmpnt:
        tip_pieces = ['Wellbore', 'well', 'wells']
    else:
        tip_pieces = ['Source', 'source', 'sources']

    well_loc_label = ttk.Label(
        frame, text="{} locations:".format(tip_pieces[0]),
        width=PARAMETER_LABEL_WIDTH)
    well_loc_label.grid(row=0, column=0, columnspan=2, sticky='w', padx=5)

    known_well_frame = ttk.Frame(frame)
    known_well_frame.grid(row=1, column=0, sticky='w', padx=20)

    componentVars[cmpnt_nm]['xCoordinates'] = StringVar()
    componentVars[cmpnt_nm]['xCoordinates'].set("100")
    x_coords_label = ttk.Label(known_well_frame, text="x-coordinates [m]:",
                               width=LABEL_WIDTH)
    x_coords_field = tk.Entry(
        known_well_frame, width=2*DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
        textvariable=componentVars[cmpnt_nm]['xCoordinates'])
    x_coords_label.grid(row=0, column=0, pady=5, padx=5, sticky='w')
    x_coords_field.grid(row=0, column=1, pady=5, padx=10, sticky='w')
    tool_tip.bind(x_coords_field, ''.join([
        'Enter x-coordinates of known {} locations, ',
        'separated by comma.']).format(tip_pieces[1]))

    componentVars[cmpnt_nm]['yCoordinates'] = StringVar()
    componentVars[cmpnt_nm]['yCoordinates'].set("100")
    y_coords_label = ttk.Label(known_well_frame, text="y-coordinates [m]:",
                               width=LABEL_WIDTH)
    y_coords_field = tk.Entry(
        known_well_frame, width=2*DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
        textvariable=componentVars[cmpnt_nm]['yCoordinates'])
    y_coords_label.grid(row=1, column=0, pady=5, padx=5, sticky='w')
    y_coords_field.grid(row=1, column=1,
                        pady=5, padx=10, sticky='w')
    tool_tip.bind(y_coords_field, ''.join([
        'Enter y-coordinates of known {} locations, ',
        'separated by comma.']).format(tip_pieces[1]))

    rand_well_label_frame = ttk.Frame(frame)
    rand_well_label_frame.grid(row=2, column=0, sticky='w', padx=0)

    rand_well_domain_label = ttk.Label(
        rand_well_label_frame,
        text="Use random {} domain:".format(tip_pieces[2]),
        width=PARAMETER_LABEL_WIDTH)
    rand_well_domain_label.grid(row=0, column=0, sticky='w', padx=5)

    frame.rand_locs_subframe = ttk.Frame(frame)
    frame.rand_locs_subframe.grid(row=3, column=0, sticky='w', padx=20)

    domain_checkbox = tk.Checkbutton(
        rand_well_label_frame,
        variable=componentVars[cmpnt_nm]['useRandomLocDomain'],
        command=lambda: disable_rand_locs_frame_widgets(
            componentVars[cmpnt_nm], frame.rand_locs_subframe))
    domain_checkbox.grid(row=0, column=1, sticky='w', padx=5)
    tool_tip.bind(domain_checkbox,
                  ''.join(['Check to generate additional random {} ',
                           'locations over the domain.']).format(tip_pieces[1]))

    # produces a bounding box for all random wellbores to be placed within
    componentVars[cmpnt_nm]['RandomLocDomain'] = {}
    setup_dict = {'xmin': 40, 'ymin': 50, 'xmax': 100, 'ymax': 140, 'seed': 345}
    for key in setup_dict:
        if key != 'seed':
            componentVars[cmpnt_nm]['RandomLocDomain'][key] = DoubleVar()
        else:
            componentVars[cmpnt_nm]['RandomLocDomain'][key] = IntVar()
        componentVars[cmpnt_nm]['RandomLocDomain'][key].set(setup_dict[key])

    seed_label = ttk.Label(
        frame.rand_locs_subframe, text="Seed:", width=LABEL_WIDTH)
    seed_txtField = tk.Entry(
        frame.rand_locs_subframe, width=DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
        textvariable=componentVars[cmpnt_nm]['RandomLocDomain']['seed'])

    randomWellXMin_label = ttk.Label(
        frame.rand_locs_subframe, text="x-minimum [m]:", width=LABEL_WIDTH)
    randomWellXMin_txtField = tk.Entry(
        frame.rand_locs_subframe, width=DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
        textvariable=componentVars[cmpnt_nm]['RandomLocDomain']['xmin'])
    randomWellYMin_label = ttk.Label(
        frame.rand_locs_subframe, text="y-minimum [m]:", width=LABEL_WIDTH)
    randomWellYMin_txtField = tk.Entry(
        frame.rand_locs_subframe, width=DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
        textvariable=componentVars[cmpnt_nm]['RandomLocDomain']['ymin'])

    randomWellXMax_label = ttk.Label(
        frame.rand_locs_subframe, text="x-maximum [m]:", width=LABEL_WIDTH)
    randomWellXMax_txtField = tk.Entry(
        frame.rand_locs_subframe, width=DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
        textvariable=componentVars[cmpnt_nm]['RandomLocDomain']['xmax'])
    randomWellYMax_label = ttk.Label(
        frame.rand_locs_subframe, text="y-maximum [m]:", width=LABEL_WIDTH)
    randomWellYMax_txtField = tk.Entry(
        frame.rand_locs_subframe, width=DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
        textvariable=componentVars[cmpnt_nm]['RandomLocDomain']['ymax'])

    seed_label.grid(row=0, column=0, pady=5, padx=5, sticky='w')
    seed_txtField.grid(row=0, column=1, pady=5, padx=5, sticky='w')
    tool_tip.bind(seed_txtField,
                  ''.join(['Enter seed (starting point) for generation ',
                           'of random {} locations.']).format(tip_pieces[1]))

    randomWellXMin_label.grid(row=1, column=0, pady=5, padx=5, sticky='w')
    randomWellXMin_txtField.grid(row=1, column=1, pady=5, padx=5, sticky='w')
    tool_tip.bind(randomWellXMin_txtField,
                  'Enter minimum x-value for random {} domain.'.format(
                      tip_pieces[2]))

    randomWellYMin_label.grid(row=1, column=2, pady=5, padx=5, sticky='w')
    randomWellYMin_txtField.grid(row=1, column=3, pady=5, padx=5, sticky='w')
    tool_tip.bind(randomWellYMin_txtField,
                  'Enter minimum y-value for random {} domain.'.format(
                      tip_pieces[2]))

    randomWellXMax_label.grid(row=2, column=0, pady=5, padx=5, sticky='w')
    randomWellXMax_txtField.grid(row=2, column=1, pady=5, padx=5, sticky='w')
    tool_tip.bind(randomWellXMax_txtField,
                  'Enter maximum x-value for random {} domain.'.format(
                      tip_pieces[2]))

    randomWellYMax_label.grid(row=2, column=2, pady=5, padx=5, sticky='w')
    randomWellYMax_txtField.grid(row=2, column=3, pady=5, padx=5, sticky='w')
    tool_tip.bind(randomWellYMax_txtField,
                  'Enter maximum y-value for random {} domain.'.format(
                      tip_pieces[2]))
    disable_rand_locs_frame_widgets(
            componentVars[cmpnt_nm], frame.rand_locs_subframe)

    # Add reference to the locations frame to the controller to use
    # when the saved simulation is loaded
    controller.setvar(name=cmpnt_nm+'_locs', value=frame)


def disable_rand_locs_frame_widgets(variable, frame):
    """ Disable/enable widgets located on random sources/wells frame."""
    par_frames_state = {0: 'normal', 1: 'disabled'}

    check_button_state = variable['useRandomLocDomain'].get()

    for widget in frame.winfo_children():
        widget.configure(state=par_frames_state[1-check_button_state])


def add_obs_locs_frame_widgets(controller, cmpnt_nm, frame, tool_tip):
    obs_loc_label = ttk.Label(
        frame, text="Observation locations:", width=PARAMETER_LABEL_WIDTH)
    obs_loc_label.grid(row=0, column=0, columnspan=2, sticky='w', padx=5)

    coords_frame = ttk.Frame(frame)
    coords_frame.grid(row=1, column=0, sticky='w', padx=20)

    coords = ['x', 'y']
    arg_names = ['xCoordinates', 'yCoordinates']
    arg_labels = ['x-coordinates [m]:', 'y-coordinates [m]:']

    coords_labels = []
    coords_fields = []
    # Create and place label and entry widgets
    for ind in range(2):
        # Create variable to keep value of the coordinate
        componentVars[cmpnt_nm][arg_names[ind]] = StringVar()
        componentVars[cmpnt_nm][arg_names[ind]].set("")

        coords_labels.append(ttk.Label(
            coords_frame, text=arg_labels[ind], width=LABEL_WIDTH))

        coords_fields.append(tk.Entry(
            coords_frame,
            width=2*DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            textvariable=componentVars[cmpnt_nm][arg_names[ind]]))
        coords_labels[-1].grid(row=ind, column=0, pady=5, padx=5, sticky='w')
        coords_fields[-1].grid(row=ind, column=1, pady=5, padx=10, sticky='w')
        tool_tip.bind(coords_fields[-1], ''.join([
            'Enter {}-coordinates of observation locations, ',
            'separated by comma.\nIf the field is left empty no additional ',
            'observations will be generated.']).format(coords[ind]))


def add_inj_well_frame_widgets(controller, cmpnt_nm, cmpnt_type, frame, tool_tip):

    if cmpnt_type in ['SimpleReservoir', 'AnalyticalReservoir', 'GenericReservoir']:
       inj_well_label_text = "Injection well location:"
       arg_labels = ['x-coordinate [m]:', 'y-coordinate [m]:']
       tool_tip_text =  'Enter {}-coordinate of injection well.'

    elif cmpnt_type =='TheisReservoir':
       inj_well_label_text = "Injection well location(s):"
       arg_labels = ['x-coordinate(s) [m]:', 'y-coordinate(s) [m]:']
       tool_tip_text = ''.join([
          'Enter {}-coordinate of each injection well, separated by commas.\n',
          'The number of the provided x- and y-coordinates must be the same.'])

    inj_well_label = ttk.Label(frame, text=inj_well_label_text)
    inj_well_label.grid(row=0, column=0, columnspan=2, sticky='w', padx=5)

    inj_well_coords_frame = ttk.Frame(frame)
    inj_well_coords_frame.grid(row=1, column=0, sticky='w', padx=20)

    coords = ['x', 'y']
    arg_names = ['injX', 'injY']

    coords_labels = []
    coords_fields = []
    # Create and place label and entry widgets
    for ind in range(2):
        # Create variable to keep value of the coordinate
        componentVars[cmpnt_nm][arg_names[ind]] = StringVar()
        componentVars[cmpnt_nm][arg_names[ind]].set("0")

        coords_labels.append(ttk.Label(
            inj_well_coords_frame, text=arg_labels[ind], width=LABEL_WIDTH))

        coords_fields.append(tk.Entry(
            inj_well_coords_frame,
            width=2*DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            textvariable=componentVars[cmpnt_nm][arg_names[ind]]))
        coords_labels[-1].grid(row=ind, column=0, pady=5, padx=5, sticky='w')
        coords_fields[-1].grid(row=ind, column=1, pady=5, padx=10, sticky='w')
        tool_tip.bind(coords_fields[-1],
                      tool_tip_text.format(coords[ind]))


def add_theis_inj_times_rates_frame_widgets(controller, cmpnt_nm, frame, tool_tip):
    inj_well_label = ttk.Label(frame, text="Injection well rates over time:")
    inj_well_label.grid(row=0, column=0, columnspan=2, sticky='w', padx=5)

    inj_times_rates_subframe = ttk.Frame(frame)
    inj_times_rates_subframe.grid(row=1, column=0, sticky='w', padx=20)

    arg_names = ['injTimes', 'injRates']
    arg_labels = ['Time(s) of injection rate change(s) [years]:',
                  'Injection rate(s) [m^3/s]:']
    # The rows in the file with the injection times correspond to different
    # injection wells in the order their x- and y-coordinates are specified.
    # The rows in the file with the injection rates correspond to different
    # injection wells in the order their x- and y-coordinates are specified;
    # the columns correspond to the time points of injection times.
    entry_tip_text = [
        ''.join(['Enter time points at which injection rates change manually\n',
                 'or provide path to the file containing the input.']),
        ''.join(['Enter injection rates manually (positive for injection, ',
                 'negative for extraction)\n',
                 'or provide path to the file containing the input.'])]
    button_tip_text = ['Select file containing time points.',
                       'Select file containing injection rates.']
    dialog_titles = ['Choose file containing time points',
                     'Choose file containing injection rates']

    inputs_labels = []
    # Create and place label and entry widgets
    for ind in range(len(arg_names)):
        # Create variable to keep value of the coordinate
        componentVars[cmpnt_nm][arg_names[ind]] = StringVar()
        componentVars[cmpnt_nm][arg_names[ind]].set("")

        inputs_labels.append(ttk.Label(
            inj_times_rates_subframe, text=arg_labels[ind], width=LABEL_WIDTH_THEIS))
        inputs_labels[-1].grid(row=ind, column=0, pady=5, padx=5, sticky='w')

        add_file_input_widgets(controller, inj_times_rates_subframe,
                               tool_tip, componentVars[cmpnt_nm][arg_names[ind]],
                               entry_tip_text[ind], button_tip_text[ind],
                               dialog_titles[ind], row_ind=ind, col_ind=1,
                               entry_width=FILE_ENTRY_WIDTH)

def read_obs_locations_data(data, cmpnt_nm):
    """ Read data located in the observation location and injection well
    (if applicable) frames."""
    # Check if the observation locations are provided
    if componentVars[cmpnt_nm]['xCoordinates'].get():
        if componentVars[cmpnt_nm]['yCoordinates'].get():

            data['Locations'] = {}
            data['Locations']['coordx'] = []
            data['Locations']['coordy'] = []
            for numx in componentVars[cmpnt_nm]['xCoordinates'].get().split(','):
                data['Locations']['coordx'].append(float(numx.strip()))

            for numy in componentVars[cmpnt_nm]['yCoordinates'].get().split(','):
                data['Locations']['coordy'].append(float(numy.strip()))

            if len(data['Locations']['coordx']) != len(data['Locations']['coordy']):
                messagebox.showerror('Error', ''.join([
                    'Lengths of the provided lists of x-coordinates and y-coordinates ',
                    'are not the same.']))
                return
        else:
            messagebox.showerror('Error', ''.join([
                'The y-coordinates of the locations for component {} ',
                'are not provided.']).format(cmpnt_nm))
            return
    else:
        if componentVars[cmpnt_nm]['yCoordinates'].get():
            messagebox.showerror('Error', ''.join([
            'The x-coordinates of the locations for component {} ',
            'are not provided.']).format(cmpnt_nm))
            return

    # Check the type of the reservoir component
    if data['type'] in ['SimpleReservoir', 'AnalyticalReservoir', 'GenericReservoir']:
        if componentVars[cmpnt_nm]['injX'].get():
            if componentVars[cmpnt_nm]['injY'].get():
                data['InjectionWell'] = {}
                data['InjectionWell']['coordx'] = float(
                    componentVars[cmpnt_nm]['injX'].get())
                data['InjectionWell']['coordy'] = float(
                    componentVars[cmpnt_nm]['injY'].get())
            else:
                messagebox.showerror('Error', ''.join([
                    'The y-coordinate of the injection well location for component {} ',
                    'is not provided.']).format(cmpnt_nm))
                return
        else:
            if componentVars[cmpnt_nm]['injY'].get():
                messagebox.showerror('Error', ''.join([
                    'The x-coordinate of the injection well location for component {} ',
                    'is not provided.']).format(cmpnt_nm))
                return
            else:
                messagebox.showinfo('Information', ''.join([
                    "Default location of injection well at (0, 0) ",
                    "for component {} will be used."]).format(cmpnt_nm))
                data['InjectionWell'] = {}
                componentVars[cmpnt_nm]['injX'].set("0")
                componentVars[cmpnt_nm]['injY'].set("0")
                data['InjectionWell']['coordx'] = 0.0
                data['InjectionWell']['coordy'] = 0.0

    # Check the type of the reservoir component
    if data['type'] == 'TheisReservoir':
        if componentVars[cmpnt_nm]['injX'].get():
            if componentVars[cmpnt_nm]['injY'].get():
                data['InjectionWell'] = {}
                try:
                    data['InjectionWell']['coordx'] = float(
                        componentVars[cmpnt_nm]['injX'].get())
                    data['InjectionWell']['coordy'] = float(
                        componentVars[cmpnt_nm]['injY'].get())
                except ValueError:
                    data['InjectionWell']['coordx'] = [float(val) for val in
                        componentVars[cmpnt_nm]['injX'].get().split(',')]
                    data['InjectionWell']['coordy'] = [float(val) for val in
                        componentVars[cmpnt_nm]['injY'].get().split(',')]

            else:
                messagebox.showerror('Error', ''.join([
                    'The y-coordinate(s) of the injection well location(s) for component {} ',
                    'is/are not provided.']).format(cmpnt_nm))
                return
        else:
            if componentVars[cmpnt_nm]['injY'].get():
                messagebox.showerror('Error', ''.join([
                    'The x-coordinate(s) of the injection well location(s) for component {} ',
                    'is/are not provided.']).format(cmpnt_nm))
                return
            else:
                messagebox.showinfo('Information', ''.join([
                    "Default location of injection well at (0, 0) ",
                    "for component {} will be used."]).format(cmpnt_nm))
                data['InjectionWell'] = {}
                componentVars[cmpnt_nm]['injX'].set("0")
                componentVars[cmpnt_nm]['injY'].set("0")
                data['InjectionWell']['coordx'] = 0.0
                data['InjectionWell']['coordy'] = 0.0

        if componentVars[cmpnt_nm]['injRates'].get():
            if componentVars[cmpnt_nm]['injTimes'].get():
                input_data = {'injTimes': componentVars[cmpnt_nm]['injTimes'].get(),
                              'injRates': componentVars[cmpnt_nm]['injRates'].get()}
                for input_key in ['injTimes', 'injRates']:
                    if ',' in input_data[input_key]: # list of values
                        data[input_key] = [float(val) for val in input_data[input_key].split(',')]
                    else:
                        try:
                            data[input_key] = float(input_data[input_key])
                        except ValueError:
                            data[input_key] = input_data[input_key]
            else:
                messagebox.showerror('Error', ''.join([
                    'The injTimes for component {} are not provided.']).format(cmpnt_nm))
                return
        else:
            if componentVars[cmpnt_nm]['injTimes'].get():
                messagebox.showerror('Error', ''.join([
                    'The injRates for component {} are not provided.']).format(cmpnt_nm))
                return
            else:
                messagebox.showinfo('Information', ''.join([
                    "Default injRates of 0.1 m^3/s ",
                    "for component {} will be used."]).format(cmpnt_nm))
                data['InjectionWell'] = {}
                for input_key in ['injTimes', 'injRates']:
                    componentVars[cmpnt_nm][input_key].set(str(THEIS_DEFAULTS[input_key]))
                    data[input_key] = THEIS_DEFAULTS[input_key]

    return data


def load_theis_inj_times_rates_data(comp_data, cmpnt_nm):
    """
    Load input data for Theis Reservoir component related to injection times and
    injection rates.
    """
    for input_key in ['injTimes', 'injRates']:
        if input_key in comp_data:
            input_data = comp_data[input_key]
            if isinstance(input_data, list):
                componentVars[cmpnt_nm][input_key].set(
                    ', '.join([str(val) for val in input_data]))
            elif isinstance(input_data, (float, int)):
                componentVars[cmpnt_nm][input_key].set(str(input_data))
            elif isinstance(input_data, str):  # possibly path to the file
                componentVars[cmpnt_nm][input_key].set(input_data)
        else:
            componentVars[cmpnt_nm][input_key].set(str(THEIS_DEFAULTS[input_key]))


def load_obs_locations_data(comp_data, cmpnt_nm):
    # If locations are provided
    if 'Locations' in comp_data:
        xcoords = ", ".join(
            str(item) for item in comp_data['Locations']['coordx'])
        ycoords = ", ".join(
            str(item) for item in comp_data['Locations']['coordy'])
        componentVars[cmpnt_nm]['xCoordinates'].set(xcoords)
        componentVars[cmpnt_nm]['yCoordinates'].set(ycoords)
    else:
        componentVars[cmpnt_nm]['xCoordinates'].set('')
        componentVars[cmpnt_nm]['yCoordinates'].set('')

    if comp_data['type'] in ['SimpleReservoir', 'AnalyticalReservoir', 'GenericReservoir']:
        if 'InjectionWell' in comp_data:
            componentVars[cmpnt_nm]['injX'].set(str(
                comp_data['InjectionWell']['coordx']))
            componentVars[cmpnt_nm]['injY'].set(str(
                comp_data['InjectionWell']['coordy']))
        else:
            componentVars[cmpnt_nm]['injX'].set('0')
            componentVars[cmpnt_nm]['injY'].set('0')
    elif comp_data['type'] == 'TheisReservoir':
        if 'InjectionWell' in comp_data:
            coordx_data = comp_data['InjectionWell']['coordx']
            coordy_data = comp_data['InjectionWell']['coordy']
            if isinstance(coordx_data, list) and isinstance(coordy_data, list):
                componentVars[cmpnt_nm]['injX'].set(', '.join([
                    str(val) for val in coordx_data]))
                componentVars[cmpnt_nm]['injY'].set(', '.join([
                    str(val) for val in coordy_data]))
            elif isinstance(coordx_data, (float, int)) and \
                    isinstance(coordy_data, (float, int)):
                componentVars[cmpnt_nm]['injX'].set(str(
                    comp_data['InjectionWell']['coordx']))
                componentVars[cmpnt_nm]['injY'].set(str(
                    comp_data['InjectionWell']['coordy']))
            else:
                messagebox.showinfo('Information', ''.join([
                    "x- and y-coordinates of injection well(s) for component {}",
                    "are of incompatible types."]).format(cmpnt_nm))
        else:
            componentVars[cmpnt_nm]['injX'].set('0')
            componentVars[cmpnt_nm]['injY'].set('0')


def add_cell_locs_frame_widgets(controller, cmpnt_nm, frame, variables, tool_tip):
    obs_loc_label = ttk.Label(
        frame, text="Cell locations:", width=PARAMETER_LABEL_WIDTH)
    obs_loc_label.grid(row=0, column=0, columnspan=2, sticky='w', padx=5)

    frame.coords_frame = ttk.Frame(frame)
    frame.coords_frame.grid(row=1, column=0, sticky='w', padx=20)

    coords = ['x', 'y']
    arg_names = ['coordx', 'coordy']
    arg_labels = ['x-coordinates [m]:', 'y-coordinates [m]:']

    coords_labels = []
    coords_fields = []

    # Create and place label and entry widgets
    for ind in range(2):
        coords_labels.append(ttk.Label(
            frame.coords_frame, text=arg_labels[ind], width=LABEL_WIDTH))

        coords_fields.append(tk.Entry(
            frame.coords_frame,
            width=2*DISTRIBUTION_ARG_TEXTFIELD_WIDTH,
            textvariable=variables[arg_names[ind]]))
        coords_labels[-1].grid(row=ind, column=0, pady=5, padx=5, sticky='w')
        coords_fields[-1].grid(row=ind, column=1, pady=5, padx=10, sticky='w')
        tool_tip.bind(
            coords_fields[-1],
            'Enter {}-coordinates of cell locations separated by comma.'.format(
                coords[ind]))

    frame.coords_frame.coords_labels = coords_labels
    frame.coords_frame.coords_fields = coords_fields

    file_input_label = ttk.Label(
        frame.coords_frame, text="Use file input for cell centers:",
        width=PARAMETER_LABEL_WIDTH)
    file_input_checkbox = tk.Checkbutton(
        frame.coords_frame, variable=variables['fileinput'],
        command=lambda: disable_cell_locs_frame_widgets(frame.coords_frame))
    file_input_label.grid(
        row=2, column=0, pady=5, padx=5, sticky='w')
    file_input_checkbox.grid(
        row=2, column=1, pady=5, padx=5, sticky='w')
    tool_tip.bind(file_input_checkbox,
                  'Check to use file input for cell centers.')

    frame.coords_frame.checkbox_variable = variables['fileinput']

    filename_label = ttk.Label(
        frame.coords_frame, text="Cell coordinates file:", width=PARAMETER_LABEL_WIDTH)
    filename_label.grid(
        row=3, column=0, pady=5, padx=5, sticky='w')

    add_file_input_widgets(
        controller, frame.coords_frame, tool_tip,
        componentVars[cmpnt_nm]['Cells']['Locations']['filename'],
        'Provide path to the file containing cell coordinates data.',
        'Select file containing cell coordinates data.',
        'Choose file containing cell coodinates data',
        row_ind=3, col_ind=1, entry_width=FILE_ENTRY_WIDTH)

    # Disable entry and browse button
    frame.coords_frame.filename_entry.configure(state='disabled')
    frame.coords_frame.browse_button.configure(state='disabled')


def disable_cell_locs_frame_widgets(frame):
    frame_state = {0: 'normal', 1: 'disabled'}

    check_button_state = frame.checkbox_variable.get()

    for ind in range(2):
        frame.coords_labels[ind].configure(state=frame_state[check_button_state])
        frame.coords_fields[ind].configure(state=frame_state[check_button_state])

    frame.filename_entry.configure(state=frame_state[1-check_button_state])
    frame.browse_button.configure(state=frame_state[1-check_button_state])


def add_file_input_widgets(controller, frame, tool_tip, variable,
                           entry_tooltip, button_tooltip, dialog_title,
                           row_ind=0, col_ind=2, entry_width=SETUP_ENTRY_WIDTH):
    """ Add widgets related to the parameters setup through file input."""

    filename_entry = tk.Entry(frame, width=entry_width, textvariable=variable)

    browse_button = ttk.Button(
        frame, text="Browse",
        command=lambda: controller.choose_file(variable, dialog_title),
        width=BUTTON_WIDTH)

    # Place widgets
    filename_entry.grid(row=row_ind, column=col_ind, padx=5)
    browse_button.grid(row=row_ind, column=col_ind+1, padx=5, sticky='w')
    tool_tip.bind(filename_entry, entry_tooltip)
    tool_tip.bind(browse_button, button_tooltip)

    frame.filename_entry = filename_entry
    frame.browse_button = browse_button
