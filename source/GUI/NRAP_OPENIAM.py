# Pmw is imported to enable hover text
# pickle is used to generate and read binary files
import os
import sys
import pickle
import random
import logging
import warnings
import yaml
logging.basicConfig(level=logging.WARNING)

import tkinter as tk
from tkinter import ttk, StringVar, IntVar, DoubleVar, BooleanVar, messagebox

import matplotlib
debug_msg = 'Available matplotlib backends: {}'.format(matplotlib.rcsetup.all_backends)
logging.debug(debug_msg)

with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    # matplotlib.use('TkAgg')

import numpy as np
import Pmw

from Disclaimer import Disclaimer_Page
from Dashboard import Dashboard_Page
from OpenIAM_Page import OpenIAM_Page, disable_time_frame_widgets
from PostProcessor_Page import PostProcessor_Page

from dictionarydata import (d, APP_SIZE, TAB_SIZE, componentVars, componentChoices,
                            componentTypeDictionary, connectionsDictionary,
                            DISTRIBUTION_OPTIONS, connections, connectionTypes,
                            COMPONENT_TYPES, ANALYSIS_TYPES,
                            DISTRIBUTION_MENU_WIDTH, DISTRIBUTION_ARG_LABEL_WIDTH,
                            DISTRIBUTION_ARG_TEXTFIELD_WIDTH, PARAMETER_LABEL_WIDTH,
                            GFR_PARAMETER_LABEL_WIDTH, STRATA_PARAMETER_LABEL_WIDTH,
                            FL_PARAMETER_LABEL_WIDTH,
                            MODEL_TAB_LABEL_WIDTH2, MODEL_TAB_LABEL_WIDTH3,
                            MODEL_TAB_ENTRY_WIDTH, MODEL_TAB_MENU_WIDTH,
                            DISTRIBUTION_PARS_LABELS, DISTRIBUTION_PARS_SETUPS)

from cmpnts_tabs import (src_tab, arc_tab, grc_tab, trc_tab, msw_tab, lutr_tab, cw_tab,
                         cwwr_tab, ow_tab, gfr_tab, ff_tab, fl_tab, hcl_tab, sh_tab, ca_tab,
                         aalf_tab, daa_tab, daaml_tab, fgaq_tab, fgaz_tab,
                         ga_tab, atm_tab, psa_tab, strata_tab, cws_tab, locations)
from cmpnts_tabs.parameter_entry import ParameterEntry


# Save location of source folder in the top level folder
SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(SOURCE_DIR)
CODE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
USER_DIR = os.sep.join([CODE_DIR, 'examples', 'user'])


class NRAPOpenIAM(tk.Tk):
    """
    Class for NRAP-Open-IAM GUI controller. """
    def __init__(self, *args, **kwargs):
        """ Constructor method for NRAPOpenIAM class. """
        tk.Tk.__init__(self, *args, **kwargs)

        container = tk.Frame(self)
        container.grid(row=0, column=0, sticky='nsew')
        componentVars['simName'] = StringVar()

        # Define folder where user files will be saved
        # If None it means no "save" requests were made
        self.user_dir = None
        self.sim_file = None

        self.frames = {}

        for F in (Disclaimer_Page, Dashboard_Page, OpenIAM_Page, PostProcessor_Page):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, rowspan=20, columnspan=20, sticky="nsew")

        self.show_frame(Disclaimer_Page)

    def populate_dictionary(self, ask_to_save=True):
        """
        Save the current simulation out to a binary file to allow
        later editing with the load_simulation method.
        """
        from tkinter.filedialog import asksaveasfilename

        global d

        d = {}
        try:
            fileName = componentVars['simName'].get()
            d['simName'] = componentVars['simName'].get()
        except:
            fileName = d['simName']

        if self.user_dir is None:
            self.user_dir = USER_DIR

        if ask_to_save:
            dialog_answer = asksaveasfilename(
                initialdir=self.user_dir,
                initialfile=fileName,
                title="Save simulation file",
                filetypes=[("Open IAM Files", "*.OpenIAM")])

            if dialog_answer:
                if dialog_answer[-8:] != '.OpenIAM':
                    fileName = dialog_answer+'.OpenIAM'
                else:
                    fileName = dialog_answer
                self.user_dir = os.path.dirname(os.path.abspath(fileName))
        else:
            fileName = os.sep.join([USER_DIR, fileName+'.OpenIAM'])

        self.sim_file = fileName

        d['ModelParams'] = {}
        if componentVars['outputDirectoryGenerate'].get():
            d['ModelParams']['OutputDirectory'] = componentVars[
                'outputDirectory'].get()+'{datetime}'
        else:
            d['ModelParams']['OutputDirectory'] = componentVars[
                'outputDirectory'].get()

        for i, choice in enumerate(componentChoices):
            read_tab_pars = self.get_read_vars_method(componentTypeDictionary[i])
            d[choice] = read_tab_pars(choice)

            if componentTypeDictionary[i] == 'LookupTableReservoir':
                continue

            if componentTypeDictionary[i] == 'PlumeStability':
                continue

            self.process_parameter_vars(componentVars[choice]['Params'],
                                        d[choice]['Parameters'],
                                        componentVars[choice]['Params'].keys())

            if 'Controls' in componentVars[choice]:
                d[choice]['Controls'] = {}
                self.process_control_vars(componentVars[choice]['Controls'],
                                          d[choice]['Controls'],
                                          componentVars[choice]['Controls'].keys())

            if componentTypeDictionary[i] == 'SealHorizon':
                self.process_parameter_vars(
                    componentVars[choice]['Cells'],
                    d[choice]['Cells'],
                    sh_tab.SH_CELL_PARAMETERS,
                    exclusion=True)
                sh_tab.clean_data(d, choice)
                continue

            if componentTypeDictionary[i] == 'FaultFlow':
                ff_tab.clean_data(d, choice)
                continue

        d['Stratigraphy'] = strata_tab.read_tab_vars(self)

        # Setup analysis type and corresponding parameters
        d['ModelParams']['Analysis'] = {}
        analysis_type = componentVars['analysis']['type'].get()
        d['ModelParams']['Analysis']['type'] = analysis_type

        if analysis_type == 'LHS':
            d['ModelParams']['Analysis']['siz'] = (
                componentVars['analysis']['size'].get())
            d['ModelParams']['Analysis']['seed'] = (
                componentVars['analysis']['seed'].get())

        if analysis_type == 'Parstudy':
            try:
                vals = int(componentVars['analysis']['nvals'].get())
            except ValueError:
                vals = []
                for value in componentVars['analysis']['nvals'].get().split(','):
                    vals.append(int(value.strip()))

            d['ModelParams']['Analysis']['nvals'] = vals

        d = self.process_time_frame(d)
        d['ModelParams']['Logging'] = componentVars['logging'].get()
        d['ModelParams']['Components'] = componentChoices
        # Column or row-wise orientation of resulting arrays in output files:
        # True is column-wise; False is row-wise
        d['ModelParams']['OutputType'] = self.OutputType.get()
        d['ModelParams']['GenerateOutputFiles'] = self.GenerateOutputFiles.get()
        d['ModelParams']['GenerateCombOutputFile'] = self.GenerateCombOutputFile.get()
        d['ModelParams']['GenerateStatFiles'] = self.GenerateStatFiles.get()

        try:
            with open(fileName, 'wb') as outPutFile:
                pickle.dump(d, outPutFile, pickle.HIGHEST_PROTOCOL)

            yaml_filename = fileName[0:-8]+'.yaml'
            with open(yaml_filename, 'w') as yaml_file:
                yaml.dump(d, yaml_file, default_flow_style=False, explicit_start=True)
        except:
            return

    @staticmethod
    def process_time_frame(d):
        if not componentVars['timePointsInput'].get():
            flag_var = d['ModelParams'].pop('TimePoints', 'Not used')
            d['ModelParams']['EndTime'] = componentVars['endTime'].get()
            d['ModelParams']['TimeStep'] = componentVars['timeStep'].get()

        else: # if file or manual input are chosen
            flag_var1 = d['ModelParams'].pop('EndTime', 'Not used')
            flag_var2 = d['ModelParams'].pop('TimeStep', 'Not used')
            inp_data = componentVars['timePoints'].get()

            if ',' in inp_data:
                data = inp_data.split(',')
                d['ModelParams']['TimePoints'] = [float(val.strip()) for val in data]
            else:
                inp_data_file_path = os.path.join(CODE_DIR, inp_data)
                if os.path.isfile(inp_data_file_path):
                    d['ModelParams']['TimePoints'] = inp_data
                else:
                    raise FileNotFoundError(
                        'File {} is not found.'.format(inp_data_file_path))

        return d

    def process_parameter_vars(self, var_dict_from, data_dict_to, var_keys, exclusion=None):
        """ Read information from componentVars kept in var_dict_from to copy to
            data_dict_to.

        Read information from componentVars corresponding to the parameters
        of the component.
        """
        for key in var_keys:
            distr_type = var_dict_from[key][
                'distribution'].get()
            if distr_type == 'Fixed Value':
                if exclusion is None:  # special case of fixed value parameters
                    data_dict_to[key]['value'] = var_dict_from[key]['value'].get()
                    data_dict_to[key]['vary'] = False
                else:
                    data_dict_to[key] = var_dict_from[key]['value'].get()

            if distr_type == 'List':
                values = []
                for value in var_dict_from[key]['ordvalues'].get().split(','):
                    values.append(float(value.strip()))
                data_dict_to[key] = values

            if distr_type == 'File Input':
                data_dict_to[key] = var_dict_from[key]['filename'].get()

            if distr_type == 'Uniform':
                data_dict_to[key]['dist'] = 'uniform'
                data_dict_to[key]['min'] = (
                    var_dict_from[key]['min'].get())
                data_dict_to[key]['max'] = (
                    var_dict_from[key]['max'].get())

            if distr_type == 'Normal':
                data_dict_to[key]['dist'] = 'norm'
                data_dict_to[key]['mean'] = (
                    var_dict_from[key]['mean'].get())
                data_dict_to[key]['std'] = (
                    var_dict_from[key]['std'].get())

            if distr_type == 'Lognormal':
                data_dict_to[key]['dist'] = 'lognorm'
                kwargs = {'mean': var_dict_from[key]['mean'].get(),
                          'std': var_dict_from[key]['std'].get()}
                dist_pars = self.reparametrize_lognorm_distribution(kwargs)
                data_dict_to[key]['dist_pars'] = [
                    dist_pars['s'], dist_pars['loc'], dist_pars['scale']]

            # TODO As of now scipy.stats does not sample effectively when
            # the values are at the tails.
            # See issue here: https://github.com/scipy/scipy/issues/10092
            # We leave it as it is for now, i.e. truncated
            # distribution is suggested to be replaced with other distributions
            if distr_type == 'Truncated':
                data_dict_to[key]['dist'] = 'truncnorm'
                kwargs = {
                    'min': var_dict_from[key]['min'].get(),
                    'max': var_dict_from[key]['max'].get(),
                    'mean': var_dict_from[key]['mean'].get(),
                    'std': var_dict_from[key]['std'].get()}
                dist_pars = self.reparametrize_truncated_distribution(kwargs)
                data_dict_to[key]['dist_pars'] = [
                    dist_pars['a'], dist_pars['b'],
                    dist_pars['loc'], dist_pars['scale']]

            if distr_type == 'Triangular':
                data_dict_to[key]['dist'] = 'triang'
                kwargs = {
                    'min': var_dict_from[key]['min'].get(),
                    'max': var_dict_from[key]['max'].get(),
                    'mode': var_dict_from[key]['mode'].get()}
                dist_pars = self.reparametrize_triang_distribution(kwargs)
                data_dict_to[key]['dist_pars'] = [
                    dist_pars['c'], dist_pars['loc'], dist_pars['scale']]

            if distr_type == 'Discrete':
                data_dict_to[key]['discrete_vals'] = []
                discrete_values = []
                discrete_weights = []

                for value in var_dict_from[key]['values'].get().split(','):
                    discrete_values.append(float(value.strip()))

                for weight in var_dict_from[key]['weights'].get().split(','):
                    discrete_weights.append(float(weight.strip()))

                data_dict_to[key]['discrete_vals'].append(discrete_values)
                data_dict_to[key]['discrete_vals'].append(discrete_weights)


    def process_control_vars(self, var_dict_from, data_dict_to, var_keys):
        """ Read information from componentVars[component_name]['Controls'] kept
            in var_dict_from to copy to data_dict_to[component_name]['Controls'].

        Read information from componentVars corresponding to the controls
        of the component.
        """
        for key in var_keys:
            data_dict_to[key] = var_dict_from[key].get()


    @staticmethod
    def reformat_list_presentation(val_list):
        """ Reformat list representation for tooltip hints."""

        quotient, remainder = divmod(len(val_list), 10)

        S = '['
        if len(val_list) >= 10:
            for ind in range(quotient):
                S = S + ', '.join([str(val) for val in val_list[ind*10:(ind+1)*10]]) + ',\n'

        if remainder == 0:
            S = S[0:-2] + ']'
        else:
            S = S + ', '.join([str(val) for val in val_list[quotient*10:]]) + ']'

        return S


    def add_remaining_widgets(self, par_name, dist_type, frame, toolTip, vars_dict):
        """ Add widgets related to the parameters setup. """
        # Determine type of distribution of the component parameter
        before_tooltip = ''
        # Modify content of tooltip if distribution is lognormal
        if dist_type == 'Lognormal':
            before_tooltip = 'logarithm of '

        if dist_type == 'File Input':
            locations.add_file_input_widgets(
                self, frame, toolTip, vars_dict['filename'],
                'Provide path to the file containing {} data.'.format(
                    frame.toolTipText),
                'Select file containing {} data.'.format(frame.toolTipText),
                'Choose file containing values for parameter {}'.format(
                    frame.toolTipText))
            return

        # Determine number of parameters
        num_pars = len(DISTRIBUTION_PARS_LABELS[dist_type])
        # Initialize list of Label and Entry widgets
        pars_Labels = []
        pars_Entries = []
        for ind in range(num_pars):
            # Add Label widget for each distribution parameter
            # par_text_label is something like 'Value:', 'Minimum:', 'Maximum:'
            par_text_label = DISTRIBUTION_PARS_LABELS[dist_type][ind]
            pars_Labels.append(ttk.Label(
                frame, text=par_text_label, width=DISTRIBUTION_ARG_LABEL_WIDTH))

            # key is something like 'value', 'min', 'max', etc.
            key = DISTRIBUTION_PARS_SETUPS[par_text_label][0]

            if frame.par_bounds['discrete_bounds'] is not None:
                if key in ['value', 'values']:
                    par_bounds_kwargs = {
                        'discrete_bounds': frame.par_bounds['discrete_bounds'],
                        'to_validate': True}
                    entry_tooltip_text =''.join([
                        DISTRIBUTION_PARS_SETUPS[par_text_label][1],
                        '\nPossible values are \n{1}.'])
                    standard_tooltip_text = entry_tooltip_text.format(
                        before_tooltip+par_name,
                        self.reformat_list_presentation(
                            frame.par_bounds['discrete_bounds']))

                elif key == 'weights':
                    par_bounds_kwargs = {'lower_bound': 0,
                                         'upper_bound': np.inf,
                                         'to_validate': True}
                    entry_tooltip_text =''.join([
                        DISTRIBUTION_PARS_SETUPS[par_text_label][1],
                        '\nPossible values are between 0 and inf.'])
                    standard_tooltip_text = entry_tooltip_text.format(
                        before_tooltip+par_name)

            else:

                if key in ['value', 'min', 'max', 'mode', 'mean', 'values', 'ordvalues']:
                    lower_bound = frame.par_bounds['lower_bound']
                    upper_bound = frame.par_bounds['upper_bound']
                    to_validate = True
                elif key == 'std':
                    lower_bound = 0
                    upper_bound = (
                        frame.par_bounds['upper_bound']-frame.par_bounds['lower_bound'])/2
                    to_validate = True
                else:  # for weights
                    lower_bound = 0
                    upper_bound = np.inf
                    to_validate = True

                entry_tooltip_text =''.join([
                    DISTRIBUTION_PARS_SETUPS[par_text_label][1],
                    '\nPossible values are between {1} and {2}.'])
                standard_tooltip_text = entry_tooltip_text.format(
                        before_tooltip+par_name, lower_bound, upper_bound)
                par_bounds_kwargs = {'lower_bound': lower_bound,
                                     'upper_bound': upper_bound,
                                     'to_validate': to_validate}

            # Add Entry widget for each distribution parameter
            pars_Entries.append(
                ParameterEntry(
                    frame, par_name, vars_dict[key],
                    DISTRIBUTION_ARG_TEXTFIELD_WIDTH, toolTip,
                    standard_tooltip_text=standard_tooltip_text,
                    **par_bounds_kwargs))

            # Place Label and Entry widget onto frame
            pars_Labels[-1].grid(
                row=0, column=2*(ind+1), padx=5)
            pars_Entries[-1].grid(
                row=0, column=2*(ind+1)+1, padx=5)
            pars_Entries[-1].config(width=DISTRIBUTION_ARG_TEXTFIELD_WIDTH)


    def change_distribution(self, frame, pars_label_width=PARAMETER_LABEL_WIDTH):
        """
        Change the distribution for a parameter based on user inputs.
        """
        toolTip = Pmw.Balloon(self)
        for widget in frame.winfo_children():
            toolTip.unbind(widget)
            widget.destroy()

        # Label with parameter name and units
        par_name_label = ttk.Label(frame, width=pars_label_width,
                                   text=frame.labelText)
        # Distribution menu
        distr_menu = tk.OptionMenu(
            frame, frame.distType, *frame.distr_options,
            command=lambda _: self.change_distribution(
                frame, pars_label_width=pars_label_width))

        # Configure widgets
        distr_menu.config(width=DISTRIBUTION_MENU_WIDTH)
        toolTip.bind(distr_menu, 'Select distribution for {}.'.format(
            frame.toolTipText))
        # Place widgets
        par_name_label.grid(row=0, column=0, sticky='w', padx=5)
        distr_menu.grid(row=0, column=1, padx=5)

        self.add_remaining_widgets(frame.toolTipText, frame.distType.get(),
                                   frame, toolTip, frame.paramVars)

    def load_simulation(self, data):
        """
        Load selected simulation into a global variable.

        After loading the simulation the add_component method is called to add
        any existing component models.
        """

        try:
            componentVars['simName'].set(data['simName'])
        except:
            componentVars['simName'] = StringVar()
            componentVars['simName'].set(data['simName'])

        global connectionTypes

        for component in componentChoices:
            del componentVars[component]

            self.tabControl.connection_menu.children['menu'].delete(0, 'end')

            for c in connectionsDictionary:
                self.tabControl.connection_menu.children['menu'].add_command(
                    label=c, command=lambda con=c: \
                        self.tabControl.connection_menu.connection.set(con))

            for tab in self.tabControl.tabs():
                if self.tabControl.tab(tab, option="text") == component:
                    self.tabControl.forget(tab)

        i = len(componentChoices)-1
        while i >= 0:
            componentChoices.pop(i)
            componentTypeDictionary.pop(i)
            connectionsDictionary.pop(i)
            connections.pop(i)
            i = i - 1

        connectionTypes = []
        connections[0] = 'Dynamic Parameters'

        # Setup Model Tab variables
        self.process_model_data(data)

        # Setup stratigraphy tab initial parameters
        componentVars['strata']['Params']['datumPressure'].set(
            data['Stratigraphy']['datumPressure'])
        componentVars['strata']['Params']['numberOfShaleLayers'].set(
            data['Stratigraphy']['numberOfShaleLayers']['value'])

        strata_tab.add_stratigraphy_layers(
            data['Stratigraphy']['numberOfShaleLayers']['value'], self)

        # Process loaded stratigraphy parameters
        self.process_strata_data(data)

        # Get names of components in the saved model setup
        ckeys = data['ModelParams']['Components']

        # Set connection menu command
        # tabControl.connection_menu is defined in OpenIAM_Page.py file
        # Its purpose is to handle connections between components within the same
        # system model
        for c in ckeys:
            self.tabControl.connection_menu.children['menu'].add_command(
                label=c, command=lambda con=c: \
                    self.tabControl.connection_menu.connection.set(con))

        # Process connections and dynamic input
        for key in ckeys:
            connection_name = data[key]['connection']
            try:
                aquiferName = data[key]['AquiferName']
            except KeyError:
                try:
                    aquiferName = data[key]['LeakTo']
                except KeyError:
                    aquiferName = 'none'

            try:
                controls = data[key]['Controls']
            except KeyError:
                controls = {}

            if connection_name == 'Dynamic Parameters':
                if (data[key]['type'].find('Wellbore') != -1) or (
                        data[key]['type'].find('FaultFlow') != -1) or (
                            data[key]['type'].find('Seal') != -1):
                    dp_keys = ['pressure', 'CO2saturation']
                if data[key]['type'].find('Aquifer') != -1:
                    if data[key]['type'].find('Generic') != -1:
                        dp_keys = ['brine_mass', 'co2_mass']
                    else:
                        dp_keys = ['brine_rate', 'co2_rate', 'brine_mass', 'co2_mass']
                if data[key]['type'].find('Atm') != -1:
                    dp_keys = ['co2_leakrate']

                dyn_data = []
                for dp_key in dp_keys:
                    # Check whether list was provided
                    inp_data = data[key]['DynamicParameters'][dp_key]
                    if isinstance(inp_data, list):
                        dyn_data.append(", ".join(
                            str(item) for item in inp_data))
                    # Check whether string was provided
                    elif isinstance(inp_data, str):
                        inp_data_file_path = os.path.join(CODE_DIR, inp_data)
                        if os.path.isfile(inp_data_file_path):
                            dyn_data.append(inp_data)
                        else:
                            msg = ''.join([
                                'Path to the file provided as ',
                                'a source of dynamic data is not valid.'])
                            raise Warning(msg)
            else:
                dyn_data = ['1, ']

            # Add component, component tab with widgets and call initial setup
            self.add_component(
                connection_name, aquiferName, self.tabControl, key,
                data[key]['type'], self.tabControl.connection_menu,
                self.tabControl.componentsSetupFrame, self, dyn_data, controls)

            # Call additional widgets setup for selected components
            if data[key]['type'] in ['SimpleReservoir', 'AnalyticalReservoir',
                                     'GenericReservoir']:
                locations.load_obs_locations_data(data[key], key)

            if data[key]['type'] in ['TheisReservoir']:
                locations.load_obs_locations_data(data[key], key)
                locations.load_theis_inj_times_rates_data(data[key], key)

            if data[key]['type'] == 'AtmosphericROM':
                for key_arg in ['x_receptor', 'y_receptor']:
                    kwarg_data = ", ".join(str(item) for item in data[key][key_arg])
                    componentVars[key][key_arg].set(kwarg_data)

            if data[key]['type'] in ['MultisegmentedWellbore', 'OpenWellbore',
                                     'CementedWellbore', 'CementedWellboreWR',
                                     'GeneralizedFlowRate']:
                locations.load_locations_data(self, data[key], key)

            if data[key]['type'] == 'LookupTableReservoir':
                locations.load_obs_locations_data(data[key], key)

                code, msg = lutr_tab.finish_load_setup(self, data[key], key)
                if code == 0:
                    self.frames[OpenIAM_Page].tkraise()
                    messagebox.showerror("Error", msg)
                    break
                continue

            if data[key]['type'] == 'PlumeStability':
                code, msg = psa_tab.finish_load_setup(self, data[key], key)
                if code == 0:
                    self.frames[OpenIAM_Page].tkraise()
                    messagebox.showerror("Error", msg)
                    break
                continue

            for output in data[key]['Outputs']:
                if data[key]['type'] == 'MultisegmentedWellbore':
                    if not output in componentVars[key]['outputs']:
                        # Add missing variable
                        componentVars[key]['outputs'][output] = BooleanVar()
                    componentVars[key]['outputs'][output].set(1)
                elif data[key]['type'] == 'GeneralizedFlowRate':
                    if 'CO2_aquifer' in output:
                        componentVars[key]['outputs']['CO2_aquifer'].set(1)
                    if 'brine_aquifer' in output:
                        componentVars[key]['outputs']['brine_aquifer'].set(1)
                else:
                    componentVars[key][output].set(1)

            # For compatibility with older versions of GUI files before change of the
            # parameter name brineResSaturation to aquBrineResSaturation
            if data[key]['type'] == 'MultisegmentedWellbore':
                if 'brineResSaturation' in data[key]['Parameters']:
                    if isinstance(data[key]['Parameters']['brineResSaturation'], dict):
                        data[key]['Parameters']['aquBrineResSaturation'] = \
                            data[key]['Parameters']['brineResSaturation'].copy()
                    else:
                        data[key]['Parameters']['aquBrineResSaturation'] = \
                            data[key]['Parameters']['brineResSaturation']
                    data[key]['Parameters'].pop('brineResSaturation')

            # Get parameters data
            pkeys = data[key]['Parameters'].keys()
            for pkey in pkeys:
                if pkey not in ['pressure', 'CO2saturation', 'brine_rate',
                                'co2_rate', 'ithresh', 'logf',
                                'relativeModel', 'clayType', 'influenceModel']:
                    try:
                        componentVars[key]['Params'][pkey]['value'].set(
                            data[key]['Parameters'][pkey]['value'])
                    except:
                        try:
                            distr_type = data[key]['Parameters'][pkey]['dist']
                            if distr_type == 'uniform':
                                componentVars[key]['Params'][pkey]['distribution'].set('Uniform')
                                componentVars[key]['Params'][pkey]['min'].set(
                                    data[key]['Parameters'][pkey]['min'])
                                componentVars[key]['Params'][pkey]['max'].set(
                                    data[key]['Parameters'][pkey]['max'])

                            if distr_type == 'lognorm':
                                componentVars[key]['Params'][pkey]['distribution'].set('Lognormal')
                                kwargs = {
                                    's': data[key]['Parameters'][pkey]['dist_pars'][0],
                                    'scale': data[key]['Parameters'][pkey]['dist_pars'][2]}
                                dist_pars = self.reparametrize_lognorm_distribution(
                                    kwargs)
                                componentVars[key]['Params'][pkey]['mean'].set(
                                    dist_pars['mean'])
                                componentVars[key]['Params'][pkey]['std'].set(
                                    dist_pars['std'])

                            if distr_type == 'triang':
                                componentVars[key]['Params'][pkey]['distribution'].set(
                                    'Triangular')
                                kwargs = {
                                    'c': data[key]['Parameters'][pkey]['dist_pars'][0],
                                    'loc': data[key]['Parameters'][pkey]['dist_pars'][1],
                                    'scale': data[key]['Parameters'][pkey]['dist_pars'][2]}
                                dist_pars = self.reparametrize_triang_distribution(
                                    kwargs)
                                componentVars[key]['Params'][pkey]['min'].set(
                                    dist_pars['min'])
                                componentVars[key]['Params'][pkey]['max'].set(
                                    dist_pars['max'])
                                componentVars[key]['Params'][pkey]['mode'].set(
                                    dist_pars['mode'])

                            if distr_type == 'norm':
                                componentVars[key]['Params'][pkey]['distribution'].set(
                                    'Normal')
                                componentVars[key]['Params'][pkey]['mean'].set(
                                    data[key]['Parameters'][pkey]['mean'])
                                componentVars[key]['Params'][pkey]['std'].set(
                                    data[key]['Parameters'][pkey]['std'])

                            if distr_type == 'truncnorm':
                                componentVars[key]['Params'][pkey]['distribution'].set(
                                    'Truncated')
                                kwargs = {
                                    'a': data[key]['Parameters'][pkey]['dist_pars'][0],
                                    'b': data[key]['Parameters'][pkey]['dist_pars'][1],
                                    'loc': data[key]['Parameters'][pkey]['dist_pars'][2],
                                    'scale': data[key]['Parameters'][pkey]['dist_pars'][3]}
                                dist_pars = self.reparametrize_truncated_distribution(kwargs)
                                componentVars[key]['Params'][pkey]['min'].set(
                                    dist_pars['min'])
                                componentVars[key]['Params'][pkey]['max'].set(
                                    dist_pars['max'])
                                componentVars[key]['Params'][pkey]['mean'].set(
                                    dist_pars['mean'])
                                componentVars[key]['Params'][pkey]['std'].set(
                                    dist_pars['std'])

                        except:
                            componentVars[key]['Params'][pkey]['distribution'].set('Discrete')
                            values_list = data[key]['Parameters'][pkey]['discrete_vals'][0]
                            weights_list = data[key]['Parameters'][pkey]['discrete_vals'][1]
                            componentVars[key]['Params'][pkey]['values'].set(
                                ', '.join([str(val) for val in values_list]))
                            componentVars[key]['Params'][pkey]['weights'].set(
                                ', '.join([str(val) for val in weights_list]))

                # Seal Horizon keyword parameters
                if pkey in ['relativeModel', 'clayType']:
                    componentVars[key][pkey].set(data[key]['Parameters'][pkey])
                    continue

                # Seal Horizon special parameters
                if pkey == 'influenceModel':
                    model_value = data[key]['Parameters'][pkey]['value']
                    componentVars[key]['Params']['influenceModel']['value'].set(
                        model_value)
                    continue

                if pkey in sh_tab.SH_SPECIAL_PARAMETERS:
                    # The parameters and corresponding componentVars are processed
                    # but we need separate handling of widgets
                    continue

                if pkey == 'ithresh':
                    if data[key]['Parameters'][pkey] == 2:
                        componentVars[key][pkey].set('No Impact')
                    if data[key]['Parameters'][pkey] == 1:
                        componentVars[key][pkey].set('MCL')
                    continue

                if pkey == 'logf':
                    if data[key]['Parameters'][pkey] == 0:
                        componentVars[key][pkey].set('Linear')
                    if data[key]['Parameters'][pkey] == 1:
                        componentVars[key][pkey].set('Log')
                    continue

                # Fault flow component parameters
                if pkey in ['xStart', 'yStart', 'length', 'nSegments', 'strike', 'dip']:
                    continue

                # Get frame name containing the parameter
                par_frame_name = '.'.join([key, pkey, 'frame'])
                if data[key]['type'] == 'GeneralizedFlowRate':
                    label_width = GFR_PARAMETER_LABEL_WIDTH
                elif data[key]['type'] == 'FaultLeakage':
                    label_width = FL_PARAMETER_LABEL_WIDTH
                else:
                    label_width = PARAMETER_LABEL_WIDTH
                self.change_distribution(
                        self.nametowidget(self.getvar(par_frame_name)),
                        pars_label_width=label_width)

            # Some parameters of Seal Horizon component require processing
            # after every other parameters are processed
            if data[key]['type'] == 'SealHorizon':
                # Load cell parameters values
                sh_tab.load_additional_parameters(self, data[key], key)

            # Some parameters of Open Wellbore require processing after initial
            # preprocessing depending on the controls
            if data[key]['type'] == 'OpenWellbore':
                ow_tab.process_crit_pressure_approach_pars(self, data[key], key)

        frame = self.frames[OpenIAM_Page]
        frame.tkraise()

    def process_model_data(self, data):

        """
        Process data and variables associated with Model Tab setup.
        """
        try:  # trying whether time points are provided
            time_data = data['ModelParams']['TimePoints']
        except KeyError:
            componentVars['timeStep'].set(data['ModelParams']['TimeStep'])
            componentVars['endTime'].set(data['ModelParams']['EndTime'])
            componentVars['timePointsInput'].set(0)
            componentVars['timePoints'].set('')
            disable_time_frame_widgets(self.timeFrame)
        else:
            if isinstance(time_data, list):
                componentVars['timePoints'].set(", ".join([
                    str(val) for val in time_data]))
            # Check whether string was provided
            elif isinstance(time_data, str):
                componentVars['timePoints'].set(time_data)
            componentVars['timeStep'].set(1.0)
            componentVars['endTime'].set(50.0)
            componentVars['timePointsInput'].set(1)
            disable_time_frame_widgets(self.timeFrame)

        componentVars['outputDirectory'].set(
            data['ModelParams']['OutputDirectory'])
        componentVars['analysis']['type'].set(
            data['ModelParams']['Analysis']['type'])
        componentVars['logging'].set(data['ModelParams']['Logging'])

        # Setup variables related to outputs
        self.OutputType.set(data['ModelParams'].get('OutputType', True))
        self.GenerateOutputFiles.set(
            data['ModelParams'].get('GenerateOutputFiles', True))
        self.GenerateCombOutputFile.set(
            data['ModelParams'].get('GenerateCombOutputFile', True))
        self.GenerateStatFiles.set(
            data['ModelParams'].get('GenerateStatFiles', True))


        if data['ModelParams']['Analysis']['type'] == 'LHS':
            componentVars['analysis']['size'] = IntVar()
            componentVars['analysis']['seed'] = IntVar()
            componentVars['analysis']['size'].set(
                data['ModelParams']['Analysis']['siz'])
            componentVars['analysis']['seed'].set(
                data['ModelParams']['Analysis']['seed'])

        if data['ModelParams']['Analysis']['type'] == 'Parstudy':
            componentVars['analysis']['nvals'] = StringVar()
            vals = data['ModelParams']['Analysis']['nvals']

            if isinstance(vals, int):
                componentVars['analysis']['nvals'].set(str(vals))
            else:
                componentVars['analysis']['nvals'].set(
                    ', '.join([str(val) for val in vals]))

        self.set_analysis_type(data['ModelParams']['Analysis']['type'])

    def process_strata_data(self, data):
        """
        Process parameters of stratigraphy that should be loaded to
        stratigraphy tab.
        """
        # Setup stratigraphy tab parameters
        num_shale_layers = data['Stratigraphy']['numberOfShaleLayers']['value']
        strata_par_names = ['reservoirThickness'] + [
            'shale{}Thickness'.format(ind) for ind in range(1, num_shale_layers+1)] + [
                'aquifer{}Thickness'.format(ind) for ind in range(1, num_shale_layers)]

        for par_nm in strata_par_names:
            try:
                componentVars['strata']['Params'][par_nm]['value'].set(
                    data['Stratigraphy'][par_nm]['value'])
                componentVars['strata']['Params'][par_nm]['distribution'].set(
                    'Fixed Value')
            except:
                try:
                    distr_type = data['Stratigraphy'][par_nm]['dist']
                    if distr_type == 'uniform':
                        componentVars['strata']['Params'][par_nm]['distribution'].set(
                            'Uniform')
                        componentVars['strata']['Params'][par_nm]['min'].set(
                            data['Stratigraphy'][par_nm]['min'])
                        componentVars['strata']['Params'][par_nm]['max'].set(
                            data['Stratigraphy'][par_nm]['max'])

                    if distr_type == 'norm':
                        componentVars['strata']['Params'][par_nm]['distribution'].set(
                            'Normal')
                        componentVars['strata']['Params'][par_nm]['mean'].set(
                            data['Stratigraphy'][par_nm]['mean'])
                        componentVars['strata']['Params'][par_nm]['std'].set(
                            data['Stratigraphy'][par_nm]['std'])

                    if distr_type == 'lognorm':
                        componentVars['strata']['Params'][par_nm]['distribution'].set(
                            'Lognormal')
                        kwargs = {
                            's': data['Stratigraphy'][par_nm]['dist_pars'][0],
                            'scale': data['Stratigraphy'][par_nm]['dist_pars'][2]}
                        dist_pars = self.reparametrize_lognorm_distribution(kwargs)
                        componentVars['strata']['Params'][par_nm]['mean'].set(
                            dist_pars['mean'])
                        componentVars['strata']['Params'][par_nm]['std'].set(
                            dist_pars['std'])

                    if distr_type == 'truncnorm':
                        componentVars['strata']['Params'][par_nm]['distribution'].set(
                            'Truncated')
                        kwargs = {
                            'a': data['Stratigraphy'][par_nm]['dist_pars'][0],
                            'b': data['Stratigraphy'][par_nm]['dist_pars'][1],
                            'loc': data['Stratigraphy'][par_nm]['dist_pars'][2],
                            'scale': data['Stratigraphy'][par_nm]['dist_pars'][3]}
                        dist_pars = self.reparametrize_truncated_distribution(kwargs)
                        componentVars['strata']['Params'][par_nm]['min'].set(
                            dist_pars['min'])
                        componentVars['strata']['Params'][par_nm]['max'].set(
                            dist_pars['max'])
                        componentVars['strata']['Params'][par_nm]['mean'].set(
                            dist_pars['mean'])
                        componentVars['strata']['Params'][par_nm]['std'].set(
                            dist_pars['std'])

                    if distr_type == 'triang':
                        componentVars['strata']['Params'][par_nm]['distribution'].set(
                            'Triangular')
                        kwargs = {
                            'c': data['Stratigraphy'][par_nm]['dist_pars'][0],
                            'loc': data['Stratigraphy'][par_nm]['dist_pars'][1],
                            'scale': data['Stratigraphy'][par_nm]['dist_pars'][2]}
                        dist_pars = self.reparametrize_triang_distribution(kwargs)
                        componentVars['strata']['Params'][par_nm]['mode'].set(
                            dist_pars['mode'])
                        componentVars['strata']['Params'][par_nm]['min'].set(
                            dist_pars['min'])
                        componentVars['strata']['Params'][par_nm]['max'].set(
                            dist_pars['max'])
                except:
                    if 'discrete_vals' in data['Stratigraphy'][par_nm]:
                        componentVars['strata']['Params'][par_nm]['distribution'].set(
                            'Discrete')
                        values_list = data['Stratigraphy'][par_nm]['discrete_vals'][0]
                        weights_list = data['Stratigraphy'][par_nm]['discrete_vals'][1]
                        componentVars['strata']['Params'][par_nm]['values'].set(
                            ', '.join([str(val) for val in values_list]))
                        componentVars['strata']['Params'][par_nm]['weights'].set(
                            ', '.join([str(val) for val in weights_list]))

            self.change_distribution(self.strata_par_frames[par_nm],
                                     pars_label_width=STRATA_PARAMETER_LABEL_WIDTH)

    def show_frame(self, cont):
        """
        Navigate through all frames that can be imported into the class
        being navigated from.
        """
        frame = self.frames[cont]
        frame.tkraise()

    def show_dashboard(self):
        """
        Method that allows to navigate to and from dashboard.
        """
        frame = self.frames[Dashboard_Page]
        frame.tkraise()

    @staticmethod
    def choose_output_dir(outputDir):
        """
        Set the output directory to be saved in the output dictionary.
        """
        from tkinter.filedialog import askdirectory
        fileDialog = tk.Tk()
        fileDialog.withdraw()

        try:
            dirname = askdirectory(
                initialdir=outputDir.get(),
                title="Choose directory to save outputs")
        except:
            fileDialog.destroy()
        else:
            if dirname != '':
                outputDir.set(dirname)

            fileDialog.destroy()

    @staticmethod
    def choose_file(text_field_var, dialog_title):
        """
        Select file containing specified data.
        """
        from tkinter.filedialog import askopenfilename
        fileDialog = tk.Tk()

        fileDialog.withdraw()
        try:
            filename = askopenfilename(
                initialdir=os.path.dirname(os.path.dirname(os.path.dirname(
                    os.path.abspath(__file__)))),
                title=dialog_title)
        except:
            fileDialog.destroy()
        else:
            if filename != '':
                text_field_var.set(filename)

            fileDialog.destroy()

    @staticmethod
    def Disclaimer_Page():
        """ Create Disclaimer page to allow change frame to work properly.  """
        Disclaimer_Page(tk.Frame)

    @staticmethod
    def Dashboard_Page():
        """ Create Dashboard page. """
        Dashboard_Page(tk.Frame)

    @staticmethod
    def PostProcessor_Page():
        """ Create Postprocessor page. """

        PostProcessor_Page(tk.Frame)

    @staticmethod
    def OpenIAM_Page():
        """ Create Open-IAM page. """
        OpenIAM_Page(tk.Frame)

    def remove_component(self, tab, tabControl, connection_menu):
        """
        Remove component page that is currently displayed on screen.
        """
        MsgBox = messagebox.askquestion(
            "Confirm Removal",
            "Click Yes to confirm removal and No to keep component.")

        if MsgBox == 'yes':
            del componentVars[tabControl.tab(tab)['text']]
            index = componentChoices.index(tabControl.tab(tab)['text'])

            connectionTypes = []

            cmpnt_to_be_removed = componentChoices.pop(index)
            componentTypeDictionary.pop(index)
            connectionsDictionary.pop(index)
            connections.pop(index+1)

            # Remove component data from dictionary d if data is present
            d.pop(cmpnt_to_be_removed, 'No component found')

            connection_menu.children['menu'].delete(0, 'end')

            for c in connectionsDictionary:
                connectionTypes.append(c)
                connection_menu.children['menu'].add_command(
                    label=c,
                    command=lambda con=c: connection_menu.connection.set(con))

            tabControl.select('.!frame.!openiam_page.!notebook.!frame3')
            self.connection.set(connections[0])
            self.componentType.set(COMPONENT_TYPES[0])

            for widget in tabControl.nametowidget(
                    '.!frame.!openiam_page.!notebook.!frame3.!frame.!frame2').winfo_children():
                widget.destroy()

            tabControl.forget(tab)
        else:
            return

    @staticmethod
    def populate_params_dict(cmpnt_par_values, distr_options=None, options=None):
        """ Create dictionary containing parameter values of the component.

        :param cmpnt_par_values: dictionary to contain details of the parameter
            setup
        :type cmpnt_par_values: dict

        :param distr_options: keys of dictionary cmpnt_par_values that need
            to be defined for a given parameter
        :type distr_options: list

        :param options: dictionary with keys 1 and 2; options[1] and options[2]
            are lists; options[1] is a list of keys to be defined as string variable;
            options[2] is a list of keys to be defined as double variable;
        :type options: dict

        """
        if options is None:
            options = {1: ['distribution', 'values', 'weights'],
                       2: ['value', 'mode', 'min', 'max', 'mean', 'std']}
        if distr_options is None:
            distr_options = DISTRIBUTION_OPTIONS

        cmpnt_params_dict = {}
        for par_nm, par_vals in cmpnt_par_values.items():
            cmpnt_params_dict[par_nm] = {}
            # Setup variables
            setup_dict1 = {'distribution': distr_options[0],
                           'values': '{}, {}'.format(par_vals[0], par_vals[3]),
                           'weights': '0.5, 0.5',
                           'ordvalues': str(par_vals[0]),
                           'filename': ''}
            for key in options[1]:
                cmpnt_params_dict[par_nm][key] = StringVar()
                cmpnt_params_dict[par_nm][key].set(setup_dict1[key])

            # Setup variables
            setup_dict2 = {'value': par_vals[0], 'mode': par_vals[4],
                           'min': par_vals[1], 'max': par_vals[2],
                           'mean': par_vals[4], 'std': par_vals[5]}
            for key in options[2]:
                cmpnt_params_dict[par_nm][key] = DoubleVar()
                cmpnt_params_dict[par_nm][key].set(setup_dict2[key])

        return cmpnt_params_dict


    def add_component(self, conn, aqName, tabControl, compName, compType,
                      connection_menu, componentsSetupFrame, controller,
                      dyn_data, controls):
        """
        Add component model.

        Based on component type different parameters are added and
        corresponding values are set.
        """
        toolTip = Pmw.Balloon(self)
        newTab = ttk.Frame(tabControl, padding=10)
        componentName = StringVar()
        componentType = StringVar()
        aquiferName = StringVar()

        try:
            dyn_data_vars = []
            tabControl.add(newTab, text=compName.get())
            componentName.set(compName.get())
            componentType.set(compType.get())
            connection_menu.connection.set(conn.get())
            for _, dyn_data_el in enumerate(dyn_data):
                dyn_data_vars.append(StringVar())
                dyn_data_vars[-1].set(dyn_data_el.get())
            aquiferName.set(aqName.get())
        except:
            dyn_data_vars = []
            tabControl.add(newTab, text=compName)
            componentName.set(compName)
            componentType.set(compType)
            connection_menu.connection.set(conn)

            for _, dyn_data_el in enumerate(dyn_data):
                dyn_data_vars.append(StringVar())
                inp_data_file_path = os.path.join(CODE_DIR, dyn_data_el)
                if os.path.isfile(inp_data_file_path):
                    dyn_data_vars[-1].set(dyn_data_el)
                else:
                    # TODO I'm not sure what kind of dynamic input this covers
                    dyn_data_vars[-1].set(dyn_data_el[:-2])
            aquiferName.set(aqName)

        scanv = tk.Canvas(newTab, relief=tk.SUNKEN)
        scanv.config(width=TAB_SIZE[0], height=TAB_SIZE[1])
        scanv.config(scrollregion=self.get_scroll_region(componentType.get()))
        scanv.config(highlightthickness=0)

        sybar = tk.Scrollbar(newTab, orient='vertical')

        sybar.config(command=scanv.yview)

        scanv.config(yscrollcommand=sybar.set)
        sybar.pack(side=tk.RIGHT, fill=tk.Y)
        scanv.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

        tabType = tk.Frame(scanv)
        tabControl.pack(expand=1, fill="both")
        tabType.grid(row=0, column=0, columnspan=10)

        scanv.create_window((10, 0), window=tabType, anchor='nw')

        # Some components have too many components and other widgets
        # so the placement of the Add/Remove component button has to go after
        # all of them
        last_row = 60

        buttonsFrame = tk.Frame(tabType)
        buttonsFrame.grid(row=last_row, column=0, columnspan=4, pady=20, sticky='w')

        removeComponentButton = ttk.Button(
            buttonsFrame, text='Remove this Component',
            command=lambda: controller.remove_component(
                newTab, tabControl, connection_menu))
        removeComponentButton.grid(row=last_row, column=0)
        toolTip.bind(
            removeComponentButton,
            'Remove current component model and return to Add Components tab.')

        addNextComponentButton = ttk.Button(
            buttonsFrame, text='Add another Component',
            command=lambda: tabControl.select(
                '.!frame.!openiam_page.!notebook.!frame3'))
        addNextComponentButton.grid(row=last_row, column=2, padx=30)
        toolTip.bind(addNextComponentButton, 'Return to Add Components tab.')

        if componentName.get() == '':
            tabControl.forget(newTab)
            messagebox.showerror(
                "Error",
                "You must enter a unique name for each component model.")
            return

        try:
            # Check whether the component with the same name was already added
            componentChoices.index(componentName.get())
        except ValueError:
            # Get add_widgets method corresponding to the type of component
            # All add_widgets methods have to have the same signature
            # and accept arguments in the same order
            dyn_data_list = [var.get() for var in dyn_data_vars]
            add_widgets = self.get_add_widgets_method(componentType.get())
            add_widgets(controller, tabType, componentName.get(),
                        componentType.get(), toolTip, connection_menu.connection.get(),
                        dyn_data_list, aquiferName.get(), controls)

            connection_menu.children['menu'].delete(0, 'end')
            connections.append(componentName.get())

            for c in connections:
                connectionTypes.append(c)
                connection_menu.children['menu'].add_command(
                    label=c, command=lambda con=c: connection_menu.connection.set(con))

            if tabControl.index(tk.END) == 11:
                tabControl.hide('.!frame.!openiam_page.!notebook.!frame3')

            for widget in componentsSetupFrame.winfo_children():
                widget.destroy()

            self.connection.set(connections[0])
            self.componentType.set(COMPONENT_TYPES[0])
            self.componentName_textField.delete(0, tk.END)
            tabControl.select(newTab)
        else:
            # If there is a component with the same name, show error message
            tabControl.forget(newTab)
            componentName.set('')
            messagebox.showerror(
                "Error", "You must enter a unique name for each component model.")
            return

    @staticmethod
    def get_scroll_region(cmpnt_type):
        """Return tuple with a size of scroll region for a component tab."""

        scrollregion_dict = {
            'SimpleReservoir': (0, 0, 0, 1300),
            'AnalyticalReservoir': (0, 0, 0, 1300),
            'GenericReservoir': (0, 0, 0, 1300),
            'TheisReservoir': (0, 0, 0, 1300),
            'LookupTableReservoir': (0, 0, 0, 2000),
            'MultisegmentedWellbore': (0, 0, 0, 1300),
            'CementedWellbore': (0, 0, 0, 800),
            'CementedWellboreWR': (0, 0, 0, 800),
            'OpenWellbore': (0, 0, 0, 1200),
            'GeneralizedFlowRate': (0, 0, 0, 1300),
            'FaultFlow': (0, 0, 0, 2000),
            'FaultLeakage': (0, 0, 0, 800),
            'HydrocarbonLeakage': (0, 0, 0, 800),
            'SealHorizon': (0, 0, 0, 2000),
            'CarbonateAquifer': (0, 0, 0, 2000),
            'AlluviumAquiferLF': (0, 0, 0, 1500),
            'DeepAlluviumAquifer': (0, 0, 0, 1500),
            'DeepAlluviumAquiferML': (0, 0, 0, 1500),
            'FutureGen2Aquifer': (0, 0, 0, 1500),
            'FutureGen2AZMI': (0, 0, 0, 1500),
            'GenericAquifer': (0, 0, 0, 1500),
            'AtmosphericROM': (0, 0, 0, 800),
            'PlumeStability': (0, 0, 0, 1300),
            'ChemicalWellSealing': (0, 0, 0, 800),
            }

        # The component type can contain spaces and parentheses so we remove
        # those if any are present
        return scrollregion_dict[
            cmpnt_type.replace(' ', '').replace('(', '').replace(')', '')]

    @staticmethod
    def get_read_vars_method(cmpnt_type):
        """ Return method which reads data from variables associated with widgets
        on the tab corresponding to the given type of component."""

        method_dict = {
            'SimpleReservoir': src_tab.read_tab_vars,
            'AnalyticalReservoir': arc_tab.read_tab_vars,
            'GenericReservoir': grc_tab.read_tab_vars,
            'TheisReservoir': trc_tab.read_tab_vars,
            'LookupTableReservoir': lutr_tab.read_tab_vars,
            'MultisegmentedWellbore': msw_tab.read_tab_vars,
            'CementedWellbore': cw_tab.read_tab_vars,
            'CementedWellboreWR': cwwr_tab.read_tab_vars,
            'OpenWellbore': ow_tab.read_tab_vars,
            'GeneralizedFlowRate': gfr_tab.read_tab_vars,
            'FaultFlow': ff_tab.read_tab_vars,
            'FaultLeakage': fl_tab.read_tab_vars,
            'HydrocarbonLeakage': hcl_tab.read_tab_vars,
            'SealHorizon': sh_tab.read_tab_vars,
            'CarbonateAquifer': ca_tab.read_tab_vars,
            'AlluviumAquiferLF': aalf_tab.read_tab_vars,
            'DeepAlluviumAquifer': daa_tab.read_tab_vars,
            'DeepAlluviumAquiferML': daaml_tab.read_tab_vars,
            'FutureGen2Aquifer': fgaq_tab.read_tab_vars,
            'FutureGen2AZMI': fgaz_tab.read_tab_vars,
            'GenericAquifer': ga_tab.read_tab_vars,
            'AtmosphericROM': atm_tab.read_tab_vars,
            'PlumeStability': psa_tab.read_tab_vars,
            'ChemicalWellSealing': cws_tab.read_tab_vars,
            }

        # The component type can contain spaces and parentheses so we remove
        # those if any are present
        return method_dict[
            cmpnt_type.replace(' ', '').replace('(', '').replace(')', '')]

    @staticmethod
    def get_add_widgets_method(cmpnt_type):
        """ Return method which adds widgets to the component's tab."""

        # All add_widgets methods have to have the same signature
        # and accept arguments in the same order.
        method_dict = {
            'SimpleReservoir': src_tab.add_widgets,
            'AnalyticalReservoir': arc_tab.add_widgets,
            'GenericReservoir': grc_tab.add_widgets,
            'TheisReservoir': trc_tab.add_widgets,
            'LookupTableReservoir': lutr_tab.add_widgets,
            'MultisegmentedWellbore': msw_tab.add_widgets,
            'CementedWellbore': cw_tab.add_widgets,
            'CementedWellboreWR': cwwr_tab.add_widgets,
            'OpenWellbore': ow_tab.add_widgets,
            'GeneralizedFlowRate': gfr_tab.add_widgets,
            'FaultFlow': ff_tab.add_widgets,
            'FaultLeakage': fl_tab.add_widgets,
            'HydrocarbonLeakage': hcl_tab.add_widgets,
            'SealHorizon': sh_tab.add_widgets,
            'CarbonateAquifer': ca_tab.add_widgets,
            'AlluviumAquiferLF': aalf_tab.add_widgets,
            'DeepAlluviumAquifer': daa_tab.add_widgets,
            'DeepAlluviumAquiferML': daaml_tab.add_widgets,
            'FutureGen2Aquifer': fgaq_tab.add_widgets,
            'FutureGen2AZMI': fgaz_tab.add_widgets,
            'GenericAquifer': ga_tab.add_widgets,
            'AtmosphericROM': atm_tab.add_widgets,
            'PlumeStability': psa_tab.add_widgets,
            'ChemicalWellSealing': cws_tab.add_widgets,
            }

        # The component type can contain spaces and parentheses so we remove
        # those if any are present
        return method_dict[
            cmpnt_type.replace(' ', '').replace('(', '').replace(')', '')]

    def run_simulation(self):
        """ Save and run simulation. """
        self.runSim_button.config(text="Simulation Running")
        self.runSim_button.config(state="disabled")
        self.populate_dictionary(ask_to_save=False)

        filename = self.sim_file
        run_file = os.path.join(SOURCE_DIR, 'openiam', 'openiam_cf.py')
        # Add quotation marks around path to the control file to avoid problems
        # with spaces in paths
        run_command = 'python "{0}" --file "{1}" --binary True'.format(
            run_file, filename)
        os.system(run_command)

        self.runSim_button.config(state="enabled")
        self.runSim_button.config(text="RUN SIMULATION")

    def set_analysis_type(self, analysis):
        """ Change the analysis type of the simulation. """
        for widget in self.analysisFrame.winfo_children():
            widget.destroy()

        # Create and configure common widgets
        analysis_label = ttk.Label(
            self.analysisFrame, text="Analysis:", width=MODEL_TAB_LABEL_WIDTH2)
        analysis_menu = tk.OptionMenu(
            self.analysisFrame, componentVars['analysis']['type'],
            *ANALYSIS_TYPES, command=lambda _: self.set_analysis_type(
                componentVars['analysis']['type'].get()))
        analysis_menu.config(width=MODEL_TAB_MENU_WIDTH)
        analysis_label.grid(row=0, column=0, pady=5, padx=5, sticky='w')
        analysis_menu.grid(row=0, column=1, pady=5, padx=5, sticky='w')

        setup_dict = {'size': 30, 'nvals': '5', 'seed': random.randrange(1000)}

        if analysis == 'LHS':
            for key in ['size', 'seed']:
                if not key in componentVars['analysis']:
                    componentVars['analysis'][key] = IntVar()
                    componentVars['analysis'][key].set(setup_dict[key])
            analysis_size_label = ttk.Label(
                self.analysisFrame, text="Size:", width=MODEL_TAB_LABEL_WIDTH3)
            analysis_size_txtField = tk.Entry(
                self.analysisFrame, width=MODEL_TAB_ENTRY_WIDTH,
                textvariable=componentVars['analysis']['size'])
            analysis_seed_label = ttk.Label(
                self.analysisFrame, text="Seed:", width=MODEL_TAB_LABEL_WIDTH3)
            analysis_seed_txtField = tk.Entry(
                self.analysisFrame, width=MODEL_TAB_ENTRY_WIDTH,
                textvariable=componentVars['analysis']['seed'])

            analysis_size_label.grid(
                row=0, column=2, pady=5, padx=5, sticky='w')
            analysis_size_txtField.grid(
                row=0, column=3, pady=5, padx=5, sticky='w')
            analysis_seed_label.grid(
                row=0, column=4, pady=5, padx=5, sticky='w')
            analysis_seed_txtField.grid(
                row=0, column=5, pady=5, padx=5, sticky='w')

        if analysis == 'Parstudy':
            if not 'nvals' in componentVars['analysis']:
                componentVars['analysis']['nvals'] = StringVar()
                componentVars['analysis']['nvals'].set(setup_dict['nvals'])
            analysis_size_label = ttk.Label(
                self.analysisFrame, text="nvals:", width=MODEL_TAB_LABEL_WIDTH3)
            analysis_size_txtField = tk.Entry(
                self.analysisFrame, width=MODEL_TAB_ENTRY_WIDTH,
                textvariable=componentVars['analysis']['nvals'])

            analysis_size_label.grid(
                row=0, column=2, pady=5, padx=5, sticky='w')
            analysis_size_txtField.grid(
                row=0, column=3, pady=5, padx=5, sticky='w')

    @staticmethod
    def reparametrize_triang_distribution(kwargs):
        """ Change parameters of the triangular distribution."""
        output = {}
        if ('mode' in kwargs) and ('min' in kwargs) and ('max' in kwargs):
            # We need to find c and scale
            output['loc'] = kwargs['min']
            output['scale'] = kwargs['max'] - kwargs['min']
            output['c'] = (kwargs['mode'] - kwargs['min'])/(
                kwargs['max'] - kwargs['min'])
        elif ('c' in kwargs) and ('scale' in kwargs) and ('loc' in kwargs):
            # We need to find maximum and mode values
            output['min'] = kwargs['loc']
            output['max'] = kwargs['loc'] + kwargs['scale']
            output['mode'] = kwargs['loc'] + kwargs['c']*kwargs['scale']

        return output

    @staticmethod
    def reparametrize_lognorm_distribution(kwargs):
        """ Change parameters of the lognormal distribution. """
        output = {}
        if ('mean' in kwargs) and ('std' in kwargs):
            # We need to find s and scale
            output['s'] = kwargs['std']
            output['loc'] = 0
            output['scale'] = np.exp(kwargs['mean'])
        elif ('s' in kwargs) and ('scale' in kwargs):
            # We need to find mean and standard deviation
            output['mean'] = np.log(kwargs['scale'])
            output['std'] = kwargs['s']

        return output

    @staticmethod
    def reparametrize_truncated_distribution(kwargs):
        """ Change parameters of the truncated normal distribution. """
        output = {}
        if ('min' in kwargs) and ('max' in kwargs):
            # We need to find a, b
            output['a'] = (kwargs['min'] - kwargs['mean'])/kwargs['std']
            output['b'] = (kwargs['max'] - kwargs['mean'])/kwargs['std']
            output['loc'] = kwargs['mean']
            output['scale'] = kwargs['std']
        elif ('a' in kwargs) and ('b' in kwargs):
            # We need to find min and max
            output['min'] = kwargs['a']*kwargs['scale'] + kwargs['loc']
            output['max'] = kwargs['b']*kwargs['scale'] + kwargs['loc']
            output['mean'] = kwargs['loc']
            output['std'] = kwargs['scale']

        return output

    def setup_parameter_frame(self, frame, par_name, par_bounds,
                              par_label_text, tool_tip_text,
                              par_name_label_width, distr_menu_width,
                              distr_arg_label_width, text_field_width,
                              distr_options, par_args_variables, cmpnt_nm, tool_tip):
        """ Create, configure and place widgets for a given parameter frame."""
        # Update attributes of frame
        frame.labelText = par_label_text
        frame.component = cmpnt_nm
        frame.distType = par_args_variables['distribution']
        frame.paramVars = par_args_variables
        frame.toolTipText = tool_tip_text
        frame.text = par_name
        frame.par_bounds = {'lower_bound': None,
                            'upper_bound': None,
                            'discrete_bounds': None}
        frame.distr_options = distr_options

        # Copy provided values
        for key in par_bounds:
            frame.par_bounds[key] = par_bounds[key]

        # Create common widgets
        # Label with parameter name and units
        par_name_label = ttk.Label(frame, width=par_name_label_width,
                                   text=par_label_text)
        # Distribution menu
        distr_menu = tk.OptionMenu(
            frame, par_args_variables['distribution'],
            *distr_options, command=lambda _: self.change_distribution(frame,
               pars_label_width=par_name_label_width))
        # Value argument label
        value_label = ttk.Label(frame, text='Value:',
                                width=distr_arg_label_width)

        if frame.par_bounds['discrete_bounds'] is not None:
            standard_tooltip_text = (
                'Set value of {}.\nPossible values are \n{}.'.format(
                    frame.toolTipText,
                    self.reformat_list_presentation(
                        frame.par_bounds['discrete_bounds'])))
        else:
            standard_tooltip_text=(
                'Set value of {}.\nPossible values are between {} and {}.'.format(
                    frame.toolTipText,
                    frame.par_bounds['lower_bound'],
                    frame.par_bounds['upper_bound']))

        # Value argument text field
        value_entry = ParameterEntry(
            frame, par_name, par_args_variables['value'],
            text_field_width, tool_tip,
            standard_tooltip_text=standard_tooltip_text,
            **frame.par_bounds)

        # Configure widgets on the frame
        distr_menu.config(width=distr_menu_width)
        tool_tip.bind(distr_menu, 'Select distribution for {}.'.format(
            frame.toolTipText))

        # Place widgets on the frame
        par_name_label.grid(row=0, column=0, padx=5, sticky='w')
        distr_menu.grid(row=0, column=1, padx=5, sticky='w')
        value_label.grid(row=0, column=2, padx=5, sticky='w')
        value_entry.grid(row=0, column=3, padx=5, sticky='w')

        # Save frame
        self.setvar(name='.'.join([cmpnt_nm, par_name, 'frame']), value=frame)


def ask_exit():
    if messagebox.askokcancel("Exit", "Do you want to close the application?"):
        app.destroy()

if __name__ == "__main__":
    app = NRAPOpenIAM()
    app.wm_title('NRAP-Open-IAM')
    app.geometry("{}x{}".format(APP_SIZE[0], APP_SIZE[1]))
    # app.protocol("WM_DELETE_WINDOW", ask_exit)
    app.mainloop()
