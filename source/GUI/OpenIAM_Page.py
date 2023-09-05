import os
import tkinter as tk
from tkinter import ttk
from tkinter import StringVar, DoubleVar, BooleanVar
from tkinter import messagebox
import Pmw

from dictionarydata import componentVars, APP_SIZE, TAB_SIZE
from dictionarydata import aquifers
from dictionarydata import connections
from dictionarydata import ANALYSIS_TYPES
from dictionarydata import LOGGING_TYPES
from dictionarydata import COMPONENT_TYPES
from dictionarydata import savedDictionary
from dictionarydata import LABEL_FONT
from dictionarydata import INSTRUCTIONS_FONT
from dictionarydata import BUTTON_WIDTH, FILE_ENTRY_WIDTH
from dictionarydata import (MODEL_TAB_LABEL_WIDTH1, MODEL_TAB_LABEL_WIDTH2,
                            MODEL_TAB_LARGE_LABEL_WIDTH, PARAMETER_LABEL_WIDTH,
                            MODEL_TAB_ENTRY_WIDTH, MODEL_TAB_MENU_WIDTH)
from dictionarydata import (SETUP_LABEL_WIDTH, SETUP_ENTRY_WIDTH, SETUP_MENU_WIDTH)

from cmpnts_tabs import strata_tab
from cmpnts_tabs.locations import add_file_input_widgets


WELL_DYN_INPUT_SETUP_TEXTS = {
    0: ["Pressure [Pa]:",
        ''.join(['Enter pressure data manually or provide path ',
                 'to the file containing data.']),
        'Select file containing pressure data.'],
    1: ["CO{} saturation [-]:".format(u'\u2082'),
        ''.join(
            ['Enter CO{} saturation data manually or provide path ',
             'to the file containing data.']).format(u'\u2082'),
        ''.join(['Select file containing CO{} ',
                 'saturation data.']).format(u'\u2082')]}

AQUIFER_DYN_INPUT_SETUP_TEXTS = {
    0: ["Brine flow rate [kg/s]:",
        ''.join(
            ['Enter brine flow rate data manually or ',
             'provide path to the file containing data.']),
        'Select file containing brine flow rate data.'],
    1: ["CO{} flow rate [kg/s]:".format(u'\u2082'),
        ''.join(
            ['Enter CO{} flow rate data manually or provide path to ',
             'the file containing data.']).format(u'\u2082'),
        'Select file containing CO{} flow rate data.'.format(u'\u2082')],
    2: ["Brine mass [kg]:",
        ''.join(
            ['Enter mass of brine leaked to the aquifer manually ',
             'or provide path to the file containing data.']),
        'Select file containing brine mass data.'],
    3: ["CO{} mass [kg]:".format(u'\u2082'),
        ''.join(
            ['Enter mass of CO{} leaked to the aquifer manually or ',
             'provide path to the file containing data.']).format(
                 u'\u2082'),
        'Select file containing CO{} mass data.'.format(u'\u2082')]}

ATM_DYN_INPUT_SETUP_TEXTS = {
    0: ["CO{} flow rate [kg/s]:".format(u'\u2082'),
        ''.join(
            ['Enter flow rate of CO{} leaked to the atmosphere ',
             'manually or provide path to the file containing ',
             'data.']).format(u'\u2082'),
        'Select file containing CO{} flow rate data.'.format(u'\u2082')]}


# Define default values for the dynamic kwargs for different components
WELL_DYN_INPUT_DEFAULTS = [1.0e+6, 0.0]
AQUIFER_DYN_INPUT_DEFAULTS = [1.0e-7, 1.0e-7, 100.0, 100.0]
ATM_DYN_INPUT_DEFAULTS = [1.0e-5]

# # TODO Commented out for now as it's not working quite as needed
# class Autoresized_Notebook(ttk.Notebook):
#   def __init__(self, master=None, **kw):

#     ttk.Notebook.__init__(self, master, **kw)
#     self.bind("<<NotebookTabChanged>>", self._on_tab_changed)

#   def _on_tab_changed(self, event):
#     event.widget.update_idletasks()

#     tab = event.widget.nametowidget(event.widget.select())
#     event.widget.configure(height=tab.winfo_reqheight(), width=APP_SIZE[0]-70)


def setup_dyn_pars_widgets(
        controller, frame, texts, variables, tool_tip, row_ind=2,
        dialog_title="Choose file with data for dynamic parameter"):
    """
    Create label, entry and button for each dynamic input type.

    :param controller: application
    :type controller: object of class NRAPOpenIAM

    :param frame: frame to which the widgets will be added
    :type frame: ttk.Frame

    :param texts: dictionary with keys (0-(num_of_variables-1)) where
    each item is a list of 1. text for label; 2. tool tip text for field;
    3. tool tip text for button; texts[-1] contains dialog title
    :type texts: dict

    :param variables: variables to keep value in the entry
    :type variables: list of tkinter StringVar type variables

    :param tool_tip: tooltip to which hints for entry field and button will be
        added to
    :type tool_tip: object of class Pmw.Balloon

    :param row_ind: index of row in frame that the first set of widgets
        will be placed
    :type row_ind: int
    """
    num_of_variables = len(texts.keys())
    labels = []
    entry_fields = []
    buttons = []

    # Create widgets
    labels = [ttk.Label(
        frame, width=SETUP_LABEL_WIDTH,
        text=texts[ind][0]) for ind in range(num_of_variables)]

    entry_fields = [tk.Entry(
        frame, width=SETUP_ENTRY_WIDTH,
        textvariable=variables[ind]) for ind in range(num_of_variables)]

    buttons = [ttk.Button(
        frame, text="Browse",
        command=lambda ind=ind: controller.choose_file(
            variables[ind], dialog_title),
        width=BUTTON_WIDTH) for ind in range(num_of_variables)]

    for ind in range(num_of_variables):
        # Place widgets
        labels[ind].grid(row=row_ind+ind, column=0, pady=5, padx=5, sticky='w')
        entry_fields[ind].grid(row=row_ind+ind, column=1, pady=5, padx=5, sticky='w')
        buttons[ind].grid(row=row_ind+ind, column=2, pady=5, padx=25, sticky='w')

        # Add tool tips
        tool_tip.bind(entry_fields[ind], texts[ind][1])
        tool_tip.bind(buttons[ind], texts[ind][2])


def disable_time_frame_widgets(frame):
    """
    Disable/unable widgets on the time frame of Model Tab depending on the chosen option.
    """
    frame_state = {0: 'normal', 1: 'disabled'}
    check_button_state = frame.checkbox_variable.get()

    # Widgets corresponding to end time and time step
    frame.endtime_label.configure(state=frame_state[check_button_state])
    frame.endtime_entry.configure(state=frame_state[check_button_state])
    frame.timestep_label.configure(state=frame_state[check_button_state])
    frame.timestep_entry.configure(state=frame_state[check_button_state])

    # Widgets corresponding to file input entry and button
    frame.timepoints_label.configure(state=frame_state[1-check_button_state])
    frame.filename_entry.configure(state=frame_state[1-check_button_state])
    frame.browse_button.configure(state=frame_state[1-check_button_state])


class OpenIAM_Page(tk.Frame):
    def __init__(self, parent, controller):
        """
        Constructor method for OpenIAMPage.
        """
        tk.Frame.__init__(self, parent)
        self.toolTip = Pmw.Balloon(parent)

        # Setup all variables needed for this frame
        self.controller = controller
        self.controller.connection = StringVar()
        self.controller.connection.set(connections[0])

        # dyn_data_vars are used to pass
        # pressure/saturation or brine/CO2 leakage rates, masses data
        self.dyn_data_vars = []

        componentVars['analysis'] = {}
        componentVars['strata'] = {}

        def show_dashboard():
            """
            Create and show dashboard page.
            """
            currentDictionary = {}

            for key in componentVars:
                try:
                    currentDictionary[key] = componentVars[key].get()
                except:
                    currentDictionary[key] = {}

                    for objt in componentVars[key]:
                        try:
                            currentDictionary[key][objt] = (
                                componentVars[key][objt].get())
                        except:
                            currentDictionary[key][objt] = (
                                componentVars[key][objt])

            if savedDictionary != currentDictionary:
                MsgBox = messagebox.askquestion(
                    'Save data', ''.join(['Would you like to save the data ',
                                          'before switching to the Dashboard?']),
                    default='yes', icon='question', type='yesnocancel')

                if MsgBox == 'cancel':
                    pass
                else:
                    if MsgBox == 'yes':
                        self.controller.populate_dictionary()
                    self.controller.show_dashboard()
            else:
                self.controller.show_dashboard()
            # End of method


        def add_widgets_for_dyn_pars(cmpnt_type, **kwargs):
            """
            Add relevant widgets for each dynamic input type.
            """
            if cmpnt_type == 'well_fault_seal': # works for wells, fault flow and seal components
                texts = WELL_DYN_INPUT_SETUP_TEXTS
            elif cmpnt_type == 'aquifer':
                if kwargs['num_of_variables'] == 4:
                    texts = AQUIFER_DYN_INPUT_SETUP_TEXTS
                else:
                    texts = {0: AQUIFER_DYN_INPUT_SETUP_TEXTS[2],
                             1: AQUIFER_DYN_INPUT_SETUP_TEXTS[3]}

            else: # atmospheric impact component
                texts = ATM_DYN_INPUT_SETUP_TEXTS

            setup_dyn_pars_widgets(
                controller, tabControl.componentsSetupFrame, texts,
                self.dyn_data_vars, self.toolTip, row_ind=2,
                dialog_title="Choose file with data for dynamic parameter")
            # End of method


        def change_connection(conn):
            """
            Set the component model that the new component model will
            have a connection to.
            """
            for widget in componentsSetupFrame.winfo_children():
                self.toolTip.unbind(widget)
                widget.destroy()

            connection_label = ttk.Label(
                componentsSetupFrame, width=SETUP_LABEL_WIDTH, text="Connection:")
            self.controller.connection.set(connections[0])
            connection_menu = tk.OptionMenu(
                componentsSetupFrame, self.controller.connection,
                *connections, command=change_connection)

            connection_menu.config(width=SETUP_MENU_WIDTH)
            connection_label.grid(row=0, column=0, pady=5, padx=5, sticky='w')
            connection_menu.grid(row=0, column=1, pady=5, padx=5, sticky='w')

            if (self.controller.componentType.get().find('Aquifer') != -1) or (
                    self.controller.componentType.get().find('AZMI') != -1):
                aquiferName.set(aquifers[0])
                label = ttk.Label(componentsSetupFrame, width=SETUP_LABEL_WIDTH,
                                  text="Aquifer name:")
                menu = tk.OptionMenu(componentsSetupFrame, aquiferName, *aquifers)

                menu.config(width=SETUP_MENU_WIDTH)
                label.grid(row=1, column=0, pady=5, padx=5, sticky='w')
                menu.grid(row=1, column=1, pady=5, padx=5, sticky='w')

            # # TODO Do not delete
            # # This part will replace lines above if we decide that
            # # fault flow component can provide input for any of the aquifer components
            # if (self.controller.componentType.get().find('Aquifer') != -1) or (
            #         self.controller.componentType.get().find('AZMI') != -1) or (
            #             self.controller.componentType.get().find('FaultFlow') != -1):
            #     aquiferName.set(aquifers[0])
            #     label = ttk.Label(componentsSetupFrame, width=SETUP_LABEL_WIDTH,
            #                       text="Aquifer name:")
            #     menu = tk.OptionMenu(componentsSetupFrame, aquiferName, *aquifers)

            #     menu.config(width=SETUP_MENU_WIDTH)
            #     label.grid(row=1, column=0, pady=5, padx=5, sticky='w')
            #     menu.grid(row=1, column=1, pady=5, padx=5, sticky='w')

            if conn != 'Dynamic Parameters':
                self.controller.connection.set(conn)

            if conn == 'Dynamic Parameters':
                if (self.controller.componentType.get().find('Wellbore') != -1) or (
                        self.controller.componentType.get().find('Fault Flow') != -1) or (
                            self.controller.componentType.get().find('Seal') != -1):
                    self.dyn_data_vars = [StringVar() for val in range(2)]
                    for ind in range(2):
                        self.dyn_data_vars[ind].set(WELL_DYN_INPUT_DEFAULTS[ind])

                    add_widgets_for_dyn_pars('well_fault_seal')

                if (self.controller.componentType.get().find('Aquifer') != -1) or (
                        self.controller.componentType.get().find('AZMI') != -1):
                    if self.controller.componentType.get().find('Generic Aquifer') != -1:
                        self.dyn_data_vars = [StringVar() for val in range(2)]
                        for ind in range(2):
                            self.dyn_data_vars[ind].set(AQUIFER_DYN_INPUT_DEFAULTS[ind+2])
                    else:
                        self.dyn_data_vars = [StringVar() for val in [0, 1, 2, 3]]
                        for ind in range(4):
                            self.dyn_data_vars[ind].set(AQUIFER_DYN_INPUT_DEFAULTS[ind])

                    add_widgets_for_dyn_pars(
                        'aquifer', num_of_variables=len(self.dyn_data_vars))

                if self.controller.componentType.get().find('Atmo') != -1:
                    self.dyn_data_vars = [StringVar()]
                    for ind in range(1):
                        self.dyn_data_vars[ind].set(ATM_DYN_INPUT_DEFAULTS[ind])

                    add_widgets_for_dyn_pars('atmosphere')
            # End of method


        def change_component_type_setup(componentType):
            """
            Update the add component model frame.

            Update the add component model frame to match what is needed
            for each different type of component model.
            """
            global aquifers

            for widget in componentsSetupFrame.winfo_children():
                self.toolTip.unbind(widget)
                widget.destroy()

            aquifers = ['aquifer{}'.format(ind) for ind in range(
                1, componentVars['strata']['Params']['numberOfShaleLayers'].get())]

            if componentType.find('Reservoir') != -1:
                return

            if componentType.find('Leakage') != -1:
                return

            if componentType.find('Plume') != -1:
                return

            if componentType.find('Chemical') != -1:
                return

            if (componentType.find('Generalized') != -1) or (
                    componentType.find('Aquifer') != -1) or (
                        componentType.find('AZMI') != -1):
                aquiferName.set(aquifers[0])
                label = ttk.Label(componentsSetupFrame, width=SETUP_LABEL_WIDTH,
                                  text="Aquifer name:")
                menu = tk.OptionMenu(componentsSetupFrame, aquiferName, *aquifers)

                menu.config(width=SETUP_MENU_WIDTH)
                self.toolTip.bind(menu,
                                  ''.join(['Select which of the aquifer layers ',
                                           'this component model will represent.']))
                label.grid(row=1, column=0, pady=5, padx=5, sticky='w')
                menu.grid(row=1, column=1, pady=5, padx=5, sticky='w')
                if componentType.find('Generalized') != -1:
                    # No further processing is required for this type of components
                    # since generalized flow rate component does not need connections
                    return

            # # TODO This part will replace lines above 345-363 above if we decide that
            # # fault flow component can provide input for any of the aquifer components
            # if (componentType.find('Generalized') != -1) or (
            #         componentType.find('Fault') != -1) or (
            #             componentType.find('Aquifer') != -1) or (
            #                 componentType.find('AZMI') != -1):
            #     aquiferName.set(aquifers[0])
            #     label = ttk.Label(componentsSetupFrame, width=SETUP_LABEL_WIDTH,
            #                       text="Aquifer name:")
            #     menu = tk.OptionMenu(componentsSetupFrame, aquiferName, *aquifers)

            #     menu.config(width=SETUP_MENU_WIDTH)
            #     self.toolTip.bind(menu,
            #                       ''.join(['Select which of the aquifer layers ',
            #                                'this component model will represent.']))
            #     label.grid(row=1, column=0, pady=5, padx=5, sticky='w')
            #     menu.grid(row=1, column=1, pady=5, padx=5, sticky='w')
            #     if (componentType.find('Generalized') != -1):
            #         # No further processing is required for this type of components
            #         # since generalized flow rate component does not need connections
            #         return

            connection_label = ttk.Label(
                componentsSetupFrame, width=SETUP_LABEL_WIDTH, text="Connection:")
            self.controller.connection.set(connections[0])
            connection_menu = tk.OptionMenu(
                componentsSetupFrame, self.controller.connection,
                *connections, command=change_connection)

            connection_menu.config(width=SETUP_MENU_WIDTH)
            self.toolTip.bind(connection_menu,
                              'Select the connection for your component model.')
            connection_label.grid(row=0, column=0, pady=5, padx=5, sticky='w')
            connection_menu.grid(row=0, column=1, pady=5, padx=5, sticky='w')

            if (componentType.find('Wellbore') != -1) or (
                    componentType.find('Seal') != -1) or (
                        componentType.find('Fault Flow') != -1):
                self.dyn_data_vars = [StringVar() for val in range(2)]
                for ind in range(2):
                    self.dyn_data_vars[ind].set(WELL_DYN_INPUT_DEFAULTS[ind])

                add_widgets_for_dyn_pars('well_fault_seal')

            if (componentType.find('Aquifer') != -1) or (
                    componentType.find('AZMI') != -1):
                if componentType.find('Generic Aquifer') != -1:
                    self.dyn_data_vars = [StringVar() for val in range(2)]
                    for ind in range(2):
                        self.dyn_data_vars[ind].set(AQUIFER_DYN_INPUT_DEFAULTS[ind+2])
                else:
                    self.dyn_data_vars = [StringVar() for val in range(4)]
                    for ind in range(4):
                        self.dyn_data_vars[ind].set(AQUIFER_DYN_INPUT_DEFAULTS[ind])

                add_widgets_for_dyn_pars('aquifer', num_of_variables=len(self.dyn_data_vars))

            if componentType.find('Atmo') != -1:
                self.dyn_data_vars = [StringVar()]
                self.dyn_data_vars[0].set(ATM_DYN_INPUT_DEFAULTS[0])

                add_widgets_for_dyn_pars('atmosphere')
            # End of method

        tabControl = ttk.Notebook(self, width=APP_SIZE[0]-70)
        # # TODO Commented out for now as it's not working quite as needed
        # tabControl = Autoresized_Notebook(self)
        style = ttk.Style()
        current_theme = style.theme_use()
        style.theme_settings(
            current_theme, {"TNotebook.Tab": {"configure": {"padding": [20, 5]}}})

        self.controller.tabControl = tabControl
        modelTab = ttk.Frame(tabControl, padding=10)
        modelFrame = ttk.Frame(modelTab, padding=10)
        tabControl.add(modelTab, text="Model")
        tabControl.pack(expand=1, fill="both", padx=10, pady=5)

        if componentVars['simName'].get() == '':
            componentVars['simName'].set('Default')

        def test_val(inStr, acttyp):
            if acttyp == '1': # insert
                if not list(inStr)[-1].isalnum():
                    return False
            return True
            # End of method

        # Setup model frame and content of Model tab
        modelFrame.grid(row=0, column=0, columnspan=4)

        nameFrame = ttk.Frame(modelFrame)
        nameFrame.grid(row=0, column=0, columnspan=6, sticky='w')

        simName_label = ttk.Label(nameFrame, text='Simulation name:',
                                  width=MODEL_TAB_LABEL_WIDTH2)
        simName_label.grid(row=0, column=0, padx=5, pady=(5, 10), sticky='w')

        simName_txtField = tk.Entry(nameFrame, textvariable=componentVars['simName'],
                                    validate="key", width=FILE_ENTRY_WIDTH)
        simName_txtField.grid(row=0, column=1, padx=5, pady=5, sticky='w')
        simName_txtField['validatecommand'] = (
            simName_txtField.register(test_val), '%P', '%d')
        self.toolTip.bind(simName_txtField,
                          'Enter a unique name for the simulation.')

        modelParams_label = ttk.Label(
            modelFrame, text="Model Parameters", font=LABEL_FONT,
            width=MODEL_TAB_LARGE_LABEL_WIDTH)
        modelParams_label.grid(row=1, column=0, columnspan=2, sticky='w')

        timeFrame = ttk.Frame(modelFrame)
        timeFrame.grid(row=2, column=0, columnspan=6, sticky='w')

        componentVars['endTime'] = DoubleVar()
        componentVars['endTime'].set(50.0)
        endTime_label = ttk.Label(timeFrame, text="End time [years]:",
                                  width=MODEL_TAB_LABEL_WIDTH2)
        endTime_txtField = tk.Entry(timeFrame, width=MODEL_TAB_ENTRY_WIDTH,
                                    textvariable=componentVars['endTime'])
        endTime_label.grid(row=0, column=0, pady=5, padx=5, sticky='w')
        endTime_txtField.grid(row=0, column=1, padx=5, pady=5, sticky='w')
        self.toolTip.bind(endTime_txtField, "Enter simulation time.")
        timeFrame.endtime_label = endTime_label
        timeFrame.endtime_entry = endTime_txtField

        componentVars['timeStep'] = DoubleVar()
        componentVars['timeStep'].set(1.0)
        timeStep_label = ttk.Label(timeFrame, text="Time step [years]:",
                                   width=MODEL_TAB_LABEL_WIDTH2)
        timeStep_txtField = tk.Entry(timeFrame, width=MODEL_TAB_ENTRY_WIDTH,
                                     textvariable=componentVars['timeStep'])
        timeStep_label.grid(row=1, column=0, pady=5, padx=5, sticky='w')
        timeStep_txtField.grid(row=1, column=1, padx=5, pady=5, sticky='w')
        self.toolTip.bind(timeStep_txtField, "Enter time step.")
        timeFrame.timestep_label = timeStep_label
        timeFrame.timestep_entry = timeStep_txtField

        componentVars['timePointsInput'] = BooleanVar()
        componentVars['timePointsInput'].set(0)
        file_input_label = ttk.Label(
            timeFrame, text="Use manual or file input for time points:",
            width=PARAMETER_LABEL_WIDTH+10)
        file_input_checkbox = tk.Checkbutton(
            timeFrame, variable=componentVars['timePointsInput'],
            command=lambda: disable_time_frame_widgets(timeFrame))
        file_input_label.grid(
            row=2, column=0, columnspan=2, pady=5, padx=5, sticky='w')
        file_input_checkbox.grid(
            row=2, column=2, pady=5, padx=5, sticky='w')
        self.toolTip.bind(
            file_input_checkbox,
            'Check to use manual or file input for time points.')
        timeFrame.fileinput_label = file_input_label
        timeFrame.fileinput_checkbox = file_input_checkbox
        timeFrame.checkbox_variable = componentVars['timePointsInput']

        componentVars['timePoints'] = StringVar()
        componentVars['timePoints'].set('')
        timePoints_label = ttk.Label(timeFrame, text="Time points [years]:",
                                     width=MODEL_TAB_LABEL_WIDTH2)
        timePoints_label.grid(row=3, column=0, pady=5, padx=5, sticky='w')
        timeFrame.timepoints_label = timePoints_label
        add_file_input_widgets(
            self.controller, timeFrame, self.toolTip, componentVars['timePoints'],
            ''.join(['Enter time points data manually (separated by comma) or ',
                     'provide path to the file containing data.']),
            'Select file containing time points data.',
            'Choose file containing time points data',
            row_ind=3, col_ind=1)
        disable_time_frame_widgets(timeFrame)
        self.controller.timeFrame = timeFrame

        self.controller.analysisFrame = tk.Frame(modelFrame)
        self.controller.analysisFrame.grid(row=3, column=0, columnspan=6, sticky='we')

        componentVars['analysis']['type'] = StringVar()
        componentVars['analysis']['type'].set(ANALYSIS_TYPES[0])
        analysis_label = ttk.Label(
            self.controller.analysisFrame, text="Analysis:", width=MODEL_TAB_LABEL_WIDTH2)
        analysis_menu = tk.OptionMenu(
            self.controller.analysisFrame,
            componentVars['analysis']['type'], *ANALYSIS_TYPES,
            command=lambda _: self.controller.set_analysis_type(
                componentVars['analysis']['type'].get()))
        analysis_menu.config(width=MODEL_TAB_MENU_WIDTH)
        analysis_label.grid(row=0, column=0, pady=5, padx=5, sticky='w')
        analysis_menu.grid(row=0, column=1, pady=5, padx=5, sticky='w')
        self.toolTip.bind(analysis_menu,
                          'Select type of analysis to be used for simulation.')

        loggingFrame = ttk.Frame(modelFrame)
        loggingFrame.grid(row=4, column=0, columnspan=6, sticky='w')

        componentVars['logging'] = StringVar()
        componentVars['logging'].set(LOGGING_TYPES[1])
        logging_label = ttk.Label(
            loggingFrame, text="Logging:", width=MODEL_TAB_LABEL_WIDTH2)
        logging_Menu = tk.OptionMenu(
            loggingFrame, componentVars['logging'], *LOGGING_TYPES)
        logging_Menu.config(width=MODEL_TAB_MENU_WIDTH)
        logging_label.grid(row=0, column=0, pady=5, padx=5, sticky='w')
        logging_Menu.grid(row=0, column=1, pady=5, padx=5, sticky='w')
        self.toolTip.bind(logging_Menu,
                          'Select type of logging for simulation.')

        componentVars['outputDirectory'] = StringVar()
        try:
            componentVars['outputDirectory'].set(os.path.join(
                os.path.dirname(os.path.dirname(
                    os.path.dirname(os.path.abspath(__file__)))), 'Output'))
        except:
            componentVars['outputDirectory'].set('~Documents')

        outputFrame1 = ttk.Frame(modelFrame)
        outputFrame1.grid(row=5, column=0, columnspan=6, sticky='w')

        outputDirectory_label = ttk.Label(
            outputFrame1, text="Output directory:", width=MODEL_TAB_LABEL_WIDTH2)
        outputDirectory_txtField = tk.Entry(outputFrame1, width=FILE_ENTRY_WIDTH,
                                            textvariable=componentVars['outputDirectory'])
        outputDirectory_browse = ttk.Button(
            outputFrame1, width=BUTTON_WIDTH, text="Browse",
            command=lambda: self.controller.choose_output_dir(
                componentVars['outputDirectory']))
        outputDirectory_label.grid(row=1, column=0, pady=5, padx=5, sticky='w')
        outputDirectory_txtField.grid(
            row=1, column=1, columnspan=3, pady=5, padx=5, sticky='w')
        outputDirectory_browse.grid(
            row=1, column=4, pady=5, padx=5, sticky='w')
        self.toolTip.bind(outputDirectory_txtField,
                          'Enter or select a location to save simulation outputs.')
        self.toolTip.bind(outputDirectory_browse,
                          'Open file browser and select output directory.')

        # Create frame containing setup of output files
        outputFrame2 = ttk.Frame(modelFrame)
        outputFrame2.grid(row=6, column=0, columnspan=6, sticky='w')

        # Create variables and widgets relevant to generating output files
        componentVars['outputDirectoryGenerate'] = BooleanVar()
        componentVars['outputDirectoryGenerate'].set(0)
        outputDirectoryGenerate_label = ttk.Label(
            outputFrame2, text="Generate output directory:",
            width=MODEL_TAB_LABEL_WIDTH1)
        outputDirectoryGenerate_checkbox = tk.Checkbutton(
            outputFrame2, variable=componentVars['outputDirectoryGenerate'])
        outputDirectoryGenerate_label.grid(
            row=2, column=0, pady=5, padx=5, sticky='w')
        outputDirectoryGenerate_checkbox.grid(
            row=2, column=1, pady=5, padx=5, sticky='w')
        self.toolTip.bind(outputDirectoryGenerate_checkbox,
                          ''.join(['Check to generate a directory name augmented with ',
                                   'timestamp for outputs to be saved.']))

        self.controller.OutputType = BooleanVar()
        self.controller.OutputType.set(True)
        outputType_label = ttk.Label(
            outputFrame2, text="Output orientation:")
        outputType_Selection1 = tk.Radiobutton(
            outputFrame2, variable=self.controller.OutputType,
            text="Column-wise", value=True)
        outputType_Selection2 = tk.Radiobutton(
            outputFrame2, variable=self.controller.OutputType,
            text="Row-wise", value=False)
        self.toolTip.bind(outputType_Selection1, 'Select format of output files.')
        self.toolTip.bind(outputType_Selection2, 'Select format of output files.')

        outputType_label.grid(row=3, column=0, pady=5, padx=5, sticky='w')
        outputType_Selection1.grid(row=3, column=1, pady=5, padx=5, sticky='w')
        outputType_Selection2.grid(row=3, column=2, pady=5, padx=5, sticky='w')

        def buttons_state():
            """
            Sets all relevant checkboxes to be checked if output files
            are to be generated.

            """
            s = self.controller.GenerateOutputFiles.get()

            if s:  # if output files are to be generated
                for each in self.button_set:
                    each.config(state='normal')
            else:
                for each in self.button_set:
                    each.config(state='disable')
                self.controller.GenerateCombOutputFile.set(False)
                self.controller.GenerateStatFiles.set(False)
            # End of method


        self.controller.GenerateOutputFiles = BooleanVar()
        self.controller.GenerateOutputFiles.set(True)
        GenerateOutputFiles_label = ttk.Label(
            outputFrame2, text="Generate Output Files?")
        GenerateOutputFiles_Selection1 = tk.Radiobutton(
            outputFrame2, variable=self.controller.GenerateOutputFiles,
            text="Yes", value=True, command=buttons_state)
        GenerateOutputFiles_Selection2 = tk.Radiobutton(
            outputFrame2, variable=self.controller.GenerateOutputFiles,
            text="No", value=False, command=buttons_state)
        self.toolTip.bind(GenerateOutputFiles_Selection1, 'Generate output files.')
        self.toolTip.bind(GenerateOutputFiles_Selection2, 'Do not generate output files.')

        GenerateOutputFiles_label.grid(row=4, column=0, pady=5, padx=5, sticky='w')
        GenerateOutputFiles_Selection1.grid(row=4, column=1, pady=5, padx=5, sticky='w')
        GenerateOutputFiles_Selection2.grid(row=4, column=2, pady=5, padx=5, sticky='w')

        # GenerateCombOutputFile
        self.controller.GenerateCombOutputFile = BooleanVar()
        self.controller.GenerateCombOutputFile.set(True)
        GenerateCombOutputFile_label = ttk.Label(
            outputFrame2, text="Generate a Combined Output File?")
        GenerateCombOutputFile_Selection1 = tk.Radiobutton(
            outputFrame2, variable=self.controller.GenerateCombOutputFile,
            text="Yes", value=True, state='active')
        GenerateCombOutputFile_Selection2 = tk.Radiobutton(
            outputFrame2, variable=self.controller.GenerateCombOutputFile,
            text="No", value=False, state='active')
        self.toolTip.bind(GenerateCombOutputFile_Selection1,
                          'Generate a combined output file.')
        self.toolTip.bind(GenerateCombOutputFile_Selection2,
                          'Do not generate a combined output file.')

        GenerateCombOutputFile_label.grid(row=5, column=0, pady=5, padx=5, sticky='w')
        GenerateCombOutputFile_Selection1.grid(row=5, column=1, pady=5, padx=5, sticky='w')
        GenerateCombOutputFile_Selection2.grid(row=5, column=2, pady=5, padx=5, sticky='w')

        # GenerateStatFiles
        self.controller.GenerateStatFiles = BooleanVar()
        self.controller.GenerateStatFiles.set(True)
        GenerateStatFiles_label = ttk.Label(
            outputFrame2, text="Generate a Statistics File?")
        GenerateStatFiles_Selection1 = tk.Radiobutton(
            outputFrame2, variable=self.controller.GenerateStatFiles,
            text="Yes", value=True, state='active')
        GenerateStatFiles_Selection2 = tk.Radiobutton(
            outputFrame2, variable=self.controller.GenerateStatFiles,
            text="No", value=False, state='active')
        self.toolTip.bind(GenerateStatFiles_Selection1, 'Generate a statistics file.')
        self.toolTip.bind(GenerateStatFiles_Selection2, 'Do not generate a statistics file.')

        GenerateStatFiles_label.grid(row=6, column=0, pady=5, padx=5, sticky='w')
        GenerateStatFiles_Selection1.grid(row=6, column=1, pady=5, padx=5, sticky='w')
        GenerateStatFiles_Selection2.grid(row=6, column=2, pady=5, padx=5, sticky='w')

        self.button_set = [GenerateCombOutputFile_Selection1,
                           GenerateCombOutputFile_Selection2,
                           GenerateStatFiles_Selection1,
                           GenerateStatFiles_Selection2]

        descriptionFrame = ttk.Frame(modelFrame)
        descriptionFrame.grid(row=7, column=0, columnspan=7,
                              sticky='w', pady=(20, 5))

        descriptionLabel = ttk.Label(
            descriptionFrame,
            text=''.join(['After entering system model parameters',
                          ' proceed to Stratigraphy.']),
            font=INSTRUCTIONS_FONT)
        descriptionLabel.grid(row=0, column=0, columnspan=6,
                              padx=5, pady=15, sticky='w')

        nextButton = ttk.Button(
            descriptionFrame, text='Stratigraphy',
            width=BUTTON_WIDTH, command=lambda: tabControl.select(
                '.!frame.!openiam_page.!notebook.!frame2'))
        nextButton.grid(row=0, column=6, padx=5, pady=15, sticky='w')
        self.toolTip.bind(
            nextButton,
            'Switch to the Stratigraphy tab, the second step of the setup.')

        # Setup stratigraphy tab
        new_tab = ttk.Frame(tabControl, padding=10)

        self.controller.strata_scanv = tk.Canvas(new_tab, relief=tk.SUNKEN)
        self.controller.strata_scanv.config(width=TAB_SIZE[0], height=TAB_SIZE[1])
        self.controller.strata_scanv.config(scrollregion=(0, 0, 0, 0))
        self.controller.strata_scanv.config(highlightthickness=0)

        sybar = tk.Scrollbar(new_tab, orient='vertical')
        sybar.config(command=self.controller.strata_scanv.yview)

        self.controller.strata_scanv.config(yscrollcommand=sybar.set)
        sybar.pack(side=tk.RIGHT, fill=tk.Y)
        self.controller.strata_scanv.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

        stratigraphyTab = tk.Frame(self.controller.strata_scanv)
        stratigraphyTab.grid(row=0, column=0, columnspan=10)

        self.controller.strata_scanv.create_window((10, 0), window=stratigraphyTab, anchor='nw')

        tabControl.add(new_tab, text="Stratigraphy")
        tabControl.pack(expand=1, fill="both")

        strata_tab.add_widgets(self.controller, stratigraphyTab, self.toolTip)

        addComponentTab = ttk.Frame(tabControl, padding=10)
        tabControl.add(addComponentTab, text="Add Components")
        tabControl.pack(expand=1, fill="both")
        addComponentFrame = ttk.Frame(addComponentTab, padding=10)
        addComponentFrame.pack(expand=1, fill="both", anchor=tk.NW)

        addComponent_label = ttk.Label(
            addComponentFrame, text="Add Component", font=LABEL_FONT)
        addComponent_label.grid(row=0, column=0, sticky='w')

        componentName = StringVar()
        aquiferName = StringVar()
        componentName.set('')

        componentNameFrame = ttk.Frame(addComponentFrame)
        componentNameFrame.grid(row=1, column=0, columnspan=2, sticky='w')
        componentName_label = ttk.Label(
            componentNameFrame, width=SETUP_LABEL_WIDTH, text="Component name:")
        self.controller.componentName_textField = tk.Entry(
            componentNameFrame, textvariable=componentName,
            width=SETUP_ENTRY_WIDTH, validate="key")
        self.controller.componentName_textField['validatecommand'] = (
            self.controller.componentName_textField.register(test_val), '%P', '%d')

        componentName_label.grid(row=1, column=0, pady=5, padx=5, sticky='w')
        self.controller.componentName_textField.grid(
            row=1, column=1, pady=5, padx=5, sticky='w')
        self.toolTip.bind(self.controller.componentName_textField,
                          'Assure that the component has a unique name.')

        self.controller.componentType = StringVar()
        self.controller.componentType.set(COMPONENT_TYPES[0])
        componentType_label = ttk.Label(
            componentNameFrame, width=SETUP_LABEL_WIDTH, text="Model type:")
        componentType_Menu = tk.OptionMenu(
            componentNameFrame, self.controller.componentType,
            *COMPONENT_TYPES, command=change_component_type_setup)

        componentType_Menu.config(width=SETUP_MENU_WIDTH)
        self.toolTip.bind(
            componentType_Menu,
            'Select the type of component model to be added to the simulation.')
        componentType_label.grid(row=2, column=0, pady=5, padx=5, sticky='w')
        componentType_Menu.grid(row=2, column=1, pady=5, padx=5, sticky='w')

        connection_menu = tk.OptionMenu(
            addComponentFrame, self.controller.connection, *connections)
        tabControl.connection_menu = connection_menu
        tabControl.connection_menu.connection = self.controller.connection

        componentsSetupFrame = ttk.Frame(addComponentFrame)
        tabControl.componentsSetupFrame = componentsSetupFrame
        componentsSetupFrame.grid(
            row=3, column=0, columnspan=3, sticky='w')

        addComponentButton = ttk.Button(
            addComponentFrame, text="Add Component", width=BUTTON_WIDTH,
            command=lambda: self.controller.add_component(
                self.controller.connection, aquiferName, tabControl,
                componentName, self.controller.componentType, connection_menu,
                componentsSetupFrame, self.controller, self.dyn_data_vars, {}))
        addComponentButton.grid(row=1, column=2, pady=2, padx=25, sticky='nw')
        self.toolTip.bind(
            addComponentButton,
            ''.join(['After configuration is finished click Add Component',
                     ' and switch to the corresponding tab for its setup.']))

        textDescription = ttk.Label(
            addComponentTab, font=INSTRUCTIONS_FONT,
            text=''.join([
                'Add each component model for the system to be ',
                'simulated. After specifying all component models,',
                '\nsave the model and return to Dashboard to run the ',
                'simulation. ']))
        textDescription.pack(anchor=tk.SW, pady=10, padx=5)
#        textDescription.grid(row=6, column=0, pady=(50,5),
#                             padx=5, columnspan=10, sticky='w')

        saveButton = ttk.Button(
            self, text="Save", width=BUTTON_WIDTH,
            command=lambda: self.controller.populate_dictionary())
        saveButton.pack(side='left', padx=10, pady=5)

        cancelButton = ttk.Button(self, text="Return to Dashboard",
                                  command=show_dashboard, width=BUTTON_WIDTH)
        cancelButton.pack(side='right', padx=10, pady=5)

        for key in componentVars:
            if key != 'strata':
                try:
                    savedDictionary[key] = componentVars[key].get()
                except:
                    savedDictionary[key] = {}

                    for objt in componentVars[key]:
                        savedDictionary[key][objt] = componentVars[
                            key][objt].get()
