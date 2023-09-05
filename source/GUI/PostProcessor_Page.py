import sys
import os
import tkinter as tk
from tkinter import ttk
from tkinter import StringVar
from tkinter import BooleanVar
from tkinter import IntVar
from tkinter import messagebox

# Pwm is imported to allow hover text
import Pmw

from dictionarydata import (
    POSTPROCESSOR_LABEL_WIDTH, BUTTON_WIDTH,
    POSTPROCESSOR_SPINBOX_WIDTH, POSTPROCESSOR_ENTRY_WIDTH,
    POSTPROCESSOR_FILE_ENTRY_WIDTH, POSTPROCESSOR_CHECKBOX_WIDTH,
    POSTPROCESSOR_MENU_WIDTH)
from dictionarydata import LABEL_FONT
from dictionarydata import plotTypes
from dictionarydata import processingTypes

sys.path.append(os.path.join('..', '..', 'source'))
from openiam.visualize.iam_post_processor import IAM_Post_Processor


class PostProcessor_Page(tk.Frame):
    """ Class to create a postprocessor page for OpenIAM GUI. """
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        # Setup all variables needed for the post processing interface
        processingType = StringVar()
        plotType = StringVar()
        titleText = StringVar()
        filename = StringVar()
        subplot = BooleanVar()
        ncols = IntVar()
        simType = StringVar()
        num_sensitivities = IntVar()
        capturePoint = IntVar()
        sampleNumber = IntVar()
        self.outputVars = []

        processingType.set(processingTypes[0])
        plotType.set(plotTypes[0])
        # Define dict with default figure filenames and title
        PLOTTING_OPTIONS = {
            'Plotting': {
                'savefig': '{ob}_time_series',
                'title': 'Results of Simulation'},
            'Map View (AtmROM Only)': {
                'savefig': 'atmplumes_samples_{time}_or_{time_index}',
                'title': 'Map View of $CO_2$ Release at {time} years'},
            'Correlation Coefficients': {
                'savefig': '{ctype}_correlation_coefficients_time_index_{cp}',
                'title': '{ctype} Correlation Coefficients at {ct} years'},
            'Sensitivity Coefficients': {
                'savefig': '{ob}_sensitivity_time_index_{cp}',
                'title': '{ob} Sensitivity Coefficients at {ct} years'},
            'Time Series Sensitivity': {
                'savefig': '{ob}_time_series_sensitivity',
                'title': 'Time Series Sensitivity Coefficients for {ob}'},
            'Multiple Sensitivity Coefficients': {
                'savefig': 'multisensitivities_time_index_{cp}.png',
                'title': 'Sensitivity Coefficients at {ct} years'}}

        titleText.set(PLOTTING_OPTIONS['Plotting']['title'])
        filename.set(PLOTTING_OPTIONS['Plotting']['savefig'])
        subplot.set(False)
        ncols.set(1)
        simType.set('Pearson')
        num_sensitivities.set(5)

        toolTip = Pmw.Balloon(parent)
        setupFrame = ttk.Frame(self)
        browseFrame = ttk.Frame(self)
        nameFrame = ttk.Frame(self)

        def chck_all(chckbox):
            for each in chckbox:
                if each.instate(['selected']) == False:
                    each.invoke()
            # End of method

        def update_name_frame(*args):
            """Assign default title of the plot and filename for Plotting option.
            """
            if plotType.get() == 'Atm Plume Single':
                titleText.set('Map view of $CO_2$ release at {time} years')
                filename.set('atmplumes_samples_{time}_or_{time_index}')
            elif plotType.get() == 'Atm Plume Ensemble':
                titleText.set('Probability of critical distance for {time} years')
                filename.set('atmplumes_ensemble_{time}_or_{time_index}')
            else:
                titleText.set(PLOTTING_OPTIONS['Plotting']['title'])
                filename.set(PLOTTING_OPTIONS['Plotting']['savefig'])
            # End of method

        plotType.trace("w", update_name_frame)

        # Based on the processing type chosen the panel will be filled with
        # different input parts
        def change_processing_type(processing_type):
            """ Change processing type of the simulation data.

            Change processing type of the simulation data and update
            the postprocessor page correspondingly.
            """
            widget_list = setupFrame.winfo_children()
            for widget in widget_list:
                widget.destroy()

            chckbox = []

            processing_label = ttk.Label(browseFrame, text='Processing:',
                                         width=POSTPROCESSOR_LABEL_WIDTH)
            processing_label.grid(row=1, column=0, padx=5, pady=5, sticky='w')

            processing_menu = tk.OptionMenu(browseFrame, processingType,
                                            *processingTypes,
                                            command=change_processing_type)
            processing_menu.grid(row=1, column=1, padx=5, pady=5, sticky='w')
            processing_menu.config(width=POSTPROCESSOR_MENU_WIDTH)
            toolTip.bind(processing_menu,
                         'The only processing type available is Plotting.')

            title_label = ttk.Label(
                nameFrame, text='Title:', width=POSTPROCESSOR_LABEL_WIDTH)
            title_label.grid(row=0, column=0, padx=5, pady=5, sticky='w')

            title_txtField = tk.Entry(nameFrame, textvariable=titleText,
                                      width=POSTPROCESSOR_ENTRY_WIDTH)
            title_txtField.grid(row=0, column=1, padx=5, pady=5, sticky='w')
            toolTip.bind(title_txtField, 'Enter the title of the plot here.')

            filename_label = ttk.Label(nameFrame, text='Filename:',
                                       width=POSTPROCESSOR_LABEL_WIDTH)
            filename_label.grid(row=1, column=0, padx=5, pady=5, sticky='w')

            filename_txtField = tk.Entry(nameFrame, textvariable=filename,
                                         width=POSTPROCESSOR_ENTRY_WIDTH)
            filename_txtField.grid(row=1, column=1, padx=5, pady=5, sticky='w')
            toolTip.bind(filename_txtField, ''.join([
                'Enter the file name you wish to use to save the plot. This ',
                'file will be saved in the current output directory.']))

            observations_frame = ttk.Frame(setupFrame)
            observations_frame.grid(
                row=3, column=1, padx=5, pady=5, columnspan=2)

            post_processing_plot_button = ttk.Button(
                self, text='Plot', width=BUTTON_WIDTH,
                command=lambda: plot_simulation(chckbox))
            post_processing_plot_button.grid(
                row=4, column=1, pady=15, padx=5, sticky='w')
            toolTip.bind(post_processing_plot_button,
                         'Create plots for selected simulation.')

            # Update filename and title variables dependent on the processing type
            filename.set(PLOTTING_OPTIONS[processing_type]['savefig'])
            titleText.set(PLOTTING_OPTIONS[processing_type]['title'])


            def add_observation_frame_widgets():
                """ Add labels and checkboxes to the observations_frame.
                """
                i = 1
                index = 0
                row = 1
                self.outputVars = []
                for output in ipp.filtered_outputs:
                    outputVar = BooleanVar()
                    outputVar.set(0)
                    self.outputVars.append(outputVar)

                    checkbox = ttk.Checkbutton(observations_frame, text=output,
                                               variable=self.outputVars[index],
                                               width=POSTPROCESSOR_CHECKBOX_WIDTH,
                                               onvalue=1,
                                               offvalue=0)
                    checkbox.grid(
                        row=row, column=i, sticky='w', pady=5, padx=5)
                    chckbox.append(checkbox)

                    i = i + 1
                    index = index + 1

                    if i >= 4:
                        row = row + 1
                        i = 1

                check_all = ttk.Button(
                    observations_frame, text='Check All', width=BUTTON_WIDTH,
                    command=lambda: chck_all(chckbox))
                check_all.grid(row=row+1, column=1, sticky='w', pady=5, padx=5)
                # End of add_observation_frame_widgets method


            if processing_type == 'Plotting':
                plotType.set('Time Series')
                observations = ttk.Label(setupFrame, text='Observations:',
                                         width=POSTPROCESSOR_LABEL_WIDTH)
                observations.grid(row=2, column=0, padx=5, pady=5, sticky='w')

                add_observation_frame_widgets()

                plot_label = ttk.Label(setupFrame, text='Plot type:',
                                       width=POSTPROCESSOR_LABEL_WIDTH)
                plot_label.grid(row=0, column=0, padx=5, pady=5, sticky='w')

                plot_menu = tk.OptionMenu(setupFrame, plotType, *plotTypes)
                plot_menu.grid(row=0, column=1, padx=5, pady=5, sticky='w')
                plot_menu.config(width=POSTPROCESSOR_MENU_WIDTH)
                toolTip.bind(plot_menu, 'Choose type of plot to be created.')

                subPlots_label = ttk.Label(setupFrame, text='Use subplots:',
                                           width=POSTPROCESSOR_LABEL_WIDTH)
                subPlots_label.grid(
                    row=4, column=0, padx=5, pady=5, sticky='w')

                subPlots_chckBox = ttk.Checkbutton(setupFrame, variable=subplot)
                subPlots_chckBox.grid(row=4, column=1, padx=5, pady=5, sticky='w')
                toolTip.bind(subPlots_chckBox, ''.join([
                    'Check if you would like to see more than one plot per file.']))

                ncols_label = ttk.Label(setupFrame, text='Number of columns:',
                                        width=POSTPROCESSOR_LABEL_WIDTH)
                ncols_label.grid(row=5, column=0, padx=5, pady=5, sticky='w')

                ncols_spin = tk.Spinbox(setupFrame, from_=1, to=10,
                                        textvariable=ncols,
                                        width=POSTPROCESSOR_SPINBOX_WIDTH)
                ncols_spin.grid(row=5, column=1, padx=5, pady=5, sticky='w')
                toolTip.bind(ncols_spin, ''.join([
                    'Enter the number of columns you would like ',
                    'your plots to be separated into.']))

                if (ipp.analysis_type == 'LHS') or (ipp.analysis_type == 'lhs'):
                    toolTip.bind(processing_menu,
                                 'Select the type of processing you wish to do.')
                    return

            if processing_type == 'Map View (AtmROM Only)':
                plotType.set('Atm Plume Single')
                plot_label = ttk.Label(setupFrame, text='Plot type:',
                                       width=POSTPROCESSOR_LABEL_WIDTH)
                plot_label.grid(row=0, column=0, padx=5, pady=5, sticky='w')

                plot_menu = tk.OptionMenu(setupFrame, plotType,
                                          *['Atm Plume Single', 'Atm Plume Ensemble'])
                plot_menu.grid(row=0, column=1, padx=5, pady=5, sticky='w')
                plot_menu.config(width=POSTPROCESSOR_MENU_WIDTH)
                toolTip.bind(plot_menu, 'Choose type of plot to be created.')
                toolTip.bind(processing_menu,
                             'Select the type of processing you wish to do.')
                sampleNumber.set(1)
                sampleNumber_label = ttk.Label(setupFrame, text="Sample number:",
                                               width=POSTPROCESSOR_LABEL_WIDTH)
                sampleNumber_label.grid(
                    row=1, column=0, padx=5, pady=5, sticky='w')
                sampleNumber_txtField = tk.Entry(
                    setupFrame, textvariable=sampleNumber,
                    width=POSTPROCESSOR_ENTRY_WIDTH//2)
                sampleNumber_txtField.grid(row=1, column=1, padx=5, pady=5, sticky='w')
                toolTip.bind(sampleNumber_txtField,
                             'Enter the sample number for Atm Plume Single map view plot.')

            elif processing_type != 'Plotting':
                toolTip.bind(processing_menu,
                             'Select the type of processing you wish to do.')
                capturePoint.set(len(ipp.time_array)-1)
                capturepoint_label = ttk.Label(setupFrame, text="Capture point:",
                                               width=POSTPROCESSOR_LABEL_WIDTH)
                capturepoint_label.grid(
                    row=0, column=0, padx=5, pady=5, sticky='w')

                capturepoint_spin = tk.Spinbox(setupFrame,
                                               from_=0, to=len(ipp.time_array)-1,
                                               textvariable=capturePoint,
                                               width=POSTPROCESSOR_SPINBOX_WIDTH)
                capturepoint_spin.grid(row=0, column=1, padx=5, pady=5, sticky='w')
                toolTip.bind(capturepoint_spin,
                             'Select index of time step at which results should be plotted.')

            if processing_type == 'Correlation Coefficients':

                capturepoint_label = ttk.Label(setupFrame, text="Type:",
                                               width=POSTPROCESSOR_LABEL_WIDTH)
                capturepoint_label.grid(
                    row=1, column=0, padx=5, pady=5, sticky='w')

                type_menu = tk.OptionMenu(setupFrame, simType,
                                          *['Pearson', 'Spearman'])
                type_menu.grid(row=1, column=1, padx=5, pady=5, sticky='w')
                type_menu.config(width=POSTPROCESSOR_MENU_WIDTH)
                toolTip.bind(type_menu, 'Select the coefficients type.')

                excludes = ttk.Label(setupFrame, text='Excludes:',
                                     width=POSTPROCESSOR_LABEL_WIDTH)
                excludes.grid(row=2, column=0, padx=5, pady=5, sticky='w')

                add_observation_frame_widgets()

            if processing_type.find('Sensitivity') != -1:
                outputs = ttk.Label(setupFrame, text='Observations:',
                                    width=POSTPROCESSOR_LABEL_WIDTH)
                outputs.grid(row=2, column=0, sticky='w', padx=5, pady=5)

                add_observation_frame_widgets()

            if processing_type == 'Time Series Sensitivity':
                sensitivities_label = ttk.Label(setupFrame,
                                                text="Number of sensitivities:",
                                                width=POSTPROCESSOR_LABEL_WIDTH)
                sensitivities_label.grid(
                    row=1, column=0, padx=5, pady=5, sticky='w')

                sensitivities_spin = tk.Spinbox(setupFrame, from_=1, to=100,
                                                textvariable=num_sensitivities,
                                                width=POSTPROCESSOR_SPINBOX_WIDTH)
                sensitivities_spin.grid(
                    row=1, column=1, padx=5, pady=5, sticky='w')
                toolTip.bind(sensitivities_spin,
                             'Select the number of sensitivities you wish to plot.')
            # End of change_processing_type method

        def open_simulation():
            """ Open simulation.

            Allow user to select results of a previously run simulation and get back
            to a list of possible post processing plots as well as statistics.
            """
            widget_list = setupFrame.winfo_children()
            for widget in widget_list:
                widget.destroy()

            from tkinter.filedialog import askdirectory
            fileDialog = tk.Tk()
            fileDialog.withdraw()

            try:
                dirname = askdirectory(initialdir=simulation.get(),
                                       title="Choose directory containing outputs")
            except:
                fileDialog.destroy()
            else:
                if dirname != '':
                    simulation.set(dirname)

                fileDialog.destroy()

            global ipp, plotTypes, processingTypes

            plotTypes = ['Time Series']
            processingTypes = ['Plotting']
            plotType.set(plotTypes[0])
            processingType.set(processingTypes[0])

            try:
                ipp = IAM_Post_Processor(dirname)
            except:
                messagebox.showerror(
                    'Error', ''.join([
                        'The directory you have selected is not',
                        ' a valid NRAP-Open-IAM output directory.']))
                return

            if ipp.analysis_type.lower() in ['parstudy', 'lhs']:
                toolTipText = 'Choose plot type you wish to use.'
                if len(plotTypes) == 1:
                    plotTypes.append('Time Series Stats')
                    plotTypes.append('Time Series and Stats')

#                if ipp.spec_comps['AtmosphericROM']:
#                    if len(plotTypes) == 3:
#                        plotTypes.append('Atm Plume Single')
#                        plotTypes.append('Atm Plume Ensemble')

            add_post_processing_types()
            change_processing_type('Plotting')
            # End of method

        def add_post_processing_types():
            """
            Update postprocessor page based on the opened simulation.
            """
            toolTipText = 'The only plot type available is Time Series.'
            setupFrame.grid(
                row=2, column=0, padx=5, pady=5, columnspan=3, sticky='w')
            nameFrame.grid(
                row=3, column=0, padx=5, pady=5, columnspan=3, sticky='w')

            processing_label = ttk.Label(browseFrame, text='Processing:',
                                         width=POSTPROCESSOR_LABEL_WIDTH)
            processing_label.grid(row=1, column=0, padx=5, pady=5, sticky='w')

            processing_menu = tk.OptionMenu(browseFrame, processingType,
                                            *processingTypes,
                                            command=change_processing_type)
            processing_menu.grid(
                row=1, column=1, padx=5, pady=5, sticky='w')
            processing_menu.config(width=POSTPROCESSOR_MENU_WIDTH)
            toolTip.bind(processing_menu,
                         'The only processing type available is Plotting.')

            plot_label = ttk.Label(setupFrame, text='Plot type:',
                                   width=POSTPROCESSOR_LABEL_WIDTH)
            plot_label.grid(row=0, column=0, padx=5, pady=5, sticky='w')

            plot_menu = tk.OptionMenu(setupFrame, plotType, *plotTypes)
            plot_menu.grid(row=0, column=1, padx=5, pady=5, sticky='w')
            plot_menu.config(width=POSTPROCESSOR_MENU_WIDTH)
            toolTip.bind(plot_menu, toolTipText)

            title_label = ttk.Label(nameFrame, text='Title:',
                                    width=POSTPROCESSOR_LABEL_WIDTH)
            title_label.grid(row=0, column=0, padx=5, pady=5, sticky='w')

            title_txtField = tk.Entry(nameFrame, textvariable=titleText,
                                      width=POSTPROCESSOR_ENTRY_WIDTH)
            title_txtField.grid(row=0, column=1, padx=5, pady=5, sticky='w')
            toolTip.bind(title_txtField, 'Enter the title of the plot here.')

            filename_label = ttk.Label(nameFrame, text='Filename:',
                                       width=POSTPROCESSOR_LABEL_WIDTH)
            filename_label.grid(row=1, column=0, padx=5, pady=5, sticky='w')

            filename_txtField = tk.Entry(nameFrame, textvariable=filename,
                                         width=POSTPROCESSOR_ENTRY_WIDTH)
            filename_txtField.grid(row=1, column=1, padx=5, pady=5, sticky='w')
            toolTip.bind(filename_txtField, ''.join([
                'Enter the file name you wish to use to save the plot. ',
                'This file will be saved in the current output directory.']))

            observations = ttk.Label(setupFrame, text='Observations:',
                                     width=POSTPROCESSOR_LABEL_WIDTH)
            observations.grid(row=2, column=0, sticky='w', pady=5, padx=5)

            observations_frame = ttk.Frame(setupFrame)
            observations_frame.grid(
                row=3, column=1, padx=5, pady=5, columnspan=2)

            subPlots_label = ttk.Label(setupFrame, text='Use subplots:',
                                       width=POSTPROCESSOR_LABEL_WIDTH)
            subPlots_label.grid(row=4, column=0, padx=5, pady=5, sticky='w')

            subPlots_chckBox = ttk.Checkbutton(setupFrame, variable=subplot)
            subPlots_chckBox.grid(row=4, column=1, padx=5, pady=5, sticky='w')
            toolTip.bind(
                subPlots_chckBox,
                'Check if you would like to see more than one plot per file.')

            ncols_label = ttk.Label(setupFrame, text='Number of columns:',
                                    width=POSTPROCESSOR_LABEL_WIDTH)
            ncols_label.grid(row=5, column=0, padx=5, pady=5, sticky='w')

            ncols_spin = tk.Spinbox(setupFrame,
                                    from_=1, to=10, textvariable=ncols,
                                    width=POSTPROCESSOR_SPINBOX_WIDTH)
            ncols_spin.grid(row=5, column=1, padx=5, pady=5, sticky='w')
            toolTip.bind(ncols_spin, ''.join([
                'Enter the number of columns you would like ',
                'your plots to be separated into.']))

            post_processing_plot_button = ttk.Button(
                self, text='Plot', width=BUTTON_WIDTH,
                command=lambda: plot_simulation())
            post_processing_plot_button.grid(
                row=4, column=1, pady=15, padx=5, sticky='w')
            toolTip.bind(post_processing_plot_button,
                         'Create plots for selected simulation.')

            if (ipp.analysis_type == 'LHS') or (
                    ipp.analysis_type == 'lhs') and (
                        len(processingTypes) <= 1):
                if ipp.spec_comps['AtmosphericROM']:
                    processingTypes.append('Map View (AtmROM Only)')
                processingTypes.append('Correlation Coefficients')
                processingTypes.append('Sensitivity Coefficients')
                processingTypes.append('Multiple Sensitivity Coefficients')
                processingTypes.append('Time Series Sensitivity')
                toolTip.bind(processing_menu,
                             'Select the type of processing you wish to do.')
                return
            # End of method

        def plot_simulation(chckbox):
            """
            Plot results according to the made selections.
            """
            outputs = []

            if titleText.get() == '':
                messagebox.showerror(
                    "Error", "You must enter a title for the plot.")
                return
            if filename.get() == '':
                messagebox.showinfo("warning", ''.join([
                    "You will need to enter a file name ",
                    "to save the image directly."]))
                return
            for box in chckbox:
                if box.instate(['selected']):
                    if (ipp.spec_comps['AtmosphericROM']) and (
                            box['text'] in ['x_new', 'y_new',
                                            'critical_distance', 'outflag']):
                        outputs = outputs + ipp.atm_cmp_outputs[box['text']]
                    else:
                        outputs.append(box['text'])

            if processingType.get() == 'Plotting':
                if plotType.get() == 'Time Series':
                    ipp.plotter(outputs, 'TimeSeries', title=titleText.get(),
                                subplot={'use': subplot.get(),
                                         'ncols': ncols.get()},
                                savefig=filename.get())
                if plotType.get() == 'Time Series Stats':
                    ipp.plotter(outputs, 'TimeSeriesStats',
                                title=titleText.get(),
                                subplot={'use': subplot.get(),
                                         'ncols': ncols.get()},
                                savefig=filename.get())
                if plotType.get() == 'Time Series and Stats':
                    ipp.plotter(outputs, 'TimeSeriesAndStats',
                                title=titleText.get(),
                                subplot={'use': subplot.get(),
                                         'ncols': ncols.get()},
                                savefig=filename.get())
            if processingType.get() == 'Map View (AtmROM Only)':
                if plotType.get() == 'Atm Plume Single':
                    ipp.plotter(outputs, 'AtmPlumeSingle',
                                savefig=filename.get(), title=titleText.get(),
                                sample_num=sampleNumber.get()-1)
                if plotType.get() == 'Atm Plume Ensemble':
                    ipp.plotter(outputs, 'AtmPlumeEnsemble',
                                savefig=filename.get(), title=titleText.get())

            if processingType.get() == 'Correlation Coefficients':
                ipp.sensitivity_analysis(
                    'CorrelationCoeff', outputs,
                    {'savefig': filename.get().format(
                        ctype=simType.get(), cp=capturePoint.get(),
                        ct=ipp.time_array[capturePoint.get()]/365.25),
                     'capture_point': capturePoint.get(),
                     'ctype': simType.get().lower(),
                     'title': titleText.get().format(
                         ctype=simType.get(), cp=capturePoint.get(),
                         ct=ipp.time_array[capturePoint.get()]/365.25)})
            if processingType.get() == 'Sensitivity Coefficients':
                ipp.sensitivity_analysis('SensitivityCoeff', outputs,
                                         {'savefig': filename.get(),
                                          'title': titleText.get(),
                                          'capture_point': capturePoint.get()})
            if processingType.get() == 'Multiple Sensitivity Coefficients':
                ipp.sensitivity_analysis('MultiSensitivities', outputs,
                                         {'savefig': filename.get().format(
                                             cp=capturePoint.get()),
                                          'title': titleText.get().format(
                                              cp=capturePoint.get(),
                                              ct=ipp.time_array[capturePoint.get()]/365.25),
                                          'capture_point': capturePoint.get()})
            if processingType.get() == 'Time Series Sensitivity':
                ipp.sensitivity_analysis('TimeSeriesSensitivity', outputs,
                                         {'savefig': filename.get(),
                                          'title': titleText.get(),
                                          'capture_point': capturePoint.get(),
                                          'num_sensitivities': num_sensitivities.get()})
            # End of method

        header = ttk.Label(self, text='Post Processor', font=LABEL_FONT)
        header.grid(row=0, column=0, padx=5, pady=5, sticky='w')

        selection_label = ttk.Label(browseFrame, text='Folder:',
                                    width=POSTPROCESSOR_LABEL_WIDTH)
        selection_label.grid(row=0, column=0, padx=5, pady=5, sticky='w')

        simulation = StringVar()
        simulation.set('')

        selection_txtField = tk.Entry(browseFrame, textvariable=simulation,
                                      width=POSTPROCESSOR_FILE_ENTRY_WIDTH)
        selection_txtField.grid(row=0, column=1, padx=5, pady=5, sticky='w')
        toolTip.bind(selection_txtField, ''.join([
            'Enter the output directory of the simulation ',
            'to be used for post processing.']))

        selection_button = ttk.Button(browseFrame, text='Browse',
                                      command=open_simulation,
                                      width=BUTTON_WIDTH)
        selection_button.grid(row=0, column=2, padx=5, pady=5, sticky='w')
        toolTip.bind(selection_button, ''.join([
            'Click to search for directory containing results of simulation ',
            'you wish to use for post processing.']))

        browseFrame.grid(
            row=1, column=0, columnspan=3, padx=5, pady=5, sticky='w')

        returnButton = ttk.Button(self, text="Return to Dashboard",
                                  command=controller.show_dashboard,
                                  width=BUTTON_WIDTH)
        returnButton.grid(row=4, column=0, pady=15, padx=5, sticky='w')
        toolTip.bind(returnButton, 'Click to return to the main page.')
