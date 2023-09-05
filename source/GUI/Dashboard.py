# Pmw allows for hover text
# webbrowser allows me to open a webpage link
# pickle is used to write binary files and read binary files
import os
import sys
import pickle
import webbrowser
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

import Pmw

from OpenIAM_Page import OpenIAM_Page
from PostProcessor_Page import PostProcessor_Page
from dictionarydata import LABEL_FONT


# Save location of source folder in the top level folder
SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(SOURCE_DIR)

from openiam import IAM_DIR


class Dashboard_Page(tk.Frame):
    def __init__(self, parent, controller):
        self.root = parent
        tk.Frame.__init__(self, parent)
        toolTip = Pmw.Balloon(parent)

        def open_NRAP_page():
            """ Open browser with NRAP page."""
            webbrowser.open('https://edx.netl.doe.gov/nrap/', new=2)

        def open_forum_page():
            """ Open browser with forum page."""
            webbrowser.open(
                'https://edx.netl.doe.gov/workspace/forum/nrap-tools', new=2)

        def show_user_guide():
            """ Link to the user guide."""

            os.startfile(os.path.join('..', '..', 'User_Guide.pdf'))

        def show_dev_guide():
            """ Link to the developer guide."""

            os.startfile(os.path.join('..', '..', 'Developer_Guide.pdf'))

        def show_references():
            """ Create a popup window to display all references. """
            messagebox.showinfo("References",
                ''.join(['L. H. Pan, S. W. Webb, and C. M. Oldenburg. ',
                        'Analytical solution for two-phase flow in a ',
                        'wellbore using the drift-flux model. Advances ',
                        'in Water Resources, 34(12):1656–1665, 2011. ',
                        'doi:10.1016/j.advwatres.2011.08.009.\n\n\n',
                        'D. R. Harp, R Pawar, J. W. Carey, and C. W. Gable. ',
                        'Reduced order models of transient CO₂ and brine ',
                        'leakage along abandoned wellbores from geologic ',
                        'carbon sequestration reservoirs. International ',
                        'Journal of Greenhouse Gas Control, ',
                        '45:150 – 162, 2016. ',
                        'URL: http://www.sciencedirect.com/science/article/pii/S1750583615301493, ',
                        'doi:https://doi.org/10.1016/j.ijggc.2015.12.001.\n\n\n',
                        'M. A. Celia, J. M. Nordbotten, B. Court, ',
                        'M. Dobossy, and S. Bachu. Field-scale application ',
                        'of a semi-analytical model for estimation of CO₂ ',
                        'and brine leakage along old wells. International ',
                        'Journal of Greenhouse Gas Control, 5(2):257–269,',
                        ' 2011.\n\n\n',
                        'D. H.  Bacon, et al., Geochemical impacts of carbon ',
                        'dioxide, brine, trace metal and organic leakage ',
                        'into an unconfined, oxidizing limestone aquifer. ',
                        'Energy Procedia, 2014.\n\n\n',
                        'Y. Zhang, et al., Fast estimation of dense gas ',
                        'dispersion from multiple continuous CO₂ surface ',
                        'leakage sources for risk assessment, International ',
                        'Journal of Greenhouse Gas Control, 2016.']))

        def show_acknowledgements():
            """ Create a popup window to display all acknowledgements. """
            messagebox.showinfo("Acknowledgements",
                ''.join(['This work was completed as part of the ',
                         'National Risk Assessment Partnership (NRAP) ',
                         'project. Support for this project came from ',
                         'the U.S. Department of Energy’s (DOE) Office ',
                         'of Fossil Energy’s Carbon Storage Program.  ',
                         'The authors wish to acknowledge Mark Ackiewicz ',
                         '(Division of CCS Research Program Manager), ',
                         'Traci Rodosta (Carbon Storage Technology Manager), ',
                         'Kanwal Mahajan (Carbon Storage Division Director), ',
                         'and M. Kylee Underwood (NRAP Project Manager). ',
                         'We also acknowledge the technical guidance ',
                         'of the NRAP Executive Committee, and the comments ',
                         'and feedback of beta tool testers from the ',
                         'international CCUS community.']))

        def open_simulation():
            """
            Reload a saved simulation from a binary file utilizing the pickle package.
            """
            from tkinter.filedialog import askopenfilename
            fileDialog = tk.Tk()
            fileDialog.withdraw()

            # First the open file dialog will appear and show only open iam files.
            # The next lines will try to catch any error with the opening of the file.
            # If an error occurs the file dialog is destroyed else the file is
            # read into a dictionary and passed to the data file to be
            # a globally accessible file.
            try:
                sim_filename = askopenfilename(
                    initialdir=os.path.join(IAM_DIR, 'examples', 'GUI_Files'),
                    title="Choose simulation file",
                    filetypes=[("Open IAM Files", "*.OpenIAM")])
            except:
                fileDialog.destroy()
            else:
                if sim_filename:
                    # This try only happens if the file name selected is not
                    # blank if any error with the file read or the dictionary
                    # load happens then the error message is displayed and the
                    # file dialog is destroyed
                    try:
                        with open(sim_filename, 'rb') as file:
                            data_dict = pickle.load(file)
                    except:
                        messagebox.showerror(
                            'Error with file',
                            'The file you are attempting to open has a problem.')

                    controller.load_simulation(data_dict)

                fileDialog.destroy()

        tool_label = tk.Label(self, padx=5, pady=5, text=(
            'NRAP-Open-IAM - Main Page'), font=LABEL_FONT, bg="light blue")
        tool_label.grid(row=0, column=0, columnspan=2, sticky='ew')

        left_frame = tk.Frame(self, borderwidth=2)
        left_frame.grid(row=1, column=0, sticky='w')

        model_frame = tk.Frame(left_frame, borderwidth=2)
        model_frame.grid(row=1, column=0)

        buttons_frame = ttk.Frame(model_frame, borderwidth=2)
        buttons_frame.pack(side='top')

        enter_params_button = ttk.Button(
            buttons_frame, text="Enter Parameters",
            command=lambda: controller.show_frame(OpenIAM_Page))
        enter_params_button.grid(row=1, column=0, padx=5)
        toolTip.bind(enter_params_button,
                     'Click to start entering parameters.')

        open_saved_sim_button = ttk.Button(
            buttons_frame, text="Load Simulation", command=open_simulation)
        open_saved_sim_button.grid(row=1, column=1, padx=5)
        toolTip.bind(open_saved_sim_button,
                     'Click to open existing simulation file.')

        post_processing_button = ttk.Button(
            buttons_frame, text='Post Processing',
            command=lambda: controller.show_frame(PostProcessor_Page))
        post_processing_button.grid(row=1, column=2, padx=5)
        toolTip.bind(post_processing_button,
                     'Click to post-process and/or plot results.')

        description_label = tk.Label(
            model_frame, padx=5, pady=5,
            text=''.join(
                ['NRAP-Open-IAM is an open-source ',
                 'Integrated Assessment Model (IAM)\n',
                 'for Phase III of the National Risk Assessment ',
                 'Partnership (NRAP). The goal\n',
                 'of this software is to go beyond risk ',
                 'assessment into risk management\n',
                 'and containment assurance.']), justify=tk.LEFT)
        description_label.pack(side='bottom')

        runStyle = ttk.Style()
        runStyle.configure('Bold.TButton', font=('Times', '20', 'bold'),
                           background='light blue')

        controller.runSim_button = ttk.Button(
            left_frame, text='RUN SIMULATION', width=30, style='Bold.TButton',
            command=controller.run_simulation)
        controller.runSim_button.grid(row=2, column=0)

        logo_frame = tk.Frame(left_frame, borderwidth=2)
        logo_frame.grid(row=3, column=0)

        try:
            self.nrap_logo = tk.PhotoImage(file=os.path.join('images', 'NRAP_logo.png'))
        except:
            self.nrap_logo = tk.PhotoImage(file=os.path.join('images','NRAP_logo.gif'))
        nrap_canvas = tk.Canvas(logo_frame)
        nrap_canvas.configure(width=525, height=200)
        nrap_canvas.create_image(250, 100, image=self.nrap_logo)
        nrap_canvas.grid(row=0, column=0, columnspan=5)

        self.pacificnorthwest_logo = tk.PhotoImage(
            file=os.path.join('images', 'DOE-LABS_NRAP.gif'))
        pacificnorthwest_canvas = tk.Canvas(logo_frame)
        pacificnorthwest_canvas.configure(width=900, height=90)
        pacificnorthwest_canvas.create_image(
            435, 50, image=self.pacificnorthwest_logo)
        pacificnorthwest_canvas.grid(row=1, column=0)

        right_frame = tk.Frame(self, borderwidth=2)
        right_frame.grid(row=1, column=1)

        self.well_logo = tk.PhotoImage(file=os.path.join('images', 'wellbore.gif'))
        well_canvas = tk.Canvas(right_frame)
        well_canvas.configure(width=250, height=250)
        well_canvas.create_image(0, 0, anchor='nw', image=self.well_logo)
        well_canvas.grid(row=0, column=0)

        label_frame = tk.Frame(right_frame)
        label_frame.grid(row=1, column=0)

        version_label = ttk.Label(
            label_frame, text=('Version: 2023-08-a2.7.2'))
        version_label.grid(row=0, column=0, padx=5, pady=5, sticky='w')
        toolTip.bind(version_label, 'Software version')

        contact_label = tk.Label(
            label_frame, padx=5, pady=5, text=('Developer: NRAP'))
        contact_label.grid(row=1, column=0, sticky='w')
        toolTip.bind(contact_label, 'Name of project developer')

        nrap_button = ttk.Button(label_frame, text='Home Page',
            command=open_NRAP_page)
        nrap_button.grid(row=2, column=0, sticky='ew')
        toolTip.bind(nrap_button, 'Open browser at NRAP website')

        forum_button = ttk.Button(label_frame, text='NRAP Tools Forum',
            command=open_forum_page)
        forum_button.grid(row=3, column=0, sticky='ew')
        toolTip.bind(forum_button, 'Open browser at forum page')

        user_guide_button = ttk.Button(
            label_frame, text="User's Guide", command=show_user_guide)
        user_guide_button.grid(row=4, column=0, sticky='ew')
        toolTip.bind(user_guide_button, "Open user's guide")

        dev_guide_button = ttk.Button(
            label_frame, text="Developer's Guide", command=show_dev_guide)
        dev_guide_button.grid(row=5, column=0, sticky='ew')
        toolTip.bind(dev_guide_button, "Open developer's guide")

        reference_button = ttk.Button(
            label_frame, text='References', command=show_references)
        reference_button.grid(row=6, column=0, sticky='ew')
        toolTip.bind(reference_button, 'Link to references')

        acknowledge_button = ttk.Button(label_frame, text='Acknowledgements',
            command=show_acknowledgements)
        acknowledge_button.grid(row=7, column=0, sticky='ew')
        toolTip.bind(acknowledge_button, 'Link to acknowledgements')
