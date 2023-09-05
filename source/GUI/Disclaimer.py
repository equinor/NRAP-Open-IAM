# Pmw allows for hover text
import Pmw

import tkinter as tk
from tkinter import ttk

from Dashboard import Dashboard_Page
from dictionarydata import LABEL_FONT


class Disclaimer_Page(tk.Frame):
    """ Disclaimer_Page class"""
    def __init__(self, parent, controller):
        """ Constructor method of Disclaimer page. """
        tk.Frame.__init__(self, parent)
        toolTip = Pmw.Balloon(parent)

        self.grid_columnconfigure(0, weight=1)

        notice_label = ttk.Label(self, text=(
            'NOTICE TO USERS:'), font=LABEL_FONT)
        notice_label.grid(row=0, column=0, pady=(70, 10))

        disclaimer_label = ttk.Label(self, text=(
            """Neither the United States Government, nor any of their """+
            """employees, make any warranty, express or implied, or assumes any \n"""+
            """legal liability or responsibility for the accuracy, """+
            """completeness, or usefulness of any information, apparatus, product, or process \n"""+
            """disclosed, or represents that its use would not infringe privately """+
            """owned rights. Reference herein to any specific commercial product, \n"""+
            """process, or service by trade name, trademark, manufacturer, or """+
            """otherwise does not necessarily constitute or imply its endorsement, \n"""+
            """recommendation, or favoring by the United States Government or any """+
            """agency thereof. The views and opinions of authors expressed \n"""+
            """herein do not necessarily state or reflect those of the United States """+
            """Government or any agency thereof."""), font=('Helvetica', 10))
        disclaimer_label.grid(row=1, column=0, pady=5)
        disclaimer_label.configure(anchor="center")

        okButton = ttk.Button(
            self, text="OK", command=lambda: controller.show_frame(Dashboard_Page))
        okButton.grid(row=2, column=0, pady=15)
        toolTip.bind(
            okButton,
            'Clicking this acknowledges that you understand the Disclaimer and accept it.')
