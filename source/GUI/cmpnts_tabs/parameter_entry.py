"""
This module contains class for NRAP-Open-IAM parameters entries.
"""
import tkinter as tk
import numpy as np

class ParameterEntry(tk.Entry):
    def __init__(self, master, name, text_variable, width, tool_tip_link,
                 standard_tooltip_text,
                 lower_bound=None, upper_bound=None,
                 discrete_bounds=None,
                 to_validate=True, **kwargs):
        """
        discrete_bounds - list of values that parameter has to be in

        **kwargs: optional dictionary of additional arguments
            discrete_bounds_msg: list of 4-elements array of alternative messages
            for violation of discrete bounds;
        """

        self.master = master
        self.name = name
        self.to_validate = to_validate
        self.red_flag = 0
        self.standard_tool_tip = standard_tooltip_text
        self.kwargs = kwargs

        # Tooltip
        self.tool_tip = tool_tip_link

        # Validation parameters
        if discrete_bounds is not None:
            self.discrete_bounds = discrete_bounds
        elif (lower_bound is not None) and (upper_bound is not None):
            self.lower_bound = lower_bound
            self.upper_bound = upper_bound
            self.discrete_bounds = None

        # Validation command
        vcmd = (self.master.register(self.validate), "%P")

        # Initialize entry
        tk.Entry.__init__(self, master,
                          textvariable=text_variable, width=width,
                          validate="key", validatecommand=vcmd)

    def turn_on_red(self):
        self.configure(bg="#ffcccc")
        self.red_flag = 1

    def turn_on_white(self):
        self.configure(bg="white")
        self.red_flag = 0

    def validate(self, proposed_value):
        # If entry does not require validation return True
        if not self.to_validate:
            self.tool_tip.bind(self, self.standard_tool_tip)
            return True

        value_is_valid = 1

        if proposed_value.strip() == "":
            value_is_valid = 0
            self.tool_tip.bind(self, "No value specified.")
        else:
            values = []
            try:
                values = [float(proposed_value)]
            except:
                str_values = proposed_value.split(',')
                try:
                    # Strip spaces and convert strings to numbers
                    values = [float(val.strip()) for val in str_values]

                except:
                    value_is_valid = 0
                    self.tool_tip.bind(self, "This is not a correct input.")

            if values:
                if self.discrete_bounds is None:
                    for value in values:
                        if value < self.lower_bound:
                            value_is_valid = 0
                            if ('discrete_bounds_msg' in self.kwargs) and (
                                    self.kwargs['discrete_bounds_msg'][0] != ''):
                                self.tool_tip.bind(
                                    self, self.kwargs['discrete_bounds_msg'][0])
                            else:
                                self.tool_tip.bind(self, ''.join([
                                    "The parameter value(s) must not be less ",
                                    "than {}."]).format(self.lower_bound))
                            break
                        if value > self.upper_bound:
                            value_is_valid = 0
                            if ('discrete_bounds_msg' in self.kwargs) and (
                                    self.kwargs['discrete_bounds_msg'][1] != ''):
                                self.tool_tip.bind(
                                    self, self.kwargs['discrete_bounds_msg'][1])
                            else:
                                self.tool_tip.bind(self, ''.join([
                                    "The parameter value(s) must not be greater ",
                                    "than {}."]).format(self.upper_bound))
                            break
                else:
                    if len(np.unique(values)) < len(values):
                        value_is_valid = 0
                        if ('discrete_bounds_msg' in self.kwargs) and (
                                    self.kwargs['discrete_bounds_msg'][2] != ''):
                                self.tool_tip.bind(
                                    self, self.kwargs['discrete_bounds_msg'][2])
                        else:
                            self.tool_tip.bind(
                                self, "There are duplicate values in the provided list.")
                    else:
                        for value in values:
                            if not value in self.discrete_bounds:
                                value_is_valid = 0
                                if ('discrete_bounds_msg' in self.kwargs) and (
                                    self.kwargs['discrete_bounds_msg'][3] != ''):
                                    self.tool_tip.bind(
                                        self, self.kwargs['discrete_bounds_msg'][3])
                                else:
                                    self.tool_tip.bind(self, ''.join([
                                        "The parameter value must belong to the list ",
                                        "{}."]).format(self.discrete_bounds))
                                break

        if value_is_valid:
            self.tool_tip.bind(self, self.standard_tool_tip)
            self.turn_on_white()
        else:
            self.turn_on_red()
        return True
