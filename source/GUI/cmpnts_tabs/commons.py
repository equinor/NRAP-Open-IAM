# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:34:16 2022

"""
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
from openiam import IAM_DIR
from dictionarydata import componentVars

def commons_read_tab_vars(cmpnt_nm, cmpnt_type, parameter_names=None,
                          dynamic_kwarg_names=None, observation_names=None,
                          add_aquifer=False, add_number=False,
                          add_leak_to=False):
    """
    Read values of tkinter variables associated with a given component tab.

    """
    d = {}
    d['type'] = cmpnt_type
    d['connection'] = componentVars[cmpnt_nm]['connection'].get()
    if add_aquifer:
        d['AquiferName'] = componentVars[cmpnt_nm]['aquiferName'].get()
    if add_number:
        d['number'] = abs(componentVars[cmpnt_nm]['number'].get())
    if add_leak_to:
        d['LeakTo'] = componentVars[cmpnt_nm]['LeakTo'].get()

    # Get regular parameters
    if parameter_names is not None:
        d['Parameters'] = {}
        for par_nm in parameter_names:
            d['Parameters'][par_nm] = {}

    if dynamic_kwarg_names is not None:
        # Dynamic keyword arguments
        if "Dynamic" in componentVars[cmpnt_nm]['connection'].get():
            d['connection'] = 'Dynamic Parameters'
            d['DynamicParameters'] = {}
            for key in dynamic_kwarg_names:
                inp_data = componentVars[cmpnt_nm][key]
                # Check whether list of numbers or one number (starting with digits 0-9)
                # are provided
                if ',' in inp_data or inp_data[0] in list(range(10)):
                    # Note: There is no way to enter dynamic input for several
                    # locations in this way; for multiple locations the solution
                    # is to use an input file
                    data = inp_data.split(',')
                    d['DynamicParameters'][key] = []
                    for val in data:
                        d['DynamicParameters'][key].append(float(val.strip()))
                else:
                    inp_data_file_path = os.path.join(IAM_DIR, inp_data)
                    if os.path.isfile(inp_data_file_path):
                        d['DynamicParameters'][key] = inp_data
                    else:
                        raise FileNotFoundError(
                            'File {} is not found.'.format(inp_data_file_path))

    if observation_names is not None:
        model_outputs = []
        for obs_nm in observation_names:
            if componentVars[cmpnt_nm][obs_nm].get():
                model_outputs.append(obs_nm)

        d['Outputs'] = model_outputs
    return d
