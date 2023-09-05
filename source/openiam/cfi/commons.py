import os
import sys
import re
import logging
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

try:
    from openiam import IAM_DIR
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

import openiam as iam
import openiam.cfi.strata as iam_strata


def process_parameters(component, component_data, name2obj_dict=None):
    """
    Check whether any parameters are provided and process their data.

    Process parameters of the component based on the information
    provided in the component data.

    One can provide input to a parameter in the .yaml file as a string.
    These strings must correspond to the bottom depth of a unit, such as
    shale1Depth or aquifer2Depth. For example, one might want to set the wellTop
    or reservoirDepth parameters for an OpenWellbore in this way. Note that
    reservoirDepth can be represented by shale1Depth. Using these string inputs
    is especially useful when the strike and dip option is used for the
    stratigraphy, such that there are spatial variations in unit depths. The
    code obtains the stratigraphy component created for each component, so the
    depths can vary as needed.
    """
    if ('Parameters' in component_data) and (component_data['Parameters']):
        for key in component_data['Parameters']:
            # The case for "par_name: par_value"
            if not isinstance(component_data['Parameters'][key], dict):
                component_data['Parameters'][key] = {
                    'value': component_data['Parameters'][key], 'vary': False}

            if name2obj_dict is not None and 'value' in component_data['Parameters'][key]:
                if isinstance(component_data['Parameters'][key]['value'], str):
                    str_input = component_data['Parameters'][key]['value']

                    # Get the stratigraphy component
                    strata_comp = name2obj_dict['strata']

                    strata_dict = iam_strata.get_strata_info_from_component(strata_comp)

                    numShaleLayers = strata_dict['numberOfShaleLayers']
                    shaleThicknesses = strata_dict['shaleThicknesses']
                    aquiferThicknesses = strata_dict['aquiferThicknesses']

                    unitType, unitNum, metricType = handle_str_input(
                        numShaleLayers, str_input, component.name, key)

                    if metricType == 'Thickness':
                        if unitType == 'shale':
                            thicknessValue = shaleThicknesses[unitNum - 1]

                        elif unitType == 'aquifer':
                            thicknessValue = aquiferThicknesses[unitNum - 1]

                        component_data['Parameters'][key]['value'] = thicknessValue

                    elif metricType == 'Depth':
                        depthValue = iam_strata.get_unit_depth_from_component(
                            numShaleLayers, strata_comp, unitNumber=unitNum,
                            unitType=unitType, top_or_bottom='bottom')

                        component_data['Parameters'][key]['value'] = depthValue

                    component.add_par(key, **component_data['Parameters'][key])

                    # Reset the value so that other locations for the same type
                    # of component will be handled properly.
                    component_data['Parameters'][key]['value'] = str_input

                else:
                    component.add_par(key, **component_data['Parameters'][key])

            else:
                component.add_par(key, **component_data['Parameters'][key])


def process_dynamic_inputs(component, component_data, array_like=False,
                           check_second_dim=False, **kwargs):
    """
    Check whether any dynamic inputs are provided and process their data.

    Process dynamic inputs of the component based on the information
    provided in the component data.

    component - object of ComponentModel class
    component_data - dictionary containing data about setup of the component
    array_like - flag showing whether slice of the data provided through
    dynamic input setup of the component should be array like

    if array_like is True and check_second_dim is True as well, it allows to check whether
    the slice of the data (second dimension of the data) has the right size.

    In kwargs:
        size
        quantity_to_compare
        For size, quantity_to_compare see documentation for
        check_size method
    """
    if ('DynamicParameters' in component_data) and (component_data['DynamicParameters']):
        for key in component_data['DynamicParameters']:
            if key != 'structure':
                inp_data = component_data['DynamicParameters'][key]
                if isinstance(inp_data, str):  # if filename is provided
                    inp_data_file_path = os.path.join(IAM_DIR, inp_data)
                    if os.path.isfile(inp_data_file_path):
                        data = np.genfromtxt(inp_data_file_path,
                                             delimiter=",", dtype='f8',
                                             comments='#')
                        if data.ndim == 2:
                            data = data.T
                        else:
                            if array_like:
                                data = data.reshape(data.shape[0], -1)

                    else:
                        err_msg = '{} is not a valid file name.'.format(inp_data)
                        logging.error(err_msg)
                        raise FileNotFoundError(err_msg)
                else:  # if numerical data is provided
                    if array_like:
                        data = np.array(inp_data)
                        data = data.reshape(data.shape[0], -1)
                    else:
                        data = inp_data

                if check_second_dim:
                    check_size(component, data, kwargs['size'],
                               kwargs['quantity_to_compare'], key)
                component.add_dynamic_kwarg(key, data)


def check_size(component, data, size, quantity_to_compare, dyn_input_name):
    """
    data is a matrix like np.array whose size needs to be checked
    size is a required size of a slice
    quantity_to_compare is a string indicating name of a quantity (if applicable)
    whose size should be equal to the second dimension of the data
    dyn_input_name is a string indicating name of a dynamic input (usually,
    the one used in the model method of component, e.g.,
    pressure/CO2saturation/co2rate, etc.)
    """
    if data.shape[1] != size:
        err_msg = ''.join([
            '{} does not coincide with the second dimension of the provided ',
            'dynamic keyword argument {} for component {}']).format(
                quantity_to_compare, dyn_input_name, component.name)
        logging.error(err_msg)
        raise IndexError(err_msg)


def handle_str_input(num_shale_layers, str_input, cmpnt_name, par_name):
    """
    This function takes string input provided for a parameter in a .yaml file
    and returns the corresponding unitType ('shale' or 'aquifer'), unitNum
    (unit number, 1, 2, 3, etc.), and metricType ('Thickness' or 'Depth', where
    'Depth' corresponds with the bottom depth of the unit).

    :param num_shale_layers: number of shale layers
    :type num_shale_layers: int

    :param str_input: input for parameter provided instead of value
    :type str_input: str

    :param cmpnt_name: name of the component object for which the given parameter
        is processed
    :type cmpnt_name: str

    :param par_name: name of the parameter for which the string input is
        provided
    :type par_name: str
    """
    # Create a list of possible stratigraphy parameters
    potential_str_references = [
        'shale{}Depth'.format(ind) for ind in range(1, num_shale_layers+1)]+[
            'shale{}Thickness'.format(ind) for ind in range(1, num_shale_layers+1)]+[
                'aquifer{}Depth'.format(ind) for ind in range(1, num_shale_layers)]+[
                    'aquifer{}Thickness'.format(ind) for ind in range(1, num_shale_layers)]

    # If string provided as parameter value is among possible stratigraphy parameters
    # determine whether it's for shale or aquifer and whether it's depth or thickness
    if str_input in potential_str_references:
        # Match the str_input to the pattern
        out = re.match(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", str_input)
        unitType = out.group(1)
        unitNum = int(out.group(2))
        metricType = out.group(3)

    else:
        err_msg = ''.join([
            'A string was used as input for the parameter {} in the component {}. ',
            'The string did not match any of the expected values, however. ',
            'Expected values are the depths to the bottom of a unit (e.g., '
            'shale1Depth or aquifer2Depth). Check the input ',
            'in the .yaml file.']).format(par_name, cmpnt_name)
        raise KeyError(err_msg)

    return unitType, unitNum, metricType


def get_parameter_val(comp, par_name, sm=None, yaml_data=None):
    """
    Function that checks if a parameter (par_name) is in the stochastic,
    deterministic, composite, or default parameters of a component. If so, the
    parameter value is returned.
    """
    par_val = None

    # If the parameter is a composite parameter that returns as zero (which
    # can happen for lhs or parstudy simulations) or linked to another parameter,
    # the function check_parameter_types() updates run_again_keys and then runs again.
    # When it runs again, it is directed to another component from which to obtain
    # the parameter value (e.g., getting 'aquifer2Depth' from a stratigraphy
    # component in lieu of 'wellTop').
    run_again_keys = [
        'comp', 'par_name', 'numShaleLayers', 'unit_number', 'top_or_bottom']

    run_again_info = {key: None for key in run_again_keys}

    par_val, run_again, run_again_info = check_parameter_types(
        comp, par_name, run_again_info, sm=sm, yaml_data=yaml_data)

    if run_again:
        par_val, run_again, run_again_info = check_parameter_types(
            run_again_info['comp'], run_again_info['par_name'], run_again_info,
            sm=sm, yaml_data=yaml_data)

    return par_val


def check_parameter_types(comp, par_name, run_again_info, sm=None, yaml_data=None):
    """
    This function is used by get_parameter_val() to check if the parameter in
    question (par_name) is in the stochastic (.pars), deterministic, composite,
    linked, or default parameters. It will return the value once the parameter
    is found. If a composite parameter returns as 0 (which happens in lhs and parstudy
    simulations) or if the parameter is linked to another parameter, this function
    provides the component and parameter name that will provide the correct value
    (e.g., 'aquifer2Depth' for the 'wellTop' parameter).
    """
    run_again = False

    par_val = None

    # This list is only meant to include composite parameters, not parameters
    # linked to other parameters (e.g., aqu_thick or top_depth for a GenericAquifer).
    strat_related_pars = ['wellTop', 'reservoirDepth', 'wellDepth']

    if par_name in comp.pars:
        par_val = comp.pars[par_name].value

    elif par_name in comp.deterministic_pars:
        par_val = comp.deterministic_pars[par_name].value

    elif par_name in comp.composite_pars:
        # In cases where a parameter is set as a composite parameter (e.g.,
        # setting the wellTop or reservoirDepth parameter of an OpenWellore
        # only by having 'LeakTo: aquifer2') and the simulation uses lhs or
        # parstudy analysis types, the composite_pars value will be returned
        # as zero. The 'if' statements below address such situations -
        # additional cases may need to be added for new component parameters.
        par_val = comp.composite_pars[par_name].value

        if par_val == 0:
            if 'Depth' in par_name and isinstance(comp, iam.Stratigraphy):
                if 'aquifer' in par_name:
                    unitType = 'aquifer'
                elif 'shale' in par_name:
                    unitType = 'shale'
                elif 'reservoir' in par_name:
                    unitType = 'reservoir'
                    run_again_info['top_or_bottom'] = 'top'

                par_val = iam_strata.get_unit_depth_from_component(
                    run_again_info['numShaleLayers'], comp,
                    unitNumber=run_again_info['unit_number'], unitType=unitType,
                    top_or_bottom=run_again_info['top_or_bottom'])

            elif par_name in strat_related_pars and not yaml_data is None:
                run_again_info['numShaleLayers'] = yaml_data['Stratigraphy'][
                    'numberOfShaleLayers']['value']

                try:
                    run_again_info['comp'] = sm.component_models['strata']
                except KeyError:
                    # Case for spatially variable stratigraphy
                    run_again_info['comp'] = sm.component_models['strata' + comp.name]

                if par_name == 'wellTop':
                    LeakTo = yaml_data[comp.name]['LeakTo']

                    # If the unit number is less than 10, there will be 8 characters
                    # (e.g., 'aquifer2'). If the unit number is >= 10, there will be
                    # 9 characters (e.g., 'aquifer12').
                    run_again_info['unit_number'] = int(LeakTo[7:])

                    run_again_info['top_or_bottom'] = 'bottom'

                    run_again_info['par_name'] = LeakTo + 'Depth'

                elif par_name in ['reservoirDepth', 'wellDepth']:
                    run_again_info['par_name'] = 'shale1Depth'

                    run_again_info['unit_number'] = 1
                    run_again_info['top_or_bottom'] = 'bottom'

                run_again = True

            else:
                # This is where cases for new stratigraphy-related parameters should be entered
                pass

    elif par_name in comp.parlinked_pars and sm is not None:
        par_connection = comp.parlinked_pars[par_name]

        linked_comp_name = par_connection[0:par_connection.index('.')]
        run_again_info['par_name'] = par_connection[(par_connection.index('.') + 1):None]

        run_again_info['comp'] = sm.component_models[linked_comp_name]

        run_again = True

    else:
        if 'Thickness' in par_name and ('shale' in par_name or 'aquifer' in par_name):
            # Parameters like 'shale1Thickness' will not be in the default parameters,
            # but 'shaleThickness' will be.
            if 'shale' in par_name:
                par_name = 'shaleThickness'
            elif 'aquifer' in par_name:
                par_name = 'aquiferThickness'

        par_val = comp.default_pars[par_name].value

    return par_val, run_again, run_again_info
