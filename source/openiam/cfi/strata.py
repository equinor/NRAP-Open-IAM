"""
Created: September 9th, 2022
Last modified: October 21th, 2022

Authors: Nate Mitchell, Veronika Vasylkivska
Leidos supporting NETL
"""

import sys
import os
import logging
import numpy as np
from matplotlib.colors import is_color_like
import pandas as pd


SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(SOURCE_DIR)

import openiam as iam

# Functions in this file are used to generate the colors, alpha values, and
# labels for different units.
# Color used to plot the unit
UNIT_COLOR_DICT = {
    'ReservoirColor': [0.5, 0.5, 0.5],
    'ShaleColor': [1, 0, 0],
    'AquiferColor': [0, 0, 1],
    }

# Default well color and alpha values
WELL_COLOR = [0, 0, 1]
WELL_ALPHA = 1
WELL_ALPHA_FILL = 0.75

# This is used to scale how much higher the alpha of labels / lines is than
# the alpha for filled areas:
# alpha = alphaFill + ((1 - alphaFill) * ALPHA_SCALE_FACTOR)
ALPHA_SCALE_FACTOR = 0.8

# Alpha used when plotting unit labels and the lines on the edges of units
UNIT_ALPHA_DICT = {
    'ReservoirAlpha': 0.9,
    'ShaleAlpha': 0.85,
    'AquiferAlpha': 0.85,
    }

# Alpha used when plotting a filled in area
UNIT_ALPHA_FILL_DICT = {
    'ReservoirAlphaFill': 0.5,
    'ShaleAlphaFill': 0.25,
    'AquiferAlphaFill': 0.25,
    }

# Labels for each unit. This is used in stratigraphy_plot.py (thickness labels
# are too long for that plot type).
UNIT_LABEL_DICT = {
    'ReservoirLabel': 'Reservoir',
    'ShaleLabel': 'Shale {:.0f}',
    'AquiferLabel': 'Aquifer {:.0f}',
    }

# Labels that include the unit thickness. This is used in stratigraphic_column.py.
UNIT_LABEL_THICKNESS_DICT = {
    'ReservoirLabel': 'Reservoir Thickness: {:.2f} m',
    'ShaleLabel': 'Shale {:.0f} Thickness: {:.2f} m',
    'AquiferLabel': 'Aquifer {:.0f} Thickness: {:.2f} m',
    }

DIP_DIRECTION_DEGREE_OPTIONS = [0, 90, 180, 270, 360]
DIP_DIRECTION_OPTIONS = ['N', 'S', 'E', 'W', 'NE', 'NW', 'SE', 'SW']

CHECK_CONDITIONS_MSG = ''.join([
    'Check your input. Alternatively, calculating unit thicknesses with these ',
    'strike and dip values may not be appropriate for these locations and/or ',
    'conditions.'])


def initialize_strata(yaml_data, sm):
    """
    This function handles the processing of stratigraphy input from the .yaml
    file. If spatiallyVariable is not included, then a single stratigraphy
    component is made and the parameters from the .yaml file are added. If
    spatiallyVariable is included in the .yaml file, then the function creates
    a list that will contain all of the stratigraphy components.

    :param yaml_data: dictionary of input values from the .yaml file
    :type yaml_data: dict

    :param sm: system model object
    :type sm: openiam.iam_base_classes.SystemModel

    :returns: strata, sm, spatially_variable_strata
    """

    if 'Stratigraphy' in yaml_data:
        # If spatiallyVariable is in yaml_data, set up strata to be a list
        # containing multiple stratigraphy components. If spatiallyVariable is
        # not is yaml_data, then treat the stratigraphy as flat-lying.
        if 'spatiallyVariable' in yaml_data['Stratigraphy'] or \
                'spatiallyvariable' in yaml_data['Stratigraphy']:
            spatially_variable_strata = True
            strata = []

            # Make one stratigraphy component for the reference point - the
            # location at which the stratigraphy specified in the .yaml file applies
            strata.append(sm.add_component_model_object(iam.Stratigraphy(
                name='strataRefPoint', parent=sm)))

            strat_params = yaml_data['Stratigraphy']
            for parameter in strat_params:
                if parameter not in ['spatiallyVariable', 'spatiallyvariable']:
                    if not isinstance(strat_params[parameter], dict):
                        strat_params[parameter] = {
                            'value': strat_params[parameter],
                            'vary': False}
                    strata[-1].add_par(parameter, **(strat_params[parameter]))

        else:
            spatially_variable_strata = False
            # Add stratigraphy component
            strata = sm.add_component_model_object(iam.Stratigraphy(
                name='strata', parent=sm))

            strat_params = yaml_data['Stratigraphy']
            for parameter in strat_params:
                if not isinstance(strat_params[parameter], dict):
                    strat_params[parameter] = {
                        'value': strat_params[parameter],
                        'vary': False}
                strata.add_par(parameter, **(strat_params[parameter]))

    else:
        spatially_variable_strata = False
        # If 'Stratigraphy' is not in yaml_data, still add stratigraphy comp.
        strata = sm.add_component_model_object(iam.Stratigraphy(
            name='strata', parent=sm))

    return strata, sm, spatially_variable_strata


def process_spatially_variable_strata(strata, component_name, yaml_data,
                                      locations, sm):
    """
    This function handles the creation of stratigraphy components in
    openiam_cf.py in cases where the stratigraphy is spatially variable.

    Note that when unit depths decrease, the required changes in depth are
    accomplished by making the top unit (a shale) thicker. When unit depths
    decrease so much that the unit piches out (i.e., thickness in the
    subsurface goes to zero), the thickness will only go down to the minimum
    value of 1 m.

    :param component_name: name of the component model being handled. Note that
        the component name is used to extract the location data from locations
    :type component_name: str

    :param yaml_data: dictionary of input values from the .yaml file
    :type yaml_data: dict

    :param locations: dictionary linking specific components to the
        corresponding coordinates, number of locations, and component type
        (e.g., 'well': True)
    :type locations: dict

    :param sm: system model object
    :type sm: openiam.iam_base_classes.SystemModel

    :returns: strata, sm
    """

    if 'strikeAndDip' in yaml_data['Stratigraphy']['spatiallyVariable']:
        var_type = 'strikeAndDip'

        strata_var_info = get_strata_var_info_from_yaml(yaml_data)

        strike = strata_var_info['strike']
        dip = strata_var_info['dip']
        dipDirection = strata_var_info['dipDirection']
        coordxRefPoint = strata_var_info['coordxReferencePoint']
        coordyRefPoint = strata_var_info['coordyReferencePoint']

    elif 'LookupTableStratigraphy' in yaml_data['Stratigraphy']['spatiallyVariable']:
        var_type = 'LookupTable'

    # First, find the base name of the component (e.g., excluding '_001') and
    # the number within the name ('001'), which refers to the location.
    base_name = component_name[0:component_name.index('_')]
    comp_number = component_name[(component_name.index('_') + 1):None]
    comp_number = int(comp_number)

    if base_name in locations:
        # Get the coordinates for this component model
        comp_x_val = locations[base_name]['coordx'][comp_number]
        comp_y_val = locations[base_name]['coordy'][comp_number]

    elif 'DistributedTo' in yaml_data[base_name] and len(yaml_data[
            base_name]['DistributedTo']) != 0:
        # Get the coordinates from the component model connected to this one
        comp_x_val = locations[yaml_data[base_name]['DistributedTo'][0]]['coordx'][comp_number]
        comp_y_val = locations[yaml_data[base_name]['DistributedTo'][0]]['coordy'][comp_number]

    elif 'Connection' in yaml_data[base_name]:
        # Get the coordinates from the component model connected to this one
        comp_x_val = locations[yaml_data[base_name]['Connection']]['coordx'][comp_number]
        comp_y_val = locations[yaml_data[base_name]['Connection']]['coordy'][comp_number]

    # Make a reference stratigraphy object for the current component
    strata.append(sm.add_component_model_object(iam.Stratigraphy(
        name='strata' + component_name, parent=sm)))

    if var_type == 'LookupTable':
        lut_strata_data = yaml_data['Stratigraphy']['spatiallyVariable'][
            'LookupTableStratigraphy']
        file_name = lut_strata_data['FileName']
        file_directory = lut_strata_data['FileDirectory']
        if 'MaxPointDistance' in lut_strata_data:
            max_point_distance = float(lut_strata_data['MaxPointDistance'])
        else:
            warning_msg = ''.join([
                'MaxPointDistance was not present in the LookupTableStratigraphy ',
                'section of the .yaml file. This value sets how close a point ',
                'in the LookupStratigraphy table must be to a component for ',
                'the stratigraphy information to be used. The default value ',
                'of 100 m will be used.'])
            logging.debug(warning_msg)
            max_point_distance = 100.0

        LUTStrat_dict = get_lut_stratigraphy_dict(
            file_name, file_directory, comp_x_val, comp_y_val,
            max_point_distance=max_point_distance)

        numShaleLayers = LUTStrat_dict['numberOfShaleLayers']

        shaleThickness_List = LUTStrat_dict['shaleThickness_List']
        shaleThickness_min_List = LUTStrat_dict['shaleThickness_min_List']
        shaleThickness_max_List = LUTStrat_dict['shaleThickness_max_List']

        aquiferThickness_List = LUTStrat_dict['aquiferThickness_List']
        aquiferThickness_min_List = LUTStrat_dict['aquiferThickness_min_List']
        aquiferThickness_max_List = LUTStrat_dict['aquiferThickness_max_List']

        resThickness = LUTStrat_dict['resThickness']
        resThickness_min = LUTStrat_dict['resThickness']
        resThickness_max = LUTStrat_dict['resThickness']

        resDepth = LUTStrat_dict['depth']
        resDepth_min = LUTStrat_dict['depth_min']
        resDepth_max = LUTStrat_dict['depth_max']

        datumPressure = LUTStrat_dict['datumPressure']
        datumPressure_min = LUTStrat_dict['datumPressure_min']
        datumPressure_max = LUTStrat_dict['datumPressure_max']

        if numShaleLayers < strata[0].pars_bounds['numberOfShaleLayers'][0]:
            debug_msg = ''.join([
                'The number of shale layers went below the minimum value ',
                'allowed at x = {} m and y = {} m. Setting the parameter ',
                'to the minimum value.']).format(comp_x_val, comp_y_val)
            logging.debug(debug_msg)
            shaleThickness_List = strata[0].pars_bounds['numberOfShaleLayers'][0]

        if numShaleLayers > strata[0].pars_bounds['numberOfShaleLayers'][1]:
            debug_msg = ''.join([
                'The number of shale layers went above the maximum value ',
                'allowed at x = {} m and y = {} m. Setting the parameter ',
                'to the maximum value.']).format(comp_x_val, comp_y_val)
            logging.debug(debug_msg)
            numShaleLayers = strata[0].pars_bounds['numberOfShaleLayers'][1]

        strata[-1].add_par('numberOfShaleLayers', value=numShaleLayers,
                           vary=False)

        # Assign the unit thicknesses
        for shaleRef in range(numShaleLayers):
            shale_par_name = 'shale{}Thickness'.format(shaleRef + 1)
            if not shaleThickness_List[shaleRef]:
                shaleThickness_List[shaleRef] = strata[0].default_pars[shale_par_name].value

            if shaleThickness_List[shaleRef] < strata[0].pars_bounds['shaleThickness'][0]:
                debug_msg = layer_thickness_bound_debug_message(
                    'shale {}'.format(shaleRef+1), 'minimum',
                    comp_x_val, comp_y_val, msg_option=2)
                logging.debug(debug_msg)
                shaleThickness_List[shaleRef] = strata[0].pars_bounds['shaleThickness'][0]

            if shaleThickness_List[shaleRef] > strata[0].pars_bounds['shaleThickness'][1]:
                debug_msg = layer_thickness_bound_debug_message(
                    'shale {}'.format(shaleRef+1), 'maximum',
                    comp_x_val, comp_y_val, msg_option=2)
                logging.debug(debug_msg)
                shaleThickness_List[shaleRef] = strata[0].pars_bounds['shaleThickness'][1]

            if shaleThickness_min_List[shaleRef] and shaleThickness_max_List[shaleRef]:
                strata[-1].add_par(shale_par_name,
                                   value=shaleThickness_List[shaleRef],
                                   min=shaleThickness_min_List[shaleRef],
                                   max=shaleThickness_max_List[shaleRef])
            else:
                strata[-1].add_par(shale_par_name,
                                   value=shaleThickness_List[shaleRef],
                                   vary=False)

        for shaleRef in range(numShaleLayers-1):
            aq_par_name = 'aquifer{}Thickness'.format(shaleRef + 1)
            if not aquiferThickness_List[shaleRef]:
                aquiferThickness_List[shaleRef] = strata[0].default_pars[aq_par_name].value

            if aquiferThickness_List[shaleRef] < strata[0].pars_bounds['aquiferThickness'][0]:
                debug_msg = layer_thickness_bound_debug_message(
                    'aquifer {}'.format(shaleRef+1), 'minimum',
                    comp_x_val, comp_y_val, msg_option=2)
                logging.debug(debug_msg)
                shaleThickness_List[shaleRef] = strata[0].pars_bounds['aquiferThickness'][0]

            if aquiferThickness_List[shaleRef] > strata[0].pars_bounds['aquiferThickness'][1]:
                debug_msg = layer_thickness_bound_debug_message(
                    'aquifer {}'.format(shaleRef+1), 'maximum',
                    comp_x_val, comp_y_val, msg_option=2)
                logging.debug(debug_msg)
                shaleThickness_List[shaleRef] = strata[0].pars_bounds['aquiferThickness'][1]

            if aquiferThickness_min_List[shaleRef] and aquiferThickness_max_List[shaleRef]:
                strata[-1].add_par(aq_par_name,
                                   value=aquiferThickness_List[shaleRef],
                                   min=aquiferThickness_min_List[shaleRef],
                                   max=aquiferThickness_max_List[shaleRef])
            else:
                strata[-1].add_par(aq_par_name,
                                   value=aquiferThickness_List[shaleRef],
                                   vary=False)

        if not resThickness:
            resThickness = strata[0].default_pars['reservoirThickness'].value

        if resThickness:
            if resThickness_min and resThickness_max:
                strata[-1].add_par('reservoirThickness',
                                   value=resThickness,
                                   min=resThickness_min,
                                   max=resThickness_max)
            else:
                strata[-1].add_par('reservoirThickness',
                                   value=resThickness,
                                   vary=False)

        if resDepth:
            if resDepth != (sum(shaleThickness_List) + sum(aquiferThickness_List)):
                debug_msg = ''.join([
                    'In the file {}, the entry for depth at x = {} m and ',
                    'y = {} m was not equal to the sum of the unit thicknesses at ',
                    'that location.']).format(file_name, comp_x_val, comp_y_val)
                logging.debug(debug_msg)
                resDepth = (sum(shaleThickness_List) + sum(aquiferThickness_List))
        else:
            resDepth = (sum(shaleThickness_List) + sum(aquiferThickness_List))

        if resDepth_min and resDepth_max:
            strata[-1].add_par('depth', value=resDepth,
                               min=resDepth_min, max=resDepth_max)
        else:
            strata[-1].add_par('depth', value=resDepth, vary=False)

        if not datumPressure:
            datumPressure = strata[0].default_pars['datumPressure'].value

        if datumPressure_min and datumPressure_max:
            strata[-1].add_par('datumPressure', value=datumPressure,
                               min=datumPressure_min, max=datumPressure_max)
        else:
            strata[-1].add_par('datumPressure', value=datumPressure, vary=False)

    elif var_type == 'strikeAndDip':
        # depth_increase is positive when the depth increases
        depth_increase = depth_increases_from_strike_dip(
            strike, dip, dipDirection, comp_x_val, comp_y_val,
            coordxRefPoint=coordxRefPoint, coordyRefPoint=coordyRefPoint,
            output_type='single_point')

        strata_dict = get_strata_info_from_component(strata[0])

        numShaleLayers = strata_dict['numberOfShaleLayers']
        numShaleLayers_min = strata_dict['numberOfShaleLayers_min']
        numShaleLayers_max = strata_dict['numberOfShaleLayers_max']
        numShaleLayers_vary = strata_dict['numberOfShaleLayers_vary']

        reservoirThicknessRefPoint = strata_dict['reservoirThickness']
        reservoirThicknessRefPoint_min = strata_dict['reservoirThickness_min']
        reservoirThicknessRefPoint_max = strata_dict['reservoirThickness_max']
        reservoirThicknessRefPoint_vary = strata_dict['reservoirThickness_vary']

        resDepth = strata_dict['depth']
        resDepth_min = strata_dict['depth_min']
        resDepth_max = strata_dict['depth_max']
        resDepth_vary = strata_dict['depth_vary']
        resDepth_default = strata_dict['depth_default']

        datumPressure = strata_dict['datumPressure']
        datumPressure_min = strata_dict['datumPressure_min']
        datumPressure_max = strata_dict['datumPressure_max']
        datumPressure_vary = strata_dict['datumPressure_vary']

        shaleThicknessesRefPoint = strata_dict['shaleThicknesses']
        shaleThicknessesRefPoint_min = strata_dict['shaleThicknesses_min']
        shaleThicknessesRefPoint_max = strata_dict['shaleThicknesses_max']
        shaleThicknessesRefPoint_vary = strata_dict['shaleThicknesses_vary']

        aquiferThicknessesRefPoint = strata_dict['aquiferThicknesses']
        aquiferThicknessesRefPoint_min = strata_dict['aquiferThicknesses_min']
        aquiferThicknessesRefPoint_max = strata_dict['aquiferThicknesses_max']
        aquiferThicknessesRefPoint_vary = strata_dict['aquiferThicknesses_vary']

        # Now alter and assign the parameters to the stratigraphy component
        # created for the component's location. The depth_increase value is used to
        # make the top shale unit thicker or thinner, affecting the depths of all
        # other units.
        # Do not change the numberOfShaleLayers, however.
        if numShaleLayers_vary:
            strata[-1].add_par('numberOfShaleLayers', value=numShaleLayers,
                               min=numShaleLayers_min, max=numShaleLayers_max)
        else:
            strata[-1].add_par('numberOfShaleLayers', value=numShaleLayers, vary=False)

        if datumPressure:
            if datumPressure_vary:
                strata[-1].add_par('datumPressure', value=datumPressure,
                                   min=datumPressure_min, max=datumPressure_max)
            else:
                strata[-1].add_par('datumPressure', value=datumPressure, vary=False)

        # First, handle the shallowest shale
        shaleThickness = shaleThicknessesRefPoint[-1] + depth_increase
        shaleThickness_min = shaleThicknessesRefPoint_min[-1] + depth_increase
        shaleThickness_max = shaleThicknessesRefPoint_max[-1] + depth_increase

        # If the depth change has made the shaleThickness fall beneath the
        # minimum value, set it to the minimum value
        if shaleThickness < strata[0].pars_bounds['shaleThickness'][0]:
            debug_msg = layer_thickness_bound_debug_message(
                'shale {}'.format(numShaleLayers), 'minimum', comp_x_val, comp_y_val)
            logging.debug(debug_msg)
            extra_depth_increase = shaleThickness - strata[0].pars_bounds[
                'shaleThickness'][0]
            # Make the aquifer beneath this shale thinner to accomodate
            # the extra_depth_increase
            aquiferThicknessesRefPoint[-1] += extra_depth_increase
            aquiferThicknessesRefPoint_min[-1] += extra_depth_increase
            aquiferThicknessesRefPoint_max[-1] += extra_depth_increase
            del extra_depth_increase
            shaleThickness = strata[0].pars_bounds['shaleThickness'][0]

        # If the depth change has made the shaleThickness go above the maximum
        # value, set it to the maximum value
        if shaleThickness > strata[0].pars_bounds['shaleThickness'][1]:
            debug_msg = layer_thickness_bound_debug_message(
                'shale {}'.format(numShaleLayers), 'maximum', comp_x_val, comp_y_val)
            logging.debug(debug_msg)
            extra_depth_increase = shaleThickness - strata[0].pars_bounds[
                'shaleThickness'][1]
            # Make the aquifer beneath this shale thicker to accomodate
            # the extra_depth_increase
            aquiferThicknessesRefPoint[-1] += extra_depth_increase
            aquiferThicknessesRefPoint_min[-1] += extra_depth_increase
            aquiferThicknessesRefPoint_max[-1] += extra_depth_increase
            del extra_depth_increase
            shaleThickness = strata[0].pars_bounds['shaleThickness'][1]

        # Assign the parameter
        if shaleThicknessesRefPoint_vary[-1]:
            strata[-1].add_par('shale{}Thickness'.format(numShaleLayers),
                               value=shaleThickness,
                               min=shaleThickness_min,
                               max=shaleThickness_max)
        else:
            strata[-1].add_par('shale{}Thickness'.format(numShaleLayers),
                               value=shaleThickness, vary=False)

        if not resDepth_default:
            # I use this to double check that the depth obtained from altering resDepth
            # is the same as taking the sum of the unit thicknesses.
            resDepthVerify = 0
            resDepthVerify_min = 0
            resDepthVerify_max = 0

            resDepthVerify += shaleThickness
            resDepthVerify_min += shaleThickness_min
            resDepthVerify_max += shaleThickness_max

        # Now, go from the highest aquifer to the lowest shale so that required
        # thickness changes can be carried down through the units
        for shaleRef in range(numShaleLayers - 2, -1, -1):
            aquiferThickness = aquiferThicknessesRefPoint[shaleRef]
            aquiferThickness_min = aquiferThicknessesRefPoint_min[shaleRef]
            aquiferThickness_max = aquiferThicknessesRefPoint_max[shaleRef]

            shaleThickness = shaleThicknessesRefPoint[shaleRef]
            shaleThickness_min = shaleThicknessesRefPoint_min[shaleRef]
            shaleThickness_max = shaleThicknessesRefPoint_max[shaleRef]

            # If the depth change has made the aquiferThickness to fall beneath the
            # minimum value, set it to the minimum value
            if aquiferThickness < strata[0].pars_bounds['aquiferThickness'][0]:
                debug_msg = layer_thickness_bound_debug_message(
                    'aquifer {}'.format(shaleRef+1), 'minimum', comp_x_val, comp_y_val)
                logging.debug(debug_msg)
                extra_depth_increase = aquiferThickness \
                    - strata[0].pars_bounds['aquiferThickness'][0]
                # Make the shale beneath this aquifer thinner to accomodate
                # the extra_depth_increase
                shaleThickness += extra_depth_increase
                shaleThickness_min += extra_depth_increase
                shaleThickness_max += extra_depth_increase
                del extra_depth_increase
                aquiferThickness = strata[0].pars_bounds['aquiferThickness'][0]

            # If the depth change has made the aquiferThickness to go above the maximum
            # value, set it to the maximum value
            if aquiferThickness > strata[0].pars_bounds['aquiferThickness'][1]:
                debug_msg = layer_thickness_bound_debug_message(
                    'aquifer {}'.format(shaleRef+1), 'maximum', comp_x_val, comp_y_val)
                logging.debug(debug_msg)
                extra_depth_increase = aquiferThickness \
                    - strata[0].pars_bounds['aquiferThickness'][1]
                # Make the shale beneath this aquifer thicker to accomodate
                # the extra_depth_increase
                shaleThickness += extra_depth_increase
                shaleThickness_min += extra_depth_increase
                shaleThickness_max += extra_depth_increase
                del extra_depth_increase
                aquiferThickness = strata[0].pars_bounds['aquiferThickness'][1]

            # If the depth change has made the shaleThickness fall beneath the
            # minimum value, set it to the minimum value
            if shaleThickness < strata[0].pars_bounds['shaleThickness'][0]:
                debug_msg = layer_thickness_bound_debug_message(
                    'shale {}'.format(shaleRef+1), 'minimum', comp_x_val, comp_y_val)
                logging.debug(debug_msg)
                extra_depth_increase = shaleThickness \
                    - strata[0].pars_bounds['shaleThickness'][0]
                # Make the aquifer beneath this shale thinner to accomodate
                # the extra_depth_increase. If it's a shale 1, make the reservoir thinner.
                if (shaleRef + 1) != 1:
                    aquiferThicknessesRefPoint[shaleRef - 1] += extra_depth_increase
                    aquiferThicknessesRefPoint_min[shaleRef - 1] += extra_depth_increase
                    aquiferThicknessesRefPoint_max[shaleRef - 1] += extra_depth_increase
                elif (shaleRef + 1) == 1:
                    reservoirThicknessRefPoint += extra_depth_increase
                    if reservoirThicknessRefPoint_vary:
                        reservoirThicknessRefPoint_min += extra_depth_increase
                        reservoirThicknessRefPoint_max += extra_depth_increase
                del extra_depth_increase
                shaleThickness = strata[0].pars_bounds['shaleThickness'][0]

            # If the depth change has made the shaleThickness go above the maximum
            # value, set it to the maximum value
            if shaleThickness > strata[0].pars_bounds['shaleThickness'][1]:
                debug_msg = layer_thickness_bound_debug_message(
                    'shale {}'.format(shaleRef+1), 'maximum', comp_x_val, comp_y_val)
                logging.debug(debug_msg)
                extra_depth_increase = shaleThickness \
                    - strata[0].pars_bounds['shaleThickness'][1]
                # Make the aquifer beneath this shale thicker to accomodate
                # the extra_depth_increase. If it's shale 1, make the reservoir thicker.
                if (shaleRef + 1) != 1:
                    aquiferThicknessesRefPoint[shaleRef - 1] += extra_depth_increase
                    aquiferThicknessesRefPoint_min[shaleRef - 1] += extra_depth_increase
                    aquiferThicknessesRefPoint_max[shaleRef - 1] += extra_depth_increase
                elif (shaleRef + 1) == 1:
                    reservoirThicknessRefPoint += extra_depth_increase
                    if reservoirThicknessRefPoint_vary:
                        reservoirThicknessRefPoint_min += extra_depth_increase
                        reservoirThicknessRefPoint_max += extra_depth_increase
                del extra_depth_increase
                shaleThickness = strata[0].pars_bounds['shaleThickness'][1]

            if not resDepth_default:
                resDepthVerify += aquiferThickness
                resDepthVerify_min += aquiferThickness_min
                resDepthVerify_max += aquiferThickness_max

                resDepthVerify += shaleThickness
                resDepthVerify_min += shaleThickness_min
                resDepthVerify_max += shaleThickness_max

            # Assign the parameters
            if aquiferThicknessesRefPoint_vary[shaleRef]:
                strata[-1].add_par(
                    'aquifer{}Thickness'.format(shaleRef + 1), value=aquiferThickness,
                    min=aquiferThickness_min, max=aquiferThickness_max)
            else:
                strata[-1].add_par('aquifer{}Thickness'.format(shaleRef + 1),
                                   value=aquiferThickness, vary=False)

            if shaleThicknessesRefPoint_vary[shaleRef]:
                strata[-1].add_par(
                    'shale{}Thickness'.format(shaleRef + 1), value=shaleThickness,
                    min=shaleThickness_min, max=shaleThickness_max)
            else:
                strata[-1].add_par('shale{}Thickness'.format(shaleRef + 1),
                                   value=shaleThickness, vary=False)

        if reservoirThicknessRefPoint_vary:
            strata[-1].add_par(
                'reservoirThickness', value=reservoirThicknessRefPoint,
                min=reservoirThicknessRefPoint_min, max=reservoirThicknessRefPoint_max)
        else:
            strata[-1].add_par('reservoirThickness',
                               value=reservoirThicknessRefPoint, vary=False)

        if not resDepth_default:
            resDepth += depth_increase

            if resDepth < strata[0].pars_bounds['depth'][0]:
                resDepth = strata[0].pars_bounds['depth'][0]

            if resDepth > strata[0].pars_bounds['depth'][1]:
                resDepth = strata[0].pars_bounds['depth'][1]

            if resDepth != resDepthVerify:
                err_msg = ''.join([
                    'The reservoir top depth obtained by altering the ',
                    'reservoir top depth at the reference point did not equal ',
                    'the reservoir top depth obtained through a summation ',
                    'of unit thicknesses at x = {} m and y = {} m. ',
                    CHECK_CONDITIONS_MSG]).format(comp_x_val, comp_y_val)
                raise KeyError(err_msg)

            if resDepth_vary:
                resDepth_min += depth_increase
                resDepth_max += depth_increase

                if resDepth_min < strata[0].pars_bounds['depth'][0]:
                    resDepth_min = strata[0].pars_bounds['depth'][0]

                if resDepth_max > strata[0].pars_bounds['depth'][1]:
                    resDepth_max = strata[0].pars_bounds['depth'][1]

                if resDepth_min != resDepthVerify_min:
                    err_msg = ''.join([
                        'The minimum reservoir top depth obtained by altering the ',
                        'minimum reservoir top depth at the reference point ',
                        'did not equal the minimum reservoir top depth obtained '
                        'through a summation of unit thicknesses at ',
                        'x = {} m and y = {} m. ', CHECK_CONDITIONS_MSG]).format(
                            comp_x_val, comp_y_val)
                    raise KeyError(err_msg)

                if resDepth_max != resDepthVerify_max:
                    err_msg = ''.join([
                        'The maximum reservoir top depth obtained by altering the ',
                        'maximum reservoir top depth at the reference point ',
                        'did not equal the maximum reservoir top depth obtained '
                        'through a summation of unit thicknesses at ',
                        'x = {} m and y = {} m. ', CHECK_CONDITIONS_MSG]).format(
                            comp_x_val, comp_y_val)
                    raise KeyError(err_msg)

                strata[-1].add_par('depth', value=resDepth,
                                   min=resDepth_min, max=resDepth_max)

            else:
                strata[-1].add_par('depth', value=resDepth, vary=False)

    return strata, sm


def depth_increases_from_strike_dip(strike, dip, dipDirection,
                                    x_locations, y_locations,
                                    output_type='single_point',
                                    coordxRefPoint=0, coordyRefPoint=0):
    """
    This function calculates the changes in unit depths due to a prescribed
    strike and dip.

    :param strike: Unit strike in degrees, where 0 is north, 90 is east, 180 is
        south, and 270 is west.
    :type strike: int or float

    :param dip: Unit dip in degrees, where positive means dipping down in
        the dipDirection provided
    :type dip: int or float

    :param dipDirection: Direction of the units dip, expressed as N, E, S, W,
        NE, SE, SW, or SW. The exact direction of dip in degrees 'clockwise'
        from North is calculated from the strike and dipDirection.
    :type dipDirection: str

    :param x_locations: x value or a list of x values for the location(s) at which
        the function estimates changes in depth due to strike and dip. Note
        that x is assumed to increase to the east.
    :type x_locations: int, float, or list

    :param y_locations: y value or a list of y values for the location(s) at which
        the function estimates changes in depth due to strike and dip. Note
        that y is assumed to increase to the north.
    :type x_locations: int, float, or list

    :param output_type: Option specifying the type of analysis. The first
        option is 'single_point', where the function is given the x and y
        coordinates for a single location and the function then provides the
        change in unit depth at that location (due to the strike and dip data
        provided). The second option is 'point_list', where the function is
        given lists of x and y values for multiple locations. The lists are
        assumed to be aligned, so that the first x value and first y value
        represent one location (etc.). With 'point_list', the function returns
        the changes in elevation at each of the locations provided. The third
        option is 'grid'. While 'point_list' produces a 1-dimensional list of
        depth changes for a 1-dimensional list of x and y coordinates, the grid
        option produces a 2-dimensional array of depth changes for 2-D grid of
        x and y points. The x and y grid can be made with the 'grid' option in
        get_known_locations() within file locations.py in the folder cfi.
    :type output_type: str

    :param coordxRefPoint: x value for the reference point. The reference point
        is assumed to have a change in depth value of zero, and it is used
        in the calculation of depth changes (i.e., using three points to make
        a plane).
    :type coordxRefPoint: int or float

    :param coordyRefPoint: y value for the reference point. The reference point
        is assumed to have a change in depth value of zero, and it is used
        in the calculation of depth changes (i.e., using three points to make
        a plane).
    :type coordyRefPoint: int or float

    :returns: depth_increase
    """
    if isinstance(x_locations, list) and isinstance(y_locations, list):
        if len(x_locations) != len(y_locations):
            err_msg = ''.join([
                'The x_locations and y_locations lists provided to the ',
                'function depth_increases_from_strike_dip() in strata.py ',
                'do not have equal lengths.'])
            raise ValueError(err_msg)

    if isinstance(x_locations, list) and not isinstance(y_locations, list):
        err_msg = ''.join([
            'The x_locations provided to the function depth_increases_from_strike_dip() ',
            'in strata.py were a list, but the y_locations were not a list.',
            'Check your input.'])
        raise TypeError(err_msg)

    if not isinstance(x_locations, list) and isinstance(y_locations, list):
        err_msg = ''.join([
            'The y_locations provided to the function depth_increases_from_strike_dip() ',
            'in strata.py were a list, but the x_locations were not a list.',
            'Check your input.'])
        raise TypeError(err_msg)

    if strike > 360:
        err_msg = ''.join(['The strike provided was > 360 degrees, ',
                           'check the input values.'])
        raise ValueError(err_msg)

    if strike < 0:
        err_msg = ''.join(['The strike provided was < 0, but a value ',
                           'between 0 and 360 degrees was expected.'])
        raise ValueError(err_msg)

    if dip == 0:
        dip = np.abs(dip)
        err_msg = ''.join(['The dip provided was 0 degrees, in which case ',
                           ' calculation of depth changes due to dip should ',
                           'not be performed.'])
        raise ValueError(err_msg)

    if dip < 0:
        dip = np.abs(dip)
        err_msg = ''.join(['The dip provided was < 0 but a positive ',
                           'value was expected.'])
        raise ValueError(err_msg)

    if dip >= 90:
        dip = np.abs(dip)
        err_msg = ''.join(['The dip provided was >= 90 degrees, but a ',
                           'value between 0 and 90 was expected.'])
        raise ValueError(err_msg)

    dipDirectionDegrees = obtain_dip_direction_degrees(strike, dipDirection)

    # Create a 3-dimensional plane representing the changes in the depth of
    # each unit's base. To create the plane, use the change in depths for 3
    # locations. One location is the reference point (x = coordxRefPoint, y =
    # coordyRefPoint, and z = 0). The two other locations will have z values
    # calculated with the strike and dip. Point p1 is in the dip direction,
    # while point p2 is 90 degrees counter clockwise from the dip direction
    # (in one of the directions of strike). Points 3 and 4 aren't needed here.

    x_points, y_points, z_points = make_xyz_points_from_strike_and_dip(
        dip, dipDirectionDegrees, L1=100, L2=100, L3=100,
        L4=100, point0_xyz=[coordxRefPoint, coordyRefPoint, 0])

    x_p0 = x_points[0]
    x_p1 = x_points[1]
    x_p2 = x_points[2]

    y_p0 = y_points[0]
    y_p1 = y_points[1]
    y_p2 = y_points[2]

    z_p0 = z_points[0]
    z_p1 = z_points[1]
    z_p2 = z_points[2]

    # Set of known x, y, and z points
    points = [[x_p0, y_p0, z_p0],
              [x_p1, y_p1, z_p1],
              [x_p2, y_p2, z_p2]]

    p0, p1, p2 = points
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    x2, y2, z2 = p2

    ux, uy, uz = [x1 - x0, y1 - y0, z1 - z0]
    vx, vy, vz = [x2 - x0, y2 - y0, z2 - z0]

    u_cross_v = [(uy * vz) - (uz * vy),
                 (uz * vx) - (ux * vz),
                 (ux * vy) - (uy * vx)]

    point  = np.array(p0)
    normal = np.array(u_cross_v)

    d = -point.dot(normal)

    # Solve for the change in depth at the location - depth_increase is positive
    # when the depth increases and is negative when it decreases.
    if output_type == 'single_point':
        depth_increase = -(-(normal[0] * x_locations)
                        - (normal[1] * y_locations) - d) / normal[2]

    elif output_type == 'point_list':
        depth_increase = []

        for loc_ref, xval in enumerate(x_locations):
            depth_increase_val = -(-(normal[0] * xval)
                     - (normal[1] * y_locations[loc_ref]) - d) / normal[2]
            depth_increase.append(depth_increase_val)

    elif output_type == 'grid':
        depth_increase = np.zeros((x_locations.shape[0], x_locations.shape[1]))

        for x_ref in range(x_locations.shape[0]):
            for y_ref in range(x_locations.shape[1]):
                depth_increase[x_ref, y_ref] = -(-(normal[0] * x_locations[x_ref, y_ref])
                         - (normal[1] * y_locations[x_ref, y_ref]) - d) / normal[2]

    return depth_increase


def layer_thickness_bound_debug_message(unit, bound, locx, locy, msg_option=1):
    """
    Returns string delivering debug message regarding setting unit thickness to
    either minimum or maximum value.

    unit: string, 'shale' or 'aquifer'
    index: integer, index of the unit: e.g., 2 for aquifer 2
    bound: string, 'minimum' or 'maximum'
    locx: float, x-cordinate of location of interest
    locy: float, y-cordinate of location of interest

    """
    side = {'minimum': 'below', 'maximum': 'above'}

    if msg_option == 1:
        msg = ''.join([
            'The thickness of {} went {} the {} value allowed at ',
            'x = {} m and y = {} m. Setting the parameter to the {} value ',
            'and accomplishing the required depth changes by altering ',
            'the thicknesses of the underlying units.']).format(
                unit, side[bound], bound, locx, locy, bound)
    else:
        msg = ''.join([
            'The thickness of {} went {} the {} value allowed at ',
            'x = {} m and y = {} m. Setting the parameter to the {} value.']).format(
                unit, side[bound], bound, locx, locy, bound)

    return msg


def update_stratigraphy_by_strike_and_dip(numberOfShaleLayers, shaleThicknessList,
                                          aquiferThicknessList, reservoirThickness,
                                          strike, dip, dipDirection,
                                          coordxRefPoint, coordyRefPoint,
                                          location_x, location_y, strataRefPoint):
    """
    This function provides an updated thicknesses for shales, aquifers, and the
    reservoir based on strike and dip information as well as the location at
    which the estimates are to be made. If you are using this function within a
    loop (e.g., see stratigraphy_plot.py), be sure to provide shaleThicknessList
    and aquiferThicknessList to the function as:
        shaleThicknessList = shaleThicknessList[:]
        aquiferThicknessList = aquiferThicknessList[:]
    Otherwise, the repeated use of the function within the loop may allow the
    changes to the lists to carry over, potentially resulting in negative unit
    thicknesses.

    :param numberOfShaleLayers: number of shale layers
    :type numberOfShaleLayers: int

    :param shaleThicknessList: List of the thickness for each shale layer, with
        the first entry representing the deepest shale and the last entry
        representing the shallowest shale
    :type shaleThicknessList: list

    :param aquiferThicknessList: List of the thickness [|m|] for each aquifer,
        with the first entry representing the deepest aquifer and the last entry
        representing the shallowest aquifer. These thicknesses are those at the
        reference point located at (coordxRefPoint, coordyRefPoint).
    :type aquiferThicknessList: list

    :param reservoirThickness: thickness of the reservoir [|m|] at the
        reference point.
    :type reservoirThickness: int or float

    :param strike: Unit strike in degrees, where 0 is north, 90 is east, 180 is
        south, and 270 is west.
    :type strike: int or float

    :param dip: Unit dip in degrees, where positive means dipping down in
        the dipDirection provided
    :type dip: int or float

    :param dipDirection: Direction of the units dip, expressed as N, E, S, W,
        NE, SE, SW, or SW. Note that dipDirection must be compatible with the
        strike provided.
    :type dipDirection: str

    :param coordxRefPoint: x value for the reference point. The reference point
        is assumed to have a change in depth value of zero, and it is used
        in the calculation of depth changes (i.e., using three points to make
        a plane).
    :type coordxRefPoint: int or float

    :param coordyRefPoint: y value for the reference point. The reference point
        is assumed to have a change in depth value of zero, and it is used
        in the calculation of depth changes (i.e., using three points to make
        a plane).
    :type coordyRefPoint: int or float

    :param location_x: x value [|m|] for the location at which the function is
        estimating unit thicknesses.
    :type location_x: int or float

    :param location_y: y value [|m|] for the location at which the function is
        estimating unit thicknesses.
    :type location_y: int or float

    :param strataRefPoint: stratigraphy component already created for the
        reference point. This component is used to check the parameter bounds
        for thickness values.
    :type strataRefPoint: openiam.stratigraphy_component.Stratigraphy

    :returns: updatedStratigraphy, a dictionary containing the updated
        thicknesses and top depths for all shales, aquifers, and the reservoir:

        updatedStratigraphy['shaleThicknessList'] = shaleThicknessListUpdated
        updatedStratigraphy['aquiferThicknessList'] = aquiferThicknessListUpdated
        updatedStratigraphy['reservoirThickness'] = reservoirThicknessUpdated
        updatedStratigraphy['shaleTopDepthList'] = shaleTopDepthListUpdated
        updatedStratigraphy['aquiferTopDepthList'] = aquiferTopDepthListUpdated
        updatedStratigraphy['reservoirTopDepth'] = resTopDepthUpdated
        updatedStratigraphy['shaleBottomDepthList'] = shaleBottomDepthListUpdated
        updatedStratigraphy['aquiferBottomDepthList'] = aquiferBottomDepthListUpdated
        updatedStratigraphy['reservoirBottomDepth'] = resBottomDepthUpdated

        Each of the lists within the dictionary updatedStratigraphy are listed
        in ascending order, so the first entry represents the deepest shale or
        aquifer while the last entry represents the shallowest. Note that all
        top depths represent the top of each unit, while the bottom depths
        represent the bottom of each unit.

    """
    depth_increase = depth_increases_from_strike_dip(strike, dip, dipDirection,
                                                     location_x, location_y,
                                                     output_type='single_point',
                                                     coordxRefPoint=coordxRefPoint,
                                                     coordyRefPoint=coordyRefPoint)

    reservoirThicknessUpdated = reservoirThickness

    shaleThicknessListUpdated = [None] * numberOfShaleLayers
    aquiferThicknessListUpdated = [None] * (numberOfShaleLayers - 1)

    # First, handle the highest shale
    shaleThickness = shaleThicknessList[-1] + depth_increase

    # If the depth change has made the shaleThickness fall beneath the
    # minimum value, set it to the minimum value
    if shaleThickness < strataRefPoint.pars_bounds['shaleThickness'][0]:
        debug_msg = layer_thickness_bound_debug_message(
            'shale {}'.format(numberOfShaleLayers), 'minimum', location_x, location_y)
        logging.debug(debug_msg)
        additional_depth_increase = shaleThickness - strataRefPoint.pars_bounds[
            'shaleThickness'][0]
        # Make the aquifer beneath this shale thinner to accomodate
        # the additional_depth_increase
        aquiferThicknessList[-1] += additional_depth_increase
        del additional_depth_increase
        shaleThickness = strataRefPoint.pars_bounds['shaleThickness'][0]

    # If the depth change has made the shaleThickness go above the maximum
    # value, set it to the maximum value
    if shaleThickness > strataRefPoint.pars_bounds['shaleThickness'][1]:
        debug_msg = layer_thickness_bound_debug_message(
            'shale {}'.format(numberOfShaleLayers), 'maximum', location_x, location_y)
        logging.debug(debug_msg)
        additional_depth_increase = shaleThickness - strataRefPoint.pars_bounds[
            'shaleThickness'][1]
        # Make the aquifer beneath this shale thicker to accomodate
        # the additional_depth_increase
        aquiferThicknessList[-1] += additional_depth_increase
        del additional_depth_increase
        shaleThickness = strataRefPoint.pars_bounds['shaleThickness'][1]

    shaleThicknessListUpdated[-1] = shaleThickness

    # Now, go from the highest aquifer to the lowest shale so that required
    # thickness changes can be carried down through the units
    for shaleRef in range(numberOfShaleLayers - 2, -1, -1):
        aquiferThickness = aquiferThicknessList[shaleRef]
        shaleThickness = shaleThicknessList[shaleRef]

        # If the depth change has made the aquiferThickness fall beneath the
        # minimum value, set it to the minimum value
        if aquiferThickness < strataRefPoint.pars_bounds['aquiferThickness'][0]:
            debug_msg = layer_thickness_bound_debug_message(
                'aquifer {}'.format(shaleRef+1), 'minimum', location_x, location_y)
            logging.debug(debug_msg)
            additional_depth_increase = aquiferThickness \
                - strataRefPoint.pars_bounds['aquiferThickness'][0]
            # Make the shale beneath this aquifer thinner to accomodate
            # the additional_depth_increase
            shaleThickness += additional_depth_increase
            del additional_depth_increase
            aquiferThickness = strataRefPoint.pars_bounds['aquiferThickness'][0]

        # If the depth change has made the aquiferThickness go above the maximum
        # value, set it to the maximum value
        if aquiferThickness > strataRefPoint.pars_bounds['aquiferThickness'][1]:
            debug_msg = layer_thickness_bound_debug_message(
                'aquifer {}'.format(shaleRef+1), 'maximum', location_x, location_y)
            logging.debug(debug_msg)
            additional_depth_increase = aquiferThickness \
                - strataRefPoint.pars_bounds['aquiferThickness'][1]
            # Make the shale beneath this aquifer thicker to accomodate
            # the additional_depth_increase
            shaleThickness += additional_depth_increase
            del additional_depth_increase
            aquiferThickness = strataRefPoint.pars_bounds['aquiferThickness'][1]

        aquiferThicknessListUpdated[shaleRef] = aquiferThickness

        # If the depth change has made the shaleThickness fall beneath the
        # minimum value, set it to the minimum value
        if shaleThickness < strataRefPoint.pars_bounds['shaleThickness'][0]:
            debug_msg = layer_thickness_bound_debug_message(
                'shale {}'.format(shaleRef+1), 'minimum', location_x, location_y)
            logging.debug(debug_msg)
            additional_depth_increase = shaleThickness \
                - strataRefPoint.pars_bounds['shaleThickness'][0]
            # Make the aquifer beneath this shale thinner to accomodate
            # the additional_depth_increase. If it's shale 1, make the reservoir thinner.
            if (shaleRef + 1) != 1:
                aquiferThicknessList[shaleRef - 1] += additional_depth_increase
            elif (shaleRef + 1) == 1:
                reservoirThicknessUpdated += additional_depth_increase
            del additional_depth_increase
            shaleThickness = strataRefPoint.pars_bounds['shaleThickness'][0]

        # If the depth change has made the shaleThickness go above the maximum
        # value, set it to the maximum value
        if shaleThickness > strataRefPoint.pars_bounds['shaleThickness'][1]:
            debug_msg = layer_thickness_bound_debug_message(
                'shale {}'.format(shaleRef+1), 'maximum', location_x, location_y)
            logging.debug(debug_msg)
            additional_depth_increase = shaleThickness \
                - strataRefPoint.pars_bounds['shaleThickness'][1]
            # Make the aquifer beneath this shale thicker to accomodate
            # the additional_depth_increase. If it's shale 1, make the reservoir thicker.
            if (shaleRef + 1) != 1:
                aquiferThicknessList[shaleRef - 1] += additional_depth_increase
            elif (shaleRef + 1) == 1:
                reservoirThicknessUpdated += additional_depth_increase
            del additional_depth_increase
            shaleThickness = strataRefPoint.pars_bounds['shaleThickness'][1]

        shaleThicknessListUpdated[shaleRef] = shaleThickness
    # End of loop through shales and aquifers
    del shaleThicknessList, aquiferThicknessList

    if reservoirThicknessUpdated < strataRefPoint.pars_bounds['reservoirThickness'][0]:
        debug_msg = layer_thickness_bound_debug_message(
            'reservoir', 'minimum', location_x, location_y, msg_option=2)
        logging.debug(debug_msg)
        reservoirThicknessUpdated = strataRefPoint.pars_bounds[
            'reservoirThickness'][0]

    if reservoirThicknessUpdated > strataRefPoint.pars_bounds['reservoirThickness'][1]:
        debug_msg = layer_thickness_bound_debug_message(
            'reservoir', 'maximum', location_x, location_y, msg_option=2)
        logging.debug(debug_msg)
        reservoirThicknessUpdated = strataRefPoint.pars_bounds[
            'reservoirThickness'][1]

    resTopDepthUpdated = np.sum([shaleThicknessListUpdated])
    resTopDepthUpdated += np.sum([aquiferThicknessListUpdated])

    resBottomDepthUpdated = resTopDepthUpdated + reservoirThicknessUpdated

    if resTopDepthUpdated < strataRefPoint.pars_bounds['depth'][0]:
        err_msg = ''.join([
            'The depth of the reservoir top went below the minimum value ',
            'allowed at x = {} m and y = {} m. ',
            CHECK_CONDITIONS_MSG]).format(location_x, location_y)
        raise ValueError(err_msg)

    if resTopDepthUpdated > strataRefPoint.pars_bounds['depth'][1]:
        err_msg = ''.join([
            'The depth of the reservoir top went above the maximum value ',
            'allowed at x = {} m and y = {} m. ',
            CHECK_CONDITIONS_MSG]).format(location_x, location_y)
        raise ValueError(err_msg)

    shaleTopDepthListUpdated = [None] * numberOfShaleLayers
    aquiferTopDepthListUpdated = [None] * (numberOfShaleLayers - 1)

    shaleBottomDepthListUpdated = [None] * numberOfShaleLayers
    aquiferBottomDepthListUpdated = [None] * (numberOfShaleLayers - 1)

    # The top depth of the highest shale is 0
    shaleTopDepthListUpdated[-1] = 0
    shaleBottomDepthListUpdated[-1] = shaleThicknessListUpdated[-1]

    for shaleRef in range(numberOfShaleLayers - 2, -1, -1):
        aquiferTopDepthListUpdated[shaleRef] = shaleTopDepthListUpdated[
            shaleRef + 1] + shaleThicknessListUpdated[shaleRef + 1]

        if aquiferTopDepthListUpdated[shaleRef] <= 0:
            err_msg = ''.join([
                'The top depth of aquifer {} was <= 0 m at x = {} m and ',
                'y = {} m. Only shale {} should reach the surface. ',
                CHECK_CONDITIONS_MSG]).format(
                    shaleRef + 1, location_x, location_y, numberOfShaleLayers)
            raise ValueError(err_msg)

        if aquiferTopDepthListUpdated[shaleRef] <= shaleTopDepthListUpdated[shaleRef + 1]:
            err_msg = ''.join([
                'The top depth of aquifer {} was <= the top depth of shale {} ',
                'at x = {} m and y = {} m. This aquifer should always be ',
                'beneath this shale. ', CHECK_CONDITIONS_MSG]).format(
                    shaleRef + 1, shaleRef + 2, location_x, location_y)
            raise ValueError(err_msg)

        aquiferBottomDepthListUpdated[shaleRef] = aquiferTopDepthListUpdated[
            shaleRef] + aquiferThicknessListUpdated[shaleRef]

        shaleTopDepthListUpdated[shaleRef] = aquiferTopDepthListUpdated[
            shaleRef] + aquiferThicknessListUpdated[shaleRef]

        if shaleTopDepthListUpdated[shaleRef] <= 0:
            err_msg = ''.join([
                'The top depth of shale {} was <= 0 m at x = {} m and ',
                'y = {} m. Only shale {} should reach the surface. ',
                CHECK_CONDITIONS_MSG]).format(
                    shaleRef + 1, location_x, location_y, numberOfShaleLayers)
            raise ValueError(err_msg)

        if shaleTopDepthListUpdated[shaleRef] <= aquiferTopDepthListUpdated[shaleRef]:
            err_msg = ''.join([
                'The top depth of shale {} was <= the top depth of aquifer {} ',
                'at x = {} m and y = {} m. This shale should always be beneath this ',
                'aquifer. ', CHECK_CONDITIONS_MSG]).format(
                    shaleRef + 1, shaleRef + 1, location_x, location_y)
            raise ValueError(err_msg)

        shaleBottomDepthListUpdated[shaleRef] = shaleTopDepthListUpdated[
            shaleRef] + shaleThicknessListUpdated[shaleRef]

    updatedStratigraphy = dict()

    updatedStratigraphy['shaleThicknessList'] = shaleThicknessListUpdated
    updatedStratigraphy['aquiferThicknessList'] = aquiferThicknessListUpdated
    updatedStratigraphy['reservoirThickness'] = reservoirThicknessUpdated

    updatedStratigraphy['shaleTopDepthList'] = shaleTopDepthListUpdated
    updatedStratigraphy['aquiferTopDepthList'] = aquiferTopDepthListUpdated
    updatedStratigraphy['reservoirTopDepth'] = resTopDepthUpdated

    updatedStratigraphy['shaleBottomDepthList'] = shaleBottomDepthListUpdated
    updatedStratigraphy['aquiferBottomDepthList'] = aquiferBottomDepthListUpdated
    updatedStratigraphy['reservoirBottomDepth'] = resBottomDepthUpdated

    return updatedStratigraphy


def obtain_dip_direction_degrees(strike, dipDirection):
    """
    Function that provides the dip direction in degrees clockwise from north,
    so that 90 is east, 180 is south, and 270 is west.

    :param strike: strike of the units in degrees clockwise from north, where
        0 is north, 90 is east, 180 is south, and 270 is west.
    :type strike: int or float

    :param dipDirection: Direction of dip, provided in a cardinal direction -
        N, E, S, W, NE, SE, SW, or NW.
    :type dipDirection: str

    :returns: dipDirectionDegrees
    """
    if dipDirection not in DIP_DIRECTION_OPTIONS:
        err_msg = ''.join(['Dip direction provided is not E, S, ',
                           'W, N, NE, SE, SW, or NW.'])
        raise KeyError(err_msg)

    dipDirectionDegrees = None
    if 0 < strike < 90:
        if dipDirection in ['S', 'E', 'SE']:
            dipDirectionDegrees = strike + 90
        elif dipDirection in ['N', 'W', 'NW']:
            dipDirectionDegrees = strike + 270
    elif 90 < strike < 180:
        if dipDirection in ['N', 'E', 'NE']:
            dipDirectionDegrees = strike - 90
        elif dipDirection in ['S', 'W', 'SW']:
            dipDirectionDegrees = strike + 90
    elif 180 < strike < 270:
        if dipDirection in ['N', 'W', 'NW']:
            dipDirectionDegrees = strike + 90
        elif dipDirection in ['S', 'E', 'SE']:
            dipDirectionDegrees = strike - 90
    elif 270 < strike < 360:
        if dipDirection in ['N', 'E', 'NE']:
            dipDirectionDegrees = (strike + 90) - 360
        elif dipDirection in ['S', 'W', 'SW']:
            dipDirectionDegrees = strike - 90
    elif strike in [0, 180, 360]:
        if dipDirection == 'E':
            dipDirectionDegrees = 90
        elif dipDirection == 'W':
            dipDirectionDegrees = 270
    elif strike in [90, 270]:
        if dipDirection == 'N':
            dipDirectionDegrees = 0
        elif dipDirection == 'S':
            dipDirectionDegrees = 180

    if dipDirectionDegrees is None:
        err_msg = ''.join(['The dip direction provided does not work ',
                           'with the strike provided. Dip must be ',
                           'orthogonal to strike. For example, the strike ',
                           'and dip cannot both be to the NE. A unit '
                           'striking to the NE can have dips of E, S, ',
                           'or SE (all for dip to the SE) or N, W, or NW',
                           '(all for dip to the NW). A unit striking ',
                           'N or S (0 or 180 degrees) can only dip E or W. ',
                           'A unit striking E or W (90 or 270 degrees) ',
                           'can only dip N or S.'])
        raise ValueError(err_msg)

    return dipDirectionDegrees


def get_strata_var_info_from_yaml(yaml_data):
    """
    Function that obtains strike and dip information from a .yaml file.

    :param yaml_data: The dictionary created from an input .yaml_file
    :type yaml_data: dict

    :returns: strata_var_info
    """
    # Setup default values of output dictionary
    strata_var_info = {'var_type': 'noVariation',
                       'strike': None,
                       'dip': None,
                       'dipDirection': None,
                       'coordxReferencePoint': None,
                       'coordyReferencePoint': None}

    if 'spatiallyVariable' in yaml_data['Stratigraphy']:
        data = yaml_data['Stratigraphy']['spatiallyVariable']
        if 'strikeAndDip' in data:
            strata_var_info['var_type'] = 'strikeAndDip'

            try:
                for key in ['strike', 'dip', 'dipDirection']:
                    strata_var_info[key] = data['strikeAndDip'][key]
            except KeyError:
                err_msg = ''.join(['Strike, dip, and/or dipDirection inputs were ',
                                   'not provided correctly in the .yaml file.'])
                raise KeyError(err_msg) from None

            if 'coordxRefPoint' in data['strikeAndDip']:
                strata_var_info['coordxReferencePoint'] = data['strikeAndDip']['coordxRefPoint']
            else:
                debug_msg = ''.join(['No coordxReferencePoint value was provided in the ',
                                     '.yaml file. Using the default value of ',
                                     'coordxReferencePoint = 0.'])
                logging.debug(debug_msg)
                strata_var_info['coordxReferencePoint'] = 0

            if 'coordyRefPoint' in data['strikeAndDip']:
                strata_var_info['coordyReferencePoint'] = data['strikeAndDip']['coordyRefPoint']
            else:
                debug_msg = ''.join(['No coordyReferencePoint value was provided in the ',
                                     '.yaml file. Using the default value of ',
                                     'coordyReferencePoint = 0.'])
                logging.debug(debug_msg)
                strata_var_info['coordyReferencePoint'] = 0

        elif 'LookupTableStratigraphy' in data:
            strata_var_info['var_type'] = 'LookupTable'

    return strata_var_info


def get_strata_info_from_component(strata_comp):
    """
    Function that obtains stratigraphy information from an already created
    stratigraphy component. The information obtained includes the number of
    shale layers, the thickness of each shale layer, the thickness of each
    aquifer, and the thickness of the reservoir. Minimum and maximum values for
    each value are also obtained, although those values are zero if the
    stratigraphy component does not have minimum and mazimum values.

    Note that this function deals with the input provided directly to
    stratigraphy components, while the function get_unit_depth_from_component()
    provides the unit depths (top or bottom) that are derived from the
    simulation input.

    :param strata_comp: Stratigraphy component used to obtain information
    :type strata_comp: instance of Stratigraphy component class

    """
    # Initialize output dictionary
    strata_dict = dict()
    # Defaults for parameter numberOfShaleLayers
    strata_dict['numberOfShaleLayers_min'] = 0
    strata_dict['numberOfShaleLayers_max'] = 0
    strata_dict['numberOfShaleLayers_vary'] = False

    # The min, max, and numShaleLayers_vary are not used right now, but they
    # are kept here because they might be used in future versions.
    numShaleLayers = iam.cfi.commons.get_parameter_val(
        strata_comp, 'numberOfShaleLayers')

    if 'numberOfShaleLayers' in strata_comp.pars:
        strata_dict['numberOfShaleLayers_min'] = \
            strata_comp.pars['numberOfShaleLayers'].min
        strata_dict['numberOfShaleLayers_max'] = \
            strata_comp.pars['numberOfShaleLayers'].max
        strata_dict['numberOfShaleLayers_vary'] = True

    strata_dict['numberOfShaleLayers'] = numShaleLayers

    # Initialize additional keys
    for key in ['shaleThicknesses', 'aquiferThicknesses']:
        strata_dict[key] = []

    for key in ['shaleThicknesses_min', 'shaleThicknesses_max',
                'aquiferThicknesses_min', 'aquiferThicknesses_max']:
        strata_dict[key] = numShaleLayers*[0]

    for key in ['shaleThicknesses_vary', 'aquiferThicknesses_vary']:
        strata_dict[key] = numShaleLayers*[False]

    # Get the unit thicknesses
    # Shale layers
    for shaleRef in range(numShaleLayers):
        shale_par_nm = 'shale{}Thickness'.format(shaleRef + 1)

        strata_dict['shaleThicknesses'].append(
            iam.cfi.commons.get_parameter_val(strata_comp, shale_par_nm))

        # The min and max are not used currently, but they are kept for potential
        # updates in the future.
        if shale_par_nm in strata_comp.pars:
            strata_dict['shaleThicknesses_min'][shaleRef] = \
                strata_comp.pars[shale_par_nm].min
            strata_dict['shaleThicknesses_max'][shaleRef] = \
                strata_comp.pars[shale_par_nm].max
            strata_dict['shaleThicknesses_vary'][shaleRef] = True

    # Aquifer layers
    for shaleRef in range(numShaleLayers-1):
        aq_par_nm = 'aquifer{}Thickness'.format(shaleRef + 1)

        strata_dict['aquiferThicknesses'].append(
            iam.cfi.commons.get_parameter_val(strata_comp, aq_par_nm))

        if aq_par_nm in strata_comp.pars:
            strata_dict['aquiferThicknesses_min'][shaleRef] = \
                strata_comp.pars[aq_par_nm].min
            strata_dict['aquiferThicknesses_max'][shaleRef] = \
                strata_comp.pars[aq_par_nm].max
            strata_dict['aquiferThicknesses_vary'][shaleRef] = True

    # ReservoirThickness isn't needed for this analysis, but it is needed in the
    # update_stratigraphy() function that is called.
    # Add additional keys
    strata_dict['depth_default'] = False
    for par_name in ['reservoirThickness', 'depth', 'datumPressure']:
        for key in ['_min', '_max']:
            # Default values
            strata_dict[par_name+key] = 0  # e.g. 'reservoirThickness_min'

        strata_dict[par_name+'_vary'] = False   # e.g. 'reservoirThickness_vary'

        strata_dict[par_name] = iam.cfi.commons.get_parameter_val(
            strata_comp, par_name)

        # The min and max vlaues are not used currently, but they are kept for
        # potential updates in the future.
        if par_name in strata_comp.pars:
            strata_dict[par_name+'_min'] = strata_comp.pars[par_name].min
            strata_dict[par_name+'_max'] = strata_comp.pars[par_name].max
            strata_dict[par_name+'_vary'] = True
        elif par_name in strata_comp.default_pars and par_name == 'depth':
            strata_dict['depth_default'] = True

    return strata_dict


def make_xyz_points_from_strike_and_dip(dip, dipDirectionDegrees, L1=100,
                                        L2=100, L3=100, L4=100,
                                        point0_xyz=None):
    """
    This function calculates a set of x, y, and z values based on a reference
    location, and dip direction. Here, x and y are horizontal distances and z
    is elevation relative to the point0_xy (negative means downwards).
    The point0_xy is assumed to be at 0 m above the surface.

    The dx and dy values used to scale each distance are calculated based on
    the dip direction (in degrees) and the input values for the corresponding
    length scale (L0, L1, L2, L3, L4). Only three points are needed to create
    a three-dimensional plane in the depth_change() function (p0, p1, and p2)
    but all five points are needed in the stratigraphy_plot() function in
    stratigraphy_plot.py.

    Point (p0) is at the location in point0_xy. The 2nd is at a distance of dx
    and dy from p0 in the direction of dip. The 3rd and 4th points (p2 and p3)
    are at distances dx and dy away from p0 in the directions of strike (one
    in each direction). The 5th point (p4) is at a distances of dx and dy from
    p0 in the direction of dip, and is used for the location of the number in
    the dip symbol. The dx and dy values are calculated with the length scales
    (L1, L2, L3, and L4) in a way that depends on dip direction.

    Note that L4 is doubled for certain values of dipDirectionDegrees because
    in those cases the dip number would often look like it is right on top of
    the strike and dip symbol (for the default perspectives set by view_elev
    and view_azimuth).

    :param dip: Dip angle in degrees
    :type dip: int or float

    :param dipDirectionDegrees: Direction of dip in degrees clockwise from
        north, so that 90 is east, 180 is south, and 270 is west.
    :type dipDirectionDegrees: int or float

    :param L1: Length scale (m) used to calculate the position of point 1,
        which is described by x_p1, y_p1, and z_p1. Point 1 is in the direction
        of dip.
    :type L1: int or float

    :param L2: Length scale (m) used to calculate the position of point 2,
        which is described by x_p2, y_p2, and z_p2. Point 2 is in the direction
        of strike, 90 degrees counter clockwise from the dip direction.
    :type L2: int or float

    :param L3: Length scale (m) used to calculate the position of point 3,
        which is described by x_p3, y_p3, and z_p3. Point 3 is in the direction
        of strike, 90 degrees clockwise from the dip direction.
    :type L3: int or float

    :param L4: Length scale (m) used to calculate the position of point 4,
        which is described by x_p4, y_p4, and z_p4. Point 4 is in the direction
        of dip, but is meant to be farther down dip than point 1. This point is
        used for the placement of dip numbers in the stratigraphy_plot()
        function.
    :type L3: int or float

    :param point0_xyz: Array of length 2 containing the x, y, and z values for
        point 0 (x_p0, y_p0, and z_p0). In the stratigraphy _plot() function,
        point 0 is taken as the location for the strike and dip symbol
        (SandD_Locaiton).
    :type point0_xy: array

    :returns: x_points, y_points, and z_points
    """
    if point0_xyz is None:
        x_p0, y_p0, z_p0 = 0, 0, 0
    else:
        x_p0 = point0_xyz[0]
        y_p0 = point0_xyz[1]
        z_p0 = point0_xyz[2]

    if dipDirectionDegrees not in DIP_DIRECTION_DEGREE_OPTIONS:
        if 0 < dipDirectionDegrees < 90:
            dx = L1 * np.sin(np.radians(dipDirectionDegrees))
            dy = L1 * np.cos(np.radians(dipDirectionDegrees))
            x_p1 = x_p0 + dx
            y_p1 = y_p0 + dy
            z_p1 = (-L1 * np.tan(np.radians(dip))) + z_p0

            dx = L2 * np.cos(np.radians(dipDirectionDegrees))
            dy = L2 * np.sin(np.radians(dipDirectionDegrees))
            x_p2 = x_p0 - dx
            y_p2 = y_p0 + dy
            z_p2 = z_p0

            dx = L3 * np.cos(np.radians(dipDirectionDegrees))
            dy = L3 * np.sin(np.radians(dipDirectionDegrees))
            x_p3 = x_p0 + dx
            y_p3 = y_p0 - dy
            z_p3 = z_p0

            dx = L4 * np.sin(np.radians(dipDirectionDegrees))
            dy = L4 * np.cos(np.radians(dipDirectionDegrees))
            x_p4 = x_p0 + dx
            y_p4 = y_p0 + dy
            z_p4 = (-L4 * np.tan(np.radians(dip))) + z_p0

        elif 90 < dipDirectionDegrees < 180:
            dx = L1 * np.cos(np.radians(dipDirectionDegrees - 90))
            dy = L1 * np.sin(np.radians(dipDirectionDegrees - 90))
            x_p1 = x_p0 + dx
            y_p1 = y_p0 - dy
            z_p1 = (-L1 * np.tan(np.radians(dip))) + z_p0

            dx = L2 * np.cos(np.radians(180 - dipDirectionDegrees))
            dy = L2 * np.sin(np.radians(180 - dipDirectionDegrees))
            x_p2 = x_p0 + dx
            y_p2 = y_p0 + dy
            z_p2 = z_p0

            dx = L3 * np.cos(np.radians(180 - dipDirectionDegrees))
            dy = L3 * np.sin(np.radians(180 - dipDirectionDegrees))
            x_p3 = x_p0 - dx
            y_p3 = y_p0 - dy
            z_p3 = z_p0

            dx = L4 * 2 * np.cos(np.radians(dipDirectionDegrees - 90))
            dy = L4 * 2 * np.sin(np.radians(dipDirectionDegrees - 90))
            x_p4 = x_p0 + dx
            y_p4 = y_p0 - dy
            z_p4 = (-L4 * 2 * np.tan(np.radians(dip))) + z_p0

        elif 180 < dipDirectionDegrees < 270:
            dx = L1 * np.sin(np.radians(dipDirectionDegrees - 180))
            dy = L1 * np.cos(np.radians(dipDirectionDegrees - 180))
            x_p1 = x_p0 - dx
            y_p1 = y_p0 - dy
            z_p1 = (-L1 * np.tan(np.radians(dip))) + z_p0

            dx = L2 * np.cos(np.radians(dipDirectionDegrees - 180))
            dy = L2 * np.sin(np.radians(dipDirectionDegrees - 180))
            x_p2 = x_p0 + dx
            y_p2 = y_p0 - dy
            z_p2 = z_p0

            dx = L3 * np.cos(np.radians(dipDirectionDegrees - 180))
            dy = L3 * np.sin(np.radians(dipDirectionDegrees - 180))
            x_p3 = x_p0 - dx
            y_p3 = y_p0 + dy
            z_p3 = z_p0

            dx = L4 * np.sin(np.radians(dipDirectionDegrees - 180))
            dy = L4 * np.cos(np.radians(dipDirectionDegrees - 180))
            x_p4 = x_p0 - dx
            y_p4 = y_p0 - dy
            z_p4 = (-L4 * np.tan(np.radians(dip))) + z_p0

        elif 270 < dipDirectionDegrees < 360:
            dx = L1 * np.cos(np.radians(dipDirectionDegrees - 270))
            dy = L1 * np.sin(np.radians(dipDirectionDegrees - 270))
            x_p1 = x_p0 - dx
            y_p1 = y_p0 + dy
            z_p1 = (-L1 * np.tan(np.radians(dip))) + z_p0

            dx = L2 * np.cos(np.radians(360 - dipDirectionDegrees))
            dy = L2 * np.sin(np.radians(360 - dipDirectionDegrees))
            x_p2 = x_p0 - dx
            y_p2 = y_p0 - dy
            z_p2 = z_p0

            dx = L3 * np.cos(np.radians(360 - dipDirectionDegrees))
            dy = L3 * np.sin(np.radians(360 - dipDirectionDegrees))
            x_p3 = x_p0 + dx
            y_p3 = y_p0 + dy
            z_p3 = z_p0

            dx = L4 * 2 * np.cos(np.radians(dipDirectionDegrees - 270))
            dy = L4 * 2 * np.sin(np.radians(dipDirectionDegrees - 270))
            x_p4 = x_p0 - dx
            y_p4 = y_p0 + dy
            z_p4 = (-L4 * 2 * np.tan(np.radians(dip))) + z_p0

    elif dipDirectionDegrees in [0, 360]:
        dx = 0
        dy = L1
        x_p1 = x_p0 + dx
        y_p1 = y_p0 + dy
        z_p1 = (-np.tan(np.radians(dip)) * L1) + z_p0

        dx = L2
        dy = 0
        x_p2 = x_p0 + dx
        y_p2 = y_p0 + dy
        z_p2 = z_p0

        dx = L3
        dy = 0
        x_p3 = x_p0 - dx
        y_p3 = y_p0 - dy
        z_p3 = z_p0

        dx = 0
        dy = L4
        x_p4 = x_p0 + dx
        y_p4 = y_p0 + dy
        z_p4 = (-np.tan(np.radians(dip)) * L4) + z_p0

    elif dipDirectionDegrees == 180:
        dx = 0
        dy = L1
        x_p1 = x_p0 + dx
        y_p1 = y_p0 - dy
        z_p1 = (-np.tan(np.radians(dip)) * L1) + z_p0

        dx = L2
        dy = 0
        x_p2 = x_p0 + dx
        y_p2 = y_p0 + dy
        z_p2 = z_p0

        dx = L3
        dy = 0
        x_p3 = x_p0 - dx
        y_p3 = y_p0 - dy
        z_p3 = z_p0

        dx = 0
        dy = L4 * 2
        x_p4 = x_p0 + dx
        y_p4 = y_p0 - dy
        z_p4 = (-np.tan(np.radians(dip)) * L4 * 2) + z_p0

    elif dipDirectionDegrees == 90:
        dx = L1
        dy = 0
        x_p1 = x_p0 + dx
        y_p1 = y_p0 + dy
        z_p1 = (-np.tan(np.radians(dip)) * L1) + z_p0

        dx = 0
        dy = L2
        x_p2 = x_p0 + dx
        y_p2 = y_p0 + dy
        z_p2 = z_p0

        dx = 0
        dy = L3
        x_p3 = x_p0 - dx
        y_p3 = y_p0 - dy
        z_p3 = z_p0

        dx = L4
        dy = 0
        x_p4 = x_p0 + dx
        y_p4 = y_p0 + dy
        z_p4 = (-np.tan(np.radians(dip)) * L4) + z_p0

    elif dipDirectionDegrees == 270:
        dx = L1
        dy = 0
        x_p1 = x_p0 - dx
        y_p1 = y_p0 + dy
        z_p1 = (-np.tan(np.radians(dip)) * L1) + z_p0

        dx = 0
        dy = L2
        x_p2 = x_p0 + dx
        y_p2 = y_p0 + dy
        z_p2 = z_p0

        dx = 0
        dy = L3
        x_p3 = x_p0 - dx
        y_p3 = y_p0 - dy
        z_p3 = z_p0

        dx = L4
        dy = 0
        x_p4 = x_p0 - dx
        y_p4 = y_p0 + dy
        z_p4 = (-np.tan(np.radians(dip)) * L4) + z_p0

    x_points = [x_p0, x_p1, x_p2, x_p3, x_p4]
    y_points = [y_p0, y_p1, y_p2, y_p3, y_p4]
    z_points = [z_p0, z_p1, z_p2, z_p3, z_p4]

    return x_points, y_points, z_points


def get_lut_stratigraphy_dict(file_name, file_directory, comp_locX, comp_locY,
                              max_point_distance=100):
    """
    Function that reads stratigraphy information from a csv file. The file is
    expected to have x and y columns representing distances east and north for
    each point to be used by components. The unit thicknesses in each row are
    expected to be representative of the corresponding x and y values. The
    minimum and maximum values can be provided (e.g., aquifer2Thickness_min, or
    shale3Thickness_max, or reservoirThickness_min). Parameters can also be
    entered as single values that are representative of the entire domain. The
    parameters for a given comp_locX and comp_locY values are taken as those
    from the closest point from the file. To prevent an unrepresentative value
    from being used, however, points must be within max_point_distance of the
    comp_locX and comp_locY location.
    """
    file_path_name = os.path.join(iam.IAM_DIR, file_directory, file_name)

    data = pd.read_csv(file_path_name)

    try:
        numShaleLayers = data.numberOfShaleLayers
        numShaleLayers = numShaleLayers.values
        numShaleLayers = numShaleLayers[np.isnan(numShaleLayers) != True]
        numShaleLayers = int(numShaleLayers)
    except (AttributeError, TypeError):
        numShaleLayers = None
        warning_msg = ''.join([
            'The LookupTableStratigraphy file ', file_name, ' did not have a '
            'valid entry for numberOfShaleLayers. Enter the parameter in its ',
            'own column with only one row. The default value for ',
            'numberOfShaleLayers will be used.'])
        logging.warning(warning_msg)

    x_values = data.x
    x_values = x_values.values

    y_values = data.y
    y_values = y_values.values

    # Find the closest point
    for LUTS_index, xval in enumerate(x_values):
        dx, dy = (comp_locX - xval), (comp_locY - y_values[LUTS_index])
        point_distance = (np.abs(dx)**2 + np.abs(dy)**2)**0.5
        if LUTS_index == 0:
            selected_index = LUTS_index
            last_point_distance = point_distance
        elif point_distance < last_point_distance:
            selected_index = LUTS_index
            last_point_distance = point_distance

    point_distance = (((np.abs((comp_locX - x_values[selected_index])))
                      ** 2) + ((np.abs((comp_locY - y_values[selected_index])))
                               ** 2)) ** 0.5

    if point_distance > max_point_distance:
        err_msg = ''.join([
            'The closest location in the LookupTableStratigraphy file {} ',
            '(x = {} m and y = {} m) was still {} m away from the component ',
            'at x = {} m and y = {} m. Because the MaxPointDistance is set ',
            'to {} m, no stratigraphy information could be used ',
            'at this location and the simulation could not proceed.']).format(
                file_name, x_values[selected_index], y_values[selected_index],
                point_distance, comp_locX, comp_locY, max_point_distance)
        raise ValueError(err_msg)

    par_list_keys = ['shaleThickness_List', 'aquiferThickness_List',
                     'shaleThickness_min_List', 'shaleThickness_max_List',
                     'aquiferThickness_min_List', 'aquiferThickness_max_List']

    # Initialize output dictionary
    LUTStrat_dict = {}
    for key in par_list_keys:
        LUTStrat_dict[key] = []

    LUTStrat_dict['numberOfShaleLayers'] = numShaleLayers

    # Get the unit thicknesses
    # Shale layers
    for shaleRef in range(1, numShaleLayers+1):
        par_nm_start = 'shale{}'.format(shaleRef)
        for par_nm_end in ['Thickness', 'Thickness_min', 'Thickness_max']:
            parameter_name = par_nm_start + par_nm_end
            unitThickness_temp = get_parameter_from_lut_stratigraphy(
                data, parameter_name, selected_index)
            LUTStrat_dict['shale{}_List'.format(par_nm_end)].append(unitThickness_temp)

    # Aquifer layers
    for shaleRef in range(1, numShaleLayers):
        par_nm_start = 'aquifer{}'.format(shaleRef)
        for par_nm_end in ['Thickness', 'Thickness_min', 'Thickness_max']:
            parameter_name = par_nm_start + par_nm_end
            unitThickness_temp = get_parameter_from_lut_stratigraphy(
                data, parameter_name, selected_index)
            LUTStrat_dict['aquifer{}_List'.format(par_nm_end)].append(unitThickness_temp)

    # Get other parameters
    par_nm_start = 'reservoir'
    for par_nm_end in ['Thickness', 'Thickness_min', 'Thickness_max']:
        parameter_name = par_nm_start + par_nm_end
        unitThickness_temp = get_parameter_from_lut_stratigraphy(
            data, parameter_name, selected_index)
        LUTStrat_dict['res{}'.format(par_nm_end)] = unitThickness_temp

    for parameter_name in ['depth', 'depth_min', 'depth_max',
                           'datumPressure', 'datumPressure_min', 'datumPressure_max']:
        val_temp = get_parameter_from_lut_stratigraphy(
            data, parameter_name, selected_index)
        LUTStrat_dict[parameter_name] = val_temp

    return LUTStrat_dict


def get_parameter_from_lut_stratigraphy(data, parameter_name, selected_index):
    """
    Function that takes .csv data read with pandas (data), finds the column
    corresponding with the parameter_name provided (e.g., shale1Thickness), and
    checks if the parameter is given as a list or single value. If it is a list,
    then only the value in the position of selected_index is returned. Otherwise,
    the single value is returned.
    """
    try:
        LUTS_value = data.eval(parameter_name)
        LUTS_value = LUTS_value.values
        LUTS_value = LUTS_value[np.isnan(LUTS_value) != True]
    except NameError:
        LUTS_value = None

    if not LUTS_value is None:
        try:
            LUTS_value = LUTS_value[selected_index]
        except (IndexError, TypeError):
            LUTS_value = float(LUTS_value)

    return LUTS_value


def get_unit_depth_from_component(numShaleLayers, stratigraphyComponent,
                                  unitNumber=1, unitType='shale',
                                  top_or_bottom='top'):
    """
    This function  provides the top or bottom depth of a unit.

    Note that this function provides the unit depths that are derived from the
    simulation input, such as unit thicknesses, strike, and dip. In contrast,
    the function get_strata_info_from_component() deals only with the input provided
    directly to stratigraphy components (e.g., unit thicknesses and datumPressure).

    Unit depths can also be obtained with the composite parameters made during
    the connect_with_system() method of the Stratigraphy component. Those
    parameters work in forward simulations but are 0 in LHS and parstudy
    simulations, and this function solves that problem.

    :param numShaleLayers: Number of shale layers
    :type numShaleLayers: int

    :param stratigraphyComponent: Stratigraphy component used to obtain information
    :type stratigraphyComponent: instance of Stratigraphy component class

    :param unitNumber: The number for the unit targeted (e.g., 2 for aquifer 2)
    :type unitNumber: int

    :param unitType: The type of unit. This entry can only be 'shale', 'aquifer',
        or 'reservoir.'
    :type unitType: str

    :param top_or_bottom: Option to obtain the 'top' or 'bottom' depth of a unit
    :type top_or_bottom: str

    :returns: unitDepth
    """
    res_bottom_check = False

    # If the user wants the top depth of the reservoir, obtain the bottom of
    # shale 1 (deepest shale) instead.
    if unitType == 'reservoir':
        unitType = 'shale'
        unitNumber = 1

        # If the user wants the bottom depth of the reservoir, add the reservoir
        # thickness at the end.
        if top_or_bottom == 'bottom':
            res_bottom_check = True
        elif top_or_bottom == 'top':
            # The user wanted the top of the reservoir, which is the bottom of
            # shale 1
            top_or_bottom = 'bottom'

    unitDepth = 0

    # First, handle the shallowest shale
    shale_par_nm = 'shale{}Thickness'.format(numShaleLayers)

    shaleThickness = iam.cfi.commons.get_parameter_val(
        stratigraphyComponent, shale_par_nm)

    if unitNumber < numShaleLayers:
        unitDepth += shaleThickness  # start with thickness of the top shale
    else:  # unitNumber == numShaleLayers
        if unitType == 'shale' and top_or_bottom == 'bottom':
            unitDepth += shaleThickness

    if unitNumber != numShaleLayers:
        for shaleRef in range(numShaleLayers - 1, unitNumber - 1, -1):
            aq_par_nm = 'aquifer{}Thickness'.format(shaleRef)
            shale_par_nm = 'shale{}Thickness'.format(shaleRef)

            aquiferThickness = iam.cfi.commons.get_parameter_val(
                stratigraphyComponent, aq_par_nm)

            shaleThickness = iam.cfi.commons.get_parameter_val(
                stratigraphyComponent, shale_par_nm)

            if shaleRef > unitNumber:
                unitDepth += aquiferThickness + shaleThickness
            elif shaleRef == unitNumber and unitType == 'shale':
                unitDepth += aquiferThickness
                if top_or_bottom == 'bottom':
                    unitDepth += shaleThickness
            elif shaleRef == unitNumber and unitType == 'aquifer':
                if top_or_bottom == 'bottom':
                    unitDepth += aquiferThickness

    if res_bottom_check:
        par_name = 'reservoirThickness'

        reservoirThickness = iam.cfi.commons.get_parameter_val(
            stratigraphyComponent, par_name)

        unitDepth += reservoirThickness

    return unitDepth


def check_color_alpha_label_yaml_input(yaml_input, strata_plot_data, name):
    """
    Function that checks if color, alpha, or labels provided for stratigraphic
    units in a .yaml file are acceptable. If the input is recognized and deemed
    acceptable, the input is added to yaml_input. If not, the input is excluded
    from yaml_input and a warning message is logged and printed.

    :param yaml_input: Dictionary of keys provided for the plotting function
        used (e.g., stratigraphy_plot() or stratigraphic_column).
    :type yaml_input: dict

    :param strata_plot_data: Input values provided for the plot in the Plots
        section of a .yaml file.
    :type strata_plot_data: dict

    :param name: Name of the plot
    :type name: str
    """
    alphaNames = ['ReservoirAlpha', 'ShaleAlpha', 'AquiferAlpha', 'WellAlpha']

    alphaNames += ['Shale' + str(num) + 'Alpha' for num in range(1, 31)]
    alphaNames += ['Aquifer' + str(num) + 'Alpha' for num in range(1, 30)]

    colorNames = ['ReservoirColor', 'ShaleColor', 'AquiferColor', 'WellColor']

    colorNames += ['Shale' + str(num) + 'Color' for num in range(1, 31)]
    colorNames += ['Aquifer' + str(num) + 'Color' for num in range(1, 30)]

    labelNames = ['ReservoirLabel', 'WellLabel']

    labelNames += ['Shale' + str(num) + 'Label' for num in range(1, 31)]
    labelNames += ['Aquifer' + str(num) + 'Label' for num in range(1, 30)]

    for input_type in alphaNames:
        if input_type in strata_plot_data:
            try:
                if 0 < strata_plot_data[input_type] <= 1:
                    # The alpha input corresponds with the alpha used for filled areas
                    yaml_input[input_type + 'Fill'] = strata_plot_data[input_type]

                    # The alpha value used for lines is higher
                    yaml_input[input_type] = strata_plot_data[input_type] + (
                        (1 - strata_plot_data[input_type]) * ALPHA_SCALE_FACTOR)
                else:
                    color_alpha_label_error_msg(
                        input_type, strata_plot_data[input_type], 'alpha', name)
            except:
                color_alpha_label_error_msg(
                    input_type, strata_plot_data[input_type], 'alpha', name)

    for input_type in colorNames:
        if input_type in strata_plot_data:
            if isinstance(strata_plot_data[input_type], str):
                if is_color_like(strata_plot_data[input_type]):
                    yaml_input[input_type] = strata_plot_data[input_type]
                else:
                    color_alpha_label_error_msg(
                        input_type, strata_plot_data[input_type], 'str', name)

            elif isinstance(strata_plot_data[input_type], list):
                if not is_color_like(strata_plot_data[input_type]):
                    color_alpha_label_error_msg(
                        input_type, strata_plot_data[input_type], 'list', name)
                else:
                    yaml_input[input_type] = strata_plot_data[input_type]
            else:
                color_alpha_label_error_msg(
                    input_type, strata_plot_data[input_type], 'neither', name)

    for input_type in labelNames:
        if input_type in strata_plot_data:
            if isinstance(strata_plot_data[input_type], str):
                yaml_input[input_type] = strata_plot_data[input_type]
            else:
                color_alpha_label_error_msg(
                    input_type, strata_plot_data[input_type], 'label', name)

    return yaml_input


def color_alpha_label_error_msg(input_type, input_value, error_type, name):
    """
    Function that logs and prints an error message when invalid input is
    provided for a unit color, unit alpha, or unit label.
    """
    err_msg_pt1 = ''.join(['The ', input_type, ' input provided for the figure ',
                           name, '(', input_type, ':', ' ', input_value, ')'])
    if error_type == 'alpha':
        err_msg_pt2 = ''.join([' was not a value between 0 and 1.'])
    elif error_type == 'str':
        err_msg_pt2 = ''.join([' was a string input, but not one that matplotlib ',
                               'recognized as a color (e.g., "r", "teal", ',
                               'or "#030764").'])
    elif error_type == 'list':
        err_msg_pt2 = ''.join([' was a list, but not a list of length three ',
                               'that matplotlib recognized as a color (e.g., ',
                               '"[1, 0, 0]" or "[0.67, 0.33, 0]"). The three ',
                               'values in the list must be between 0 and 1 ',
                               '(for fractions of [Red, Green, Blue]).'])
    elif error_type == 'neither':
        err_msg_pt2 = ''.join([' was neither a string nor a list that matplotlib ',
                               'recognized as a color (e.g., "r", "teal", "#030764", ',
                               '"[1, 0, 0]," or "[0.67, 0.33, 0]"). If given as ',
                               'a list, the three values in the list must be ',
                               'between 0 and 1 (for fractions of [Red, Green, Blue]).'])
    elif error_type == 'label':
        err_msg_pt2 = ''.join([' was not a string. Only string inputs (e.g., ',
                               '"AZMI", "BGWP", or "Freshwater Aquifer") can be ',
                               'provided for unit labels. '])

    err_msg_pt3 = ' The figure {} will use the default approach for setting {}.'.format(
        name, input_type)

    err_msg = err_msg_pt1 + err_msg_pt2 + err_msg_pt3

    logging.debug(err_msg)


def get_plotting_setup_for_units(strat_plot_yaml_input, numShaleLayers,
                                 reservoirThickness=None, shaleThicknessList=None,
                                 aquiferThicknessList=None, include_thickness=False):
    """
    Function that provides the colors, alphas, and labels for each unit. The
    input provided in strat_plot_yaml_input is used if present. Otherwise,
    default values are used.

    :param strata_plot_data: Dictionary containing input values for a plot.
        This dictionary should already have been updated with the function
        check_color_alpha_label_yaml_input().
    :type strata_plot_data: dict

    :param numShaleLayers: Number of shale layers.
    :type numShaleLayers: int

    :param reservoirThickness: Thickness of the reservoir (m)
    :type reservoirThickness: int or float

    :param shaleThicknessList: List of shale thicknesses (m)
    :type shaleThicknessList: list

    :param aquiferThicknessList: List of aquifer thicknesses (m)
    :type aquiferThicknessList: list
    """
    reservoirColor = strat_plot_yaml_input.get(
        'ReservoirColor', UNIT_COLOR_DICT['ReservoirColor'])
    reservoirAlpha = strat_plot_yaml_input.get(
        'ReservoirAlpha', UNIT_ALPHA_DICT['ReservoirAlpha'])
    reservoirAlphaFill = strat_plot_yaml_input.get(
        'ReservoirAlphaFill', UNIT_ALPHA_FILL_DICT['ReservoirAlphaFill'])
    if include_thickness and not reservoirThickness is None:
        reservoirLabel = strat_plot_yaml_input.get(
            'ReservoirLabel', UNIT_LABEL_THICKNESS_DICT['ReservoirLabel'].format(
                reservoirThickness))
    else:
        reservoirLabel = strat_plot_yaml_input.get(
            'ReservoirLabel', UNIT_LABEL_DICT['ReservoirLabel'])

    wellColor = strat_plot_yaml_input.get('WellColor', WELL_COLOR)
    wellAlpha = strat_plot_yaml_input.get('WellAlpha', WELL_ALPHA)
    wellAlphaFill = strat_plot_yaml_input.get('WellAlphaFill', WELL_ALPHA_FILL)
    wellLabel = strat_plot_yaml_input.get('WellLabel', None)

    shaleColor = [None] * numShaleLayers
    shaleAlpha = [None] * numShaleLayers
    shaleAlphaFill = [None] * numShaleLayers
    shaleLabel = [None] * numShaleLayers

    aquiferColor = [None] * (numShaleLayers - 1)
    aquiferAlpha = [None] * (numShaleLayers - 1)
    aquiferAlphaFill = [None] * (numShaleLayers - 1)
    aquiferLabel = [None] * (numShaleLayers - 1)

    for shaleRef in range(0, numShaleLayers):
        shaleColor[shaleRef] = strat_plot_yaml_input.get(
            'ShaleColor', UNIT_COLOR_DICT['ShaleColor'])

        # More specific color input overrides the more general color input
        if 'Shale{:.0f}Color'.format(shaleRef + 1) in strat_plot_yaml_input:
            shaleColor[shaleRef] = strat_plot_yaml_input[
                'Shale{:.0f}Color'.format(shaleRef + 1)]

        shaleAlpha[shaleRef] = strat_plot_yaml_input.get(
            'ShaleAlpha', UNIT_ALPHA_DICT['ShaleAlpha'])

        shaleAlphaFill[shaleRef] = strat_plot_yaml_input.get(
            'ShaleAlphaFill', UNIT_ALPHA_FILL_DICT['ShaleAlphaFill'])

        # More specific alpha input overrides the more general color input
        if 'Shale{:.0f}Alpha'.format(shaleRef + 1) in strat_plot_yaml_input:
            shaleAlpha[shaleRef] = strat_plot_yaml_input[
                'Shale{:.0f}Alpha'.format(shaleRef + 1)]

        if 'Shale{:.0f}AlphaFill'.format(shaleRef + 1) in strat_plot_yaml_input:
            shaleAlphaFill[shaleRef] = strat_plot_yaml_input[
                'Shale{:.0f}AlphaFill'.format(shaleRef + 1)]

        if include_thickness and not shaleThicknessList is None:
            shaleLabel[shaleRef] = strat_plot_yaml_input.get(
                'Shale{:.0f}Label'.format(shaleRef + 1),
                UNIT_LABEL_THICKNESS_DICT['ShaleLabel'].format(
                    shaleRef + 1, shaleThicknessList[shaleRef]))
        else:
            shaleLabel[shaleRef] = strat_plot_yaml_input.get(
                'Shale{:.0f}Label'.format(shaleRef + 1),
                UNIT_LABEL_DICT['ShaleLabel'].format(shaleRef + 1))

        # Get input for aquifers
        if shaleRef != (numShaleLayers - 1):
            aquiferColor[shaleRef] = strat_plot_yaml_input.get(
                'Aquifer{:.0f}Color'.format(shaleRef + 1), UNIT_COLOR_DICT['AquiferColor'])

            # More specific color input overrides the more general color input
            if 'Aquifer{:.0f}Color'.format(shaleRef + 1) in strat_plot_yaml_input:
                aquiferColor[shaleRef] = strat_plot_yaml_input[
                    'Aquifer{:.0f}Color'.format(shaleRef + 1)]

            aquiferAlpha[shaleRef] = strat_plot_yaml_input.get(
                'AquiferAlpha', UNIT_ALPHA_DICT['AquiferAlpha'])

            aquiferAlphaFill[shaleRef] = strat_plot_yaml_input.get(
                'AquiferAlphaFill', UNIT_ALPHA_FILL_DICT['AquiferAlphaFill'])

            # More specific alpha input overrides the more general color input
            if 'Aquifer{:.0f}Alpha'.format(shaleRef + 1) in strat_plot_yaml_input:
                aquiferAlpha[shaleRef] = strat_plot_yaml_input[
                    'Aquifer{:.0f}Alpha'.format(shaleRef + 1)]

            if 'Aquifer{:.0f}AlphaFill'.format(shaleRef + 1) in strat_plot_yaml_input:
                aquiferAlphaFill[shaleRef] = strat_plot_yaml_input[
                    'Aquifer{:.0f}AlphaFill'.format(shaleRef + 1)]

            if include_thickness and not aquiferThicknessList is None:
                aquiferLabel[shaleRef] = strat_plot_yaml_input.get(
                    'Aquifer{:.0f}Label'.format(shaleRef + 1),
                    UNIT_LABEL_THICKNESS_DICT['AquiferLabel'].format(
                        shaleRef + 1, aquiferThicknessList[shaleRef]))
            else:
                aquiferLabel[shaleRef] = strat_plot_yaml_input.get(
                    'Aquifer{:.0f}Label'.format(shaleRef + 1),
                    UNIT_LABEL_DICT['AquiferLabel'].format(shaleRef + 1))

    return reservoirColor, reservoirAlpha, reservoirAlphaFill, reservoirLabel, \
        shaleColor, shaleAlpha, shaleAlphaFill, shaleLabel, \
            aquiferColor, aquiferAlpha, aquiferAlphaFill, aquiferLabel, \
                wellColor, wellAlpha, wellAlphaFill, wellLabel
