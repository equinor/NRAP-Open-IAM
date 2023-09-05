"""
Last modified: September 6th, 2022

Authors: Veronika Vasylkivska, Nate Mitchell
"""
import os
import sys
import logging
import math
import random
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

try:
    from openiam import IAM_DIR
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))


def get_known_locations(cmpnt_name, cmpnt_data):
    """
    Check in which form the locations for the component are provided.

    Returns dictionary containing known coordinates for placeable components:
    e.g., reservoir, wellbore, seal horizon or fault flow components.

    :param cmpnt_name: name of component model
    :type name: str

    :param cmpnt_data: component data providing its setup parameters
    :type name: dict

    Returns
    -------
    locations: dictionary with coordinates
    """
    loc_data = cmpnt_data['Locations']  # this assumes that key 'Locations' is present in cmpnt_data

    if 'file' in loc_data:
        return file_locations(cmpnt_name, loc_data)

    if 'linspace' in loc_data:
        return linspace_locations(cmpnt_name, loc_data)

    if 'grid' in loc_data:
        return grid_locations(cmpnt_name, loc_data)

    if 'range' in loc_data:
        return range_locations(cmpnt_name, loc_data)

    if 'coordx' in loc_data or 'coordy' in loc_data:
        return coord_locations(cmpnt_name, loc_data)
    else:
        debug_msg = 'No location information is provided for component {}.'.format(
            cmpnt_name)
        logging.debug(debug_msg)
        locations = {}  # empty dictionary
        return locations


def get_random_locations(cmpnt_name, cmpnt_data, num_rand_locs):
    """
    Generate required number of random locations based on the parameters provided by user.

    Parameters
    ----------
    Returns dictionary containing randomly generated coordinates for wellbore components.

    :param cmpnt_name: name of component model
    :type name: str

    :param cmpnt_data: component data providing its setup parameters
    :type name: dict

    :param num_rand_locs: number of locations to be generated
    :type num_rand_locs: int

    Returns
    -------
    locations: dictionary with coordinates

    """
    # Assuming that key 'RandomLocDomain' is present in cmpnt_data
    rand_loc_data = cmpnt_data['RandomLocDomain']
    locations = {}
    if 'seed' in rand_loc_data:
        loc_seed = rand_loc_data['seed']
    elif 'Seed' in rand_loc_data:
        loc_seed = rand_loc_data['Seed']
    else:
        loc_seed = random.randint(1, 10000)

    # Get minimums and maximums of the uniform distribution
    mins = {'x': rand_loc_data['xmin'], 'y': rand_loc_data['ymin']}
    maxs = {'x': rand_loc_data['xmax'], 'y': rand_loc_data['ymax']}

    # Generate additional random locations
    rng = np.random.default_rng(loc_seed)

    rand_xloc = rng.uniform(mins['x'], maxs['x'], num_rand_locs).tolist()
    rand_yloc = rng.uniform(mins['y'], maxs['y'], num_rand_locs).tolist()
    locations['coordx'] = rand_xloc
    locations['coordy'] = rand_yloc

    # Check if boundaries for z-coordinates are also provided
    try:
        mins['z'] = rand_loc_data['zmin']
        maxs['z'] = rand_loc_data['zmax']
    except KeyError: # if zmin and/or zmax not in comp_data['RandomLocDomain']
        for cx, cy in zip(rand_xloc, rand_yloc):
            debug_msg = 'Creating random location for component {} at ({}, {}).'.format(
                cmpnt_name, cx, cy)
            logging.debug(debug_msg)
    else:
        rand_zloc = rng.uniform(mins['z'], maxs['z'], num_rand_locs).tolist()
        locations['coordz'] = rand_zloc
        for cx, cy, cz in zip(rand_xloc, rand_yloc, rand_zloc):
            debug_msg = 'Creating random location for component {} at ({}, {}, {}).'.format(
                cmpnt_name, cx, cy, cz)
            logging.debug(debug_msg)

    locations['number'] = num_rand_locs

    return locations


def file_locations(cmpnt_name, loc_data):
    """
    Get locations provided through file.

    Parameters
    ----------
    :param cmpnt_name: name of component model
    :type name: str

    :param loc_data: locations data
    :type name: dict

    Raises
    ------
    FileNotFoundError
        The error is raised if the file path is not valid.

    Returns
    -------
    dictionary locations with coordinates

    """
    # Define output
    locations = {}
    inp_data_file_path = os.path.join(IAM_DIR, loc_data['file'])
    if os.path.isfile(inp_data_file_path):
        data = np.genfromtxt(
            inp_data_file_path, delimiter=",", dtype='f8', comments='#')

        # If the first row doesn't start with '#' and the first row has labels,
        # that row will be read as NaN values.
        if np.isnan(data[0, 0]):
            data = data[1:None, :]

        locations['coordx'] = data[:, 0].tolist()
        locations['coordy'] = data[:, 1].tolist()
        locations['number'] = data.shape[0]

        # By default, z-values are not provided and should not be considered
        read_z_values = loc_data.get('read_z_values', False)

        if read_z_values: # this will serve as flag of presence of z-values
            try:
                locations['coordz'] = data[:, 2].tolist()
            except IndexError:
                err_msg = ''.join(['Option "read_z_values" is set to True ',
                                   'but z-coordinates are not available.'])
                logging.error(err_msg)
                raise IndexError(err_msg) from None

        # If you are reading points from a lookup reservoir table, the point
        # densities can become high near the injection wells. Handling a large
        # number of points is computationally taxing, and you might not need
        # a high number of points (e.g., for the AoR_plot() function). This
        # section lets you thin out the points. Specifically, this code
        # prevents the inclusion of points that are within min_x_spacing and
        # min_y_spacing from a previously included point (for the x and y
        # directions, respectively).
        if 'thin_point_density' in loc_data:
            if loc_data['thin_point_density']:
                # Get minimum x- and y-spacing if present, otherwise set to np.Inf
                min_x_spacing = loc_data.get('min_x_spacing', np.Inf)
                min_y_spacing = loc_data.get('min_y_spacing', np.Inf)

                # Determine unique x- and y-values
                unique_x_vals = np.unique(locations['coordx'])
                unique_y_vals = np.unique(locations['coordy'])
                chosen_x_vals = np.array([])
                chosen_y_vals = np.array([])

                last_point = unique_x_vals[0]
                for x_val in unique_x_vals[1:]:
                    if x_val - last_point >= min_x_spacing:
                        chosen_x_vals = np.append(chosen_x_vals, x_val)
                        last_point = x_val

                last_point = unique_y_vals[0]
                for y_val in unique_y_vals[1:]:
                    if y_val - last_point >= min_y_spacing:
                        chosen_y_vals = np.append(chosen_y_vals, y_val)
                        last_point = y_val

                # Go through the x- and y-values and add only the points with values
                # that are in the chosen_x_vals and chosen_y_vals
                x_loc_redo = []
                y_loc_redo = []
                ref_indices = []

                for loc_ref, (x_val, y_val) in enumerate(zip(locations['coordx'],
                                                             locations['coordy'])):
                    if x_val in chosen_x_vals and y_val in chosen_y_vals:
                        x_loc_redo.append(x_val)
                        y_loc_redo.append(y_val)
                        ref_indices.append(loc_ref)

                locations['coordx'] = x_loc_redo
                locations['coordy'] = y_loc_redo
                locations['number'] = len(x_loc_redo)
                if read_z_values:
                    z_loc_redo = [locations['coordz'][ind] for ind in ref_indices]
                    locations['coordz'] = z_loc_redo

    else:
        err_msg = ''.join(['File {} does not exist. Please check {} ',
                           'component setup.']).format(loc_data['file'],
                                                     cmpnt_name)
        logging.error(err_msg)
        raise FileNotFoundError(err_msg)

    debug_msg = 'Known locations for component {} are {}.'.format(
        cmpnt_name, locations)
    logging.debug(debug_msg)

    return locations


def linspace_locations(cmpnt_name, loc_data):
    """
    Generate locations using numpy linspace with parameters provided by user.

    Parameters
    ----------
    :param cmpnt_name: name of component model
    :type name: str

    :param loc_data: locations data
    :type name: dict

    Raises
    ------
    KeyError
        The error is raised if some of the required keys are missing.

    Returns
    -------
    dictionary locations with coordinates

    """
    # Define output
    locations = {}
    # Check whether all the needed parameters are present
    required_keys = ['xmin', 'xmax', 'ymin', 'ymax', 'size']
    missing_keys = check_missing_keys(loc_data['linspace'], required_keys)
    if missing_keys:  # if the list is not empty
        raise_missing_keys_error(cmpnt_name, 'linspace', missing_keys)
    else:
        # Explicitly stating that the endpoint is included
        locations['coordx'] = np.linspace(loc_data['linspace']['xmin'],
                                          loc_data['linspace']['xmax'],
                                          num=loc_data['linspace']['size'],
                                          endpoint=True).tolist()
        locations['coordy'] = np.linspace(loc_data['linspace']['ymin'],
                                          loc_data['linspace']['ymax'],
                                          num=loc_data['linspace']['size'],
                                          endpoint=True).tolist()
        locations['number'] = loc_data['linspace']['size']

        try: # get z coordinates if setup is provided
            locations['coordz'] = np.linspace(loc_data['linspace']['zmin'],
                                              loc_data['linspace']['zmax'],
                                              num=loc_data['linspace']['size'],
                                              endpoint=True).tolist()
        except KeyError:
            pass

    debug_msg = 'Known locations for component {} are {}.'.format(
        cmpnt_name, locations)
    logging.debug(debug_msg)

    return locations


def grid_locations(cmpnt_name, loc_data):
    """
    Generate locations using numpy.meshgrid method with parameters provided by user.

    The produced x and y arrays are flattened after that.
    y coordinate is assumed to change faster than x coordinates.
    If z coordinates are produced they are the ones changing the fastest, i.e,
    from one index to another.

    Parameters
    ----------
    :param cmpnt_name: name of component model
    :type name: str

    :param loc_data: locations data
    :type name: dict

    Raises
    ------
    KeyError
        The error is raised if some of the required keys are missing.

    Returns
    -------
    dictionary locations with coordinates

    """
    # Define output
    locations = {}
    # Start with checking first option: whether x and y are provided
    required_keys = ['x', 'y']
    missing_keys = check_missing_keys(loc_data['grid'], required_keys)
    if not missing_keys:  # if both keys 'x' and 'y' are present, missing_keys is empty list
        # Get x and y values
        uniq_x = loc_data['grid']['x']
        uniq_y = loc_data['grid']['y']
        try: # let's see if z array is also provided
            uniq_z = loc_data['grid']['z']
        except KeyError: # no z, only x and y arrays are provided
            yy, xx = np.meshgrid(uniq_y, uniq_x)
            # Flatten arrays and assign to the corresponding keys
            locations['coordx'] = xx.flatten().tolist()
            locations['coordy'] = yy.flatten().tolist()
            locations['number'] = len(locations['coordx'])
        else:  # x, y, and z arrays are provided
            yy, xx, zz = np.meshgrid(uniq_y, uniq_x, uniq_z)
            # Flatten arrays and assign to the corresponding keys
            locations['coordx'] = xx.flatten().tolist()
            locations['coordy'] = yy.flatten().tolist()
            locations['coordz'] = zz.flatten().tolist()
        locations['number'] = len(locations['coordx'])
    elif len(missing_keys) == 1:
        raise_missing_keys_error(cmpnt_name, 'grid', missing_keys)
    else:  # both keys 'x' and 'y' are missing
        # Check for the second option
        required_keys = ['xmin', 'xmax', 'xsize', 'ymin', 'ymax', 'ysize']
        missing_keys = check_missing_keys(loc_data['grid'], required_keys)
        if missing_keys:  # if list is not empty, then something is missing
            raise_missing_keys_error(cmpnt_name, 'grid', missing_keys)
        else:
            uniq_x = np.linspace(loc_data['grid']['xmin'],
                                 loc_data['grid']['xmax'],
                                 num=loc_data['grid']['xsize'],
                                 endpoint=True)
            uniq_y = np.linspace(loc_data['grid']['ymin'],
                                 loc_data['grid']['ymax'],
                                 num=loc_data['grid']['ysize'],
                                 endpoint=True)
            try:
                uniq_z = np.linspace(loc_data['grid']['zmin'],
                                     loc_data['grid']['zmax'],
                                     num=loc_data['grid']['zsize'],
                                     endpoint=True)
            except KeyError:  # no z data provided
                yy, xx = np.meshgrid(uniq_y, uniq_x)
                # Flatten arrays and assign to the corresponding keys
                locations['coordx'] = xx.flatten().tolist()
                locations['coordy'] = yy.flatten().tolist()
            else:  # x, y, and z data are provided
                yy, xx, zz = np.meshgrid(uniq_y, uniq_x, uniq_z)
                # Flatten arrays and assign to the corresponding keys
                locations['coordx'] = xx.flatten().tolist()
                locations['coordy'] = yy.flatten().tolist()
                locations['coordz'] = zz.flatten().tolist()
            locations['number'] = len(locations['coordx'])

    debug_msg = 'Known locations for component {} are {}.'.format(
        cmpnt_name, locations)
    logging.debug(debug_msg)

    return locations


def range_locations(cmpnt_name, loc_data):
    """
    Generate locations using numpy.arange method with parameters provided by user.

    Parameters
    ----------
    :param cmpnt_name: name of component model
    :type name: str

    :param loc_data: locations data
    :type name: dict

    Returns
    -------
    dictionary locations with coordinates

    """
    locations = {}
    # Check whether all the needed parameters are present
    required_keys = ['xmin', 'xmax', 'dx', 'ymin', 'ymax', 'dy']
    missing_keys = check_missing_keys(loc_data['range'], required_keys)
    if missing_keys:  # if some keys are missing
        raise_missing_keys_error(cmpnt_name, 'range', missing_keys)
    else:
        uniq_x = np.arange(loc_data['range']['xmin'],
                           loc_data['range']['xmax']+loc_data['range']['dx'],
                           loc_data['range']['dx']).tolist()
        uniq_y = np.arange(loc_data['range']['ymin'],
                           loc_data['range']['ymax']+loc_data['range']['dy'],
                           loc_data['range']['dy']).tolist()
        try:
            uniq_z = np.arange(loc_data['range']['zmin'],
                               loc_data['range']['zmax']+loc_data['range']['dz'],
                               loc_data['range']['dz']).tolist()
        except KeyError:
            if len(uniq_x) == len(uniq_y):
                locations['coordx'] = uniq_x
                locations['coordy'] = uniq_y
                locations['number'] = len(uniq_x)
            else:
                min_len = min(uniq_x, uniq_y)
                warn_msg = ''.join(['Lengths of the resulting location ',
                                    'arrays for x- and y-coordinates for ',
                                    'component {} are not the same. The arrays ',
                                    'will be truncated to the shorter ',
                                    'length of {}.']). format(cmpnt_name, min_len)
                logging.warning(warn_msg)
                locations['coordx'] = uniq_x[0:min_len]
                locations['coordy'] = uniq_y[0:min_len]
                locations['number'] = min_len
        else:
            if len(uniq_x) == len(uniq_y) == len(uniq_z):
                locations['coordx'] = uniq_x
                locations['coordy'] = uniq_y
                locations['coordz'] = uniq_z
                locations['number'] = len(uniq_x)
            else:
                min_len = min(uniq_x, uniq_y, uniq_z)
                warn_msg = ''.join(['Lengths of the resulting location ',
                                    'arrays for x-, y- and z-coordinates for ',
                                    'component {} are not the same. The arrays ',
                                    'will be truncated to the shorter ',
                                    'length of {}.']). format(cmpnt_name, min_len)
                logging.warning(warn_msg)
                locations['coordx'] = uniq_x[0:min_len]
                locations['coordy'] = uniq_y[0:min_len]
                locations['coordz'] = uniq_z[0:min_len]
                locations['number'] = min_len

    debug_msg = 'Known locations for component {} are {}.'.format(
        cmpnt_name, locations)
    logging.debug(debug_msg)

    return locations


def coord_locations(cmpnt_name, loc_data):
    """
    Get locations provided through coordinates.

    Parameters
    ----------
    :param cmpnt_name: name of component model
    :type name: str

    :param loc_data: locations data
    :type name: dict

    Returns
    -------
    dictionary locations with coordinates

    """
    locations = {}
    # Check whether all the needed parameters are present
    required_keys = ['coordx', 'coordy']
    key_present = {'coordx': 'coordy', 'coordy': 'coordx'}
    missing_keys = check_missing_keys(loc_data, required_keys)
    if len(missing_keys) == 1:  # if one of the keys are missing
        err_msg = "".join(["Parameter '{}' is provided but parameter '{}' is ",
                           "missing in the Locations setup of component {}."]).format(
                               key_present[missing_keys[0]], missing_keys[0],
                               cmpnt_name)
        logging.error(err_msg)
        raise KeyError(err_msg)
    if len(missing_keys) == 2:
        debug_msg = 'No x- and y-coordinates are provided for component {}.'.format(cmpnt_name)
        logging.debug(debug_msg)
    else:
        locations['coordx'] = loc_data['coordx']
        locations['coordy'] = loc_data['coordy']

        try:
            locations['coordz'] = loc_data['coordz']
        except KeyError:
            if len(locations['coordx']) == len(locations['coordy']):
                locations['number'] = len(locations['coordx'])
            else:
                min_len = min(locations['coordx'], locations['coordy'])
                locations['coordx'] = locations['coordx'][0:min_len]
                locations['coordy'] = locations['coordy'][0:min_len]
                locations['number'] = min_len
                warn_msg = ''.join([
                    'Lengths of the provided lists of x- and ',
                    'y-coordinates for component {} are not the same. ',
                    'The lists will be truncated to the shorter length ',
                    'of {}.']).format(cmpnt_name, min_len)
                logging.warning(warn_msg)
        else:
            if len(locations['coordx']) == len(locations['coordy']) == len(locations['coordz']):
                locations['number'] = len(locations['coordx'])
            else:
                min_len = min(locations['coordx'], locations['coordy'], len(locations['coordz']))
                locations['coordx'] = locations['coordx'][0:min_len]
                locations['coordy'] = locations['coordy'][0:min_len]
                locations['coordz'] = locations['coordz'][0:min_len]
                locations['number'] = min_len
                warn_msg = ''.join([
                    'Lengths of the provided lists of x-, y- and ',
                    'z-coordinates for component {} are not the same. ',
                    'The lists will be truncated to the shorter length ',
                    'of {}.']).format(cmpnt_name, min_len)
                logging.warning(warn_msg)

    debug_msg = 'Known locations for component {} are {}.'.format(
        cmpnt_name, locations)
    logging.debug(debug_msg)

    return locations


def raise_missing_keys_error(cmpnt_name, option, missing_pars):
    """
    Raise error if some or all required parameters of the setup of Locations
    are missing.

    Parameters
    ----------
    cmpnt_name : str
        name of component for which setup of locations
    option : str
        one of the possible ways to setup locations ('range', 'grid', etc.)
    missing_pars : list
        parameter names missing in the intended setup of a chosen option

    Raises
    ------
    KeyError
        Provides information what component, and what parameters are missing
        for the chosen option

    Returns
    -------
    None.

    """
    err_msg = "".join(["Required parameters {} are not provided ",
                       "in the Locations '{}' setup of the component {}."]).format(
                           missing_pars, option, cmpnt_name)
    logging.error(err_msg)
    raise KeyError(err_msg)


def check_missing_keys(data_dict, required_keys):
    """
    Check whether dictionary contains all required keys.

    Returns list of missing keys or empty list if all are present.
    """
    missing_keys = []

    for key in required_keys:
        if key not in data_dict:
            missing_keys.append(key)

    return missing_keys


def process_cell_centers(comp_data, comp_model, locations):
    """
    Analyze data provided in the setup of Seal Horizon component and
    obtain information about coordinates of cell centers.
    """
    try:
        cell_data = comp_data['Cells']
    except KeyError:
        msg = "".join(["Required keyword 'Cells' is missing in the setup ",
                       "of the SealHorizon component {}."]).format(comp_model)
        logging.error(msg)
        raise KeyError(msg) from None

    # Determine number of cells
    try:
        num_cells = cell_data['Number']
    except KeyError:
        msg = "".join(["Keyword 'Number' is missing in the setup of ",
                       "'Cells' for SealHorizon component {}."]).format(
            comp_model)
        logging.error(msg)
        raise KeyError(msg) from None

    # Determine locations of cells
    try:
        _ = cell_data['Locations']
    except KeyError:
        msg = "".join(["Keyword 'Locations' is missing in the setup of ",
                       "'Cells' for SealHorizon component {}."]).format(
            comp_model)
        logging.error(msg)
        raise KeyError(msg) from None
    else:
        cell_locations = get_known_locations(comp_model, cell_data)
        if cell_locations:
            locX = cell_locations['coordx']
            locY = cell_locations['coordy']
        else:
            msg = "".join(["Cell centers coordinates are not provided ",
                           "in the setup of SealHorizon component {}.",
                           "For example, keywords 'coordx'/'coordy' or 'file' ",
                           "can be provided"]).format(comp_model)
            logging.error(msg)
            raise KeyError(msg)

    # Setup cell centers
    if num_cells == len(locX):
        locations[comp_model] = {'number': len(locX),
                                 'coordx': locX,
                                 'coordy': locY,
                                 # Indicate whether component is a wellbore component
                                 'well': False}
    else:
        msg = "".join(["Length of array with cell centers coordinates is not equal ",
                       "to the value provided in the keyword 'Number'."])
        logging.error(msg)
        raise ValueError(msg)

    return locations


def process_fault_segment_centers(comp_data, comp_model, locations):
    """ Obtain information about coordinates of the fault segment centers."""
    try:
        segm_data = comp_data['Segments']
    except KeyError:
        msg = "".join(["Required keyword 'Segments' is missing in the setup ",
                       "of the FaultFlow component {}."]).format(comp_model)
        logging.error(msg)
        raise KeyError(msg) from None

    # Determine number of segments
    try:
        num_segments = segm_data['Number']
    except KeyError:
        msg = "".join(["Keyword 'Number' is missing in the setup of ",
                       "'Segments' for FaultFlow component {}."]).format(
            comp_model)
        logging.error(msg)
        raise KeyError(msg) from None

    # Determine locations of fault segments centers
    try:
        _ = segm_data['Locations']
    except KeyError:
        if 'DynamicParameters' not in comp_data:
            msg = "".join(["Keyword 'Locations' is missing in the setup of ",
                           "'Segments' for FaultFlow component {}. ",
                           "Locations will be calculated based on values ",
                           "of strike and dip."]).format(comp_model)
            logging.info(msg)
        # Get parameters data
        par_data = comp_data['Parameters']
        # Get strike, dip and length
        par_vals = {'strike': 30.0, 'length': 1000.0, 'xStart': 500.0, 'yStart': 300.0}
        for nm in par_vals:
            if nm in par_data:
                if not isinstance(par_data[nm], dict):
                    par_vals[nm] = par_data[nm]
                else:
                    par_vals[nm] = par_data[nm]['value']

        # Initialize list of locations
        locX = [par_vals['xStart']]
        locY = [par_vals['yStart']]
        if num_segments > 1:
            seg_length = par_vals['length']/num_segments
            for _ in range(1, num_segments):
                new_x = locX[-1] + seg_length*math.sin(math.radians(par_vals['strike']))
                new_y = locY[-1] + seg_length*math.cos(math.radians(par_vals['strike']))
                locX.append(new_x)
                locY.append(new_y)
    else:
        segm_locations = get_known_locations(comp_model, segm_data)
        if segm_locations:
            locX = segm_locations['coordx']
            locY = segm_locations['coordy']
        else:
            msg = "".join(["Coordinates of the fault segment centers are not provided ",
                           "in the setup of FaultFlow component {}. ",
                           "For example, keywords 'coordx'/'coordy' or 'file' ",
                           "can be provided."]).format(comp_model)
            logging.error(msg)
            raise KeyError(msg)


    # Setup fault segments centers
    if num_segments == len(locX):
        locations[comp_model] = {'number': len(locX),
                                 'coordx': locX,
                                 'coordy': locY,
                                 # Indicate whether component is a wellbore component
                                 'well': False}
    else:
        msg = "".join(["Length of array(s) with fault segment centers coordinates ",
                       "provided through 'Locations' setup for component {} ",
                       "is not equal to the value provided in the keyword ",
                       "'Number'."]).format(comp_model)
        logging.error(msg)
        raise ValueError(msg)

    return locations


def process_reservoir_locations(comp_data, comp_model,
                                locations, inj_well_locations, comp_type):
    """
    Analyze and process reservoir location data.
    """
    # Сheck if locations are provided
    if 'Locations' in comp_data:
        if comp_data['Locations']:
            point_locations = get_known_locations(comp_model, comp_data)
            locations[comp_model] = {'number': point_locations['number'],
                                     'coordx': point_locations['coordx'],
                                     'coordy': point_locations['coordy'],
                                     # Indicate whether component is a wellbore component
                                     'well': False}
            try:
                locations[comp_model]['coordz'] = point_locations['coordz']
            except KeyError:
                pass

        else:
            err_msg = "".join([
                "Parameter 'Locations' is used but no information ",
                "about coordinates is provided."])
            logging.error(err_msg)
            raise KeyError(err_msg)

    if comp_type in ['SimpleReservoir', 'AnalyticalReservoir', 'GenericReservoir',
                     'TheisReservoir']:
        if 'InjectionWell' in comp_data:
            if comp_data['InjectionWell']:
                inj_well_locations[comp_model] = [
                    comp_data['InjectionWell']['coordx'],
                    comp_data['InjectionWell']['coordy']]
            else:
                warn_msg = "".join([
                "Parameter 'InjectionWell' is used but no information ",
                "about coordinates is provided. Default location of injection ",
                "well at (0, 0) will be used."])
                inj_well_locations[comp_model] = [0.0, 0.0]
                logging.warning(warn_msg)
                raise Warning(warn_msg)

    return locations, inj_well_locations


def process_wellbore_locations(comp_data, comp_model, locations):
    """
    Analyze and process wellbore location data.
    """
    # Initialize dictionary key corresponding to the given component
    # Setup coordinates to empty lists in case we can get wellbore locations
    # from at least randomly generating them
    locations[comp_model] = {'coordx': [],
                             'coordy': [],
                             # we'll remove coordz key at the end if the list stays empty
                             'coordz': [],
                             # Indicate whether component is a wellbore component
                             'well': True}

    # Check if number of locations is provided
    if 'Locations' in comp_data:
        known_locations = get_known_locations(comp_model, comp_data)
        if known_locations:
            locations[comp_model]['coordx'] = known_locations['coordx']
            locations[comp_model]['coordy'] = known_locations['coordy']
            locations[comp_model]['number'] = known_locations['number']
            try:
                locations[comp_model]['coordz'] = known_locations['coordz']
            except KeyError:
                pass

    # Check if number of locations is provided
    if 'number' in comp_data:
        locations[comp_model]['number'] = comp_data['number']
    elif 'Number' in comp_data:
        locations[comp_model]['number'] = comp_data['Number']
    else:
        if not locations[comp_model]['coordx']:  # still empty list of x-coordinates
            err_msg = ''.join([
                "It is not possible to determine the total number ",
                "of wellbore locations: parameter 'Number' is not provided",
                " for component {}. Locations are not provided either."]).format(
                    comp_model)
            logging.error(err_msg)
            raise KeyError(err_msg)
        # Even when there are known locations provided but no total number of wells
        # is given we cannot determine how many random locations we need.
        # So we ignore data provided in RandomLocDomain setup and let user know about it
        if 'RandomLocDomain' in comp_data:
            warn_msg = ''.join([
                'Information provided in RandomLocDomain section ',
                'of the wellbore component {} is not used since ',
                'the total number of wells is not specified.']).format(comp_model)
            logging.warning(warn_msg)

    # Now two cases are possible:
    #    1) total number of wells is determined by 'number' parameter
    #    2) total number of wells is determined through length of known well locations list
    # In the second case we don't need to do anything.

    # Save number of locations
    num_locs = locations[comp_model]['number']

    # Check the total needed number of locations against the ones we get from
    # the known locations
    # If the total number of wells is less the number of coordinates provided
    # we use only the requested number of locations
    if num_locs < len(locations[comp_model]['coordx']):
        warn_msg = ''.join([
            "Number of provided locations (through x- and y-coordinates) ",
            "exceeds value specified in parameter 'number'. Only the ",
            "first {} locations will be used for simulation."]).format(num_locs)
        logging.warning(warn_msg)
        locations[comp_model]['coordx'] = locations[comp_model][
            'coordx'][0:num_locs]
        locations[comp_model]['coordy'] = locations[comp_model][
            'coordy'][0:num_locs]

        # Check presence of z-coordinates
        try:
            locations[comp_model]['coordz'] = locations[comp_model][
                'coordz'][0:num_locs]
        except KeyError:
            pass

    # If the total number of needed wells exceeds number of known wells we can try
    # to use random domain
    elif num_locs > len(locations[comp_model]['coordx']):
        if 'RandomLocDomain' in comp_data:
            num_rand_locs = num_locs-len(locations[comp_model]['coordx'])
            rand_locations = get_random_locations(comp_model, comp_data, num_rand_locs)
            locations[comp_model]['coordx']+= rand_locations['coordx']
            locations[comp_model]['coordy']+= rand_locations['coordy']
            try:
                locations[comp_model]['coordz']+= rand_locations['coordz']
            except KeyError:
                pass
        elif locations[comp_model]['coordx']:
            warn_msg = ''.join([
                "Number of provided locations (through x- and ",
                "y-coordinates) is less than the specified value ",
                "in parameter 'Number'. Only provided locations ",
                "will be used."])
            locations[comp_model]['number'] = len(
                locations[comp_model]['coordx'])
        else:
            err_msg = 'No locations data is provided for component {}'.format(
                comp_model)
            logging.error(err_msg)
            raise KeyError(err_msg)

    # If list with z-coordinates stayed empty we delete it from the dictionary of locations
    if not locations[comp_model]['coordz']:
        locations[comp_model].pop('coordz')

    return locations
