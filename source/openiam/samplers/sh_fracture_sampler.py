# -*- coding: utf-8 -*-
# @author: Veronika Vasylkivska, Ernest Lindner
# July, 2022
import os
import sys
import logging
import numpy as np

source_folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(source_folder)
sys.path.append(os.path.join(source_folder, 'components', 'seal'))


try:
    from openiam import SystemModel, SamplerModel
except ImportError as err:
    print('Unable to load IAM class module: '+ str(err))

try:
    import components.seal.frac_decipher as fde
    import components.seal.frac_origin as fog
    import components.seal.frac_stochastic as fstoch
    import components.seal.frac_user as fuser
    import components.seal.seal_config as sconf
    import components.seal.seal_intro as sintro
    import components.seal.seal_model as smod
    import components.seal.seal_units as sunit
    from components.seal.frac_setup import FRAC_SETUP_DICT
    from components.seal.seal_setup import SEAL_SETUP_DICT
except ImportError:
    print('\nERROR: Unable to load Fracture Sampler for Seal Horizon\n')
    sys.exit()


class SHFractureSampler(SamplerModel):
    """
    Class for model providing fractured seal permeability for Seal Horizon component.

    The class can be used as provider of permeability values for any other component.

    Setup parameters for fractured seal permeability

    * **connectFactor** [-] (0 to 1.01) - a reduction factor to consider vertical
      connectivity of fractures as well as the effects of fracture roughness
      and tortuosity with cubic law (default: 1)

    * **aveDensity** [number of fractures per |m^2|] (1.0e-10 to 100) - mean
      of the triangular distribution for fracture areal density (default: 2.0e-4)

    * **minDensity** [number of fractures per |m^2|] (1.0e-12 to 50) - minimum
      of the triangular distribution for fracture areal density (default: 1.0e-4)

    * **maxDensity** [number of fractures per |m^2|] (1.0e-10 to 500) - maximum
      of the triangular distribution for fracture areal density (default: 3.0e-4)

    * **aveTrend** [|deg|] (-180 to 180) - mean orientation/strike from North
      of the von Mises distribution of fractures (default: 30)

    * **spreadTrend** [|deg|] (0 to 180) - approximate 2-sigma spread in strike
      orientation of the random fractures from the mean, using a von Mises
      distribution (default: 20)

    * **lengthFunction** [-] (lognormal, powerlaw) - random distribution used
      to generate length value. Either of two types are allowed: lognormal,
      power law.

    * **aveLength** [|m|] (0.1 to 1000) - mean of the log-normal distribution
      for length. Parameter is used only if lengthFunction is lognormal (default: 50)

    * **stdDevLength** [|m|] (1.0e-3 to 50) - standard deviation of the log-normal
      distribution for length. Parameter is used only if lengthFunction is
      lognormal (default: 20)

    * **expLength** [-] (-5 to 5) - exponent used in power law representation
      length. Parameter is used only if lengthFunction is powerlaw (default: -0.8)

    * **minLength** [|m|] (0.1 to 50) - minimum length; this value censors the
      distribution and limits lower values (default: 0.1)

    * **maxLength** [|m|] (1 to 5000) - maximum length; this value censors the
      distribution and limits higher values (default: 500)

    * **aveBaseDepth** [|m|] (800 to 9500) - average depth to base of
      cells (default: 1100)

    * **aveBasePressure** [|Pa|] (1.0e+6 to 6.0e+7) - average pressure
      at cells base (default: 3.3e+7)

    * **staticDepth** [|m|] (800 to 9500) - reference depth of the top of cell
      (default: 1000)

    * **staticPressure** [|Pa|] (1.0e+6 to 6.0e+7) - pressure at static
      reference depth at the cell top (default: 1.0e+7).

    The following four parameters are used only if aperture_length_correlation
    attribute is set to False.

    * **aveAperture** [|m|] (1.0e-7 to 5.0e-2) - mean of the log-normal
      distribution for fracture aperture. This parameter is also used
      as a reference for threshold pressure (default: 5.0e-7)

    * **stdDevAperture** [|m|] (0 to 1.0e-2) - standard deviation of the log-normal
      distributionc for fracture aperture (default: 1.0e-7)

    * **minAperture** [|m|] (0 to 1.0e-2) - minimum aperture; this value censors the
      aperture distribution and limits lower values (default: 1.0e-9)

    * **maxAperture** [|m|] (1.0e-6 to 1.0e-1) - maximum aperture; this value
      censors the aperture distribution and limits higher values (default: 5.0e-6)

    The following two parameters are used only if aperture_length_correlation
    attribute is set to True.

    * **alpha** [-] (0.4 to 1) - exponent for aperture correlation (default: 0.6)

    * **beta** [|m|] (1.0e-8 to 1.0e-3) - factor in power-law length correlation;
      in effect, beta is equal to the aperture of 1 meter long fracture
      in the correlation relation (default: 1.0e-7)

    * **refEntryPressure** [|Pa|] (1 to 5.0e+6) - reference threshold pressure
      used in a defined relationship with average aperture to compute local threshold
      pressures (default: 5000)

    * **resHydrAperture** [|m|] (1.0e-7 to 1.0e-4) - residual hydraulic
      aperture (default: 5.0e-6)

    * **maxHydrAperture** [|m|] (1.0e-6 to 1.0e-2) - maximum hydraulic aperture
      (default: 5.0e-5)

    * **stressLimit** [|Pa|] (1.0e+6 to 5.0e+7) - stress limit for nonlinear response
      (default: 2.5e+7)

    * **theta** [-] (0.1 to 10) - stiffness history factor (default: 2.5)

    The following two parameters are used only if user_fractures argument is
    provided (i.e., not None) and define additional input for user fractures.

    * **refAperture** [|m|] (1.0e-7 to 5.0e-2) - matrix reference aperture
      (default: 1.0e-5)

    * **refPressure** [|Pa|] (1 to 5.0e+6) - matrix reference pressure
      (default: 1.0e+6)

    The following parameters define the rock matrix permeability properties.

    * **aveRockPerm** [|m^2|] (1.0e-26 to 1.0e-14) - mean of the log-normal
    distribution for matrix permeability (default: 1.0e-26)

    * **stdDevRockPerm** [|m^2|] (0 to 1.0e-16) - standard deviation of the log-normal
    distribution for matrix permeability (default: 1.0e-28)

    * **minRockPerm** [|m^2|] (1.0e-28 to 1.0e-16) - minimum matrix permeability;
    this value truncates the results of the random distribution and limits
    lower values (default: 1.0e-27)

    * **maxRockPerm** [|m^2|] (1.0e-24 to 1.0e-15) - maximum matrix permeability;
    this value truncates the results of the random distribution and limits higher
    values (default: 3.0e-24)

    * **refRockPerm** [|m^2|] (1.0e-28 to 1.0e-14) - permeability reference
    value used in defined relationship with reference threshold pressure
    for local threshold pressure (default: 1.0e-26)

    * **refRockPressure** [|Pa|] (1 to 5.0e+7) - threshold pressure reference
    value used in defined relationship for local threshold pressure (default: 1.0e+6)

    The possible outputs from the SH Fracture Sampler component are permeability and
    correlated entry pressure. The names of the observations are of the form:

    * **permeability** [|m^2|] - equivalent permeability for defined area or cells

    * **entryPressure** [|Pa|] - threshold pressure values for defined area or cells
    """
    def __init__(self, name, parent, grid_shape, geometry, coordinates,
                 constr_type='array', reproducible=True, **kwargs):
        """
        Constructor method of SHFractureSampler class

        :param name: name of component model
        :type name: [str]

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel [object]

        :param grid_shape: shape of grid for which the output is produced;
            for constr_type='matrix' shape of output is grid_shape, i.e.,
            (grid_shape[0], grid_shape[1]);
            for constr_type='array' shape of output is (grid_shape[0]*grid_shape[1],)
        :type shape: tuple of one or two integers

        :param geometry: dictionary with one ('area') or 2 keys ('height' and 'width')
            defining size of the cell for which permeability is generated. In the
            case of 'area' key provided height is assumed to be equal to width, and
            is calculated as square root of area.
        :type geometry: dict

        :param constr_type: type of the gridded observation. Possible values:
            'array', 'matrix'. The parameter indicates what type of gridded
            observation, a given component will provide. 'array' option means that
            the returned observation obs will be an array-like object and
            every element can be accessed as obs[ind]; 'matrix' means that the returned
            observation is a matrix-like object and every element can be accessed
            by providing 2 (or 3) indices of elements: either like
            obs[ind1][ind2] or obs[ind1,ind2] for 2d case, or
            obs[ind1][ind2][ind3] or obs[ind1,ind2,ind3] for 3d case.
        :type constr_type: str

        :param coordinates: x, y (and z) coordinates of the grid associated
            with the observation. coordinates is be a dictionary
            with keys 'x', 'y' (and 'z'). coordinates[key] should be an array-like object
            for any key in coordinates of the same shape as grid_shape
            if constr_type is 'matrix' and of length equal to the produce of
            dimensions if constr_type is 'array'.
        :type coordinates: dict of array-like objects

        :param kwargs: dictionary of additional arguments. Possible keys:
            'random_fractures' - indicates whether random fractures are to be
            generated; possible values are True and False; default value is True
            'user_fractures' - user provided fractures defined as array-like object
            of shape (n, 5) where n is a total number of fractures. The elements
            of each row are x_start, y_start, x_end, y_end, aperture_value in m,
            i.e, coordinates of the start and end of fracture,and its aperture.
            Default value of argument user_fractures is None
            'aperture_length_correlation' - indicates whether there is a aperture
            length correlation; possible values are True and False;
            default value is False
            'aperture_pressure_correlation' - indicates whether there is a pressure
            aperture correlation; possible values are True and False;
            default value is False
        :type kwargs: dict

        :returns: SHFractureSampler class object
        """
        super().__init__(name=name, parent=parent, reproducible=reproducible)

        # Add type attribute
        self.class_type = 'SHFractureSampler'

        self.model_kwargs = {}

        # Define gridded observations names
        self.grid_obs_keys = ['permeability', 'entryPressure']

        # Set default parameters
        self.setup_default_pars()

        # Return updated dictionary with default fault setup
        # using user provided keyword arguments
        # Attribute user_fractures is also setup in the process_additional_input method
        fs_setup_dict = self.process_additional_input(**kwargs)

        # Define default controls
        self.controls = {}
        self.controls = fde.interpret_frac_parameters(fs_setup_dict, self.controls)

        # Define parameters boundaries dictionary self.pars_bounds
        self.define_pars_bounds()

        # Logging output of default parameters and their values
        debug_msg = ''.join([par_nm+': '+str(par_obj.value) \
            for par_nm, par_obj in self.default_pars.items()])
        logging.debug(debug_msg)

        # Define additional controls
        self.process_geometry(grid_shape, geometry)
        self.check_coordinates(coordinates, constr_type)

        # Create grid of cells for compatibility with standalone code
        self.create_grid()

        # Add gridded observation to be linked
        # It makes sense to add it here since the model has only one observation
        # and it's assumed that if the object of SHFractureSampler class is
        # created then it's observation will be used by some other component.
        # Otherwise, beyond test purposes there are no other reasons to create
        # the object
        for obs_nm in self.grid_obs_keys:
            self.add_obs_to_be_linked(
                obs_nm, obs_type='grid', constr_type=self.controls['constr_type'],
                coordinates=coordinates, index=[0])

    def process_additional_input(self, **kwargs):
        """
        Process additional keyword arguments provided to the constructor method.

        :param kwargs: dictionary of additional arguments. Possible keys:
            'random_fractures' - indicates whether random fractures are to be
            generated; possible values are True and False; default value is True
            'user_fractures' - user provided fractures defined as array-like object
            of shape (n, 5) where n is a total number of fractures. The elements
            of each row are x_start, y_start, x_end, y_end, aperture_value in m,
            i.e, coordinates of the start and end of fracture,and its aperture.
            Default value of argument user_fractures is None
            'aperture_length_correlation' - indicates whether there is a aperture
            length correlation; possible values are True and False;
            default value is False
            'aperture_pressure_correlation' - indicates whether there is a pressure
            aperture correlation; possible values are True and False;
            default value is False
        :type kwargs: dict
        """
        # Set attribute for user fractures
        self.user_fractures = kwargs.get('user_fractures', None)

        # Get default setup
        setup_dict = FRAC_SETUP_DICT.copy()

        # Update some of the setup based on user input
        if self.user_fractures is not None:
            setup_dict['Controls']['userFracApproach'] = True

        # Setup controls
        setup_dict['Controls']['randomFracApproach'] = kwargs.get(
            'random_fractures', True)
        setup_dict['Controls']['apertureLengthApproach'] = kwargs.get(
            'aperture_length_correlation', False)
        setup_dict['Controls']['pressureApproach'] = kwargs.get(
            'aperture_pressure_correlation', False)

        return setup_dict

    def check_coordinates(self, coordinates, constr_type):
        """
        Check whether coordinates argument of constructor method satisfies the
        required conditions.

        :param coordinates: dictionary with keys 'x', 'y', and possibly 'z';
            each item is an array-like object
        :type coordinates: dict

        :param constr_type: type of structure of requested output; one of the two
            values are possible: 'array' or 'matrix'
        :type constr_type: str
        """
        if constr_type == 'array':
            for key, coord_arr in coordinates.items():
                if len(coord_arr) != self.controls['num_cells']:
                    err_msg = ''.join(
                        ['Length of the provided {}-coordinates ({}) is not equal ',
                         'to the total number of cells ({}).']).format(
                             key, len(coord_arr), self.controls['num_cells'])
                    logging.error(err_msg)
                    raise IndexError(err_msg)
        elif constr_type == 'matrix':
            for key, coord_matr in coordinates.items():
                if coord_matr.shape[0] != self.controls['grid_rows'] or \
                        coord_matr.shape[1] != self.controls['grid_cols']:
                    err_msg = ''.join(
                        ['Shape of the provided {}-coordinates is ({}, {}) ',
                         'but should be ({}, {}).']).format(
                             key, coord_matr.shape[0], coord_matr.shape[1],
                             self.controls['grid_rows'], self.controls['grid_cols'])
                    logging.error(err_msg)
                    raise IndexError(err_msg)
        else:
            err_msg = 'Value {} of argument constr_type is not valid'.format(
                constr_type)
            logging.error(err_msg)
            raise ValueError(err_msg)

        # Save constructio type and coordinates of the cell centers
        self.controls['constr_type'] = constr_type
        self.coordinates = {key: coords.flatten() for key, coords in coordinates.items()}

    def process_geometry(self, grid_shape, geometry):
        """
        Process constructor method input related to the area of cells and
        desired shape of output

        :param grid_shape: shape of grid for which the output is produced;
            for constr_type='matrix' shape of output is grid_shape, i.e.,
            (grid_shape[0], grid_shape[1]);
            for constr_type='array' shape of output is (grid_shape[0]*grid_shape[1],)
        :type shape: tuple of one or two integers

        :param geometry: dictionary with one ('area') or 2 keys ('height' and 'width')
            defining size of the cell for which permeability is generated. In the
            case of 'area' key provided height is assumed to be equal to width, and
            is calculated as square root of area.
        :type geometry: dict
        """
        if len(grid_shape) == 2: # length is 2
            self.controls['grid_approach'] = True
            self.controls['grid_rows'] = grid_shape[0]
            self.controls['grid_cols'] = grid_shape[1]
            self.controls['num_cells'] = grid_shape[0]*grid_shape[1]
        else:
            err_msg = ''.join([
                'Argument grid_shape is of wrong length. It should have ',
                'length 2 but has length {}.']).format(len(grid_shape))
            logging.error(err_msg)
            raise ValueError(err_msg)

        if 'area' in geometry:
            self.controls['area'] = geometry['area']
            self.controls['cell_height'] = self.controls['area']**0.5
            self.controls['cell_width'] = self.controls['cell_height']
        elif 'height' in geometry and 'width' in geometry:
            self.controls['cell_height'] = geometry['height']
            self.controls['cell_width'] = geometry['width']
            self.controls['area'] = self.controls['cell_height']*self.controls['cell_width']
        else:
            err_msg = ''.join(['Not enough data is provided for the geometry ',
                               '(area or height and width) of cells ',
                               'for fracture sampler {}']).format(self.name)
            logging.error(err_msg)
            raise KeyError(err_msg)

    def create_grid(self):
        """
        Create grid of cells and assign corresponding attributes.
        """
        # Create grid of cells
        self.grid = []
        for ind in range(self.controls['num_cells']):
            # Create a cell
            self.grid.append(
                smod.Cell(x_center=self.coordinates['x'][ind],
                          y_center=self.coordinates['y'][ind]))
            self.grid[-1].area = self.controls['area']

    def setup_default_pars(self):
        """
        Add parameters of the component with default values.
        """
        # Add default parameters specific to this sampler
        # Note that seed parameter is already parameter of the sampler
        # by definition of sampler
        # Connectivity and roughness factor for site
        self.add_default_par(
            'connectFactor', value=FRAC_SETUP_DICT['Connectivity']['geometric'])

        # Distribution parameters for random fractures
        self.add_default_par(
            'aveDensity',
            value=FRAC_SETUP_DICT['RandomFracs']['Density']['aveDensity']['valu'])
        self.add_default_par(
            'minDensity',
            value=FRAC_SETUP_DICT['RandomFracs']['Density']['minDensity']['valu'])
        self.add_default_par(
            'maxDensity',
            value=FRAC_SETUP_DICT['RandomFracs']['Density']['maxDensity']['valu'])

        # Reference parameters
        self.add_default_par(
            'aveBaseDepth',
            value=SEAL_SETUP_DICT['Description']['aveBaseDepth']['valu'])
        self.add_default_par(
            'aveBasePressure',
            value=SEAL_SETUP_DICT['Description']['aveBasePressure']['valu'])
        self.add_default_par(
            'staticDepth',
            value=SEAL_SETUP_DICT['Description']['staticDepth']['valu'])
        self.add_default_par(
            'staticPressure',
            value=SEAL_SETUP_DICT['Description']['staticPressure']['valu'])

        # Parameters for random fractures orientation
        self.add_default_par(
            'aveTrend',
            value=FRAC_SETUP_DICT['RandomFracs']['Orientation']['aveTrend'])
        self.add_default_par(
            'spreadTrend',
            value=FRAC_SETUP_DICT['RandomFracs']['Orientation']['spreadTrend'])

        # Parameters for random fractures length
        self.model_kwargs['lengthFunction'] = \
            FRAC_SETUP_DICT['RandomFracs']['Length']['function']
        self.add_default_par(
            'aveLength',
            value=FRAC_SETUP_DICT['RandomFracs']['Length']['aveLength']['valu'])
        self.add_default_par(
            'stdDevLength',
            value=FRAC_SETUP_DICT['RandomFracs']['Length']['stdDevLength']['valu'])
        self.add_default_par(
            'expLength',
            value=FRAC_SETUP_DICT['RandomFracs']['Length']['expLength']['valu'])
        self.add_default_par(
            'minLength',
            value=FRAC_SETUP_DICT['RandomFracs']['Length']['minLength']['valu'])
        self.add_default_par(
            'maxLength',
            value=FRAC_SETUP_DICT['RandomFracs']['Length']['maxLength']['valu'])

        # Parameters for random fractures apertures
        self.add_default_par(
            'aveAperture',
            value=FRAC_SETUP_DICT['RandomFracs']['VaryAperture']['aveAperture']['valu'])
        self.add_default_par(
            'stdDevAperture',
            value=FRAC_SETUP_DICT['RandomFracs']['VaryAperture']['stdDevAperture']['valu'])
        self.add_default_par(
            'minAperture',
            value=FRAC_SETUP_DICT['RandomFracs']['VaryAperture']['minAperture']['valu'])
        self.add_default_par(
            'maxAperture',
            value=FRAC_SETUP_DICT['RandomFracs']['VaryAperture']['maxAperture']['valu'])

        # Parameters for correlation
        self.add_default_par(
            'alpha',
            value=FRAC_SETUP_DICT['RandomFracs']['CorrelateAperture']['alpha'])
        self.add_default_par(
            'beta',
            value=FRAC_SETUP_DICT['RandomFracs']['CorrelateAperture']['beta']['valu'])

        # Pressure related parameters
        self.add_default_par(
            'refEntryPressure',
            value=FRAC_SETUP_DICT['RandomFracs']['Threshold']['refEntry']['valu'])
        self.add_default_par(
            'resHydrAperture',
            value=FRAC_SETUP_DICT['RandomFracs']['PressureCorrection']['resAperture']['valu'])
        self.add_default_par(
            'maxHydrAperture',
            value=FRAC_SETUP_DICT['RandomFracs']['PressureCorrection']['wideAperture']['valu'])
        self.add_default_par(
            'stressLimit',
            value=FRAC_SETUP_DICT['RandomFracs']['PressureCorrection']['stressLimit']['valu'])
        self.add_default_par(
            'theta',
            value=FRAC_SETUP_DICT['RandomFracs']['PressureCorrection']['theta'])

        # Parameters related to user input fractures
        self.add_default_par(
            'refAperture',
            value=FRAC_SETUP_DICT['InputFracs']['Threshold']['refAperture']['valu'])
        self.add_default_par(
            'refPressure',
            value=FRAC_SETUP_DICT['InputFracs']['Threshold']['refPressure']['valu'])

        # Rock matrix parameters
        self.add_default_par(
            'aveRockPerm',
            value=FRAC_SETUP_DICT['RockMatrix']['Permeability']['aveMatrix']['valu'])
        self.add_default_par(
            'stdDevRockPerm',
            value=FRAC_SETUP_DICT['RockMatrix']['Permeability']['stdDevMatrix']['valu'])
        self.add_default_par(
            'minRockPerm',
            value=FRAC_SETUP_DICT['RockMatrix']['Permeability']['minMatrix']['valu'])
        self.add_default_par(
            'maxRockPerm',
            value=FRAC_SETUP_DICT['RockMatrix']['Permeability']['maxMatrix']['valu'])
        self.add_default_par(
            'refRockPerm',
            value=FRAC_SETUP_DICT['RockMatrix']['Threshold']['matrixPerm']['valu'])
        self.add_default_par(
            'refRockPressure',
            value=FRAC_SETUP_DICT['RockMatrix']['Threshold']['matrixPress']['valu'])

    def define_pars_bounds(self):
        """
        Define dictionary of parameter boundaries.

        The method uses boundaries imposed in the standalone Seal Flux code.
        """
        # Get boundaries from frac_origin module for
        params_bounds = fog.define_param_limits()

        sh_params_bounds = sintro.define_input_limits()

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        # np.inf is used to ephasize that no upper limit is imposed
        self.pars_bounds['seed'] = [0, np.inf]

        # Connectivity
        self.pars_bounds['connectFactor'] = params_bounds['connect_factor']  # 'geometric'

        # Random Fracs - Density (per m2):
        self.pars_bounds['aveDensity'] = params_bounds['density_ave']
        self.pars_bounds['minDensity'] = params_bounds['density_min']
        self.pars_bounds['maxDensity'] = params_bounds['density_max']

        self.pars_bounds['aveBaseDepth'] = sh_params_bounds['ave_base_depth']
        self.pars_bounds['aveBasePressure'] = sh_params_bounds['ave_base_pressure']
        self.pars_bounds['staticDepth'] = sh_params_bounds['static_depth']
        self.pars_bounds['staticPressure'] = sh_params_bounds['static_pressure']

        # Random Fracs - Orientation (deg):
        self.pars_bounds['aveTrend'] = params_bounds['orient_mu']
        self.pars_bounds['spreadTrend'] = params_bounds['orient_sigma']

        # Random Fracs - Length (m):
        self.pars_bounds['aveLength'] = params_bounds['length_ave']
        self.pars_bounds['stdDevLength'] = params_bounds['length_dev']
        self.pars_bounds['expLength'] = params_bounds['length_eta']
        self.pars_bounds['minLength'] = params_bounds['length_min']
        self.pars_bounds['maxLength'] = params_bounds['length_max']

        # Random Fracs - Vary Aperture (mm):
        # Convert boundaries from mm to m
        mm_to_m = sunit.mm_to_m()
        self.pars_bounds['aveAperture'] = [
            val*mm_to_m for val in params_bounds['aperture_ave']]
        self.pars_bounds['stdDevAperture'] = [
            val*mm_to_m for val in params_bounds['aperture_dev']]
        self.pars_bounds['minAperture'] = [
            val*mm_to_m for val in params_bounds['aperture_min']]
        self.pars_bounds['maxAperture'] = [
            val*mm_to_m for val in params_bounds['aperture_max']]

        # Random Fracs - Correlate Aperture:
        self.pars_bounds['alpha'] = params_bounds['aperture_alpha']
        # Transform from mm to m
        self.pars_bounds['beta'] = [val*mm_to_m for val in params_bounds['aperture_beta']]

        # Random Threshold (Pa):
        self.pars_bounds['refEntryPressure'] = params_bounds['entry_pressure'] # 'refEntry'

        # Pressure-Aperture Correction:
        self.pars_bounds['resHydrAperture'] = [
            val*mm_to_m for val in params_bounds['residual_aperture']] # 'resAperture'
        self.pars_bounds['maxHydrAperture'] = [
            val*mm_to_m for val in params_bounds['wide_aperture']]     # 'wideAperture'
        self.pars_bounds['stressLimit'] = params_bounds['stress_limit']
        self.pars_bounds['theta'] = params_bounds['theta_aperture']

        # Reference aperture:
        self.pars_bounds['refAperture'] = [
            val*mm_to_m for val in params_bounds['ref_user_aperture']]
        self.pars_bounds['refPressure'] = params_bounds['ref_user_pressure']

        # Matrix Permeability:
        mcD_to_m2 = sunit.microd_to_metersq()
        self.pars_bounds['aveRockPerm'] = [
            val*mcD_to_m2 for val in params_bounds['rock_perm_ave']]    # 'aveMatrix'
        self.pars_bounds['stdDevRockPerm'] = [
            val*mcD_to_m2 for val in params_bounds['rock_perm_dev']]    # 'stdDevMatrix'
        self.pars_bounds['minRockPerm'] = [
            val*mcD_to_m2 for val in params_bounds['rock_perm_min']]    # 'minMatrix'
        self.pars_bounds['maxRockPerm'] = [
            val*mcD_to_m2 for val in params_bounds['rock_perm_max']]    # 'maxMatrix'

        # Matrix Reference:
        self.pars_bounds['refRockPerm'] = [
            val*mcD_to_m2 for val in params_bounds['ref_matrix_perm']]  # 'matrixPerm'
        self.pars_bounds['refRockPressure'] = params_bounds['ref_matrix_threshold'] # 'matrixPress'

    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        debug_msg = 'Input parameters of component {} are {}'.format(self.name, p)
        logging.debug(debug_msg)

        # Check the rest of the sampler parameters
        for key in self.default_pars.keys():
            if key in p:
                val = p[key]
                if ((val < self.pars_bounds[key][0]) or
                        (val > self.pars_bounds[key][1])):
                    warn_msg = ''.join([
                        'Parameter {} of value {} of SH Fracture Sampler {} ',
                        'is out of boundaries [{}, {}].']).format(
                            key, val, self.name, self.pars_bounds[key][0],
                            self.pars_bounds[key][1])
                    logging.warning(warn_msg)

    def update_controls(self, p, **kwargs):
        """
        Update self.controls attribute with the user provided data.
        """
        # Update controls related to user parameters
        self.controls['connect_factor'] = p['connectFactor']
        self.controls['density_ave'] = p['aveDensity']
        self.controls['density_min'] = p['minDensity']
        self.controls['density_max'] = p['maxDensity']
        self.controls['orient_mu'] = p['aveTrend']
        self.controls['orient_sigma'] = p['spreadTrend']

        check_val = kwargs.get('lengthFunction', 'lognormal').upper()
        if check_val == 'LOGNORMAL':
            self.controls['length_approach'] = 'LOGNORM'
            self.controls['length_ave'] = p['aveLength']
            self.controls['length_dev'] = p['stdDevLength']
        else:
            self.controls['length_approach'] = 'POWER'
            self.controls['length_eta'] = p['expLength']

        self.controls['length_min'] = p['minLength']
        self.controls['length_max'] = p['maxLength']

        m_to_mm = sunit.m_to_mm()
        self.controls['aperture_ave'] = p['aveAperture']*m_to_mm
        self.controls['aperture_dev'] = p['stdDevAperture']*m_to_mm
        self.controls['aperture_min'] = p['minAperture']*m_to_mm
        self.controls['aperture_max'] = p['maxAperture']*m_to_mm
        self.controls['aperture_alpha'] = p['alpha']
        self.controls['aperture_beta'] = p['beta']*m_to_mm
        self.controls['entry_pressure'] = p['refEntryPressure']

        self.controls['residual_aperture'] = p['resHydrAperture']
        self.controls['wide_aperture'] = p['maxHydrAperture']
        self.controls['stress_limit'] = p['stressLimit']
        self.controls['theta_aperture'] = p['theta']
        self.controls['static_depth'] = p['staticDepth']
        self.controls['static_pressure'] = p['staticPressure']
        self.controls['ave_base_depth'] = p['aveBaseDepth']
        self.controls['ave_base_pressure'] = p['aveBasePressure']

        self.controls['ref_user_aperture'] = p['refAperture']
        self.controls['ref_user_pressure'] = p['refPressure']

        m2_to_mcD = sunit.metersq_to_microd()
        self.controls['rock_perm_ave'] = p['aveRockPerm']*m2_to_mcD
        self.controls['rock_perm_dev'] = p['stdDevRockPerm']*m2_to_mcD
        self.controls['rock_perm_min'] = p['minRockPerm']*m2_to_mcD
        self.controls['rock_perm_max'] = p['maxRockPerm']*m2_to_mcD
        self.controls['ref_matrix_perm'] = p['refRockPerm']*m2_to_mcD
        self.controls['ref_matrix_threshold'] = p['refRockPressure']

        # For compatibility with standalone
        self.controls['plot_fractures'] = False

    def check_pressure(self):
        """Check if input for pressure is OK.

        Method is based on the same method from standalone Seal Flux code.
        """
        # Get minimum lithostatic pressure (in Pa) & input pressure.
        depth = self.controls['ave_base_depth']
        min_pressure = sunit.brine_pressure(depth)
        input_pressure = self.controls['ave_base_pressure']

        # Check pressure and report if in error.
        if input_pressure < min_pressure:
            # Error message.
            err_msg = ''.join(['Parameter aveBasePressure of value {} is less than ',
                           'the estimate minimum {} from gradient-based ',
                           'calculations.']).format(input_pressure, min_pressure)
            logging.error(err_msg)
            raise ValueError(err_msg)

    def process_random_approach(self):
        """
        Process additional parameters if random approach is chosen.

        Code based on the part of standalone code of Seal Horizon
        in frac_origin.py module for method frac_launch
        """
        # Define other parameters only if random approach
        if self.controls['random_approach']:
            # Setup log-normal parameters
            self.controls = fog.setup_aperture_lognormal(self.controls)

            # Convert sigma to kappa for von Mises distribution
            self.controls = fog.setup_kappa(self.controls)

            # Define orientation of second fracture set
            self.controls['orient_set_2'] = self.controls['orient_mu'] + 90.0
            if self.controls['orient_set_2'] > 360.0:
                self.controls['orient_set_2'] -= 360.0

    def evaluate_fracs(self):
        """
        Evaluate the effect of fracturing effect on permeability and pressure.

        Method is based on the standalone code Seal Horizon method evaluate_fracs
        in the module frac_origin.py
        """
        # Setup zero permeability arrays
        random_perm, user_perm, matrix_perm = fog.create_perm_storage(self.controls)

        # Setup zero threshold pressure arrays
        random_threshold, user_threshold, matrix_threshold = \
            fog.create_entry_storage(self.controls)

        # Setup list for combined output
        combo_list = [[], [], []]

        # Evaluate random fractures
        if self.controls['random_approach']:
            # Generate random fractures and calculate permeabilities for grid
            random_perm, random_threshold, combo_list[0] = \
                fstoch.compute_random_permeability(self.controls, self.grid, False)

        # Evaluate user fractures
        if self.controls['user_approach']:
            # Input user fractures.
            user_data = self.user_fractures
            # Transform the last column of the data containing aperture values
            # from meters to millimeters units
            user_data[:, 4] = user_data[:, 4]*sunit.m_to_mm()

            # Process user fractures to obtain permeability
            user_perm, user_threshold, combo_list[1] = \
                fuser.process_user_lines(self.controls, self.grid, user_data)

        # Evaluate matrix characteristics
        # Use probability distribution to get matrix values
        matrix_perm, matrix_threshold, combo_list[2] = \
            fuser.define_matrix_perm(self.controls)

        # Compute final equivalent properties
        # Sum/process case arrays to get equivalent values for grid
        perm_sum = random_perm + user_perm + matrix_perm
        entry_sum = fuser.sum_threshold_arrays(self.controls, random_threshold,
                                               user_threshold, matrix_threshold)

        # Transform permeability from microDarcy to m^2
        mcD_to_m2 = sunit.microd_to_metersq()
        perm_sum = perm_sum*mcD_to_m2

        return perm_sum, entry_sum

    def sample(self, p, **kwargs):
        """
        :param p: input parameters of SHFractureSampler model
        :type p: dict

        :param kwargs: additional arguments of sample method. Possible keys:
            'lengthFunction' - either of two values: 'lognormal' and 'powerlaw';
            default value: 'lognormal'

        """
        if self.reproducible: # if repeatable results are needed
            if 'seed' in p:
                # check_seed_parameter is a default method of SamplerModel Class
                # checking whether seed is a positive integer
                self.check_seed_parameter(p['seed'])
                # Set seed for random numbers on each simulation
                rng = np.random.default_rng(p['seed'])
            else:
                rng = np.random.default_rng(sconf.SEEDX)
        else:
            rng = np.random.default_rng()

        # Save random number generator reference
        self.controls['rng'] = rng

        # Obtain the default values of the parameters from dictionary of
        # default parameters.
        actual_p = {k: v.value for k, v in self.default_pars.items()}

        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Set default if not provided
        lengthFunction = kwargs.get('lengthFunction', 'lognormal')
        if lengthFunction not in ['lognormal', 'powerlaw']:
            warn_msg = ''.join(['Parameter lengthFunction has a wrong value ',
                                'of {}. Default value of "lognormal" will ',
                                'be used.']).format(lengthFunction)
            kwargs['lengthFunction'] = 'lognormal'
            logging.warning(warn_msg)

        # Update control values
        self.update_controls(actual_p, **kwargs)

        # Check whether provided pressure parameters satisfy required conditions
        self.check_pressure()

        # Add matrix parameters for threshold & lognormal
        self.controls = fog.setup_threshold(self.controls)
        self.controls = fog.setup_matrix_lognormal(self.controls)

        # Define other parameters only if random approach
        self.process_random_approach()

        # Evaluate fractures to get permeability and entry pressure
        permeability, entry_pressure = self.evaluate_fracs()

        # Check whether matrix-type output is needed
        if self.controls['constr_type'] == 'matrix':
            permeability = permeability.reshape(self.controls['grid_rows'],
                                                self.controls['grid_cols'])
            entry_pressure = entry_pressure.reshape(self.controls['grid_rows'],
                                                self.controls['grid_cols'])

        # Define output
        out = {'permeability': permeability,
               'entryPressure': entry_pressure}

        # Return output
        return out


def test_scenario1():
    """ Run script example illustrating work of SH Permeability Sampler model.

    SH Permeability Sampler uses only randomly generated fractures.

    """
    # Setup logging and constants
    logging.basicConfig(level=logging.DEBUG)

    # Define keyword arguments of the system model
    model_kwargs = {'time_point': 0.0} # time is given in days

    output_dir = '../../../output/samplers/'
    # Create system model
    sm = SystemModel(model_kwargs=model_kwargs)
    nx = 10
    ny = 12
    dx = 30
    dy = 30
    # The domain size is approximately 300x350 m
    x = dx/2. + dx*np.array(range(nx))
    y = dy/2. + dy*np.array(range(ny))
    yy, xx = np.meshgrid(y, x)
    frac_sampler = sm.add_component_model_object(
        SHFractureSampler(name='ps', parent=sm, grid_shape=(nx, ny),
                          reproducible=False,
                          constr_type='matrix',
                          geometry={'width': dx, 'height': dy},
                          coordinates={'x': xx, 'y': yy}))
    frac_sampler.add_par('seed', value=234, vary=False)
    frac_sampler.add_par('connectFactor', value=1, vary=False)

    # Distribution parameters for random fractures
    frac_sampler.add_par('aveDensity', value=2.0e-4, vary=False)
    frac_sampler.add_par('minDensity', value=1.0e-4, vary=False)
    frac_sampler.add_par('maxDensity', value=3.0e-4, vary=False)

    # Reference parameters
    base_depth = 1500.0
    brine_density = 1070.0
    top_depth = 1200.0
    frac_sampler.add_par('aveBaseDepth', value=base_depth, vary=False)
    frac_sampler.add_par('aveBasePressure', vary=False,
                         value=101325+9.8*brine_density*base_depth)
    frac_sampler.add_par('staticDepth', value=top_depth, vary=False)
    frac_sampler.add_par('staticPressure',vary=False,
                         value=101325+9.8*brine_density*top_depth)

    # Parameters for random fractures orientation
    frac_sampler.add_par('aveTrend', value=45, vary=False)
    frac_sampler.add_par('spreadTrend', value=15, vary=False)

    # Parameters for random fractures length
    frac_sampler.model_kwargs['lengthFunction'] = 'lognormal'
    frac_sampler.add_par('aveLength', value=15, vary=False)
    frac_sampler.add_par('stdDevLength', value=5, vary=False)
    frac_sampler.add_par('minLength', value=10, vary=False)
    frac_sampler.add_par('maxLength', value=24, vary=False)

    # Parameters for random fractures apertures
    frac_sampler.add_par('aveAperture', value=2.0e-5, vary=False)
    frac_sampler.add_par('stdDevAperture', value=1.0e-5, vary=False)
    frac_sampler.add_par('minAperture', value=1.0e-6, vary=False)
    frac_sampler.add_par('maxAperture', value=5.0e-5, vary=False)

    # Pressure related parameters
    frac_sampler.add_par('refEntryPressure', value=1.0e+4, vary=False)
    frac_sampler.add_par('resHydrAperture', value=5.0e-5, vary=False)
    frac_sampler.add_par('maxHydrAperture', value=1.0e-5, vary=False)
    frac_sampler.add_par('stressLimit', value=2.0e+7, vary=False)
    frac_sampler.add_par('theta', value=2.0, vary=False)

    # Rock matrix parameters
    frac_sampler.add_par('aveRockPerm', value=1.0e-22, vary=False)
    frac_sampler.add_par('stdDevRockPerm', value=1.0e-21, vary=False)
    frac_sampler.add_par('minRockPerm', value=5.0e-23, vary=False)
    frac_sampler.add_par('maxRockPerm', value=1.0e-21, vary=False)
    frac_sampler.add_par('refRockPerm', value=1.0e-22, vary=False)
    frac_sampler.add_par('refRockPressure', value=1.0e+6)
    frac_sampler.add_grid_obs('permeability', constr_type='array', index=[0],
                              output_dir=output_dir)
    frac_sampler.add_grid_obs('entryPressure', constr_type='array', index=[0],
                              output_dir=output_dir)

    sm.forward()

    # Collect saved observation
    perm_data = sm.collect_gridded_observations_as_time_series(
        frac_sampler, 'permeability', output_dir,
        indices=[0], rlzn_number=0)

    print(perm_data[0].shape)
    # print(perm_data[0])

    pressure_data = sm.collect_gridded_observations_as_time_series(
        frac_sampler, 'entryPressure', output_dir,
        indices=[0], rlzn_number=0)

    print(pressure_data[0].shape)
    # print(pressure_data[0])

    from matplotlib import cm
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d

    zz = np.zeros(nx*ny,)

    dxx = dx*np.ones(nx*ny,)
    dyy = dy*np.ones(nx*ny,)
    dzz = perm_data[0].reshape(nx*ny,)
    cmap = cm.get_cmap('Spectral')
    max_height = np.max(dzz)
    min_height = np.min(dzz)
    colors = [cmap((k-min_height)/max_height) for k in dzz]
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.bar3d(xx.reshape(nx*ny,), yy.reshape(nx*ny,), zz,
              dxx, dyy, dzz, shade=True, color=colors)
    ax2.set_xlabel('x [m]')
    ax2.set_ylabel('y [m]')
    ax2.set_zlabel(r'k [$m^2$]')


def test_scenario2():
    """ Run script example illustrating work of SH Permeability Sampler model.

    SH Permeability Sampler uses only user provided fractures.
    """
    # Setup logging and constants
    logging.basicConfig(level=logging.DEBUG)

    # Define keyword arguments of the system model
    model_kwargs = {'time_point': 0.0} # time is given in days

    output_dir = '../../../output/samplers/'
    # Create system model
    sm = SystemModel(model_kwargs=model_kwargs)
    nx = 10
    ny = 12
    dx = 30
    dy = 30
    # The domain size is approximately 300x350 m
    x = dx/2. + dx*np.array(range(nx))
    y = dy/2. + dy*np.array(range(ny))
    yy, xx = np.meshgrid(y, x)
    known_fractures = np.array([[25., 50., 250., 50., 2.0e-5],      # fracture 1: x1, y1, x2, y2, aperture
                                [250., 100., 250., 300., 2.0e-5]])  # fracture 2: x1, y1, x2, y2, aperture
    frac_sampler = sm.add_component_model_object(
        SHFractureSampler(name='ps', parent=sm, grid_shape=(nx, ny),
                          reproducible=False,
                          constr_type='matrix',
                          geometry={'width': dx, 'height': dy},
                          coordinates={'x': xx, 'y': yy},
                          random_fractures=False,
                          user_fractures=known_fractures))
    frac_sampler.add_par('seed', value=234, vary=False)
    frac_sampler.add_par('connectFactor', value=1, vary=False)

    # Reference parameters
    base_depth = 1500.0
    brine_density = 1070.0
    top_depth = 1200.0
    frac_sampler.add_par('aveBaseDepth', value=base_depth, vary=False)
    frac_sampler.add_par('aveBasePressure', vary=False,
                         value=101325+9.8*brine_density*base_depth)
    frac_sampler.add_par('staticDepth', value=top_depth, vary=False)
    frac_sampler.add_par('staticPressure',vary=False,
                         value=101325+9.8*brine_density*top_depth)

    # Pressure related parameters
    frac_sampler.add_par('refEntryPressure', value=1.0e+4, vary=False)
    frac_sampler.add_par('resHydrAperture', value=5.0e-5, vary=False)
    frac_sampler.add_par('maxHydrAperture', value=1.0e-5, vary=False)
    frac_sampler.add_par('stressLimit', value=2.0e+7, vary=False)
    frac_sampler.add_par('theta', value=2.0, vary=False)

    frac_sampler.add_par('refAperture', value=1.0e-6, vary=False)
    frac_sampler.add_par('refPressure', value=1.0e+5, vary=False)

    # Rock matrix parameters
    frac_sampler.add_par('aveRockPerm', value=1.0e-22, vary=False)
    frac_sampler.add_par('stdDevRockPerm', value=1.0e-21, vary=False)
    frac_sampler.add_par('minRockPerm', value=5.0e-23, vary=False)
    frac_sampler.add_par('maxRockPerm', value=1.0e-21, vary=False)
    frac_sampler.add_par('refRockPerm', value=1.0e-22, vary=False)
    frac_sampler.add_par('refRockPressure', value=1.0e+6)
    frac_sampler.add_grid_obs('permeability', constr_type='array', index=[0],
                              output_dir=output_dir)
    frac_sampler.add_grid_obs('entryPressure', constr_type='array', index=[0],
                              output_dir=output_dir)

    sm.forward()

    # Collect saved observation
    perm_data = sm.collect_gridded_observations_as_time_series(
        frac_sampler, 'permeability', output_dir,
        indices=[0], rlzn_number=0)

    print(perm_data[0].shape)
    # print(perm_data[0])

    pressure_data = sm.collect_gridded_observations_as_time_series(
        frac_sampler, 'entryPressure', output_dir,
        indices=[0], rlzn_number=0)

    print(pressure_data[0].shape)
    # print(pressure_data[0])

    from matplotlib import cm
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d

    zz = np.zeros(nx*ny,)
    dxx = dx*np.ones(nx*ny,)
    dyy = dy*np.ones(nx*ny,)
    dzz = perm_data[0].reshape(nx*ny,)
    cmap = cm.get_cmap('Spectral')
    max_height = np.max(dzz)
    min_height = np.min(dzz)
    colors = [cmap((k-min_height)/max_height) for k in dzz]

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.bar3d(xx.reshape(nx*ny,), yy.reshape(nx*ny,), zz,
              dxx, dyy, dzz, shade=True, color=colors)
    ax2.set_xlabel('x [m]')
    ax2.set_ylabel('y [m]')
    ax2.set_zlabel(r'k [$m^2$]')

if __name__ == "__main__":

    test = 1
    if test == 1:
        test_scenario1()
    elif test == 2:
        test_scenario2()
