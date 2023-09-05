"""
This module contains several classes allowing to setup permeability of the wellbore
components randomly.

Created: July, 2018
Last modified: June, 2023

@author: Veronika Vasylkivska, Greg Lackey
"""
import sys
import os
import logging
from math import ceil
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from openiam import SamplerModel
except ImportError as err:
    print('Unable to load IAM class module: '+str(err))


PERM_GROUPS_DISTRIBUTIONS = {'small': [-20.0, -18.0],
                            'medium': [-18.0, -10.0],
                            'large': [-10.0, -9.0]}

def generate_seq(num_elems, seed=None, array_to_be_returned=2):
    """
    Generate two arrays: array representing permutation of indices and/or
    array of uniformly distributed random values in [0, 1]

    Parameters
    ----------
    num_elems : int
        Size of array(s) to be produced.
    seed : int, optional
        The default value is None. The value should be set if the output is
        to be reproducible.
    array_to_be_returned : int, optional. Can be 0, 1, or 2.
        The default value is 2: both arrays are returned. 0: only permutation
        array is returned. 1: only array with uniformly distributed values are
        returned.

    Returns
    -------
    numpy array or list of two numpy arrays.

    """
    if seed is None:
        prng = np.random.default_rng()  # pseudo-random number generator
    else:
        prng = np.random.default_rng(ceil(seed))  # seed has to be an integer

    # Create a random permutation of components indices
    indices_permuted = prng.permutation(range(1, num_elems+1))

    # Generate uniform random variables for permeability transformation
    rand_val = prng.random(num_elems)

    if array_to_be_returned == 0:
        return indices_permuted

    if array_to_be_returned == 1:
        return rand_val

    # If none of the checked conditions hold return both
    return indices_permuted, rand_val


class ParameterSetup1(SamplerModel):
    """
    ParameterSetup1 component generates permeability values for wellbore
    components of interest. Risk scores, if assigned to wellbores, can be used
    to generate permeability values. Two possible configurations for wellbore
    permeability distributions can be used.

    Configuration 1 generates at most two groups of wells: leaking
    ('exist_path') and non-leaking ('fixed'). If scores list
    is provided then the permeability distribution is defined by the component scores.

    Configuration 2 generates at most three groups of wells: leaking ('exist_path'),
    non-leaking ('fixed') and with paths to be possibly developed ('TBD').
    If scores are not provided, then configuration 2 can generate scores
    of 1 or 0.5 uniformly distributed within all three groups.
    """
    def __init__(self, name, parent, cmpnt_nms=None, scores=None,
                 par_name='logWellPerm', configuration=1, reproducible=True,
                 to_produce_scores=False):
        """
        Constructor method of ParameterSetup1 class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param cmpnt_nms: names of components whose parameters/keyword
            arguments/attributes will be controlled by the configurer
        :type cmpnt_nms: list()

        :param scores: scores based on which components arguments/attributes
            will be determined.
        :type scores: list()

        :param par_name: name of the parameter of the components that
            will be configured by the configurer
        :type par_name: str

        :param configuration: number defining parameter distribution
            for all components. Possible values are 1 and 2. Value of 1 means
            that it is possible to have at most two groups of wells: leaking
            ('exist_path') and non-leaking ('fixed'). If scores list
            is not None, then the distribution is defined by the component scores.
            If the value of configuration is 2 then it is possible to have
            three groups of wells: leaking ('exist_path'), non-leaking ('fixed')
            and with paths to be possibly developed ('TBD').
            If scores are not provided, then configuration 2 can generate scores
            of 1 or 0.5 approximately equally within all three groups.
        :type configuration: int

        :param reproducible: flag variable indicating whether the produced setup
            should be reproducible, i.e., if seed parameter of the model method
            should be used to generate the sequence of random numbers. By default,
            the value is True.
        :type reproducible: boolean

        :param to_produce_scores: flag indicating whether the scores
            should be produced for the components; it works only with configuration 2.

        :returns: ParameterSetup1 class object
        """
        # SamplerModel components have default parameter 'seed' assigned value of 1,
        # attribute 'run_frequency' set to 1 (component run only once), and
        # attribute 'reproducible' set to True (by default)
        super().__init__(name=name, parent=parent, reproducible=reproducible)

        # Add type attribute
        self.class_type = 'ParameterSetup1'

        # Setup the list of names of components which can be configured during the run
        if cmpnt_nms is None:
            self.cmpnts_to_configure = []
        else:
            self.cmpnts_to_configure = cmpnt_nms

        self.scores = scores
        self.parameter_name = par_name

        # Save configuration
        self.configuration = configuration
        self.to_produce_scores = to_produce_scores

        # Set default parameters of the component model
        self.add_default_par('seed', value=1)
        self.add_default_par('num_intervals', value=3)
        # Define base percentage of wells
        # 3 is maximum number of groups for score based perm distribution
        self.add_default_par('percentage1', value=0.95)  # small leak probability
        self.add_default_par('percentage2', value=0.04)  # medium leak probability
        self.add_default_par('percentage3', value=0.01)  # large leak probability

        # Define percentage adjustment factor for the second and third group
        # The default value of the factor is defined in a way so that
        # the percentages for the components with score 0.5 are not modified
        # Whether the percentages are modified in the model method depends on
        # whether the scores are provided by the contructor method.
        self.add_default_par('adj_factor', value=2.0)

        # Define percentage of wells with paths to be developed in each group
        self.add_default_par('dev_path_percentage1', value=0.0)
        self.add_default_par('dev_path_percentage2', value=0.5)
        self.add_default_par('dev_path_percentage3', value=0.5)

        # Define the value that defines the threshold, i.e. components with
        # the parameter below that threshold are not run.
        self.add_default_par('threshold', value=-np.inf)

        # Define min and max of the uniform distribution for each percentage group
        self.add_default_par('low_bound1', value=-20.0)
        self.add_default_par('high_bound1', value=-18.0)
        self.add_default_par('low_bound2', value=-18.0)
        self.add_default_par('high_bound2', value=-10.0)
        self.add_default_par('low_bound3', value=-10.0)
        self.add_default_par('high_bound3', value=-9.0)

        msg = 'ParameterSetup1 created with name {name}'.format(name=self.name)
        logging.debug(msg)

    @property
    def num_cmpnts(self):
        return len(self.cmpnts_to_configure)

    def sample(self, p):
        """
        Return initial value of the parameters for the configured components.

        :param p: dictionary of input parameters
        :type p: dict

        :returns: dictionary of observations
        """
        if len(self.cmpnts_to_configure) == 0:
            logging.error('No components were linked for configuration.')
            raise ValueError(
                ''.join(['The model method of component {} ',
                         'cannot create an initial distribution ',
                         'configuration: the list of components ',
                         'to be configured is empty']).format(self.name))

        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Check that the sum of percentages is 1, otherwise normalize
        S = sum([actual_p['percentage{}'.format(ind+1)] for ind in range(
            actual_p['num_intervals'])])
        if S != 1:
            for ind in range(actual_p['num_intervals']):
                # Divide each percentage by S
                actual_p['percentage{}'.format(ind+1)] = (
                    actual_p['percentage{}'.format(ind+1)]/S)

        # Determine the boundaries of the intervals
        boundaries = np.zeros((len(self.cmpnts_to_configure), actual_p['num_intervals']+1))
        boundaries[:, -1] = 1

        # Check whether the scores were provided
        if self.scores is None:
            for ind in range(1, actual_p['num_intervals']):
                boundaries[:, ind] = boundaries[:, ind-1] + actual_p['percentage{}'.format(ind)]
        else:
            # Check whether the length of the scores list and number of components are the same
            if self.num_cmpnts != len(self.scores):
                err_msg = ''.join([
                    'Number of scores provided does not coincide ',
                    'with the number of linked components.'])
                logging.error(err_msg)
                raise ValueError(err_msg)

            boundaries[:, 2] = np.maximum(
                np.zeros(self.num_cmpnts),
                1 - actual_p['percentage3']*actual_p['adj_factor']*self.scores)
            boundaries[:, 1] = np.maximum(
                np.zeros(self.num_cmpnts),
                boundaries[:, 2] - actual_p['percentage2']*actual_p['adj_factor']*self.scores)

        # Initialize the output dictionary of the model
        out = {}
        out[self.parameter_name] = np.zeros(self.num_cmpnts)
        out['run'] = 2*np.ones(self.num_cmpnts)

        # The consecutive random number generation will be the same
        # for the runs with the same seed and on the same machine
        if self.reproducible:    # if reproducible results are needed
            indices_permuted, rand_val = generate_seq(
                self.num_cmpnts, seed=actual_p['seed'])
        else:
            indices_permuted, rand_val = generate_seq(self.num_cmpnts)

        # Setup a dictionary to keep information about
        # in which group each of the components ends up
        self.details['ind_groups'] = {ind+1: [] for ind in range(actual_p['num_intervals'])}

        # Cycle over all components to set the needed parameter
        for cind in range(self.num_cmpnts):
            for ind in range(1, actual_p['num_intervals']+1):
                ind_ratio = indices_permuted[cind]/self.num_cmpnts
                if (boundaries[cind, ind-1] < ind_ratio <= boundaries[cind, ind]):
                    # Calculate parameter value
                    out[self.parameter_name][cind] = (actual_p['low_bound{}'.format(ind)]+(
                        actual_p['high_bound{}'.format(ind)]
                        -actual_p['low_bound{}'.format(ind)])*rand_val[cind])
                    # Save component index in the corresponding group
                    self.details['ind_groups'][ind].append(cind)

        # Cycle over all components to check which of the components
        # have parameters assigned below the threshold
        for cind in range(self.num_cmpnts):
            # Get the component
            nm = self.cmpnts_to_configure[cind]
            cmpnt = self._parent.component_models[nm]
            # If generated value of the parameter is below the threshold
            # the component should not be run
            if out[self.parameter_name][cind] <= actual_p['threshold']:
                cmpnt.run_frequency = 0
                out['run'][cind] = 0
                cmpnt.details['risk_type'] = 'fixed'
            else:
                cmpnt.run_frequency = 2
                cmpnt.details['risk_type'] = 'exist_path'

            if 'normalized_leak_score' in cmpnt.details:
                if cmpnt.details['normalized_leak_score'] is None:
                    cmpnt.details['normalized_leak_score'] = 0.0

        # Determine length of each group
        group_length = {
            ind: len(self.details['ind_groups'][ind]) for ind in range(
                1, actual_p['num_intervals']+1)}

        # Produce scores if requested
        if self.to_produce_scores:
            pos_scores_options = {1: 1, 0: 0.5}
            scores = {ind: [pos_scores_options[i%2] for i in range(
                group_length[ind])] for ind in range(1, actual_p['num_intervals']+1)}

        # Check whether the second configuration is required:
        # Ni% of wells in each of the i interval
        # with existing leaking paths are replaced with paths
        # to be possibly developed in the future
        if self.configuration == 2:
            # Cycle over all groups of wells
            for ind in range(1, actual_p['num_intervals']+1):
                # Cycle over indices of components in the current group
                for j, cmpnt_ind in enumerate(self.details['ind_groups'][ind]):
                    # Some of components in the current group have to change
                    # their properties if their relative index satisfies the condition
                    nm = self.cmpnts_to_configure[cmpnt_ind]
                    cmpnt = self._parent.component_models[nm]
                    if (j+1)/group_length[ind] < actual_p['dev_path_percentage{}'.format(ind)]:
                        cmpnt.run_frequency = 0
                        cmpnt.details['risk_type'] = 'TBD'
                    if self.to_produce_scores:
                        cmpnt.details['normalized_leak_score'] = scores[ind][j]

#            # Turn off old wells
#            for cind in range(self.num_cmpnts):
#                # Get the component
#                nm = self.cmpnts_to_configure[cind]
#                cmpnt = self._parent.component_models[nm]
#                if cmpnt.details['normalized_leak_score'] == 1:
#                    cmpnt.details['normalized_leak_score'] = 0.0
#                    cmpnt.details['risk_type'] = 'fixed'
#                    cmpnt.run_frequency = 0
#                    out['run'][cind] = 0

        return out


class ParameterSetup2(SamplerModel):
    """
    ParameterSetup2 component generates permeability values for wellbore
    components of interest. Different permeability distribution can be specified
    for three groups of wells: (1) not leaking, (2) leaking, and (3) potentially
    leaking. Permeability is classified as being 'small', 'medium', or 'large'.
    The boundaries are defined and can be changed by the component setup.
    Each of the three wellbore groups can be assigned a different percent
    for each group of permeability.
    """
    def __init__(self, name, parent, cmpnt_nms=None, par_name='logWellPerm',
                 reproducible=True, well_status_perc=None, perm_bounds=None,
                 perm_perc=None):
        """
        Constructor method of ParameterSetup2 class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param cmpnt_nms: names of components whose parameters/keyword
            arguments/attributes will be controlled by the configurer
        :type cmpnt_nms: list()

        :param par_name: name of the parameter of the components that
            will be configured by the configurer
        :type par_name: str

        :param reproducible: flag variable indicating whether the produced setup
            should be reproducible, i.e. if seed parameter of the model method
            should be used to generate the sequence of random numbers. By default,
            the value is True.
        :type reproducible: boolean

        :param well_status_perc: percentages of the wells belonging to
            not leaking (1), leaking (2) and potentially leaking (3) groups.
            By default, {1: 0.0, 2: 1.0, 3: 0.0}.
        :type well_status_perc: dict()

        :param perm_bounds: dictionary with keys ['small', 'medium', 'large']
            defining boundaries for each permeability group. For example,
            dictionary
                {'small': [-20.0, -18.0],
                 'medium': [-18.0, -10.0],
                 'large': [-10.0, -9.0]}}
            defines 'small' permeability as being between 10^(-20) and 10^(-18).
            By default, {'small': [-20.0, -18.0], 'medium': [-18.0, -10.0],
                         'large': [-10.0, -9.0]}
        :type perm_bounds: dict()

        :param perm_perc: dictionary with keys [1, 2, 3] defining percentages of
            wellbores with different permeability within each group of wellbores.
            For example, dictionary
                {1: [1.0, 0.0, 0.0],
                 2: [0.75, 0.2, 0.05],
                 3: [0.75, 0.2, 0.05]}
            defines that (not leaking) wells within group 1 can have only 'small'
            permeability (value of 1 corresponds to 100%); (leaking) wells within
            group2 have higher chance (75%) to have 'small' permeability, 20% have
            'medium' permeability and 5% have 'large' permeability.
            By default, {1: [1.0, 0.0, 0.0],  2: [0.75, 0.2, 0.05],
                         3: [0.75, 0.2, 0.05]}
        :type param perm_perc: dict()

        :returns: ParameterSetup2 class object
        """
        # SamplerModel components have default parameter 'seed' assigned value of 1,
        # attribute 'run_frequency' set to 1 (component run only once), and
        # attribute 'reproducible' set to True (by default)
        super().__init__(name=name, parent=parent, reproducible=reproducible)

        # Add type attribute
        self.class_type = 'ParameterSetup2'

        # Setup the list of names of components which can be configured during the run
        if cmpnt_nms is None:
            self.cmpnts_to_configure = []
        else:
            self.cmpnts_to_configure = cmpnt_nms

        # Save parameter name
        self.parameter_name = par_name

        # Save distribution properties
        # well_status_perc: 1 - not leaking, 2 - leaking, 3 - potentially leaking
        if well_status_perc is None:
            self.well_status_perc = {1: 0.0, 2: 1.0, 3: 0.0} # 100% of wells are leaking
        else:
            self.well_status_perc = well_status_perc

        if perm_bounds is None:
            self.perm_bounds = PERM_GROUPS_DISTRIBUTIONS
        else:
            self.perm_bounds = perm_bounds

        if perm_perc is None:
            self.perm_perc = {1: [1.0, 0.0, 0.0],    # 1 - perm dist for not leaking wells
                              2: [0.75, 0.2, 0.05],  # 2 - for leaking wells
                              3: [0.75, 0.2, 0.05]}  # 3 - for potentially leaking wells
        else:
            self.perm_perc = perm_perc

        msg = 'ParameterSetup2 created with name {name}'.format(name=self.name)
        logging.debug(msg)

    def sample(self, p):
        """
        Return initial value of the parameters for the configured components.

        :param p: dictionary of input parameters
        :type p: dict

        :returns: dictionary of observations
        """
        if len(self.cmpnts_to_configure) == 0:
            logging.error('No components were linked for configuration.')
            raise ValueError(' '.join([
                'The model method of component {}',
                'cannot create an initial distribution configuration: the list of components',
                'to be configured is empty']).format(self.name))

        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Calculate number of components
        num_cmpnts = len(self.cmpnts_to_configure)

        # The consecutive random number generation will be the same
        # for the runs with the same seed and on the same machine
        if self.reproducible:    # if reproducible results are needed
            indices_permuted1 = generate_seq(
                num_cmpnts, seed=actual_p['seed'], array_to_be_returned=0)
        else:
            indices_permuted1 = generate_seq(
                num_cmpnts, array_to_be_returned=0)

        # Initialize the output dictionary of the model
        out = {}
        out[self.parameter_name] = np.zeros(num_cmpnts)
        out['run'] = 2*np.ones(num_cmpnts)
        self.well_status_groups = {k: [] for k in range(1, 4)}

        # Choose a status for all wells using provided desired distribution
        # Cycle over all components
        for cind in range(num_cmpnts):
            nm = self.cmpnts_to_configure[cind]
            cmpnt = self._parent.component_models[nm]
            ind_ratio = indices_permuted1[cind]/num_cmpnts
            if ind_ratio <= self.well_status_perc[1]:
                cmpnt.run_frequency = 0
                out['run'][cind] = 0
                cmpnt.details['risk_type'] = 'fixed'
                self.well_status_groups[1].append(cind)
            elif ind_ratio <= self.well_status_perc[1]+self.well_status_perc[2]:
                cmpnt.run_frequency = 2
                out['run'][cind] = 2
                cmpnt.details['risk_type'] = 'exist_path'
                self.well_status_groups[2].append(cind)
            else:
                cmpnt.run_frequency = 0
                out['run'][cind] = 0
                cmpnt.details['risk_type'] = 'TBD'
                self.well_status_groups[3].append(cind)

            # Determine length of each group
            groups_len = {k: len(self.well_status_groups[k]) for k in range(1, 4)}

            if self.reproducible:    # if reproducible results are needed
                seeds = [actual_p['seed']+g_ind for g_ind in range(1, 4)]
            else:
                seeds = 3*[None]

            self.details['well_status_groups'] = self.well_status_groups
            self.details['perm_groups'] = {'small': [], 'medium': [], 'large': []}

            # Assign permeability to the wells in each group
            for g_ind in range(1, 4):
                indices_permuted2, rand_val = generate_seq(
                    groups_len[g_ind],
                    seed=seeds[g_ind-1],
                    array_to_be_returned=2)

                for ind, sind in enumerate(self.well_status_groups[g_ind]):
                    nm = self.cmpnts_to_configure[sind]
                    cmpnt = self._parent.component_models[nm]
                    ind_ratio = indices_permuted2[ind]/groups_len[g_ind]

                    if ind_ratio <= self.perm_perc[g_ind][0]:
                        self.details['perm_groups']['small'].append(sind)
                        out[self.parameter_name][sind] = self.perm_bounds['small'][0]+(
                            self.perm_bounds['small'][1]
                            -self.perm_bounds['small'][0])*rand_val[ind]
                    elif ind_ratio <= (self.perm_perc[g_ind][0]+self.perm_perc[g_ind][1]):
                        self.details['perm_groups']['medium'].append(sind)
                        out[self.parameter_name][sind] = self.perm_bounds['medium'][0]+(
                            self.perm_bounds['medium'][1]
                            -self.perm_bounds['medium'][0])*rand_val[ind]
                    else:
                        self.details['perm_groups']['large'].append(sind)
                        out[self.parameter_name][sind] = self.perm_bounds['large'][0]+(
                            self.perm_bounds['large'][1]
                            -self.perm_bounds['large'][0])*rand_val[ind]

        return out


class ParameterSetup3(SamplerModel):
    """
    ParameterSetup3 component generates permeability values for wellbore
    components of interest. Different permeability distribution can be specified
    for multiple groups of wells separated by their age. Permeability
    is classified as being 'small', 'medium', or 'large'.
    The boundaries for 'small', 'mediuam', and 'large' are defined
    and can be changed by the component setup.
    Each of the wellbore groups (based on age) can be assigned a different percent
    for each group of permeability.
    """
    def __init__(self, name, parent, cmpnt_nms=None, ages=None,
                 par_name='logWellPerm', reproducible=True,
                 well_status_perc=None, well_status_perc_adj_factor=None,
                 perm_bounds=None, perm_perc=None):
        """
        Constructor method of ParameterSetup3 class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param cmpnt_nms: names of components whose parameters/keyword
            arguments/attributes will be controlled by the configurer
        :type cmpnt_nms: list()

        :param ages: ages based on which components arguments/attributes
            will be determined.
        :type scores: list()

        :param par_name: name of the parameter of the components that
            will be configured by the configurer
        :type par_name: str

        :param reproducible: flag variable indicating whether the produced setup
            should be reproducible, i.e. if seed parameter of the model method
            should be used to generate the sequence of random numbers. By default,
            the value is True.
        :type reproducible: boolean

        :param well_status_perc: percentages of the wells belonging to
            not leaking (1), leaking (2) and potentially leaking (3) groups.
            By default, {1: 0.0, 2: 1.0, 3: 0.0}.
        :type well_status_perc: dict()

        :param well_status_perc_adj_factor: parameter adjusting percentage of
            the three well groups described in the description of well_status_perc
            parameter. By default, {1: 1.5, 2: 1.0, 3: 0.5, 4: 0.1}.
        :type well_status_perc_adj_factor: dict()

        :param perm_bounds: dictionary with keys ['small', 'medium', 'large']
            defining boundaries for each permeability group. For example,
            dictionary
                {'small': [-20.0, -18.0],
                 'medium': [-18.0, -10.0],
                 'large': [-10.0, -9.0]}}
            defines 'small' permeability as being between 10^(-20) and 10^(-18).
            By default, {'small': [-20.0, -18.0], 'medium': [-18.0, -10.0],
                         'large': [-10.0, -9.0]}
        :type perm_bounds: dict()

        :param perm_perc: dictionary with keys [1, 2, 3, 4] defining percentages of
            wellbores with different permeability within each age group of wellbores.
            For example, dictionary
                {1: [0.0, 0.8, 0.2],
                 2: [0.0, 0.8, 0.2],
                 3: [0.0, 0.8, 0.2],
                 4: [0.0, 0.8, 0.2]}
            defines that all age groups have the same permeability
            distribution: 80% have 'medium' permeability and 20% have 'large'
            permeability. By default, {1: [0.0, 0.8, 0.2], 2: [0.0, 0.8, 0.2],
                                       3: [0.0, 0.8, 0.2], 4: [0.0, 0.8, 0.2]}
        :type param perm_perc: dict()


        :returns: ParameterSetup3 class object
        """
        # SamplerModel components have default parameter 'seed' assigned value of 1,
        # attribute 'run_frequency' set to 1 (component run only once), and
        # attribute 'reproducible' set to True (by default)
        super().__init__(name=name, parent=parent, reproducible=reproducible)

        # Add type attribute
        self.class_type = 'ParameterSetup3'

        if cmpnt_nms is None:
            self.cmpnts_to_configure = []
        else:
            self.cmpnts_to_configure = cmpnt_nms

        # Save distribution properties
        # well_status_perc: 1 - not leaking, 2 - leaking, 3 - potentially leaking
        if well_status_perc is None:
            self.well_status_perc = {1: 0.75, 2: 0.125, 3: 0.125}
        else:
            self.well_status_perc = well_status_perc

        if well_status_perc_adj_factor is None:
            well_status_perc_adj_factor = {1: 1.5, 2: 1.0, 3: 0.5, 4: 0.1},

        if perm_bounds is None:
            self.perm_bounds = PERM_GROUPS_DISTRIBUTIONS
        else:
            self.perm_bounds = perm_bounds

        if perm_perc is None:
            self.perm_perc = {1: [0.0, 0.8, 0.2],   # 1 - perm dist for wells of age group 1
                              2: [0.0, 0.8, 0.2],   # 2 - perm dist for wells of age group 2
                              3: [0.0, 0.8, 0.2],   # 3 - perm dist for wells of age group 3
                              4: [0.0, 0.8, 0.2]}   # 4 - perm dist for wells of age group 4
        else:
            self.perm_perc = perm_perc

        # Save wellbore ages
        self.ages = ages

        # Determine which of the wells in which group
        self.age_groups = {}
        self.well_status_perc_adj = {}

        for age in np.unique(self.ages):
            self.age_groups[age] = np.where(self.ages == age)[0]
            # Adjust percentage of the wells of each age in each event category
            self.well_status_perc_adj[age] = {}
            if age in well_status_perc_adj_factor:
                scalar = well_status_perc_adj_factor[age]

                self.well_status_perc_adj[age][3] = min(
                    1.0, self.well_status_perc[3]*scalar)

                self.well_status_perc_adj[age][2] = min(
                    1.0-self.well_status_perc_adj[age][3],
                    self.well_status_perc[2]*scalar)

                self.well_status_perc_adj[age][1] = max(
                    0.0,
                    1.0-self.well_status_perc_adj[age][2]-self.well_status_perc_adj[age][3])
            else:
                self.well_status_perc_adj[age] = self.well_status_perc

        # Save name of parameter that is getting changed
        self.parameter_name = par_name

        msg = 'ParameterSetup3 created with name {name}'.format(name=self.name)
        logging.debug(msg)

    def sample(self, p):
        """
        Return initial value of the parameters for the configured components.

        :param p: dictionary of input parameters
        :type p: dict

        :returns: dictionary of observations
        """
        if len(self.cmpnts_to_configure) == 0:
            logging.error('No components were linked for configuration.')
            raise ValueError(' '.join([
                'The model method of component {}',
                'cannot create an initial distribution configuration: the list of components',
                'to be configured is empty']).format(self.name))

        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Calculate number of components
        num_cmpnts = len(self.cmpnts_to_configure)

        # The consecutive random number generation will be the same
        # for the runs with the same seed and on the same machine
        if self.reproducible:    # if reproducible results are needed
            indices_permuted1 = generate_seq(
                num_cmpnts, seed=actual_p['seed'], array_to_be_returned=0)
        else:
            indices_permuted1 = generate_seq(
                num_cmpnts, array_to_be_returned=0)

        # Initialize the output dictionary of the model
        out = {}
        out[self.parameter_name] = np.zeros(num_cmpnts)
        out['run'] = 2*np.ones(num_cmpnts)
        self.well_status_groups = {k: [] for k in range(1, 4)}

        # Choose a status for all wells using provided desired distribution
        # Cycle over all components
        for cind in range(num_cmpnts):
            nm = self.cmpnts_to_configure[cind]
            cmpnt = self._parent.component_models[nm]
            ind_ratio = indices_permuted1[cind]/num_cmpnts
            c_age = self.ages[cind]
            if ind_ratio <= self.well_status_perc_adj[c_age][1]:
                cmpnt.run_frequency = 0
                out['run'][cind] = 0
                cmpnt.details['risk_type'] = 'fixed'
                self.well_status_groups[1].append(cind)
            elif ind_ratio <= (
                    self.well_status_perc_adj[c_age][1]
                    +self.well_status_perc_adj[c_age][2]):
                cmpnt.run_frequency = 2
                out['run'][cind] = 2
                cmpnt.details['risk_type'] = 'exist_path'
                self.well_status_groups[2].append(cind)
            else:
                cmpnt.run_frequency = 0
                out['run'][cind] = 0
                cmpnt.details['risk_type'] = 'TBD'
                self.well_status_groups[3].append(cind)

            if self.reproducible:    # if reproducible results are needed
                seeds = [actual_p['seed']+ag_ind+1 for ag_ind in range(len(self.age_groups))]
            else:
                seeds = [None for ag_ind in range(len(self.age_groups))]

            self.details['well_status_groups'] = self.well_status_groups
            self.details['age_groups'] = self.age_groups
            self.details['perm_groups'] = {'small': [], 'medium': [], 'large': []}

            # Assign permeability to the wells in each group
            for g_ind, age_gr in enumerate(self.age_groups.keys()):
                gr_indices = self.age_groups[age_gr]
                indices_permuted2, rand_val = generate_seq(
                    len(gr_indices), seed=seeds[g_ind-1], array_to_be_returned=2)

                for ind, sind in enumerate(gr_indices):
                    nm = self.cmpnts_to_configure[sind]
                    cmpnt = self._parent.component_models[nm]
                    ind_ratio = indices_permuted2[ind]/len(gr_indices)

                    if ind_ratio <= self.perm_perc[age_gr][0]:
                        self.details['perm_groups']['small'].append(sind)
                        out[self.parameter_name][sind] = self.perm_bounds['small'][0]+(
                            self.perm_bounds['small'][1]
                            -self.perm_bounds['small'][0])*rand_val[ind]

                    elif ind_ratio <= (self.perm_perc[age_gr][0]+self.perm_perc[age_gr][1]):
                        self.details['perm_groups']['medium'].append(sind)
                        out[self.parameter_name][sind] = self.perm_bounds['medium'][0]+(
                            self.perm_bounds['medium'][1]
                            -self.perm_bounds['medium'][0])*rand_val[ind]

                    else:
                        self.details['perm_groups']['large'].append(sind)
                        out[self.parameter_name][sind] = self.perm_bounds['large'][0]+(
                            self.perm_bounds['large'][1]
                            -self.perm_bounds['large'][0])*rand_val[ind]

        return out


class ParameterSetup4(SamplerModel):
    """
    ParameterSetup4 component generates permeability values for wellbore
    components of interest. Different permeability distribution can be specified
    for three groups of wells: (1) not leaking, (2) leaking, and (3) potentially
    leaking. Permeability is classified as being 'small', 'medium', or 'large'.
    The boundaries for 'small', 'mediuam', and 'large' are defined
    and can be changed by the component setup. Percentage of the wells assigned
    to different groups can be controlled further by setup of the
    well_status_perc_adj_factor argument.
    """
    def __init__(self, name, parent, cmpnt_nms=None, scores=None,
                 par_name='logWellPerm', reproducible=True,
                 well_status_perc=None, well_status_perc_adj_factor=None,
                 perm_bounds=None, perm_perc=None):
        """
        Constructor method of ParameterSetup4 class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param cmpnt_nms: names of components whose parameters/keyword
            arguments/attributes will be controlled by the configurer
        :type cmpnt_nms: list()

        :param scores: scores are used to determine components arguments/attributes
        :type scores: list()

        :param par_name: name of the parameter of the components that
            will be configured by the configurer
        :type par_name: str

        :param reproducible: flag variable indicating whether the produced setup
            should be reproducible, i.e. if seed parameter of the model method
            should be used to generate the sequence of random numbers. By default,
            the value is True.
        :type reproducible: boolean

        :param well_status_perc: base percentages of the wells belonging to
            not leaking (1), leaking (2) and potentially leaking (3) groups.
            By default, {1: 0.0, 2: 1.0, 3: 0.0}.
        :type well_status_perc: dict()

        :param well_status_perc_adj_factor: the values of scores and of corresponding
            adjustment factors which would modify the base percentages. Argument
            is a matrix of size (N, 2) where N is a number of known scores
            and associated adjustment factors. The first column represents scores,
            the second column represents corresponding adjustment factors.
            By default, np.array(
                [[1, 1.5], [0.8, 1.], [0.5, 0.5], [0.2, 0.1], [0.0, 0.0]])
        :type ws_perc_adj_factor: np.array

        :param perm_bounds: dictionary with keys ['small', 'medium', 'large']
            defining boundaries for each permeability group. For example,
            dictionary
                {'small': [-20.0, -18.0],
                 'medium': [-18.0, -10.0],
                 'large': [-10.0, -9.0]}}
            defines 'small' permeability as being between 10^(-20) and 10^(-18).
            By default, {'small': [-20.0, -18.0], 'medium': [-18.0, -10.0],
                         'large': [-10.0, -9.0]}
        :type perm_bounds: dict()

        :param perm_perc: dictionary which contain the information
            of the permeability distribution for each group of wells based
            on the type of leakage. By default, {1: [1.0, 0.0, 0.0],
                                                 2: [0.75, 0.2, 0.05],
                                                 3: [0.75, 0.2, 0.05]}
        :type perm_perc: dict()

        :returns: ParameterSetup4 class object
        """
        # SamplerModel components have default parameter 'seed' assigned value of 1,
        # attribute 'run_frequency' set to 1 (component run only once), and
        # attribute 'reproducible' set to True (by default)
        super().__init__(name=name, parent=parent, reproducible=reproducible)

        # Add type attribute
        self.class_type = 'ParameterSetup4'

        # Setup the list of names of components which can be configured during the run
        if cmpnt_nms is None:
            self.cmpnts_to_configure = []
        else:
            self.cmpnts_to_configure = cmpnt_nms

        # Save distribution properties
        # well_status_perc: 1 - not leaking, 2 - leaking, 3 - potentially leaking
        if well_status_perc is None:
            self.well_status_perc = {1: 0.75, 2: 0.125, 3: 0.125}
        else:
            self.well_status_perc = well_status_perc

        if well_status_perc_adj_factor is None:
            well_status_perc_adj_factor = np.array(
                [[1, 1.5], [0.8, 1.], [0.5, 0.5], [0.2, 0.1], [0.0, 0.0]])

        if perm_bounds is None:
            self.perm_bounds = PERM_GROUPS_DISTRIBUTIONS
        else:
            self.perm_bounds = perm_bounds

        if perm_perc is None:
            self.perm_perc = {1: [1.0, 0.0, 0.0],    # 1 - perm dist for not leaking wells
                              2: [0.75, 0.2, 0.05],  # 2 - for leaking wells
                              3: [0.75, 0.2, 0.05]}  # 3 - for potentially leaking wells
        else:
            self.perm_perc = perm_perc


        # Save wellbore ages
        self.scores = scores

        # Determine which of the wells in which group
        self.well_status_perc_adj = {}
        sort_ind = np.argsort(well_status_perc_adj_factor[:, 0])

        self.adj_factors = np.interp(
            scores,
            well_status_perc_adj_factor[:, 0][sort_ind],
            well_status_perc_adj_factor[:, 1][sort_ind])

        self.well_status_perc_adj = [
            {1: self.well_status_perc[1],
             2: self.well_status_perc[2],
             3: self.well_status_perc[3]} for score in self.scores]

        for ind in range(len(self.scores)):
            self.well_status_perc_adj[ind][3] = min(
                1.0, self.well_status_perc[3]*self.adj_factors[ind])
            self.well_status_perc_adj[ind][2] = min(
                1.0-self.well_status_perc_adj[ind][3],
                self.well_status_perc[2]*self.adj_factors[ind])
            self.well_status_perc_adj[ind][1] = max(
                0.0,
                1.0-self.well_status_perc_adj[ind][2]-self.well_status_perc_adj[ind][3])

        # Save name of parameter that is getting changed
        self.parameter_name = par_name

        msg = 'ParameterSetup4 created with name {name}'.format(name=self.name)
        logging.debug(msg)

    def sample(self, p):
        """
        Return initial value of the parameters for the configured components.

        :param p: dictionary of input parameters
        :type p: dict

        :returns: dictionary of observations
        """
        if len(self.cmpnts_to_configure) == 0:
            logging.error('No components were linked for configuration.')
            raise ValueError(' '.join([
                'The model method of component {}',
                'cannot create an initial distribution configuration: the list of components',
                'to be configured is empty']).format(self.name))

        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Calculate number of components
        num_cmpnts = len(self.cmpnts_to_configure)

        # The consecutive random number generation will be the same
        # for the runs with the same seed and on the same machine
        if self.reproducible:    # if reproducible results are needed
            indices_permuted1 = generate_seq(
                num_cmpnts, seed=actual_p['seed'], array_to_be_returned=0)
        else:
            indices_permuted1 = generate_seq(
                num_cmpnts, array_to_be_returned=0)

        # Initialize the output dictionary of the model
        out = {}
        out[self.parameter_name] = np.zeros(num_cmpnts)
        out['run'] = 2*np.ones(num_cmpnts)
        self.well_status_groups = {k: [] for k in range(1, 4)}

        # Choose a status for all wells using provided desired distribution
        # Cycle over all components
        for cind in range(num_cmpnts):
            nm = self.cmpnts_to_configure[cind]
            cmpnt = self._parent.component_models[nm]
            ind_ratio = indices_permuted1[cind]/num_cmpnts

            if ind_ratio <= self.well_status_perc_adj[cind][1]:
                cmpnt.run_frequency = 0
                out['run'][cind] = 0
                cmpnt.details['risk_type'] = 'fixed'
                self.well_status_groups[1].append(cind)
            elif ind_ratio <= (
                    self.well_status_perc_adj[cind][1]
                    +self.well_status_perc_adj[cind][3]):
                cmpnt.run_frequency = 0
                out['run'][cind] = 1
                cmpnt.details['risk_type'] = 'TBD'
                self.well_status_groups[3].append(cind)
            else:
                cmpnt.run_frequency = 2
                out['run'][cind] = 2
                cmpnt.details['risk_type'] = 'exist_path'
                self.well_status_groups[2].append(cind)

            if self.reproducible:    # if reproducible results are needed
                seeds = [actual_p['seed']+ind for ind in [23, 87]]  # just arbitrary numbers
            else:
                seeds = 2*[None]

            self.details['well_status_groups'] = self.well_status_groups
            self.details['perm_groups'] = {'small': [],
                                           'medium': [],
                                           'large': []}

            # Assign permeability to the wells in each group
            rand_val1 = generate_seq(
                num_cmpnts, seed=seeds[0], array_to_be_returned=1)
            rand_val2 = generate_seq(
                num_cmpnts, seed=seeds[1], array_to_be_returned=1)

            for g_ind in range(1, 4):
                for sind in self.well_status_groups[g_ind]:
                    nm = self.cmpnts_to_configure[sind]
                    cmpnt = self._parent.component_models[nm]
                    ind_ratio = rand_val1[sind]

                    if ind_ratio <= self.perm_perc[g_ind][0]:
                        self.details['perm_groups']['small'].append(sind)
                        out[self.parameter_name][sind] = self.perm_bounds['small'][0]+(
                            self.perm_bounds['small'][1]
                            -self.perm_bounds['small'][0])*rand_val2[sind]

                    elif ind_ratio <= (self.perm_perc[g_ind][0]+self.perm_perc[g_ind][1]):
                        self.details['perm_groups']['medium'].append(sind)
                        out[self.parameter_name][sind] = self.perm_bounds['medium'][0]+(
                            self.perm_bounds['medium'][1]
                            -self.perm_bounds['medium'][0])*rand_val2[sind]

                    else:
                        self.details['perm_groups']['large'].append(sind)
                        out[self.parameter_name][sind] = self.perm_bounds['large'][0]+(
                            self.perm_bounds['large'][1]
                            -self.perm_bounds['large'][0])*rand_val2[sind]

        return out


class ParameterSetup5(SamplerModel):
    """
    ParameterSetup5 component generates permeability values for wellbore
    components of interest. The wells are separated into groups based on
    the leakage types ('small', 'medium', 'large' leaks). Each of the leakage
    type group is assigned different permeability ranges.
    """
    def __init__(self, name, parent, cmpnt_nms=None, scores=None,
                 par_name='logWellPerm', reproducible=True,
                 # leakage_type_perc: keys 1 - small leak, 2 - medium leak, 3 - large leak
                 leakage_type_perc=None, perm_bounds=None, configuration=1):
        """
        Constructor method of ParameterSetup5 class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param cmpnt_nms: names of components whose parameters/keyword
            arguments/attributes will be controlled by the configurer
        :type cmpnt_nms: list()

        :param scores: scores are used to determine components arguments/attributes
        :type scores: list()

        :param par_name: name of the parameter of the components that
            will be configured by the configurer
        :type par_name: str

        :param reproducible: flag variable indicating whether the produced setup
            should be reproducible, i.e. if seed parameter of the model method
            should be used to generate the sequence of random numbers. By default,
            the value is True. By default, the value is True.
        :type reproducible: boolean

        :param leakage_type_perc: base percentages of the wells having 'small'
            'medium', amd 'large' permeability. By default,
            {1: 0.954, 2: 0.044, 3: 0.002}
        :type leakage_type_perc: dict()

        :param perm_bounds: dictionary which contain the information
            of the permeability distribution for each group of wells based
            on the type of leakage. By default, {'small': [-19.0, -19.0],
                                                 'medium': [-17.0, -14.0],
                                                 'large': [-13.0, -12.0]}.
        :type perm_perc: dict()

        :param configuration: additional parameter defining distribution
            for all components. Possible values are 1 (default), 2, and 3.
            Value of 1 means that permeability is not correlated with scores.
            Value of 2 means that permeability is correlated with scores but only
            in terms of distribution: there is no relationship between score and
            assumed permeability.
            Value of 3 means that permeability is correlated with scores: wells
            with larger scores will be assigned larger permeability; permeability
            values are sorted after randomization and assigned to the wells
            in the order the scores rank the wells.
        :type configuration: int

        :returns: ParameterSetup5 class object
        """
        # SamplerModel components have default parameter 'seed' assigned value of 1,
        # attribute 'run_frequency' set to 1 (component run only once), and
        # attribute 'reproducible' set to True (by default)
        super().__init__(name=name, parent=parent, reproducible=reproducible)

        # Add type attribute
        self.class_type = 'ParameterSetup5'

        # Setup the list of names of components which can be configured during the run
        if cmpnt_nms is None:
            self.cmpnts_to_configure = []
        else:
            self.cmpnts_to_configure = cmpnt_nms

        # Save distribution properties
        if leakage_type_perc is None:
            # 1 - small leak, 2 - medium leak, 3 - large leak
            self.leakage_type_perc = {1: 0.954, 2: 0.044, 3: 0.002}
        else:
            self.leakage_type_perc = leakage_type_perc

        if perm_bounds is None:
            self.perm_bounds = {'small': [-19.0, -19.0],
                                'medium': [-17.0, -14.0],
                                'large': [-13.0, -12.0]}
        else:
            self.perm_bounds = perm_bounds

        # Save wellbore ages
        self.scores = scores

        self.configuration = configuration

        # Save name of parameter that is getting changed
        self.parameter_name = par_name

        msg = 'ParameterSetup5 created with name {name}'.format(name=self.name)
        logging.debug(msg)

    def sample(self, p):
        """
        Return initial value of the parameters for the configured components.

        :param p: dictionary of input parameters
        :type p: dict

        :returns: dictionary of observations
        """
        if len(self.cmpnts_to_configure) == 0:
            logging.error('No components were linked for configuration.')
            raise ValueError(' '.join([
                'The model method of component {}',
                'cannot create an initial distribution configuration: the list of components',
                'to be configured is empty']).format(self.name))

        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Calculate number of components
        num_cmpnts = len(self.cmpnts_to_configure)

        # The consecutive random number generation will be the same
        # for the runs with the same seed and on the same machine
        if self.reproducible:    # if reproducible results are needed
            indices_permuted1, rand_val = generate_seq(
                num_cmpnts, seed=actual_p['seed'], array_to_be_returned=2)
        else:
            indices_permuted1, rand_val = generate_seq(
                num_cmpnts, array_to_be_returned=2)

        # Initialize the output dictionary of the model
        out = {}
        out[self.parameter_name] = np.zeros(num_cmpnts)
        out['run'] = 2*np.ones(num_cmpnts)
        self.well_status_groups = {k: [] for k in range(1, 4)}

        # Calculate number of wells in each leak group
        leak_groups_counts = {
            k: int(np.floor(self.leakage_type_perc[k]*num_cmpnts)) for k in range(1, 3)}
        leak_groups_counts[3] = num_cmpnts - leak_groups_counts[1] - leak_groups_counts[2]

        # Generate permeability values for all wells
        par_values = np.zeros(num_cmpnts)

        par_values[0:leak_groups_counts[1]] = self.perm_bounds['small'][0]+(
            self.perm_bounds['small'][1]
            -self.perm_bounds['small'][0])*rand_val[0:leak_groups_counts[1]]

        par_values[leak_groups_counts[1]:(leak_groups_counts[1]+leak_groups_counts[2])] = (
            self.perm_bounds['medium'][0]+(
                self.perm_bounds['medium'][1]
                -self.perm_bounds['medium'][0])*rand_val[
                    leak_groups_counts[1]:(leak_groups_counts[1]+leak_groups_counts[2])])

        par_values[(leak_groups_counts[1]+leak_groups_counts[2]):] = self.perm_bounds['large'][0]+(
            self.perm_bounds['large'][1]
            -self.perm_bounds['large'][0])*rand_val[
                (leak_groups_counts[1]+leak_groups_counts[2]):]

        if self.configuration == 1:  # permeability does not correlate with scores
            out[self.parameter_name] = par_values[indices_permuted1-1]

        elif self.configuration == 2:  # permeability correlates with scores
            # but only in terms of distribution
            # Return indices where permeability is greater
            # than the higher end of small leak permeability (1.0e-19)
            perm_med_or_large_ind = np.where(
                par_values > self.perm_bounds['small'][1])[0]

            # Assign all wells permeability from the smallest leak permeability
            out[self.parameter_name] = self.perm_bounds['small'][0]+(
                self.perm_bounds['small'][1]-self.perm_bounds['small'][0])*rand_val

            # Return indices that would sort scores
            scores_sort_ind = np.argsort(-self.scores)[0:len(perm_med_or_large_ind)]

            # Assign medium or large permeability to the wells with larger scores
            out[self.parameter_name][scores_sort_ind] = par_values[perm_med_or_large_ind]

        elif self.configuration == 3:  # permeability correlates with scores
            # Sort permeability
            par_values[::-1].sort()

            # Return indices that would sort scores
            scores_sort_ind = np.argsort(-self.scores)

            # Assign the wells with larger leak scores the larger permeability
            out[self.parameter_name][scores_sort_ind] = par_values

#            # Check whether the configurer produced what is needed
#            # Return indices that would sort scores
#            scores_sort_ind = np.argsort(-self.scores)[0:60]
#            # Print scores and assigned perm
#            print(self.scores[scores_sort_ind], out[self.parameter_name][scores_sort_ind])

        # Save what wells are in what group
        self.details['perm_groups'] = {}
        self.details['perm_groups']['small'] = np.intersect1d(
            np.where(out[self.parameter_name] >= self.perm_bounds['small'][0])[0],
            np.where(out[self.parameter_name] <= self.perm_bounds['small'][1])[0])

        self.details['perm_groups']['medium'] = np.intersect1d(
            np.where(out[self.parameter_name] >= self.perm_bounds['medium'][0])[0],
            np.where(out[self.parameter_name] <= self.perm_bounds['medium'][1])[0])

        self.details['perm_groups']['large'] = np.intersect1d(
            np.where(out[self.parameter_name] >= self.perm_bounds['large'][0])[0],
            np.where(out[self.parameter_name] <= self.perm_bounds['large'][1])[0])

        self.well_status_groups[1] = len(self.details['perm_groups']['small'])
        self.well_status_groups[2] = len(self.details['perm_groups']['medium'])
        self.well_status_groups[3] = len(self.details['perm_groups']['large'])

#            # Turn off wells with small permeability
#            for ind, cind in enumerate(self.details['perm_groups']['small']):
#                nm = self.cmpnts_to_configure[cind]
#                cmpnt = self._parent.component_models[nm]
#                cmpnt.run_frequency = 0
#                out['run'][cind] = 0

        return out


def test_case_parameter_setup1():
    # Import needed packages and classes
    from matk import pyDOE
    from openiam import SystemModel, AnalyticalReservoir, MultisegmentedWellbore

    # Define output directory
    output_directory = os.path.join('..', '..', 'output',
                                    'test_parameter_setup1')
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # Define keyword arguments of the system model
    num_years = 10
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}   # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Create randomly located leaky well locations
    # within box defined by xmin,xmax,ymin,ymax
    xymins = np.array([100., 450.])
    xymaxs = np.array([300., 500.])
    num_wells = 10
    well_xys = xymins + pyDOE.lhs(2, samples=num_wells)*(xymaxs-xymins)

    ress_nms = ['res'+str(ind+1) for ind in range(num_wells)]
    mss_nms = ['ms'+str(ind+1) for ind in range(num_wells)]

    param_stp = sm.add_component_model_object(
            ParameterSetup1(name='stp1', parent=sm, cmpnt_nms=mss_nms))
    param_stp.add_grid_obs('logWellPerm', constr_type='array',
                           output_dir=output_directory, index=[0])
    param_stp.add_obs_to_be_linked('logWellPerm', obs_type='grid', constr_type='array')

    ress = []
    mss = []
    for ind, crds in enumerate(well_xys):

        # Add reservoir components
        ress.append(sm.add_component_model_object(
                AnalyticalReservoir(name=ress_nms[ind], parent=sm,
                injX=0., injY=0., locX=crds[0], locY=crds[1])))

        # Add parameters of reservoir component model
        ress[-1].add_par('numberOfShaleLayers', value=3, vary=False)
        ress[-1].add_par('injRate', value=0.8, vary=False)
        ress[-1].add_par('shale1Thickness', min=30.0, max=50., value=40.0)
        ress[-1].add_par('aquifer1Thickness', min=20.0, max=60., value=50.0)
        ress[-1].add_par('aquifer2Thickness', min=30.0, max=50., value=45.0)
        ress[-1].add_par('logResPerm', min=-13.,max=-11., value=-12.)

        # Add observations of reservoir component model to be used by the next component
        ress[-1].add_obs_to_be_linked('pressure')
        ress[-1].add_obs_to_be_linked('CO2saturation')

        # Add observations of reservoir component model
        ress[-1].add_obs('pressure')
        ress[-1].add_obs('CO2saturation')

        # Add multisegmented wellbore components
        mss.append(sm.add_component_model_object(
            MultisegmentedWellbore(name=mss_nms[ind], parent=sm)))

        # A lot of parameters of multisegmented wellbore component
        # are the same as for the reservoir component
        # Add parameters linked to the same parameters from reservoir model
        mss[-1].add_par_linked_to_par(
            'numberOfShaleLayers', ress[-1].deterministic_pars['numberOfShaleLayers'])
        mss[-1].add_par_linked_to_par(
            'shale1Thickness', ress[-1].pars['shale1Thickness'])
        mss[-1].add_par_linked_to_par(
            'shale2Thickness', ress[-1].default_pars['shaleThickness'])
        mss[-1].add_par_linked_to_par(
            'shale3Thickness', ress[-1].default_pars['shaleThickness'])
        mss[-1].add_par_linked_to_par(
            'aquifer1Thickness', ress[-1].pars['aquifer1Thickness'])
        mss[-1].add_par_linked_to_par(
            'aquifer2Thickness', ress[-1].pars['aquifer2Thickness'])
        mss[-1].add_par_linked_to_obs(
            'logWellPerm', param_stp.linkobs['logWellPerm'],
            obs_type='grid', loc_ind=ind, constr_type='array')
        # Add keyword arguments linked to the output provided by reservoir model
        mss[-1].add_kwarg_linked_to_obs(
            'pressure', ress[-1].linkobs['pressure'])
        mss[-1].add_kwarg_linked_to_obs(
            'CO2saturation', ress[-1].linkobs['CO2saturation'])
        mss[-1].add_obs('brine_aquifer1')
        mss[-1].add_obs('CO2_aquifer1')
        mss[-1].details['risk_type'] = 'exist_path'
        mss[-1].details['normalized_leak_score'] = None

    sm.forward()

    # Extract parameter setup component outputs
    permeability = sm.collect_gridded_observations_as_time_series(
        param_stp, 'logWellPerm', output_directory, indices=[0], rlzn_number=0)[0]

    print('Generated permeability \n', permeability)


def test_case_parameter_setup2():
    pass


def test_case_parameter_setup3():
    pass

def test_case_parameter_setup4():
    pass


def test_case_parameter_setup5():
    pass



if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)

    test_case = 3

    test_to_run = {1: test_case_parameter_setup1,
                   2: test_case_parameter_setup2,
                   3: test_case_parameter_setup3,
                   4: test_case_parameter_setup4,
                   5: test_case_parameter_setup5}

    # Run corresponding test example
    test_to_run[test_case]()
