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
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from openiam import ComponentModel
    from openiam.parameter_setup_component import generate_seq
except ImportError as err:
    print('Unable to load IAM class module: '+str(err))


class PressureBasedRiskConfigurer(ComponentModel):
    """
    Component that uses pressure at the bottom of the well to decide whether
    a particular well is going to leak or not.
    """
    def __init__(self, name, parent, cmpnt_nms=None, reproducible=True):
        """
        Constructor method of PressureBasedRiskConfigurer class.
        This configurer does not change the permeability of the components.
        It only controls when the leakage does start.

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to; it's the same parent that components to be configured
            are linked to.
        :type parent: SystemModel object

        :param cmpnt_nms: list of components whose parameters/keyword
            arguments/attributes will be controlled by the configurer
        type cmpnt_nms: list()

        :returns: PressureBasedRiskConfigurer class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = dict()
        model_kwargs['time_point'] = 365.25   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'PressureBasedRiskConfigurer'

        # Setup the list of names of components which can be configured during the run
        if cmpnt_nms is None:
            self.cmpnts_to_configure = []
        else:
            self.cmpnts_to_configure = cmpnt_nms

        # Set a flag indicating whether for each new simulation
        # the results will be random or reproducible.
        # reproducible value of True indicates that the model method will use
        # either the default or provided seed parameter for random number generator
        # (so the results can be reproduced if the same seed is used).
        self.reproducible = reproducible

        # Set default parameters of the component model
        self.add_default_par('pressure_threshold', value=36.0e+6)
        self.add_default_par('seed', value=1)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['pressure_threshold'] = [20.0e+6, 80.0e+6]

        msg = 'PressureBasedRiskConfigurer created with name {name}'.format(
            name=self.name)
        logging.debug(msg)

    @property
    def num_cmpnts_on(self):
        if len(self.cmpnts_to_configure) > 0:
            return sum([1 for nm in self.cmpnts_to_configure
                        if self._parent.component_models[nm].run_frequency == 2])

        logging.error('No components were linked for configuration.')
        raise ValueError(' '.join(['List of components to be configured',
                                   'is empty.']))
    @property
    def num_cmpnts(self):
        return len(self.cmpnts_to_configure)

    def create_details(self):
        """
        Create details of the wellbore components needed for the work of the configurer.
        """
        details = {'leak_prob': None,
                   'leak_prob_kwargs': None,
                   'leak_start_time': None,
                   'detection_prob': None,
                   'detection_prob_kwargs': {},
                   'leak_detected_time': None,
                   'pressure_threshold': None,
                   'risk_type': None,
                   'location': None,
                   'normalized_leak_score': None}
        return details

    def simulation_model(self, p, time_point=365.25, pressure=0.0):
        """
        Initiate leakage of the well based on pressure exceeding the threshold.
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Initialize the output dictionary of the model
        out = {}
        if len(self.cmpnts_to_configure) == 0:
            logging.error('No components were linked for configuration.')
            raise ValueError(
                ' '.join(['The model method of component {}',
                          'cannot create a leakage configuration: the list of components',
                          'to be configured is empty']).format(self.name))

        if time_point == 0.0:
            # Zero value of status means nothing changed in the status of the wells
            self.details['status'] = np.zeros(len(self.cmpnts_to_configure))
            out['status'] = self.details['status']
            for ind, nm in enumerate(self.cmpnts_to_configure):
                cmpnt = self._parent.component_models[nm]

                # Check which of the components have no pressure thresholds assigned
                if cmpnt.details['pressure_threshold'] is None:
                    cmpnt.details['pressure_threshold'] = actual_p['pressure_threshold']

                # Check which of the components have no scores assigned
                if cmpnt.details['normalized_leak_score'] is None:
                    cmpnt.details['normalized_leak_score'] = 0.0
        else:
            if self.reproducible:    # if reproducible results are needed
                # Generate sequences with defined seeds
                leak_val = generate_seq(
                    self.num_cmpnts,
                    seed=actual_p['seed']+time_point+1,
                    array_to_be_returned=1)

                detect_val = generate_seq(
                    self.num_cmpnts,
                    seed=actual_p['seed']+time_point+10,
                    array_to_be_returned=1)
            else:
                # Otherwise, generate sequences randomly
                leak_val = generate_seq(self.num_cmpnts, array_to_be_returned=1)
                detect_val = generate_seq(self.num_cmpnts, array_to_be_returned=1)

            for ind, nm in enumerate(self.cmpnts_to_configure):
                cmpnt = self._parent.component_models[nm]
                if (cmpnt.run_frequency == 0) and (cmpnt.details['risk_type'] == 'TBD'):
                    # If condition is met
                    if pressure[ind] >= cmpnt.details['pressure_threshold']:
                        # Check whether the event indeed has occurred
                        if leak_val[ind] < cmpnt.details['leak_prob'](
                                cmpnt, time_point, pressure[ind],
                                **cmpnt.details['leak_prob_kwargs']):
                            self.details['status'][ind] = 2
                            cmpnt.run_frequency = 2
                            cmpnt.details['risk_type'] = 'exist_path'

                elif cmpnt.run_frequency == 2:
                    # Check whether the leak was detected
                    if (cmpnt.details['leak_detected_time'] is not None) or (
                            detect_val[ind] < cmpnt.details['detection_prob'](
                                cmpnt, time_point, **cmpnt.details['detection_prob_kwargs'])):
                        self.details['status'][ind] = 3
                        cmpnt.run_frequency = 0
                        cmpnt.details['risk_type'] = 'fixed'

        out['status'] = self.details['status']

        out['num_cmpnts_on'] = self.num_cmpnts_on
        return out


class DataBasedRiskConfigurer(ComponentModel):
    """
    Component that uses pressure at the bottom of the well and its risk score
    to decide whether the well is going to leak or not.
    """
    def __init__(self, name, parent, cmpnt_nms=None, reproducible=True,
                 perm_bounds=None):
        """
        Constructor method of DataBasedRiskConfigurer class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to; it's the same parent that components to be configured
            are linked to.
        :type parent: SystemModel object

        :param cmpnt_nms: list of components whose parameters/keyword
            arguments/attributes will be controlled by the configurer
        type cmpnt_nms: list()

        :param reproducible: flag variable indicating whether the produced configuration
            should be reproducible, i.e. if seed parameter of the model method
            should be used to generate the sequence of random numbers. By default,
            the value is True. By default, the value is True.
        :type reproducible: boolean

        :param perm_bounds: dictionary which contain the information
            of the permeability distribution for each group of wells based
            on the type of leakage. By default, {'small': [-19.0, -19.0],
                                                 'medium': [-17.0, -14.0],
                                                 'large': [-13.0, -12.0]}.
        :type perm_perc: dict()

        :returns: DataBasedRiskConfigurer class object

        This configurer does not change the permeability of the components.
        It only controls when the leakage does start.
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = dict()
        model_kwargs['time_point'] = 365.25   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'DataBasedRiskConfigurer'

        # Setup the list of names of components which can be configured during the run
        if cmpnt_nms is None:
            self.cmpnts_to_configure = []
        else:
            self.cmpnts_to_configure = cmpnt_nms

        # Save distribution properties
        if perm_bounds is None:
            self.perm_bounds = {'small': [-19.0, -19.0],
                                'medium': [-17.0, -14.0],
                                'large': [-13.0, -12.0]}
        else:
            self.perm_bounds = perm_bounds

        # Set a flag indicating whether for each new simulation
        # the results will be random or reproducible.
        # reproducible value of True indicates that the model method will use
        # either the default or provided seed parameter for random number generator
        # (so the results can be reproduced if the same seed is used).
        self.reproducible = reproducible

        # Set default parameters of the component model
        self.add_default_par('pressure_threshold', value=36.0e+6)
        self.add_default_par('score_threshold', value=0.0)
        self.add_default_par('seed', value=1)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['pressure_threshold'] = [20.0e+6, 80.0e+6]
        self.pars_bounds['score_threshold'] = [0.0, 1.0]

        msg = 'DataBasedRiskConfigurer created with name {name}'.format(
            name=self.name)
        logging.debug(msg)

    @property
    def num_cmpnts_on(self):
        if not self.cmpnts_to_configure:
            logging.error('No components were linked for configuration.')
            raise ValueError(' '.join(['List of components to be configured',
                                       'is empty.']))

        return sum([1 for nm in self.cmpnts_to_configure
                        if self._parent.component_models[nm].run_frequency == 2])

    @property
    def num_cmpnts(self):
        return len(self.cmpnts_to_configure)

    def create_details(self):
        """
        Create details of the wellbore components needed for the work of the configurer.
        """
        details = {'leak_prob': None,
                   'leak_prob_kwargs': {'score': None},
                   'leak_start_time': None,
                   'detection_prob': None,
                   'detection_prob_kwargs': {},
                   'leak_detected_time': None,
                   'leak_started_time': None,
                   'pressure_threshold': None,
                   'risk_type': None,
                   'location': None,
                   'logWellPerm': None,
                   'normalized_leak_score': None}
        return details

    def simulation_model(self, p, time_point=365.25, pressure=0.0, logWellPerm=0.0):
        """
        Initiate leakage of the well based on pressure exceeding the threshold.
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Initialize the output dictionary of the model
        out = {}
        if len(self.cmpnts_to_configure) == 0:
            logging.error('No components were linked for configuration.')
            raise ValueError(
                ' '.join(['The model method of component {}',
                          'cannot create a leakage configuration: the list of components',
                          'to be configured is empty']).format(self.name))

        if time_point == 0.0:
            # zero means nothing changed in the status of the wells
            self.details['status'] = np.zeros(len(self.cmpnts_to_configure))
            out['status'] = self.details['status']
            num_leak_events = 0

            # Check whether initial permeability data were provided
            if len(logWellPerm) == 0:
                out['logWellPerm'] = -19.0 + 7*np.random.rand(self.num_cmpnts)
            else:
                out['logWellPerm'] = np.asarray(logWellPerm)

            for ind, nm in enumerate(self.cmpnts_to_configure):
                cmpnt = self._parent.component_models[nm]

                # Check which of the components have no pressure thresholds assigned
                if cmpnt.details['pressure_threshold'] is None:
                    cmpnt.details['pressure_threshold'] = actual_p['pressure_threshold']

                if cmpnt.details['logWellPerm'] is None:
                    cmpnt.details['logWellPerm'] = out['logWellPerm'][ind]
                else:
                    if (len(logWellPerm) == 0) and (
                            cmpnt.details['logWellPerm'] != out['logWellPerm'][ind]):
                        logging.warning(
                            ''.join(['Permeability value {} passed',
                                     ' to the pressure based leakage configurer',
                                     ' does not coincide with the value {} setup',
                                     ' in the wellbore component with index {}.']).format(
                                         out['logWellPerm'][ind],
                                         cmpnt.details['logWellPerm'],
                                         ind))
                        cmpnt.details['logWellPerm'] = out['logWellPerm'][ind]

        else:
            num_leak_events = 0
            out['logWellPerm'] = np.array(
                [self._parent.component_models[nm].details[
                    'logWellPerm'] for ind, nm in enumerate(self.cmpnts_to_configure)])
            if self.reproducible:    # if reproducible results are needed
                # Generate sequences with defined seeds
                leak_val = generate_seq(
                    self.num_cmpnts,
                    seed=actual_p['seed']+time_point+1,
                    array_to_be_returned=1)
                detect_val = generate_seq(
                    self.num_cmpnts,
                    seed=actual_p['seed']+time_point+10,
                    array_to_be_returned=1)
                generate_val = generate_seq(
                    self.num_cmpnts,
                    seed=actual_p['seed']+time_point+100,
                    array_to_be_returned=1)
            else:
                # Otherwise, generate sequences randomly
                leak_val = generate_seq(self.num_cmpnts, array_to_be_returned=1)
                detect_val = generate_seq(self.num_cmpnts, array_to_be_returned=1)
                generate_val = generate_seq(self.num_cmpnts, array_to_be_returned=1)

            for ind, nm in enumerate(self.cmpnts_to_configure):
                cmpnt = self._parent.component_models[nm]

                if (cmpnt.details['leak_detected_time'] is not None) or (
                        detect_val[ind] < cmpnt.details['detection_prob'](
                            cmpnt, time_point, **cmpnt.details['detection_prob_kwargs'])):
                    cmpnt.details['leak_started_time'] = None
                    self.details['status'][ind] = 3    # leak is detected
                    out['logWellPerm'][ind] = -19.0    # remediate well
                    cmpnt.details['logWellPerm'] = out['logWellPerm'][ind]
                    cmpnt.details['leak_prob_kwargs']['score'] = 0.2
                    cmpnt.details['normalized_leak_score'] = 0.2

                # If the leakage was not detected but the pressure exceeds the threshold
                elif ((cmpnt.details['normalized_leak_score'] >= actual_p['score_threshold']) and (
                        pressure[ind] >= cmpnt.details['pressure_threshold'])):
                    # Check whether the event indeed has occurred
                    if leak_val[ind] < cmpnt.details['leak_prob'](
                            cmpnt, time_point, pressure[ind],
                            **cmpnt.details['leak_prob_kwargs']):

                        if cmpnt.details['leak_started_time'] is None:
                            self.details['status'][ind] = 2  # leak has occurred
                            cmpnt.details['leak_detected_time'] = None
                            # Record time
                            cmpnt.details['leak_started_time'] = time_point
                            # Increase permeability sampling from larger permeability group
                            if cmpnt.details['logWellPerm'] <= self.perm_bounds['small'][1]:
                                out['logWellPerm'][ind] = self.perm_bounds['medium'][0] + (
                                    self.perm_bounds['medium'][1]
                                    -self.perm_bounds['medium'][0])*generate_val[ind]
                            elif cmpnt.details['logWellPerm'] <= self.perm_bounds['medium'][1]:
                                out['logWellPerm'][ind] = self.perm_bounds['large'][0] + (
                                    self.perm_bounds['large'][1]
                                    -self.perm_bounds['large'][0])*generate_val[ind]
                            else:
                                val = self.perm_bounds['large'][0] + (
                                    self.perm_bounds['large'][1]
                                    -self.perm_bounds['large'][0])*generate_val[ind]
                                out['logWellPerm'][ind] = max(
                                    [val, cmpnt.details['logWellPerm']])
                            num_leak_events = num_leak_events + 1

                        cmpnt.details['logWellPerm'] = out['logWellPerm'][ind]

                # maxp = np.max(
                #     np.array([
                #         self._parent.component_models[nm].details[
                #             'logWellPerm'] for ind, nm in enumerate(
                #                 self.cmpnts_to_configure)]))
                # print(np.max(out['logWellPerm']), maxp)

            out['status'] = self.details['status']

        out['num_leak_events'] = num_leak_events
        out['num_cmpnts_on'] = self.num_cmpnts_on
        return out


class WellDepthRiskConfigurer(ComponentModel):
    """
    Component that compares the depth of a well to the depth of the reservoir
    to decide whether the well can be considered for leak scenario or not.
    """
    def __init__(self, name, parent, cmpnt_nms=None):
        """
        Constructor method of WellDepthRiskConfigurer class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to; it's the same parent that components to be configured
            are linked to.
        :type parent: SystemModel object

        :param cmpnt_nms: list of components whose parameters/keyword
            arguments/attributes will be controlled by the configurer
        type cmpnt_nms: list() or list of list() for multiple types of components
            being handled
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = dict()
        model_kwargs['time_point'] = 365.25   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'WellDepthRiskConfigurer'

        # Setup the list of names of components which can be configured during the run
        if cmpnt_nms is None:
            self.cmpnts_to_configure = np.array([])
        else:
            self.cmpnts_to_configure = np.array(cmpnt_nms)
            if len(self.cmpnts_to_configure.shape) == 1:
                self.cmpnts_to_configure = self.cmpnts_to_configure.reshape((1, -1))
            self.num_comp_types = self.cmpnts_to_configure.shape[0]

        # Set default parameters of the component model
        self.add_default_par('reservoirDepth', value=2000)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['reservoirDepth'] = [5.0, 30000.0]

        msg = 'WellDepthRiskConfigurer created with name {}'.format(self.name)
        logging.debug(msg)

        # The model only needs to be run once, as output for all times can be
        # obtained at once.
        self.default_run_frequency = 1
        self.run_frequency = 1

    @property
    def num_cmpnts_on(self):
        if self.cmpnts_to_configure.size == 0:
            logging.error('No components were linked for configuration.')
            raise ValueError(' '.join(['List of components to be configured',
                                       'is empty.']))

        return sum([1 for nm in self.cmpnts_to_configure.flatten()
                        if self._parent.component_models[nm].run_frequency != 0])

    @property
    def num_cmpnts(self):
        return self.cmpnts_to_configure.size

    def simulation_model(self, p, time_point=365.25, wellDepth=None, **kwargs):
        """
        Initiate leakage of the well based on wellDepth exceeding the reservoirDepth.
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        if wellDepth is not None:
            # The locZ values from a LocationGenerator component will be negative.
            # To be consistent with the conventions of other components, however,
            # the reservoirDepth parameter is taken as positive (e.g., the
            # OpenWellbore component's reservoirDepth parameter and the CementedWellbore's
            # wellDepth parameter are taken as positive). If wellDepth is negative,
            # it is made positive - we never focus on elevations above the Earth's surface.
            well_depth = np.abs(wellDepth)
        else:
            err_msg = "".join([
                "Argument 'wellDepth' is not defined (None). ",
                "WellDepthRiskConfigurer {} cannot proceed with the ",
                "setup."]).format(self.name)
            logging.error(err_msg)
            raise ValueError(err_msg)

        # Get reservoir depth
        reservoir_depth = kwargs.get('reservoirDepth', actual_p['reservoirDepth'])
        if isinstance(reservoir_depth, (int, float)):  # scalar for a reservoir depth
            reservoir_depth = reservoir_depth*np.ones(self.num_cmpnts//self.num_comp_types)
        else:
            if len(reservoir_depth) != self.num_cmpnts//self.num_comp_types:
                err_msg = "".join([
                    "Number of elements in the array keyword argument 'reservoirDepth' ",
                    "is not equal to the number of components to be configured. ",
                    "Check your input."])
                logging.error(err_msg)
                raise ValueError(err_msg)

        # Initialize the output dictionary of the model
        out = {}

        if self.cmpnts_to_configure.size == 0:
            logging.error('No components were linked for configuration.')
            raise ValueError(
                ' '.join(['The model method of component {}',
                          'cannot create a leakage configuration: the list of components',
                          'to be configured is empty.']).format(self.name))

        less_indices = np.where(well_depth < reservoir_depth)[0]
        for type_ind in range(self.num_comp_types):
            for ind in less_indices:
                cmpnt = self._parent.component_models[self.cmpnts_to_configure[type_ind, ind]]
                cmpnt.run_frequency = 0

        out['num_cmpnts_on'] = self.num_cmpnts_on

        return out


def test_case_pressure_based_risk_configurer():
    pass


def test_case_data_based_risk_configurer():
    pass

def test_case_well_depth_risk_configurer():
    pass


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)

    test_case = 2

    test_to_run = {1: test_case_pressure_based_risk_configurer,
                   2: test_case_data_based_risk_configurer}

    # Run corresponding test example
    test_to_run[test_case]()
