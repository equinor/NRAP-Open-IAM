"""
This module contains several monitoring tools classes.

Created: September, 2018
Last modified: June, 2023

@author: Veronika Vasylkivska
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


class MonitoringScheduler1(ComponentModel):
    """
    Monitoring Scheduler 1 component produces indices of the components to be
    monitored based on a prescribed schedule and whether the wellbores have
    associated pressure values exceeding the defined threshold.

    Monitoring scheduler component should precede the monitoring tool,
    wellbore components, or any other components it provides input
    for in the setup of the system model.
    In the case the component is created after the to be scheduled components
    one of the system model methods (reorder_component_models, swap_components,
    or put_in_front_of_component) should be used to correct the order of all
    components.
    """
    def __init__(self, name, parent, cmpnt_nms=None, reproducible=True):
        """
        Constructor method of MonitoringScheduler1 class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to; it's the same parent that components to be configured
            are linked to.
        :type parent: SystemModel object

        :param cmpnt_nms: list of components whose parameters/keyword
            arguments/attributes will be checked by the monitoring tool
        :type cmpnt_nms: list()

        :param reproducible: flag variable indicating whether the produced setup
            should be reproducible, i.e. if seed parameter of the model method
            should be used to generate the sequence of random numbers.
        :type reproducible: boolean

        :returns: MonitoringScheduler1 class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25,   # default value of 365.25 days
                        'time_index': 0}        # default value of 0, single point

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'MonitoringScheduler1'

        # Setup the list of names of components that can be monitored
        if cmpnt_nms is None:
            self.cmpnts_to_monitor = []
        else:
            self.cmpnts_to_monitor = cmpnt_nms

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

        # Add accumulators that would keep an updated list of indices of components
        # that still can be monitored
        self.add_accumulator('current_cmpnts_to_monitor', sim=0.0)

        msg = 'MonitoringScheduler1 created with name {name}'.format(name=self.name)
        logging.debug(msg)

    @property
    def max_num_cmpnts_to_monitor(self):
        try:
            val = len(self.cmpnts_to_monitor)
            return val
        except:
            raise TypeError('Attribute cmpnts_to_monitor is of wrong type.')

    def check_input_parameters(self, p):
        """
        Check setup of the monitoring component and input parameters.

        """
        # Check whether the component was linked to any components to be monitored
        if self.max_num_cmpnts_to_monitor == 0:
            warn_msg = ''.join([
                'No components are linked to the monitoring scheduler ',
                'component {}'.format(self.name)])
            logging.warning(warn_msg)

    def simulation_model(self, p, time_point=365.25, time_index=0, pressure=0.0,
                         num_cmpnts_to_check=5):
        """
        :param p: dictionary of input parameters
        :type p: dict

        :param time_point: time point (in days) for which the monitoring is
            done; by default, its value is 365.25 days
        :type time_point: float

        :param time_index: index of time point defined in the system time_array attribute
        :type time_index: int

        :param pressure: list of pressure values at the monitored wells
            to be compared with the pressure threshold
        :type pressure: np.array of floats

        :param num_cmpnts_to_check: parameter defines frequency of testing
        :type num_cmpnts_to_check: int

        Setups schedule of testing the components.
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        out = {}
        if time_point == 0.0:
            out['cmpnts_indices_to_check'] = np.array([])
            self.accumulators['current_cmpnts_to_monitor'].sim = np.arange(
                len(self.cmpnts_to_monitor))
        else:
            # Check whether the components to be monitored are linked
            if self.max_num_cmpnts_to_monitor > 0:
                if num_cmpnts_to_check > 0:
                    # Determine which of the components have pressure exceeding threshold
                    init_indices = np.where(pressure >= actual_p['pressure_threshold'])[0]
                    # If there is something to check according to the criteria
                    # Compare the generated list with the indices of the components
                    # that need to be monitored
                    if len(init_indices) > 0:
                        common_indices = np.array(
                            [ind for ind in init_indices if ind in self.accumulators[
                                'current_cmpnts_to_monitor'].sim])
                    else: # if not, then all eligible indices can be used
                        common_indices = self.accumulators['current_cmpnts_to_monitor'].sim

                    # Determine maximum number of components that can be tested
                    max_num_cmpnts_to_monitor = len(common_indices)

                    # Generate random indices of components to check
                    num_indices = min(num_cmpnts_to_check, max_num_cmpnts_to_monitor)
                    if self.reproducible:
                        indices = generate_seq(
                            max_num_cmpnts_to_monitor,
                            seed=actual_p['seed'] + 8.7654*time_point + 1, # random constant
                            array_to_be_returned=0)[0:num_indices]
                    else:
                        indices = generate_seq(max_num_cmpnts_to_monitor,
                                               array_to_be_returned=0)[0:num_indices]
                    indices = indices-1
#                    print(common_indices)
                    cmpnts_indices = common_indices[indices]

                    # Update accumulator variable
                    self.accumulators['current_cmpnts_to_monitor'].sim = np.setdiff1d(
                        self.accumulators['current_cmpnts_to_monitor'].sim,
                        cmpnts_indices)
                else:
                    cmpnts_indices = np.array([])  # no components are checked


            out['cmpnts_indices_to_check'] = cmpnts_indices

        return out


class MonitoringScheduler2(ComponentModel):
    """
    Monitoring Scheduler 2 component produces indices of the components to be
    monitored based on a prescribed schedule. In selection of the next wellbore
    components to monitor the preference is given to the wellbores with higher
    risk scores.

    Monitoring scheduler component should precede the monitoring tool,
    wellbore components, or any other components it provides input
    for in the setup of the system model.
    In the case the component is created after the to be scheduled components
    one of the system model methods (reorder_component_models, swap_components,
    or put_in_front_of_component) should be used to correct the order of all
    components.
    """
    def __init__(self, name, parent, cmpnt_nms=None):
        """
        Constructor method of MonitoringScheduler2 class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to; it's the same parent that components to be configured
            are linked to.
        :type parent: SystemModel object

        :param cmpnt_nms: list of components whose parameters/keyword
            arguments/attributes will be checked by the monitoring tool
        :type cmpnt_nms: list()

        :returns: MonitoringScheduler2 class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25,   # default value of 365.25 days
                        'time_index': 0}        # default value of 0, single point

        super().__init__(name, parent, model=self.simulation_model, model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'MonitoringScheduler2'

        # Setup the list of names of components that can be monitored
        if cmpnt_nms is None:
            self.cmpnts_to_monitor = []
        else:
            self.cmpnts_to_monitor = cmpnt_nms

        # Set default parameters of the component model
        self.add_default_par('pressure_threshold', value=36.0e+6)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['pressure_threshold'] = [20.0e+6, 80.0e+6]

        # Add accumulators that would keep an updated list of indices
        # of components that still can be monitored
        self.add_accumulator('current_cmpnts_to_monitor', sim=0.0)
        self.add_accumulator('initial_pressure', sim=0.0)

        msg = 'MonitoringScheduler2 created with name {name}'.format(name=self.name)
        logging.debug(msg)

    @property
    def max_num_cmpnts_to_monitor(self):
        try:
            val = len(self.cmpnts_to_monitor)
            return val
        except:
            raise TypeError('Attribute cmpnts_to_monitor is of wrong type.')

    def check_input_parameters(self, p):
        """
        Check setup of the monitoring component and input parameters.

        """
        # Check whether the component was linked to any components to be monitored
        if self.max_num_cmpnts_to_monitor == 0:
            warn_msg = ''.join([
                'No components are linked to the monitoring scheduler ',
                'component {}'.format(self.name)])
            logging.warning(warn_msg)

    def simulation_model(self, p, time_point=365.25, time_index=0, pressure=0.0,
                         num_cmpnts_to_check=5):
        """
        :param p: dictionary of input parameters
        :type p: dict

        :param time_point: time point (in days) for which the monitoring is
            done; by default, its value is 365.25 days
        :type time_point: float

        :param time_index: index of time point defined in the system time_array attribute
        :type time_index: int

        :param num_cmpnts_to_check: parameter defines frequency of testing
        :type num_cmpnts_to_check: int

        Setups schedule of testing the components.
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        out = {}
        if time_point == 0.0:
            out['cmpnts_indices_to_check'] = np.array([])
            self.accumulators['current_cmpnts_to_monitor'].sim = np.arange(
                len(self.cmpnts_to_monitor))
            self.accumulators['initial_pressure'].sim = pressure  #  save initial pressure

        else:
            # Check whether the components to be monitored are linked
            if self.max_num_cmpnts_to_monitor > 0:
                if num_cmpnts_to_check > 0:

                    # Calculate number of components that still can be tested
                    num_cmpnts = len(self.accumulators['current_cmpnts_to_monitor'].sim)

                    # Get (possibly updated) scores
                    scores = np.zeros(num_cmpnts)
                    for ind, cind in enumerate(self.accumulators['current_cmpnts_to_monitor'].sim):
                        nm = self.cmpnts_to_monitor[cind]
                        cmpnt = self._parent.component_models[nm]
                        scores[ind] = cmpnt.details['normalized_leak_score']

                    # Return the indices that would sort scores in decreasing order
                    sorted_scores_ind = np.argsort(-scores)

                    # Determine the maximum number of components that can be tested
                    num_indices = min(num_cmpnts_to_check, num_cmpnts)

                    cmpnts_indices = self.accumulators[
                        'current_cmpnts_to_monitor'].sim[sorted_scores_ind][0:num_indices]

                    # Update accumulator variable
                    self.accumulators['current_cmpnts_to_monitor'].sim = np.setdiff1d(
                        self.accumulators['current_cmpnts_to_monitor'].sim,
                        cmpnts_indices)
                else:
                    cmpnts_indices = np.array([])  # no components are checked


            out['cmpnts_indices_to_check'] = cmpnts_indices

        return out


class MonitoringScheduler3(ComponentModel):
    """
    Monitoring Scheduler 3 component produces indices of the components to be
    monitored based on a prescribed schedule. The wellbores to be monitored
    at the current time step are chosen randomly.

    Monitoring scheduler component should precede the monitoring tool,
    wellbore components, or any other components it provides input
    for in the setup of the system model.
    In the case the component is created after the to be scheduled components
    one of the system model methods (reorder_component_models, swap_components,
    or put_in_front_of_component) should be used to correct the order of all
    components.
    """
    def __init__(self, name, parent, cmpnt_nms=None, cmpnt_inds=None,
                 reproducible=True, num_cmpnts_to_remediate=0):
        """
        Constructor method of MonitoringScheduler3 class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to; it's the same parent that components to be configured
            are linked to.
        :type parent: SystemModel object

        :param cmpnt_nms: list of components whose parameters/keyword
            arguments/attributes will be checked by the monitoring tool
        :type cmpnt_nms: list()

        :param reproducible: flag variable indicating whether the produced setup
            should be reproducible, i.e. if seed parameter of the model method
            should be used to generate the sequence of random numbers.
        :type reproducible: boolean

        :returns: MonitoringScheduler3 class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25,   # default value of 365.25 days
                        'time_index': 0}        # default value of 0, single point

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'MonitoringScheduler3'

        # Setup the list of names of components that can be monitored
        if cmpnt_nms is None:
            self.cmpnts_to_monitor = []
        else:
            self.cmpnts_to_monitor = cmpnt_nms

        if not cmpnt_inds:  # cover cases when cmpnt_inds is empty list or None
            self.cmpnts_to_monitor_inds = list(range(len(cmpnt_nms)))
        else:
            self.cmpnts_to_monitor_inds = cmpnt_inds
            if len(cmpnt_inds) != len(cmpnt_nms):
                raise TypeError(''.join(['Length of array with components names',
                                         ' does not coincide with the length of',
                                         ' array with components indices.']))

        self.num_cmpnts_to_remediate = num_cmpnts_to_remediate

        # Set a flag indicating whether for each new simulation
        # the results will be random or reproducible.
        # reproducible value of True indicates that the model method will use
        # either the default or provided seed parameter for random number generator
        # (so the results can be reproduced if the same seed is used).
        self.reproducible = reproducible

        # Set default parameters of the component model
        self.add_default_par('seed', value=1)

        # Add accumulators that would keep an updated list of indices
        # of components that still can be monitored
        self.add_accumulator('current_cmpnts_to_monitor', sim=0.0)

        msg = 'MonitoringScheduler3 created with name {name}'.format(name=self.name)
        logging.debug(msg)

    @property
    def max_num_cmpnts_to_monitor(self):
        try:
            val = len(self.cmpnts_to_monitor)
            return val
        except:
            raise TypeError('Attribute cmpnts_to_monitor is of wrong type.')

    def check_input_parameters(self, p):
        """
        Check setup of the monitoring component and input parameters.

        """
        # Check whether the component was linked to any components to be monitored
        if self.max_num_cmpnts_to_monitor == 0:
            warn_msg = ''.join([
                'No components are linked to the monitoring scheduler ',
                'component {}'.format(self.name)])
            logging.warning(warn_msg)

    def simulation_model(self, p, time_point=365.25, time_index=0,
                         num_cmpnts_to_check=5):
        """
        :param p: dictionary of input parameters
        :type p: dict

        :param time_point: time point (in days) for which the monitoring is
            done; by default, its value is 365.25 days
        :type time_point: float

        :param time_index: index of time point defined in the system time_array attribute
        :type time_index: int

        :param num_cmpnts_to_check: parameter defines frequency of testing
        :type num_cmpnts_to_check: int

        Setups schedule of testing the components.
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        out = {}
        out['cmpnts_indices_to_remediate'] = np.array([]) # for other cases
        if time_point == 0.0:
            out['cmpnts_indices_to_check'] = np.array([])
            self.accumulators['current_cmpnts_to_monitor'].sim = np.array(
                self.cmpnts_to_monitor_inds)

        else:

            if time_index == 1:
                if self.num_cmpnts_to_remediate > 0:
                    if self.reproducible:
                        indices = generate_seq(
                            self.max_num_cmpnts_to_monitor,
                            seed=actual_p['seed'] + 7.6543*time_point + 1, # random constant
                            array_to_be_returned=0)[0:self.num_cmpnts_to_remediate]-1
                    else:
                        indices = generate_seq(
                            self.max_num_cmpnts_to_monitor,
                            array_to_be_returned=0)[0:self.num_cmpnts_to_remediate]-1
                out['cmpnts_indices_to_remediate'] = np.array(
                    self.cmpnts_to_monitor_inds)[indices]
                self.accumulators['current_cmpnts_to_monitor'].sim = np.setdiff1d(
                    np.array(self.cmpnts_to_monitor_inds),
                    out['cmpnts_indices_to_remediate'])

            # Check whether the components to be monitored are linked
            if self.max_num_cmpnts_to_monitor > 0:
                if  num_cmpnts_to_check > 0:
                    avail_cmpnts_indices = self.accumulators['current_cmpnts_to_monitor'].sim
                    num_avail_cmpnts = len(avail_cmpnts_indices)
                    if num_avail_cmpnts > num_cmpnts_to_check:
                        # Randomize selection of components if more available than needed
                        if self.reproducible:
                            indices = generate_seq(
                                num_avail_cmpnts,
                                seed=actual_p['seed'] + 6.5432*time_point + 1, # random constant
                                array_to_be_returned=0)[0:num_cmpnts_to_check]-1
                        else:
                            indices = generate_seq(
                                num_avail_cmpnts,
                                array_to_be_returned=0)[0:num_cmpnts_to_check]-1

                        chosen_cmpnts_indices = avail_cmpnts_indices[indices]

                        # Update accumulator variable
                        self.accumulators['current_cmpnts_to_monitor'].sim = np.setdiff1d(
                            self.accumulators['current_cmpnts_to_monitor'].sim,
                            chosen_cmpnts_indices)

                    elif num_avail_cmpnts == num_cmpnts_to_check:
                        chosen_cmpnts_indices = avail_cmpnts_indices

                        # Update accumulator variable
                        self.accumulators['current_cmpnts_to_monitor'].sim = np.setdiff1d(
                            np.array(self.cmpnts_to_monitor_inds), chosen_cmpnts_indices)
                    else:
                        chosen_cmpnts_indices = np.zeros(num_cmpnts_to_check, dtype=int)
                        chosen_cmpnts_indices[0:num_avail_cmpnts] = avail_cmpnts_indices

                        # Randomize selection of components
                        if self.reproducible:
                            indices = generate_seq(
                                self.max_num_cmpnts_to_monitor-num_avail_cmpnts,
                                seed=actual_p['seed'] + 5.4321*time_point + 1, # random constant
                                array_to_be_returned=0)[0:(num_cmpnts_to_check-num_avail_cmpnts)]-1
                        else:
                            indices = generate_seq(
                                self.max_num_cmpnts_to_monitor-num_avail_cmpnts,
                                array_to_be_returned=0)[0:(num_cmpnts_to_check-num_avail_cmpnts)]-1

                        chosen_cmpnts_indices[num_avail_cmpnts:] = np.setdiff1d(
                            np.array(self.cmpnts_to_monitor_inds),
                            avail_cmpnts_indices)[indices]

                        # Update accumulator variable
                        self.accumulators['current_cmpnts_to_monitor'].sim = np.setdiff1d(
                            np.array(self.cmpnts_to_monitor_inds), chosen_cmpnts_indices)

                else:
                    chosen_cmpnts_indices = np.array([])  # no components are to be checked

            out['cmpnts_indices_to_check'] = chosen_cmpnts_indices

        return out


def test_case_monitoring_scheduler1():
    pass


def test_case_monitoring_scheduler2():
    pass


def test_case_monitoring_scheduler3():
    pass


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)

    test_case = 3

    test_to_run = {1: test_case_monitoring_scheduler1,
                   2: test_case_monitoring_scheduler2,
                   3: test_case_monitoring_scheduler3}

    # Run corresponding test example
    test_to_run[test_case]()
