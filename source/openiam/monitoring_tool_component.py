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


class MonitoringTool1(ComponentModel):
    """
    Monitoring Tool 1 component compares the distance between its associated
    location and a wellbore component and under consideration of the amount
    of leaked CO2 mass through the wellbore component makes a prediction
    about whether the leak can be detected or not.

    Monitoring tool component should precede the wellbore components or any other
    component it monitors in the setup of the system model.
    In the case the component is created after the to be monitored components
    one of the system model methods (reorder_component_models, swap_components,
    or put_in_front_of_component) should be used to correct the order of all
    components.
    """
    def __init__(self, name, parent, locX, locY, cmpnt_nms=None):
        """
        Constructor method of MonitoringTool1 class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to; it's the same parent that components to be configured
            are linked to.
        :type parent: SystemModel object

        :param locX: the x-coordinate of the monitoring well location
        :type locX: float

        :param locY: the y-coordinate of the monitoring well location
        :type locY: float

        :param cmpnt_nms: list of components whose parameters/keyword
            arguments/attributes will be checked by the monitoring tool
        type cmpnt_nms: list()

        :returns: MonitoringTool1 class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25}   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'MonitoringTool1'

        # Setup the list of names of components that can be monitored
        if cmpnt_nms is None:
            self.cmpnts_to_monitor = []
        else:
            self.cmpnts_to_monitor = cmpnt_nms

        # Setup monitoring tool location
        self.locX = locX
        self.locY = locY

        # Distances between monitoring well and wellbore locations
        self.distances = None

        # Set default parameters of the component model
        self.add_default_par('aquiferThickness', value=100.0)
        self.add_default_par('aquiferPorosity', value=1.0)
        self.add_default_par('brineDensity', value=1000.0)
        self.add_default_par('CO2Density', value=479.0)
        self.add_default_par('brineViscosity', value=2.535e-3)
        self.add_default_par('CO2Viscosity', value=3.95e-5)
        self.add_default_par('brineResSaturation', value=0.1)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['aquiferThickness'] = [1.0, 1000.0]
        self.pars_bounds['aquiferPorosity'] = [0.01, 1.0]
        self.pars_bounds['brineDensity'] = [900.0, 1500.0]
        self.pars_bounds['CO2Density'] = [100.0, 1500.0]
        self.pars_bounds['brineViscosity'] = [1.0e-4, 5.0e-3]
        self.pars_bounds['CO2Viscosity'] = [1.0e-6, 1.0e-4]
        self.pars_bounds['brineResSaturation'] = [0.0, 0.7]

    def simulation_model(self, p, time_point=365.25, CO2_mass=0.0):
        """
        Return number of wells identified as leaking.

        :param p: dictionary of input parameters
        :type p: dict

        :param time_point: time point (in days) for which the monitoring is
            done; by default, its value is 365.25 days
        :type time_point: float

        :param CO2_mass: mass of CO2 leaked to the monitored aquifer, [kg]
        :type CO2_mass: [float]

        :returns: dictionary of observations
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # Initialize the output dictionary of the model
        # Output can be number of wells detected to be leaking
        out = {}

        if time_point == 0.0:   # all wells have zero leakage detection probability at time 0
            # Find components locations and calculate distances between tool and wells
            self.distances = np.zeros(len(self.cmpnts_to_monitor))
            for ind, nm in enumerate(self.cmpnts_to_monitor):
                cmpnt = self._parent.component_models[nm]
                cmpnt_loc = cmpnt.details['location']
                self.distances[ind] = np.sqrt(
                    (cmpnt_loc[0]-self.locX)**2+(cmpnt_loc[1]-self.locY)**2)

            num_leak_wells = 0

            for ind, nm in enumerate(self.cmpnts_to_monitor):
                cmpnt = self._parent.component_models[nm]
                if cmpnt.run_frequency == 2:
                    num_leak_wells = num_leak_wells + 1

            out['num_ident_leak_wells'] = 0
            out['num_leak_wells'] = num_leak_wells
        else:
            scalar_const = np.sqrt(actual_p['CO2Viscosity']/(
                np.pi*actual_p['CO2Density']*actual_p['aquiferThickness']/3.*
                actual_p['aquiferPorosity']*(1-actual_p['brineResSaturation'])*
                actual_p['brineViscosity']))

            num_ident_leak_wells = 0
            num_leak_wells = 0

            for ind, nm in enumerate(self.cmpnts_to_monitor):
                cmpnt = self._parent.component_models[nm]
                effect_radius = scalar_const*np.sqrt(CO2_mass[ind])

                if cmpnt.run_frequency == 2:
                    num_leak_wells = num_leak_wells + 1
                    if effect_radius > self.distances[ind]:
                        num_ident_leak_wells = num_ident_leak_wells + 1
                        cmpnt.details['leak_detected_time'] = time_point

            out['num_ident_leak_wells'] = num_ident_leak_wells
            out['num_leak_wells'] = num_leak_wells

        # Return dictionary of outputs
        return out


class MonitoringTool2(ComponentModel):
    """
    Monitoring Tool 2 component checks whether any of the linked wellbore
    components have any integrity issues.

    Monitoring tool component should precede the wellbore components or any other
    component it monitors in the setup of the system model.
    In the case the component is created after the to be monitored components
    one of the system model methods (reorder_component_models, swap_components,
    or put_in_front_of_component) should be used to correct the order of all
    components.
    """
    def __init__(self, name, parent, cmpnt_nms=None, reproducible=True):
        """
        Constructor method of MonitoringTool2 class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to; it's the same parent that components to be configured
            are linked to.
        :type parent: SystemModel object

        :param cmpnt_nms: list of components whose parameters/keyword
            arguments/attributes will be checked by the monitoring tool
        type cmpnt_nms: list()

        :param reproducible: flag variable indicating whether the check of the wells
            should be reproducible, i.e. if seed parameter of the model method
            should be used to generate the sequence of random numbers.
        :type reproducible: boolean

        :returns: MonitoringTool2 class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25,   # default value of 365.25 days
                        'time_index': 0}        # default value of 0, single point

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'MonitoringTool2'

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
        self.add_default_par('seed', value=1)

        # Define attribute that will save the history of what components were checked
        self.details['cmpnts_indices'] = {}

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
                'No components are linked to the monitoring ',
                'component {}'.format(self.name)])
            logging.warning(warn_msg)

    def simulation_model(self, p, time_point=365.25, time_index=0,
                         num_cmpnts_to_check=0, cmpnts_indices_to_check=None):
        """
        Check whether a given well has any integrity issues.

        :param p: dictionary of input parameters
        :type p: dict

        :param time_point: time point (in days) for which the monitoring is
            done; by default, its value is 365.25 days
        :type time_point: float

        :param time_index: index of time point defined in the system time_array attribute
        :type time_index: int

        :param cmpnts_indices_to_check: indices of components
            that will be tested at the current time point
        type cmpnts_indices_to_check: array-like

        :returns: dictionary of observations
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        out = {}
        if time_point == 0.0:
            out['num_ident_leak_wells'] = 0
            self.details['cmpnts_indices'][time_index] = None
        else:
            num_ident_leak_wells = 0
            self.details['cmpnts_indices'][time_index] = None

            # Check whether the components to be monitored are linked
            if self.max_num_cmpnts_to_monitor > 0:
                # If no indices were provided
                if not cmpnts_indices_to_check:
                    if num_cmpnts_to_check > 0:
                        # Generate random indices of components to check
                        num_indices = min(num_cmpnts_to_check,
                                          self.max_num_cmpnts_to_monitor)
                        if self.reproducible:
                            cmpnts_indices = generate_seq(
                                self.max_num_cmpnts_to_monitor,
                                seed=actual_p['seed'] + 10.9876*time_point + 1, # random constant
                                array_to_be_returned=0)[0:num_indices]
                        else:
                            cmpnts_indices = generate_seq(
                                self.max_num_cmpnts_to_monitor,
                                array_to_be_returned=0)[0:num_indices]
                        cmpnts_indices = cmpnts_indices-1
                    else:
                        cmpnts_indices = np.array([])  # no components are checked
                else:
                    # If indices are provided save them
                    cmpnts_indices = np.asarray(cmpnts_indices_to_check)

                # Cycle over all provided indices and check whether the components are leaking
                for ind in cmpnts_indices:
                    nm = self.cmpnts_to_monitor[ind]
                    cmpnt = self._parent.component_models[nm]
                    if cmpnt.run_frequency == 2:
                        num_ident_leak_wells = num_ident_leak_wells + 1
                        cmpnt.details['leak_detected_time'] = time_point
                    else:
                        cmpnt.details['risk_type'] = 'fixed'

                # Save which components were checked
                self.details['cmpnts_indices'][time_index] = cmpnts_indices

            out['num_ident_leak_wells'] = num_ident_leak_wells

        # Return dictionary of outputs
        return out


class MonitoringTool3(ComponentModel):
    """
    Monitoring Tool 3 component checks whether mass of CO2 or brine leaked from
    any of the leaked wellbore components exceeds the predefined thresholds.
    If that is the case the leak is considered to be detected.

    Monitoring tool component should precede the wellbore components or any other
    component it monitors in the setup of the system model.
    In the case the component is created after the to be monitored components
    one of the system model methods (reorder_component_models, swap_components,
    or put_in_front_of_component) should be used to correct the order of all
    components.
    """
    def __init__(self, name, parent, cmpnt_nms=None, reproducible=True):
        """
        Constructor method of MonitoringTool3 class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to; it's the same parent that components to be configured
            are linked to.
        :type parent: SystemModel object

        :param cmpnt_nms: list of components whose parameters/keyword
            arguments/attributes will be checked by the monitoring tool
        :type cmpnt_nms: list()

        :param reproducible: flag variable indicating whether the check of the wells
            should be reproducible, i.e. if seed parameter of the model method
            should be used to generate the sequence of random numbers.
        :type reproducible: boolean

        :param CO2_mass_threshold: threshold value used to compare mass of leaked
            CO2 to classify the leak as being detected
        type CO2_mass_threshold: float

        :param brine_mass_threshold: threshold value used to compare mass of leaked
            brine to classify the leak as being detected
        type brine_mass_threshold: float

        :returns: MonitoringTool3 class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25,   # default value of 365.25 days
                        'time_index': 0}        # default value of 0, single point

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)

        print(self.model_kwargs)

        # Add type attribute
        self.class_type = 'MonitoringTool3'

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
        self.add_default_par('seed', value=1)
        self.add_default_par('CO2_mass_threshold', value=1.0e+3)
        self.add_default_par('brine_mass_threshold', value=1.0e+3)
        self.add_default_par('CO2_rate_threshold', value=1.0e-6)
        self.add_default_par('brine_rate_threshold', value=1.0e-6)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['CO2_mass_threshold'] = [1.0e+2, 1.0e+7]
        self.pars_bounds['brine_mass_threshold'] = [1.0e+2, 1.0e+7]
        self.pars_bounds['CO2_rate_threshold'] = [1.0e-7, 1.0e+3]
        self.pars_bounds['brine_rate_threshold'] = [1.0e-7, 1.0e+3]

        # Define attribute that will save the history of what components were checked
        self.details['cmpnts_indices'] = {}

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
                'No components are linked to the monitoring ',
                'component {}'.format(self.name)])
            logging.warning(warn_msg)

    def simulation_model(self, p, time_point=365.25, time_index=0,
                         num_cmpnts_to_check=0, cmpnts_indices_to_remediate=None,
                         cmpnts_indices_to_check=None,
                         # CO2_rate=0.0, brine_rate=0.0,
                         CO2_mass=0.0, brine_mass=0.0):
        """
        Check whether a given well has any integrity issues.

        :param p: dictionary of input parameters
        :type p: dict

        :param time_point: time point (in days) for which the monitoring is
            done; by default, its value is 365.25 days
        :type time_point: float

        :param time_index: index of time point defined in the system time_array attribute
        :type time_index: int

        :param num_cmpnts_to_check: number of components that will be checked
        :type num_cmpnts_to_check: int

        :param cmpnts_indices_to_remediate: indices of components that will be
            remediated at the current time point
        :type cmpnts_indices_to_remediate: array-like

        :param cmpnts_indices_to_check: indices of components
            that will be tested at the current time point
        :type cmpnts_indices_to_check: array-like

        :param CO2_mass: mass of CO2 leaked along the linked wellbore components
        :type CO2_mass: array-like

        :param brine_mass: mass of CO2 leaked along the linked wellbore components
        :type brine_mass: array-like

        :returns: dictionary of observations
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)

        out = {}
        if time_point == 0.0:
            out['num_ident_leak_wells'] = 0
            self.details['cmpnts_indices'][time_index] = None
        else:
            if cmpnts_indices_to_remediate:
                for ind in cmpnts_indices_to_remediate:
                    nm = self.cmpnts_to_monitor[ind]
                    cmpnt = self._parent.component_models[nm]
                    cmpnt.details['leak_detected_time'] = time_point

            num_ident_leak_wells = 0
            self.details['cmpnts_indices'][time_index] = None

            # Check whether the components to be monitored are linked
            if self.max_num_cmpnts_to_monitor > 0:
                # If no indices were provided
                if not cmpnts_indices_to_check:  # None or empty list
                    if num_cmpnts_to_check > 0:
                        # Generate random indices of components to check
                        num_indices = min(num_cmpnts_to_check,
                                          self.max_num_cmpnts_to_monitor)
                        if self.reproducible:
                            cmpnts_indices = generate_seq(
                                self.max_num_cmpnts_to_monitor,
                                actual_p['seed'] + 9.8765*time_point + 1, # random constant
                                array_to_be_returned=0)[0:num_indices]
                        else:
                            cmpnts_indices = generate_seq(
                                self.max_num_cmpnts_to_monitor,
                                array_to_be_returned=0)[0:num_indices]
                        cmpnts_indices = cmpnts_indices-1
                    else:
                        cmpnts_indices = np.array([])  # no components are checked
                else:
                    # If indices are provided save them
                    cmpnts_indices = np.asarray(cmpnts_indices_to_check)

                # Cycle over all provided indices and check whether the components are leaking
                for ind in cmpnts_indices:
                    nm = self.cmpnts_to_monitor[ind]
                    cmpnt = self._parent.component_models[nm]
#                    if (CO2_rate[ind] >= actual_p['CO2_rate_threshold']) or (
#                        brine_rate[ind] >= actual_p['brine_rate_threshold']):
                    if (CO2_mass[ind] >= actual_p['CO2_mass_threshold']) or (
                            brine_mass[ind] >= actual_p['brine_mass_threshold']):
                        num_ident_leak_wells = num_ident_leak_wells + 1
                        cmpnt.details['leak_detected_time'] = time_point

                # Save which components were checked
                self.details['cmpnts_indices'][time_index] = cmpnts_indices

            out['num_ident_leak_wells'] = num_ident_leak_wells

        # Return dictionary of outputs
        return out



def reservoir_data():
    pressure_data = np.array([[27175000., 29385000., 31743000.,
                               32337000., 32931000., 33525000.,
                               33675200., 33825400., 33975600.,
                               34125800., 34276000.],
                              [27048522.44574737, 29118406.7102302, 31037215.67475864,
                               31583376.16990262, 32129536.6650466, 32675697.16019057,
                               32842093.84428151, 33008490.52837246, 33174887.21246339,
                               33341283.89655434, 33507680.58064527],
                              [27309256.65893723, 30208998.91070913, 32124966.74351997,
                               32648529.89465057, 33172093.04578118, 33695656.19691178,
                               33829912.61599338, 33964169.03507498, 34098425.45415658,
                               34232681.87323818, 34366938.29231977],
                              [27236298.06957644, 29945784.11194711, 31993456.91025444,
                               32551093.53539488, 33108730.16053533, 33666366.78567577,
                               33811835.18328377, 33957303.58089177, 34102771.97849976,
                               34248240.37610775, 34393708.77371575]])


    saturation_data = np.array([[0., 0.097242, 0.36106, 0.37590667,
                                 0.39075333, 0.4056, 0.411114, 0.416628,
                                 0.422142, 0.427656, 0.43317],
                                [0., 0.07658825, 0.28437252, 0.31036439,
                                 0.33635627, 0.36234815, 0.37052551, 0.37870287,
                                 0.38688023, 0.39505759, 0.40323496],
                                [0., 0.13858407, 0.35395587, 0.36760335,
                                 0.38125082, 0.3948983, 0.40007711, 0.40525593,
                                 0.41043475, 0.41561356, 0.42079238],
                                [0., 0.12657617, 0.35644597, 0.37346205,
                                 0.39047814, 0.40749422, 0.41309076, 0.4186873,
                                 0.42428384, 0.42988038, 0.43547692]])

    return pressure_data, saturation_data


def test_case_monitoring_tool1():
    from openiam import SystemModel, MultisegmentedWellbore

    num_years = 10
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    pressure_data, saturation_data = reservoir_data()

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Create randomly located leaky well locations
    # within box defined by xmin,xmax,ymin,ymax
    well_xys = np.array([[37538.96951585, 48450.91705453],
                         [37331.27665035, 48297.67494295],
                         [37375.00857811, 48356.92539461],
                         [37461.00669181, 48254.02433701]])
    num_wells = len(well_xys)

    # Add monitoring tool component: the wellbore components are not yet linked
    # but the component is created
    mss_nms = ['ms{}'.format(ind+1) for ind in range(0, num_wells)]
    mntr = sm.add_component_model_object(
        MonitoringTool1(name='mntr', parent=sm, locX=37540.0, locY=48300.0,
                        cmpnt_nms=mss_nms))

    # Initialize list of wellbore components
    mss = []
    for ind, crds in enumerate(well_xys):
        # Add multisegmented wellbore component
        mss.append(sm.add_component_model_object(
            MultisegmentedWellbore(name=mss_nms[ind], parent=sm)))

        mss[-1].add_par('numberOfShaleLayers', value=4, vary=False)
        mss[-1].add_par('shale1Thickness', value=250., vary=False)
        mss[-1].add_par('shale2Thickness', value=550., vary=False)
        mss[-1].add_par('shale3Thickness', value=400., vary=False)
        mss[-1].add_par('shale4Thickness', value=400., vary=False)
        mss[-1].add_par('aquifer1Thickness', value=150., vary=False)
        mss[-1].add_par('aquifer2Thickness', value=720., vary=False)
        mss[-1].add_par('aquifer3Thickness', value=400., vary=False)
        mss[-1].add_par('reservoirThickness', value=400., vary=False)
        mss[-1].add_par('logWellPerm', value=-13.0, vary=False)

        # Add keyword arguments linked to the output provided by reservoir model
        mss[-1].add_dynamic_kwarg('pressure', pressure_data[ind, :])
        mss[-1].add_dynamic_kwarg('CO2saturation', saturation_data[ind, :])

        mss[-1].add_obs_to_be_linked('mass_CO2_aquifer1')
        mss[-1].details['location'] =  well_xys[ind, :]

    # Create collection of leaked masses
    mass_collection = [mss[ind].linkobs['mass_CO2_aquifer1'] for ind in range(num_wells)]

    # Connect collection to the input of monitoring tool
    mntr.add_kwarg_linked_to_collection('CO2_mass', mass_collection)

    # Add monitoring tool outputs
    for obs_nm in ['num_ident_leak_wells', 'num_leak_wells']:
        mntr.add_obs(obs_nm)

    # Run forward model
    sm.forward()

    # Collect outputs
    outputs = {}
    for obs_nm in ['num_ident_leak_wells', 'num_leak_wells']:
        outputs[obs_nm] = sm.collect_observations_as_time_series(mntr, obs_nm)

    # Print outputs
    for obs_nm in ['num_ident_leak_wells', 'num_leak_wells']:
        print(obs_nm, '\n', outputs[obs_nm])


def test_case_monitoring_tool2():
    from openiam import SystemModel, MultisegmentedWellbore

    num_years = 10
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    pressure_data, saturation_data = reservoir_data()

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Create randomly located leaky well locations
    # within box defined by xmin,xmax,ymin,ymax
    well_xys = np.array([[37538.96951585, 48450.91705453],
                         [37331.27665035, 48297.67494295],
                         [37375.00857811, 48356.92539461],
                         [37461.00669181, 48254.02433701]])
    num_wells = len(well_xys)

    # Add monitoring tool component: the wellbore components are not yet linked
    # but the component is created
    mss_nms = ['ms{}'.format(ind+1) for ind in range(0, num_wells)]
    mntr = sm.add_component_model_object(
        MonitoringTool2(name='mntr', parent=sm, cmpnt_nms=mss_nms))

    # Initialize list of wellbore components
    mss = []
    for ind, crds in enumerate(well_xys):
        # Add multisegmented wellbore component
        mss.append(sm.add_component_model_object(
            MultisegmentedWellbore(name=mss_nms[ind], parent=sm)))

        mss[-1].add_par('numberOfShaleLayers', value=4, vary=False)
        mss[-1].add_par('shale1Thickness', value=250., vary=False)
        mss[-1].add_par('shale2Thickness', value=550., vary=False)
        mss[-1].add_par('shale3Thickness', value=400., vary=False)
        mss[-1].add_par('shale4Thickness', value=400., vary=False)
        mss[-1].add_par('aquifer1Thickness', value=150., vary=False)
        mss[-1].add_par('aquifer2Thickness', value=720., vary=False)
        mss[-1].add_par('aquifer3Thickness', value=400., vary=False)
        mss[-1].add_par('reservoirThickness', value=400., vary=False)
        mss[-1].add_par('logWellPerm', value=-13.0, vary=False)

        # Add keyword arguments linked to the output provided by reservoir model
        mss[-1].add_dynamic_kwarg('pressure', pressure_data[ind, :])
        mss[-1].add_dynamic_kwarg('CO2saturation', saturation_data[ind, :])

        mss[-1].details['location'] =  well_xys[ind, :]
        mss[-1].details['leak_detected_time'] = None
        mss[-1].details['risk_type'] = 'exist_path'

    # Add monitoring tool dynamic inputs
    mntr.add_dynamic_kwarg('num_cmpnts_to_check', [0]+10*[1])

    # Add monitoring tool outputs
    mntr.add_obs('num_ident_leak_wells')

    # Component details before the simulation is run
    for ind in range(num_wells):
        print('Component', mss[ind].name, 'details:')
        print(mss[ind].details)
    # Run forward model
    sm.forward()

    # Collect outputs
    outputs = {}
    for obs_nm in ['num_ident_leak_wells']:
        outputs[obs_nm] = sm.collect_observations_as_time_series(
            mntr, obs_nm)

    # Print outputs
    for obs_nm in ['num_ident_leak_wells']:
        print(obs_nm, '\n', outputs[obs_nm])

    # Component details after the simulation is run
    for ind in range(num_wells):
        print('Component', mss[ind].name, 'details:')
        print(mss[ind].details)

    print('Components indices that were checked at each time point')
    print(mntr.details['cmpnts_indices'])

def test_case_monitoring_tool3():
    from openiam import SystemModel, MultisegmentedWellbore

    num_years = 10
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    pressure_data, saturation_data = reservoir_data()

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Create randomly located leaky well locations
    # within box defined by xmin,xmax,ymin,ymax
    well_xys = np.array([[37538.96951585, 48450.91705453],
                         [37331.27665035, 48297.67494295],
                         [37375.00857811, 48356.92539461],
                         [37461.00669181, 48254.02433701]])
    num_wells = len(well_xys)

    # Add monitoring tool component: the wellbore components are not yet linked
    # but the component is created
    mss_nms = ['ms{}'.format(ind+1) for ind in range(0, num_wells)]
    mntr = sm.add_component_model_object(
        MonitoringTool3(name='mntr', parent=sm, cmpnt_nms=mss_nms))
    # Add parameters of monitoring tool component
    mntr.add_par('CO2_mass_threshold', value=1.0e+3, vary=False)
    mntr.add_par('brine_mass_threshold', value=1.0e+3, vary=False)

    # Initialize list of wellbore components
    mss = []
    for ind, crds in enumerate(well_xys):
        # Add multisegmented wellbore component
        mss.append(sm.add_component_model_object(
            MultisegmentedWellbore(name=mss_nms[ind], parent=sm)))

        mss[-1].add_par('numberOfShaleLayers', value=4, vary=False)
        mss[-1].add_par('shale1Thickness', value=250., vary=False)
        mss[-1].add_par('shale2Thickness', value=550., vary=False)
        mss[-1].add_par('shale3Thickness', value=400., vary=False)
        mss[-1].add_par('shale4Thickness', value=400., vary=False)
        mss[-1].add_par('aquifer1Thickness', value=150., vary=False)
        mss[-1].add_par('aquifer2Thickness', value=720., vary=False)
        mss[-1].add_par('aquifer3Thickness', value=400., vary=False)
        mss[-1].add_par('reservoirThickness', value=400., vary=False)
        mss[-1].add_par('logWellPerm', value=-13.0, vary=False)

        # Add keyword arguments linked to the output provided by reservoir model
        mss[-1].add_dynamic_kwarg('pressure', pressure_data[ind, :])
        mss[-1].add_dynamic_kwarg('CO2saturation', saturation_data[ind, :])

        mss[-1].add_obs_to_be_linked('mass_CO2_aquifer1')
        mss[-1].details['location'] =  well_xys[ind, :]
        mss[-1].details['leak_detected_time'] = None
        mss[-1].details['risk_type'] = 'exist_path'

    # Create collection of leaked masses
    mass_collection = [mss[ind].linkobs['mass_CO2_aquifer1'] for ind in range(num_wells)]

    # Connect collection to the input of monitoring tool
    # Mass of brine won't be analized
    mntr.add_kwarg_linked_to_collection('CO2_mass', mass_collection)
    mntr.model_kwargs['brine_mass'] = 4*[0]

    # Add monitoring tool dynamic inputs
    mntr.add_dynamic_kwarg('num_cmpnts_to_check', [0]+10*[2])

    # Add monitoring tool outputs
    for obs_nm in ['num_ident_leak_wells']:
        mntr.add_obs(obs_nm)

    # Run forward model
    sm.forward()

    # Collect outputs
    outputs = {}
    for obs_nm in ['num_ident_leak_wells']:
        outputs[obs_nm] = sm.collect_observations_as_time_series(mntr, obs_nm)

    # Print outputs
    for obs_nm in ['num_ident_leak_wells']:
        print(obs_nm, '\n', outputs[obs_nm])

    print('Components indices that were checked at each time point')
    print(mntr.details['cmpnts_indices'])


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)

    test_case = 3

    test_to_run = {1: test_case_monitoring_tool1,
                   2: test_case_monitoring_tool2,
                   3: test_case_monitoring_tool3}

    # Run corresponding test example
    test_to_run[test_case]()
