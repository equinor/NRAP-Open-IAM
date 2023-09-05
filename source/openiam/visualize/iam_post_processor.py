# -*- coding: utf-8 -*-
"""
Post processing analysis and plots for GUI in NRAP-Open-IAM
"""
import collections
import itertools
import os
import sys
from re import split
import pickle
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import openiam

try:
    from .time_series import time_series_plot
    from .atmosphere_plot import map_plume_plot_single, map_plume_plot_ensemble
    from .sensitivity_analysis import (
        correlations_at_time, simple_sensitivities_barplot,
        multi_sensitivities_barplot, time_series_sensitivities)
except ModuleNotFoundError:
    from time_series import time_series_plot
    from atmosphere_plot import map_plume_plot_single, map_plume_plot_ensemble
    from sensitivity_analysis import (
        correlations_at_time, simple_sensitivities_barplot,
        multi_sensitivities_barplot, time_series_sensitivities)


class IAM_Post_Processor():
    """
    Post processor to build plots of OpenIAM simulation results and
    compute sensitivity analysis and correlation coeffiecents.

    Attributes:
    analysis_type - 'forward', 'lhs', or 'parstudy' openiam simulation type
    outputs - list of outputs valid for plotting
    atm_cmp - boolean if atmospheric components exist for atm map plots
    """
    def __init__(self, folder):
        """
        Constructor method for Post processor class

        :param folder: Path to output folder for OpenIAM simulation results
        :type folder: str

        :returns: IAM_Post_Processor class object
        """
        self.folder = folder
        self.load_folder(folder)
        self.get_outputs()
        self.analysis_type = self.sim_data['ModelParams']['Analysis_type']
        self.components = self.sim_data['components']
        self.time_array = self.sim_data['time_array']

        # Check the presence of components with either special type of outputs
        # or special naming conventions of outputs
        self.spec_comps = {
            'AtmosphericROM': bool(self.find_atm_comp()),
            'GenericAquifer': bool(self.find_ga_comp()),
            'SealHorizon': bool(self.find_sh_comp()),
            'FaultFlow': bool(self.find_ff_comp())}

        # Filter outputs if AtmosphericROM component is present in the system model
        if self.spec_comps['AtmosphericROM']:
            self.filter_outputs()
        else:
            self.filtered_outputs = self.outputs.copy()

    def filter_outputs(self):
        """ Method to separate original form of the outputs from the modified
        observations related to the AtmROM component.
        """
        self.filtered_outputs = []
        self.atm_cmp_outputs = {}
        for key in self.outputs:
            test_res1 = split("_s[0-9]+$", key)  # e.g., x_new_s###
            test_res2 = split("_r[0-9]+$", key)  # e.g., outflag_r###

            if len(test_res1) == 2:  # for x_news###, etc.
                orig_obs_nm = test_res1[0]
            else:
                orig_obs_nm = test_res2[0]  # for outflag_r### and other observations

            if orig_obs_nm not in self.filtered_outputs:
                self.filtered_outputs.append(orig_obs_nm)

            if orig_obs_nm in ['x_new', 'y_new', 'critical_distance', 'outflag']:
                if orig_obs_nm not in self.atm_cmp_outputs:
                    self.atm_cmp_outputs[orig_obs_nm] = []
                self.atm_cmp_outputs[orig_obs_nm].append(key)

    def load_folder(self, folder):
        """
        Method to load OpenIAM results into post processor object
        """
        output_filename1 = os.path.join(folder, 'combined_output.pkl')
        # For compatibility with previously saved results in file named output_dump.pkl
        output_filename2 = os.path.join(folder, 'output_dump.pkl')
        if not os.path.exists(output_filename1):
            if not os.path.exists(output_filename2):
                raise ValueError('IAM Output file not found in specified folder')
            else:
                with open(output_filename2, 'rb') as fi:
                    self.sim_data = pickle.load(fi)
        else:
            with open(output_filename1, 'rb') as fi:
                self.sim_data = pickle.load(fi)

    def get_outputs(self):
        """
        Populate self.outputs attribute with a list of outputs from
        simulation that can be plotted.
        """
        output_list = self.sim_data['output_list']
        self.outputs = []
        for ov in itertools.chain(*output_list.values()):
            if ov not in self.outputs:
                self.outputs.append(ov)

    def plotter(self, outputs, plottype, name=None, title=None,
                savefig=None, subplot=None, sample_num=0):
        """
        Method for plotting simulation results

        :param outputs: List of output names to plot
        :type outputs: list

        :param plottype: String of plottype. Values can be: TimeSeries,
            TimeSeriesStats, TimeSeriesAndStats, AtmPlumeSingle,
            or AtmPlumeEnsemble
        :type plottype: str

        :param name: Name of the figure to create
        :type name: str

        :param title: Title to print on the plot
        :type title: str

        :param savefig: Filename to save the file as. Name only, saved to output folder.
        :type savefig: str

        :param subplot: Dictionary for subplot controls,
            use=True will use subplotting (boolean default=False),
            ncols=n will use n columns of subplots (positive integer default 1)
            comp.obs_name = sub_title will title the subplots of
            comp.obs_name subtitle (string default=comp.var_name)
        :type subplot: dict
        """
        if subplot is None:
            subplot = {'use': False}
        if not name:
            name = ' '.join(outputs)

        savefig_base_name = savefig
        if savefig:
            savefig_base_name = os.path.join(self.folder, savefig)
            if not os.path.exists(os.path.dirname(savefig_base_name)):
                os.mkdir(os.path.dirname(savefig_base_name))

        ts_plot_options = {'TimeSeries': ['real'],
                           'TimeSeriesStats': ['stats'],
                           'TimeSeriesAndStats': ['real', 'stats']}

        if plottype in ['TimeSeries', 'TimeSeriesStats', 'TimeSeriesAndStats']:
            plot_data = {plottype: None,
                         'UseMarkers': False,       # default is False
                         'UseLines': True,          # default is True
                         'VaryLineStyles': False,   # default is False
                         'FigureDPI': 100,
                         'subplot': subplot}
            if not subplot['use']:
                for obs_nm in outputs:
                    res_components = self.resolve_components(obs_nm)
                    for comp in res_components:
                        filename = None
                        if savefig:
                            filename = savefig_base_name.format(
                                ob='.'.join([comp.name, obs_nm]))
                        # Setup figure name
                        name = '_'.join([comp.name, obs_nm])
                        # Update plot_data dictionary
                        plot_data[plottype] = [obs_nm]
                        # Create plots
                        time_series_plot(
                            [obs_nm], self.sim_data['sm'], self.sim_data['s'],
                            plot_data, {comp: res_components[comp]}, name=name,
                            analysis=self.analysis_type, savefig=filename,
                            title=title, subplot=subplot,
                            plot_type=ts_plot_options[plottype])
            else:
                # Update plot_data dictionary
                plot_data[plottype] = outputs
                # Create plots
                time_series_plot(
                    outputs, self.sim_data['sm'], self.sim_data['s'], plot_data,
                    self.sim_data['output_list'], name=name,
                    analysis=self.analysis_type,
                    savefig=savefig_base_name.format(ob='_'.join(outputs)),
                    title=title, subplot=subplot,
                    plot_type=ts_plot_options[plottype])
            plt.show()

        else:
            if plottype == 'AtmPlumeSingle':
                plot_data = {'AtmPlumeSingle': {
                    'Realization': sample_num,
                    'FigureDPI': 100,
                    'PlotInjectionSites': False,
                    'PlotReceptors': False,
                    'SaveCSVFiles': False}}
                satm = self.find_atm_comp()
                map_plume_plot_single(
                    plot_data, 'atm_plume_single', self.sim_data['sm'],
                    self.sim_data['s'], satm, self.time_array,
                    self.folder, analysis=self.analysis_type,
                    savefig=savefig_base_name, title=title)

            elif plottype == 'AtmPlumeEnsemble':
                plot_data = {'AtmPlumeEnsemble': {
                    'FigureDPI': 100,
                    'PlotInjectionSites': False,
                    'PlotReceptors': False,
                    'SaveCSVFiles': False}}
                satm = self.find_atm_comp()
                map_plume_plot_ensemble(
                    plot_data, 'atm_plume_ensemble', self.sim_data['sm'],
                    self.sim_data['s'], satm, self.time_array,
                    self.folder, analysis=self.analysis_type,
                    savefig=savefig_base_name, title=title)

            return

    def find_atm_comp(self):
        """
        Method to find AtmosphericROM component of a system model
        """
        for comp in self.components[::-1]:
            if isinstance(comp, openiam.AtmosphericROM):
                return comp

        return False

    def find_ga_comp(self):
        """
        Method to find Generic Aquifer component(s) of a system model.
        """
        comp_list = []
        for comp in self.components[::-1]:
            if isinstance(comp, openiam.GenericAquifer):
                comp_list.append(comp)

        if comp_list:
            if len(comp_list) == 1:
                return comp_list[0]
            else:
                return comp_list
        else:
            return False

    def find_sh_comp(self):
        """
        Method to find Seal Horizon component(s) of a system model.
        """
        comp_list = []
        for comp in self.components[::-1]:
            if isinstance(comp, openiam.SealHorizon):
                comp_list.append(comp)

        if comp_list:
            if len(comp_list) == 1:
                return comp_list[0]
            else:
                return comp_list
        else:
            return False

    def find_ff_comp(self):
        """
        Method to find Fault Flow component(s) of a system model.
        """
        comp_list = []
        for comp in self.components[::-1]:
            if isinstance(comp, openiam.FaultFlow):
                comp_list.append(comp)

        if comp_list:
            if len(comp_list) == 1:
                return comp_list[0]
            else:
                return comp_list
        else:
            return False

    def sensitivity_analysis(self, sensitivity_type, outputs, sensitivity_dict):
        """
        Method for computing and plotting correlation or sensitivity coeffiecents

        :param sensitivity_type: Type of Sensitivity Analysis to perform.
            Valid inputs are CorrelationCoeff, SensitivityCoeff,
            MultiSensitivities, TimeSeriesSensitivity
        :type sensitivity_type: str

        :param outputs: List of output names to perform sensitivity analysis on.
        :type outputs: list

        :param sensitivity_dict: Dictionary for sensitivity analysis options
        :type sensitivity_dict: dict
        """
        if sensitivity_type == 'CorrelationCoeff':
            if outputs:
                outputs_to_exclude = self.resolve_obs_names(outputs)
            else:
                outputs_to_exclude = []
            corrcoeff_dict = {'capture_point': len(self.time_array)-1,
                              'excludes': outputs_to_exclude,
                              'ctype': 'pearson',
                              'plot': True,
                              'printout': False,
                              'plotvals': True,
                              'figsize': (12, 12),
                              'title': 'Pearson Correlation Coefficients at {ct} years',
                              'xrotation': 90,
                              'savefig': 'correlation_coefficients_time_index_{cp}.png',
                              'outfile': 'correlation_coefficients_time_index_{cp}.csv',
                              'GUI_output': True}
            corrcoeff_dict.update(sensitivity_dict)
            if corrcoeff_dict['savefig']:
                if corrcoeff_dict['savefig'][-4:] != '.png':
                    corrcoeff_dict['savefig'] = corrcoeff_dict['savefig'] + '.png'
                corrcoeff_dict['savefig'] = os.path.join(
                    self.folder, corrcoeff_dict['savefig'])
            if corrcoeff_dict['outfile']:
                corrcoeff_dict['outfile'] = os.path.join(
                    self.folder, corrcoeff_dict['outfile'])
            correlations_at_time(self.sim_data['s'], self.time_array, **corrcoeff_dict)

        elif sensitivity_type == 'SensitivityCoeff':
            if isinstance(outputs, str):
                outputs = [outputs]
            res_obs = self.resolve_obs_names(outputs)
            if 'capture_point' in sensitivity_dict:
                capture_point = sensitivity_dict.pop('capture_point')
            else:
                capture_point = len(self.time_array)-1
            if not isinstance(capture_point, collections.Iterable):
                capture_points = [capture_point]
            else:
                capture_points = capture_point
            for capture_point in capture_points:
                cp_obs = ['{ob}_{cp}'.format(ob=ob, cp=capture_point)
                          for ob in res_obs]
                sens_dict_full = {'title': '{ob} Sensitivity Coefficients',
                                  'ylabel': None,
                                  'savefig': '{ob}_sensitivity.png',
                                  'outfile': '{ob}_sensitivity.csv',
                                  'name': '{ob}_{cp}',
                                  'GUI_output': True}
                sens_dict_full.update(sensitivity_dict)
                for ob in cp_obs:
                    if sens_dict_full['savefig']:
                        if sens_dict_full['savefig'][-4:] != '.png':
                            sens_dict_full['savefig'] = sens_dict_full['savefig'] + '.png'
                        sens_dict_full['savefig'] = os.path.join(
                            self.folder, sens_dict_full['savefig'])
                    if sens_dict_full['outfile']:
                        sens_dict_full['outfile'] = os.path.join(
                            self.folder, sens_dict_full['outfile'])
                    sens_dict2 = {}
                    for key, value in list(sens_dict_full.items()):
                        if isinstance(value, str):
                            v = value.format(ob=ob, cp=capture_point,
                                             ct=self.time_array[capture_point]/365.25)
                        else:
                            v = value
                        sens_dict2[key] = v
                    sensitivities = self.sim_data['s'].rbd_fast(
                        obsname=ob, print_to_console=False)
                    simple_sensitivities_barplot(
                        sensitivities, self.sim_data['sm'], **sens_dict2)

        elif sensitivity_type == 'MultiSensitivities':
            if isinstance(outputs, str):
                outputs = [outputs]
            res_obs = self.resolve_obs_names(outputs)
            if 'capture_point' in sensitivity_dict:
                capture_point = sensitivity_dict.pop('capture_point')
            else:
                capture_point = len(self.time_array)-1
            if not isinstance(capture_point, collections.Iterable):
                capture_points = [capture_point]
            else:
                capture_points = capture_point
            for capture_point in capture_points:
                cp_obs = ['{ob}_{cp}'.format(ob=ob, cp=capture_point)
                          for ob in res_obs]
                sens_dict_full = {
                    'title': 'Sensitivity Coefficients at {ct} years'.format(
                        ct=self.time_array[capture_point]/365.25),
                    'ylabel': None,
                    'savefig': 'multisensitivities_time_{ct}.png'.format(
                        ct=self.time_array[capture_point]/365.25),
                    'outfile': 'multisensitivities_time_{ct}.csv'.format(
                        ct=self.time_array[capture_point]/365.25),
                    'GUI_output': True}
                sens_dict_full.update(sensitivity_dict)
                if sens_dict_full['savefig']:
                    if sens_dict_full['savefig'][-4:] != '.png':
                        sens_dict_full['savefig'] = sens_dict_full['savefig'] + '.png'
                    sens_dict_full['savefig'] = os.path.join(
                        self.folder, sens_dict_full['savefig'])
                if sens_dict_full['outfile']:
                    sens_dict_full['outfile'] = os.path.join(
                        self.folder, sens_dict_full['outfile'])
                multi_sensitivities_barplot(cp_obs,
                                            self.sim_data['sm'],
                                            self.sim_data['s'],
                                            **sens_dict_full)

        elif sensitivity_type == 'TimeSeriesSensitivity':
            if isinstance(outputs, str):
                outputs = [outputs]
            res_obs = self.resolve_obs_names(outputs)
            if 'capture_point' in sensitivity_dict:
                capture_point = sensitivity_dict.pop('capture_point')
            else:
                capture_point = len(self.time_array)-1
            sens_dict_full = {'title': 'Time Series Sensitivity Coefficients for {ob}',
                              'ylabel': None,
                              'num_sensitivities': 5,
                              'savefig': '{ob}_time_series_sensitivity.png',
                              'outfile': '{ob}_time_series_sensitivity.csv',
                              'GUI_output': True}
            sens_dict_full.update(sensitivity_dict)
            for ob in res_obs:
                if sens_dict_full['savefig']:
                    if sens_dict_full['savefig'][-4:] != '.png':
                        sens_dict_full['savefig'] = sens_dict_full['savefig'] + '.png'
                    sens_dict_full['savefig'] = os.path.join(
                        self.folder, sens_dict_full['savefig'])
                if sens_dict_full['outfile']:
                    sens_dict_full['outfile'] = os.path.join(
                        self.folder, sens_dict_full['outfile'])
                sens_dict = {}
                for key, value in list(sens_dict_full.items()):
                    if isinstance(value, str):
                        v = value.format(ob=ob)
                    else:
                        v = value
                    sens_dict[key] = v

                time_series_sensitivities(ob, self.sim_data['sm'],
                                          self.sim_data['s'],
                                          self.time_array,
                                          capture_point=capture_point,
                                          **sens_dict)
        else:
            raise IOError('sensitivity_type not recognized')

        plt.show()

    def resolve_components(self, obs):
        """
        Returns dictionary with keys being the components which have observations
        with name obs.
        """
        output_list = self.sim_data['output_list']
        resolved_components = {}
        for output_component in list(output_list.keys()):
            if obs in output_list[output_component]:
                resolved_components[output_component] = output_list[output_component]

        return resolved_components


    def resolve_obs_names(self, obs):
        '''
        Takes in base observation name and returns list of appended cm.obs for
        all component model names that have observations matching the base name.
        '''
        output_list = self.sim_data['output_list']
        resolved_names = []
        for ob in obs:
            for output_component in list(output_list.keys()):
                if ob in output_list[output_component]:
                    resolved_names.append('.'.join([output_component.name, ob]))
        return resolved_names
