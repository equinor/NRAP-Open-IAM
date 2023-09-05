# -*- coding: utf-8 -*-
"""
Code to create map-view figures showing the evolution of contaminant plumes and
the estimated time to first detection (TTFD) for a given monitoring network.

Examples illustrating applications of ttfd_plot method:
    iam_sys_reservoir_mswell_futuregen_ttfdplot_dipping_strata.py
    iam_sys_reservoir_mswell_futuregen_ttfdplot_no_dip.py
    iam_sys_reservoir_mswell_futuregen_ttfdplot_no_dip_lhs.py
    ControlFile_ex39a.yaml
    ControlFile_ex39b.yaml
    ControlFile_ex40.yaml
    ControlFile_ex41.yaml
    ControlFile_ex42.yaml
    ControlFile_ex43.yaml

Last Modified: July, 2023

@author: Nate Mitchell (Nathaniel.Mitchell@NETL.DOE.GOV)
@contributor: Veronika Vasylkivska (Veronika.Vasylkivska@NETL.DOE.GOV)
LRST (Battelle|Leidos) supporting NETL
"""
import os
import sys
import logging
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import ticker


SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(SOURCE_DIR)

import openiam.cfi.commons as iamcommons
import openiam.cfi.strata as strata

TTFD_RESERVOIR_COMPONENTS = ['LookupTableReservoir',
                             'SimpleReservoir',
                             'AnalyticalReservoir']

RC_FONT = {'family': 'Arial', 'weight': 'normal', 'size': None}

TTFD_AQUIFER_COMPONENTS = ['CarbonateAquifer', 'FutureGen2Aquifer',
                           'FutureGen2AZMI', 'GenericAquifer',
                           'DeepAlluviumAquifer', 'DeepAlluviumAquiferML']

TITLE_OPTIONS = {
    'plumeTimings': {
        'Pressure': 'Spatiotemporal Evolution of Pressure Plumes{}',
        'pH': 'Spatiotemporal Evolution of pH Plumes{}',
        'TDS': 'Spatiotemporal Evolution of TDS Plumes{}',
        'Dissolved_CO2': 'Spatiotemporal Evolution of Dissolved CO$_2$ Plumes{}',
        'Dissolved_salt': 'Spatiotemporal Evolution of Dissolved Salt Plumes{}',
        'CarbonateAquifer':
            'Spatiotemporal Evolution of Plumes (Carbonate Aquifer Component){}'},
    'monitoringTTFD': {
        'Pressure': 'Time to First Detection (TTFD) for Pressure Plumes{}',
        'pH': 'Time to First Detection (TTFD) for pH Plumes{}',
        'TDS': 'Time to First Detection (TTFD) for TDS Plumes{}',
        'Dissolved_CO2': 'Time to First Detection (TTFD) for Dissolved CO$_2$ Plumes{}',
        'Dissolved_salt': 'Time to First Detection (TTFD) for Dissolved Salt Plumes{}',
        'CarbonateAquifer':
            'Time to First Detection (TTFD) for Plumes (Carbonate Aquifer Component){}'},
    'plumeProbabilities': {
        'Pressure': 'Probabilities of Pressure Plume Occurence{}',
        'pH': 'Probabilities of pH Plume Occurence{}',
        'TDS': 'Probabilities of TDS Plume Occurence{}',
        'Dissolved_CO2': 'Probabilities of Dissolved CO$_2$ Plume Occurence{}',
        'Dissolved_salt': 'Probabilities of Dissolved Salt Plume Occurence{}',
        'CarbonateAquifer':
            'Probabilities of Plume Occurence (Carbonate Aquifer Component){}'}}

AX_TITLE_PARTIAL = 3*[dict()]  # list of 3 dictionaries with the same keys
# The first dictionary in the list
AX_TITLE_PARTIAL[0] = {
    'plumeTimings': 'No Plumes',
    'monitoringTTFD': 'No TTFD From Monitoring Network',
    'plumeProbabilities': 'All Plume Probabilities are 0%'}
# The second dictionary in the list
AX_TITLE_PARTIAL[1] = {
    'plumeTimings': 'Plume Timing Range: ',
    'monitoringTTFD': 'Range in TTFD from Monitoring Network: ',
    'plumeProbabilities': 'Probability Range: '}
# The third dictionary in the list
AX_TITLE_PARTIAL[2] = {
    'plumeTimings': 'Singular Plume Timing: ',
    'monitoringTTFD': 'TTFD from Monitoring Network: ',
    'plumeProbabilities': 'Singular Probability: '}

AX_TITLE_PARTIAL_2 = 3*[dict()]  # list of 3 dictionaries with the same keys
# The first dictionary in the list
AX_TITLE_PARTIAL_2[0] = {
    'plumeTimings': '{}{}{}',
    'monitoringTTFD': '{}{}{}',
    'plumeProbabilities': '{}{}{}'}
# The second dictionary in the list
AX_TITLE_PARTIAL_2[1] = {
    'plumeTimings': '{:.0f} year{} to {:.0f} years',
    'monitoringTTFD': '{:.0f} year{} to {:.0f} years',
    'plumeProbabilities': '{:.2f} %{} to {:.2f} %'}
# The third dictionary in the list
AX_TITLE_PARTIAL_2[2] = {
    'plumeTimings': '{:.0f} year{}{}',
    'monitoringTTFD': '{:.0f} year{}',
    'plumeProbabilities': '{:.2f} %{}{}'}

BACKGROUND_COLOR = [0.67, 0.67, 0.67]

# The min. plume timing results are obtained by making a grid with very high
# values (MAX_TIME), calculating the plume timings, and then taking the minimum
# values across the grid. Cases with valid plume timings are found as
# plumeTiming < THRESHOLD_TIME.
MAX_TIME = 1.0e30
THRESHOLD_TIME = 1.0e29

MIN_PROBABILITY = 1.0e-6

THRESHOLDS = {'pH': 0.2, 'Pressure': 0.00065, 'Temperature': 0.0003,
             'TDS': 0.1, 'Dissolved_CO2': 0.2, 'Dissolved_salt': 0.1}
INDICATOR = {'pH': 'absolute', 'Pressure': 'relative',
             'Temperature': 'relative', 'TDS': 'relative',
             'Dissolved_CO2': 'relative', 'Dissolved_salt': 'relative',
             'CarbonateAquifer': 'absolute'}


def ttfd_plot(yaml_data, model_data, sm, s,
              output_dir, name='TTFD_Figure1', analysis='lhs',
              figsize=(10, 8), genfontsize=12,
              axislabelfontsize=14, titlefontsize=14, boldlabels=True):
    """
    Makes a map-view figure showing the estimated time to first detection (TTFD)
    for a given monitoring network.

    Note that aq_name_list (used below) is a list of the names for Aquifer
    components producing the plume metrics. Note that these names are those
    typically shown in the .yaml files (e.g., GenericAquifer1, GenericAquifer2)
    and not the names used within the control file interface (e.g.,
    GenericAquifer1_000 and GenericAquifer2_000 for GenericAquifer1 and
    GenericAquifer 2 at location 1).

    :param yaml_data: Dictionary of input values
    :type yaml_data: dict

    :param model_data: Input from the 'ModelParams' section of the .yaml file
    :type model_data: dict

    :param locations: dictionary of locations assigned to each applicable component
    :type locations: dict

    :param sm: OpenIAM System model for which plots are created
    :type sm: openiam.SystemModel object

    :param s: SampleSet with results of analysis performed on the system model.
        None for forward analysis.
    :type s: matk.SampleSet object

    :param name: Figure Name to be used/created.
    :type name: str

    :param analysis: Type of OpenIAM system analysis performed ('forward',
        'lhs', or 'parstudy')
    :type analysis: str

    :param savefig: Filename to save resulting figure to. No file saved if None.
    :type savefig: str

    :param figsize: width and height of the figure (width, height), in inches.
        Default value is (10, 8).
    :type figsize: tuple

    :param genfontsize: fontsize for tick labels, etc.
    :type genfontsize: float or int

    :param axislabelfontsize: fontsize for x and y axis labels
    :type axislabelfontsize: float or int

    :param titlefontsize: fontsize for the title
    :type titlefontsize: float or int

    :param boldlabels: option to use bold x and y labels and bold titles. Set to
        True for bold labels, False for normal labels.
    :type boldlabels: bool

    :return: None
    """
    plume_metric_abbrev = yaml_data['Plots'][name]['TTFD']['PlumeType']
    aq_name_list = yaml_data['Plots'][name]['TTFD']['ComponentNameList']

    if not isinstance(aq_name_list, list):
        aq_name_list = [aq_name_list]

    TTFD_yaml_input_dict = get_ttfd_plot_yaml_input(yaml_data, name)

    if TTFD_yaml_input_dict['monitoringTTFD'] is not None:
        monitoringTTFD = TTFD_yaml_input_dict['monitoringTTFD']
        monitoringCoordX = TTFD_yaml_input_dict['monitoringCoordX']
        monitoringCoordY = TTFD_yaml_input_dict['monitoringCoordY']
        monitorHorizontalWindow = TTFD_yaml_input_dict['HorizontalWindow']
        monitorVerticalWindow = TTFD_yaml_input_dict['VerticalWindow']

        if plume_metric_abbrev != 'CarbonateAquifer':
            monitoringCoordZ = TTFD_yaml_input_dict['monitoringCoordZ']
        else:
            monitoringCoordZ = None

    if not monitoringTTFD:
        monitoringCoordX = None
        monitoringCoordY = None
        monitoringCoordZ = None

    numPointsInAq = TTFD_yaml_input_dict['numPointsInAq']
    numPointsInShales = TTFD_yaml_input_dict['numPointsInShales']

    x_grid_spacing = TTFD_yaml_input_dict['x_grid_spacing']
    y_grid_spacing = TTFD_yaml_input_dict['y_grid_spacing']

    EnforceGridXandYLims = TTFD_yaml_input_dict['EnforceGridXandYLims']
    gridXLims = TTFD_yaml_input_dict['gridXLims']
    gridYLims = TTFD_yaml_input_dict['gridYLims']

    save_results = TTFD_yaml_input_dict['saveCSVFiles']
    write_DREAM_output = TTFD_yaml_input_dict['write_DREAM_output']

    if boldlabels:
        selected_labelfontweight = 'bold'
    else:
        selected_labelfontweight = 'normal'

    time_array = sm.time_array

    aq_components, aq_component_types, aq_component_indices, \
        aq_component_xvals, aq_component_yvals, aq_component_ithresh, \
            res_comp_injX, res_comp_injY, x_range, y_range = \
                get_aq_comp_lists_and_xy_grids(sm, yaml_data,
                                               name, aq_name_list)

    x_grid, y_grid = make_xandy_grids(
        x_range, y_range, x_grid_spacing=x_grid_spacing,
        y_grid_spacing=y_grid_spacing, EnforceGridXandYLims=EnforceGridXandYLims,
        gridXLims=gridXLims, gridYLims=gridYLims,
        monitoringCoordX=monitoringCoordX, monitoringCoordY=monitoringCoordY)

    strata_var_info = strata.get_strata_var_info_from_yaml(yaml_data)

    var_type = strata_var_info['var_type']
    strike = strata_var_info['strike']
    dip = strata_var_info['dip']
    dipDirection = strata_var_info['dipDirection']
    coordxReferencePoint = strata_var_info['coordxReferencePoint']
    coordyReferencePoint = strata_var_info['coordyReferencePoint']

    # Get a stratigraphy component
    if var_type == 'noVariation':
        strata_comp = sm.component_models['strata']
    elif var_type == 'strikeAndDip':
        strata_comp = sm.component_models['strataRefPoint']

    strata_dict = strata.get_strata_info_from_component(strata_comp)

    numShaleLayers = strata_dict['numberOfShaleLayers']
    shaleThicknesses = strata_dict['shaleThicknesses']
    aquiferThicknesses = strata_dict['aquiferThicknesses']
    reservoirThickness = strata_dict['reservoirThickness']

    # If unit thicknesses change over space, these are updated within the loop
    aquiferBottomDepths = get_aquifer_bottom_depths(
        numShaleLayers, shaleThicknesses, aquiferThicknesses)

    z_grid = get_z_values(
        numShaleLayers, shaleThicknesses, aquiferThicknesses, reservoirThickness,
        aquiferBottomDepths, x_grid, y_grid, aq_component_types,
        var_type=var_type, strike=strike, dip=dip, dipDirection=dipDirection,
        coordxReferencePoint=coordxReferencePoint,
        coordyReferencePoint=coordyReferencePoint,
        stratigraphy_comp=strata_comp, numPointsInAq=numPointsInAq,
        numPointsInShales=numPointsInShales, monitoringCoordZ=monitoringCoordZ)

    if output_dir and write_DREAM_output:
        write_dream_grid(x_grid, y_grid, z_grid, output_dir, var_type=var_type)

    checkCarbAq = 'CarbonateAquifer' in aq_component_types

    if save_results:
        save_grid_to_csv(x_grid, y_grid, z_grid, output_dir, var_type=var_type,
                         checkCarbAq=checkCarbAq)

    if analysis in ['lhs', 'parstudy']:
        num_samples =  model_data['Analysis']['siz']

        if checkCarbAq:
            plumeProb = np.zeros((len(x_grid), len(x_grid)))
        else:
            plumeProb = np.zeros((z_grid.shape[0], len(y_grid), len(x_grid)))
    else:
        num_samples = 1

    for sample in range(num_samples):
        realization = str(sample)

        plotType = 'plumeTimings'
        plumeTimings = get_plume_timings(
            sm, s, time_array, sample, plume_metric_abbrev,
            x_grid, y_grid, z_grid, aq_components, aq_component_types,
            aq_component_xvals, aq_component_yvals, aq_component_indices,
            numShaleLayers, shaleThicknesses, aquiferThicknesses,
            reservoirThickness, strata_comp, analysis=analysis,
            strike=strike, dip=dip, dipDirection=dipDirection,
            coordxReferencePoint=coordxReferencePoint,
            coordyReferencePoint=coordyReferencePoint, var_type=var_type)

        if save_results:
            save_results_to_csv(
                plumeTimings, x_grid, y_grid, z_grid, output_dir, plume_metric_abbrev,
                plotType, ttfd_list=None, ttfd_x_list=None, ttfd_y_list=None,
                ttfd_z_list=None, analysis=analysis, realization=realization,
                num_samples=num_samples, var_type=var_type,
                checkCarbAq=checkCarbAq)

        if analysis in ['lhs', 'parstudy']:
            plumeProb[plumeTimings < THRESHOLD_TIME] += 1

        if write_DREAM_output:
            write_dream_output(
                plumeTimings, x_grid, y_grid, z_grid, plume_metric_abbrev,
                output_dir, realization=realization, var_type=var_type,
                aq_component_ithresh=aq_component_ithresh, checkCarbAq=checkCarbAq)

        plot_plume_metric(
            plumeTimings, plotType, yaml_data, num_samples, time_array, x_grid,
            y_grid, z_grid, plume_metric_abbrev, aq_component_types,
            aq_component_xvals, aq_component_yvals, ttfd_list=None,
            ttfd_x_list=None, ttfd_y_list=None, name=name,
            analysis=analysis, realization=realization, output_dir=output_dir,
            genfontsize=genfontsize, axislabelfontsize=axislabelfontsize,
            titlefontsize=titlefontsize, labelfontweight=selected_labelfontweight,
            colormap='viridis', figsize=figsize, res_comp_injX=res_comp_injX,
            res_comp_injY=res_comp_injY, var_type=var_type)

        if monitoringTTFD:
            monitoringXIndex, monitoringYIndex = get_xandy_monitoring_indices(
                monitoringCoordX, monitoringCoordY, x_grid, y_grid)

            plotType = 'monitoringTTFD'
            ttfd_list, ttfd_x_list, ttfd_y_list, ttfd_z_list = \
                get_monitoring_location_ttfd(
                    plumeTimings, x_grid, y_grid, z_grid,
                    monitoringCoordX, monitoringCoordY, monitoringCoordZ,
                    monitoringXIndex, monitoringYIndex,
                    monitorHorizontalWindow, monitorVerticalWindow,
                    aq_component_types, var_type=var_type)

            if save_results:
                save_results_to_csv(
                    None, x_grid, y_grid, z_grid, output_dir, plume_metric_abbrev,
                    plotType, ttfd_list=ttfd_list, ttfd_x_list=ttfd_x_list,
                    ttfd_y_list=ttfd_y_list, ttfd_z_list=ttfd_z_list,
                    analysis=analysis, realization=realization,
                    num_samples=num_samples, var_type=var_type,
                    checkCarbAq=checkCarbAq)

            plot_plume_metric(
                plumeTimings, plotType, yaml_data, num_samples, time_array, x_grid,
                y_grid, z_grid, plume_metric_abbrev, aq_component_types,
                aq_component_xvals, aq_component_yvals, ttfd_list=ttfd_list,
                ttfd_x_list=ttfd_x_list, ttfd_y_list=ttfd_y_list,
                name=name, analysis=analysis,
                realization=realization, output_dir=output_dir,
                genfontsize=genfontsize, axislabelfontsize=axislabelfontsize,
                titlefontsize=titlefontsize, labelfontweight=selected_labelfontweight,
                colormap='viridis', figsize=figsize, res_comp_injX=res_comp_injX,
                res_comp_injY=res_comp_injY, var_type=var_type)
    # End of loop through LHS samples (or a singular forward simulation)

    if analysis in ['lhs', 'parstudy']:
        plumeProb /= num_samples
        plumeProb *= 100

        plotType = 'plumeProbabilities'
        if save_results:
            save_results_to_csv(
                plumeProb, x_grid, y_grid, z_grid, output_dir, plume_metric_abbrev,
                plotType, ttfd_list=ttfd_list, ttfd_x_list=ttfd_x_list,
                ttfd_y_list=ttfd_y_list, ttfd_z_list=ttfd_z_list,
                analysis=analysis, realization=realization,
                num_samples=num_samples, var_type=var_type,
                checkCarbAq=checkCarbAq)

        plot_plume_metric(
            plumeProb, plotType, yaml_data, num_samples, time_array, x_grid,
            y_grid, z_grid, plume_metric_abbrev, aq_component_types,
            aq_component_xvals, aq_component_yvals, ttfd_list=None,
            ttfd_x_list=None, ttfd_y_list=None, name=name,
            analysis=analysis, realization=realization, output_dir=output_dir,
            genfontsize=genfontsize, axislabelfontsize=axislabelfontsize,
            titlefontsize=titlefontsize, labelfontweight=selected_labelfontweight,
            colormap='plasma', figsize=figsize, res_comp_injX=res_comp_injX,
            res_comp_injY=res_comp_injY, var_type=var_type)


def get_aquifer_bottom_depths(numShaleLayers, shaleThicknesses, aquiferThicknesses):
    """
    This function is used to create aquifer bottom depths from unit thicknesses
    """
    aquiferBottomDepths = [None] * (numShaleLayers - 1)

    for aquiferRef in range(numShaleLayers - 2, -1, -1):
        aquiferBottomDepths[aquiferRef] = sum(shaleThicknesses[
            aquiferRef + 1:None]) + sum(aquiferThicknesses[aquiferRef:None])

    return aquiferBottomDepths


def get_z_values(numShaleLayers, shaleThicknesses, aquiferThicknesses,
                 reservoirThickness, aquiferBottomDepths, x_grid, y_grid,
                 aq_component_types, var_type='noVariation', strike=None,
                 dip=None, dipDirection=None, coordxReferencePoint=None,
                 coordyReferencePoint=None, stratigraphy_comp=None,
                 numPointsInAq=10, numPointsInShales=3,
                 monitoringCoordZ=None, min_z_spacing=0.01, lowest_depth=0):
    """
    This function is used to create the z values used in the calculation of TTFD.
    Any monitoringCoordZ values provided are included.
    """
    # This is used to keep track of whether each monitoringCoordZ value
    # has been added (1 for unused, 0 for already inserted).
    if isinstance(monitoringCoordZ, list):
        monitoringCoordZCheck = np.ones(len(monitoringCoordZ))
        # For 'strikeAndDip', the ZCheck needs to be reset within the loop
        monitoringCoordZCheckOrig = monitoringCoordZCheck[:]

    elif isinstance(monitoringCoordZ, (int, float)):
        monitoringCoordZCheck = 1
        monitoringCoordZCheckOrig = 1

    if var_type == 'strikeAndDip' and not 'CarbonateAquifer' in aq_component_types:
        z = None

        for xRef, xVal in enumerate(x_grid):
            for yRef, yVal in enumerate(y_grid[:, 0]):
                # Reset monitoringCoordZCheck, if it is being used
                if monitoringCoordZ is not None:
                    monitoringCoordZCheck = monitoringCoordZCheckOrig[:]

                updatedStratigraphy = strata.update_stratigraphy_by_strike_and_dip(
                    numberOfShaleLayers=numShaleLayers,
                    shaleThicknessList=shaleThicknesses[:],
                    aquiferThicknessList=aquiferThicknesses[:],
                    reservoirThickness=reservoirThickness,
                    strike=strike, dip=dip, dipDirection=dipDirection,
                    coordxRefPoint=coordxReferencePoint,
                    coordyRefPoint=coordyReferencePoint,
                    location_x=float(xVal), location_y=float(yVal),
                    strataRefPoint=stratigraphy_comp)

                shaleThicknessesUpdated = updatedStratigraphy['shaleThicknessList']
                aquiferThicknessesUpdated = updatedStratigraphy['aquiferThicknessList']

                aquiferBottomDepths = get_aquifer_bottom_depths(
                    numShaleLayers, shaleThicknessesUpdated,
                    aquiferThicknessesUpdated)

                for shaleRef in range(numShaleLayers - 1):
                    z_aq = np.linspace(-aquiferBottomDepths[shaleRef],
                                       -aquiferBottomDepths[shaleRef]
                                       + aquiferThicknessesUpdated[shaleRef],
                                       numPointsInAq)
                    z_sh = np.linspace(z_aq[-1], z_aq[-1]
                                       + shaleThicknessesUpdated[shaleRef + 1],
                                       numPointsInShales)

                    if shaleRef == 0:
                        z_sh_lowest = np.linspace(
                            z_aq[0] - shaleThicknessesUpdated[shaleRef], z_aq[0],
                            numPointsInShales)

                    # Check if a monitoringCoordZ value fits within the current
                    # aquifer or shale. If so, add it.
                    if monitoringCoordZ is not None:
                        if isinstance(monitoringCoordZ, list):
                            for monitorRef, monitoringCoordZVal in enumerate(monitoringCoordZ):
                                if monitoringCoordZCheck[monitorRef]:
                                    # Check if the point is within this aquifer or shale
                                    if np.min(z_aq) <= monitoringCoordZVal < np.max(z_aq):
                                        # If there is already a point at that
                                        # depth, or close enough, then do
                                        # not add it.
                                        min_diff = np.min(np.abs(
                                            z_aq - monitoringCoordZVal))
                                        if min_diff > min_z_spacing:
                                            z_aq_list = z_aq.tolist()
                                            z_aq_list.append(monitoringCoordZVal)
                                            # This puts it in the correct order
                                            z_aq = np.unique(z_aq_list)
                                        monitoringCoordZCheck[monitorRef] = 0

                                    elif np.min(z_sh) <= monitoringCoordZVal < np.max(z_sh):
                                        min_diff = np.min(np.abs(
                                            z_sh -monitoringCoordZVal))
                                        if min_diff > min_z_spacing:
                                            z_sh_list = z_sh.tolist()
                                            z_sh_list.append(monitoringCoordZVal)
                                            # This puts it in the correct order
                                            z_sh = np.unique(z_sh_list)
                                        monitoringCoordZCheck[monitorRef] = 0

                        else:
                            # Only one monitoring location
                            if monitoringCoordZCheck:
                                if np.min(z_aq) <= monitoringCoordZ < np.max(z_aq):
                                    min_diff = np.min(np.abs(
                                        z_aq - monitoringCoordZ))
                                    if min_diff > min_z_spacing:
                                        z_aq_list = z_aq.tolist()
                                        z_aq_list.append(monitoringCoordZ)
                                        z_aq = np.unique(z_aq_list)
                                    monitoringCoordZCheck = 0

                                elif np.min(z_sh) <= monitoringCoordZ < np.max(z_sh):
                                    min_diff = np.min(np.abs(
                                        z_sh - monitoringCoordZ))
                                    if min_diff > min_z_spacing:
                                        z_sh_list = z_sh.tolist()
                                        z_sh_list.append(monitoringCoordZ)
                                        z_sh = np.unique(z_sh_list)
                                    monitoringCoordZCheck = 0

                    if shaleRef == 0:
                        z_temp = np.concatenate(
                            (z_sh_lowest[:-1], z_aq[:-1], z_sh[:-1]), axis=0)

                    elif shaleRef not in [0, (numShaleLayers - 2)]:
                        z_temp = np.concatenate((z_temp[:-1], z_aq[:-1],
                                                 z_sh[:-1]), axis=0)

                    elif shaleRef == (numShaleLayers - 2):
                        z_temp = np.concatenate((z_temp[:-1], z_aq[:-1], z_sh[:]),
                                           axis=0)

                        if z is None:
                            z = np.zeros((len(z_temp), len(y_grid), len(x_grid)))

                        if len(z_temp) == z.shape[0]:
                            z[:, yRef, xRef] = z_temp[:]
                        elif len(z_temp) > z.shape[0]:
                            # If the changes in depth due to strike and dip make
                            # the number of z grid points higher than the number
                            # of points it had for the first x and y values
                            # considered, do not use the extra points.
                            z[:, yRef, xRef] = z_temp[0:z.shape[0]]
                        elif len(z_temp) < z.shape[0]:
                            num_extra_points = z.shape[0] - len(z_temp)
                            z_extra = np.linspace(z_temp[-1], lowest_depth,
                                               num_extra_points)

                            # If the changes in depth due to strike and dip make
                            # the number of z grid points lower than the number
                            # of points it had for the first x and y values
                            # considered, add extra points to compensate.
                            z_temp = np.concatenate((z_temp[:], z_extra[:]), axis=0)

                            z[:, yRef, xRef] = z_temp[:]

    elif var_type == 'noVariation' and not 'CarbonateAquifer' in aq_component_types:
        for shaleRef in range(numShaleLayers - 1):
            z_aq = np.linspace(-aquiferBottomDepths[shaleRef],
                               -aquiferBottomDepths[shaleRef]
                               + aquiferThicknesses[shaleRef], numPointsInAq)
            z_sh = np.linspace(z_aq[-1], z_aq[-1]
                               + shaleThicknesses[shaleRef + 1],
                               numPointsInShales)

            if shaleRef == 0:
                z_sh_lowest = np.linspace(
                    z_aq[0] - shaleThicknesses[shaleRef], z_aq[0],
                    numPointsInShales)

            # Check if a monitoringCoordZ value fits within the current
            # aquifer or shale. If so, add it.
            if monitoringCoordZ is not None:
                if isinstance(monitoringCoordZ, list):
                    for monitorRef, monitoringCoordZVal in enumerate(monitoringCoordZ):
                        if monitoringCoordZCheck[monitorRef]:
                            if np.min(z_aq) <= monitoringCoordZVal < np.max(z_aq):
                                min_diff = np.min(
                                    np.abs(z_aq - monitoringCoordZVal))
                                if min_diff > min_z_spacing:
                                    z_aq_list = z_aq.tolist()
                                    z_aq_list.append(monitoringCoordZVal)
                                    z_aq = np.unique(z_aq_list)
                                monitoringCoordZCheck[monitorRef] = 0

                            elif np.min(z_sh) <= monitoringCoordZVal < np.max(z_sh):
                                min_diff = np.min(
                                    np.abs(z_sh - monitoringCoordZVal))
                                if min_diff > min_z_spacing:
                                    z_sh_list = z_sh.tolist()
                                    z_sh_list.append(monitoringCoordZVal)
                                    z_sh = np.unique(z_sh_list)
                                monitoringCoordZCheck[monitorRef] = 0

                else:
                    # Only one monitoring location
                    if monitoringCoordZCheck:
                        if np.min(z_aq) <= monitoringCoordZ < np.max(z_aq):
                            min_diff = np.min(np.abs(z_aq - monitoringCoordZ))
                            if min_diff > min_z_spacing:
                                z_aq_list = z_aq.tolist()
                                z_aq_list.append(monitoringCoordZ)
                                z_aq = np.unique(z_aq_list)
                            monitoringCoordZCheck = 0

                        elif np.min(z_sh) <= monitoringCoordZ < np.max(z_sh):
                            min_diff = np.min(np.abs(
                                z_sh_list - monitoringCoordZ))
                            if min_diff > min_z_spacing:
                                z_sh_list = z_sh.tolist()
                                z_sh_list.append(monitoringCoordZ)
                                z_sh = np.unique(z_sh_list)
                            monitoringCoordZCheck = 0

            if shaleRef == 0:
                z = np.concatenate(
                    (z_sh_lowest[:-1], z_aq[:-1], z_sh[:-1]), axis=0)

            elif shaleRef not in [0, (numShaleLayers - 2)]:
                z = np.concatenate((z[:-1], z_aq[:-1], z_sh[:-1]), axis=0)

            elif shaleRef == (numShaleLayers - 2):
                z = np.concatenate((z[:-1], z_aq[:-1], z_sh[:]), axis=0)[:, None, None]

    elif 'CarbonateAquifer' in aq_component_types:
        z = None

    return z


def get_plume_timings(sm, s, time_array, sample, plume_metric_abbrev,
                      x_grid, y_grid, z_grid, aq_components, aq_component_types,
                      aq_component_xvals, aq_component_yvals, aq_component_indices,
                      numShaleLayers, shaleThicknesses, aquiferThicknesses,
                      reservoirThickness, strata_comp, analysis='lhs',
                      strike=None, dip=None, dipDirection=None,
                      coordxReferencePoint=None, coordyReferencePoint=None,
                      var_type='noVariation'):
    """
    Function to calculate earliest plume timings based on contaminant
    plume dimensions over time.

    This function can handle the plume metrics from the components
    FutureGen2Aquifer, FutureGen2AZMI, GenericAquifer, CarbonateAquifer,
    DeepAlluviumAquifer, DeepAlluviumAquiferML. Such metrics include TDS_dx,
    TDS_dy, and TDS_dz from FutureGen2Aquifer, Dissolved_salt_dr and
    Dissolved_salt_dz from GenericAquifer, and dx and dx from Carbonate Aquifer.
    The metrics to be used must be produced by (an) aquifer component(s).
    The metric to be used is specified with PlumeType in the TTFD section of
    the .yaml file. The PlumeType can be TDS, pH, Dissolved_salt,
    Dissolved_CO2, or CarbonateAquifer. CarbonateAquifer components produce
    plume dimensions representing the impacted region, and this output does not
    isolate the influences of TDS and pH (for example).

    :param sm: OpenIAM System model for which plots are created
    :type sm: openiam.SystemModel object

    :param s: SampleSet with results of analysis performed on the system model.
        None for forward analysis.
    :type s: matk.SampleSet object

    :param time_array: Array of the times assessed in the simulation(s)
    :type time_array: numpy.ndarray

    :param sample: index of the LHS sample (sample number - 1), or 0 for forward
    :type sample: int

    :param plume_metric_abbrev: prefix of the output plume metrics for the type
        of contaminant plume assessed (e.g., "TDS" for TDS_dx, TDS_dy, and
        TDS_dz or "Dissolved_CO2" for Dissolved_CO2_dr and Dissolved_CO2_dz).
    :type plume_metric_abbrev: str

    :param x_grid: Array of x values to be evaluated. The array has the shape
        (number of x values,).
    :type x_grid: numpy.ndarray

    :param y_grid: Array of y values to be evaluated. The array has the shape
        (number of y values, 1).
    :type y_grid: numpy.ndarray

    :param z_grid: Array of z values to be evaluated. The array has the shape
        (number of z values, 1, 1).
    :type z_grid: numpy.ndarray

    :param aq_components: List containing each of the aquifer components
    :type aq_components: list

    :param aq_component_xvals: List of the x values for aquifer components.
        Values should be ordered to correspond with each component in aq_components.
    :type aq_component_xvals: list

    :param aq_component_yvals: List of the y values for aquifer components.
        Values should be ordered to correspond with each component in aq_components.
    :type aq_component_yvals: list

    :param aq_component_indices: List of the aquifer indices for each of the
        aquifer components. For example, consider two aquifer components
        representing aquifer 1 and aquifer 2. If these components are linked
        to one well at one location, the list would be [0, 1]. If the aquifer
        components are linked to two wells at two locations, there are really
        four aquifer components being used and the list would be [0, 1, 0, 1]
        (correspondign with locations [#1, #1, #2, #2]). Values should be
        ordered to correspond with each component in aq_components.
    :type aq_component_indices: list

    :param numShaleLayers: number of shale layers
    :type numShaleLayers: int

    :param shaleThicknesses: List of the thickness for each shale layer, with
        the first entry representing the lowest shale and the last entry
        representing the highest shale
    :type shaleThicknesses: list

    :param aquiferThicknesses: List of the thickness [|m|] for each aquifer,
        with the first entry representing the lowest aquifer and the last entry
        representing the highest aquifer. These thicknesses are those at the
        reference point located at (coordxRefPoint, coordyRefPoint).
    :type aquiferThicknesses: list

    :param reservoirThickness: thickness of the reservoir [|m|] at the
        reference point.
    :type reservoirThickness: int or float

    :param strata_comp: stratigraphy component. If spatially variable
        stratigraphy is being used, this is the component created for the
        reference point.
    :type strata_comp: openiam.stratigraphy_component.Stratigraphy

    :param strike: Unit strike in degrees, where 0 is north, 90 is east, 180 is
        south, and 270 is west.
    :type strike: int, float, or None

    :param dip: Unit dip in degrees, where positive means dipping down in
        the dipDirection provided
    :type dip: int, float, or None

    :param dipDirection: Direction of the units dip, expressed as N, E, S, W,
        NE, SE, SW, or SW. Note that dipDirection must be compatible with the
        strike provided.
    :type dipDirection: str or None

    :param coordxReferencePoint: x value for the reference point. Only used if
    :type coordxReferencePoint: int, float, or None

    :param coordyReferencePoint: y value for the reference point. Only used if
        spatially variable stratigraphy is used.
    :type coordyReferencePoint: int, float, or None

    :param var_type: Option used for the domain stratigraphy. 'noVariation'
        means the strata are flat and unchanging over space while
        'strikeAndDip' means that unit thicknesses change according to strike,
        dip, and dipDirection values.
    :type var_type: str

    """
    # Go through the different types of aquifer components
    if 'GenericAquifer' in aq_component_types:
        # Initialize time to first detection grid to large number
        plumeTiming = np.ones((z_grid.shape[0], len(y_grid),
                               len(x_grid))) * MAX_TIME

        for well, (x0, y0) in enumerate(zip(aq_component_xvals, aq_component_yvals)):
            z0 = get_z0(x0, y0, strata_comp, numShaleLayers, shaleThicknesses[:],
                        aquiferThicknesses[:], reservoirThickness,
                        aquiferReference=aq_component_indices[well],
                        var_type=var_type, strike=strike, dip=dip,
                        dipDirection=dipDirection,
                        coordxReferencePoint=coordxReferencePoint,
                        coordyReferencePoint=coordyReferencePoint)

            # Get plume dimensions
            for indd, tt in enumerate(time_array):
                if analysis == 'forward':
                    r = sm.collect_observations_as_time_series(
                        aq_components[well], plume_metric_abbrev+'_dr')[indd]
                    h = sm.collect_observations_as_time_series(
                        aq_components[well], plume_metric_abbrev+'_dz')[indd]
                elif analysis in ['lhs', 'parstudy']:
                    r = s.recarray[aq_components[well].name + '.'
                                   + plume_metric_abbrev
                                   + '_dr_' + str(indd)][sample]
                    h = s.recarray[aq_components[well].name + '.'
                                   + plume_metric_abbrev
                                   + '_dz_' + str(indd)][sample]
                # Plume diameter
                d = r * 2

                # Calculate time to first detection
                if (d > 0 and h > 0):
                    if var_type == 'noVariation':
                        ellipsoid = ((x_grid - x0) / (d * 0.5)) ** 2 + (
                            (y_grid - y0) / (d * 0.5)) ** 2 + (
                                (z_grid - z0) / (h * 0.5)) ** 2 <= 1
                        in_aquifer = z_grid - z0 < 0
                        half_ellipsoid = np.logical_and(ellipsoid, in_aquifer)
                        mask = np.logical_not(half_ellipsoid) * MAX_TIME
                        plume = half_ellipsoid * tt
                        plumeTiming = np.minimum(plumeTiming, mask + plume)

                    elif var_type == 'strikeAndDip':
                        ellipsoid = (
                            (x_grid - x0) / (d * 0.5)) ** 2 + (
                                (y_grid - y0) / (d * 0.5)) ** 2 + (
                                    (z_grid- z0) / (h * 0.5)) ** 2 <= 1
                        in_aquifer = z_grid - z0 < 0
                        half_ellipsoid = np.logical_and(
                            ellipsoid, in_aquifer)
                        mask = np.logical_not(half_ellipsoid) * MAX_TIME
                        plume = half_ellipsoid * tt
                        plumeTiming = np.minimum(plumeTiming, mask + plume)

    elif 'CarbonateAquifer' in aq_component_types:
        # Initialize time to first detection grid to large number
        plumeTiming = np.ones((len(x_grid), len(y_grid))) * MAX_TIME

        for well, (x0, y0) in enumerate(zip(aq_component_xvals, aq_component_yvals)):
            # These are remade for each sample, but getting the values
            # from the right component is critical.
            ind_list = list(range(len(time_array)))

            if analysis == 'forward':
                dxs = sm.collect_observations_as_time_series(
                    aq_components[well], 'dx')
                dys = sm.collect_observations_as_time_series(
                    aq_components[well], 'dy')

                dx = dxs[:]
                dy = dys[:]

            elif analysis in ['lhs', 'parstudy']:
                dxs = np.array([s.recarray[aq_components[well].name + '.dx_'
                                           + str(indd)] for indd in ind_list])
                dys = np.array([s.recarray[aq_components[well].name + '.dy_'
                                           + str(indd)] for indd in ind_list])
                dx = dxs[:, sample]
                dy = dys[:, sample]

            # Calculate earliest plume timing
            for tt, a, b in zip(time_array, dx, dy):
                if (a > 0 and b > 0):
                    # True for points inside the ellipse
                    ellipse = (((x_grid - x0) / (a * 0.5)) ** 2) + (
                        ((y_grid - y0) / (b * 0.5)) ** 2) <= 1
                    mask = np.logical_not(ellipse) * MAX_TIME
                    plume = ellipse * tt
                    plumeTiming = np.minimum(plumeTiming, mask + plume)

    else:
        # Get the plume timings from FutureGen2Aquifer, FutureGen2AZMI,
        # DeepAlluviumAquifer, and DeepAlluviumAquiferML components
        # Initialize timings grid to large number
        plumeTiming = np.ones((z_grid.shape[0], len(y_grid), len(x_grid))) * MAX_TIME

        for well, (x0, y0) in enumerate(zip(aq_component_xvals, aq_component_yvals)):
            z0 = get_z0(x0, y0, strata_comp, numShaleLayers, shaleThicknesses[:],
                        aquiferThicknesses[:], reservoirThickness,
                        aquiferReference=aq_component_indices[well],
                        var_type=var_type, strike=strike,
                        dip=dip, dipDirection=dipDirection,
                        coordxReferencePoint=coordxReferencePoint,
                        coordyReferencePoint=coordyReferencePoint)

            # Get plume dimensions
            for indd, tt in enumerate(time_array):
                if analysis == 'forward':
                    a = sm.collect_observations_as_time_series(
                        aq_components[well], plume_metric_abbrev + '_dx')[indd]
                    b = sm.collect_observations_as_time_series(
                        aq_components[well], plume_metric_abbrev + '_dy')[indd]
                    c = sm.collect_observations_as_time_series(
                        aq_components[well], plume_metric_abbrev + '_dz')[indd]
                elif analysis in ['lhs', 'parstudy']:
                    a = s.recarray[aq_components[well].name + '.'
                                   + plume_metric_abbrev
                                   + '_dx_' + str(indd)][sample]
                    b = s.recarray[aq_components[well].name + '.'
                                   + plume_metric_abbrev
                                   + '_dy_' + str(indd)][sample]
                    c = s.recarray[aq_components[well].name + '.'
                                   + plume_metric_abbrev
                                   + '_dz_' + str(indd)][sample]

                # Calculate earliest plume timing
                if (a > 0 and b > 0 and c > 0):
                    if var_type == 'noVariation':
                        ellipsoid = ((x_grid - x0) / (a * 0.5)) ** 2 + (
                            (y_grid - y0)/(b * 0.5)) ** 2 + (
                                (z_grid - z0) / (c * 0.5)) ** 2 <= 1
                        in_aquifer = z_grid - z0 < 0
                        half_ellipsoid = np.logical_and(ellipsoid, in_aquifer)
                        mask = np.logical_not(half_ellipsoid) * MAX_TIME
                        plume = half_ellipsoid * tt
                        plumeTiming = np.minimum(plumeTiming, mask + plume)

                    elif var_type == 'strikeAndDip':
                        ellipsoid = (
                            (x_grid - x0) / (a * 0.5)) ** 2 + (
                                (y_grid - y0) / (b * 0.5)) ** 2 + (
                                    (z_grid- z0) / (c * 0.5)) ** 2 <= 1
                        in_aquifer = z_grid - z0 < 0
                        half_ellipsoid = np.logical_and(
                            ellipsoid, in_aquifer)
                        mask = np.logical_not(half_ellipsoid) * MAX_TIME
                        plume = half_ellipsoid * tt
                        plumeTiming = np.minimum(plumeTiming, mask + plume)

    return plumeTiming


def get_monitoring_location_ttfd(plumeTimings, x_grid, y_grid, z_grid,
                                 monitoringCoordX, monitoringCoordY,
                                 monitoringCoordZ, monitoringXIndex,
                                 monitoringYIndex, monitorHorizontalWindow,
                                 monitorVerticalWindow, aq_component_types,
                                 var_type='noVariation'):
    """
    Function that takes plume timing results and monitoring location data and
    returns lists of time to first detection (ttfd) at each monitoring location.
    The lists returned contain the actual ttfd as well as the corresponding
    sensors' x, y, and z values.
    """
    # List of TTFD values at or sufficiently close to monitoring location(s)
    ttfd_list = []

    # List of x, y, and z values corresponding with the ttfd_list values
    ttfd_x_list = []
    ttfd_y_list = []

    if 'CarbonateAquifer' in aq_component_types:
        checkCarbAq = True
        ttfd_z_list = None
    else:
        checkCarbAq = False
        ttfd_z_list = []

    if isinstance(monitoringCoordX, list) and isinstance(monitoringCoordY, list):
        for monitorRef, monitoringXIndexVal in enumerate(monitoringXIndex):
            if checkCarbAq:
                # Find any TTFD within HorizontalWindow of the Monitoring Location
                for xRef, xVal in enumerate(x_grid):
                    for yRef, yVal in enumerate(y_grid[:, 0]):
                        distance_temp = (
                            ((monitoringCoordX[monitorRef] - xVal) ** 2) \
                                + ((monitoringCoordY[monitorRef] - yVal) ** 2)) ** 0.5

                        if distance_temp <= monitorHorizontalWindow:
                            ttfdVal = plumeTimings[xRef, yRef]

                            if ttfdVal < THRESHOLD_TIME:
                                ttfd_x_list.append(xVal)
                                ttfd_y_list.append(yVal)
                                ttfd_list.append(ttfdVal)

            else:
                # Output from any aquifer components besides CarbonateAquifer
                # Find any TTFD within HorizontalWindow and VerticalWindow
                # of the monitoring location.
                for xRef, xVal in enumerate(x_grid):
                    for yRef, yVal in enumerate(y_grid[:, 0]):
                        distance_temp = (
                            ((monitoringCoordX[monitorRef] - xVal) ** 2) \
                                + ((monitoringCoordY[monitorRef] - yVal) ** 2)) ** 0.5

                        if distance_temp <= monitorHorizontalWindow:
                            for zRef in range(z_grid.shape[0]):
                                if var_type == 'noVariation':
                                    z_grid_val = z_grid[zRef, 0, 0]
                                elif var_type == 'strikeAndDip':
                                    z_grid_val = z_grid[
                                        zRef, monitoringYIndex[monitorRef],
                                        monitoringXIndexVal]

                                # These should both be negative, so take the abs()
                                vertical_distance_temp = np.abs(
                                    np.abs(z_grid_val)
                                    - np.abs(monitoringCoordZ[monitorRef]))

                                if vertical_distance_temp <= monitorVerticalWindow:
                                    ttfdVal = plumeTimings[zRef, yRef, xRef]

                                    if ttfdVal < THRESHOLD_TIME:
                                        ttfd_list.append(ttfdVal)

                                        ttfd_x_list.append(xVal)
                                        ttfd_y_list.append(yVal)
                                        ttfd_z_list.append(z_grid_val)
    else:
        # Single monitoring location, not a list
        if checkCarbAq:
            # Find any TTFD within HorizontalWindow of the Monitoring Location
            for xRef, xVal in enumerate(x_grid):
                for yRef, yVal in enumerate(y_grid[:, 0]):
                    distance_temp = (
                        ((monitoringCoordX - xVal) ** 2) + (
                            (monitoringCoordY - yVal) ** 2)) ** 0.5

                    if distance_temp <= monitorHorizontalWindow:
                        ttfdVal = plumeTimings[xRef, yRef]

                        if ttfdVal < THRESHOLD_TIME:
                            ttfd_list.append(ttfdVal)
                            ttfd_x_list.append(x_grid[monitoringXIndex])
                            ttfd_y_list.append(y_grid[monitoringYIndex, 0])
        else:
            # Output from any aquifer compoent besides CarbonateAquifer
            # Find any TTFD within HorizontalWindow and VerticalWindow
            # of the Monitoring Location
            for xRef, xVal in enumerate(x_grid):
                for yRef, yVal in enumerate(y_grid[:, 0]):
                    distance_temp = np.abs(
                        (((monitoringCoordX - xVal) ** 2) + (
                            (monitoringCoordY - yVal) ** 2)) ** 0.5)

                    if distance_temp <= monitorHorizontalWindow:
                        for zRef in range(z_grid.shape[0]):
                            if var_type == 'noVariation':
                                z_grid_val = z_grid[zRef, 0, 0]
                            elif var_type == 'strikeAndDip':
                                z_grid_val = z_grid[zRef, monitoringYIndex,
                                                    monitoringXIndex]

                            # These should both be negative, so take the abs()
                            vertical_distance_temp = np.abs(
                                np.abs(z_grid_val) - np.abs(monitoringCoordZ))

                            if vertical_distance_temp <= monitorVerticalWindow:
                                ttfdVal = plumeTimings[zRef, yRef, xRef]

                                if ttfdVal < THRESHOLD_TIME:
                                    ttfd_list.append(ttfdVal)

                                    ttfd_x_list.append(xVal)
                                    ttfd_y_list.append(yVal)
                                    ttfd_z_list.append(z_grid_val)

    return ttfd_list, ttfd_x_list, ttfd_y_list, ttfd_z_list


def get_monitors_within_z_range(monitoringCoordX, monitoringCoordY, monitoringCoordZ,
                                min_z_subplot, max_z_subplot):
    """
    Function that returns the monitoring locations that are within a certain
    range of z values.
    """
    monitoringCoordX_temp = []
    monitoringCoordY_temp = []
    monitoringCoordZ_temp = []

    if isinstance(monitoringCoordZ, (int, float)):
        if min_z_subplot <= monitoringCoordZ < max_z_subplot:
            monitoringCoordX_temp = [monitoringCoordX]
            monitoringCoordY_temp = [monitoringCoordY]
            monitoringCoordZ_temp = [monitoringCoordZ]

    else:
        for xVal, yVal, zVal in zip(monitoringCoordX, monitoringCoordY, monitoringCoordZ):
            if min_z_subplot <= zVal < max_z_subplot:
                monitoringCoordX_temp.append(xVal)
                monitoringCoordY_temp.append(yVal)
                monitoringCoordZ_temp.append(zVal)

    return  monitoringCoordX_temp, monitoringCoordY_temp, monitoringCoordZ_temp


def get_xandy_monitoring_indices(monitoringCoordX, monitoringCoordY, x_grid,
                                 y_grid):
    """
    Function that returns the indices for the monitoring location(s) given the
    x_grid and y_grid values.
    """
    if isinstance(monitoringCoordX, (int, float)):
        xMinProximity = None
        for xRef, xVal in enumerate(x_grid):
            xProximity_temp = np.abs(xVal - monitoringCoordX)

            if xMinProximity is None:
                xMinProximity = xProximity_temp
                monitoringXIndex = xRef
            elif xProximity_temp < xMinProximity:
                xMinProximity = xProximity_temp
                monitoringXIndex = xRef

    else:
        monitoringXIndex = []
        for monitoringCoordXVal in monitoringCoordX:
            xMinProximity = None
            for xRef, xVal in enumerate(x_grid):
                xProximity_temp = np.abs(xVal - monitoringCoordXVal)

                if xMinProximity is None:
                    xMinProximity = xProximity_temp
                    monitoringXIndex_current = xRef
                elif xProximity_temp < xMinProximity:
                    xMinProximity = xProximity_temp
                    monitoringXIndex_current = xRef

            monitoringXIndex.append(monitoringXIndex_current)

    if isinstance(monitoringCoordY, (int, float)):
        yMinProximity = None
        for yRef, yVal in enumerate(y_grid[:, 0]):
            yProximity_temp = np.abs(yVal - monitoringCoordY)

            if yMinProximity is None:
                yMinProximity = yProximity_temp
                monitoringYIndex = yRef
            elif yProximity_temp < yMinProximity:
                yMinProximity = yProximity_temp
                monitoringYIndex = yRef

    else:
        monitoringYIndex = []
        for _, monitoringCoordYVal in enumerate(monitoringCoordY):
            yMinProximity = None
            for yRef, yVal in enumerate(y_grid[:, 0]):
                yProximity_temp = np.abs(yVal - monitoringCoordYVal)

                if yMinProximity is None:
                    yMinProximity = yProximity_temp
                    monitoringYIndex_current = yRef
                elif yProximity_temp < yMinProximity:
                    yMinProximity = yProximity_temp
                    monitoringYIndex_current = yRef

            monitoringYIndex.append(monitoringYIndex_current)

    return monitoringXIndex, monitoringYIndex


def get_z0(x_val, y_val, stratigraphy_comp, numShaleLayers,
           shaleThicknessList, aquiferThicknessList, reservoirThickness,
           aquiferReference=0, var_type='noVariation', strike=None, dip=None,
           dipDirection=None, coordxReferencePoint=0,
           coordyReferencePoint=0):
    """
    Function that returns the z0 value used in get_plume_timings().
    """
    if var_type == 'noVariation':
        aquiferBottomDepths = get_aquifer_bottom_depths(numShaleLayers,
                                                        shaleThicknessList,
                                                        aquiferThicknessList)

        z0 = -aquiferBottomDepths[aquiferReference] \
            + aquiferThicknessList[aquiferReference]

    elif var_type == 'strikeAndDip':
        updatedStratigraphy = strata.update_stratigraphy_by_strike_and_dip(
            numberOfShaleLayers=numShaleLayers,
            shaleThicknessList=shaleThicknessList[:],
            aquiferThicknessList=aquiferThicknessList[:],
            reservoirThickness=reservoirThickness,
            strike=strike, dip=dip, dipDirection=dipDirection,
            coordxRefPoint=coordxReferencePoint,
            coordyRefPoint=coordyReferencePoint,
            location_x=x_val, location_y=y_val,
            strataRefPoint=stratigraphy_comp)

        shaleThicknessesUpdated = updatedStratigraphy['shaleThicknessList']
        aquiferThicknessesUpdated = updatedStratigraphy['aquiferThicknessList']

        aquiferBottomDepths = get_aquifer_bottom_depths(numShaleLayers,
                                                        shaleThicknessesUpdated,
                                                        aquiferThicknessesUpdated)

        z0 = -aquiferBottomDepths[aquiferReference] \
            + aquiferThicknessesUpdated[aquiferReference]

    return z0


def write_dream_grid(x_grid, y_grid, z_grid, output_dir, var_type='noVariation'):
    """
    Function that writes x_grid, y_grid, and z_grids for DREAM. Note that
    lengths are converted from meters to feet.
    """
    # Write grid in feet for DREAM
    filename = os.sep.join([output_dir, 'iam.grid'])

    with open(filename, "w+") as f:

        if z_grid is None:
            f.write('x_feet, y_feet, \n')

            for xVal in x_grid:
                f.write('{},'.format(xVal / 0.3048))
                for yVal in y_grid[:, 0]:
                    f.write('{},{},\n'.format(
                        xVal / 0.3048, yVal / 0.3048))

        else:
            f.write('x_feet, y_feet, z_feet, \n')

            for xRef, xVal in enumerate(x_grid):
                f.write('{},'.format(xVal / 0.3048))
                for yRef, yVal in enumerate(y_grid[:, 0]):
                    for zRef in range(z_grid.shape[0]):
                        if var_type == 'noVariation':
                            f.write('{},{},{},\n'.format(
                                xVal / 0.3048, yVal / 0.3048,
                                z_grid[zRef, 0, 0] / 0.3048))
                        elif var_type == 'strikeAndDip':
                            f.write('{},{},{},\n'.format(
                                xVal / 0.3048, yVal / 0.3048,
                                z_grid[zRef, yRef, xRef] / 0.3048))


def write_dream_output(plumeTimings, x_grid, y_grid, z_grid, plume_metric_abbrev,
                       output_dir, realization=0, var_type='noVariation',
                       aq_component_ithresh=None, checkCarbAq=False):
    """
    Function to write output for DREAM. Note that lengths are assumed to
    initially be in meters. The output for DREAM then include lengths converted
    to feet.
    """
    # Write output (in feet) for DREAM
    filename = os.sep.join([output_dir, 'ttfd_{}_{}.iam'.format(
        plume_metric_abbrev, realization)])

    if checkCarbAq:
        with open(filename, "w+") as f:
            if np.max(aq_component_ithresh) == 1:
                threshold_CarbAq = 6.5
            elif np.min(aq_component_ithresh) == 2:
                threshold_CarbAq = 6.6
            else:
                threshold_CarbAq = 'both_ithresh_values_were_used_by_Carbonate_Aquifer_components'

            f.write('IAM,{},{},{},{},\n'.format(
                realization, plume_metric_abbrev, INDICATOR[plume_metric_abbrev],
                threshold_CarbAq))

            f.write('x_feet, y_feet, ttfd_days,\n')

            ix = -1
            for xx in x_grid:
                ix = ix + 1
                iy = -1
                for yy in y_grid:
                    iy = iy + 1
                    checkValid = check_metric_validity(
                        plumeTimings[ix][iy], 'plumeTimings', valType = 'single')
                    if checkValid:
                        f.write('{},{},{}\n'.format(xx / 0.3048, yy[0] / 0.3048,
                                                    plumeTimings[ix][iy]))

    else:
        with open(filename, "w+") as f:
            f.write('IAM,{},{},{},{},\n'.format(
                realization, plume_metric_abbrev, INDICATOR[plume_metric_abbrev],
                THRESHOLDS[plume_metric_abbrev]))

            if var_type == 'noVariation':
                f.write('x_feet, y_feet, z_feet, ttfd_days,\n')
            elif var_type == 'strikeAndDip':
                f.write('x_feet, y_feet, z_feet, ttfd_days, Note: the strikeAndDip ',
                        'option was used for stratigraphy so z values vary with x and y,\n')

            for xRef, xVal in enumerate(x_grid):
                for yRef, yVal in enumerate(y_grid[:, 0]):
                    for zRef in range(z_grid.shape[0]):
                        checkValid = check_metric_validity(
                            plumeTimings[zRef][yRef][xRef], 'plumeTimings',
                            valType = 'single')

                        if checkValid:
                            if var_type == 'noVariation':
                                f.write('{},{},{},{}\n'.format(
                                    xVal / 0.3048, yVal / 0.3048,
                                    z_grid[zRef, 0, 0] / 0.3048,
                                    plumeTimings[zRef][yRef][xRef]))

                            elif var_type == 'strikeAndDip':
                                f.write('{},{},{},{}\n'.format(
                                    xVal / 0.3048, yVal / 0.3048,
                                    z_grid[zRef, yRef, xRef] / 0.3048,
                                    plumeTimings[zRef][yRef][xRef]))


def get_ttfd_plot_yaml_input(yaml_data, name):
    """
    Function that reads the TTFD plot input provided in the .yaml file and
    returns a dictionary containing the input.

    :param yaml_data: Dictionary of input values from the entire .yaml file
    :type yaml_data: dict

    :param name: Name of the TTFD plot provided in the Plots section of the
        .yaml file
    :type name: str

    :returns: TTFD_yaml_input_dict
    """
    plume_type = yaml_data['Plots'][name]['TTFD']['PlumeType']

    defaultFigureDPI = 100
    FigureDPI = defaultFigureDPI

    defaultHorizontalWindow = 1
    HorizontalWindow = defaultHorizontalWindow
    defaultVerticalWindow = 1
    VerticalWindow = defaultVerticalWindow

    defaultNumPointsInAq = 10
    numPointsInAq = defaultNumPointsInAq
    defaultNumPointsInShales = 3
    numPointsInShales = defaultNumPointsInShales

    NumZPointInAqsWarning = False
    NumZPointInShalesWarning = False

    gridXLimsWarning = False
    gridYLimsWarning = False

    defaultSaveCSVFiles = True
    saveCSVFiles = defaultSaveCSVFiles
    defaultWriteDreamOutput = False
    write_DREAM_output = defaultWriteDreamOutput

    windowType_debug_msg = ''.join([
        'The {}Window provided for the monitoring locations in the TTFD plot ',
        name, ' was not of type int or float. This ',
        'parameter will be set to the default value of {} m.'])

    num_zPoints_debug_msg = ''.join([
        'The NumZPointsWithin{} provided in the TTFD plot ', name, ' was not ',
        'of type int. This parameter will be set to the default value of {}.'])

    grid_spacing_debug_msg = ''.join([
        'The {}GridSpacing provided in the TTFD plot ', name, ' was not of type ',
        'int or type float. The {}_grid will be created in the default manner. ',
        'Check your inputs in the .yaml file.'])

    grid_limits_debug_msg = ''.join([
        'The limits provided for the {}_grid used for the TTFD plot ', name,
        ' is not a list of length 2 (SpecifyXandYGridLims: grid{}Lims). The ',
        '{}_grid will be created in the default manner. Check your inputs in ',
        'the .yaml file.'])

    InjectionCoord_debug_msg = ''.join([
        'InjectionCoord{} was provided for the TTFD plot ', name, ', but not ',
        'InjectionCoord{}. Check your input. Injection sites will not be ',
        'displayed.'])

    axis_lims_debug_msg = ''.join([
        'The {}-axis limits provided for the TTFD plot ', name,
        ' (SpecifyXandYLims: {}Lims) are not of length 2 and will not be used. ',
        'Check your inputs in the .yaml file.'])

    bool_option_debug_msg = ''.join([
        'The input provided for {} in the TTFD plot ', name,
        ' was not of type boolean. {} will be set to the default value of {}.'])

    defaultMonitoringTTFD = False
    monitoringTTFD = defaultMonitoringTTFD
    monitoringCoordX = None
    monitoringCoordY = None
    monitoringCoordZ = None

    defaultPlotInjectionSites = False
    plot_injection_sites = defaultPlotInjectionSites

    EnforceXandYLims = False
    xLims = None
    yLims = None

    EnforceGridXandYLims = False
    gridXLims = None
    gridYLims = None

    defaultPlotInjectionSites = False
    plot_injection_sites = defaultPlotInjectionSites

    if 'FigureDPI' in yaml_data['Plots'][name]['TTFD']:
        FigureDPI = yaml_data['Plots'][name]['TTFD']['FigureDPI']

        if not isinstance(FigureDPI, (int, float)):
            debug_msg = ''.join([
                'The FigureDPI provided for the TTFD plot ', name, ' was not of ',
                'type int or float. The dpi (dots-per-inch) will be set to the ',
                'default value of {}']).format(str(defaultFigureDPI))
            logging.debug(debug_msg)
            FigureDPI = defaultFigureDPI

    if 'MonitoringLocations' in yaml_data['Plots'][name]['TTFD']:
        monitoringTTFD = True
        monitoringCoordX = yaml_data['Plots'][name]['TTFD'][
            'MonitoringLocations']['coordx']
        monitoringCoordY = yaml_data['Plots'][name]['TTFD'][
            'MonitoringLocations']['coordy']

        if plume_type != 'CarbonateAquifer':
            monitoringCoordZ = yaml_data['Plots'][name]['TTFD'][
                'MonitoringLocations']['coordz']

        if 'HorizontalWindow' in yaml_data['Plots'][name]['TTFD'][
                'MonitoringLocations']:
            HorizontalWindow = yaml_data['Plots'][name]['TTFD'][
                'MonitoringLocations']['HorizontalWindow']
            if not isinstance(HorizontalWindow, (float, int)):
                debug_msg = windowType_debug_msg.format(
                    'Horizontal', str(defaultHorizontalWindow))
                logging.debug(debug_msg)
                HorizontalWindow = defaultHorizontalWindow

        if plume_type != 'CarbonateAquifer':
            if 'VerticalWindow' in yaml_data['Plots'][name]['TTFD'][
                    'MonitoringLocations']:
                VerticalWindow = yaml_data['Plots'][name]['TTFD'][
                    'MonitoringLocations']['VerticalWindow']

                if not isinstance(VerticalWindow, (float, int)):
                    debug_msg = windowType_debug_msg.format(
                        'Vertical', str(defaultVerticalWindow))
                    logging.debug(debug_msg)
                    VerticalWindow = defaultVerticalWindow

        if not isinstance(monitoringCoordX, (float, int, list)):
            debug_msg = gen_monitor_loc_type_warning_msg('coordx', name)
            logging.debug(debug_msg)
            monitoringTTFD = defaultMonitoringTTFD

        if not isinstance(monitoringCoordY, (float, int, list)):
            debug_msg = gen_monitor_loc_type_warning_msg('coordy', name)
            logging.debug(debug_msg)
            monitoringTTFD = defaultMonitoringTTFD

        if plume_type != 'CarbonateAquifer':
            if not isinstance(monitoringCoordZ, (float, int, list)):
                debug_msg = gen_monitor_loc_type_warning_msg('coordz', name)
                logging.debug(debug_msg)
                monitoringTTFD = defaultMonitoringTTFD

            if np.max(monitoringCoordZ) > 0:
                debug_msg = ''.join([
                    'The coordz values provided for MonitoringLocations in the ',
                    'TTFD plot ', name, ' had (a) positive value(s). The depth ',
                    'values represented by coordz are taken as negative when ',
                    'beneath the surface. Check your input. The TTFD at the actual ',
                    'monitoring locations will not be plotted.'])
                logging.debug(debug_msg)
                monitoringTTFD = defaultMonitoringTTFD

        if plume_type != 'CarbonateAquifer':
            if isinstance(monitoringCoordX, list) and isinstance(monitoringCoordY, list) \
                    and isinstance(monitoringCoordZ, list):
                if len(monitoringCoordX) != len(monitoringCoordY) \
                        or len(monitoringCoordX) != len(monitoringCoordZ):
                    debug_msg = gen_monitoring_loc_len_warning(plume_type, name)
                    logging.debug(debug_msg)
                    monitoringTTFD = defaultMonitoringTTFD
        else:
            if isinstance(monitoringCoordX, list) and isinstance(monitoringCoordY, list):
                if len(monitoringCoordX) != len(monitoringCoordY):
                    debug_msg = gen_monitoring_loc_len_warning(plume_type, name)
                    logging.debug(debug_msg)
                    monitoringTTFD = defaultMonitoringTTFD

    if 'NumZPointsWithinAquifers' in yaml_data['Plots'][name]['TTFD']:
        numPointsInAq = yaml_data['Plots'][name]['TTFD']['NumZPointsWithinAquifers']

        if isinstance(numPointsInAq, list):
            if len(numPointsInAq) == 1:
                numPointsInAq = numPointsInAq[0]

                if not isinstance(numPointsInAq, int):
                    NumZPointInAqsWarning = True
            else:
                NumZPointInAqsWarning = True

        elif not isinstance(numPointsInAq, int):
            NumZPointInAqsWarning = True

    if NumZPointInAqsWarning:
        debug_msg = num_zPoints_debug_msg.format(
            'Aquifers', str(defaultNumPointsInAq))
        logging.debug(debug_msg)
        numPointsInAq = defaultNumPointsInAq

    if 'NumZPointsWithinShales' in yaml_data['Plots'][name]['TTFD']:
        numPointsInShales = yaml_data['Plots'][name]['TTFD']['NumZPointsWithinShales']

        if isinstance(numPointsInShales, list):
            if len(numPointsInShales) == 1:
                numPointsInShales =  numPointsInShales[0]

                if not isinstance(numPointsInShales, int):
                    NumZPointInShalesWarning = True
            else:
                NumZPointInShalesWarning = True

        elif not isinstance(numPointsInShales, int):
            NumZPointInShalesWarning = True

    if NumZPointInShalesWarning:
        debug_msg = num_zPoints_debug_msg.format(
            'Shales', str(defaultNumPointsInShales))
        logging.debug(debug_msg)
        numPointsInShales = defaultNumPointsInShales

    # These are used to make the x and y grids
    if 'xGridSpacing' in yaml_data['Plots'][name]['TTFD']:
        x_grid_spacing = yaml_data['Plots'][name]['TTFD']['xGridSpacing']

        if not isinstance(x_grid_spacing, (int, float)):
            debug_msg = grid_spacing_debug_msg.format('x', 'x')
            logging.debug(debug_msg)
            x_grid_spacing = None
    else:
        x_grid_spacing = None

    if 'yGridSpacing' in yaml_data['Plots'][name]['TTFD']:
        y_grid_spacing = yaml_data['Plots'][name]['TTFD']['yGridSpacing']

        if not isinstance(x_grid_spacing, (int, float)):
            debug_msg = grid_spacing_debug_msg.format('y', 'y')
            logging.debug(debug_msg)
            y_grid_spacing = None
    else:
        y_grid_spacing = None

    if 'SpecifyXandYGridLims' in yaml_data['Plots'][name]['TTFD']:
        EnforceGridXandYLims = True
        gridXLims = yaml_data['Plots'][name]['TTFD']['SpecifyXandYGridLims'][
            'gridXLims']
        gridYLims = yaml_data['Plots'][name]['TTFD']['SpecifyXandYGridLims'][
            'gridYLims']

        if isinstance(gridXLims, list):
            if len(gridXLims) != 2:
                gridXLimsWarning = True

        if gridXLimsWarning:
            debug_msg = grid_limits_debug_msg.format('x', 'X', 'x')
            logging.debug(debug_msg)
            EnforceGridXandYLims = False
            gridXLims = None

        if isinstance(gridYLims, list):
            if len(gridYLims) != 2:
                gridYLimsWarning = True

        if gridYLimsWarning:
            debug_msg = grid_limits_debug_msg.format('y', 'Y', 'y')
            logging.debug(debug_msg)
            EnforceGridXandYLims = False
            gridYLims = None

    if 'PlotInjectionSites' in yaml_data['Plots'][name]['TTFD']:
        plot_injection_sites = yaml_data['Plots'][name]['TTFD']['PlotInjectionSites']
        if not isinstance(plot_injection_sites, bool):
            debug_msg = ''.join(['The input provided for PlotInjectionSites ',
                                 'in  TTFD section of the .yaml file was ',
                                 'not of type bool. PlotInjectionSites will ',
                                 'be set to the default value of False.'])
            logging.debug(debug_msg)
            plot_injection_sites = defaultPlotInjectionSites

    InjectionCoordx = None
    InjectionCoordy = None
    if 'InjectionCoordx' in yaml_data['Plots'][name]['TTFD']:
        InjectionCoordx = yaml_data['Plots'][name]['TTFD']['InjectionCoordx']

    if 'InjectionCoordy' in yaml_data['Plots'][name]['TTFD']:
        InjectionCoordy = yaml_data['Plots'][name]['TTFD']['InjectionCoordy']
        if InjectionCoordx is None:
            debug_msg = InjectionCoord_debug_msg.format('y', 'x')
            logging.debug(debug_msg)
            plot_injection_sites = defaultPlotInjectionSites

        elif isinstance(InjectionCoordx, list) and isinstance(InjectionCoordy, list):

            if len(InjectionCoordx) != len(InjectionCoordy):
                debug_msg = ''.join(['The InjectionCoordy provided was not of ',
                                     'the same length as the InjectionCoordx ',
                                     'provided. Check your input. Injection ',
                                     'sites will not be displayed.'])
                logging.debug(debug_msg)
                plot_injection_sites = False

    if InjectionCoordx is not None and InjectionCoordy is None:
        debug_msg = InjectionCoord_debug_msg.format('x', 'y')
        logging.debug(debug_msg)
        plot_injection_sites = defaultPlotInjectionSites

    if plot_injection_sites and InjectionCoordx is not None \
            and InjectionCoordy is not None:
        try:
            InjectionCoordx = float(InjectionCoordx)
        except TypeError:
            InjectionCoordx = np.array(list(InjectionCoordx), dtype=float)

        try:
            InjectionCoordy = float(InjectionCoordy)
        except TypeError:
            InjectionCoordy = np.array(list(InjectionCoordy), dtype=float)

    elif not plot_injection_sites:
        InjectionCoordx = None
        InjectionCoordy = None

    if 'SpecifyXandYLims' in yaml_data['Plots'][name]['TTFD']:
        EnforceXandYLims = True
        xLims = yaml_data['Plots'][name]['TTFD']['SpecifyXandYLims']['xLims']
        yLims = yaml_data['Plots'][name]['TTFD']['SpecifyXandYLims']['yLims']

        if len(xLims) != 2:
            debug_msg = axis_lims_debug_msg.format('x', 'x')
            logging.debug(debug_msg)
            EnforceXandYLims = False

        if len(yLims) != 2:
            debug_msg = axis_lims_debug_msg.format('y', 'y')
            logging.debug(debug_msg)
            EnforceXandYLims = False

    if 'SaveCSVFiles' in yaml_data['Plots'][name]['TTFD']:
        saveCSVFiles = yaml_data['Plots'][name]['TTFD']['SaveCSVFiles']

        if not isinstance(saveCSVFiles, bool):
            debug_msg = bool_option_debug_msg.format(
                'SaveCSVFiles', 'SaveCSVFiles', str(defaultSaveCSVFiles))
            logging.debug(debug_msg)
            saveCSVFiles = defaultSaveCSVFiles

    if 'WriteDreamOutput' in yaml_data['Plots'][name]['TTFD']:
        write_DREAM_output = yaml_data['Plots'][name]['TTFD']['WriteDreamOutput']

        if not isinstance(write_DREAM_output, bool):
            debug_msg = bool_option_debug_msg.format(
                'WriteDreamOutput', 'WriteDreamOutput', str(defaultWriteDreamOutput))
            logging.debug(debug_msg)
            write_DREAM_output = defaultWriteDreamOutput

    TTFD_yaml_input_dict = dict()

    TTFD_yaml_input_dict['FigureDPI'] = FigureDPI

    TTFD_yaml_input_dict['monitoringTTFD'] = monitoringTTFD
    TTFD_yaml_input_dict['monitoringCoordX'] = monitoringCoordX
    TTFD_yaml_input_dict['monitoringCoordY'] = monitoringCoordY
    TTFD_yaml_input_dict['monitoringCoordZ'] = monitoringCoordZ
    TTFD_yaml_input_dict['HorizontalWindow'] = HorizontalWindow
    TTFD_yaml_input_dict['VerticalWindow'] = VerticalWindow

    TTFD_yaml_input_dict['plot_injection_sites'] = plot_injection_sites

    TTFD_yaml_input_dict['InjectionCoordx'] = InjectionCoordx
    TTFD_yaml_input_dict['InjectionCoordy'] = InjectionCoordy

    TTFD_yaml_input_dict['EnforceXandYLims'] = EnforceXandYLims
    TTFD_yaml_input_dict['xLims'] = xLims
    TTFD_yaml_input_dict['yLims'] = yLims

    TTFD_yaml_input_dict['EnforceGridXandYLims'] = EnforceGridXandYLims
    TTFD_yaml_input_dict['gridXLims'] = gridXLims
    TTFD_yaml_input_dict['gridYLims'] = gridYLims

    TTFD_yaml_input_dict['x_grid_spacing'] = x_grid_spacing
    TTFD_yaml_input_dict['y_grid_spacing'] = y_grid_spacing

    TTFD_yaml_input_dict['numPointsInAq'] = numPointsInAq
    TTFD_yaml_input_dict['numPointsInShales'] = numPointsInShales

    TTFD_yaml_input_dict['saveCSVFiles'] = saveCSVFiles

    TTFD_yaml_input_dict['write_DREAM_output'] = write_DREAM_output

    return TTFD_yaml_input_dict


def gen_monitor_loc_type_warning_msg(loc_type, name):
    """
    Function generating the warning message used in get_ttfd_plot_yaml_input()
    for the input types of monitoring coordx, coordy, and coordz values.
    """
    debug_msg = ''.join([
        'The {} values provided for monitoring locations in ',
        'the TTFD plot ', name, ' were not of type int, float, or list. ',
        'Check your input in the .yaml file. The TTFD at the ',
        'actual monitoring locations will not be plotted.']).format(loc_type)

    return debug_msg


def gen_monitoring_loc_len_warning(plume_type, name):
    """
    Function generating the warning message used in get_ttfd_plot_yaml_input()
    for the lengths monitoring coordx, coordy, and coordz values.
    """
    if plume_type == 'CarbonateAquifer':
        debug_msg = ''.join([
            'The coordx and coordy lists provided ',
            'for monitoring locations in the TTFD plot ', name,
            ' did not have the same length. Check your input in ',
            'the .yaml file. The TTFD at the actual monitoring ',
            'locations will not be plotted.'])
    else:
        debug_msg = ''.join([
            'The coordx, coordy, and coordz lists provided ',
            'for monitoring locations in the TTFD plot ', name,
            ' did not have the same length. Check your input in ',
            'the .yaml file. The TTFD at the actual monitoring ',
            'locations will not be plotted.'])

    return debug_msg


def get_aq_comp_lists_and_xy_grids(sm, yaml_data, name, aq_name_list):
    """
    Function that produces lists of aquifer components, their characteristics
    (e.g., x and y values and the index for the corresponding aquifer), reservoir
    components, and the injection locations. By going through these components,
    it also obtains and provides the minimum x and y values used in the
    simulation(s). These boundaries are then used to make the x_grid and y_grid
    values required for the function get_plume_timings().
    """
    # Define grid coordinates. These are overwritten within the loop below
    min_x_val = 9e99
    max_x_val = 0
    min_y_val = 9e99
    max_y_val = 0

    aq_components = []
    aq_component_types = []
    aq_component_indices = []
    aq_component_xvals = []
    aq_component_yvals = []

    # This is only used for Carbonate Aquifers
    aq_component_ithresh = []

    res_comp_injX = []
    res_comp_injY = []

    components = list(sm.component_models.values())
    for comp in components:
        if comp.class_type != 'Stratigraphy':
            comp_data = yaml_data[comp.name]
        else:
            continue

        if comp.class_type in TTFD_AQUIFER_COMPONENTS:
            for aqRef, aq_name in enumerate(aq_name_list):

                if aq_name in comp.name:
                    if 'AquiferName' in comp_data:
                        aq_name_from_yaml = comp_data['AquiferName']
                        # See if the aq. number has one digit or two, then get the number
                        try:
                            aq_number_from_yaml = int(aq_name_from_yaml[7:])
                        except ValueError:
                            err_msg = ''.join([
                                'The aquifer component {} has an entry of {} for ',
                                'the AquiferName field. This entry does not ',
                                'match the expected format of aquifer#, where ',
                                '# is an integer between 1 and 29. Check your ',
                                'input in the .yaml file.']).format(
                                    comp.name, comp_data['AquiferName'])
                            raise ValueError(err_msg) from None
                    else:
                        err_msg = ''.join([
                            'The aquifer component {} does not have an entry ',
                            'for AquiferName in the .yaml file (e.g., aquifer2). ',
                            'This entry is required for the TTFD plot type. ',
                            'Check your input.']).format(comp.name)
                        raise KeyError(err_msg)

                    if comp.class_type != 'CarbonateAquifer':
                        aq_components.append(comp)
                        aq_component_types.append(comp.class_type)

                        try:
                            # If reservoir data to wellbore components come from
                            # reservoir component get locations from it
                            res_comp = sm.component_models[yaml_data[yaml_data[
                                comp.name]['Connection']]['Connection']]

                            aq_component_xvals.append(res_comp.locX)
                            aq_component_yvals.append(res_comp.locY)
                        except KeyError:
                            # Get location from the wellbore component itself
                            well_comp = yaml_data[comp.name]['Connection']
                            aq_component_xvals.append(yaml_data[well_comp]['locX'])
                            aq_component_yvals.append(yaml_data[well_comp]['locY'])

                        aq_component_indices.append(aq_number_from_yaml - 1)

                    elif comp.class_type == 'CarbonateAquifer':
                        ithresh_val = iamcommons.get_parameter_val(comp, 'ithresh')

                        try:
                            carbaq_x_vals = float(comp_data['locX'])
                            carbaq_y_vals = float(comp_data['locY'])

                        except TypeError:
                            carbaq_x_vals = np.array(list(comp_data['locX']),
                                                     dtype = float)
                            carbaq_y_vals = np.array(list(comp_data['locY']),
                                                     dtype = float)

                        if isinstance(carbaq_x_vals, float):
                            aq_components.append(comp)
                            aq_component_types.append(comp_data['Type'])

                            aq_component_xvals.append(carbaq_x_vals)
                            aq_component_yvals.append(carbaq_y_vals)

                            aq_component_indices.append(aq_number_from_yaml - 1)

                            aq_component_ithresh.append(ithresh_val)

                        else:
                            for carbaq_x_value in carbaq_x_vals:
                                aq_components.append(comp)
                                aq_component_types.append(comp_data['Type'])

                                aq_component_indices.append(aqRef)

                                aq_component_xvals.append(carbaq_x_value)

                                # These lists needs to have the same length as
                                # the x and y lists. If the Carb. Aq. has a list
                                # of locX and locY values, I append the aq. #s
                                # and ithresh values here to keep the list
                                # lengths the same.
                                aq_component_indices.append(aq_number_from_yaml - 1)

                                aq_component_ithresh.append(ithresh_val)

                            for carbaq_y_value in carbaq_y_vals:
                                aq_components.append(comp)
                                aq_component_types.append(comp_data['Type'])

                                aq_component_yvals.append(carbaq_y_value)

        if comp.class_type in TTFD_RESERVOIR_COMPONENTS:
            # Get the injection sites
            if comp.class_type != 'LookupTableReservoir':
                res_comp_injX.append(comp.injX)
                res_comp_injY.append(comp.injY)

            else:
                res_comp_injX.append(None)
                res_comp_injY.append(None)

            TTFD_yaml_input_dict = get_ttfd_plot_yaml_input(yaml_data, name)

            if TTFD_yaml_input_dict['plot_injection_sites']:
                if comp.class_type != 'LookupTableReservoir':
                    if comp.injX < min_x_val:
                        min_x_val = comp.injX

                    if comp.injX > max_x_val:
                        max_x_val = comp.injX

                    if comp.injY < min_y_val:
                        min_y_val = comp.injY

                    if comp.injY > max_y_val:
                        max_y_val = comp.injY
                else:
                    InjectionCoordx = TTFD_yaml_input_dict['InjectionCoordx']
                    InjectionCoordy = TTFD_yaml_input_dict['InjectionCoordy']

                    if np.min(InjectionCoordx) < min_x_val:
                        min_x_val = np.min(InjectionCoordx)

                    if np.min(InjectionCoordx) > max_x_val:
                        max_x_val = np.min(InjectionCoordx)

                    if np.min(InjectionCoordy) < min_y_val:
                        min_y_val = np.min(InjectionCoordy)

                    if np.min(InjectionCoordy) > max_y_val:
                        max_y_val = np.min(InjectionCoordy)

            # Get the well / aquifer component locations from the res. comp.
            x_vals = comp.locX
            y_vals = comp.locY

            if np.min(x_vals) < min_x_val:
                min_x_val = np.min(x_vals)

            if np.max(x_vals) > max_x_val:
                max_x_val = np.max(x_vals)

            if np.min(y_vals) < min_y_val:
                min_y_val = np.min(y_vals)

            if np.max(y_vals) > max_y_val:
                max_y_val = np.max(y_vals)

    x_range = [min_x_val, max_x_val]
    y_range = [min_y_val, max_y_val]

    return aq_components, aq_component_types, aq_component_indices, \
        aq_component_xvals, aq_component_yvals, aq_component_ithresh, \
            res_comp_injX, res_comp_injY, x_range, y_range


def make_xandy_grids(x_range, y_range, x_grid_spacing=None,
                     y_grid_spacing=None, EnforceGridXandYLims=False,
                     gridXLims=None, gridYLims=None,
                     monitoringCoordX=None, monitoringCoordY=None,
                     min_x_range=200, min_y_range=200, use_buffer=True,
                     min_grid_spacing=0.01):
    """
    Function that creates the x and y grids used by get_plume_timings() and other
    functions.
    """

    min_x_val = x_range[0]
    max_x_val = x_range[1]

    min_y_val = y_range[0]
    max_y_val = y_range[1]

    # If the ranges in x and y values are not greater than min_x_range and
    # min_y_range, then increase the ranges. For example, even if a model uses
    # a single well the grid should still cover a larger area.
    if (max_x_val - min_x_val) < min_x_range:
        max_x_val += (min_x_range / 2)
        min_x_val -= (min_x_range / 2)

    if (max_y_val - min_y_val) < min_y_range:
        max_y_val += (min_y_range / 2)
        min_y_val -= (min_y_range / 2)

    # One might not want a point of interest (e.g., a well) to be right at the
    # boundary of the grid. The buffer makes it so the location is situated
    # within the grid.
    if use_buffer:
        x_buffer = (max_x_val - min_x_val) / 10
        y_buffer = (max_y_val - min_y_val) / 10
    else:
        x_buffer = 0
        y_buffer = 0

    if EnforceGridXandYLims:
        if x_grid_spacing:
            x_grid = np.arange(gridXLims[0], gridXLims[1] + x_grid_spacing,
                               x_grid_spacing)
        else:
            x_grid = np.linspace(gridXLims[0], gridXLims[1],
                                 101, endpoint=True)

        # y values of interest, as a "column" array
        if y_grid_spacing:
            y_grid = np.arange(gridYLims[0], gridYLims[1] + y_grid_spacing,
                               y_grid_spacing)[:, None]
        else:
            y_grid = np.linspace(gridYLims[0], gridYLims[1],
                                 101, endpoint=True)[:, None]
    else:
        if x_grid_spacing:
            x_grid = np.arange(min_x_val - x_grid_spacing,
                               max_x_val + (2 * x_grid_spacing),
                               x_grid_spacing)
        else:
            x_grid = np.linspace(min_x_val - x_buffer, max_x_val + x_buffer,
                                 101, endpoint=True)

        # y values of interest, as a "column" array
        if y_grid_spacing:
            y_grid = np.arange(min_y_val - y_grid_spacing,
                               max_y_val + (2 * y_grid_spacing),
                               y_grid_spacing)[:, None]
        else:
            y_grid = np.linspace(min_y_val - y_buffer, max_y_val + y_buffer,
                                 101, endpoint=True)[:, None]

    if monitoringCoordX is not None and monitoringCoordY is not None:
        if isinstance(monitoringCoordX, list) and isinstance(monitoringCoordY, list):
            for monitorCoord in monitoringCoordX:
                # If there is already an x value there, or close enough, then
                # do not add the monitoringCoordX.
                min_diff = np.min(np.abs(x_grid - monitorCoord))
                if min_diff > min_grid_spacing:
                    x_grid_list = x_grid.tolist()
                    x_grid_list.append(monitorCoord)
                    # This puts it in the correct order
                    x_grid = np.unique(x_grid_list)

            for monitorCoord in monitoringCoordY:
                min_diff = np.min(np.abs(y_grid - monitorCoord))
                if min_diff > min_grid_spacing:
                    y_grid_list = y_grid[:, 0].tolist()
                    y_grid_list.append(monitorCoord)
                    y_grid = np.unique(y_grid_list)[:, None]

        else:
            min_diff = np.min(np.abs(x_grid - monitoringCoordX))
            if min_diff > min_grid_spacing:
                x_grid_list = x_grid.tolist()
                x_grid_list.append(monitoringCoordX)
                x_grid = np.unique(x_grid_list)

            min_diff = np.min(np.abs(y_grid - monitoringCoordY))
            if min_diff > min_grid_spacing:
                y_grid_list = y_grid[:, 0].tolist()
                y_grid_list.append(monitoringCoordY)
                y_grid = np.unique(y_grid_list)[:, None]

    return x_grid, y_grid


def plot_plume_metric(plumeMetric, plotType, yaml_data, num_samples,
                      time_array, x_grid, y_grid, z_grid, plume_metric_abbrev,
                      aq_component_types, aq_component_xvals, aq_component_yvals,
                      ttfd_list=None, ttfd_x_list=None, ttfd_y_list=None,
                      name='TTFD_Figure1', analysis='lhs', realization=0, output_dir=None,
                      genfontsize=12, axislabelfontsize=14, titlefontsize=14,
                      labelfontweight='bold', colormap='plasma', figsize=(10, 8),
                      res_comp_injX=None, res_comp_injY=None, min_num_points=25,
                      var_type='noVariation'):
    """
    Function that creates plots of different plume metrics: earliest plume
    timings over space, earliest plume timings at monitoring locations
    (time to first-detection, TTFD), and plume occurence probabilities. Note
    that plume probabilities can only be made in LHS simulations and are
    calculated as the fraction of LHS realizations in which a plume occurred at
    each location within the domain.
    """
    TTFD_yaml_input_dict = get_ttfd_plot_yaml_input(yaml_data, name)

    figureDPI = TTFD_yaml_input_dict['FigureDPI']

    plot_injection_sites = TTFD_yaml_input_dict['plot_injection_sites']

    InjectionCoordx = TTFD_yaml_input_dict['InjectionCoordx']
    InjectionCoordy = TTFD_yaml_input_dict['InjectionCoordy']

    # If (a) LookupTableReservoir(s) is/are being used, InjectionCoordx and
    # InjectionCoordy must be provided in the .yaml file. Otherwise, the
    # injection locations are obtained through res_comp_injX and res_comp_injy.
    if plot_injection_sites and InjectionCoordx is None:
        InjectionCoordx = res_comp_injX
        InjectionCoordy = res_comp_injY

    monitoringTTFD = TTFD_yaml_input_dict['monitoringTTFD']
    monitoringCoordX = TTFD_yaml_input_dict['monitoringCoordX']
    monitoringCoordY = TTFD_yaml_input_dict['monitoringCoordY']
    monitoringCoordZ = TTFD_yaml_input_dict['monitoringCoordZ']

    EnforceXandYLims = TTFD_yaml_input_dict['EnforceXandYLims']
    xLims = TTFD_yaml_input_dict['xLims']
    yLims = TTFD_yaml_input_dict['yLims']

    plumeMetricMarkerSize = 5
    ttfdMarkerSize = 14

    if name is None:
        name = plotType

    # Check the file name for an extension like .tiff or .eps
    if '.' in name:
        name_extension = name[name.index('.'):None]
    else:
        name_extension = '.png'

    if plotType == 'plumeTimings':
        file_name = '{}_Plume_Timings{}{}'
        colorBarLabel = 'Time (years)'
    elif plotType == 'monitoringTTFD':
        file_name = '{}_Monitoring_TTFD{}{}'
        colorBarLabel = 'Time (years)'
    elif plotType == 'plumeProbabilities':
        file_name = '{}_Plume_Probabilities{}{}'
        colorBarLabel = 'Probability (%)'

    if plotType in ['plumeTimings', 'monitoringTTFD']:
        min_level = 0
        max_level = max(time_array) / 365.25
        interval = (max_level - min_level) / 100
        levels = np.arange(
            min_level, max_level + interval, interval)
        contourf_norm_factor = 365.25
    elif plotType == 'plumeProbabilities':
        levels = np.arange(MIN_PROBABILITY, 101, 1)
        contourf_norm_factor = 1

    X, Y = np.meshgrid(x_grid, y_grid)

    checkCarbAq = 'CarbonateAquifer' in aq_component_types

    # This is used to keep track of whether a monitoring sensor has been shown
    # in a legend yet. It only has to be shown once, not multiple times.
    sensor_lgnd_check = False

    # Figure formatting
    font = RC_FONT
    font['size'] = genfontsize
    plt.rc('font', **font)

    # Make the figure
    fig = plt.figure(name, figsize = figsize)
    cmap = plt.cm.get_cmap(colormap)

    ax = plt.gca()

    if plotType == 'monitoringTTFD':
        if not checkCarbAq:
            min_z = np.min(z_grid)
            max_z = np.max(z_grid)
        else:
            min_z = None
            max_z = None

        sensor_lgnd_check = plot_wells_inj_sites(
            ax, aq_component_xvals, aq_component_yvals, min_z,
            max_z, plot_injection_sites, InjectionCoordx,
            InjectionCoordy, monitoringTTFD, monitoringCoordX,
            monitoringCoordY, monitoringCoordZ, sensor_lgnd_check,
            genfontsize, checkCarbAq=checkCarbAq, labelWells=True)

        # Make the colorbar
        cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                            values=levels)
        cbar.set_label(
            colorBarLabel, rotation=90, fontsize=axislabelfontsize,
            fontweight=labelfontweight)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()

        min_metric_for_subplot = None
        max_metric_for_subplot = None

        if ttfd_list is not None:
            for listRef, ttfd_list_item in enumerate(ttfd_list):
                rgba = cmap(
                    ((ttfd_list_item / 365.25) - np.min(levels))
                    / (np.max(levels) - np.min(levels)))

                plt.plot(ttfd_x_list[listRef] / 1000, ttfd_y_list[listRef] / 1000,
                    color=rgba[0:3], marker='o', markerfacecolor=rgba[0:3],
                    linewidth=1, linestyle='none', markersize=ttfdMarkerSize,
                    zorder=1)

            if len(ttfd_list) > 0:
                min_metric_for_subplot = np.min(ttfd_list)
                max_metric_for_subplot = np.max(ttfd_list)

        plt.xlabel('Easting (km)', fontsize=axislabelfontsize,
                   fontweight=labelfontweight)
        plt.ylabel('Northing (km)', fontsize=axislabelfontsize,
                   fontweight=labelfontweight)

        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

        axtitle = gen_ax_title_for_subplot(
            plotType, np.min(z_grid), np.max(z_grid), min_metric_for_subplot,
            max_metric_for_subplot, checkCarbAq=checkCarbAq)

        ax.set_title(axtitle, fontweight=labelfontweight,
                     fontsize=genfontsize)

        if EnforceXandYLims:
            plt.xlim(xLims[0] / 1000.0, xLims[1] / 1000.0)
            plt.ylim(yLims[0] / 1000.0, yLims[1] / 1000.0)
        else:
            plt.xlim(np.min(x_grid) / 1000.0, np.max(x_grid) / 1000.0)
            plt.ylim(np.min(y_grid) / 1000.0, np.max(y_grid) / 1000.0)

    elif checkCarbAq:
        sensor_lgnd_check = plot_wells_inj_sites(
            ax, aq_component_xvals, aq_component_yvals, None,
            None, plot_injection_sites, InjectionCoordx,
            InjectionCoordy, monitoringTTFD, monitoringCoordX,
            monitoringCoordY, monitoringCoordZ, sensor_lgnd_check,
            genfontsize, checkCarbAq=checkCarbAq, labelWells=True)

        for _, spine in ax.spines.items():
            spine.set_zorder(3e6)

        _ = plt.contourf(X / 1000.0, Y / 1000.0, plumeMetric / contourf_norm_factor,
                          levels=levels, cmap=colormap, zorder=3)

        # Make the colorbar
        cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                            values=levels)
        cbar.set_label(colorBarLabel, rotation=90, fontsize=axislabelfontsize,
                       fontweight=labelfontweight)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()

        plt.xlabel('Easting (km)', fontsize=axislabelfontsize,
                   fontweight=labelfontweight)
        plt.ylabel('Northing (km)', fontsize=axislabelfontsize,
                   fontweight=labelfontweight)

        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

        min_metric_for_subplot = None
        max_metric_for_subplot = None

        checkValid = check_metric_validity(plumeMetric, plotType, valType='multiple')

        if checkValid:
            if plotType == 'plumeProbabilities':
                plumeMetricValid = plumeMetric[plumeMetric > MIN_PROBABILITY]
            elif plotType == 'plumeTimings':
                min_metric_for_subplot = np.min(plumeMetric)
                plumeMetricValid = plumeMetric[plumeMetric < THRESHOLD_TIME]

            # If there aren't enough points, contourf won't work
            if len(plumeMetricValid) < min_num_points:
                for plumeMetricVal in plumeMetricValid:
                    X_temp = X[plumeMetric[:, :] == plumeMetricVal]
                    Y_temp = Y[plumeMetric[:, :] == plumeMetricVal]

                    rgba = cmap(((plumeMetricVal / contourf_norm_factor
                                  ) - np.min(levels)) / (
                                      np.max(levels) - np.min(levels)))

                    plt.plot(X_temp / 1000, Y_temp / 1000, color=rgba[0:3],
                             marker='o', markerfacecolor=rgba[0:3],
                             linewidth=1, linestyle='none',
                             markersize=plumeMetricMarkerSize, zorder=1)

            if plotType == 'plumeProbabilities':
                min_metric_for_subplot = np.min(plumeMetric[
                    plumeMetric > MIN_PROBABILITY])
                max_metric_for_subplot = np.max(plumeMetric)
            elif plotType == 'plumeTimings':
                min_metric_for_subplot = np.min(plumeMetric)
                max_metric_for_subplot = np.max(plumeMetric[
                    plumeMetric < THRESHOLD_TIME])

        axtitle = gen_ax_title_for_subplot(
                plotType, None, None, min_metric_for_subplot,
                max_metric_for_subplot, checkCarbAq=checkCarbAq)

        ax.set_title(axtitle, fontweight=labelfontweight,
                     fontsize=genfontsize)

        if EnforceXandYLims:
            plt.xlim(xLims[0] / 1000.0, xLims[1] / 1000.0)
            plt.ylim(yLims[0] / 1000.0, yLims[1] / 1000.0)
        else:
            plt.xlim(np.min(x_grid) / 1000.0, np.max(x_grid) / 1000.0)
            plt.ylim(np.min(y_grid) / 1000.0, np.max(y_grid) / 1000.0)

    else:
        # These are the thresholds used to switch between the 4 subplots.
        min_z_subplot1 = np.min(z_grid)
        min_z_subplot2 = ((np.min(z_grid) - np.max(z_grid)) * 0.75) + (
            np.max(z_grid))
        min_z_subplot3 = ((np.min(z_grid) - np.max(z_grid)) * 0.5) + (
            np.max(z_grid))
        min_z_subplot4 = ((np.min(z_grid) - np.max(z_grid)) * 0.25) + (
            np.max(z_grid))
        max_z_subplot4 = np.max(z_grid)

        subplot_min_z_vals = [min_z_subplot1, min_z_subplot2, min_z_subplot3, min_z_subplot4]
        subplot_max_z_vals = [min_z_subplot2, min_z_subplot3, min_z_subplot4, max_z_subplot4]

        subplot_use_x_label = [False, False, True, True]
        subplot_use_y_label = [True, False, True, False]

        subplot_label_wells = [True, False, False, False]

        for subplotRef, (subplot_min_z_val, subplot_max_z_val)  in enumerate(
                zip(subplot_min_z_vals, subplot_max_z_vals)):
            # This is used to make certain results plot on top of other results
            # (lowest plume timings or highest plume probabilities).
            zorder_val = 1e6

            min_metric_for_subplot = None
            max_metric_for_subplot = None
            mean_metric_for_subplot = None

            ax = plt.subplot(2, 2, subplotRef + 1)

            sensor_lgnd_check = plot_wells_inj_sites(
                ax, aq_component_xvals, aq_component_yvals,
                subplot_min_z_val, subplot_max_z_val,
                plot_injection_sites, InjectionCoordx, InjectionCoordy,
                monitoringTTFD, monitoringCoordX, monitoringCoordY,
                monitoringCoordZ, sensor_lgnd_check, genfontsize,
                labelWells=subplot_label_wells[subplotRef])

            if subplot_use_x_label[subplotRef]:
                plt.xlabel('Easting (km)', fontsize=axislabelfontsize,
                           fontweight=labelfontweight)

            if subplot_use_y_label[subplotRef]:
                plt.ylabel('Northing (km)', fontsize=axislabelfontsize,
                           fontweight=labelfontweight)

            # Makes sure the outer edge of the plot stays on top of other elements
            for _, spine in ax.spines.items():
                spine.set_zorder(3e6)

            # Make the colorbar
            cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                                values=levels)
            cbar.set_label(
                colorBarLabel, rotation=90, fontsize=axislabelfontsize,
                fontweight=labelfontweight)
            tick_locator = ticker.MaxNLocator(nbins=5)
            cbar.locator = tick_locator
            cbar.update_ticks()

            if EnforceXandYLims:
                ax.set_xlim(xLims[0] / 1000.0, xLims[1] / 1000.0)
                ax.set_ylim(yLims[0] / 1000.0, yLims[1] / 1000.0)
            else:
                ax.set_xlim(np.min(x_grid) / 1000.0, np.max(x_grid) / 1000.0)
                ax.set_ylim(np.min(y_grid) / 1000.0, np.max(y_grid) / 1000.0)

            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
            ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

            if var_type == 'noVariation':
                zIndices = np.array(range(len(z_grid)))
                zMask = np.ma.masked_inside(
                    z_grid[:, 0, 0], subplot_min_z_val, subplot_max_z_val)
                zIndices = zIndices[zMask.mask]
            elif var_type == 'strikeAndDip':
                zIndices = range(len(z_grid))

            # Loop through the z_grid values and plume metrics
            for zRef in zIndices:
                z_temp = z_grid[zRef, :, :]

                # Get the plume metric and z values for the current zRef
                plumeMetric_temp = plumeMetric[zRef, :, :].copy()

                # Use only the metrics with z values in the current range. Each
                # layer of z values (z_grid[zRef, :, :]) is flat-lying / the same
                # when spatially uniform stratigraphy is used, but each layer dips
                # when dipping stratigraphy is used. This function excludes results
                # that fall outside the current depth range.
                if var_type == 'strikeAndDip':
                    plumeMetric_temp = clip_results_outside_of_z_range(
                        plumeMetric_temp, z_temp, plotType,
                        subplot_min_z_val, subplot_max_z_val)

                # These are used to make the lower plumeTimings or higher
                # plumeProbabilities plot on top of the other results. They are
                # also used the first time any valid results occur. Using only
                # the minimum or maximum values does not consistently put the
                # appropriate layer on top, so these are adjusted based on mean values.
                checkMinPlumeTimings = False
                checkMaxPlumeProbabilities = False

                checkValid = check_metric_validity(plumeMetric_temp, plotType,
                                                   valType='multiple')

                if checkValid:
                    if plotType == 'plumeTimings':
                        plumeMetricValid = plumeMetric_temp[
                            plumeMetric_temp < THRESHOLD_TIME]

                        if not mean_metric_for_subplot:
                            checkMinPlumeTimings = True
                        elif np.mean(plumeMetricValid) < mean_metric_for_subplot:
                            checkMinPlumeTimings = True

                    elif plotType == 'plumeProbabilities':
                        plumeMetricValid = plumeMetric_temp[
                            plumeMetric_temp > MIN_PROBABILITY]

                        if not mean_metric_for_subplot:
                            checkMaxPlumeProbabilities = True
                        elif np.mean(plumeMetricValid) > mean_metric_for_subplot:
                            checkMaxPlumeProbabilities = True

                    if not min_metric_for_subplot:
                        min_metric_for_subplot = np.min(plumeMetricValid)
                    elif np.min(plumeMetricValid) < min_metric_for_subplot:
                        min_metric_for_subplot = np.min(plumeMetricValid)

                    if not max_metric_for_subplot:
                        max_metric_for_subplot = np.max(plumeMetricValid)
                    elif np.max(plumeMetricValid) > max_metric_for_subplot:
                        max_metric_for_subplot = np.max(plumeMetricValid)

                    if not mean_metric_for_subplot:
                        mean_metric_for_subplot = np.mean(plumeMetricValid)
                    elif np.max(plumeMetricValid) > max_metric_for_subplot:
                        mean_metric_for_subplot = np.mean(plumeMetricValid)

                    # Plot the results with contourf
                    if checkMinPlumeTimings or checkMaxPlumeProbabilities:
                        zorder_val += 1
                        # Make the lowest plumeTiming or highest plumeProbability
                        # plot on top of other values.
                        ax.contourf(X / 1000.0, Y / 1000.0,
                                    plumeMetric_temp / contourf_norm_factor,
                                    levels=levels, cmap=colormap,
                                    zorder=zorder_val)
                    else:
                        ax.contourf(X / 1000.0, Y / 1000.0,
                                    plumeMetric_temp / contourf_norm_factor,
                                    levels=levels, cmap=colormap, zorder=2)

                    # If there are very few points with results, the contourf
                    # function won't work and you need to plot the points manually.
                    if len(plumeMetricValid) < min_num_points:
                        cmap = plt.cm.get_cmap(colormap)
                        for plumeMetricVal in plumeMetricValid:
                            X_temp = X[plumeMetric_temp == plumeMetricVal]
                            Y_temp = Y[plumeMetric_temp == plumeMetricVal]

                            rgba = cmap(((plumeMetricVal / contourf_norm_factor
                                          ) - np.min(levels)) / (
                                              np.max(levels) - np.min(levels)))

                            ax.plot(X_temp / 1000, Y_temp / 1000, color=rgba[0:3],
                                    marker='o', markerfacecolor=rgba[0:3],
                                    linewidth=1, linestyle='none',
                                    markersize=plumeMetricMarkerSize, zorder=1)

            # End of the z loop - now that the min and max metric values are
            # known, add the subplot title
            axtitle = gen_ax_title_for_subplot(
                plotType, subplot_min_z_val, subplot_max_z_val,
                min_metric_for_subplot, max_metric_for_subplot, checkCarbAq=checkCarbAq)

            ax.set_title(axtitle, fontweight=labelfontweight, fontsize=genfontsize)

    output_title = TITLE_OPTIONS[plotType].get(plume_metric_abbrev, plume_metric_abbrev)

    if checkCarbAq:
        output_title = output_title.format('')
    else:
        if plotType == 'plumeTimings':
            title_addition = ',\nLayers with Lower Times Shown Above Other Layers'
        elif plotType == 'plumeProbabilities':
            title_addition = ',\nLayers with Higher Probabilities Shown Above Other Layers'
        elif plotType == 'monitoringTTFD':
            title_addition = ''
        output_title = output_title.format(title_addition)

    if plotType == 'plumeProbabilities':
        output_title += ',\nNumber of LHS Samples: ' + str(num_samples)

    elif analysis in ['lhs', 'parstudy']:
        if plotType == 'monitoringTTFD':
            output_title += ',\nRealization {}'.format(str(realization))
        else:
            output_title += ', Realization {}'.format(str(realization))

    plt.suptitle(output_title, fontweight=labelfontweight,
                 fontsize=titlefontsize)

    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9,
                        top=0.8, wspace=0.1, hspace=0.1)
    plt.tight_layout()

    if output_dir:
        if plotType == 'plumeProbabilities' or analysis == 'forward':
            file_name = file_name.format(
                plume_metric_abbrev, '', name_extension)

        elif analysis in ['lhs', 'parstudy']:
            file_name = file_name.format(
                plume_metric_abbrev, '_Realization_{}'.format(str(realization)),
                name_extension)

        plt.savefig(os.sep.join([output_dir, file_name.format(
            plume_metric_abbrev, name_extension)]), dpi=figureDPI)
        plt.close()
    else:
        plt.show()


def gen_ax_title_for_subplot(plotType, min_z_value, max_z_value,
                             min_metric_for_subplot, max_metric_for_subplot,
                             checkCarbAq=False):
    """
    Function that produces the subplot titles used in plot_plume_metric(). If
    a CarbonateAquifer is used, a single subplot is used. Otherwise, there are
    four subplots that each cover 1/4th of the depth range in z_grid.
    """
    # titleRef is used for AX_TITLE_PARTIAL(_2). For titleRef = 0, there was no
    # suitable value. For titleRef = 1, there was a range in suitable values.
    # For titleRef = 2, there was only one suitable value.
    if plotType in ['plumeTimings', 'monitoringTTFD']:
        if min_metric_for_subplot is None:
            titleRef = 0
            axTitleAddPt1 = ''
            axTitleAddPt3 = ''
        elif min_metric_for_subplot < THRESHOLD_TIME:
            if min_metric_for_subplot == max_metric_for_subplot:
                titleRef = 2
                axTitleAddPt1 = min_metric_for_subplot / 365.25
                axTitleAddPt3 = ''
            else:
                titleRef = 1
                axTitleAddPt1 = min_metric_for_subplot / 365.25
                axTitleAddPt3 = max_metric_for_subplot / 365.25
        else:
            titleRef = 0
            axTitleAddPt1 = ''
            axTitleAddPt3 = ''

        if min_metric_for_subplot is not None:
            if min_metric_for_subplot > 1:
                # year{s}, plural
                axTitleAddPt2 = 's'
        else:
            axTitleAddPt2 = ''

    elif plotType == 'plumeProbabilities':
        if min_metric_for_subplot is None:
            titleRef = 0
            axTitleAddPt1 = ''
            axTitleAddPt2 = ''
            axTitleAddPt3 = ''
        elif min_metric_for_subplot > MIN_PROBABILITY:
            if min_metric_for_subplot == max_metric_for_subplot:
                titleRef = 2
                axTitleAddPt1 = min_metric_for_subplot
                axTitleAddPt3 = ''
            else:
                titleRef = 1
                axTitleAddPt1 = min_metric_for_subplot
                axTitleAddPt3 = max_metric_for_subplot
        else:
            titleRef = 0
            axTitleAddPt1 = ''
            axTitleAddPt3 = ''

        axTitleAddPt2 = ''

    partialTitle = AX_TITLE_PARTIAL[
        titleRef].get(plotType, plotType)
    partialTitlePt2 = AX_TITLE_PARTIAL_2[
        titleRef].get(plotType, plotType)

    if not checkCarbAq:
        axTitle = ''.join([
            'Depth Range: {:.0f} m to '.format(min_z_value),
            '{:.0f} m'.format(max_z_value), ',\n{}'.format(partialTitle),
            partialTitlePt2.format(axTitleAddPt1, axTitleAddPt2, axTitleAddPt3)])
    else:
        axTitle = ''.join([
            '{}'.format(partialTitle), partialTitlePt2.format(
                axTitleAddPt1, axTitleAddPt2, axTitleAddPt3)])

    return axTitle


def plot_wells_inj_sites(ax, aq_component_xvals, aq_component_yvals, min_z,
                         max_z, plot_injection_sites, InjectionCoordx,
                         InjectionCoordy, monitoringTTFD, monitoringCoordX,
                         monitoringCoordY, monitoringCoordZ, sensor_lgnd_check,
                         genfontsize, checkCarbAq=False, labelWells=False):
    """
    Function that plots the wells, monitoring locations, and injection
    locations used in the simulation.
    """
    change_sensor_lgnd_check = False

    ax.set_facecolor(BACKGROUND_COLOR)

    if plot_injection_sites:
        if isinstance(InjectionCoordx, float):
            if labelWells:
                plt.plot(InjectionCoordx / 1000, InjectionCoordy / 1000,
                         marker='s', color='k', linestyle='none',
                         markeredgewidth=2, markersize=8,
                         markerfacecolor='none', zorder=1e7,
                         label='Injection Site')
            else:
                plt.plot(InjectionCoordx / 1000, InjectionCoordy / 1000,
                         marker='s', color='k', linestyle='none',
                         markeredgewidth=2, markersize=8,
                         markerfacecolor='none', zorder=1e7)
        else:
            for injRef, (icoordX, icoordY) in enumerate(
                    zip(InjectionCoordx, InjectionCoordy)):
                if injRef == 0 and labelWells:
                    plt.plot(icoordX / 1000, icoordY / 1000,
                             marker='s', color='k', linestyle='none',
                             markeredgewidth=2, markersize=8,
                             markerfacecolor='none', zorder=1e7,
                             label='Injection Site')
                else:
                    plt.plot(icoordX / 1000, icoordY / 1000,
                             marker='s', color='k', linestyle='none',
                             markeredgewidth=2, markersize=8,
                             markerfacecolor='none', zorder=1e7)

    if monitoringTTFD and not checkCarbAq:
        monitoringCoordX_temp, monitoringCoordY_temp, \
            monitoringCoordZ_temp = get_monitors_within_z_range(
                monitoringCoordX, monitoringCoordY, monitoringCoordZ,
                min_z, max_z)

        if len(monitoringCoordZ_temp) > 0:
            if isinstance(monitoringCoordX_temp, list):
                for monitorRef, (mcoordX, mcoordY) in enumerate(
                        zip(monitoringCoordX_temp, monitoringCoordY_temp)):
                    if monitorRef == 0 and not sensor_lgnd_check:
                        plt.plot(mcoordX / 1000, mcoordY / 1000,
                                 marker='^', color='k',
                                 linestyle='none', markeredgewidth=1.5,
                                 markersize=10, markerfacecolor='none',
                                 label='Sensor', zorder=1.5e7)
                        change_sensor_lgnd_check = True
                    else:
                        plt.plot(mcoordX / 1000, mcoordY / 1000,
                                 marker='^', color='k',
                                 linestyle='none', markeredgewidth=1.5,
                                 markersize=10,  markerfacecolor='none',
                                 zorder=1.5e7)
            elif not sensor_lgnd_check:
                plt.plot(monitoringCoordX_temp / 1000,
                         monitoringCoordY_temp / 1000,
                         marker='^', color='k', linestyle='none',
                         markeredgewidth=1.5, markersize=10,
                         markerfacecolor='none',
                         label='Sensor', zorder=1.5e7)
                change_sensor_lgnd_check = True
            else:
                plt.plot(monitoringCoordX_temp / 1000,
                         monitoringCoordY_temp / 1000,
                         marker='^', color='k', linestyle='none',
                         markeredgewidth=1.5, markersize=10,
                         markerfacecolor='none', zorder=1.5e7)

    elif monitoringTTFD and checkCarbAq:
        if isinstance(monitoringCoordX, list):
            for monitorRef, (mcoordX, mcoordY) in enumerate(zip(monitoringCoordX,
                                                                monitoringCoordY)):
                if monitorRef == 0:
                    plt.plot(mcoordX / 1000, mcoordY / 1000,
                             marker='^', color='k',
                             linestyle='none',
                             markeredgewidth=1.5,
                             markersize=12,
                             markerfacecolor='none',
                             label='Sensor', zorder=1e7)
                else:
                    plt.plot(mcoordX / 1000, mcoordY / 1000,
                             marker='^', color='k',
                             linestyle='none',
                             markeredgewidth=1.5,
                             markersize=12,
                             markerfacecolor='none', zorder=1e7)
        else:
            plt.plot(monitoringCoordX / 1000, monitoringCoordY / 1000,
                     marker='^', color='k', linestyle='none',
                     markeredgewidth=1.5, markersize=12,
                     markerfacecolor='none', label='Sensor', zorder=1e7)

    if labelWells:
        plt.plot(np.array(aq_component_xvals) / 1000,
                  np.array(aq_component_yvals) / 1000,
                  linestyle='none', marker='o', color='k',
                  markeredgewidth=1.5, markersize=6,
                  markerfacecolor='none', zorder=1e7, label='Well')
    else:
        plt.plot(np.array(aq_component_xvals) / 1000,
                  np.array(aq_component_yvals) / 1000,
                  linestyle='none', marker='o', color='k',
                  markeredgewidth=1.5, markersize=6,
                  markerfacecolor='none', zorder=1e7)

    if labelWells or change_sensor_lgnd_check or checkCarbAq:
        ax.legend(fancybox=False, fontsize=genfontsize - 4,
                  edgecolor=[0, 0, 0], loc='upper left',
                  framealpha=0.5).set_zorder(2e7)
        if change_sensor_lgnd_check:
            # The sensor only needs to be displayed in a legend once
            sensor_lgnd_check = True

    return sensor_lgnd_check


def check_metric_validity(value, resultsType, valType='single'):
    """
    Function that checks if (a) value(s) is/are valid. Invalid times will be
    equal to MAX_TIME. To avoid the inclusion of very low probabilities, valid
    plumeProbabilities must be > MIN_PROBABILITY. When resultsType ==
    'monitoringTTFD', it checks if the ttfd_list provided as value has a length
    greater than 0.
    """
    if resultsType == 'plumeTimings':
        if valType == 'single':
            checkValid = (value < THRESHOLD_TIME)
        elif valType == 'multiple':
            checkValid = (np.min(value) < THRESHOLD_TIME)

    elif resultsType == 'monitoringTTFD':
        checkValid = len(value) > 0

    elif resultsType == 'plumeProbabilities':
        if valType == 'single':
            checkValid = (value > MIN_PROBABILITY)

        elif valType == 'multiple':
            checkValid = (np.max(value) > MIN_PROBABILITY)

    return checkValid


def clip_results_outside_of_z_range(plumeMetric_temp, z_temp, plotType,
                                    subplot_min_z, subplot_max_z):
    """
    When using a strike and dip, the z grid values depend on x and y. This function
    examines z_temp (the current z_grid[zRef, :, :] within the z loop) - for
    z values outside the current subplot range, it sets the corresponding
    plumeMetric_temp to a value that will be excluded from the plot
    (MIN_PROBABILITY or THRESHOLD_TIME).
    """
    mask = np.ma.masked_outside(z_temp, subplot_min_z, subplot_max_z)

    if plotType == 'plumeProbabilities':
        plumeMetric_temp[mask.mask] = MIN_PROBABILITY
    elif plotType == 'plumeTimings':
        plumeMetric_temp[mask.mask] = MAX_TIME

    return plumeMetric_temp


def save_results_to_csv(metric, x_grid, y_grid, z_grid, output_dir,
                        plume_metric_abbrev, resultsType, ttfd_list=None,
                        ttfd_x_list=None, ttfd_y_list=None, ttfd_z_list=None,
                        analysis='lhs', realization=0, num_samples=None,
                        var_type='noVariation', checkCarbAq=False):
    """
    Function that saves plume timing output to a series of .csv files. This
    function is used in cases where a CarbonateAquifer component is not used.
    """
    checkValidResults = True

    if resultsType == 'plumeProbabilities':
        resultsNormFactor = 1
    else:
        resultsNormFactor = 365.25

    if not os.path.exists(os.path.join(output_dir, 'csv_files')):
        os.mkdir(os.path.join(output_dir, 'csv_files'))

    if analysis in ['lhs', 'parstudy']:
        name_addition = '_Realization_' + str(realization)
    else:
        name_addition = ''

    if resultsType == 'plumeTimings':
        filename = plume_metric_abbrev + '_Plume_Timings{}.csv'.format(
            name_addition)
        metricLabel = 'Earliest_Plume_Timing_years'

    elif resultsType == 'monitoringTTFD':
        filename = plume_metric_abbrev + '_Monitoring_TTFD{}.csv'.format(
            name_addition)
        metricLabel = 'TTFD_years'

    elif resultsType == 'plumeProbabilities':
        filename = plume_metric_abbrev + '_Plume_Probabilities_' + str(
            num_samples) + '_LHS_Samples.csv'
        metricLabel = 'Plume_Presence_Probability_percent'

    filename = os.path.join(output_dir, 'csv_files', filename)

    if resultsType == 'monitoringTTFD':
        checkValidResults = check_metric_validity(ttfd_list, resultsType,
                                                  valType='multiple')
    else:
        checkValidResults = check_metric_validity(metric, resultsType,
                                                  valType='multiple')

    if checkCarbAq:
        if checkValidResults:
            first_row = [
                'Easting_x_m', 'Northing_y_m', metricLabel,
                'Note: Any locations that never had results are skipped. All '
                + 'x, and y values are saved in separate .csv files. Results '
                + 'were evaluated at all combinations of these x and y values.']
        else:
            first_row = ''.join([
                'In {}, no {} plumes occured within the area ',
                'considered (the x and y grids, saved in separate .csv files).',
                '{}'])
            if resultsType == 'plumeProbabilities':
                first_row = [first_row.format(
                    'the ' + str(num_samples) + ' LHS simulations assessed',
                    plume_metric_abbrev, ' All plume probabilities were 0%.')]
            else:
                if analysis != 'forward':
                    first_row = [first_row.format(
                        'realization ' + str(realization), plume_metric_abbrev, '')]
                else:
                    first_row = [first_row.format(
                        'this simulation', plume_metric_abbrev, '')]
    else:
        if checkValidResults:
            first_row = ['Easting_x_m', 'Northing_y_m', 'Depth_m', metricLabel,
                         'Note: Any locations that never had results are skipped. '
                         + 'All x, y, and z values are saved in separate .csv files '
                         + 'Results were evaluated at all combinations of these '
                         + 'x, y, and z values.']
        else:
            first_row = ''.join([
                'In {}, no {} plumes occured within the area considered (the ',
                'x and y grids, saved in separate .csv files).{}'])
            if resultsType == 'plumeProbabilities':
                first_row = [first_row.format(
                    'the ' + str(num_samples) + ' LHS simulations assessed',
                    plume_metric_abbrev, ' All plume probabilities were 0%.')]
            else:
                if analysis != 'forward':
                    first_row = [first_row.format(
                        'realization ' + str(realization), plume_metric_abbrev, '')]
                else:
                    first_row = [first_row.format(
                        'this simulation', plume_metric_abbrev, '')]

    if not checkValidResults:
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(first_row)

    else:
        if checkCarbAq:
            with open(filename, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(first_row)

                if resultsType == 'monitoringTTFD':
                    for listRef, ttfd_list_item in enumerate(ttfd_list):
                        row_temp = [ttfd_x_list[listRef], ttfd_y_list[listRef],
                                    ttfd_list_item / resultsNormFactor]
                        writer.writerow(row_temp)

                else:
                    for xRef in range(len(list(x_grid))):
                        for yRef in range(len(list(y_grid[:, 0]))):
                            checkValid = check_metric_validity(
                                metric[xRef, yRef],
                                resultsType, valType = 'single')

                            if checkValid:
                                row_temp = [x_grid[xRef], y_grid[yRef, 0],
                                            metric[xRef, yRef] / resultsNormFactor]
                                writer.writerow(row_temp)

        else:
            with open(filename, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(first_row)

                if resultsType == 'monitoringTTFD':
                    for listRef, ttfd_list_item in enumerate(ttfd_list):
                        row_temp = [ttfd_x_list[listRef], ttfd_y_list[listRef],
                                    ttfd_z_list[listRef],
                                    ttfd_list_item / resultsNormFactor]
                        writer.writerow(row_temp)

                else:
                    for xRef in range(len(list(x_grid))):
                        for yRef in range(len(list(y_grid[:, 0]))):
                            for zRef in range(len(z_grid)):
                                checkValid = check_metric_validity(
                                    metric[zRef, yRef, xRef], resultsType,
                                    valType = 'single')

                                if checkValid:
                                    if var_type == 'noVariation':
                                        row_temp = [x_grid[xRef], y_grid[yRef, 0],
                                                    z_grid[zRef, 0, 0],
                                                    metric[zRef, yRef, xRef]
                                                    / resultsNormFactor]

                                    elif var_type == 'strikeAndDip':
                                        row_temp = [x_grid[xRef], y_grid[yRef, 0],
                                                    z_grid[zRef, yRef, xRef],
                                                    metric[zRef, yRef, xRef]
                                                    / resultsNormFactor]

                                    writer.writerow(row_temp)


def save_grid_to_csv(x_grid, y_grid, z_grid, output_dir, var_type='noVariation',
                     checkCarbAq=False):
    """
    Saves the x, y, and z grid points to a .csv file.
    """
    if not os.path.exists(os.path.join(output_dir, 'csv_files')):
        os.mkdir(os.path.join(output_dir, 'csv_files'))

    if checkCarbAq:
        coordinates = ['x', 'y']
    else:
        coordinates = ['x', 'y', 'z']

    if var_type == 'noVariation' or checkCarbAq:
        for coord in coordinates:
            filename = 'TTFD_{}_grid_points.csv'.format(coord)
            filename = os.path.join(output_dir, 'csv_files', filename)

            first_row = ['{}_grid_points_m'.format(coord)]

            if coord == 'x':
                grid = x_grid[:]
            elif coord == 'y':
                grid = y_grid[:, 0]
            elif coord == 'z':
                grid = z_grid[:, 0, 0]

            with open(filename, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(first_row)

                for coordRef in range(len(list(grid))):
                    writer.writerow([grid[coordRef]])

    elif var_type == 'strikeAndDip':
        filename = 'TTFD_xyz_grid_points.csv'
        filename = os.path.join(output_dir, 'csv_files', filename)

        first_row = ['x_grid_points_m', 'y_grid_points_m', 'z_grid_points_m']

        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(first_row)

            for xRef in range(len(list(x_grid))):
                for yRef in range(len(list(y_grid[:, 0]))):
                    for zRef in range(len(z_grid)):
                        row_temp = [x_grid[xRef], y_grid[yRef, 0],
                                    z_grid[zRef, yRef, xRef]]

                        writer.writerow(row_temp)
