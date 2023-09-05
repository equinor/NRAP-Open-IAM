
import os
import sys
import logging
import csv
import numpy as np
import pandas as pd

SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(SOURCE_DIR)

try:
    from openiam import IAM_DIR
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

import openiam as iam
import openiam.cfi.strata as iam_strata
from openiam.visualize.area_of_review import CSV_FILE_NAME_TAGS as AOR_CSV_FILE_NAME_TAGS
from openiam.visualize.area_of_review import CSV_FILE_COLUMNS as AOR_CSV_FILE_COLUMNS
import openiam.visualize.area_of_review as AoR


AOR_CRIT_PRESSURE_COLUMN = 'Critical Pressure [MPa]'

WORKFLOW_OPTIONS = ['AoR']

WORKFLOW_PLOT_NAMES = {
    'AoR': [{'Pressure_Plot': 'pressure'},
            {'CO2_Saturation_Plot': 'CO2saturation'},
            {'Aq_CO2_Impact_Plot': '{AqCO2Metric}'},
            {'Aq_Brine_Metric_Plot': '{AqSaltMetric}'}],
    }

AOR_AQUIFER_COMPONENTS = ['FutureGen2Aquifer', 'FutureGen2AZMI',
                          'GenericAquifer', 'CarbonateAquifer',
                          'DeepAlluviumAquifer', 'DeepAlluviumAquiferML']

AOR_RESERVOIR_COMPONENT_OUTPUT = ['pressure', 'CO2saturation']

AOR_WELLBORE_COMPONENT_OUTPUT = ['brine_aquifer{}', 'CO2_aquifer{}']

AOR_AQUIFER_COMPONENT_OUTPUT = {
    'FutureGen2Aquifer': ['pH_volume', 'TDS_volume'],
    'FutureGen2AZMI': ['pH_volume', 'TDS_volume'],
    'GenericAquifer': ['Dissolved_CO2_volume', 'Dissolved_salt_volume'],
    'CarbonateAquifer': ['pH_volume', 'TDS_volume'],
    'DeepAlluviumAquifer': ['pH_volume', 'TDS_volume'],
    'DeepAlluviumAquiferML': ['pH_volume', 'TDS_volume']}

DEFAULT_LUTR_FILE_DIRECTORY = os.path.join(IAM_DIR, 'source', 'components',
                                           'reservoir', 'lookuptables',
                                           'FutureGen2', '1008_sims')

DEFAULT_LUTR_ENTRIES = {
    'FileDirectory': DEFAULT_LUTR_FILE_DIRECTORY,
    'TimeFile': 'time_points.csv',
    'ParameterFilename': 'parameters_and_filenames_trunc.csv',
    'Interpolation2D': True, 'Parameters': {'index': 1}}

YAML_INPUT_WARNING_MSG = ''.join([
    'The input provided for {} under the Workflow section of the .yaml file ',
    'was not of type {}. The default setting of {} will be used.'])


def iam_workflow_setup(yaml_data, strata):
    """
    This function reads the input provided in the Workflow section of a .yaml
    file. Depending on the options provided, the function can add to the yaml_data
    dictionary which is then returned by the function. For example, The AoR
    Workflow requires the inclusion of AoR plot entries using multiple metrics
    (e.g., pressure, CO2saturation, Dissolved_salt_volume, and Dissolved_CO2_volume).
    Setting up these plots and the required components can require an experienced
    user, so this function is designed to handle much of that effort.
    """

    if 'Type' in yaml_data['Workflow']:
        if yaml_data['Workflow']['Type'] in WORKFLOW_OPTIONS:
            workflow_type = yaml_data['Workflow']['Type']
        else:
            warning_msg = ''.join([
                'A Workflow section was included in the .yaml file, but the Workflow ',
                'Type provided (', yaml_data['Workflow']['Type'], ') was not ',
                'recognized as one of the available options (', WORKFLOW_OPTIONS,
                '). Therefore, the input in the Workflow section will not be used.'])
            logging.warning(warning_msg)
    else:
        info_msg = ''.join([
            'A Workflow section was included in the .yaml file, but the section ',
            'did not include a Type entry. Therefore, the input in the Workflow ',
            'section will not be used.'])
        logging.info(info_msg)

    if 'Options' not in yaml_data['Workflow']:
        info_msg = ''.join([
            'The Workflow section of the .yaml file did not include an ',
            'Options section. The default settings will be used.'])
        logging.info(info_msg)

        yaml_data['Workflow']['Options'] = {}

    if workflow_type == 'AoR':
        yaml_data = aor_workflow_setup(yaml_data, strata)

    return yaml_data


def aor_workflow_setup(yaml_data, strata):
    """
    Sets up the components and plots required for the AoR Workflow.
    """
    automation_input = get_automation_input(yaml_data)

    res_component_type = yaml_data['Workflow']['Options'].get(
        'ReservoirComponentType', 'LookupTableReservoir')

    well_component_type = yaml_data['Workflow']['Options'].get(
        'WellboreComponentType', 'OpenWellbore')

    aquifer_component_type = yaml_data['Workflow']['Options'].get(
        'AquiferComponentType', 'GenericAquifer')

    strata_dict = iam_strata.get_strata_info_from_component(strata)
    default_aquifer_name = 'aquifer{}'.format(str(strata_dict['numberOfShaleLayers'] - 1))

    aquifer_name = yaml_data['Workflow']['Options'].get(
        'AquiferName', default_aquifer_name)

    if 'Components' not in yaml_data['ModelParams']:
        yaml_data['ModelParams']['Components'] = []

    elif not isinstance(yaml_data['ModelParams']['Components'], list):
        comps = yaml_data['ModelParams']['Components']

        yaml_data['ModelParams']['Components'] = []
        for comp in comps:
            yaml_data['ModelParams']['Components'].append(comp)

    if 'CriticalPressureMPa' not in yaml_data['Workflow']['Options']:
        yaml_data['Workflow']['Options']['CriticalPressureMPa'] = 'Calculated'

    if automation_input['ResComp']:
        yaml_data, res_component_name = add_aor_res_component_entries(
            yaml_data, res_component_type)

        yaml_data['ModelParams']['Components'].append(res_component_name)

    if automation_input['WellComp']:
        yaml_data, well_component_name = add_aor_well_component_entries(
            yaml_data, well_component_type, aquifer_name, res_component_name)

        yaml_data['ModelParams']['Components'].append(well_component_name)

    if automation_input['AqComp']:
        yaml_data, aq_component_name = add_aor_aq_component_entries(
            yaml_data, aquifer_component_type, aquifer_name, well_component_name)

        yaml_data['ModelParams']['Components'].append(aq_component_name)

    if automation_input['Plots']:
        yaml_data = set_up_aor_plots_section(yaml_data, aquifer_component_type)

    return yaml_data


def workflow_analysis(yaml_data, sm, analysis):
    """
    After the simulation has completed in openiam_cf.py, this function is called
    if a Workflow section is included. This function then checks for and runs
    a specific Workflow analysis.
    """
    if 'Type' in yaml_data['Workflow']:
        # Add future Workflow analysis types here.
        if yaml_data['Workflow']['Type'] == 'AoR':
            aor_workflow_analysis(yaml_data, sm, analysis)


def aor_workflow_analysis(yaml_data, sm, analysis):
    """
    Evaluates the analysis of the data saved by each AoR plot entry and delineates
    an overall AoR that relects all metrics considered: pressure, CO2saturation,
    and both CO2 and brine impacts on the aquifer considered (i.e., plume volumes).
    """

    TimeList_option = False
    if 'TimeList' in yaml_data['Workflow']['Options']:
        if isinstance(yaml_data['Workflow']['Options']['TimeList'], list) or \
                yaml_data['Workflow']['Options']['TimeList'] == 'All':
            TimeList_option = True

            TimeList = yaml_data['Workflow']['Options']['TimeList']

            time = sm.time_array / 365.25
            if TimeList == 'All':
                time_index_list = range(len(time))
            else:
                time_index_list = AoR.get_t_indices(TimeList, time)

    if TimeList_option:
        for tIndex in time_index_list:
            All_x_points_km, All_y_points_km, AoR_point_included, pressure_included = \
                get_points_in_aor(yaml_data, sm=sm, time_index=tIndex)

            AoR.plot_aor_workflow_results(
                yaml_data, sm, All_x_points_km, All_y_points_km,
                AoR_point_included, time_index=tIndex, analysis=analysis)
    else:
        All_x_points_km, All_y_points_km, AoR_point_included, pressure_included = \
            get_points_in_aor(yaml_data, sm=sm)

        AoR.plot_aor_workflow_results(
            yaml_data, sm, All_x_points_km, All_y_points_km,
            AoR_point_included, analysis=analysis, pressure_included=pressure_included)


def get_points_in_aor(yaml_data, sm=None, time_index=None):
    """
    Checks the AoR plot type results for multiple metrics and determines an
    AoR that reflects all metrics considered. Pressure results are only considered
    if a critical pressure is provided.
    """
    output_dir = yaml_data['ModelParams']['OutputDirectory']

    aquifer_component_type = yaml_data['Workflow']['Options'].get(
        'AquiferComponentType', 'GenericAquifer')

    # Get the stratigraphy information from the .yaml file
    strata_var_info = iam_strata.get_strata_var_info_from_yaml(yaml_data)
    var_type = strata_var_info['var_type']

    AoR_included_x_km = []
    AoR_included_y_km = []

    if time_index is None:
        AoR_filename = 'AoR_{}.csv'
    else:
        AoR_filename = 'AoR_{}' + '_tIndex_{}.csv'.format(time_index)

    pressure_included = True
    # Reservoir pressures, only use if a critical pressure was provided
    if 'CriticalPressureMPa' in yaml_data['Workflow']['Options']:
        # Reservoir pressures
        pressure_file = AoR_filename.format(AOR_CSV_FILE_NAME_TAGS['pressure'])

        AoR_pressure_results = pd.read_csv(os.path.join(IAM_DIR, output_dir,
                                                        'csv_files', pressure_file))

        pressure_x_km = AoR_pressure_results['x (km)'].values
        pressure_y_km = AoR_pressure_results['y (km)'].values

        pressure_vals_MPa = AoR_pressure_results[AOR_CSV_FILE_COLUMNS['pressure']].values

        try:
            if var_type == 'noVariation':
                critPressureVal = np.min(AoR_pressure_results[AOR_CRIT_PRESSURE_COLUMN])
            elif var_type in ['strikeAndDip', 'LookupTable']:
                critPressureVal = AoR_pressure_results[AOR_CRIT_PRESSURE_COLUMN]
        except:
            if yaml_data['Workflow']['Options']['CriticalPressureMPa'] == 'Calculated':
                # returned in MPa
                critPressureVal = get_crit_pressure_aor_analysis(
                    len(pressure_x_km), yaml_data, sm)
            else:
                critPressureVal = float(yaml_data['Workflow']['Options']['CriticalPressureMPa'])

        if var_type == 'noVariation':
            pressure_x_km_AoR = pressure_x_km[pressure_vals_MPa >= critPressureVal]
            pressure_y_km_AoR = pressure_y_km[pressure_vals_MPa >= critPressureVal]

            if len(pressure_x_km_AoR) > 0:
                AoR_included_x_km += pressure_x_km_AoR.tolist()
                AoR_included_y_km += pressure_y_km_AoR.tolist()

        elif var_type in ['strikeAndDip', 'LookupTable']:
            for loc_ref in range(len(pressure_x_km)):
                if pressure_vals_MPa[loc_ref] >= critPressureVal[loc_ref]:
                    AoR_included_x_km.append(pressure_x_km[loc_ref])
                    AoR_included_y_km.append(pressure_y_km[loc_ref])

    else:
        pressure_x_km = None
        pressure_y_km = None

        pressure_included = False

    # Reservoir CO2 saturations
    CO2saturation_file = AoR_filename.format(AOR_CSV_FILE_NAME_TAGS['CO2saturation'])

    AoR_CO2saturation_results = pd.read_csv(os.path.join(
        IAM_DIR, output_dir, 'csv_files', CO2saturation_file))

    CO2saturation_x_km = AoR_CO2saturation_results['x (km)'].values
    CO2saturation_y_km = AoR_CO2saturation_results['y (km)'].values

    CO2saturation_vals = AoR_CO2saturation_results[AOR_CSV_FILE_COLUMNS['CO2saturation']].values

    CO2saturation_x_km_AoR = CO2saturation_x_km[CO2saturation_vals > 0]
    CO2saturation_y_km_AoR = CO2saturation_y_km[CO2saturation_vals > 0]

    if len(CO2saturation_x_km_AoR) > 0:
        AoR_included_x_km += CO2saturation_x_km_AoR.tolist()
        AoR_included_y_km += CO2saturation_y_km_AoR.tolist()

    # pH or dissolved CO2 contaminant plume volumes
    CO2impact_file = AoR_filename.format(AOR_CSV_FILE_NAME_TAGS[
        AOR_AQUIFER_COMPONENT_OUTPUT[aquifer_component_type][0]])

    AoR_CO2impact_results = pd.read_csv(os.path.join(IAM_DIR, output_dir,
                                                     'csv_files', CO2impact_file))

    CO2impact_x_km = AoR_CO2impact_results['x (km)'].values
    CO2impact_y_km = AoR_CO2impact_results['y (km)'].values

    CO2impact_vals_m3 = AoR_CO2impact_results[AOR_CSV_FILE_COLUMNS[
        AOR_AQUIFER_COMPONENT_OUTPUT[aquifer_component_type][0]]].values

    CO2impact_x_km_AoR = CO2impact_x_km[CO2impact_vals_m3 > 0]
    CO2impact_y_km_AoR = CO2impact_y_km[CO2impact_vals_m3 > 0]

    if len(CO2impact_x_km_AoR) > 0:
        AoR_included_x_km += CO2impact_x_km_AoR.tolist()
        AoR_included_y_km += CO2impact_y_km_AoR.tolist()

    # TDS or dissolved salt contaminant plume volumes
    brineimpact_file = AoR_filename.format(AOR_CSV_FILE_NAME_TAGS[
        AOR_AQUIFER_COMPONENT_OUTPUT[aquifer_component_type][1]])

    AoR_brineimpact_results = pd.read_csv(os.path.join(IAM_DIR, output_dir,
                                                       'csv_files', brineimpact_file))

    brineimpact_x_km = AoR_brineimpact_results['x (km)'].values
    brineimpact_y_km = AoR_brineimpact_results['y (km)'].values

    brineimpact_vals_m3 = AoR_brineimpact_results[AOR_CSV_FILE_COLUMNS[
        AOR_AQUIFER_COMPONENT_OUTPUT[aquifer_component_type][1]]].values

    brineimpact_x_km_AoR = brineimpact_x_km[brineimpact_vals_m3 > 0]
    brineimpact_y_km_AoR = brineimpact_y_km[brineimpact_vals_m3 > 0]

    if len(brineimpact_x_km_AoR) > 0:
        AoR_included_x_km += brineimpact_x_km_AoR.tolist()
        AoR_included_y_km += brineimpact_y_km_AoR.tolist()

    # Remove redundant entries
    AoR_included_x_km, AoR_included_y_km = remove_redundant_points(
        AoR_included_x_km, AoR_included_y_km)

    # Assemble the points that are not in the AoR.
    # To do that, evaluate all points considered. The x and y points for all
    # result types should be the same, this function just verifies that.
    All_x_points_km, All_y_points_km = verify_aor_points(
        pressure_x_km, pressure_y_km, CO2saturation_x_km, CO2saturation_y_km,
        CO2impact_x_km, CO2impact_y_km, brineimpact_x_km, brineimpact_y_km)

    AoR_point_included = np.zeros(len(All_x_points_km))

    for loc_ref, (x_loc, y_loc) in enumerate(zip(All_x_points_km, All_y_points_km)):
        AoR_point_temp = (x_loc, y_loc)

        point_included = False
        for (x_loc_in, y_loc_in) in zip(AoR_included_x_km, AoR_included_y_km):
            AoR_included_loc_temp = (x_loc_in, y_loc_in)

            if AoR_point_temp == AoR_included_loc_temp:
                point_included = True

        if point_included:
            AoR_point_included[loc_ref] = 1

    results_formatted = np.empty(((len(All_x_points_km) + 1), 3), dtype=list)

    results_formatted[0, 0] = 'Considered AoR Point, x (km)'
    results_formatted[0, 1] = 'Considered AoR Point, y (km)'
    results_formatted[0, 2] = 'Point Included in AoR'

    results_formatted[1:, 0] = All_x_points_km
    results_formatted[1:, 1] = All_y_points_km
    results_formatted[1:, 2] = AoR_point_included

    if time_index is not None:
        filename = os.path.join(output_dir, 'csv_files',
                                'AoR_Workflow_Output_tIndex{}.csv'.format(time_index))
    else:
        filename = os.path.join(output_dir, 'csv_files', 'AoR_Workflow_Output.csv')

    # Save the workflow outputs
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        for row_ref in range(len(All_x_points_km) + 1):
            writer.writerow(results_formatted[row_ref, :])

    f.close()

    return All_x_points_km, All_y_points_km, AoR_point_included, pressure_included


def remove_redundant_points(AoR_x_km, AoR_y_km):
    """
    This function takes lists of x and y points and returns the lists after
    removing redundant combinations of x and y.
    """
    AoR_x_km_edit = []
    AoR_y_km_edit = []

    AoR_locs = []
    for (xLoc, yLoc) in zip(AoR_x_km, AoR_y_km):
        AoR_loc_temp = (xLoc, yLoc)

        if AoR_loc_temp not in AoR_locs:
            AoR_x_km_edit.append(xLoc)
            AoR_y_km_edit.append(yLoc)

            AoR_locs.append(AoR_loc_temp)

    AoR_x_km = AoR_x_km_edit
    AoR_y_km = AoR_y_km_edit

    return AoR_x_km, AoR_y_km


def verify_aor_points(pressure_x_km, pressure_y_km, CO2saturation_x_km,
                      CO2saturation_y_km, CO2impact_x_km, CO2impact_y_km,
                      brineimpact_x_km, brineimpact_y_km):
    """
    The x and y locations used in each of the plot types should be the same:
    this function ensures that they are indeed the same. If not, it produces an error.
    """
    x_all_same = True
    y_all_same = True

    metrics_string = '{}CO2 saturations{} and aquifer plume volumes from CO2 and brine Leakage'
    if pressure_x_km is None and pressure_y_km is None:
        x_all_same = (CO2saturation_x_km.all() == CO2impact_x_km.all() and
                      CO2impact_x_km.all() == brineimpact_x_km.all())

        y_all_same = (CO2saturation_y_km.all() == CO2impact_y_km.all() and
                      CO2impact_y_km.all() == brineimpact_y_km.all())

        metrics_string = metrics_string.format('', '')
    else:
        x_all_same = (pressure_x_km.all() == CO2saturation_x_km.all() and
                      CO2saturation_x_km.all() == CO2impact_x_km.all() and
                      CO2impact_x_km.all() == brineimpact_x_km.all())

        y_all_same = (pressure_y_km.all() == CO2saturation_y_km.all() and
                      CO2saturation_y_km.all() == CO2impact_y_km.all() and
                      CO2impact_y_km.all() == brineimpact_y_km.all())

        metrics_string = metrics_string.format('pressures, ', ',')

    if x_all_same and y_all_same:
        All_x_points_km = brineimpact_x_km
        All_y_points_km = brineimpact_y_km
    else:
        error_msg = ''.join([
            'The x and y points from the AoR .csv files for ', metrics_string,
            'were not identical. The x and y points for each AoR plot type ',
            'should be the same, if they are from one simulation with fixed ',
            'wellbore locations. Check your input for the simulation. The AoR ',
            'Workflow cannot be finished.'])
        logging.error(error_msg)

    return All_x_points_km, All_y_points_km


def add_aor_res_component_entries(yaml_data, res_component_type):
    """
    This automatically sets up the reservoir component entry for the AoR
    workflow, if AutomateResCompSetup is set to True.
    """
    compName = res_component_type + '1'

    output_list = AOR_RESERVOIR_COMPONENT_OUTPUT

    yaml_data[compName] = {'Type': res_component_type, 'Outputs': output_list}

    if res_component_type == 'LookupTableReservoir':
        lutr_entries = DEFAULT_LUTR_ENTRIES

        if 'ReservoirOptions' in yaml_data['Workflow']['Options']:
            lutr_entries.update(
                yaml_data['Workflow']['Options']['ReservoirOptions'])

        for entry in list(lutr_entries.keys()):
            yaml_data[compName][entry] = lutr_entries[entry]

    if 'ReservoirOptions' in yaml_data['Workflow']['Options'] and \
            res_component_type != 'LookupTableReservoir':
        if 'InjectionWell' in yaml_data['Workflow']['Options']['ReservoirOptions']:
            injection_x = yaml_data['Workflow']['Options']['ReservoirOptions'][
                'InjectionWell']['coordx']
            injection_y = yaml_data['Workflow']['Options']['ReservoirOptions'][
                'InjectionWell']['coordy']

            yaml_data[compName]['InjectionWell'] = {'coordx': injection_x,
                                                    'coordy': injection_y}

        if 'Parameters' in yaml_data['Workflow']['Options']['ReservoirOptions']:
            par_names = get_parameter_names(res_component_type)

            first_time_check = True
            for par in par_names:
                if par in yaml_data['Workflow']['Options']['ReservoirOptions'][
                        'Parameters']:
                    if first_time_check:
                        yaml_data[compName]['Parameters'] = {}
                        first_time_check = False

                    yaml_data[compName]['Parameters'][par] = yaml_data['Workflow'][
                        'Options']['ReservoirOptions']['Parameters'][par]

    return yaml_data, compName


def add_aor_well_component_entries(yaml_data, well_component_type, aquifer_name,
                                   res_component_name):
    """
    This automatically sets up the component entries for the AoR workflow, if
    AutomateWellCompSetup is set to True.
    """
    compName = well_component_type + '1'

    aquiferNum = aquifer_name[7:]

    output_list = [val.format(aquiferNum) for val in AOR_WELLBORE_COMPONENT_OUTPUT]

    controls = {'critPressureApproach': True, 'enforceCritPressure': False}

    loc_data_check = False
    if 'WellboreOptions' in yaml_data['Workflow']['Options']:
        if 'Controls' in yaml_data['Workflow']['Options']['WellboreOptions']:
            controls.update(yaml_data['Workflow']['Options']['WellboreOptions']['Controls'])

        if 'Locations' in yaml_data['Workflow']['Options']['WellboreOptions']:
            loc_data = yaml_data['Workflow']['Options']['WellboreOptions']['Locations']
            loc_data_check = True

    if yaml_data[res_component_name]['Type'] == 'LookupTableReservoir' and not loc_data_check:
        file_dir = yaml_data[res_component_name]['FileDirectory']

        par_file_name = yaml_data[res_component_name]['ParameterFilename']

        index = int(yaml_data[res_component_name]['Parameters']['index'])

        data = pd.read_csv(os.path.join(IAM_DIR, file_dir, par_file_name))

        filename = data['filename'][index]

        file_name_dir = os.path.join(IAM_DIR, file_dir, filename)

        read_z_values = False

        data = pd.read_csv(file_name_dir)

        if 'z' in data:
            read_z_values = True

        del data

        # If no location data are given for the OpenWelbores, use the LUTR file itself by default
        thin_point_density = yaml_data['Workflow']['Options'].get(
            'thin_point_density', True)
        min_x_spacing = yaml_data['Workflow']['Options'].get(
            'min_x_spacing', 20000)
        min_y_spacing = yaml_data['Workflow']['Options'].get(
            'min_y_spacing', 20000)

        loc_data = {'file': file_name_dir, 'read_z_values': read_z_values,
                    'thin_point_density': thin_point_density,
                    'min_x_spacing': min_x_spacing, 'min_y_spacing': min_y_spacing}
        loc_data_check = True

    # If no wellbore locations were given and the simulation does not use a
    # LookupTableReservoir (from which to draw wellbore locations), use default values
    if not loc_data_check:
        xmin = -50000
        xmax = 50000
        xsize = 6
        ymin = -50000
        ymax = 50000
        ysize = 6
        loc_data = {'grid': {'xmin': xmin, 'xmax': xmax, 'xsize': xsize,
                             'ymin': ymin, 'ymax': ymax, 'ysize': ysize}}

    yaml_data[compName] = {'Type': well_component_type, 'LeakTo': aquifer_name,
                           'Connection': res_component_name, 'Locations': loc_data,
                           'Controls': controls, 'Outputs': output_list}

    if 'WellboreOptions' in yaml_data['Workflow']['Options']:
        if 'Parameters' in yaml_data['Workflow']['Options']['WellboreOptions']:
            par_names = get_parameter_names(well_component_type)

            first_time_check = True
            for par in par_names:
                if par in yaml_data['Workflow']['Options']['WellboreOptions'][
                        'Parameters']:
                    if first_time_check:
                        yaml_data[compName]['Parameters'] = {}
                        first_time_check = False

                    yaml_data[compName]['Parameters'][par] = yaml_data['Workflow'][
                        'Options']['WellboreOptions']['Parameters'][par]

            # critPressure will not be in the default_pars.keys(), which is used
            # in get_parameter_names(). So check that parameter separately.
            if 'critPressure' in yaml_data['Workflow']['Options']['WellboreOptions'][
                    'Parameters']:
                yaml_data[compName]['Parameters']['critPressure'] = yaml_data[
                    'Workflow']['Options']['WellboreOptions']['Parameters'][
                        'critPressure']

    return yaml_data, compName


def add_aor_aq_component_entries(yaml_data, aquifer_component_type, aquifer_name,
                                 well_component_name):
    """
    This automatically sets up the component entries for the AoR workflow, if
    AutomateAqCompSetup is set to True.
    """
    compName = aquifer_component_type + '1'

    output_list = AOR_AQUIFER_COMPONENT_OUTPUT.get(aquifer_component_type,
                                                   AOR_AQUIFER_COMPONENT_OUTPUT['GenericAquifer'])

    yaml_data[compName] = {'Type': aquifer_component_type, 'AquiferName': aquifer_name,
                           'Connection': well_component_name, 'Outputs': output_list}

    if 'AquiferOptions' in yaml_data['Workflow']['Options']:
        if 'Parameters' in yaml_data['Workflow']['Options']['AquiferOptions']:
            par_names = get_parameter_names(aquifer_component_type)

            first_time_check = True
            for par in par_names:
                if par in yaml_data['Workflow']['Options']['AquiferOptions'][
                        'Parameters']:
                    if first_time_check:
                        yaml_data[compName]['Parameters'] = {}
                        first_time_check = False

                    yaml_data[compName]['Parameters'][par] = yaml_data['Workflow'][
                        'Options']['AquiferOptions']['Parameters'][par]

    return yaml_data, compName


def set_up_aor_plots_section(yaml_data, aquifer_component_type):
    """
    This automatically sets up the Plots section of the .yaml file for the AoR
    Workflow, if AutomatePlotsSetup is set to True.
    """

    if 'Plots' not in yaml_data:
        yaml_data['Plots'] = {}

    figure_dpi = yaml_data['Workflow']['Options'].get('FigureDPI', 100)

    plot_types_to_add = WORKFLOW_PLOT_NAMES['AoR']

    for plotNum, plot in enumerate(plot_types_to_add):

        for plotName in list(plot.keys()):

            if plotNum == 2:
                metric = plot[plotName].format(
                    AqCO2Metric=AOR_AQUIFER_COMPONENT_OUTPUT[
                        aquifer_component_type][0])
            elif plotNum == 3:
                metric = plot[plotName].format(
                    AqSaltMetric=AOR_AQUIFER_COMPONENT_OUTPUT[
                        aquifer_component_type][1])
            else:
                metric = plot[plotName]

            TimeList = None
            if 'TimeList' in yaml_data['Workflow']['Options']:
                if isinstance(yaml_data['Workflow']['Options']['TimeList'], list) or \
                        yaml_data['Workflow']['Options']['TimeList'] == 'All':
                    TimeList = yaml_data['Workflow']['Options']['TimeList']

                else:
                    warning_msg = ''.join(['The TimeList entry was provided under ',
                                           'the Options subsection of the Workflow ',
                                           'section in the .yaml file, but the '
                                           'TimeList entry given was neither a list ',
                                           'nor the word "All". Therefore, this input ',
                                           'will be ignored.'])
                    logging.warning(warning_msg)
                    TimeList = None

            plot_injection_sites = False
            if 'PlotInjectionSites' in yaml_data['Workflow']['Options']:
                if isinstance(yaml_data['Workflow']['Options']['PlotInjectionSites'], bool):
                    plot_injection_sites = yaml_data['Workflow']['Options']['PlotInjectionSites']
                else:
                    warning_msg = YAML_INPUT_WARNING_MSG.format(
                        'PlotInjectionSites', 'boolean', 'False')
                    logging.warning(warning_msg)

            InjectionCoordx = None
            InjectionCoordy = None
            if 'InjectionCoordx' in yaml_data['Workflow']['Options'] and \
                    'InjectionCoordy' in yaml_data['Workflow']['Options']:
                InjectionCoordx = yaml_data['Workflow']['Options']['InjectionCoordx']
                InjectionCoordy = yaml_data['Workflow']['Options']['InjectionCoordy']

            yaml_data['Plots'][plotName] = {'AoR': [metric], 'SaveCSVFiles': True,
                                            'PlotInjectionSites': plot_injection_sites,
                                            'FigureDPI': figure_dpi}

            if TimeList is not None:
                yaml_data['Plots'][plotName]['TimeList'] = TimeList

            if InjectionCoordx is not None and InjectionCoordy is not None:
                yaml_data['Plots'][plotName]['InjectionCoordx'] = InjectionCoordx
                yaml_data['Plots'][plotName]['InjectionCoordy'] = InjectionCoordy

            if plot[plotName] == 'pressure':
                critPressureInput = yaml_data['Workflow']['Options'].get(
                    'CriticalPressureMPa', 'Calculated')

                yaml_data['Plots'][plotName][
                    'CriticalPressureMPa'] = critPressureInput

    return yaml_data


def get_automation_input(yaml_data):
    """
    Checks for input related to the automatic setup of a .yaml control file for
    a Workflow.
    """
    automation_input = {}
    for key in ['AutomateResCompSetup', 'AutomateWellCompSetup',
                'AutomateAqCompSetup', 'AutomatePlotsSetup']:
        inp_key = key[8:-5]
        automation_input[inp_key] = True
        if key in yaml_data['Workflow']['Options']:
            if isinstance(yaml_data['Workflow']['Options'][key], bool):
                automation_input[inp_key] = yaml_data['Workflow']['Options'][key]
            else:
                warning_msg = YAML_INPUT_WARNING_MSG.format(key, 'boolean', True)
                logging.warning(warning_msg)

    return automation_input


def get_parameter_names(component_type):
    """
    Function that returns the names of all parameters for a given component type.
    """
    num_years = 5
    time_array = 365.25 * np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    dummy_sm = iam.SystemModel(model_kwargs=sm_model_kwargs)

    comp_type = getattr(iam, component_type)

    dummy_comp = dummy_sm.add_component_model_object(comp_type(
        name='dummy_comp', parent=dummy_sm))

    par_names = dummy_comp.default_pars.keys()

    return par_names


def get_crit_pressure_aor_analysis(num_pressure_points, yaml_data, sm):
    """
    Function that calculates the critical pressures (in MPa) for OpenWellbore
    components. When using spatially variable stratigraphy, the critical pressures
    are returned in an array where the value in each row corresponds with the
    pressure_x_km and pressure_y_km values in the same row.
    """
    # Get the stratigrapy information from the .yaml file
    strata_var_info = iam_strata.get_strata_var_info_from_yaml(yaml_data)
    var_type = strata_var_info['var_type']

    if var_type == 'noVariation':
        critPressure = None
    elif var_type in ['strikeAndDip', 'LookupTable']:
        critPressure = np.zeros((num_pressure_points, 1))

    components = list(sm.component_models.values())

    for output_component in components:

        if isinstance(output_component, iam.OpenWellbore):
            if var_type == 'noVariation' and critPressure is None:
                # If using uniform stratigraphy, only do this once
                critPressureVal = AoR.get_crit_pressure(
                    output_component, sm=sm, yaml_data=yaml_data)

                critPressure = critPressureVal / 1.0e+6

            elif var_type in ['strikeAndDip', 'LookupTable']:
                critPressureVal = AoR.get_crit_pressure(
                    output_component, sm=sm, yaml_data=yaml_data)

                # Find the location index
                loc_ref = int(output_component.name.split('_')[-1])

                critPressure[loc_ref] = critPressureVal / 1.0e+6

    return critPressure
