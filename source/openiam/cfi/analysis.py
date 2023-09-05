import os
import sys
import logging
import collections

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

try:
    import openiam.visualize as iam_vis
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))


def process_analysis(yaml_data, model_data, sm, s, output_list, analysis, time_array):
    """
    Analyse control file setup related to the analysis setup to determine
    whether additional plots need to be produced.
    """
    if analysis != 'lhs':
        logging.warning('Sensitivity analysis only available for lhs sampling')
    else:
        analysis_dict = yaml_data['Analysis']

        if 'CorrelationCoeff' in analysis_dict:
            correlation_coefficients(analysis_dict, model_data, s, time_array)

        if 'SensitivityCoeff' in analysis_dict:
            sensitivity_coefficients(
                analysis_dict, model_data, sm, s, output_list, time_array)

        if 'MultiSensitivities' in analysis_dict:
            multi_sensitivities(
                analysis_dict, model_data, sm, s, output_list, time_array)

        if 'StackedSensitivityBars' in analysis_dict:
            stacked_sensitivity_bars(
                analysis_dict, model_data, sm, s, output_list, time_array)

        if 'TimeSeriesSensitivity' in analysis_dict:
            time_series_sensitivity(
                analysis_dict, model_data, sm, s, output_list, time_array)


def adjust_yaml_input(analysis_dict, analysis_type='CorrelationCoeff'):
    """
    This function takes input provided in the .yaml file and adjusts the naming
    conventions used for the input.

    For .yaml input to plots and analysis, we use a convention of no underscores
    and capitalizing the first letter of each word. The variables given to the
    functions, however, follow the PEP convention of using underscores and lower
    case characters.
    """
    not_bool_warning = ''.join([
        'The entry for {} provided for {} in the analysis section was not of ',
        'type boolean. The default settings will be used.'])

    if analysis_dict[analysis_type] is not None:
        if 'FigureDPI' in analysis_dict[analysis_type]:
            analysis_dict[analysis_type]['figure_dpi'] = analysis_dict[
                analysis_type]['FigureDPI']

            analysis_dict[analysis_type].pop('FigureDPI')

        if 'UseFormattedLabels' in analysis_dict[analysis_type]:
            if isinstance(analysis_dict[analysis_type]['SaveCSVFiles'], bool):
                analysis_dict[analysis_type]['use_formatted_labels'] = analysis_dict[
                    analysis_type]['UseFormattedLabels']
            else:
                warning_msg = not_bool_warning.format('UseFormattedLabels',
                                                      analysis_type)
                logging.debug(warning_msg)

            analysis_dict[analysis_type].pop('UseFormattedLabels')

        if 'OutputTimeIndex' in analysis_dict[analysis_type]:
            analysis_dict[analysis_type]['capture_point'] = analysis_dict[
                analysis_type]['OutputTimeIndex']

            analysis_dict[analysis_type].pop('OutputTimeIndex')

        if 'CorrelationType' in analysis_dict[analysis_type]:
            ctype = analysis_dict[analysis_type]['CorrelationType']

            if ctype in ['pearson', 'spearman']:
                analysis_dict[analysis_type]['ctype'] = analysis_dict[
                    analysis_type]['CorrelationType']
            else:
                warning_msg = ''.join([
                    'The CorrelationType entry provided for ', analysis_type,
                    ' in the Analysis section was not pearson or spearman. ',
                    'The default setting of pearson will be used.'])

            analysis_dict[analysis_type].pop('CorrelationType')

        if 'Excludes' in analysis_dict[analysis_type]:
            analysis_dict[analysis_type]['excludes'] = analysis_dict[
                analysis_type]['Excludes']

            analysis_dict[analysis_type].pop('Excludes')

        if 'NumberIncluded' in analysis_dict[analysis_type]:
            analysis_dict[analysis_type]['num_sensitivities'] = analysis_dict[
                analysis_type]['NumberIncluded']

            analysis_dict[analysis_type].pop('NumberIncluded')

        if 'SaveCSVFiles' in analysis_dict[analysis_type]:
            if isinstance(analysis_dict[analysis_type]['SaveCSVFiles'], bool):
                if not analysis_dict[analysis_type]['SaveCSVFiles']:
                    analysis_dict[analysis_type]['outfile'] = None
            else:
                warning_msg = not_bool_warning.format('SaveCSVFiles',
                                                      analysis_type)
                logging.debug(warning_msg)

            analysis_dict[analysis_type].pop('SaveCSVFiles')

    return analysis_dict


def correlation_coefficients(analysis_dict, model_data, s, time_array):
    """
    Analyze input and produce a correlation coefficients type of plot.
    """
    corrcoeff_dict = {'capture_point': len(time_array)-1,
                      'excludes': [],
                      'ctype': 'pearson',
                      'plot': True,
                      'printout': False,
                      'plotvals': True,
                      'figsize': (15, 15),
                      'title': 'Pearson Correlation Coefficients at {ct} years',
                      'xrotation': 90,
                      'savefig': 'correlation_coefficients_time_index_{cp}.png',
                      'outfile': 'correlation_coefficients_time_index_{cp}.csv',
                      'figure_dpi': 100}

    if isinstance(analysis_dict['CorrelationCoeff'], dict):
        analysis_dict = adjust_yaml_input(
            analysis_dict, analysis_type='CorrelationCoeff')

        corrcoeff_dict.update(analysis_dict['CorrelationCoeff'])

    if 'OutputDirectory' in model_data:
        corrcoeff_dict['savefig'] = os.path.join(
            model_data['OutputDirectory'], corrcoeff_dict['savefig'])

        corrcoeff_dict['outfile'] = os.path.join(
            model_data['OutputDirectory'], 'csv_files',
            corrcoeff_dict['outfile'])

        if not os.path.exists(os.path.join(
                model_data['OutputDirectory'], 'csv_files')):
            os.mkdir(os.path.join(
                model_data['OutputDirectory'], 'csv_files'))

    iam_vis.correlations_at_time(s, time_array, **corrcoeff_dict)


def sensitivity_coefficients(analysis_dict, model_data, sm, s, output_list, time_array):
    """
    Analyze input and produce a sensitivity coefficients type of plot.
    """
    if isinstance(analysis_dict['SensitivityCoeff'], dict):
        analysis_dict = adjust_yaml_input(
            analysis_dict, analysis_type='SensitivityCoeff')

        sens_dict = analysis_dict['SensitivityCoeff']
        obs = sens_dict.pop('Outputs')

    else:
        obs = analysis_dict['SensitivityCoeff']
        sens_dict = {}

    if isinstance(obs, str):
        obs = [obs]
    elif not isinstance(obs, list):
        err_msg = ''.join([
            'Setup of SensitivityCoeff analysis is not valid. Please ',
            'check your input file.'])
        raise TypeError(err_msg)

    res_obs = resolve_obs_names(obs, output_list)

    if 'capture_point' in sens_dict:
        capture_point = sens_dict.pop('capture_point')
    else:
        capture_point = len(time_array)-1

    if not isinstance(capture_point, collections.Iterable):
        capture_points = [capture_point]
    else:
        capture_points = capture_point

    for capture_point in capture_points:
        cp_obs = ['{ob}_{cp}'.format(ob=ob, cp=capture_point)
                  for ob in res_obs]

        sens_dict_full = {'title': None, #'{ob} Sensitivity Coefficients',
                          'ylabel': None,
                          'savefig': '{ob}_sensitivity.png',
                          'outfile': '{ob}_sensitivity.csv',
                          'obs': '{ob}',
                          'use_formatted_labels': False,
                          'figure_dpi': 100}

        sens_dict_full.update(sens_dict)

        for ob in cp_obs:
            if 'OutputDirectory' in model_data:
                sens_dict_full['savefig'] = os.path.join(
                    model_data['OutputDirectory'], sens_dict_full['savefig'])

                sens_dict_full['outfile'] = os.path.join(
                    model_data['OutputDirectory'], 'csv_files',
                    sens_dict_full['outfile'])

                if not os.path.exists(os.path.join(model_data['OutputDirectory'],
                                                   'csv_files')):
                    os.mkdir(os.path.join(
                        model_data['OutputDirectory'], 'csv_files'))

            sens_dict2 = {}

            for key, value in list(sens_dict_full.items()):
                if isinstance(value, str):
                    v = value.format(ob=ob, cp=capture_point)
                else:
                    v = value
                sens_dict2[key] = v

            sensitivities = s.rbd_fast(obsname=ob, print_to_console=False)
            iam_vis.simple_sensitivities_barplot(sensitivities, sm, **sens_dict2)


def multi_sensitivities(analysis_dict, model_data, sm, s, output_list, time_array):
    """
    Analyze input and produce a multi sensitivities type of plot.
    """
    if isinstance(analysis_dict['MultiSensitivities'], dict):
        analysis_dict = adjust_yaml_input(
            analysis_dict, analysis_type='MultiSensitivities')

        sens_dict = analysis_dict['MultiSensitivities']

        obs = sens_dict.pop('Outputs')
    else:
        obs = analysis_dict['MultiSensitivities']
        sens_dict = {}

    if isinstance(obs, str):
        obs = [obs]

    res_obs = resolve_obs_names(obs, output_list)

    if 'capture_point' in sens_dict:
        capture_point = sens_dict.pop('capture_point')
    else:
        capture_point = len(time_array)-1

    if not isinstance(capture_point, collections.Iterable):
        capture_points = [capture_point]
    else:
        capture_points = capture_point

    for cp in capture_points:
        cp_obs = ['{ob}_{cp}'.format(ob=ob, cp=cp)
                  for ob in res_obs]

        sens_dict_full = {
            'title': None,
            'ylabel': None,
            'savefig': 'multisensitivities_time_index_{}.png'.format(cp),
            'outfile': 'multisensitivities_time_index_{}.csv'.format(cp),
            'use_formatted_labels': False,
            'figure_dpi': 100}

        sens_dict_full.update(sens_dict)

        if 'OutputDirectory' in model_data:
            ob = '_'.join(obs)

            sens_dict_full['savefig'] = os.path.join(
                model_data['OutputDirectory'],
                sens_dict_full['savefig'].format(ob=ob, cp=cp))

            sens_dict_full['outfile'] = os.path.join(
                model_data['OutputDirectory'], 'csv_files',
                sens_dict_full['outfile'].format(ob=ob, cp=cp))

            if not os.path.exists(os.path.join(
                    model_data['OutputDirectory'], 'csv_files')):
                os.mkdir(os.path.join(
                    model_data['OutputDirectory'], 'csv_files'))

        iam_vis.multi_sensitivities_barplot(cp_obs, sm, s, **sens_dict_full)


def stacked_sensitivity_bars(analysis_dict, model_data, sm, s, output_list, time_array):
    """
    Analyze input and produce a stacked sensitivities bar plot.
    """
    if isinstance(analysis_dict['StackedSensitivityBars'], dict):
        analysis_dict = adjust_yaml_input(
            analysis_dict, analysis_type='StackedSensitivityBars')

        sens_dict = analysis_dict['StackedSensitivityBars']

        obs = sens_dict.pop('Outputs')
    else:
        obs = analysis_dict['StackedSensitivityBars']
        sens_dict = {}

    if isinstance(obs, str):
        obs = [obs]

    res_obs = resolve_obs_names(obs, output_list)

    if 'capture_point' in sens_dict:
        capture_point = sens_dict.pop('capture_point')
    else:
        capture_point = len(time_array)-1

    if not isinstance(capture_point, collections.Iterable):
        capture_points = [capture_point]
    else:
        capture_points = capture_point

    for cp in capture_points:
        cp_obs = ['{ob}_{cp}'.format(ob=ob, cp=cp) for ob in res_obs]

        sens_dict_full = {
            'title': None,
            'ylabel': None,
            'savefig': 'StackedSensitivities_Time_Index_{}.png'.format(cp),
            'outfile': 'StackedSensitivities_Time_Index_{}.csv'.format(cp),
            'figure_dpi': 100}

        sens_dict_full.update(sens_dict)

        sensitivities_list = []
        for ob in cp_obs:
            if 'OutputDirectory' in model_data:
                sens_dict_full['savefig'] = os.path.join(
                    model_data['OutputDirectory'],
                    sens_dict_full['savefig'].format(ob=ob, cp=cp))

                sens_dict_full['outfile'] = os.path.join(
                    model_data['OutputDirectory'], 'csv_files',
                    sens_dict_full['outfile'].format(ob=ob, cp=cp))

                if not os.path.exists(os.path.join(
                        model_data['OutputDirectory'], 'csv_files')):
                    os.mkdir(os.path.join(
                        model_data['OutputDirectory'], 'csv_files'))

            sens_dict2 = {}

            for key, value in list(sens_dict_full.items()):
                if isinstance(value, str):
                    v = value.format(ob=ob, cp=cp)
                else:
                    v = value
                sens_dict2[key] = v

            sensitivities = s.rbd_fast(obsname=ob, print_to_console=False)
            sensitivities_list.append(sensitivities)

        iam_vis.stacked_sensitivities_barplot(
            sensitivities_list, cp_obs, sm, **sens_dict_full)


def time_series_sensitivity(analysis_dict, model_data, sm, s, output_list, time_array):
    """
    Analyze input and produce a time series sensitivity type of plot.
    """
    if isinstance(analysis_dict['TimeSeriesSensitivity'], dict):
        analysis_dict = adjust_yaml_input(
            analysis_dict, analysis_type='TimeSeriesSensitivity')

        sens_dict = analysis_dict['TimeSeriesSensitivity']

        obs = sens_dict.pop('Outputs')

    else:
        obs = analysis_dict['TimeSeriesSensitivity']
        sens_dict = {}

    if isinstance(obs, str):
        obs = [obs]

    res_obs = resolve_obs_names(obs, output_list)

    if 'capture_point' in sens_dict:
        capture_point = sens_dict.pop('capture_point')
    else:
        capture_point = len(time_array)-1

    sens_dict_full = {
        'title': 'Time Sensitivity Coefficients for the\nOutput {ob} from {comp}',
        'ylabel': None,
        'num_sensitivities': 5,
        'savefig': '{ob}_time_series_sensitivity.png',
        'outfile': '{ob}_time_series_sensitivity.csv',
        'use_formatted_labels': False,
        'figure_dpi': 100}

    sens_dict_full.update(sens_dict)

    for ob in res_obs:
        if 'OutputDirectory' in model_data:
            sens_dict_full['savefig'] = os.path.join(
                model_data['OutputDirectory'], sens_dict_full['savefig'])

            sens_dict_full['outfile'] = os.path.join(
                model_data['OutputDirectory'], 'csv_files',
                sens_dict_full['outfile'])

            if not os.path.exists(os.path.join(
                    model_data['OutputDirectory'], 'csv_files')):
                os.mkdir(os.path.join(
                    model_data['OutputDirectory'], 'csv_files'))

        sens_dict = {}
        for key, value in list(sens_dict_full.items()):
            if isinstance(value, str) and not key == 'title':
                v = value.format(ob=ob)
            else:
                v = value

            sens_dict[key] = v

        iam_vis.time_series_sensitivities(
            ob, sm, s, time_array,
            capture_point=capture_point, **sens_dict)


def resolve_obs_names(obs, output_list):
    """
    Takes in base observation name and returns list of appended cm.obs for
    all component model names that have observations matching the base name.
    """
    resolved_names = []
    for ob in obs:
        for output_component in list(output_list.keys()):
            if ob in output_list[output_component]:
                resolved_names.append('.'.join([output_component.name, ob]))
    return resolved_names
