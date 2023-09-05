"""
This module contains set of methods assigning score to the well based on their
properties.

Last modified: June, 2023

@author: Greg Lackey, Veronika Vasylkivska
"""
import datetime
import numpy as np
import pandas as pd


SCORES = ['deviation_score', 'scvf_gm_score', 'suspension_status_score',
          'cement_top_score', 'casing_failure_score', 'age_score', 'depth_score',
          'type_score']


def well_data_classifier(wells):
    """
    Convert data into categories according to Duguid et al., 2017.
    """
    # Well age
    wells['best_date'] = wells['spud_date'].fillna(wells['first_production_date'])
    wells['best_date'] = wells['best_date'].fillna(wells['permit_date'])
    wells['best_date'] = pd.to_datetime(wells['best_date'])
    conditions = [
        (wells['best_date'] >= datetime.date(1997, 1, 1)),
        (wells['best_date'] < datetime.date(1997, 1, 1)) and (
            wells['best_date'] >= datetime.date(1986, 1, 1)),
        (wells['best_date'] < datetime.date(1986, 1, 1)) and (
            wells['best_date'] >= datetime.date(1975, 1, 1)),
        (wells['best_date'] < datetime.date(1975, 1, 1))]
    choices = ['spud after 1997',
               'spud between 1986 and 1997',
               'spud between 1975 and 1985',
               'spud pre-1975']
    wells['age_class'] = np.select(conditions, choices, default='unknown')

    # Well depth
    wells['best_depth'] = wells['tvd'].fillna(wells['prd_depth'])
    wells['best_depth'] = wells['best_depth'].fillna(wells['dtd'])
    wells['best_depth'] = wells['best_depth'].fillna(wells['ltd'])
    wells['best_depth'] = wells['best_depth']*0.3048
    conditions = [
        (wells['best_depth'] >= 2000.),
        (wells['best_depth'] < 2000.) and (wells['best_depth'] >= 1000.),
        (wells['best_depth'] < 1000.)]
    choices = ['>=2000 mkb',
               '>1000 mkb <2000 mkb',
               '<=1000 mkb']
    wells['depth_class'] = np.select(conditions, choices, default='unknown')

    # SCVF GM
    wells['scvf_gm_class'] = 'not tested'

    # Deviation
    conditions = [
        (wells['wellbore_configuration'] == 'H'),
        (wells['wellbore_configuration'] == 'D'),
        (wells['wellbore_configuration'] == 'V')]
    choices = ['horizontal', 'deviated', 'vertical']
    wells['dev_class'] = np.select(conditions, choices, default='unknown')

    # Cement top
    wells['best_surf_depth'] = wells['surf_depth'].fillna(0.)
    wells['best_prd_depth'] = wells['prd_depth'].fillna(wells['tvd'])
    wells['best_prd_depth'] = wells['best_prd_depth'].fillna(wells['dtd'])
    wells['best_prd_depth'] = wells['best_prd_depth'].fillna(wells['ltd'])
    wells['casing_len'] = wells['best_prd_depth'] - wells['surf_depth']
    wells['len_cem'] = wells['best_prd_depth'] - wells['prd_cement_top']
    wells['perc_cem'] = (wells['len_cem']/wells['casing_len'])*100.
    conditions = [
        (wells['perc_cem'] >= 99.),
        (wells['perc_cem'] < 99.) and (wells['perc_cem'] >= 50.),
        (wells['perc_cem'] < 50.) and (wells['perc_cem'] >= 15.),
        (wells['perc_cem'] < 15.)]
    choices = ['above surface casing shoe', 'moderate', 'low', 'very low']
    wells['cem_class'] = np.select(conditions, choices, default='missing data')

    # Well type
    conditions = [
        (wells['lateral_type'] == 'WI'),
        (wells['lateral_type'] == 'SWD'),
        (wells['lateral_type'] == 'GI'),
        (wells['lateral_type'] == 'GS'),
        (wells['lateral_type'] == 'OIL'),
        (wells['lateral_type'] == 'GAS'),
        (wells['lateral_type'] == 'GC'),
        (wells['lateral_type'] == 'CM'),
        (wells['lateral_type'] == 'SHG'),
        (wells['lateral_type'] == 'WW'),
        (wells['lateral_type'] == 'WS')]

    choices = ['h2o injection well', 'h2o disposal well',
               'co2 injection well', 'co2 injection well',
               'flowing oil well', 'medium risk gas well',
               'medium risk gas well', 'medium risk gas well',
               'medium risk gas well', 'pressure observation well',
               'pressure observation well']

    wells['type_class'] = np.select(
        conditions, choices, default='flowing oil well')

    return wells



def leak_scorer(wells, deviation=None, scvf_gm=None, suspension_status=None,
                cement_top=None, casing_failure=None, age=None, depth=None,
                well_type=None):
    """
    Assign leakage scores based on the wells properties.

    Method creates scores based on the categories from Duguid 2017.
    """
    inputs = [deviation,
              scvf_gm,
              suspension_status,
              cement_top,
              casing_failure,
              age,
              depth,
              well_type]

    # Initialize arrays of scores
    for score_key in SCORES:
        wells[score_key] = 0

    if age is not None:
        wells[age] = wells[age].str.lower()
        conditions = [
            (wells[age] == 'unknown'),
            (wells[age] == 'spud pre-1975'),
            (wells[age] == 'spud between 1975 and 1985'),
            (wells[age] == 'spud between 1986 and 1997'),
            (wells[age] == 'spud after 1997')]
        choices = [5, 4, 3, 2, 1]
        wells['age_score'] = np.select(conditions, choices, default=5)

    if deviation is not None:
        wells[deviation] = wells[deviation].str.lower()
        conditions = [
            (wells[deviation] == 'slant'),
            (wells[deviation] == 'horizontal'),
            (wells[deviation] == 'deviated'),
            (wells[deviation] == 'vertical'),
            (wells[deviation] == 'pressure observation (vertical)')]
        choices = [5, 3, 3, 1, 1]
        wells['deviation_score'] = np.select(conditions, choices, default=5)

    if depth is not None:
        wells[depth] = wells[depth].str.lower()
        conditions = [
            (wells[depth] == 'unknown'),
            (wells[depth] == '>=2000 mkb'),
            (wells[depth] == '>1000 mkb <2000 mkb'),
            (wells[depth] == '<=1000 mkb')]
        choices = [5, 5, 3, 1]
        wells['depth_score'] = np.select(conditions, choices, default=5)

    if well_type is not None:
        wells[well_type] = wells[well_type].str.lower()
        conditions = [
            (wells[well_type] == 'h2o injection well'),
            (wells[well_type] == 'h2o disposal well'),
            (wells[well_type] == 'co2 injection well'),
            (wells[well_type] == 'flowing oil well'),
            (wells[well_type] == 'medium risk gas well'),
            (wells[well_type] == 'non-flowing oil well'),
            (wells[well_type] == 'other low risk wells'),
            (wells[well_type] == 'low risk oil h2s well'),
            (wells[well_type] == 'low risk gas well'),
            (wells[well_type] == 'pressure observation well')]
        choices = [5, 5, 5, 4, 3, 3, 2, 1, 1, 1]
        wells['type_score'] = np.select(conditions, choices, default=5)

    if scvf_gm is not None:
        wells[scvf_gm] = wells[scvf_gm].str.lower()
        conditions = [
            (wells[scvf_gm] == 'serious'),
            (wells[scvf_gm] == 'non-serious only'),
            (wells[scvf_gm] == 'not tested'),
            (wells[scvf_gm] == 'none')]
        choices = [5, 4, 3, 1]
        wells['scvf_gm_score'] = np.select(conditions, choices, default=3)

    if casing_failure is not None:
        wells[casing_failure] = wells[casing_failure].str.lower()
        conditions = [
            (wells[casing_failure] == 'unrepaired'),
            (wells[casing_failure] == 'yes'),
            (wells[casing_failure] == 'unknown'),
            (wells[casing_failure] == 'probable'),
            (wells[casing_failure] == 'repaired'),
            (wells[casing_failure] == 'not likely'),
            (wells[casing_failure] == 'none')]
        choices = [5, 5, 5, 4, 2, 2, 1]
        wells['casing_failure_score'] = np.select(conditions, choices, default=5)

    if cement_top is not None:
        wells[cement_top] = wells[cement_top].str.lower()
        conditions = [
            (wells[cement_top] == 'very low'),
            (wells[cement_top] == 'unknown'),
            (wells[cement_top] == 'missing data'),
            (wells[cement_top] == 'low'),
            (wells[cement_top] == 'moderate'),
            (wells[cement_top] == 'above surface casing shoe')]
        choices = [5, 5, 5, 4, 3, 1]
        wells['cement_top_score'] = np.select(conditions, choices, default=5)

    if suspension_status is not None:
        wells[suspension_status] = wells[suspension_status].str.lower()
        conditions = [
            (wells[suspension_status] == 'inactive not suspended'),
            (wells[suspension_status] == 'suspension not downhole compliant'),
            (wells[suspension_status] == 'suspension > 10 year non-compliant suspension'),
            (wells[suspension_status] == 'suspension > 10 yrs'),
            (wells[suspension_status] == 'suspension < 10 yrs')]
        choices = [5, 4, 3, 2, 1]
        wells['suspension_status_score'] = np.select(conditions, choices, default=5)


    wells['leak_score'] = (wells['deviation_score'] + wells['scvf_gm_score']
                           + wells['suspension_status_score'] + wells['cement_top_score']
                           + wells['casing_failure_score'] + wells['age_score']
                           + wells['depth_score'] + wells['type_score'])

    num_scores = sum(x is not None for x in inputs)
    max_score = num_scores*5
    wells['normalized_leak_score'] = wells['leak_score']/max_score

    return np.array(wells['leak_score']), np.array(wells['normalized_leak_score'])

if __name__ == "__main__":
    pass
    # # Example of methods application
    # # Read file with wells data
    # wells = pd.read_csv('well_details.csv')
    # classified_wells = well_data_classifier(wells)
    # leak_score, normalized_leak_score =  leak_scorer(
    #     classified_wells, age='age_class', depth='depth_class',
    #     deviation='dev_class', cement_top='cem_class', well_type='type_class',
    #     scvf_gm='scvf_gm_class')

    # classified_wells['leak_score'] = leak_score
    # classified_wells['normalized_leak_score'] = normalized_leak_score

    # cleaned_output = classified_wells[[
    #     'api_num', 'latitude', 'longitude', 'lateral_status',
    #     'lateral_type', 'x','y', 'deviation_score', 'scvf_gm_score',
    #     'suspension_status_score', 'cement_top_score', 'casing_failure_score',
    #     'age_score', 'depth_score', 'type_score', 'leak_score',
    #     'normalized_leak_score']]

    # cleaned_output.to_csv('well_data_and_scores.csv', index=False)
