'''
This example demonstrates that the results from Control Files and scripts are the same.

Created: 6-17-2022
Last modified: 11-10-2022

Author: Nate Mitchell
Leidos Research Support Team
Nathaniel.Mitchell@netl.doe.gov

Example of run:
$ python test_file1_script_vs_controlfile.py
'''
import sys
import os
import datetime
import shutil
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.sep.join(['..', '..', '..', 'source']))
from openiam import (MultisegmentedWellbore, SystemModel, Stratigraphy,
                     SimpleReservoir, openiam_cf)

if __name__ == "__main__":

    # Number of years assessed in each model run
    num_years = 100

    # time is given in days
    time_array = 365.25 * np.arange(0.0, num_years + 1)
    sm_model_kwargs = {'time_array': time_array}

    # Set up directories for saving output
    original_location = os.getcwd()
    cf_yaml_location = os.path.join(
        original_location, '..', '..', 'Control_Files', 'test_files')

    Open_IAM_Output_folder = os.path.join(original_location, '..', '..', '..', 'output')
    os.chdir(Open_IAM_Output_folder)

    master_output_dir = 'script_vs_cf_' + str(datetime.date.today())
    master_output_location = os.path.join(Open_IAM_Output_folder, master_output_dir)

    if os.path.isdir(master_output_location) is False:
        os.mkdir(master_output_dir)

    # Run the .yaml control file
    cf_output_dir_name = 'cf_output'

    # This is the initial location, BEFORE the folder is moved
    cf_output_location = os.path.join(Open_IAM_Output_folder, cf_output_dir_name)

    # Get the name of the .yaml file for this model version
    yaml_filename = 'test_control_file1.yaml'

    # Go to the Control_Files folder
    os.chdir(cf_yaml_location)

    # Run the Control File
    openiam_cf.main(yaml_filename)

    # Close the figures made by openiam_cf.py
    plt.close('Pressure')
    plt.close('CO2saturation')
    plt.close('CO2_Leakage')
    plt.close('CO2_Mass')
    plt.close('Brine_Leakage')
    plt.close('Brine_Mass')

    # If the folder for the control file output is already there (e.g., from a previous run),
    # shutil will freak out when attempting to move the folder. In that case, delete the folder.
    # If you have a figure from that folder open, this part will cause an error.
    if os.path.isdir(os.path.join(master_output_location, cf_output_dir_name)):
        shutil.rmtree(os.path.join(master_output_location, cf_output_dir_name))

    # Move the control file output folder into the master_output_location
    shutil.move(cf_output_location, master_output_location)
    # This is the final location for the control file output for this model version
    cf_output_location = os.path.join(master_output_location, cf_output_dir_name)

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add stratigraphy component
    strata = sm.add_component_model_object(Stratigraphy(name='strata', parent=sm))

    strata.add_par('numberOfShaleLayers', value=3, vary=False)
    strata.add_par('shale1Thickness', value=200.0, vary=False)
    strata.add_par('shale2Thickness', value=200.0, vary=False)
    strata.add_par('shale3Thickness', value=200.0, vary=False)
    strata.add_par('aquifer1Thickness', value=200.0, vary=False)
    strata.add_par('aquifer2Thickness', value=200.0, vary=False)
    strata.add_par('reservoirThickness', value=200.0, vary=False)

    # Add composite parameters of Stratigraphy component
    strata.add_composite_par(
        'depth', expr='+'.join(['strata.shale1Thickness',
                                'strata.shale2Thickness',
                                'strata.shale3Thickness',
                                'strata.aquifer1Thickness',
                                'strata.aquifer2Thickness']))

    # Add reservoir component
    sres = sm.add_component_model_object(
        SimpleReservoir(name='sres', parent=sm, locX=100, locY=100))

    sres.add_par('injRate', value=1.0, vary=False)
    sres.add_par('logResPerm', value=-12.0, vary=False)

    sres.add_par_linked_to_par('numberOfShaleLayers',
                               strata.deterministic_pars['numberOfShaleLayers'])
    sres.add_par_linked_to_par('shale1Thickness',
                               strata.deterministic_pars['shale1Thickness'])
    sres.add_par_linked_to_par('shale2Thickness',
                               strata.deterministic_pars['shale2Thickness'])
    sres.add_par_linked_to_par('shale3Thickness',
                               strata.deterministic_pars['shale3Thickness'])
    sres.add_par_linked_to_par('aquifer1Thickness',
                               strata.deterministic_pars['aquifer1Thickness'])
    sres.add_par_linked_to_par('aquifer2Thickness',
                               strata.deterministic_pars['aquifer2Thickness'])
    sres.add_par_linked_to_par('reservoirThickness',
                               strata.deterministic_pars['reservoirThickness'])
    sres.add_par_linked_to_par('datumPressure', strata.default_pars['datumPressure'])

    # Add observations of reservoir component model
    sres.add_obs_to_be_linked('pressure')
    sres.add_obs_to_be_linked('CO2saturation')
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')
    sres.add_obs('mass_CO2_reservoir')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))

    # The other parameters will resort to default settings
    ms.add_par('wellRadius', value=0.05, vary=False)
    ms.add_par('logWellPerm', value=-13.0, vary=False)

    # Add linked parameters: common to reservoir and wellbore components
    ms.add_par_linked_to_par('numberOfShaleLayers',
                             strata.deterministic_pars['numberOfShaleLayers'])
    ms.add_par_linked_to_par('shale1Thickness',
                             strata.deterministic_pars['shale1Thickness'])
    ms.add_par_linked_to_par('shale2Thickness',
                             strata.deterministic_pars['shale2Thickness'])
    ms.add_par_linked_to_par('shale3Thickness',
                             strata.deterministic_pars['shale3Thickness'])
    ms.add_par_linked_to_par('aquifer1Thickness',
                             strata.deterministic_pars['aquifer1Thickness'])
    ms.add_par_linked_to_par('aquifer2Thickness',
                             strata.deterministic_pars['aquifer2Thickness'])
    ms.add_par_linked_to_par('reservoirThickness',
                             strata.deterministic_pars['reservoirThickness'])
    ms.add_par_linked_to_par('datumPressure', strata.default_pars['datumPressure'])

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

    # Add observations of multisegmented wellbore component model
    ms.add_obs_to_be_linked('CO2_aquifer1')
    ms.add_obs_to_be_linked('CO2_aquifer2')
    ms.add_obs_to_be_linked('brine_aquifer1')
    ms.add_obs_to_be_linked('brine_aquifer2')
    ms.add_obs_to_be_linked('mass_CO2_aquifer1')
    ms.add_obs_to_be_linked('mass_CO2_aquifer2')
    ms.add_obs_to_be_linked('brine_atm')
    ms.add_obs_to_be_linked('CO2_atm')
    ms.add_obs('brine_aquifer1')
    ms.add_obs('CO2_aquifer1')
    ms.add_obs('brine_aquifer2')
    ms.add_obs('CO2_aquifer2')

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically

    # Get results
    os.chdir(cf_output_location)

    # CO2 leakage rates to Aquifer 1
    CO2_aquifer1 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer1')
    CO2_aquifer1_norm = np.empty(len(CO2_aquifer1))
    load = np.loadtxt(os.path.join('csv_files', 'msw1_000.CO2_aquifer1.csv'),
                      delimiter=',', skiprows=1)[:, 1]
    for norm_ref, val in enumerate(load):
        if float(val) != 0:
            CO2_aquifer1_norm[norm_ref] = CO2_aquifer1[norm_ref] / float(val)
        else:
            CO2_aquifer1_norm[norm_ref] = np.NaN

    # CO2 leakage rates to Aquifer 2
    CO2_aquifer2 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer2')
    CO2_aquifer2_norm = np.empty(len(CO2_aquifer2))

    load = np.loadtxt(os.path.join('csv_files', 'msw1_000.CO2_aquifer2.csv'),
                      delimiter=',', skiprows=1)[:, 1]
    for norm_ref, val in enumerate(load):
        if float(val) != 0:
            CO2_aquifer2_norm[norm_ref] = CO2_aquifer2[norm_ref] / float(val)
        else:
            CO2_aquifer2_norm[norm_ref] = np.NaN

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    print('CO2_aquifer1 from Script:\n', CO2_aquifer1, sep='\n')
    print('CO2_aquifer2 from Script:\n', CO2_aquifer2, sep='\n')
    print('CO2_aquifer1, Script results / CF Results: ', CO2_aquifer1_norm, sep='\n')
    print('CO2_aquifer2, Script results / CF Results:', CO2_aquifer2_norm, sep='\n')
