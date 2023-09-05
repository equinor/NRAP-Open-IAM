# -*- coding: utf-8 -*-
"""
NRAP OpenIAM Control File reader code
Reads YAML formatted control file and runs IAM analysis.
Use: python openiam_cf.py --file path_to_control_file/ControlFileName.yaml

Created: October 11, 2017
Last modified: August 5th, 2022

Authors: Seth King, Veronika Vasylkivska, Nate Mitchell
"""
import argparse
from datetime import datetime
import logging
import os
import pickle
import random
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np
import yaml

SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(SOURCE_DIR)

# import matk.ordereddict
import openiam as iam
from openiam import IAM_DIR

from openiam.cfi.commons import process_parameters, process_dynamic_inputs
from openiam.cfi.locations import (process_cell_centers,
                                          process_fault_segment_centers,
                                          process_reservoir_locations,
                                          process_wellbore_locations)
from openiam.cfi.analysis import process_analysis
from openiam.cfi.plots import process_plots
from openiam.cfi.output import process_output
from openiam.cfi.text import system_model_to_text, component_models_to_text
from openiam.cfi.strata import (initialize_strata,
                                       process_spatially_variable_strata)

from openiam.cfi.examples_data import GUI_EXAMPLES, CFI_EXAMPLES

import openiam.cfi.workflow as workflow

# The following line creates a parser to parse arguments from the command line to the code.
parser = argparse.ArgumentParser(description='YAML control file IAM reader')

# Create a group of mutually exclusive command line arguments.
# parser will make sure that only one of the arguments in the mutually
# exclusive group is present on the command line.
group = parser.add_mutually_exclusive_group(required=False)

# lhs_sample_size and parstudy_divisions are mutually exclusive.
# --lhs and --parstudy are the keys used to indicate which of the parameters
# are given in the command line.
group.add_argument('--lhs', type=int,
                   dest='lhs_sample_size', default=20,
                   help='Latin Hypercube Sampling mode: enter number of samples')
group.add_argument('--parstudy', type=int,
                   dest='parstudy_divisions', default=3,
                   help='Parameter study mode: enter number of parameter partitions')

# Parameter ncpus has default value 1, i.e., by default, all simulations are run
# sequentially. If a user wants to run simulations in parallel, then
# a different number of cpus to use must be specified
parser.add_argument('--ncpus', type=int, default=1,
                    help='Number of processors to use to run concurrent simulations')

# Parameter file sets the control file to read in.
parser.add_argument('--file', type=str, dest='yaml_cf_name',
                    # default='test_CFI',  # to test all control files examples
                    # default='test_GUI',  # to test all GUI examples
                    # default='../../test/test_control_file.yaml',
                    default='../../examples/Control_Files/ControlFile_ex9a.yaml',
                    # default='../../examples/GUI_Files/14_Forward_HCL.OpenIAM',
                    help='NRAP-Open-IAM Control File Name')
parser.add_argument('--binary', type=bool, dest='binary_file',
                    default=False, help='Set to true for binary control file')
args = parser.parse_args()

output_header = "".join(["NRAP-Open-IAM version: {iam_version}",
                         "\nRuntime: {now} \n"])

pathway_components = ['LookupTableReservoir',
                      'SimpleReservoir',
                      'AnalyticalReservoir',
                      'GenericReservoir',
                      'TheisReservoir',
                      'MultisegmentedWellbore',
                      'CementedWellbore',
                      'CementedWellboreWR',
                      'OpenWellbore',
                      'GeneralizedFlowRate',
                      'AlluviumAquifer',
                      'AlluviumAquiferLF',
                      'DeepAlluviumAquifer',
                      'DeepAlluviumAquiferML',
                      'FutureGen2Aquifer',
                      'FutureGen2AZMI',
                      'GenericAquifer']

reservoir_components = ['LookupTableReservoir',
                        'SimpleReservoir',
                        'AnalyticalReservoir',
                        'GenericReservoir',
                        'TheisReservoir']

wellbore_components = ['MultisegmentedWellbore',
                       'CementedWellbore',
                       'CementedWellboreWR',
                       'OpenWellbore',
                       'GeneralizedFlowRate']

# Components accepting dynamic inputs as scalar
single_input_aquifer_components = ['AlluviumAquifer',
                                   'AlluviumAquiferLF',
                                   'DeepAlluviumAquifer',
                                   'DeepAlluviumAquiferML',
                                   'FutureGen2Aquifer',
                                   'FutureGen2AZMI',
                                   'GenericAquifer']

# Components accepting dynamic inputs as arrays
multi_input_aquifer_components = ['CarbonateAquifer']

alluvium_aquifer_components = ['AlluviumAquifer', 'DeepAlluviumAquifer',
                               'AlluviumAquiferLF', 'DeepAlluviumAquiferML']


def main(yaml_filename, binary_file=False):
    """
    Reads in yaml data control file to create OpenIAM model and run it.

    :param yaml_filename: yaml formatted OpenIAM control file name
    :type filename: str

    :returns: None
    """
    start_time = datetime.now()
    now = start_time.strftime('%Y-%m-%d_%H.%M.%S')
    # Load yaml file data
    if binary_file:
        with open(yaml_filename, 'rb') as cf:
            yaml_data = pickle.load(cf)
    else:
        with open(yaml_filename, 'r') as yaml_cf:
            yaml_data = yaml.load(yaml_cf, Loader=yaml.SafeLoader)

    model_data = yaml_data['ModelParams']

    if 'Analysis' not in model_data:
        model_data['Analysis'] = 'forward'
    if isinstance(model_data['Analysis'], str):
        analysis = model_data['Analysis'].lower()
        analysis_dict = {}
    elif isinstance(model_data['Analysis'], dict):
        analysis_dict = model_data['Analysis']
        analysis = analysis_dict.pop('type').lower()

    if 'OutputDirectory' in model_data:
        out_dir = os.path.join(IAM_DIR, model_data['OutputDirectory'])
        out_dir = out_dir.format(datetime=now)
        if not os.path.exists(os.path.dirname(out_dir)):
            os.mkdir(os.path.dirname(out_dir))
    else:
        out_dir = os.path.join(IAM_DIR, 'output')

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    csv_files_dir = os.path.join(out_dir, 'csv_files')

    model_data['OutputDirectory'] = out_dir
    logging_file_name = os.path.join(out_dir, 'IAM_log.log')
    analysis_log = os.path.join(out_dir, 'Analysis.log')

    # Copy input YAML file
    shutil.copy2(yaml_filename, out_dir)

    # If logging level not specified, set default
    if 'Logging' not in model_data:
        model_data['Logging'] = 'Info'
    log_dict = {'All': logging.NOTSET,
                'Debug': logging.DEBUG,
                'Info': logging.INFO,
                'Warning': logging.WARNING,
                'Error': logging.ERROR,
                'Critical': logging.CRITICAL}

    log_level = log_dict[model_data['Logging']]

    logger = logging.getLogger('')
    # Remove default existing handlers
    logger.handlers.clear()
    logger.setLevel(log_level)
    # logging formatter for log files with more details
    log_formatter1 = logging.Formatter(
        fmt='%(levelname)s %(module)s.%(funcName)s: %(message)s',
        datefmt='%m-%d %H:%M')
    # logging formatter for console output
    log_formatter2 = logging.Formatter(
        fmt='%(levelname)-8s %(message)s',
        datefmt='%m-%d %H:%M')

    # Setup logging to log file
    file_handler = logging.FileHandler(filename=logging_file_name, mode='w')
    file_handler.setLevel(log_level)
    file_handler.setFormatter(log_formatter1)
    logger.addHandler(file_handler)

    # Setup logging to console
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(log_formatter2)
    logger.addHandler(console)

    iam_version_msg = output_header.format(iam_version=iam.__version__, now=now)
    logging.info('\n'+iam_version_msg)

    info_msg = '\nRunning file {}\n'.format(yaml_filename)
    logging.info(info_msg)

    # Check whether user input time points through input file or list
    if 'TimePoints' in model_data:
        time_data = model_data['TimePoints']
        if isinstance(time_data, str):  # if filename is provided
            time_data_file_path = os.path.join(IAM_DIR, time_data)
            if os.path.isfile(time_data_file_path):
                time_array = np.genfromtxt(
                    time_data_file_path, delimiter=",", dtype='f8')
            else:
                logging.debug('Wrong path is provided for time points')
        else:  # list is provided
            time_array = np.array(time_data)

        # Convert time points data to days
        time_array = 365.25 * time_array
    else:
        if 'EndTime' in model_data:
            num_years = model_data['EndTime']
        else:
            num_years = 1
        if 'TimeStep' in model_data:
            time_step = model_data['TimeStep']
        else:
            time_step = 1
        time_array = 365.25 * np.arange(0.0, num_years + time_step, time_step)

    sm_model_kwargs = {'time_array': time_array}  # time is given in days

    # Create system model
    sm = iam.SystemModel(model_kwargs=sm_model_kwargs)

    strata, sm, spatially_variable_strata = initialize_strata(yaml_data, sm)

    # Perform extra needed actions to connect it to the system model
    # For spatially variable stratigraphy, strata is a list of components
    if spatially_variable_strata:
        if not 'LookupTableStratigraphy' in yaml_data['Stratigraphy'][
                'spatiallyVariable']:
            strata[-1].connect_with_system()
    else:
        strata.connect_with_system()

    if 'Workflow' in yaml_data:
        if isinstance(yaml_data['Workflow'], dict):
            if spatially_variable_strata:
                strata_comp = strata[-1]
            else:
                strata_comp = strata

            yaml_data = workflow.iam_workflow_setup(yaml_data, strata_comp)

            # Separate path to the file from the file name itself
            filename_head, filename_tail = os.path.split(yaml_filename)

            yaml_filename_workflow = filename_tail[
                0:filename_tail.rindex('.')] + '_WorkflowSetup' + \
                    filename_tail[filename_tail.rindex('.'):]

            # Write updated input YAML file, with the Workflow entry deleted so
            # it won't be redone if this .yaml is given to openiam_cf.py.
            yaml_data_copy = yaml_data.copy()

            del yaml_data_copy['Workflow']

            # add back the full analysis dictionary
            if analysis in ['lhs', 'parstudy']:
                yaml_data_copy['ModelParams']['Analysis'] = analysis_dict.copy()
                yaml_data_copy['ModelParams']['Analysis']['type'] = analysis

            with open(os.path.join(out_dir, yaml_filename_workflow), 'w') as f:
                yaml.dump(yaml_data_copy, f)
            # f.close()

            del yaml_data_copy

    # Initialize component data that would keep the information
    # for which components the given component provides an output
    for comp_model in model_data['Components']:
        yaml_data[comp_model]['DistributedTo'] = []

    output_prov_cmpnts = []
    for comp_model in model_data['Components']:
        comp_data = yaml_data[comp_model]
        if 'connection' in comp_data:
            comp_data['Connection'] = comp_data['connection']
            comp_data.pop('connection', None)

        if 'Connection' in comp_data:
            if (comp_data['Connection'] == 'none') or (
                    comp_data['Connection'] == 'Dynamic Parameters'):
                comp_data.pop('Connection', None)
                continue

            # Works for lists and single components
            connections = np.array([comp_data['Connection']]).flatten()
            for connect in connections:
                if connect not in output_prov_cmpnts:
                    output_prov_cmpnts.append(connect)
                yaml_data[connect]['DistributedTo'].append(comp_model)

    debug_msg = ''.join([
        'List of components providing output to other components {} and ',
        'to which they provide an output {}']).format(
            output_prov_cmpnts,
            {comp_model: yaml_data[comp_model][
                'DistributedTo'] for comp_model in model_data['Components']})
    logging.debug(debug_msg)

    debug_msg = 'Initial list of components from control file {}'.format(
        model_data['Components'])
    logging.debug(debug_msg)

    # Sort components according to the order of connections
    comp_list = []
    for comp_model in model_data['Components']:
        comp_data = yaml_data[comp_model]
        if comp_model in comp_list:
            continue
        if 'Connection' in comp_data:
            # Works for lists and single components
            connections = np.array([comp_data['Connection']]).flatten()
            for connect in connections:
                if connect not in comp_list:
                    comp_list.append(connect)
            # Add component after adding its connections
            comp_list.append(comp_model)
        elif comp_model in output_prov_cmpnts:
            # If component does not have connections but provides outputs
            # for other components it has to be moved to the
            # beginning of the list to avoid being end up at the end of the list
            comp_list.insert(0, comp_model)
        else:
            # If the component does not have any connections and does not provide
            # outputs for other components
            comp_list.append(comp_model)

    debug_msg = 'List of components sorted according to connections {}'.format(comp_list)
    logging.debug(debug_msg)

    # List of components in the system from pathway_components group
    path_list = []

    # List of important locations for each component
    locations = {}
    # List of injection well locations
    inj_well_locations = {}

    # ad_connect is a list of pairs
    # where 1st element is names of wellbores (to be linked to adapters),
    # and 2nd element is name of aquifer to which the leakage rates should be output
    # adapters are components needed to link aquifer components to wellbores
    ad_connect = []
    # List of all the final components of the system model
    comp_list2 = []

    # Dictionary of system model collectors
    # collectors is a dictionary with keys being names of inputs coming from
    # observations of 'connection' components (mainly wellbore) and having
    # information about to what aquifer leakage rates are input
    sm.collectors = {}
    debug_msg = 'System component list {}'.format(comp_list)
    logging.debug(debug_msg)

    for comp_model in comp_list:
        comp_data = yaml_data[comp_model]

        if 'type' in comp_data:
            comp_data['Type'] = comp_data['type']

        if comp_data['Type'] in pathway_components:
            path_list.append(comp_model)

        if comp_data['Type'] in wellbore_components:
            locations = process_wellbore_locations(
                comp_data, comp_model, locations)

        if comp_data['Type'] in reservoir_components:
            locations, inj_well_locations = process_reservoir_locations(
                comp_data, comp_model,
                locations, inj_well_locations, comp_data['Type'])

        if comp_data['Type'] == 'SealHorizon':
            locations = process_cell_centers(
                comp_data, comp_model, locations)

        if comp_data['Type'] == 'FaultFlow':
            locations = process_fault_segment_centers(
                comp_data, comp_model, locations)

        comp_type = getattr(iam, comp_data['Type'])
        if hasattr(comp_type, 'adapters') and comp_type.adapters:
            path_list.append('adapter')

            if 'Connection' in comp_data:
                # Works for lists and single components
                connections = np.array([comp_data['Connection']]).flatten()
                ad_connect.append([connections, yaml_data[comp_model]['AquiferName']])
                comp_list2.append('adapter')

        if hasattr(comp_type, 'system_collected_inputs') and comp_type.system_collected_inputs:
            # sci dictionary with keys names of model method keyword argument and values
            # names of observations needed to be linked to the arguments
            sci = comp_type.system_collected_inputs
            if 'Connection' in comp_data:
                # For each keyword argument of model method
                for sinput in sci:
                    if not sinput in sm.collectors:
                        sm.collectors[sinput] = {}

                    sm.collectors[sinput][comp_model] = {
                        'Connection': np.array([comp_data['Connection']]).flatten(),
                        'argument': sci[sinput],
                        'data': []}
        comp_list2.append(comp_model)

    debug_msg = ''.join(["Updated component list with adapters {}.\n",
                         "List of adapters {}.\nPath list {}."]).format(
                             comp_list2, ad_connect, path_list)
    logging.debug(debug_msg)

    # Create component list with all pathways
    comps = comp_list2
    comp_list = []
    num_adapters = 0
    adapters = {}
    for comp_model in comps:  # go over all components in the components list
        if comp_model in path_list:  # reservoir, wellbore or aquifer (except carbonate)
            if comp_model == 'adapter':
                num_adapters = process_adapter(
                    yaml_data, comp_list, adapters, ad_connect,
                    locations, num_adapters)
                continue  # proceed to the next component in "comps" list

            # If comp_model is not an adapter get its data
            comp_data = yaml_data[comp_model]

            # For reservoir components providing output for other components
            if comp_data['Type'] in reservoir_components:
                process_reservoir(
                    yaml_data, comp_list, comp_data, comp_model,
                    output_prov_cmpnts, locations, inj_well_locations)
                continue

            if comp_data['Type'] in wellbore_components:
                process_wellbore(yaml_data, comp_list, comp_data,
                                 comp_model, locations)
                continue

            if comp_data['Type'] in single_input_aquifer_components:
                process_aquifer(yaml_data, comp_list, comp_data,
                                comp_model, locations)
                continue


        else:  # if comp_model is not in the path_list
            # For components not in the path_list
            # As of now these are carbonate aquifer, atmospheric ROM,
            # seal horizon, fault flow component, plume stability
            # Carbonate aquifer and atmospheric ROM would not belong to the path_list
            # but they need locations to be collected from their connections
            # Seal horizon and fault don't need locations collected from their connection
            # since they themselves require locations in their setup but they fit this flow
            comp_data = yaml_data[comp_model]

            if comp_data['Type'] in ['SealHorizon', 'FaultFlow']:
                process_connections_to_reservoir(yaml_data, comp_data,
                                                 comp_model, locations, sm)

            if comp_data['Type'] == 'CarbonateAquifer':
                process_carb_aquifer(yaml_data, comp_data, comp_model, locations)

            if comp_data['Type'] == 'AtmosphericROM':
                process_atm_rom(yaml_data, comp_data, comp_model, locations)

            comp_list.append(comp_model)

    debug_msg = "Updated component list {}".format(comp_list)
    logging.debug(debug_msg)
    if spatially_variable_strata:
        name2obj_dict = {'strata': strata[0]}
    else:
        # Create dictionary to map component names to objects
        name2obj_dict = {'strata': strata}
    # Create dictionary of component models with outputs and lists of output names
    output_list = {}
    # Create list of component models
    components = []

    for component_name in comp_list:
        component_data = yaml_data[component_name]
        comp_type = getattr(iam, component_data['Type'])

        components.append(sm.add_component_model_object(
            comp_type(name=component_name, parent=sm)))

        if spatially_variable_strata and not isinstance(
                components[-1], iam.RateToMassAdapter):
            strata, sm = process_spatially_variable_strata(strata, component_name,
                                                           yaml_data, locations, sm)
            strata[-1].connect_with_system()
            name2obj_dict['strata'] = strata[-1]

        name2obj_dict[component_name] = components[-1]

        if hasattr(components[-1], 'connect_with_system'):
            components[-1].connect_with_system(component_data,
                                               name2obj_dict,
                                               locations,
                                               adapters,
                                               output_dir=out_dir)
            # Outputs
            if 'Outputs' in component_data:
                # Due to the specific types of observations for Fault Flow
                # and Seal Horizon component they are added inside the connect_with_system
                # method and not here
                if component_data['Type'] not in [
                        'FaultFlow', 'SealHorizon', 'ChemicalWellSealing', 'GenericAquifer']:
                    comp_outputs = component_data['Outputs']
                    for output in comp_outputs:
                        components[-1].add_obs(output)
                    output_list[components[-1]] = comp_outputs
                else:
                    output_list[components[-1]] = component_data['Outputs']

            # 'continue' below is present to emphasize that connect_with_system
            # method should be able to handle parameters and dynamic inputs
            # of the component
            continue

        # Parameters of component model method
        process_parameters(components[-1], component_data, name2obj_dict)

        # Dynamic keyword arguments of component model method
        process_dynamic_inputs(components[-1], component_data)

        # Handle pressure/saturation locations for wells or reservoir
        if hasattr(components[-1], 'needs_locXY') and components[-1].needs_locXY:
            try:
                if 'locX' in component_data and 'locY' in component_data:
                    components[-1].locX = component_data['locX']
                    components[-1].locY = component_data['locY']
                elif 'coordx' in component_data and 'coordy' in component_data:
                    components[-1].locX = component_data['coordx']
                    components[-1].locY = component_data['coordy']
                else:
                    raise IndexError
            except IndexError:
                info_msg = ''.join(['Coordinates of location of interest ',
                                    'are unknown, placing at (100, 100)'])
                logging.info(info_msg)
                components[-1].locX = 100
                components[-1].locY = 100

        # Handle injection well locations
        if hasattr(components[-1], 'needs_injXY') and components[-1].needs_injXY:
            try:
                components[-1].injX = component_data['injX']
                components[-1].injY = component_data['injY']
            except KeyError:
                pass

        if hasattr(components[-1], 'needsXY'):
            if components[-1].needsXY:
                components[-1].model_kwargs['x'] = component_data['locX']
                components[-1].model_kwargs['y'] = component_data['locY']

        # Outputs
        if 'Outputs' in component_data:
            comp_outputs = component_data['Outputs']
            output_list[components[-1]] = comp_outputs
            for output in comp_outputs:
                components[-1].add_obs(output)

        # Make model connections
        if 'Connection' in component_data:
            connection = None
            try:
                connection = name2obj_dict[component_data['Connection']]
            except KeyError:
                pass

            if hasattr(components[-1], 'system_inputs'):
                for sinput in components[-1].system_inputs:
                    # TODO add adapter linked obs for system_inputs
                    # if adapters different from RateToMassAdapter
                    # will be ever created
                    connection.add_obs_to_be_linked(sinput)
                    components[-1].add_kwarg_linked_to_obs(
                        sinput, connection.linkobs[sinput])

            if (hasattr(components[-1], 'system_collected_inputs') and
                    components[-1].system_collected_inputs):
                collectors = sm.collectors
                for cdict in collectors.values():
                    argument = cdict['argument'].format(
                        aquifer_name=component_data['AquiferName'])
                    connections = cdict['Connection']

                    for connect in connections:
                        for ind in range(locations[connect]['number']):
                            cname = connect + '_{0:03}'.format(ind)
                            for ad_nm, adptr_item in adapters.items():
                                if (adptr_item['Connection'] == cname) and (
                                        adptr_item['AquiferName'] == component_data[
                                            'AquiferName']):
                                    aname = ad_nm

                            adapter = name2obj_dict[aname]
                            connector = name2obj_dict[cname]
                            if argument in adapter.linkobs:
                                cdict['data'].append(adapter.linkobs[argument])
                            else:
                                cdict['data'].append(connector.linkobs[argument])

                sci = components[-1].system_collected_inputs
                for sinput in sci:
                    components[-1].add_kwarg_linked_to_collection(
                        sinput, collectors[sinput][component_name]['data'])

            if hasattr(components[-1], 'composite_inputs'):
                for key in components[-1].composite_inputs:
                    components[-1].add_composite_par(
                        key, components[-1].composite_inputs[key].format(
                            driver=connection.name,
                            selfm=components[-1].name,
                            strata='strata'))
        # End of "if Connection in ..." statement

        if hasattr(components[-1], 'system_params'):
            if yaml_data[component_name]['Type'] in [
                    'SimpleReservoir', 'AnalyticalReservoir', 'MultisegmentedWellbore']:
                if spatially_variable_strata:
                    if 'numberOfShaleLayers' in strata[-1].pars:
                        connect = strata[-1].pars
                    elif 'numberOfShaleLayers' in strata[-1].deterministic_pars:
                        connect = strata[-1].deterministic_pars
                    else:
                        connect = strata[-1].default_pars
                else:
                    if 'numberOfShaleLayers' in strata.pars:
                        connect = strata.pars
                    elif 'numberOfShaleLayers' in strata.deterministic_pars:
                        connect = strata.deterministic_pars
                    else:
                        connect = strata.default_pars

                components[-1].add_par_linked_to_par(
                    'numberOfShaleLayers', connect['numberOfShaleLayers'])

                nSL = connect['numberOfShaleLayers'].value
                components[-1].system_params = [
                    'shale{}Thickness'.format(ind) for ind in range(1, nSL + 1)]\
                    + ['aquifer{}Thickness'.format(ind) for ind in range(1, nSL)]\
                    + ['reservoirThickness', 'datumPressure']

            if yaml_data[component_name]['Type'] == 'GenericReservoir':
                components[-1].system_params = ['reservoirThickness',
                                                'reservoirDepth']

            for sparam in components[-1].system_params:
                connect = None

                if spatially_variable_strata:
                    if sparam in strata[-1].pars:
                        connect = strata[-1].pars
                    elif sparam in strata[-1].deterministic_pars:
                        connect = strata[-1].deterministic_pars
                    elif sparam in strata[-1].default_pars:
                        connect = strata[-1].default_pars
                    else:
                        info_msg = 'Unable to find parameter {}.'.format(sparam)
                        logging.info(info_msg)

                else:
                    if sparam in strata.composite_pars:
                        connect = strata.composite_pars
                    elif sparam in strata.pars:
                        connect = strata.pars
                    elif sparam in strata.deterministic_pars:
                        connect = strata.deterministic_pars
                    elif sparam in strata.default_pars:
                        connect = strata.default_pars
                    else:
                        info_msg = 'Unable to find parameter {}.'.format(sparam)
                        logging.info(info_msg)

                components[-1].add_par_linked_to_par(sparam, connect[sparam])
    # End Component model loop

    # Setup Analysis
    if analysis == 'lhs':
        if 'siz' not in analysis_dict:
            analysis_dict['size'] = args.lhs_sample_size
        if 'seed_size' in analysis_dict:
            seed_size = [int(ss) for ss in analysis_dict.pop('seed_size').strip('()').split(',')]
            analysis_dict['seed'] = random.randint(*seed_size)
        if 'seed' not in analysis_dict:
            analysis_dict['seed'] = random.randint(500, 1100)
    elif analysis == 'parstudy':
        if 'nvals' not in analysis_dict:
            analysis_dict['nvals'] = args.parstudy_divisions
    model_data['Analysis_type'] = analysis
    model_data['Analysis_dict'] = analysis_dict
    # End setup

    calc_start_time = datetime.now()
    setup_time = calc_start_time - start_time
    debug_msg = 'Model setup time: {}'.format(setup_time)
    logging.debug(debug_msg)

    # Run Analysis
    if analysis == 'forward':
        results = sm.forward()
        s = None
    elif analysis == 'lhs':
        s = sm.lhs(**analysis_dict)
        results = s.run(cpus=args.ncpus, verbose='progress', logfile=None)
        # New line in console
        print("")
    elif analysis == 'parstudy':
        s = sm.parstudy(**analysis_dict)
        results = s.run(cpus=args.ncpus, verbose='progress', logfile=None)
        # New line in console
        print("")
    else:
        raise ValueError('Analysis type {atype} not found.\nNo analysis ran'.format(atype=analysis))

    yaml_data['Results'] = results
    yaml_data['s'] = s
    yaml_data['sm'] = sm
    yaml_data['output_list'] = output_list
    yaml_data['components'] = components
    yaml_data['time_array'] = time_array
    # End Analysis

    calc_end_time = datetime.now()
    calc_time = calc_end_time - calc_start_time
    debug_msg = 'Model analysis time: {}'.format(calc_time)
    logging.debug(debug_msg)

    if 'OutputDirectory' in model_data:
        with open (os.path.join(out_dir, 'openiam_version_info.txt'), 'w') as f:
            f.write(iam_version_msg)
        process_print_system_options(model_data, sm, out_dir, analysis)
        clean_components(yaml_data)
        process_output(yaml_data, model_data, output_list, out_dir, sm, s,
                       analysis, time_array, csv_files_dir)

    if 'Analysis' in yaml_data:
        process_analysis(yaml_data, model_data, sm, s, output_list,
                         analysis, time_array)

    if 'Plots' in yaml_data:
        process_plots(yaml_data, model_data, sm, s, output_list, analysis,
                      time_array, components, locations)

    total_time = datetime.now() - start_time
    info_msg = '\nAnalysis completed at {}.\nTotal run time: {} \n'.format(
        datetime.now().strftime('%Y-%m-%d_%H.%M.%S'), total_time)
    logging.info(info_msg)

    if binary_file:
        # For GUI produced binary control files
        info_msg = ''.join([
            '\nSimulation results can be found in the output folder: \n{}.',
            '\nProceed to the Post Processing tab to access plotting options. \n']).format(
                model_data['OutputDirectory'])
    else:
        info_msg = ''.join([
            '\nSimulation results and plots can be found ',
            'in the output folder: \n{}']).format(model_data['OutputDirectory'])
    logging.info(info_msg)

    if 'Workflow' in yaml_data:
        workflow.workflow_analysis(yaml_data, sm, analysis)

    # Remove all handlers from the logger for proper work in the consecutive runs
    logger.handlers.clear()

    # Needed for control file test in a test suite
    return True


def process_print_system_options(model_data, sm, out_dir, analysis):
    """
    Decompose, filter, and print the system model and component model objects
    into a human-readable text file for debugging and information purposes.

    Parameters
    ----------
    model_data : TYPE
        DESCRIPTION.
    sm : TYPE
        DESCRIPTION.
    out_dir : TYPE
        DESCRIPTION.
    analysis : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    to_print_system = model_data.get('GenerateSystemPrintFile', False)
    to_print_components = model_data.get('GenerateComponentsPrintFiles', False)

    # Create file containing info about system model sm
    if to_print_system:
        system_model_to_text(sm, out_dir, analysis)

    if to_print_components:
        component_models_to_text(sm.__dict__['component_models'], out_dir, analysis)


def process_adapter(yaml_data, comp_list, adapters, ad_connect, locations,
                    prev_num_adapters):
    """
    Add adapter components to yaml data based on presence of wellbore
    and aquifer components.
    """
    # ad_connect is a list of pairs
    # where 1st element is names of wellbores (to be linked to adapters),
    # and 2nd element is name of aquifer to which the leakage rates should be output
    # adapters are components needed to link aquifer components to wellbores
    # ad_connect_name_base is possibly a list of connections
    # Remove the first element from the list and return it
    ad_connect_name_base = ad_connect.pop(0)
    # Index to keep track of the total number of adapters
    ad_ind = prev_num_adapters
    # Since adapter through aquifer can be connected to
    # several groups of wells
    # connect_cmpnt represents a particular group of wells
    # ad_connect_name_base[0] is a list of wellbore components
    # which provide input for adapters
    # ad_connect_name_base[1] is a name of aquifer to which
    # leakage rates are of interest
    for connect_cmpnt in ad_connect_name_base[0]:
        # Depending on the number of wells in the given group
        for ind in range(locations[connect_cmpnt]['number']):
            ad_name = 'adapter_{0:03}'.format(ad_ind)
            # ad_connect_name is name of (wellbore) component
            # from which output will be requested
            ad_connect_name = connect_cmpnt + '_{0:03}'.format(ind)
            # Save information about adapter component in the yaml_data dictionary
            # where information about other components is saved as well
            yaml_data[ad_name] = {'Type': 'RateToMassAdapter',
                                  'Connection': ad_connect_name,
                                  'AquiferName': ad_connect_name_base[1]}
            # Save information about adapter component in adapters dictionary
            adapters[ad_name] = yaml_data[ad_name]
            # Add adapters to the component list
            comp_list.append(ad_name)
            # Increase count index for adapters
            ad_ind = ad_ind + 1

    return ad_ind


def process_reservoir(yaml_data, comp_list, comp_data, comp_model,
                      output_prov_cmpnts, locations, inj_well_locations):
    """
    Transform reservoir component data into setup of multiple components
    of the same type corresponding to different locations.
    """
    # Check whether location of injection well is provided
    if comp_model in inj_well_locations:
        comp_data_copy = comp_data.copy()
        comp_data_copy['injX'] = inj_well_locations[comp_model][0]
        comp_data_copy['injY'] = inj_well_locations[comp_model][1]
        yaml_data[comp_model] = comp_data_copy
        comp_data = yaml_data[comp_model]

    # Check if reservoir provides output for other components
    res_ind = 0
    if comp_model in output_prov_cmpnts:
        # comp_data['DistributedToIndices'] is a dictionary with keys
        # being names of components to which reservoir data is provided
        # and values being tuples of the first and the last indices
        # used to create corresponding reservoir components from which
        # the corresponding reservoir data will be used
        comp_data['DistributedToIndices'] = {}

        # Here, comp_data['DistributedTo'] is a list of well groups
        # connected to the given reservoir component
        for cmpnt in comp_data['DistributedTo']:
            # Reservoir component can provide input for
            # several groups of wells. To make sure each group
            # is connected to the right reservoir component
            # we save the reservoir component indices
            # corresponding to the well group
            comp_data['DistributedToIndices'][cmpnt] = (
                res_ind, res_ind + locations[cmpnt]['number'])

            for ind in range(locations[cmpnt]['number']):
                comp_name = comp_model + '_{0:03}'.format(res_ind)
                comp_data_copy = comp_data.copy()
                # Locations for reservoir component stored in component data
                comp_data_copy['coordx'] = locations[cmpnt]['coordx'][ind]
                comp_data_copy['coordy'] = locations[cmpnt]['coordy'][ind]
                try:
                    # Check whether z-coordinates are provided
                    # There can be IndexError and KeyError
                    # KeyError means there are no z coordinates provided
                    # IndexError means there are not enough coordinates provided
                    comp_data_copy['coordz'] = locations[cmpnt]['coordz'][ind]
                except KeyError:  # handling absence of z-coordinates
                    pass
                except IndexError:  # handling accessing index out of provided range
                    # Here are the wellbores who are provided with not enough z-coordinates
                    err_msg = ''.join(['Not enough z-coordinates are provided ',
                                       'for component {}']).format(cmpnt)
                    logging.error(err_msg)
                    raise IndexError(err_msg) from None

                yaml_data[comp_name] = comp_data_copy
                # Add reservoir component corresponding to a specific location
                # to the list of final components
                comp_list.append(comp_name)
                res_ind = res_ind + 1

    # If separate (from wellbore) locations are specified for a reservoir component
    if comp_model in locations:
        for ind in range(locations[comp_model]['number']):
            comp_name = comp_model + '_{0:03}'.format(res_ind)
            comp_data_copy = comp_data.copy()
            comp_data_copy['coordx'] = locations[comp_model]['coordx'][ind]
            comp_data_copy['coordy'] = locations[comp_model]['coordy'][ind]
            try:
                # Check whether z-coordinates are provided
                # There can be IndexError and KeyError
                # KeyError means there are no z coordinates provided
                # IndexError means there are not enough coordinates provided
                comp_data_copy['coordz'] = locations[comp_model]['coordz'][ind]
            except KeyError:  # handling absence of z-coordinates
                pass
            except IndexError:  # handling accessing index out of provided range
                # Here are the reservoirs who are provided with not enough z-coordinates
                err_msg = ''.join(['Not enough z-coordinates are provided ',
                                   'for component {}']).format(comp_model)
                logging.error(err_msg)
                raise IndexError(err_msg) from None
            yaml_data[comp_name] = comp_data_copy
            comp_list.append(comp_name)
            res_ind = res_ind + 1


def process_wellbore(yaml_data, comp_list, comp_data, comp_model, locations):
    """
    Transform wellbore component data into setup of multiple components
    of the same type corresponding to different locations.
    """
    for ind in range(locations[comp_model]['number']):
        comp_name = comp_model + '_{0:03}'.format(ind)
        comp_data_copy = comp_data.copy()

        # If the well locations are provided
        if locations[comp_model]['coordx']:
            comp_data_copy['locX'] = locations[comp_model]['coordx'][ind]
            comp_data_copy['locY'] = locations[comp_model]['coordy'][ind]
            # Since 'Connection' determines whether z-coordinates are provided and
            # the connection component (typically, reservoir) only need them,
            # then the wellbore component do not need to save them at this stage.
            # This might change in the future.

        yaml_data[comp_name] = comp_data_copy
        comp_list.append(comp_name)

    # For wells the possible connection is reservoir component
    if 'Connection' in comp_data:
        connect = comp_data['Connection']
        # This is where we use the indices of the corresponding
        # reservoir components to link well to the right reservoir
        # component
        ind_offset = yaml_data[connect]['DistributedToIndices'][comp_model][0]

        for ind in range(locations[comp_model]['number']):
            comp_name = comp_model + '_{0:03}'.format(ind)
            yaml_data[comp_name]['Connection'] = (
                connect + '_{0:03}'.format(ind + ind_offset))


def process_aquifer(yaml_data, comp_list, comp_data, comp_model, locations):
    """
    Transform aquifer component data into setup of multiple components
    of the same type corresponding to different locations.
    """
    if 'Connection' in comp_data:
        connect = comp_data['Connection']
        if connect in locations:  # connection is a wellbore component
            for ind in range(locations[connect]['number']):
                comp_name = comp_model + '_{0:03}'.format(ind)
                comp_data_copy = yaml_data[comp_model].copy()
                comp_data_copy['Connection'] = connect + '_{0:03}'.format(ind)
                comp_data_copy['locX'] = locations[connect]['coordx'][ind]
                comp_data_copy['locY'] = locations[connect]['coordy'][ind]
                yaml_data[comp_name] = comp_data_copy
                comp_list.append(comp_name)
        else:
            err_msg = ''.join([
                'Connection {} does not exist and is not possible for component {}. ',
                'Please check and update the control file.']).format(
                    connect, comp_model)
            logging.error(err_msg)
    else:
        # No connections. Possible situations: testing of single input
        # aquifer components as single components
        if 'Locations' in comp_data:
            num_pathways = len(comp_data['Locations']['coordx'])
            well_data_x = comp_data['Locations']['coordx']
            well_data_y = comp_data['Locations']['coordy']

        else:
            num_pathways = 1
            info_msg = ''.join(
                ['Locations are not specified for {} component. '.format(comp_model),
                 'Default location (100, 100) will be used.'])
            logging.info(info_msg)
            well_data_x = [100.0]
            well_data_y = [100.0]

        for ind in range(num_pathways):
            comp_name = comp_model + '_{0:03}'.format(ind)
            comp_data_copy = yaml_data[comp_model].copy()
            comp_data_copy['locX'] = well_data_x[ind]
            comp_data_copy['locY'] = well_data_y[ind]
            yaml_data[comp_name] = comp_data_copy
            comp_list.append(comp_name)


def process_carb_aquifer(yaml_data, comp_data, comp_model, locations):
    """
    Analyze carbonate aquifer component data to determine locations
    of leakage sources.
    """
    if 'Connection' in comp_data:
        connections = np.array([comp_data['Connection']]).flatten()

        yaml_data[comp_model]['locX'] = []
        yaml_data[comp_model]['locY'] = []

        for connect in connections:
            if connect in locations:
                yaml_data[comp_model]['locX'] = (
                    yaml_data[comp_model]['locX'] + locations[connect]['coordx'])
                yaml_data[comp_model]['locY'] = (
                    yaml_data[comp_model]['locY'] + locations[connect]['coordy'])
    else:
        if 'Locations' in comp_data:
            yaml_data[comp_model]['locX'] = comp_data['Locations']['coordx']
            yaml_data[comp_model]['locY'] = comp_data['Locations']['coordy']

        else:
            info_msg = ''.join([
                'Locations are not specified for {} component.',
                'Default location (100, 100) will be used.']).format(comp_model)
            logging.info(info_msg)
            yaml_data[comp_model]['locX'] = [100.0]
            yaml_data[comp_model]['locY'] = [100.0]


def process_atm_rom(yaml_data, comp_data, comp_model, locations):
    """
    Analyze atmospheric impact component data to determine locations
    of leakage sources.
    """
    if 'Connection' in comp_data:
        connections = np.array([comp_data['Connection']]).flatten()

        yaml_data[comp_model]['locX'] = []
        yaml_data[comp_model]['locY'] = []

        for connect in connections:
            if connect in locations:
                yaml_data[comp_model]['locX'] = (
                    yaml_data[comp_model]['locX'] + locations[connect]['coordx'])
                yaml_data[comp_model]['locY'] = (
                    yaml_data[comp_model]['locY'] + locations[connect]['coordy'])
    else:
        if 'Locations' in comp_data:
            yaml_data[comp_model]['locX'] = comp_data['Locations']['coordx']
            yaml_data[comp_model]['locY'] = comp_data['Locations']['coordy']

        else:
            info_msg = ''.join([
                'Locations are not specified for {} component.',
                'Default location (100, 100) will be used.']).format(comp_model)
            logging.info(info_msg)
            yaml_data[comp_model]['locX'] = [100.0]
            yaml_data[comp_model]['locY'] = [100.0]


def process_connections_to_reservoir(yaml_data, comp_data, comp_model, locations, sm):
    """
    Setup collectors required for setup of Seal Horizon and Fault Flow components.
    """
    if 'Connection' in comp_data:
        # Get name of reservoir component whose output will be used
        # by Seal or Fault component
        connect = comp_data['Connection']

        # We use the indices of the corresponding reservoir components
        # to link Seal/Fault to the right reservoir components
        ind_offset = yaml_data[connect]['DistributedToIndices'][comp_model][0]

        # Get number of reservoir components that will be linked
        num = locations[comp_model]['number']

        # We need to update comp_data['Connection'] with a list
        # of reservoir components that will be linked to the Seal or Fault
        # component
        comp_data['Connection'] = [connect + '_{0:03}'.format(
            ind + ind_offset) for ind in range(num)]

        for obs_nm in ['pressure', 'CO2saturation']:
            if obs_nm not in sm.collectors:
                sm.collectors[obs_nm] = {}
            sm.collectors[obs_nm][comp_model] = {
                'Connection': comp_data['Connection'],
                'data': []}

def clean_components(yaml_data):
    """
    Remove or reset problematic attributes of components (if present in the setup)
    before pickling the data.
    """
    for cmpnt_obj in yaml_data['components']:
        cmpnt_nm = cmpnt_obj.name
        # Reset attributes containing reference to the ML models and unnecessary data
        if yaml_data[cmpnt_nm]['Type'] in ['DeepAlluviumAquiferML', 'GenericAquifer',
                                           'AlluviumAquiferLF', 'KimberlinaWellbore',
                                           'CementedWellbore', 'CementedWellboreWR',
                                           'HydrocarbonLeakage']:
            cmpnt_obj.sol = None
            continue

        if yaml_data[cmpnt_nm]['Type'] == 'GenericReservoir':
            cmpnt_obj.sol = None
            cmpnt_obj.x_scaler_p = None
            cmpnt_obj.y_scaler_p = None
            cmpnt_obj.model_p = None
            cmpnt_obj.x_scaler_s = None
            cmpnt_obj.y_scaler_s = None
            cmpnt_obj.model_s = None
            continue


if __name__ == "__main__":

    __spec__ = None

    if args.yaml_cf_name == 'test_CFI':
        examples_description = CFI_EXAMPLES

        for key, example_feature in examples_description.items():
            print('Testing file ControlFile_ex{}.yaml:'.format(key))
            print('Example features: {}'.format(example_feature))
            try:
                main('../../examples/Control_Files/ControlFile_ex{}.yaml'.format(key), False)
                plt.close('all')
            except:
                print('An exception was raised during run.')

        print('Test of all control file interface examples is done.')

    elif args.yaml_cf_name == 'test_GUI':

        examples_description = GUI_EXAMPLES

        for key in examples_description:
            print('Testing file {}:'.format(key))
            main('../../examples/GUI_Files/{}'.format(key), True)
            plt.close('all')

        print('Test of all GUI examples is done.')

    else:
        if args.yaml_cf_name[-5:].lower() == '.yaml':
            main(args.yaml_cf_name, False)
        elif args.yaml_cf_name[-8:].lower() == '.openiam':
            main(args.yaml_cf_name, True)
