'''
Example illustrates linking of simple reservoir and multisegmented
wellbore component. The system model is run for a single time point -
initial time 0.

The code takes in several optional arguments which can be changed if needed;
the code can also be run with default values:

lhs      - number of samples for Latin Hypercube Sampling mode
parstudy - number of parameter partitions
ncpus    - number of processors to run concurrent simulations

This example also illustrates how to save all the outputs produced by the simulation.
Changing variable 'save_output' at the beginning of simulation (line 64)
from True to False allows to cancel saving of the outputs.

Examples of run:
$ python iam_sys_reservoir_mswell_lhs_parstudy.py --lhs 10 --ncpus 2
$ python iam_sys_reservoir_mswell_lhs_parstudy.py --parstudy 10 --ncpus 2
'''

import sys
import os
import argparse     # argparse module makes it easy to write command-line interfaces.

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import SystemModel, SimpleReservoir, MultisegmentedWellbore

# The following line creates a parser to parse arguments from the command line to the code.
parser = argparse.ArgumentParser(description='Run examples.')

# Create a group of mutually exclusive command line arguments.
# parser will make sure that only one of the arguments in the mutually
# exclusive group is present on the command line.
group = parser.add_mutually_exclusive_group(required=False)

# lhs_sample_size and parstudy_divisions are mutually exclusive.
# --lhs and --parstudy are the keys used to indicate which of the parameters
# are given in the command line.
group.add_argument('--lhs', type=int,
    dest='lhs_sample_size', default=None,
    help='Latin Hypercube Sampling mode: enter number of samples')
group.add_argument('--parstudy', type=int,
    dest='parstudy_divisions', default=30,
    help='Parameter study mode: enter number of parameter partitions')

# Parameter ncpus has default value 1, i.e., by default, all simulations are run
# sequentially. If a user wants to run simulations in parallel, he/she has to
# specify the number of cpus to use.
parser.add_argument('--ncpus', type=int, default=4,
    help='Number of processors to use to run concurrent simulations')

args = parser.parse_args()


if __name__=='__main__':
    # For multiprocessing in Spyder
    __spec__ = None

    # Change the variable value to False if saving the outputs is not needed.
    # By default, the results will be saved in the folder 'output/csv_files' within root
    # folder of NRAP-Open-IAM
    save_output = True

    # Simple example for single time point which is 0.0 by default
    # System model created without arguments assumes single time point 0
    sm = SystemModel()

    # Add reservoir component
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))

    # Add parameters of reservoir component model
    sres.add_par('numberOfShaleLayers', value=3, vary=False)
    sres.add_par('logResPerm', min=-13., max=-11., value=-12.)

    # Add observations of reservoir component model to be used by the next component
    sres.add_obs_to_be_linked('pressure')
    sres.add_obs_to_be_linked('CO2saturation')
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))
    # A lot of parameters of multisegmented wellbore component
    # are the same as for the reservoir component
    # Add parameters linked to the same parameters from reservoir model
    ms.add_par_linked_to_par(
        'numberOfShaleLayers', sres.deterministic_pars['numberOfShaleLayers'])
    ms.add_par_linked_to_par('shale1Thickness', sres.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('shale2Thickness', sres.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('shale3Thickness', sres.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('aquifer1Thickness', sres.default_pars['aquiferThickness'])
    ms.add_par_linked_to_par('aquifer2Thickness', sres.default_pars['aquiferThickness'])
    ms.add_par_linked_to_par('reservoirThickness', sres.default_pars['reservoirThickness'])
    ms.add_par_linked_to_par('datumPressure', sres.default_pars['datumPressure'])

    # Add parameters specific to the wellbore component
    ms.add_par('logWellPerm', min=-14., max=-11., value=-13.5)

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])
    ms.add_obs('brine_aquifer1')
    ms.add_obs('CO2_aquifer1')
    sm.forward()

    if args.lhs_sample_size is not None:
        # Draw Latin hypercube samples of parameter values
        s = sm.lhs(siz=args.lhs_sample_size, seed=1000, save_output=save_output) # create sample set
    elif args.parstudy_divisions is not None :
        # Generate parameter study samples
        s = sm.parstudy(nvals=args.parstudy_divisions, save_output=save_output)  # create sample set

    # Run model using values in samples for parameter values
    s.run(cpus=args.ncpus, verbose=False)

    # Plot histograms, scatterplots, and correlation coefficients in paired matrix
    s.panels(figsize=(20, 20))
