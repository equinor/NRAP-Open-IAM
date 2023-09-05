'''
Example illustrates system model containing a single simple component model
utilizing simple fortran code compiled into library. This example illustrates
a non-intrusive way (from ROM developer's view point) to introduce the new component model
to the IAM framework, using tools already available in iam_base_classes.py

The code takes in several optional arguments which can be changed if needed; the code
can also be run with default values:

lhs      - number of samples for Latin Hypercube Sampling mode
parstudy - number of parameter partitions
ncpus    - number of processors to run concurrent simulations

Examples of run:
$ python iam_sys_cmpnt1_fortran.py --lhs 10 --ncpus 2
$ python iam_sys_cmpnt1_fortran.py --parstudy 10 --ncpus 2
'''

import sys,os
import argparse     # argparse module makes it easy to write command-line interfaces.
import numpy as np
import random

sys.path.insert(0,os.sep.join(['..', '..', 'source']))
from openiam import SystemModel

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
    dest='lhs_sample_size', default=50,
    help='Latin Hypercube Sampling mode: enter number of samples')
group.add_argument('--parstudy', type=int,
    dest='parstudy_divisions', default=None,
    help='Parameter study mode: enter number of parameter partitions')

# Parameter ncpus has default value 1, i.e., by default, all simulations are run
# sequentially. If a user wants to run simulations in parallel, he/she has to
# specify the number of cpus to use.
parser.add_argument('--ncpus', type=int, default=1,
    help='Number of processors to use to run concurrent simulations')

args = parser.parse_args()


if __name__=='__main__':
    # For multiprocessing in Spyder
    __spec__ = None
    from iam_simple_models import fortran_based_model
    from sys import platform
    import os.path

    # Setup library and needed function names
    if platform in ("linux", "linux2"):
        # linux
        library = "example_rom.so"
    elif platform == "darwin":
        # OS X
        library = "example_rom.dylib"
    elif platform == "win32":
        # Windows...
        library = "example_rom.dll"
    if not os.path.isfile(library):
        print('\nERROR: Library {} do not exist.'.format(library)+
              ' Library has to be compiled before running this example. \n')
        sys.exit()

    # Create system model
    sm = SystemModel()

    # Add a component model (fortran_based_model)
    input_array = np.array([1.1, 2.3, 3.5, 4.7, 5.9, 6.1, 7.3, 8.5])
    model_kwargs = {'input_array': input_array}
    fb = sm.add_component_model(name='fb', model=fortran_based_model,
                                model_kwargs=model_kwargs)

    # Add parameters of the component model
    fb.add_par('var1', min=-5, max=5, value=1.5)
    fb.add_par('var2', min=2, max=10, value=4.3)

    # Add observations
    fb.add_obs('out_var1')

    # if one of the command line arguments was lhs_sample_size
    if args.lhs_sample_size is not None:
        # Draw Latin hypercube samples of parameter values
        s = sm.lhs(siz=args.lhs_sample_size,seed=random.randint(500, 1100))   # create sample set
    elif args.parstudy_divisions is not None:
        # Generate parameter study samples
        s = sm.parstudy(nvals=args.parstudy_divisions)  # create sample set

    # Run model using values in samples for parameter values
    s.run(cpus=args.ncpus, verbose=False)

    # Plot histograms, scatterplots, and correlation coefficients in paired matrix
    s.panels(figsize=(10, 10))