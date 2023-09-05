"""
Test example illustrates sampling of model parameters utilizing pyDOE.bbdesign
"""
import sys
import os
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from matk import matk, pyDOE

# Number of samples
numsamples = 480

# Create matk object
p = matk()

# Parameters
p.add_par('par1', min=30.0, max=90.0)
p.add_par('par2', min=-550.0, max=-700.0)
p.add_par('par3', min=0.02, max=0.2)
p.add_par('par4', min=-14.0, max=-11.0)
p.add_par('par4', min=0.0, max=3.0)
p.add_par('par5', min=0.0, max=1.0)
p.add_par('par6', min=-7.0, max=1.5)
p.add_par('par7', min=-7.0, max=1.5)

# Create sampleset using pyDOE and parameter mins and maxs
#s = p.parmins + pyDOE.lhs(len(p.pars), samples=numsamples) * (p.parmaxs-p.parmins)
#s = p.parmins + pyDOE.bbdesign(len(p.pars)) * (p.parmaxs-p.parmins)
s = (p.parmins + p.parmaxs)/2 + pyDOE.bbdesign(8, center=1)*(p.parmaxs-p.parmins)/2
print(s)
