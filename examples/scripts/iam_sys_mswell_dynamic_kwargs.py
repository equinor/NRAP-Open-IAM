"""
This example illustrates the use of time series for dynamic keyword arguments
of multisegmented wellbore component.

Example of run:
$ python iam_sys_mswell_dynamic_kwargs.py
"""

import sys
import os
import random
import numpy as np


sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel


if __name__ == "__main__":
    # For multiprocessing in Spyder
    __spec__ = None
    from openiam import MultisegmentedWellbore

    # Define keyword arguments of the system model
    num_years = 5
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))

    # Add parameters of multisegmented wellbore component
    ms.add_par('numberOfShaleLayers', value=3, vary=False)
    ms.add_par('shale1Thickness', min=30.0, max=50., value=40.0)
    ms.add_par('aquifer1Thickness', min=20.0, max=60., value=50.0)
    ms.add_par('aquifer2Thickness', min=30.0, max=50., value=45.0)
    ms.add_par('logWellPerm', min=-14., max=-11., value=-13.5)

    # Add time series that will be used for dynamic keyword arguments
    pressure_time_series = [27119928.57, 33166857.14, 34104428.57, 34160928.57,
                            34088642.86, 33982357.14]
    saturation_time_series = [0., 0.34310857, 0.45430571, 0.47404286,
                              0.48760786, 0.50069071]

    # Add dynamic keyword arguments
    ms.add_dynamic_kwarg('pressure', pressure_time_series)
    ms.add_dynamic_kwarg('CO2saturation', saturation_time_series)

    # Specify indices of the time points of interest
    ms.add_obs('brine_aquifer1', index=[1, 3, 5])
    ms.add_obs('CO2_aquifer1', index=[1, 3, 5])

    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=100, seed=random.randint(500, 1100))   # create sample set

    # Run model using values in samples for parameter values
    s.run(cpus=5, verbose=False)

    # Plot histograms, scatterplots, and correlation coefficients in paired matrix
    s.panels(figsize=(25, 25))
