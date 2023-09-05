'''
Example illustrates system model containing only analytical reservoir model.
System model is run for an array of time points. Parameters of the reservoir
model are defined as default, deterministic and stochastic. Plots of pressure
and |CO2| saturation are returned as output.

Examples of run:
$ python iam_analytical_reservoir.py
'''

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, AnalyticalReservoir


if __name__ == '__main__':

    # Create system model
    time_array = np.arange(0, 201)   # in days
    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    sres = sm.add_component_model_object(AnalyticalReservoir(name='sres', parent=sm))

    # Add parameters of reservoir component model
    sres.add_par('numberOfShaleLayers', value=3, vary=False)
    sres.add_par('injRate', min=0.001, max=10., value=0.0185)
    sres.add_par('shale1Thickness', min=25.0, max=80., value=50.0)
    sres.add_par('shale2Thickness', min=35.0, max=90., value=45.0)

    # Add observations of reservoir component model
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')

    # Run system model using current values of its parameters
    sm.forward()

    # Assign observations of the model to pressure and CO2saturation variables
    pressure = sm.collect_observations_as_time_series(sres, 'pressure')
    CO2saturation = sm.collect_observations_as_time_series(sres, 'CO2saturation')

    # Print pressure and saturation
    print('pressure', pressure, sep='\n')
    print('CO2 saturation', CO2saturation, sep='\n')

    # Plot pressure and saturation
    plt.figure(1)
    plt.plot(time_array, pressure/1.0e+6,
             'c-', linewidth=2, label='pressure')
    plt.xlabel('Time, t (years)')
    plt.ylabel('Pressure, P (MPa)')

    plt.figure(2)
    plt.plot(time_array, CO2saturation,
             'r-', linewidth=2, label='saturation')
    plt.xlabel('Time, t (years)')
    plt.ylabel('Saturation, S')
