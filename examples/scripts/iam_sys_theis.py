'''
Example illustrates system model containing only Theis reservoir model.
System model is run for an array of time points. Parameters of the reservoir
model are defined with a single injection well with time-varying rate
and post-injection period. Plot of injection rate and pressure vs time is
returned as output.

Examples of run:
$ python iam_sys_theis.py
'''

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, TheisReservoir


if __name__=='__main__':

    # Create system model
    time_array = np.array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.])
    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    rates = [5.0e-3, 6.0e-3, 7.0e-3, 8.0e-3, 0, 0, 0, 0, 0, 0, 0, 0,]
    tres = sm.add_component_model_object(
        TheisReservoir(name='tres', parent=sm, injX=100, injY=100,
                       locX=100, locY=45, injTimes=time_array,
                       injRates=rates))

    # Set input parameters
    tres.add_par('initialPressure', value=1.0e6)    # Pa
    tres.add_par('reservoirThickness', value=30)    # m
    tres.add_par('logResPerm', value=-10.69897)     # m^2
    tres.add_par('reservoirPorosity', value=.2)
    tres.add_par('brineDensity', value=1000)        # kg/m^3
    tres.add_par('brineViscosity', value=2.535e-3)  # Pa*s
    tres.add_par('CO2Density', value=800)           # kg/m^3
    tres.add_par('compressibility', value=2.46e-9)  # 1/Pa

    # Add observations of reservoir component model
    tres.add_obs('pressure')

    # Run system model using current values of its parameters
    sm.forward()

    # Assign observations of the model to pressure variables
    pressures = sm.collect_observations_as_time_series(tres, 'pressure')

    def twolines(x, y1, y2, label1, label2):
        _, ax = plt.subplots()

        # Plot linear sequence, and set tick labels to the same color
        ax.plot(x, y1, color='red')
        ax.tick_params(axis='y', labelcolor='red')
        ax.set_ylabel(label1, color='red')

        # Generate a new Axes instance, on the twin-X axes (same position)
        ax2 = ax.twinx()
        ax2.plot(x, y2, color='blue', linestyle='dashed')
        ax2.tick_params(axis='y', labelcolor='blue')
        ax2.set_ylabel(label2, color='blue')

        plt.tight_layout()
        plt.show()

    # plot injection rate and pressure vs. time
    twolines(time_array, rates, pressures, 'Injection Rate', 'Pressure')
