'''
Example illustrates system model containing only Theis reservoir model.
The mathematical aquifer is confined, horizontal, with constant thickness,
homogeneous, isotropic, and extends to infinity. System model is run for
an array of time points. Parameters of the reservoir model are defined
with a single injection well with constant rate. Plot of injection rate
and pressure vs time is returned as output. Results are compared with
those of the STOMP-W simulator.

Examples of run:
$ python iam_sys_theis_benchmark1.py
'''

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, TheisReservoir


if __name__=='__main__':

    time_array = np.array([
        0.00000E+00, 1.66667E-01, 4.16667E-01, 7.91667E-01,
        1.35417E+00, 2.19792E+00, 3.46354E+00, 5.36198E+00,8.20964E+00,
        1.24811E+01, 1.88883E+01, 2.84992E+01, 4.29154E+01, 6.45398E+01,
        9.69764E+01, 1.20000E+02])/24.

    stomp_pressures = np.array([
        6.00000E+05, 6.07358E+05, 6.08329E+05,
        6.08911E+05, 6.09370E+05, 6.09773E+05, 6.10145E+05, 6.10498E+05,
        6.10841E+05, 6.11176E+05, 6.11507E+05, 6.11835E+05, 6.12161E+05,
        6.12486E+05, 6.12810E+05, 6.12990E+05])

    # Create system model
    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    tres = sm.add_component_model_object(
        TheisReservoir(name='tres', parent=sm,
                       injX=0, injY=0, locX=3.131, locY=0,
                       injTimes=time_array,
                       injRates=np.ones_like(time_array)*0.0115741))

    # Set input parameters
    tres.add_par('initialPressure', value=600000.)      # Pa
    tres.add_par('reservoirThickness', value=50)        # m
    tres.add_par('logResPerm', value=-10.6271524099)    # m^2
    tres.add_par('reservoirPorosity', value=.35)
    tres.add_par('brineDensity', value=9.98664E+02)     # kg/m^3
    tres.add_par('brineViscosity', value=1.01759E-03)   # Pa*s
    # Assume we're injecting water, not CO2
    tres.add_par('CO2Density', value=9.98664E+02)       # kg/m^3
    tres.add_par('compressibility', value=4.3531E-11)   # 1/Pa

    # Add observations of reservoir component model
    tres.add_obs('pressure')

    # Run system model using current values of its parameters
    sm.forward()

    # Assign observations of the model to pressure variables
    pressures = sm.collect_observations_as_time_series(tres, 'pressure')

    # Compare pressures predicted by Theis solution and STOMP-W
    fig, ax = plt.subplots()
    ax.plot(time_array, pressures, color='red', label='Theis')
    ax.scatter(time_array, stomp_pressures, color='blue', label='STOMP-W')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Pressure (Pa)')
    ax.legend()
    plt.tight_layout()
    plt.show()
