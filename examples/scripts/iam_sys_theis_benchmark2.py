'''
Example illustrates system model containing only Theis reservoir model.
The mathematical aquifer is confined, horizontal, with constant thickness,
homogeneous, isotropic, and extends to infinity. System model is run for
an array of time points. Parameters of the reservoir model are defined
with a single injection well with time-varying rate. Plot of injection rate
and pressure vs time is returned as output. Results are compared with
those of the STOMP-W simulator.

Examples of run:
$ python iam_sys_theis_benchmark2.py
'''

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, TheisReservoir


if __name__=='__main__':

    # Benchmark data
    time_array = np.array([
        0.,   6.,   12.,  18.,  24.,  30.,  36.,  42.,  48.,  54.,  60.,
        66.,  72.,  78.,  84.,  90.,  96.,  102., 108., 114., 120., 126.,
        132., 138., 144., 150., 156., 162., 168., 174., 180., 186., 192.,
        198., 204., 210., 216., 222., 228., 234., 240.])/24.
    stomp_pressures = np.array([
        600000., 610215., 611011., 611408., 611673., 622087., 623042.,
        623572., 623950., 644679., 646518., 647525., 648241., 689667.,
        693318., 695310., 696723., 616110., 610664., 608272., 606846.,
        605875., 605163., 604615., 604178., 603820., 603520., 603266.,
        603047., 602856., 602688., 602538., 602404., 602284., 602174.,
        602074., 601982., 601898., 601820., 601747., 601679.])
    rates = np.array([
        0.0115741, 0.0115741, 0.0115741, 0.0115741, 0.0231482, 0.0231482,
        0.0231482, 0.0231482, 0.0462964, 0.0462964, 0.0462964, 0.0462964,
        0.0925928, 0.0925928, 0.0925928, 0.0925928, 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       ])

    # Create system model
    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    tres = sm.add_component_model_object(
        TheisReservoir(name='tres', parent=sm,
                       injX=0, injY=0, locX=3.131, locY=0,
                       injTimes=time_array, injRates=rates))

    # Set input parameters
    tres.add_par('initialPressure', value=600000.)        # Pa
    tres.add_par('reservoirThickness', value=50)          # m
    tres.add_par('logResPerm', value=-10.6271524099)      # m^2
    tres.add_par('reservoirPorosity', value=.35)
    tres.add_par('brineDensity', value=9.98664E+02)       # kg/m^3
    tres.add_par('brineViscosity', value=1.01759E-03)     # Pa*s
    # Assume we're injecting water, not CO2
    tres.add_par('CO2Density', value=9.98664E+02)         # kg/m^3
    tres.add_par('compressibility', value=4.3531E-11)     # 1/Pa

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
