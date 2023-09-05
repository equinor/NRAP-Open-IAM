'''
Example illustrates system model containing four Theis reservoir models,
each with a different time-varying injection rate.  System model is run
for an array of time points. Pressure changes predicted by each reservoir
model are added to calculate resulting pressure at the observation location.
This script is equivalent to the script iam_sys_theis_4inj.py with a single
reservoir component.


Examples of run:
$ python iam_sys_theis_4inj_wells.py
'''

import sys
import os
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, TheisReservoir


if __name__=='__main__':

    # Create system model
    time_array = np.array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.])
    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Specify x, y coordinates of four wells
    wx = [739, 584, 287, 519]
    wy = [423, 822, 333, 856]
    num_wells = len(wx)

    # Specify injection rates over time for four wells
    init_rate = np.array([5.0e-3, 6.0e-3, 7.0e-3, 8.0e-3, 0, 0, 0, 0, 0, 0, 0, 0])
    rates = np.zeros((4, len(init_rate)))
    rates[0, :] = init_rate
    rates[1, :] = -init_rate
    rates[2, :] = 2*init_rate
    rates[3, :] = -2*init_rate

    # Create a single Theis reservoir component model
    # Add reservoir component
    tres = sm.add_component_model_object(
        TheisReservoir(name='tres', parent=sm,
                       injX=wx, injY=wy,
                       locX=500, locY=500,
                       injTimes=time_array, injRates=rates))

    # Set input parameters
    tres.add_par('initialPressure', value=1.0e6)     # Pa
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

    # add pressure changes from all four wells
    pressure = sm.collect_observations_as_time_series(tres, 'pressure')

    # print results
    for j, time in enumerate(time_array):
        print(rates[0][j], rates[1][j], rates[2][j], rates[3][j], pressure[j])

    # Expected results
    # 0.005 -0.005 0.01 -0.01 1000000.0
    # 0.006 -0.006 0.012 -0.012 1002052.6816498041
    # 0.007 -0.007 0.014 -0.014 1002554.5023124588
    # 0.008 -0.008 0.016 -0.016 1003014.6444541777
    # 0.0 -0.0 0.0 -0.0 1003465.5585837264
    # 0.0 -0.0 0.0 -0.0 1000218.5325733812
    # 0.0 -0.0 0.0 -0.0 1000090.2336011328
    # 0.0 -0.0 0.0 -0.0 1000051.0244263352
    # 0.0 -0.0 0.0 -0.0 1000033.218697469
    # 0.0 -0.0 0.0 -0.0 1000023.4801213143
    # 0.0 -0.0 0.0 -0.0 1000017.5280396653
    # 0.0 -0.0 0.0 -0.0 1000013.6076707245
