'''
Example illustrates system model containing four Theis reservoir models,
each with a different time-varying injection rate.  System model is run
for an array of time points. Pressure changes predicted by each reservoir
model are added to calculate resulting pressure at the observation location.

Examples of run:
$ python iam_sys_theis_4inj.py
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
    rates = []
    rates.append(np.array([5.0e-3, 6.0e-3, 7.0e-3, 8.0e-3, 0, 0, 0, 0, 0, 0, 0, 0,]))
    rates.append(-rates[0])
    rates.append(2*rates[0])
    rates.append(-2*rates[0])

    # Create four Theis well component models, each predicting pressure at the same location
    tress = []
    for i in range(num_wells):

        # Add reservoir component
        tress.append(sm.add_component_model_object(
            TheisReservoir(name='tres'+str(i), parent=sm,
                           injX=wx[i], injY=wy[i],
                           locX=500, locY=500,
                           injTimes=time_array, injRates=rates[i])))

        # Set input parameters
        tress[i].add_par('initialPressure', value=1.0e6)     # Pa
        tress[i].add_par('reservoirThickness', value=30)    # m
        tress[i].add_par('logResPerm', value=-10.69897)     # m^2
        tress[i].add_par('reservoirPorosity', value=.2)
        tress[i].add_par('brineDensity', value=1000)        # kg/m^3
        tress[i].add_par('brineViscosity', value=2.535e-3)  # Pa*s
        tress[i].add_par('CO2Density', value=800)           # kg/m^3
        tress[i].add_par('compressibility', value=2.46e-9)  # 1/Pa

        # Add observations of reservoir component model
        tress[i].add_obs('pressure')

    # Run system model using current values of its parameters
    sm.forward()

    # add pressure changes from all four wells
    pressures = sum(sm.collect_observations_as_time_series(tres, 'pressure') \
                    - 1.0e6 for i, tres in enumerate(tress)) + 1.0e6

    # print results
    for j, time in enumerate(time_array):
        print(rates[0][j], rates[1][j], rates[2][j], rates[3][j], pressures[j])

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
