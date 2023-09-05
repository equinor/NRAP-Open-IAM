# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 09:04:42 2022

@author: Veronika Vasylkivska
"""
import sys
import os
import logging
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import SystemModel, GenericReservoir

if __name__ == "__main__":
    # For multiprocessing in Spyder
    __spec__ = None

    logging.basicConfig(level=logging.WARNING)

    test_case = 1  # 1 is a forward run; 2 is an lhs run

    num_years = 5
    time_array = 365.25*np.arange(0, num_years+1)

    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    well_x = 250.0
    well_y = 350.0

    # Add reservoir component
    res = sm.add_component_model_object(GenericReservoir(
        name='res', parent=sm, injX=0., injY=0., locX=well_x, locY=well_y))

    # Add parameters of reservoir component model
    res.add_par('reservoirDepth', value=1500.0, vary=False)
    res.add_par('reservoirThickness', value=50.0, vary=False)
    res.add_par('logResPerm', value=-14.0, min=-14.5, max=-13.5)
    res.add_par('resTempGradient', value=30.0, vary=False)
    res.add_par('injRate', value=100, vary=False)
    res.add_par('initialSalinity', value=0.05, vary=False)
    res.add_par('wellRadius', value=0.05, vary=False)
    res.add_par('reservoirPorosity', value=0.1, vary=False)

    # Add observations (output) from the reservoir model
    res.add_obs('pressure')
    res.add_obs('CO2saturation')

    if test_case == 1:  # forward simulation
        # print('------------------------------------------------------------------')
        # print('                  Forward method illustration ')
        # print('------------------------------------------------------------------')
        sm.forward()
        # Collect results
        pressure = sm.collect_observations_as_time_series(res, 'pressure')
        saturation = sm.collect_observations_as_time_series(res, 'CO2saturation')
        print('pressure', pressure)

    else:
        num_samples = 10
        ncpus = 2
        # Draw Latin hypercube samples of parameter values
        s = sm.lhs(siz=num_samples, seed=123)

        # Run model using values in samples for parameter values
        results = s.run(cpus=ncpus, verbose=False)
