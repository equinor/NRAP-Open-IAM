'''
This example couples the analytical reservoir, multisegmented wellbore and
generic aquifer models. The saturation/pressure output produced by analytical
reservoir model is used to drive leakage from a single multisegmented wellbore
model, which is passed to the input of an adapter that provides well
coordinates, |CO2| and brine leakage rates and cumulative mass fluxes to the
generic aquifer model.

Example of run:
$ python iam_sys_analytical_mswell_generic.py
'''

import sys
import os
import numpy as np
sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, AnalyticalReservoir,
                     MultisegmentedWellbore, GenericAquifer, RateToMassAdapter)


if __name__ == "__main__":
    # For multiprocessing in Spyder
    __spec__ = None
    # Define keyword arguments of the system model
    num_years = 20
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}   # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # legacy well location
    xloc = 200
    yloc = 200

    # Add reservoir component
    ares = sm.add_component_model_object(AnalyticalReservoir(
        name='ares', parent=sm, injX=0., injY=0., locX=xloc, locY=yloc))

    # Add parameters of reservoir component model
    ares.add_par('numberOfShaleLayers', value=3, vary=False)
    ares.add_par('shale1Thickness', value=100.0, vary=False)
    ares.add_par('aquifer1Thickness', value=100.0, vary=False)
    ares.add_par('shale2Thickness', value=100.0, vary=False)
    ares.add_par('aquifer2Thickness', value=100.0, vary=False)
    ares.add_par('shale3Thickness', value=500.0, vary=False)
    ares.add_par('injRate', value=1.0, vary=False)

    # Add observations of reservoir component model
    ares.add_obs_to_be_linked('pressure')
    ares.add_obs_to_be_linked('CO2saturation')
    ares.add_obs('pressure')
    ares.add_obs('CO2saturation')
    ares.add_obs('mass_CO2_reservoir')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms',parent=sm))
    ms.add_par('logWellPerm', min=-14.0, max=-12.0, value=-13.0)

    # Add linked parameters: common to reservoir and wellbore components
    ms.add_par_linked_to_par(
        'numberOfShaleLayers', ares.deterministic_pars['numberOfShaleLayers'])
    ms.add_par_linked_to_par(
        'shale1Thickness', ares.deterministic_pars['shale1Thickness'])
    ms.add_par_linked_to_par(
        'shale2Thickness', ares.deterministic_pars['shale2Thickness'])
    ms.add_par_linked_to_par(
        'shale3Thickness', ares.deterministic_pars['shale3Thickness'])
    ms.add_par_linked_to_par(
        'aquifer1Thickness', ares.deterministic_pars['aquifer1Thickness'])
    ms.add_par_linked_to_par(
        'aquifer2Thickness', ares.deterministic_pars['aquifer2Thickness'])
    ms.add_par_linked_to_par(
        'reservoirThickness', ares.default_pars['reservoirThickness'])
    ms.add_par_linked_to_par(
        'datumPressure', ares.default_pars['datumPressure'])

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', ares.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', ares.linkobs['CO2saturation'])

    # Add observations of multisegmented wellbore component model
    ms.add_obs_to_be_linked('CO2_aquifer1')
    ms.add_obs_to_be_linked('CO2_aquifer2')
    ms.add_obs_to_be_linked('brine_aquifer1')
    ms.add_obs_to_be_linked('brine_aquifer2')
    ms.add_obs_to_be_linked('mass_CO2_aquifer1')
    ms.add_obs_to_be_linked('mass_CO2_aquifer2')
    ms.add_obs_to_be_linked('brine_atm')
    ms.add_obs_to_be_linked('CO2_atm')
    ms.add_obs('brine_aquifer1')
    ms.add_obs('CO2_aquifer1')

    # Add adapter that transforms leakage rates to accumulated mass
    adapt = sm.add_component_model_object(RateToMassAdapter(name='adapt', parent=sm))
    adapt.add_kwarg_linked_to_collection('CO2_aquifer1',
        [ms.linkobs['CO2_aquifer1'], ms.linkobs['CO2_aquifer2']])
    adapt.add_kwarg_linked_to_collection('CO2_aquifer2',
        [ms.linkobs['CO2_aquifer2'], ms.linkobs['CO2_atm']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer1',
        [ms.linkobs['brine_aquifer1'], ms.linkobs['brine_aquifer2']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer2',
        [ms.linkobs['brine_aquifer2'], ms.linkobs['brine_atm']])
    adapt.add_obs_to_be_linked('mass_CO2_aquifer1')
    adapt.add_obs_to_be_linked('mass_CO2_aquifer2')
    adapt.add_obs_to_be_linked('mass_brine_aquifer1')
    adapt.add_obs_to_be_linked('mass_brine_aquifer2')
    adapt.add_obs('mass_CO2_aquifer1')
    adapt.add_obs('mass_brine_aquifer1')
    adapt.add_obs('mass_CO2_aquifer2')
    adapt.add_obs('mass_brine_aquifer2')

    # Add generic aquifer model object and define parameters
    ga = sm.add_component_model_object(GenericAquifer(name='ga', parent=sm))
    ga.add_par_linked_to_par('aqu_thick', ares.deterministic_pars['aquifer1Thickness'])
    ga.add_composite_par('top_depth',
        expr='ares.shale2Thickness + ares.shale3Thickness' +
        '+ ares.aquifer2Thickness')
    ga.add_par('por', value=1.965259282453879763e-01, vary=False)
    ga.add_par('log_permh', value=-1.191464515905165555e+01, vary=False)
    ga.add_par('log_aniso', value=8.046003470121247947e-01, vary=False)
    ga.add_par('aquifer_salinity', value=1.267995018132549341e-02, vary=False)
    ga.add_par('reservoir_salinity', value=4.159677791928499679e-02, vary=False)
    ga.add_par('dissolved_salt_threshold', value=0.015, vary=False)
    ga.add_par('dissolved_co2_threshold', value=0.001, vary=False)

    ga.add_kwarg_linked_to_obs('co2_mass', adapt.linkobs['mass_CO2_aquifer1'])
    ga.add_kwarg_linked_to_obs('brine_mass', adapt.linkobs['mass_brine_aquifer1'])

    # Add observations (output) from the generic aquifer model
    ga.add_obs('Dissolved_salt_volume')
    ga.add_obs('Dissolved_CO2_volume')

    # Define output folder to keep data files with gridded observations
    output_dir = 'dream_data'

    # Add gridded observations of the aquifer component
    ga.add_grid_obs('r_coordinate', constr_type='matrix', output_dir=output_dir)
    ga.add_grid_obs('z_coordinate', constr_type='matrix', output_dir=output_dir)
    ga.add_grid_obs('Dissolved_CO2_mass_fraction', constr_type='matrix', output_dir=output_dir)
    ga.add_grid_obs('Dissolved_salt_mass_fraction', constr_type='matrix', output_dir=output_dir)

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically

    print('Pressure',
          sm.collect_observations_as_time_series(ares, 'pressure'), sep='\n')
    print('------------------------------------------------------------------')
    print('CO2saturation',
          sm.collect_observations_as_time_series(ares, 'CO2saturation'), sep='\n')
    print('------------------------------------------------------------------')
    print('CO2_aquifer1',
          sm.collect_observations_as_time_series(ms, 'CO2_aquifer1'), sep='\n')
    print('------------------------------------------------------------------')
    print('brine_aquifer1',
          sm.collect_observations_as_time_series(ms, 'brine_aquifer1'), sep='\n')
    print('------------------------------------------------------------------')
    print('mass_CO2_aquifer1',
          sm.collect_observations_as_time_series(adapt, 'mass_CO2_aquifer1'), sep='\n')
    print('------------------------------------------------------------------')
    print('mass_brine_aquifer1',
          sm.collect_observations_as_time_series(adapt, 'mass_brine_aquifer1'), sep='\n')
    print('------------------------------------------------------------------')
    print('Dissolved_CO2_volume',
          sm.collect_observations_as_time_series(ga, 'Dissolved_CO2_volume'), sep='\n')
    print('------------------------------------------------------------------')
    print('Dissolved_salt_volume',
          sm.collect_observations_as_time_series(ga, 'Dissolved_salt_volume'), sep='\n')

    # Expected output:
    # Pressure
    # [9318225. 9456895. 9456895. 9456895. 9456895. 9456895. 9456895. 9456895.
    #  9456895. 9456895. 9456895. 9456895. 9456895. 9456895. 9456895. 9456895.
    #  9456895. 9456895. 9456895. 9456895. 9456895.]
    # ------------------------------------------------------------------
    # CO2saturation
    # [0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]
    # ------------------------------------------------------------------
    # CO2_aquifer1
    # [0.00e+00 1.01e-04 9.97e-05 9.97e-05 9.97e-05 9.97e-05 9.97e-05 9.97e-05
    #  9.97e-05 9.97e-05 9.97e-05 9.97e-05 9.97e-05 9.97e-05 9.97e-05 9.97e-05
    #  9.97e-05 9.97e-05 9.97e-05 9.97e-05 9.97e-05]
    # ------------------------------------------------------------------
    # brine_aquifer1
    # [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    # ------------------------------------------------------------------
    # mass_CO2_aquifer1
    # [    0.    3185.88  5183.36  7180.9   9178.49 11176.1  13173.74 15171.41
    #  17169.09 19166.78 21164.5  23162.22 25159.96 27157.71 29155.47 31153.24
    #  33151.02 35148.8  37146.6  39144.4  41142.21]
    # ------------------------------------------------------------------
    # mass_brine_aquifer1
    # [0.   1.35 0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
    #  0.   0.   0.   0.   0.   0.   0.  ]
    # ------------------------------------------------------------------
    # Dissolved_CO2_volume
    # [    0.    2638.33  2638.33  2638.33  2968.13  2968.13  4061.25  4061.25
    #   5154.38  5154.38  5154.38  6247.51  8433.76  9526.89 12654.11 12654.11
    #  13747.23 12654.11 13747.23 13747.23 16950.2 ]
    # ------------------------------------------------------------------
    # Dissolved_salt_volume
    # [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
