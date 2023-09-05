"""
Example illustrates linking of two simple reservoir, multisegmented
wellbore, and rate-to-mass adapter components. The locations of the wellbores
are random. Observations obtained from all three types of components are printed.

Example of run:
$ python iam_sys_reservoir_mswell_adapter_2locs.py
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, SimpleReservoir,
                     MultisegmentedWellbore, RateToMassAdapter)
from matk import pyDOE


if __name__=='__main__':

    # Define keyword arguments of the system model
    num_years = 5
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Create randomly located leaky well locations
    # within box defined by xmin,xmax,ymin,ymax
    xymins = np.array([100., 450.])
    xymaxs = np.array([300., 500.])
    num_wells = 2
    well_xys = xymins + pyDOE.lhs(2, samples=num_wells)*(xymaxs-xymins)

    sress = []
    mss = []
    adapts = []
    for i, crds in enumerate(well_xys):
        # Add reservoir component
        sress.append(sm.add_component_model_object(
                SimpleReservoir(name='sres'+str(i), parent=sm,
                injX=0., injY=0., locX=crds[0], locY=crds[1])))

        # Add parameters of reservoir component model
        sress[-1].add_par('numberOfShaleLayers', value=3, vary=False)
        sress[-1].add_par('injRate', value=0.8, vary=False)
        sress[-1].add_par('shale1Thickness', min=30.0, max=50., value=40.0)
        sress[-1].add_par('aquifer1Thickness', min=20.0, max=60., value=50.0)
        sress[-1].add_par('aquifer2Thickness', min=30.0, max=50., value=45.0)
        sress[-1].add_par('logResPerm', min=-13.,max=-11., value=-12.)

        # Add observations of reservoir component model to be used by the next component
        sress[-1].add_obs_to_be_linked('pressure')
        sress[-1].add_obs_to_be_linked('CO2saturation')

        # Add observations of reservoir component model
        sress[-1].add_obs('pressure')
        sress[-1].add_obs('CO2saturation')

        # Add multisegmented wellbore component
        mss.append(sm.add_component_model_object(
            MultisegmentedWellbore(name='ms'+str(i), parent=sm)))
        # A lot of parameters of multisegmented wellbore component
        # are the same as for the reservoir component
        # Add parameters linked to the same parameters from reservoir model
        mss[-1].add_par_linked_to_par('numberOfShaleLayers',
                                      sress[-1].deterministic_pars['numberOfShaleLayers'])
        mss[-1].add_par_linked_to_par('shale1Thickness', sress[-1].pars['shale1Thickness'])
        mss[-1].add_par_linked_to_par('shale2Thickness', sress[-1].default_pars['shaleThickness'])
        mss[-1].add_par_linked_to_par('shale3Thickness', sress[-1].default_pars['shaleThickness'])
        mss[-1].add_par_linked_to_par('aquifer1Thickness', sress[-1].pars['aquifer1Thickness'])
        mss[-1].add_par_linked_to_par('aquifer2Thickness', sress[-1].pars['aquifer2Thickness'])

        # Add keyword arguments linked to the output provided by reservoir model
        mss[-1].add_kwarg_linked_to_obs('pressure', sress[-1].linkobs['pressure'])
        mss[-1].add_kwarg_linked_to_obs('CO2saturation', sress[-1].linkobs['CO2saturation'])
        mss[-1].add_obs('brine_aquifer1')
        mss[-1].add_obs('CO2_aquifer1')
        mss[-1].add_obs_to_be_linked('brine_aquifer1')
        mss[-1].add_obs_to_be_linked('CO2_aquifer1')
        mss[-1].add_obs_to_be_linked('brine_aquifer2')
        mss[-1].add_obs_to_be_linked('CO2_aquifer2')
        mss[-1].add_obs_to_be_linked('brine_atm')
        mss[-1].add_obs_to_be_linked('CO2_atm')

        # Add adapter that transforms leakage rates to accumullated mass
        # Its keyword arguments are linked to a collection (list) of observations
        adapts.append(sm.add_component_model_object(
            RateToMassAdapter(name='adapt'+str(i),parent=sm)))
        adapts[-1].add_kwarg_linked_to_collection('CO2_aquifer1',
            [mss[-1].linkobs['CO2_aquifer1'], mss[-1].linkobs['CO2_aquifer2']])
        adapts[-1].add_kwarg_linked_to_collection('CO2_aquifer2',
            [mss[-1].linkobs['CO2_aquifer2'], mss[-1].linkobs['CO2_atm']])
        adapts[-1].add_kwarg_linked_to_collection('brine_aquifer1',
            [mss[-1].linkobs['brine_aquifer1'], mss[-1].linkobs['brine_aquifer2']])
        adapts[-1].add_kwarg_linked_to_collection('brine_aquifer2',
            [mss[-1].linkobs['brine_aquifer2'], mss[-1].linkobs['brine_atm']])
        adapts[-1].add_obs('mass_CO2_aquifer1')
        adapts[-1].add_obs('mass_CO2_aquifer2')
        adapts[-1].add_obs('mass_brine_aquifer1')
        adapts[-1].add_obs('mass_brine_aquifer2')

    sm.forward()

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    print('====================== PRESSURE AND SATURATION ===================')
    print('')
    # Print pressure and saturation
    for i, sres in enumerate(sress):
        print('Pressure at wellbore {}'.format(i+1),
              sm.collect_observations_as_time_series(sres, 'pressure'), sep='\n')
        print('CO2 saturation at wellbore {}'.format(i+1),
              sm.collect_observations_as_time_series(sres, 'pressure'), sep='\n')
        print('-----------------------------------------')
    print('')

    print('======================== LEAKAGE RATES ===========================')
    print('')
    for i, ms in enumerate(mss):
        print('Leakage rates at wellbore {}'.format(i+1))
        print('CO2: aquifer 1',
              sm.collect_observations_as_time_series(ms, 'CO2_aquifer1'), sep='\n')
        print('Brine: aquifer 1',
              sm.collect_observations_as_time_series(ms, 'brine_aquifer1'), sep='\n')
        print('-----------------------------------------')
    print('')

    print('========================= ACCUMULATED MASS =======================')
    print('')
    for i, adapt in enumerate(adapts):
        print('Accumulated CO2 mass from wellbore {}'.format(i+1))
        print('aquifer 1', sm.collect_observations_as_time_series(
            adapt, 'mass_CO2_aquifer1'), sep='\n')
        print('aquifer 2', sm.collect_observations_as_time_series(
            adapt, 'mass_CO2_aquifer2'), sep='\n')
        print('-----------------------------------------')
        print('Accumulated brine mass from wellbore {}'.format(i+1))
        print('aquifer 1', sm.collect_observations_as_time_series(
            adapt, 'mass_brine_aquifer1'), sep='\n')
        print('aquifer 2', sm.collect_observations_as_time_series(
            adapt, 'mass_brine_aquifer2'), sep='\n')
        print('-----------------------------------------')
