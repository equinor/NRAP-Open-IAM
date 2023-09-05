'''
This example couples the simple reservoir, multisegmented wellbore and
carbonate aquifer models. The saturation/pressure output produced by simple
reservoir model is used to drive leakage from five multisegmented wellbore
models separated into two groups according to their properties. Output from
wellbore components is passed to the input of an adapter that requires well
coordinates, |CO2| and brine leakage rates, and cumulative mass of the fluids
to the carbonate aquifer model.

Example of run:
$ python iam_sys_reservoir_mswell_2aquifers_2plus3locs.py
'''

import sys
import os
import numpy as np

sys.path.insert(0,os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, SimpleReservoir, MultisegmentedWellbore,
                     CarbonateAquifer, RateToMassAdapter)



if __name__ == '__main__':
    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0,num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    num_wells = 5
    well_xys = np.array([[100, 100],
                         [450, 125],
                         [200, 100],
                         [450, 225],
                         [340, 50]])

    sress = []
    mss = []
    adapts = []
    for i,crds in enumerate(well_xys):
        # Add reservoir component
        sress.append(sm.add_component_model_object(
                SimpleReservoir(name='sres'+str(i), parent=sm,
                injX=0., injY=0., locX=crds[0], locY=crds[1])))

        # Add parameters of reservoir component model
        sress[-1].add_par('numberOfShaleLayers', value=3, vary=False)
        sress[-1].add_par('injRate', min=0.1, max=10.0, value=1.0)
        sress[-1].add_par('logResPerm', min=-12.5, max=-11., value=-11.5)
        sress[-1].add_par('shale1Thickness', min=10.0, max=25.0, value=15.0)
        sress[-1].add_par('shale2Thickness', min=10.0, max=25.0, value=25.0)
        sress[-1].add_par('shale3Thickness', min=300.0, max=500.0, value=350.0)
        sress[-1].add_par('aquifer1Thickness', min=100.0, max=150.0, value=120.0)
        sress[-1].add_par('aquifer2Thickness', min=100.0, max=140.0, value=110.0)
        sress[-1].add_par('reservoirThickness', min=30.0, max=50.0, value=40.0)
        sress[-1].add_par('brineResSaturation', min=0.1, max=0.3, value=0.2)

        # Add observations of reservoir component model to be used by the next component
        sress[-1].add_obs_to_be_linked('pressure')
        sress[-1].add_obs_to_be_linked('CO2saturation')

        # Add observations of reservoir component model
        sress[-1].add_obs('pressure')
        sress[-1].add_obs('CO2saturation')
#        sress[-1].add_obs('mass_CO2_reservoir')

        # Add multisegmented wellbore component
        mss.append(sm.add_component_model_object(MultisegmentedWellbore(
            name='ms'+str(i), parent=sm)))
        # Add parameters specific to the multisegmented wellbore
        if i < 2:
            mss[-1].add_par('logWellPerm', min=-14.0, max=-12.0, value=-12.5)
            mss[-1].add_par('wellRadius', min=0.01, max=0.02, value=0.015)
        else:
            mss[-1].add_par('logWellPerm', min=-13.0, max=-11.5, value=-12)
            mss[-1].add_par('wellRadius', min=0.015, max=0.018, value=0.016)
        # Add parameters linked to the same parameters from reservoir model
        mss[-1].add_par_linked_to_par('numberOfShaleLayers',
                                      sress[-1].deterministic_pars['numberOfShaleLayers'])
        mss[-1].add_par_linked_to_par('shale1Thickness', sress[-1].pars['shale1Thickness'])
        mss[-1].add_par_linked_to_par('shale2Thickness', sress[-1].pars['shale2Thickness'])
        mss[-1].add_par_linked_to_par('shale3Thickness', sress[-1].pars['shale3Thickness'])
        mss[-1].add_par_linked_to_par('aquifer1Thickness', sress[-1].pars['aquifer1Thickness'])
        mss[-1].add_par_linked_to_par('aquifer2Thickness', sress[-1].pars['aquifer2Thickness'])
        mss[-1].add_par_linked_to_par('reservoirThickness', sress[-1].pars['reservoirThickness'])
        mss[-1].add_par_linked_to_par('brineResSaturation', sress[-1].pars['brineResSaturation'])
        mss[-1].add_par_linked_to_par('datumPressure', sress[-1].default_pars['datumPressure'])

        # Add keyword arguments linked to the output provided by reservoir model
        mss[-1].add_kwarg_linked_to_obs('pressure', sress[-1].linkobs['pressure'])
        mss[-1].add_kwarg_linked_to_obs('CO2saturation', sress[-1].linkobs['CO2saturation'])

        # Add observations of wellbore component model to be used by the next component
        mss[-1].add_obs_to_be_linked('brine_aquifer1')
        mss[-1].add_obs_to_be_linked('CO2_aquifer1')
        mss[-1].add_obs_to_be_linked('brine_aquifer2')
        mss[-1].add_obs_to_be_linked('CO2_aquifer2')
        mss[-1].add_obs_to_be_linked('brine_atm')
        mss[-1].add_obs_to_be_linked('CO2_atm')

        # Add observations of wellbore component
        mss[-1].add_obs('brine_aquifer1')
        mss[-1].add_obs('CO2_aquifer1')
        mss[-1].add_obs('brine_aquifer2')
        mss[-1].add_obs('CO2_aquifer2')
        mss[-1].add_obs('brine_atm')
        mss[-1].add_obs('CO2_atm')

        # Add adapter that transforms leakage rates to accumullated mass.
        # Its keyword arguments are linked to a collection (list) of observations.
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
        adapts[-1].add_obs_to_be_linked('mass_CO2_aquifer1')
        adapts[-1].add_obs_to_be_linked('mass_CO2_aquifer2')
        adapts[-1].add_obs_to_be_linked('mass_brine_aquifer1')
        adapts[-1].add_obs_to_be_linked('mass_brine_aquifer2')

    # Add carbonate aquifer model object and define parameters
    cas = []
    for ind in range(2):
        cas.append(sm.add_component_model_object(
            CarbonateAquifer(name='ca'+str(ind),parent=sm)))
        cas[-1].add_par_linked_to_par(
            'aqu_thick', sress[0].pars['aquifer{}Thickness'.format(ind+1)])
        cas[-1].model_kwargs['x'] = well_xys[:,0]
        cas[-1].model_kwargs['y'] = well_xys[:,1]

        # Create collections of references to observations containing CO2 and brine leakage rates
        CO2_rate_obs_list = []
        brine_rate_obs_list = []
        for wind in range(5):   # 5 is a number of wellbores
            CO2_rate_obs_list.append(mss[wind].linkobs['CO2_aquifer{}'.format(ind+1)])
            brine_rate_obs_list.append(mss[wind].linkobs['brine_aquifer{}'.format(ind+1)])

        # Print created collections
        print('------------------------------------------------------------------')
        print('CO2 rate list', CO2_rate_obs_list)
        print('Brine rate list', brine_rate_obs_list)

        # Add aquifer component's keyword argument co2_rate linked
        # to the collection created above
        cas[-1].add_kwarg_linked_to_collection('co2_rate', CO2_rate_obs_list)
        # Add aquifer component's keyword argument brine_rate linked
        # to the collection created above
        cas[-1].add_kwarg_linked_to_collection('brine_rate', brine_rate_obs_list)

        # Create collections of references to observations containing
        # CO2 and brine accumulated masses
        CO2_mass_obs_list = []
        brine_mass_obs_list = []
        for wind in range(5):   # 5 is a number of wellbores
            CO2_mass_obs_list.append(adapts[wind].linkobs['mass_CO2_aquifer{}'.format(ind+1)])
            brine_mass_obs_list.append(adapts[wind].linkobs['mass_brine_aquifer{}'.format(ind+1)])
        # Print created collections
        print('------------------------------------------------------------------')
        print('CO2 mass list', CO2_mass_obs_list)
        print('Brine rate list', brine_mass_obs_list)

        # Add aquifer component's keyword argument co2_mass linked to the collection created above
        cas[-1].add_kwarg_linked_to_collection('co2_mass', CO2_mass_obs_list)
        # Add aquifer component's keyword argument brine_rate linked to the collection created above
        cas[-1].add_kwarg_linked_to_collection('brine_mass', brine_mass_obs_list)

        # Print aquifer component keywrod arguments linked to collections
        print('------------------------------------------------------------------')
        print(cas[-1].collection_linked_kwargs)

        # Add observations (output) from the carbonate aquifer model
        cas[-1].add_obs('pH_volume')
        cas[-1].add_obs('TDS_volume')

    # Run system model using current values
    sm.forward()

    to_print = True
    if to_print:
        # Print some of the observations
        for i,sres in enumerate(sress):
            pressure = sm.collect_observations_as_time_series(sres, 'pressure')
            CO2saturation = sm.collect_observations_as_time_series(sres, 'CO2saturation')
            print('------------------------------------------------------------------')
            print('Pressure at wellbore {}'.format(i+1), pressure, sep='\n')
            print('CO2 saturation at wellbore {}'.format(i+1), CO2saturation, sep='\n')

        for i,ms in enumerate(mss):
            CO2_aquifer1 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer1')
            brine_aquifer1 = sm.collect_observations_as_time_series(ms, 'brine_aquifer1')
            CO2_aquifer2 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer2')
            brine_aquifer2 = sm.collect_observations_as_time_series(ms, 'brine_aquifer2')
            print('------------------------------------------------------------------')
            print('CO2 leakage at wellbore {} to aquifer 1'.format(i+1),
                  CO2_aquifer1, sep='\n')
            print('Brine leakage at wellbore {} to aquifer 1'.format(i+1),
                  brine_aquifer1, sep='\n')
            print('------------------------------------------------------------------')
            print('CO2 leakage at wellbore {} to aquifer 2'.format(i+1),
                  CO2_aquifer2, sep='\n')
            print('Brine leakage at wellbore {} to aquifer 2'.format(i+1),
                  brine_aquifer2, sep='\n')


        print('------------------------------------------------------------------')
        print('pH:',
              sm.collect_observations_as_time_series(cas[0],'pH_volume'), sep='\n')
        print('------------------------------------------------------------------')
        print('TDS:',
              sm.collect_observations_as_time_series(cas[0],'TDS_volume'), sep='\n')

        print('------------------------------------------------------------------')
        print('pH:',
              sm.collect_observations_as_time_series(cas[1],'pH_volume'), sep='\n')
        print('------------------------------------------------------------------')
        print('TDS:',
              sm.collect_observations_as_time_series(cas[1],'TDS_volume'), sep='\n')
