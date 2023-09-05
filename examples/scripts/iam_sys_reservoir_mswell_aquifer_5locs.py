'''
This example couples the simple reservoir, multisegmented wellbore and
carbonate aquifer models. The saturation/pressure output produced by simple
reservoir model is used to drive leakage from five multisegmented wellbore
models, which is passed to the input of an adapter that provides well
coordinates, |CO2| and brine leakage rates and cumulative mass fluxes to the
carbonate aquifer model.

Example of run:
$ python iam_sys_reservoir_mswell_aquifer_5locs.py
'''

import sys
import os
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, SimpleReservoir, MultisegmentedWellbore,
                     CarbonateAquifer, RateToMassAdapter)
from matk import pyDOE


if __name__=='__main__':
    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Create randomly located leaky well locations
    # within box defined by xmin,xmax,ymin,ymax
    xymins = np.array([100., 350.])
    xymaxs = np.array([200., 400.])
    num_wells = 5
    well_xys = xymins + pyDOE.lhs(2,samples=num_wells)*(xymaxs-xymins)

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
        sress[-1].add_par('injRate', value=0.8, vary=False)
        sress[-1].add_par('shale1Thickness', min=30.0, max=50., value=40.0)
        sress[-1].add_par('aquifer1Thickness', value=300., vary=False)
        sress[-1].add_par('aquifer2Thickness', min=30.0, max=50., value=45.0)
        sress[-1].add_par('shale3Thickness', min=10.0, max=50., value=35.0)
        sress[-1].add_par('logResPerm', min=-13.,max=-11., value=-12.)

        # Add observations of reservoir component model to be used by the next component
        sress[-1].add_obs_to_be_linked('pressure')
        sress[-1].add_obs_to_be_linked('CO2saturation')

        # Add observations of reservoir component model
        sress[-1].add_obs('pressure')
        sress[-1].add_obs('CO2saturation')
        sress[-1].add_obs('mass_CO2_reservoir')

        # Add multisegmented wellbore component
        mss.append(sm.add_component_model_object(
            MultisegmentedWellbore(name='ms'+str(i), parent=sm)))
        # A lot of parameters of multisegmented wellbore component
        # are the same as for the reservoir component
        # Add parameters linked to the same parameters from reservoir model
        mss[-1].add_par_linked_to_par(
            'numberOfShaleLayers', sress[-1].deterministic_pars['numberOfShaleLayers'])
        mss[-1].add_par_linked_to_par(
            'shale1Thickness', sress[-1].pars['shale1Thickness'])
        mss[-1].add_par_linked_to_par(
            'shale2Thickness', sress[-1].default_pars['shaleThickness'])
        mss[-1].add_par_linked_to_par(
            'shale3Thickness', sress[-1].pars['shale3Thickness'])
        mss[-1].add_par_linked_to_par(
            'aquifer1Thickness', sress[-1].deterministic_pars['aquifer1Thickness'])
        mss[-1].add_par_linked_to_par(
            'aquifer2Thickness', sress[-1].pars['aquifer2Thickness'])

        # Add keyword arguments linked to the output provided by reservoir model
        mss[-1].add_kwarg_linked_to_obs(
            'pressure', sress[-1].linkobs['pressure'])
        mss[-1].add_kwarg_linked_to_obs(
            'CO2saturation', sress[-1].linkobs['CO2saturation'])
        mss[-1].add_obs('brine_aquifer1')
        mss[-1].add_obs('CO2_aquifer1')
        mss[-1].add_obs_to_be_linked('brine_aquifer1')
        mss[-1].add_obs_to_be_linked('CO2_aquifer1')
        mss[-1].add_obs_to_be_linked('brine_aquifer2')
        mss[-1].add_obs_to_be_linked('CO2_aquifer2')
        mss[-1].add_obs_to_be_linked('brine_atm')
        mss[-1].add_obs_to_be_linked('CO2_atm')

        # Add adapter that transforms leakage rates to accumullated mass.
        # Its keyword arguments are linked to a collection (list) of observations.
        adapts.append(sm.add_component_model_object(
            RateToMassAdapter(name='adapt'+str(i), parent=sm)))
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
    ca = sm.add_component_model_object(CarbonateAquifer(name='ca', parent=sm))
    ca.add_par('perm_var', min=0.017, max=1.89, value=0.9535)
    ca.add_par('corr_len', min=1.0, max=3.95, value=2.475)
    ca.add_par('aniso', min=1.1, max=49.1, value=25.1)
    ca.add_par('mean_perm', min=-13.8, max=-10.3, value=-12.05)
    ca.add_par_linked_to_par('aqu_thick', sress[0].deterministic_pars['aquifer1Thickness'])
    ca.add_par('hyd_grad', min=2.88e-4, max=1.89e-2, value=9.59e-3)
    ca.add_par('calcite_ssa', min=0, max=1.e-2, value=5.5e-03)
    ca.add_par('organic_carbon', min=0, max=1.e-2, value=5.5e-03)
    ca.add_par('benzene_kd', min=1.49, max=1.73, value=1.61)
    ca.add_par('benzene_decay', min=0.15, max=2.84, value=1.5)
    ca.add_par('nap_kd', min=2.78, max=3.18, value=2.98)
    ca.add_par('nap_decay', min=-0.85, max=2.04, value=0.595)
    ca.add_par('phenol_kd', min=1.21, max=1.48, value=1.35)
    ca.add_par('phenol_decay', min=-1.22, max=2.06, value=0.42)
    ca.add_par('cl', min=0.1, max=6.025, value=0.776)
    ca.model_kwargs['x'] = well_xys[:, 0]
    ca.model_kwargs['y'] = well_xys[:, 1]

    # Create collections of references to observations containing CO2 and brine leakage rates
    CO2_rate_obs_list = []
    brine_rate_obs_list = []
    for ind in range(5):   # 5 is a number of wellbores
        CO2_rate_obs_list.append(mss[ind].linkobs['CO2_aquifer1'])
        brine_rate_obs_list.append(mss[ind].linkobs['brine_aquifer1'])
    # Print created collections
    print('------------------------------------------------------------------')
    print(CO2_rate_obs_list)
    print(brine_rate_obs_list)

    # Add aquifer component's keyword argument co2_rate linked to the collection created above
    ca.add_kwarg_linked_to_collection('co2_rate', CO2_rate_obs_list)
    # Add aquifer component's keyword argument brine_rate linked to the collection created above
    ca.add_kwarg_linked_to_collection('brine_rate', brine_rate_obs_list)

    # Create collections of references to observations containing CO2 and brine accumulated masses
    CO2_mass_obs_list = []
    brine_mass_obs_list = []
    for ind in range(5):   # 5 is a number of wellbores
        CO2_mass_obs_list.append(adapts[ind].linkobs['mass_CO2_aquifer1'])
        brine_mass_obs_list.append(adapts[ind].linkobs['mass_brine_aquifer1'])
    # Print created collections
    print('------------------------------------------------------------------')
    print(CO2_mass_obs_list)
    print(brine_mass_obs_list)

    # Add aquifer component's keyword argument co2_mass linked to the collection created above
    ca.add_kwarg_linked_to_collection('co2_mass', CO2_mass_obs_list)
    # Add aquifer component's keyword argument brine_rate linked to the collection created above
    ca.add_kwarg_linked_to_collection('brine_mass', brine_mass_obs_list)

    # Print aquifer component keywrod arguments linked to collections
    print('------------------------------------------------------------------')
    print(ca.collection_linked_kwargs)

    # Add observations (output) from the carbonate aquifer model
    ca.add_obs('pH_volume')
    ca.add_obs('TDS_volume')

    # Run system model using current values
    sm.forward()

    # Print some of the observations
    for i, sres in enumerate(sress):
        pressure = sm.collect_observations_as_time_series(sres, 'pressure')
        CO2saturation = sm.collect_observations_as_time_series(sres, 'CO2saturation')
        print('------------------------------------------------------------------')
        print('Pressure at wellbore {}'.format(i+1), pressure, sep='\n')
        print('CO2 saturation at wellbore {}'.format(i+1), CO2saturation, sep='\n')

    print('------------------------------------------------------------------')
    print('pH:',
          sm.collect_observations_as_time_series(ca,'pH_volume'), sep='\n')
    print('------------------------------------------------------------------')
    print('TDS:',
          sm.collect_observations_as_time_series(ca,'TDS_volume'), sep='\n')
