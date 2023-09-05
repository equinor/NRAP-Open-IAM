"""
This example shows setup of a system model containing simple reservoir and
multisegmented wellbore component. It also illustrates different variants
of use collect_observations_as_time_series method of the system model.
"""
import os
import sys
import logging
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import SystemModel, SimpleReservoir, MultisegmentedWellbore


if __name__ == "__main__":
    __spec__ = None

    logging.basicConfig(level=logging.WARNING)
    # Define keyword arguments of the system model
    num_years = 50
    num_obs = num_years+1
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))

    # Add parameters of reservoir component model
    sres.add_par('numberOfShaleLayers', value=3, vary=False)
    sres.add_par('injRate', value=0.1, vary=False)
    sres.add_par('reservoirThickness', value=40.0, vary=False)
    sres.add_par('shale1Thickness', min=30.0, max=150., value=55.0)
    sres.add_par('shale2Thickness', min=10.0, max=150., value=20.0)
    sres.add_par('shale3Thickness', min=10.0, max=150., value=25.0)
    sres.add_par('aquifer1Thickness', min=10.0, max=150., value=51.0)
    sres.add_par('aquifer2Thickness', min=10.0, max=150., value=28.0)

    # Add observations of reservoir component model
    sres.add_obs_to_be_linked('pressure')
    sres.add_obs_to_be_linked('CO2saturation')
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')
    sres.add_obs('mass_CO2_reservoir')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms',
                                                              parent=sm))
    ms.add_par('wellRadius', min=0.01, max=0.02, value=0.015)

    # Add linked parameters: common to both components
    ms.add_par_linked_to_par('numberOfShaleLayers',
                             sres.deterministic_pars['numberOfShaleLayers'])
    par_names = ['shale{}Thickness'.format(ind) for ind in range(1, 4)]+[
        'aquifer{}Thickness'.format(ind) for ind in range(1, 3)]
    for nm in par_names:
        ms.add_par_linked_to_par(nm, sres.pars[nm])
    ms.add_par_linked_to_par('reservoirThickness',
                             sres.deterministic_pars['reservoirThickness'])
    ms.add_par_linked_to_par('datumPressure',
                             sres.default_pars['datumPressure'])

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

    # Add observations of multisegmented wellbore component model
    obs_names = ['CO2_aquifer1', 'CO2_aquifer2',
                 'brine_aquifer1', 'brine_aquifer2',
                 'mass_CO2_aquifer1', 'mass_CO2_aquifer2']
    for nm in obs_names:
        ms.add_obs(nm)

    print('------------------------------------------------------------------')
    print('                          UQ illustration ')
    print('------------------------------------------------------------------')

    import random
    num_samples = 50
    ncpus = 1
    # Draw Latin hypercube samples of parameter values
    seed = random.randint(500, 1100)
    s = sm.lhs(siz=num_samples, seed=seed)

    # Run model using values in samples for parameter values
    s.run(cpus=ncpus, verbose=False)

    # Extract results from stochastic simulations
    outputs1 = s.collect_observations_as_time_series()

    # Another way to extract simulation results
    outputs2 = {}
    for nm in obs_names:
        outputs2['ms.'+nm] = s.collect_observations_as_time_series(
            cmpnt=ms, obs_nm=nm)['ms.'+nm]

    # Compare results
    for nm in obs_names:
        full_nm = 'ms.'+nm
        sum_diff = np.sum(np.abs(outputs1[full_nm]-outputs2[full_nm]))
        print('For observation {} the sum of differences is {}.'.format(
            nm, sum_diff))
