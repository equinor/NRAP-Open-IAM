'''
This example illustrates linking of simple reservoir, open wellbore and
atmosphere models. The saturation/pressure output produced by several simple
reservoir components is used to drive leakage from open wellbores.
|CO2| leakage rates are passed to the atmosphere model.

This example also illustrates how to save all the outputs produced by the simulation.
Changing variable 'save_output' at the beginning of simulation (line 32)
from True to False allows to cancel saving of the outputs.

Example of run:
$ python iam_sys_reservoir_openwell_atmosphere_5locs_lhs.py
'''

import sys
import os
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import SystemModel, SimpleReservoir, OpenWellbore, AtmosphericROM
from matk import pyDOE


if __name__ == '__main__':
    # For multiprocessing in Spyder
    __spec__ = None

    # Change the variable value to False if saving the outputs is not needed.
    # By default, the results will be saved in the folder 'output/csv_files' within root
    # folder of NRAP-Open-IAM
    save_output = True

    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Create randomly located leaky well locations
    # within box defined by xmin,xmax,ymin,ymax
    xymins = np.array([200., 550.])
    xymaxs = np.array([600., 700.])
    num_wells = 5
    well_xys = xymins + pyDOE.lhs(2, samples=num_wells)*(xymaxs-xymins)

    sress = []
    mss = []
    ow = []
    co2_leakrates_collector = []
    for i, crds in enumerate(well_xys):

        # Add reservoir components
        sress.append(sm.add_component_model_object(
            SimpleReservoir(name='sres'+str(i), parent=sm,
                            injX=0., injY=0., locX=crds[0], locY=crds[1])))

        # Add parameters of reservoir component model
        sress[-1].add_par('numberOfShaleLayers', value=3, vary=False)
        sress[-1].add_par('injRate', value=0.05, vary=False)
        if i == 0:
            sress[-1].add_par('shale1Thickness', min=300.0, max=500., value=400.0)
            sress[-1].add_par('aquifer1Thickness', min=20.0, max=60., value=50.0)
            sress[-1].add_par('shale2Thickness', min=700.0, max=900., value=800.0)
            sress[-1].add_par('aquifer2Thickness', min=60.0, max=80., value=75.0)
            sress[-1].add_par('shale3Thickness', min=150.0, max=250., value=200.0)
            sress[-1].add_par('logResPerm', min=-12.5, max=-11.5, value=-12.)
        else:
            sress[-1].add_par_linked_to_par(
                'shale1Thickness', sress[0].pars['shale1Thickness'])
            sress[-1].add_par_linked_to_par(
                'aquifer1Thickness', sress[0].pars['aquifer1Thickness'])
            sress[-1].add_par_linked_to_par(
                'shale2Thickness', sress[0].pars['shale2Thickness'])
            sress[-1].add_par_linked_to_par(
                'aquifer2Thickness', sress[0].pars['aquifer2Thickness'])
            sress[-1].add_par_linked_to_par(
                'shale3Thickness', sress[0].pars['shale3Thickness'])
            sress[-1].add_par_linked_to_par(
                'logResPerm', sress[0].pars['logResPerm'])

        # Add observations of reservoir component model to be used by the next component
        sress[-1].add_obs_to_be_linked('pressure')
        sress[-1].add_obs_to_be_linked('CO2saturation')

        # Add observations of reservoir component model
        sress[-1].add_obs('pressure')
        sress[-1].add_obs('CO2saturation')

        # Add open wellbore component
        ow.append(sm.add_component_model_object(
            OpenWellbore(name='ow'+str(i), parent=sm)))

        # Add parameters of open wellbore component
        ow[-1].add_par('wellRadius', min=0.025, max=0.035, value=0.03)
        ow[-1].add_par('logReservoirTransmissivity', min=-11.0, max=-9.0, value=-10.0)
        ow[-1].add_par('logAquiferTransmissivity', min=-11.0, max=-9.0, value=-10.0)
        ow[-1].add_par('brineSalinity', value=0.1, vary=False)
        ow[-1].add_par('wellTop', value=0.0, vary=False)

        # Add keyword arguments of the open wellbore component model
        ow[-1].add_kwarg_linked_to_obs('pressure', sress[-1].linkobs['pressure'])
        ow[-1].add_kwarg_linked_to_obs('CO2saturation', sress[-1].linkobs['CO2saturation'])

        # Add composite parameter of open wellbore component
        ow[-1].add_composite_par(
            'reservoirDepth', expr='+'.join(['sres0.shale1Thickness',
                                             'sres0.shale2Thickness',
                                             'sres0.shale3Thickness',
                                             'sres0.aquifer1Thickness',
                                             'sres0.aquifer2Thickness']))

        # Add observations of open wellbore component model
        ow[-1].add_obs_to_be_linked('CO2_atm')
        ow[-1].add_obs('CO2_aquifer')    # zero since well top is in atm
        ow[-1].add_obs('brine_aquifer')  # zero since well top is in aquifer
        ow[-1].add_obs('CO2_atm')
        ow[-1].add_obs('brine_atm')

        # Create collector for leakrates
        co2_leakrates_collector.append(ow[-1].linkobs['CO2_atm'])

    # Add Atm ROM component
    satm = sm.add_component_model_object(AtmosphericROM(name='satm', parent=sm))
    satm.add_par('T_amb', min=13.0, max=15, value=18.0, vary=False)
    satm.add_par('P_amb', min=0.99, max=1.02, value=1.01E+00, vary=True)
    satm.add_par('V_wind', min=3.0, max=8.0, value=5.0, vary=False)
    satm.add_par('C0_critical', value=0.01, vary=False)

    satm.model_kwargs['x_receptor'] = [200, 250, 400, 350]
    satm.model_kwargs['y_receptor'] = [580, 660, 525, 600]

    satm.model_kwargs['x_coor'] = well_xys[:, 0]
    satm.model_kwargs['y_coor'] = well_xys[:, 1]

    satm.add_kwarg_linked_to_collection('co2_leakrate', co2_leakrates_collector)

    # Add observations for receptors
    for i, r in enumerate(satm.model_kwargs['x_receptor']):
        satm.add_obs('outflag_r{0:03}'.format(i))

    n_sources = len(satm.model_kwargs['x_coor'])
    for i in range(n_sources):
        satm.add_obs('x_new_s{0:03}'.format(i))
        satm.add_obs('y_new_s{0:03}'.format(i))
        satm.add_obs('critical_distance_s{0:03}'.format(i))

    satm.add_obs('num_sources')

    # There is no forward method call as the corresponding results are not of interest
    print('Latin Hypercube sampling study started...')
    import random
    num_samples = 50
    ncpus = 4
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=random.randint(500, 1100), save_output=save_output)

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)
    print(results)
    print('------------------------------------------------------------------')
    print('                  End ')
    print('------------------------------------------------------------------')
