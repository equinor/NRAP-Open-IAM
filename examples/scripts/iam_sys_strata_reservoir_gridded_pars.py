# -*- coding: utf-8 -*-
"""
This example couples stratigraphy and simple reservoir components.
Stratigraphy component defines gridded parameter shale1Thickness which is
then used to define a deterministic parameter of the same name
for the simple reservoir component.

Example of run:
$ python iam_sys_strata_reservoir_gridded_pars.py
"""
import os
import sys
import logging
import numpy as np

sys.path.insert(0,os.sep.join(['..','..','source']))
from openiam import DataInterpolator, Stratigraphy, SimpleReservoir

try:
    from openiam import SystemModel
except ImportError as err:
    print('Unable to load IAM class module: '+str(err))


if __name__ == "__main__":
    # Define logging level
    logging.basicConfig(level=logging.WARNING)

    # Define keyword arguments of the system model
    num_years = 5
    time_array = 365.25*np.arange(0.0,num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    int1 = DataInterpolator(name='int1', parent=sm,
                            header_file_dir=os.path.join(
                                '..', '..', 'source', 'components', 'reservoir',
                                'lookuptables', 'Test_2d'),
                            data_file='shale1Thickness.csv')
    strata = sm.add_component_model_object(Stratigraphy(name='strata', parent=sm))

    # Add parameters of reservoir component model
    strata.add_par('numberOfShaleLayers', value=3, vary=False)
    strata.add_gridded_par('shale1Thickness', int1)
    strata.add_par('shale2Thickness', min=40.0, max=60., value=50.0)

    # Add reservoir component
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm,
        injX=37200., injY=48000., locX=37478.0, locY=48333.0))  # chosen injX, injY are arbitrary

    # Add parameters of reservoir component model
    sres.add_par_linked_to_par('numberOfShaleLayers',
                         strata.deterministic_pars['numberOfShaleLayers'])
    sres.add_par('injRate', min=0.4, max=0.6, value=0.5)
    sres.add_par('shale1Thickness', vary=False,
        value=strata.gridded_pars['shale1Thickness'](np.array([[sres.locX,sres.locY]])))
    sres.add_par_linked_to_par('shale2Thickness', strata.pars['shale2Thickness'])
    sres.add_par_linked_to_par('shale3Thickness', strata.default_pars['shaleThickness'])

    sres.add_par_linked_to_par('aquifer1Thickness', strata.default_pars['aquiferThickness'])
    sres.add_par_linked_to_par('aquifer2Thickness', strata.default_pars['aquiferThickness'])
    sres.add_par_linked_to_par('aquifer3Thickness', strata.default_pars['aquiferThickness'])

    sres.add_par_linked_to_par('reservoirThickness', strata.default_pars['reservoirThickness'])

    sres.add_par_linked_to_par('datumPressure', strata.default_pars['datumPressure'])

    sres.add_obs('pressure')

    # Run system model using current values of its parameters
    sm.forward()

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    # Print pressure
    print('pressure', sm.collect_observations_as_time_series(sres,'pressure'), sep='\n')
