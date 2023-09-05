'''
This example couples the stratigraphy, lookup table reservoir, open wellbore and
carbonate aquifer components. The saturation/pressure output produced by lookup table
reservoir model is used to drive leakage from a single open wellbore
model, which is passed to the input of an adapter that provides well
coordinates, |CO2| and brine leakage rates, and cumulative mass fluxes to the
carbonate aquifer model. The stratigraphy component defines the thickness of the aquifer
and shale layers above the reservoir, which are then used as input parameters
for the open wellbore and aquifer components.

This example requires the additional Kimberlina data set.
Kimberlina data set (pressure-in-pa-kimb-54-sims.zip or
Pressure_In_Pa_Kimb_54_sims.zip) can be downloaded from one of the following places:
1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam
2. https://gitlab.com/NRAP/Kimberlina_data

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/Kimb_54_sims

Example of run:
$ python iam_sys_strata_lutreservoir_openwell_aquifer.py
'''

import sys
import os
import logging
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, ReservoirDataInterpolator, LookupTableReservoir,
                     OpenWellbore, CarbonateAquifer, RateToMassAdapter, Stratigraphy)

if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)

    file_directory = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                  'lookuptables', 'Kimb_54_sims'])

    if not os.path.exists(os.sep.join([file_directory, 'Reservoir_data_sim01.csv'])):
        msg = ''.join(['\nKimberlina data set can be downloaded ',
                       'from one of the following places:\n',
                       '1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam \n',
                       '2. https://gitlab.com/NRAP/Kimberlina_data  \n',
                       'Check this example description for more information.'])
        logging.error(msg)

    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Read file with signatures of interpolators and names of files with the corresponding data
    signature_data = np.genfromtxt(
        os.path.join(file_directory, 'parameters_and_filenames.csv'),
        delimiter=",", dtype='str')

    # The first row (except the last element) of the file contains names of the parameters
    par_names = signature_data[0, 1:-1]
    num_pars = len(par_names)

    num_interpolators = signature_data.shape[0]-1  # -1 since the first line is a header

    # Create and add interpolators to the system model
    for ind in range(num_interpolators):

        signature = {par_names[j]: float(signature_data[ind+1, j+1]) for j in range(num_pars)}

        intpr = sm.add_interpolator(
            ReservoirDataInterpolator(
                name='int'+str(ind+1), parent=sm,
                header_file_dir=os.path.join('..', '..', 'source', 'components', 'reservoir',
                                             'lookuptables', 'Kimb_54_sims'),
                time_file='time_points.csv',
                data_file='Reservoir_data_sim{ind1:02}.csv'.format(ind1=ind+1),
                index=int(signature_data[ind+1, 0]),
                signature=signature,
                default_values={'salinity': 0.1, 'temperature': 50.0}),
            intr_family='reservoir')

        msg = 'Signature of the created interpolator is {}'.format(signature)
        logging.debug(msg)

    logging.debug('All interpolators are created')

    # Add reservoir component
    ltres = sm.add_component_model_object(
        LookupTableReservoir(name='ltres', parent=sm, intr_family='reservoir',
                             locX=37478.0, locY=48333.0))

    # Add parameters of reservoir component model
    for j in range(num_pars):
        # Add arbitrary line of values from signature_file
        ltres.add_par(par_names[j], value=float(signature_data[2, j+1]), vary=False)

    # Add observations of reservoir component model
    ltres.add_obs('pressure')
    ltres.add_obs('CO2saturation')
    ltres.add_obs_to_be_linked('pressure')
    ltres.add_obs_to_be_linked('CO2saturation')
    ltres.add_obs_to_be_linked('salinity')

    # Add stratigraphy component
    strata = sm.add_component_model_object(Stratigraphy(name='strata', parent=sm))

    # Add parameters of stratigraphy component model
    strata.add_par('numberOfShaleLayers', value=3, vary=False)
    strata.add_par('shale1Thickness', value=50.0, vary=False)
    strata.add_par('shale2Thickness', value=50.0, vary=False)
    strata.add_par('shale3Thickness', value=50.0, vary=False)
    strata.add_par('aquifer1Thickness', value=500.0, vary=False)
    strata.add_par('aquifer2Thickness', value=400.0, vary=False)

    # Add open wellbore component
    ow = sm.add_component_model_object(OpenWellbore(name='ow', parent=sm))

    # Add parameters of open wellbore component
    ow.add_par('wellRadius', min=0.001, max=0.002, value=0.0015)
    ow.add_par('logReservoirTransmissivity', min=-11.0, max=-9.0, value=-10.0)
    ow.add_par('logAquiferTransmissivity', min=-11.0, max=-9.0, value=-10.0)
    ow.add_par_linked_to_obs('brineSalinity', ltres.linkobs['salinity'])

    # Add keyword arguments of the open wellbore component model
    ow.add_kwarg_linked_to_obs('pressure', ltres.linkobs['pressure'])
    ow.add_kwarg_linked_to_obs('CO2saturation', ltres.linkobs['CO2saturation'])

    # Add composite parameters of open wellbore component
    ow.add_composite_par(
        'reservoirDepth', expr='+'.join(['strata.shale1Thickness',
                                         'strata.shale2Thickness',
                                         'strata.shale3Thickness',
                                         'strata.aquifer1Thickness',
                                         'strata.aquifer2Thickness']))
    ow.add_composite_par(
        'wellTop', expr='strata.shale3Thickness + strata.aquifer2Thickness')

    # Add observations of open wellbore component model
    ow.add_obs_to_be_linked('CO2_aquifer')
    ow.add_obs_to_be_linked('brine_aquifer')
    ow.add_obs_to_be_linked('brine_atm')
    ow.add_obs_to_be_linked('CO2_atm')
    ow.add_obs('CO2_aquifer')
    ow.add_obs('brine_aquifer')
    ow.add_obs('CO2_atm')     # zero since well top is in aquifer
    ow.add_obs('brine_atm')   # zero since well top is in aquifer

    # Add adapter that transforms leakage rates to accumulated mass
    adapt = sm.add_component_model_object(
        RateToMassAdapter(name='adapt', parent=sm))
    adapt.add_kwarg_linked_to_collection(
        'CO2_aquifer', [ow.linkobs['CO2_aquifer'], ow.linkobs['CO2_atm']])
    adapt.add_kwarg_linked_to_collection(
        'brine_aquifer', [ow.linkobs['brine_aquifer'], ow.linkobs['brine_atm']])
    adapt.add_obs_to_be_linked('mass_CO2_aquifer')
    adapt.add_obs_to_be_linked('mass_brine_aquifer')
    adapt.add_obs('mass_CO2_aquifer')
    adapt.add_obs('mass_brine_aquifer')

    # Add carbonate aquifer model object and define parameters
    ca = sm.add_component_model_object(
        CarbonateAquifer(name='ca', parent=sm))
    ca.add_par('perm_var', min=0.017, max=1.89, value=0.9535)
    ca.add_par('corr_len', min=1.0, max=3.95, value=2.475)
    ca.add_par('aniso', min=1.1, max=49.1, value=25.1)
    ca.add_par('mean_perm', min=-13.8, max=-10.3, value=-12.05)
    ca.add_par_linked_to_par('aqu_thick',
                             strata.deterministic_pars['aquifer2Thickness'])
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
    ca.model_kwargs['x'] = [37478.0]
    ca.model_kwargs['y'] = [37555.0]

    CO2_rate_obs_list = []
    brine_rate_obs_list = []
    CO2_mass_obs_list = []
    brine_mass_obs_list = []
    CO2_rate_obs_list.append(ow.linkobs['CO2_aquifer'])
    brine_rate_obs_list.append(ow.linkobs['brine_aquifer'])
    CO2_mass_obs_list.append(adapt.linkobs['mass_CO2_aquifer'])
    brine_mass_obs_list.append(adapt.linkobs['mass_brine_aquifer'])

    # Add aquifer component's keyword argument co2_rate linked to the collection created above
    ca.add_kwarg_linked_to_collection('co2_rate', CO2_rate_obs_list)
    ca.add_kwarg_linked_to_collection('brine_rate', brine_rate_obs_list)
    ca.add_kwarg_linked_to_collection('co2_mass', CO2_mass_obs_list)
    ca.add_kwarg_linked_to_collection('brine_mass', brine_mass_obs_list)

    # Add observations (output) from the carbonate aquifer model
    ca.add_obs('pH_volume')
    ca.add_obs('TDS_volume')

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically
    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    print('pressure',
          sm.collect_observations_as_time_series(ltres, 'pressure'), sep='\n')
    print('------------------------------------------------------------------')
    print('CO2saturation',
          sm.collect_observations_as_time_series(ltres, 'CO2saturation'), sep='\n')
    print('------------------------------------------------------------------')
    print('CO2_aquifer',
          sm.collect_observations_as_time_series(ow, 'CO2_aquifer'), sep='\n')
    print('------------------------------------------------------------------')
    print('CO2_atm',
          sm.collect_observations_as_time_series(ow, 'CO2_atm'), sep='\n')
    print('------------------------------------------------------------------')
    print('brine_aquifer',
          sm.collect_observations_as_time_series(ow, 'brine_aquifer'), sep='\n')
    print('------------------------------------------------------------------')
    print('brine_atm',
          sm.collect_observations_as_time_series(ow, 'brine_atm'), sep='\n')
    print('------------------------------------------------------------------')
    print('mass_CO2_aquifer',
          sm.collect_observations_as_time_series(adapt, 'mass_CO2_aquifer'), sep='\n')
    print('------------------------------------------------------------------')
    print('mass_brine_aquifer',
          sm.collect_observations_as_time_series(adapt, 'mass_brine_aquifer'), sep='\n')
    print('------------------------------------------------------------------')
    print('pH_volume', sm.collect_observations_as_time_series(ca, 'pH_volume'), sep='\n')
    print('------------------------------------------------------------------')
    print('TDS_volume', sm.collect_observations_as_time_series(ca, 'TDS_volume'), sep='\n')
