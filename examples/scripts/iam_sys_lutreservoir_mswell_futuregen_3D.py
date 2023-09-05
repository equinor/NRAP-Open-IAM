'''
This examples couples the stratigraphy, lookup table reservoir, open wellbore and
FutureGen aquifer components.  The saturation/pressure output produced by lookup table
reservoir model is used to drive leakage from a single open wellbore
model, which is passed to the input of an adapter that provides well
coordinates, CO2 and brine leakage rates and cumulative mass fluxes to the
futuregen aquifer model.  The stratigraphy component defines the thickness of the aquifer
and shale layers above the reservoir, which are then used as input parameters
for the open wellbore and aquifer components.

This example was created based on the code in iam_sys_lutreservoir_openwell_futuregen_aor.py.
Specifically, this example adds 3D interpolation, which was not included
in the original example. This script, however, does not incorporate
Latin Hypercube Sampling (LHS) or the determination of Area of Review (AoR).

This example has been made to be similar to the control file example in the file
ControlFile_ex_LUTR_3D_Interpolation.yaml.

This example requires the additional FutureGen 2.0 data set.
FutureGen 2.0 data set can be downloaded from the following source:
https://edx.netl.doe.gov/dataset/futuregen-2-0-1008-simulation-reservoir-lookup-table

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/FutureGen2/1008_sims

Author of iam_sys_lutreservoir_openwell_futuregen_aor.py: Diana Bacon
at diana.bacon@pnnl.gov
This example was created by editing the above file.
Author: Nate Mitchell at nathaniel.mitchell@netl.doe.gov

Last Modified: August 5th, 2022
'''


import os
import sys
import datetime
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, ReservoirDataInterpolator,
                     LookupTableReservoir, MultisegmentedWellbore,
                     FutureGen2Aquifer, FutureGen2AZMI,
                     RateToMassAdapter, Stratigraphy)

if __name__ == "__main__":

    # For multiprocessing in Spyder
    __spec__ = None

    file_directory = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                  'lookuptables', 'FutureGen2', '1008_sims'])

    if not os.path.exists(os.sep.join([file_directory, 'fg1.csv'])):
        url = ''.join([
            'https://edx.netl.doe.gov/dataset/',
            'futuregen-2-0-1008-simulation-reservoir-lookup-table \n'])
        msg = ''.join([
            '\nFutureGen 2.0 data set can be downloaded ',
            'from the following source:\n',
            url,
            'Check example description for more information.'])
        print(msg)

    # Define keyword arguments of the system model
    num_years = 20
    time_array = 365.25 * np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Define stratigraphy
    num_aquifers = 4
    aquifer_name = ['Ironton-Galesville', 'Potosi', 'New Richmond', 'St Peter']
    aqu_thick = [33.2, 84.1, 31.1, 61.6]
    por = [0.118, 0.038, 0.132, 0.18]
    log_permh = [-13.39, -11.05, -12.48, -11.92]
    log_aniso = [0.30, 1.00, 0.30, 0.30]
    rel_vol_frac_calcite = [0.1, 0.5, 0.1, 0.1]
    num_shales = 5
    shale = [198.7, 74.4, 110.3, 118.9, 530.4]

    # These locations are used in the example. The values were selected by examining
    # the distribution of z with x and y in the tables (e.g., fg1.csv).
    x_loc = [225000, 275000]
    y_loc = [4410000, 4410000]
    z_loc = [-1000.0, -1350.0]
    num_locs = len(x_loc)

    output_location = os.path.join(
        os.getcwd(), '..', '..', 'Output',
        'Script_LUTR_3D_Interp_' + str(datetime.date.today()))

    if os.path.isdir(output_location) is False:
        os.mkdir(output_location)

    # Read file with signatures of interpolators and names of files with the corresponding data
    signature_data = np.genfromtxt(
        os.sep.join([file_directory, 'parameters_and_filenames_trunc.csv']),
        delimiter=",", dtype='str')

    # The first row (except the last element) of the file contains names of the parameters
    par_names = signature_data[0, 1:-1]

    num_pars = len(par_names)

    num_interpolators = signature_data.shape[0]-1  # -1 since the first line is a header

    # Index used to select the parameters from the file parameters_and_filenames_trunc.csv
    index_ref = 5
    for i in range(num_locs):
        print('Starting {}'.format(i))
        x = x_loc[i]
        y = y_loc[i]
        z = z_loc[i]

        # Adjust strata based on model depth
        shale[4] = (- z - shale[0] - shale[1] - shale[2] - shale[3] - aqu_thick[0]
                    - aqu_thick[1] - aqu_thick[2] - aqu_thick[3])
        depth = 4*[None]  # list with four elements to be defined below
        depth[3] = shale[4] + aqu_thick[3]
        depth[2] = depth[3] + shale[3] + aqu_thick[2]
        depth[1] = depth[2] + shale[2] + aqu_thick[1]
        depth[0] = depth[1] + shale[1] + aqu_thick[0]

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Create and add interpolator to the system model
        signature = {par_names[j]: float(
            signature_data[index_ref, j+1]) for j in range(num_pars)}

        sm.add_interpolator(ReservoirDataInterpolator(
            name='int{}'.format(index_ref), parent=sm,
            header_file_dir=file_directory,
            time_file='time_points.csv', interp_2d=False,
            data_file='fg{}.csv'.format(index_ref),
            index=int(signature_data[index_ref, 0]),
            signature=signature), intr_family='reservoir')

        # Add reservoir component
        ltres = sm.add_component_model_object(LookupTableReservoir(
            name='ltres', parent=sm, intr_family='reservoir',
            locX=x, locY=y, locZ=z, interp_2d=False))

        # Add parameters of reservoir component model
        for j in range(num_pars):
            # add arbitrary line of values from signature_file
            ltres.add_par(
                par_names[j], value=float(signature_data[index_ref, j+1]), vary=False)

        # Add observations of reservoir component model
        ltres.add_obs('pressure')
        ltres.add_obs('CO2saturation')
        ltres.add_obs_to_be_linked('pressure')
        ltres.add_obs_to_be_linked('CO2saturation')

        # Add stratigraphy component
        strata = sm.add_component_model_object(
            Stratigraphy(name='strata', parent=sm))

        # Add parameters of stratigraphy component model
        strata.add_par('numberOfShaleLayers', value=5, vary=False)
        strata.add_par('shale1Thickness', value=shale[0], vary=False)
        strata.add_par('shale2Thickness', value=shale[1], vary=False)
        strata.add_par('shale3Thickness', value=shale[2], vary=False)
        strata.add_par('shale4Thickness', value=shale[3], vary=False)
        strata.add_par('shale5Thickness', value=shale[4], vary=False)
        strata.add_par('aquifer1Thickness', value=aqu_thick[0], vary=False)
        strata.add_par('aquifer2Thickness', value=aqu_thick[1], vary=False)
        strata.add_par('aquifer3Thickness', value=aqu_thick[2], vary=False)
        strata.add_par('aquifer4Thickness', value=aqu_thick[3], vary=False)
        strata.add_par('reservoirThickness', value=7, vary=False)

        # Add multisegmented wellbore component
        ms = sm.add_component_model_object(
            MultisegmentedWellbore(name = 'ms', parent = sm))

        # Add parameters of multisegmented wellbore component
        ms.add_par('wellRadius', value=0.05715, vary=False)
        ms.add_par('logAqu1Perm', value=-13.39, vary=False)
        ms.add_par('logAqu2Perm', value=-11.05, vary=False)
        ms.add_par('logAqu3Perm', value=-12.48, vary=False)
        ms.add_par('logAqu4Perm', value=-11.92, vary=False)
        ms.add_par('logWellPerm', value=-13.0, vary=False)
        ms.add_par('brineDensity', value=1030.9, vary=False)
        ms.add_par('CO2Density', value=775.0, vary=False)
        ms.add_par('brineViscosity', value=7.5e-4, vary=False)
        ms.add_par('CO2Viscosity', value=6.6e-5, vary=False)
        ms.add_par_linked_to_par('numberOfShaleLayers',
                                 strata.deterministic_pars['numberOfShaleLayers'])
        ms.add_par_linked_to_par('shale1Thickness',
                                 strata.deterministic_pars['shale1Thickness'])
        ms.add_par_linked_to_par('shale2Thickness',
                                 strata.deterministic_pars['shale2Thickness'])
        ms.add_par_linked_to_par('shale3Thickness',
                                 strata.deterministic_pars['shale3Thickness'])
        ms.add_par_linked_to_par('shale4Thickness',
                                 strata.deterministic_pars['shale4Thickness'])
        ms.add_par_linked_to_par('shale5Thickness',
                                 strata.deterministic_pars['shale5Thickness'])
        ms.add_par_linked_to_par('aquifer1Thickness',
                                 strata.deterministic_pars['aquifer1Thickness'])
        ms.add_par_linked_to_par('aquifer2Thickness',
                                 strata.deterministic_pars['aquifer2Thickness'])
        ms.add_par_linked_to_par('aquifer3Thickness',
                                 strata.deterministic_pars['aquifer3Thickness'])
        ms.add_par_linked_to_par('aquifer4Thickness',
                                 strata.deterministic_pars['aquifer4Thickness'])
        ms.add_par_linked_to_par('reservoirThickness',
                                 strata.deterministic_pars['reservoirThickness'])

        # Add keyword arguments linked to the output provided by reservoir model
        ms.add_kwarg_linked_to_obs('pressure', ltres.linkobs['pressure'])
        ms.add_kwarg_linked_to_obs('CO2saturation', ltres.linkobs['CO2saturation'])

        # Add observations of multisegmented wellbore component model
        ms.add_obs_to_be_linked('CO2_aquifer1')
        ms.add_obs_to_be_linked('brine_aquifer1')
        ms.add_obs_to_be_linked('CO2_aquifer2')
        ms.add_obs_to_be_linked('brine_aquifer2')
        ms.add_obs_to_be_linked('CO2_aquifer3')
        ms.add_obs_to_be_linked('brine_aquifer3')
        ms.add_obs_to_be_linked('CO2_aquifer4')
        ms.add_obs_to_be_linked('brine_aquifer4')
        ms.add_obs_to_be_linked('CO2_atm')
        ms.add_obs_to_be_linked('brine_atm')
        ms.add_obs('CO2_aquifer1')
        ms.add_obs('brine_aquifer1')
        ms.add_obs('CO2_aquifer2')
        ms.add_obs('brine_aquifer2')
        ms.add_obs('CO2_aquifer3')
        ms.add_obs('brine_aquifer3')
        ms.add_obs('CO2_aquifer4')
        ms.add_obs('brine_aquifer4')
        ms.add_obs('CO2_atm')
        ms.add_obs('brine_atm')

        # Add adapter that transforms leakage rates to accumulated mass
        adapt = sm.add_component_model_object(
            RateToMassAdapter(name='adapt', parent=sm))
        adapt.add_kwarg_linked_to_collection('CO2_aquifer1',
                                             [ms.linkobs['CO2_aquifer1'],
                                              ms.linkobs['CO2_aquifer2']])
        adapt.add_kwarg_linked_to_collection('CO2_aquifer2',
                                             [ms.linkobs['CO2_aquifer2'],
                                              ms.linkobs['CO2_aquifer3']])
        adapt.add_kwarg_linked_to_collection('CO2_aquifer3',
                                             [ms.linkobs['CO2_aquifer3'],
                                              ms.linkobs['CO2_aquifer4']])
        adapt.add_kwarg_linked_to_collection('CO2_aquifer4',
                                             [ms.linkobs['CO2_aquifer4'],
                                              ms.linkobs['CO2_atm']])
        adapt.add_kwarg_linked_to_collection('brine_aquifer1',
                                             [ms.linkobs['brine_aquifer1'],
                                              ms.linkobs['brine_aquifer2']])
        adapt.add_kwarg_linked_to_collection('brine_aquifer2',
                                             [ms.linkobs['brine_aquifer2'],
                                              ms.linkobs['brine_aquifer3']])
        adapt.add_kwarg_linked_to_collection('brine_aquifer3',
                                             [ms.linkobs['brine_aquifer3'],
                                              ms.linkobs['brine_aquifer4']])
        adapt.add_kwarg_linked_to_collection('brine_aquifer4',
                                             [ms.linkobs['brine_aquifer4'],
                                              ms.linkobs['brine_atm']])
        adapt.add_obs_to_be_linked('mass_CO2_aquifer4')
        adapt.add_obs_to_be_linked('mass_brine_aquifer4')
        adapt.add_obs('mass_CO2_aquifer4')
        adapt.add_obs('mass_brine_aquifer4')

        # Add futuregen aquifer model object and define parameters
        if float(depth[3]) < 700:
            fga = sm.add_component_model_object(
                FutureGen2Aquifer(name='fga', parent=sm))
        else:
            fga = sm.add_component_model_object(
                FutureGen2AZMI(name='fga', parent=sm))

        # St. Peter Sandstone
        fga.add_par('aqu_thick', value=aqu_thick[3], vary=False)
        fga.add_par('depth', value=depth[3], vary=False)
        fga.add_par('por', value=0.18, vary=False)
        fga.add_par('log_permh', value=-11.92, vary=False)
        fga.add_par('log_aniso', value=0.30, vary=False)
        fga.add_par('rel_vol_frac_calcite', value=0.1, vary=False)

        # Add aquifer component's keyword argument co2_rate
        # linked to the collection created above
        fga.add_kwarg_linked_to_obs(
            'co2_rate', ms.linkobs['CO2_aquifer4'])
        fga.add_kwarg_linked_to_obs(
            'brine_rate', ms.linkobs['brine_aquifer4'])
        fga.add_kwarg_linked_to_obs(
            'co2_mass', adapt.linkobs['mass_CO2_aquifer4'])
        fga.add_kwarg_linked_to_obs(
            'brine_mass', adapt.linkobs['mass_brine_aquifer4'])

        # Add observations (output) from the futuregen aquifer model
        fga.add_obs('pH_volume')
        fga.add_obs('TDS_volume')

        # Run system model using current values of its parameters
        sm.forward()  # system model is run deterministically

        print('Results for (x, y, z) of (' + str(x) + ' m, ' + str(y) + \
              ' m, ' + str(z) + ' m): ')

        # Get results
        # Convert pressure to MPa
        pressure = sm.collect_observations_as_time_series(ltres, 'pressure')/1.0e6
        pdiff = [p - pressure[0] for p in pressure]
        sat = sm.collect_observations_as_time_series(ltres, 'CO2saturation')
        CO2_aquifer4 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer4')
        brine_aquifer4 = sm.collect_observations_as_time_series(ms, 'brine_aquifer4')
        pH = sm.collect_observations_as_time_series(fga, 'pH_volume')
        TDS = sm.collect_observations_as_time_series(fga, 'TDS_volume')

        # Print results
        print('Pressure (MPa): ', pressure, sep='\n')
        print('Change in Pressure (MPa) over time relative to t = 0 yrs: ',
              pdiff, sep='\n')
        print('CO2_aquifer4 [kg/s]: ', CO2_aquifer4, sep='\n')
        print('brine_aquifer4 [kg/s]: ', brine_aquifer4, sep='\n')
        print('CO2saturation [-]: ', sat, sep='\n')
        print('pH_volume (m^3): ', pH, sep='\n')
        print('TDS_volume (m^3): ', TDS, sep='\n')

        # Plot results
        fig = plt.figure(1, figsize=(9, 4))
        ax = fig.add_subplot(1, 2, i + 1)
        plt.plot(time_array / 365.25, pressure, linewidth=1,
                 label='Location ' + str(i + 1))
        plt.legend(fancybox=False, framealpha=0.5)
        plt.xlabel('Time, t (years)', fontsize=14, fontweight='bold')
        plt.ylabel('Reservoir Pressure, P (MPa)', fontsize=14, fontweight='bold')
        plt.title('Location ' + str(i + 1), fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.tick_params(labelsize=12)
        plt.xlim([0, num_years])
        if i == 1:
            os.chdir(output_location)
            plt.savefig('Pressure_LUTR_3D_Interp.png', dpi=300)

        fig = plt.figure(2, figsize = (9, 4))
        ax = fig.add_subplot(1, 2, i + 1)
        plt.plot(time_array / 365.25, brine_aquifer4, linewidth=1, \
                 label='Location ' + str(i + 1))
        plt.legend(fancybox=False, framealpha=0.5)
        plt.xlabel('Time, t (years)', fontsize=14, fontweight='bold')
        plt.ylabel('Brine Leakage Rate\nto Aquifer 4, (kg/s)',
                   fontsize=14, fontweight='bold')
        plt.title('Location ' + str(i + 1), fontsize=14, fontweight='bold')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0),
                             useMathText=True)
        plt.tight_layout()
        plt.tick_params(labelsize=12)
        plt.xlim([0, num_years])
        if i == 1:
            plt.savefig('Brine_Leakage_LUTR_3D_Interp.png', dpi=300)
