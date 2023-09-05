'''
NRAP IAM AOR FutureGen Area of Review (AoR)
This examples couples the strata, lookup table reservoir, open wellbore and
FutureGen2 aquifer components.  The saturation/pressure output produced by lookup table
reservoir model is used to drive leakage from a single open wellbore
model, which is passed to the input of an adapter that provides well
coordinates, CO2 and brine leakage rates and cumulative mass fluxes to the
FutureGen2 aquifer model.  The strata component defines the thickness of the aquifer
and shale layers above the reservoir, which are then used as input parameters
for the open wellbore and aquifer components. The risk-based Area of Review (AoR) is defined
based on the probability of whether there is an impact to groundwater at each grid location.

This example requires the additional FutureGen 2.0 data set.
FutureGen 2.0 data set can be downloaded from the following source:
https://edx.netl.doe.gov/dataset/futuregen-2-0-1008-simulation-reservoir-lookup-table

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/FutureGen2/1008_sims

Usage examples:
$ python iam_sys_lutreservoir_openwell_futuregen_aor_3d.py --run 1
'''
# @author: Diana Bacon
# diana.bacon@pnnl.gov
# Edited by Nate Mitchell
# Nathaniel.Mitchell@netl.doe.gov

import os
import sys
import argparse
import datetime
import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, ReservoirDataInterpolator,
                     LookupTableReservoir, OpenWellbore,
                     FutureGen2Aquifer, FutureGen2AZMI,
                     RateToMassAdapter, Stratigraphy)

if __name__ == "__main__":
    # For multiprocessing in Spyder
    __spec__ = None

    use_lutr_locations_option = True

    if use_lutr_locations_option:
        use_location_range_option = False
    else:
        use_location_range_option = True

    Scenario_name = 'fg1'

    # These values are used if use_location_range_option is True
    xmin = 5000
    xmax = 55000
    dx = 5000

    ymin = 5000
    ymax = 55000
    dy = 5000

    file_directory = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                  'lookuptables', 'FutureGen2', '1008_sims'])

    # Input arguments
    parser = argparse.ArgumentParser(description='Calculate risk-based AoR')
    parser.add_argument('--run', default='1', help='run number to process')
    args = parser.parse_args()
    run = int(args.run)

    # Define keyword arguments of the system model
    num_years = 20
    time_array = 365*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Define stratigraphy
    num_aquifers = 4
    aquifer_name = ['Ironton-Galesville', 'Potosi', 'New Richmond', 'St Peter']
    aq_thick = [33.2, 84.1, 31.1, 61.6]
    por = [0.118, 0.038, 0.132, 0.18]
    log_permh = [-13.39, -11.05, -12.48, -11.92]
    log_aniso = [0.30, 1.00, 0.30, 0.30]
    rel_vol_frac_calcite = [0.1, 0.5, 0.1, 0.1]
    num_shales = 5
    shale = [198.7, 74.4, 110.3, 118.9, 530.4]


    if use_lutr_locations_option:
        # Get grid locations from selected data file
        locations = np.genfromtxt(
            os.path.join(file_directory, Scenario_name + '.csv'),
            delimiter=',', skip_header=1, dtype=None, usecols=(0, 1, 2))

        x_loc = locations[:, 0]
        y_loc = locations[:, 1]
        z_loc = locations[:, 2]

        x_loc_redo = np.array([])
        y_loc_redo = np.array([])
        z_loc_redo = np.array([])

        # First, trim x, y, and z by the x and y boundaries provided
        min_x_loc_m = 1
        min_y_loc_m = 1
        max_x_loc_m = 9.99e+10
        max_y_loc_m = 9.99e+10
        for loc_ref in range(0, len(x_loc)):
            if x_loc[loc_ref] >= min_x_loc_m  and x_loc[loc_ref] <= max_x_loc_m \
                and y_loc[loc_ref] >= min_y_loc_m and y_loc[loc_ref] <= max_y_loc_m:
                x_loc_redo = np.append(x_loc_redo, x_loc[loc_ref])
                y_loc_redo = np.append(y_loc_redo, y_loc[loc_ref])
                z_loc_redo = np.append(z_loc_redo, z_loc[loc_ref])

        x_loc = x_loc_redo
        y_loc = y_loc_redo
        z_loc = z_loc_redo

        # Now, trim x, y, and z by the minimum x and y spacings. THe lowest x and y
        # values are taken as the starting positions. Then, it proceeds to the next
        # x or y value and sees if it's at least min_x_spacing_m or min_y_spacing_m
        # away from the last position. If it isn't, locations with that x or y value
        # are not included.
        min_x_spacing_m = 10000
        min_y_spacing_m = 10000
        unique_x_vals = np.unique(x_loc)
        unique_y_vals = np.unique(y_loc)
        chosen_x_vals = np.array([])
        chosen_y_vals = np.array([])
        for loc_ref in range(1, len(unique_x_vals)):
            if loc_ref == 1:
                last_point = unique_x_vals[0]

            if unique_x_vals[loc_ref] - last_point >= min_x_spacing_m:
                chosen_x_vals = np.append(chosen_x_vals, unique_x_vals[loc_ref])
                last_point = unique_x_vals[loc_ref]

        for loc_ref in range(1, len(unique_y_vals)):
            if loc_ref == 1:
                last_point = unique_y_vals[0]

            if unique_y_vals[loc_ref] - last_point >= min_y_spacing_m:
                chosen_y_vals = np.append(chosen_y_vals, unique_y_vals[loc_ref])
                last_point = unique_y_vals[loc_ref]

        # Go through the x values and add only the points with values that are
        # in the chosen_x_vals
        x_loc_redo = np.array([])
        y_loc_redo = np.array([])
        z_loc_redo = np.array([])

        for loc_ref in range(0, len(x_loc)):
            if x_loc[loc_ref] in chosen_x_vals:
                x_loc_redo = np.append(x_loc_redo, x_loc[loc_ref])
                y_loc_redo = np.append(y_loc_redo, y_loc[loc_ref])
                z_loc_redo = np.append(z_loc_redo, z_loc[loc_ref])

        x_loc = x_loc_redo
        y_loc = y_loc_redo
        z_loc = z_loc_redo

        # Go through the y values and add only the points with values that are
        # in the chosen_y_vals
        x_loc_redo = np.array([])
        y_loc_redo = np.array([])
        z_loc_redo = np.array([])

        for loc_ref in range(0, len(y_loc)):
            if y_loc[loc_ref] in chosen_y_vals:
                x_loc_redo = np.append(x_loc_redo, x_loc[loc_ref])
                y_loc_redo = np.append(y_loc_redo, y_loc[loc_ref])
                z_loc_redo = np.append(z_loc_redo, z_loc[loc_ref])

        x_loc = x_loc_redo
        y_loc = y_loc_redo
        z_loc = z_loc_redo

    elif use_location_range_option:
        x_vals = np.arange(xmin, xmax + dx, dx)
        y_vals = np.arange(ymin, ymax + dy, dy)

        x_loc = np.zeros(shape = (len(x_vals) * len(y_vals)))
        y_loc = np.zeros(shape = (len(x_vals) * len(y_vals)))

        row_ref = -1
        for x_ref in range(0, len(x_vals)):
            for y_ref in range(0, len(y_vals)):
                row_ref += 1
                x_loc[row_ref] = x_vals[x_ref]
                y_loc[row_ref] = y_vals[y_ref]



    print('Number of locations', len(x_loc))

    xmin = min(x_loc)
    xmax = max(x_loc)
    ymin = min(y_loc)
    ymax = max(y_loc)
    zmin = min(z_loc)
    zmax = max(z_loc)

    output_directory = os.sep.join([
        '..', '..', 'Output', Scenario_name + '_AoR_{date_time_stamp}'.format(
            date_time_stamp=datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))])
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Read file with signatures of interpolators and names of files with the corresponding data
    signature_data = np.genfromtxt(
        os.sep.join([file_directory, 'parameters_and_filenames_trunc.csv']),
        delimiter=",", dtype='str')

    # The first row (except the last element) of the file contains names of the parameters
    par_names = signature_data[0, 1:-1]

    num_pars = len(par_names)

    num_interpolators = signature_data.shape[0]-1  # -1 since the first line is a header



    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Create and add interpolator to the system model
    ind = run-1
    signature = {par_names[j]: float(signature_data[ind+1, j+1]) for j in range(num_pars)}

    sm.add_interpolator(ReservoirDataInterpolator(
        name='int{}'.format(ind+1), parent=sm,
        header_file_dir=file_directory,
        time_file='time_points.csv', interp_2d=False,
        data_file=Scenario_name + '.csv',
        index=int(signature_data[ind+1, 0]),
        signature=signature), intr_family='reservoir')

    ltress = []
    stratas = []
    ows = []
    adapts = []
    fgas = []

    for loc_ref in range(0, len(x_loc)):

        x = x_loc[loc_ref]
        y = y_loc[loc_ref]
        z = z_loc[loc_ref]

        # Adjust strata based on model depth
        shale[4] = (- z - shale[0] - shale[1] - shale[2] - shale[3] - aq_thick[0]
                    - aq_thick[1] - aq_thick[2] - aq_thick[3])
        depth = 4*[None]  # list with four elements to be defined below
        depth[3] = shale[4] + aq_thick[3]
        depth[2] = depth[3] + shale[3] + aq_thick[2]
        depth[1] = depth[2] + shale[2] + aq_thick[1]
        depth[0] = depth[1] + shale[1] + aq_thick[0]

        # Add stratigraphy component
        stratas.append(sm.add_component_model_object(Stratigraphy(name='strata', parent=sm)))

        # Add parameters of stratigraphy component model
        stratas[-1].add_par('numberOfShaleLayers', value=5, vary=False)
        stratas[-1].add_par('shale1Thickness', value=shale[0], vary=False)
        stratas[-1].add_par('shale2Thickness', value=shale[1], vary=False)
        stratas[-1].add_par('shale3Thickness', value=shale[2], vary=False)
        stratas[-1].add_par('shale4Thickness', value=shale[3], vary=False)
        stratas[-1].add_par('shale5Thickness', value=shale[4], vary=False)
        stratas[-1].add_par('aquifer1Thickness', value=aq_thick[0], vary=False)
        stratas[-1].add_par('aquifer2Thickness', value=aq_thick[1], vary=False)
        stratas[-1].add_par('aquifer3Thickness', value=aq_thick[2], vary=False)
        stratas[-1].add_par('aquifer4Thickness', value=aq_thick[3], vary=False)
        stratas[-1].add_par('reservoirThickness', value=7.0, vary=False)

        # Add reservoir component
        ltress.append(sm.add_component_model_object(LookupTableReservoir(
            name='ltres' + str(loc_ref), parent=sm, intr_family='reservoir',
            locX=x_loc[loc_ref], locY=y_loc[loc_ref], locZ=z, interp_2d=False)))

        # Add parameters of reservoir component model
        for j in range(num_pars):
            # add arbitrary line of values from signature_file
            ltress[-1].add_par(
                par_names[j], value=float(signature_data[run, j+1]), vary=False)

        # Add observations of reservoir component model
        ltress[-1].add_obs('pressure')
        ltress[-1].add_obs('CO2saturation')
        ltress[-1].add_obs_to_be_linked('pressure')
        ltress[-1].add_obs_to_be_linked('CO2saturation')

        # Add open wellbore component
        ows.append(sm.add_component_model_object(
            OpenWellbore(name='ow' + str(loc_ref), parent=sm)))

        # Add parameters of open wellbore component
        ows[-1].add_par('wellRadius', value=0.05, vary=False)
        ows[-1].add_par('logReservoirTransmissivity', value=-10.0, vary=False)
        ows[-1].add_par('logAquiferTransmissivity', value=-10.0, vary=False)
        ows[-1].add_par('brineSalinity', value=0.0475, vary=False)

        # Add keyword arguments of the open wellbore component model
        ows[-1].add_kwarg_linked_to_obs('pressure', ltress[-1].linkobs['pressure'])
        ows[-1].add_kwarg_linked_to_obs('CO2saturation', ltress[-1].linkobs['CO2saturation'])

        res_depth = (shale[0] + shale[1] + shale[2] + shale[3] + shale[4]
                     + aq_thick[0] + aq_thick[1] + aq_thick[2] + aq_thick[3])
        # Add composite parameter of open wellbore component
        ows[-1].add_par('reservoirDepth', value=res_depth, vary=False)

        wellTop_val = shale[4] + aq_thick[3]
        ows[-1].add_par('wellTop', value=wellTop_val, vary=False)

        # Add observations of open wellbore component model
        ows[-1].add_obs_to_be_linked('CO2_aquifer')
        ows[-1].add_obs_to_be_linked('brine_aquifer')
        ows[-1].add_obs_to_be_linked('brine_atm')
        ows[-1].add_obs_to_be_linked('CO2_atm')
        ows[-1].add_obs('CO2_aquifer')
        ows[-1].add_obs('brine_aquifer')
        ows[-1].add_obs('CO2_atm') # zero since well top is in aquifer
        ows[-1].add_obs('brine_atm') # zero since well top is in aquifer

        # Add adapter that transforms leakage rates to accumulated mass
        adapts.append(sm.add_component_model_object(
            RateToMassAdapter(name='adapt' + str(loc_ref), parent=sm)))
        adapts[-1].add_kwarg_linked_to_collection(
            'CO2_aquifer', [ows[-1].linkobs['CO2_aquifer'],
                            ows[-1].linkobs['CO2_atm']])
        adapts[-1].add_kwarg_linked_to_collection('brine_aquifer',
                                             [ows[-1].linkobs['brine_aquifer'],
                                              ows[-1].linkobs['brine_atm']])
        adapts[-1].add_obs_to_be_linked('mass_CO2_aquifer')
        adapts[-1].add_obs_to_be_linked('mass_brine_aquifer')
        adapts[-1].add_obs('mass_CO2_aquifer')
        adapts[-1].add_obs('mass_brine_aquifer')

        # Add futuregen aquifer model object and define parameters
        if float(depth[3]) < 700:
            fgas.append(sm.add_component_model_object(
                FutureGen2Aquifer(name='fga' + str(loc_ref), parent=sm)))
        else:
            fgas.append(sm.add_component_model_object(
                FutureGen2AZMI(name='fga' + str(loc_ref), parent=sm)))

        # St. Peter Sandstone
        fgas[-1].add_par('aqu_thick', value=aq_thick[3], vary=False)
        fgas[-1].add_par('depth', value=depth[3], vary=False)
        fgas[-1].add_par('por', value=0.18, vary=False)
        fgas[-1].add_par('log_permh', value=-11.92, vary=False)
        fgas[-1].add_par('log_aniso', value=0.3, vary=False)
        fgas[-1].add_par('rel_vol_frac_calcite', value=0.01, vary=False)

        # Add aquifer component's keyword argument co2_rate linked to the collection created above
        fgas[-1].add_kwarg_linked_to_obs('co2_rate',
                                         ows[-1].linkobs['CO2_aquifer'])
        fgas[-1].add_kwarg_linked_to_obs('brine_rate',
                                         ows[-1].linkobs['brine_aquifer'])
        fgas[-1].add_kwarg_linked_to_obs('co2_mass',
                                         adapts[-1].linkobs['mass_CO2_aquifer'])
        fgas[-1].add_kwarg_linked_to_obs('brine_mass',
                                         adapts[-1].linkobs['mass_brine_aquifer'])

        # Add observations (output) from the futuregen aquifer model
        fgas[-1].add_obs('pH_volume')
        fgas[-1].add_obs('TDS_volume')

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically

    results = np.zeros((len(x_loc), 5))

    for loc_ref in range(0, len(x_loc)):
        # Get results
        pressure = sm.collect_observations_as_time_series(
            ltress[loc_ref], 'pressure') / 1.0e6
        results[loc_ref, 0] = max(pressure)
        results[loc_ref, 1] = max([p - pressure[0] for p in pressure])
        results[loc_ref, 2] = max(
            sm.collect_observations_as_time_series(ltress[loc_ref], 'CO2saturation'))
        results[loc_ref, 3] = max(
            sm.collect_observations_as_time_series(fgas[loc_ref], 'pH_volume'))
        results[loc_ref, 4] = max(
            sm.collect_observations_as_time_series(fgas[loc_ref], 'TDS_volume'))


    results_formatted = np.empty(((len(x_loc) + 1), 7), dtype=list)
    results_formatted[0, 0] = 'x (km)'
    results_formatted[0, 1] = 'y (km)'
    for row_ref in range(0, len(x_loc)):
        results_formatted[row_ref + 1, 0] = str('%.2e' % (x_loc[row_ref] / 1000))
        results_formatted[row_ref + 1, 1] = str('%.2e' % (y_loc[row_ref] / 1000))
    results_formatted[0, 2] = 'Max Pressure (MPa)'
    results_formatted[0, 3] = 'Max Pressure Inc. (MPa)'
    results_formatted[0, 4] = 'Max CO2 Saturation'
    results_formatted[0, 5] = 'Max pH Plume Volume (m^3)'
    results_formatted[0, 6] = 'Max TDS Plume Volume (m^3)'
    results_formatted[1:None, 2:None] = results

    os.chdir(output_directory)

    # Save the ouput for the simulation
    with open('AOR_output.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        for row_ref in range(0, len(x_loc) + 1):
            writer.writerow(results_formatted[row_ref, :])

    # Figures
    font = {'family': 'Arial',
            'weight': 'normal',
            'size': 12}
    plt.rc('font', **font)

    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    # Pressure Figure
    column = 1
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    ax.set_facecolor([0.67, 0.67, 0.67])
    plt.plot(x_loc / 1000, y_loc / 1000, linestyle='none', marker='o',
             color='k', markeredgewidth=1, markersize=6,
             markerfacecolor='none')
    if np.max(results[:, column]) != 0:
        # I use results[:, 1] * 1.01 here because having the maximum level
        # equal to the maximum value will cause the point with the maximum
        # value to be left out (i.e., there would be an uncolored area there).
        levels = np.arange(np.min(results[:, column][results[:, column] > 0]),
                           np.max(results[:, column]) * 1.01,
                           ((np.max(results[:, column]) * 1.01)
                            - np.min(results[:, column][results[:, column] > 0])) / 100)
        plt.tricontourf(x_loc / 1000.0, y_loc / 1000.0, results[:, column], levels,
                        locator=ticker.MaxNLocator(nbins=100, prune='lower'),
                        cmap="viridis")
        cbar = plt.colorbar()
        cbar.set_label('Increase in Pressure (MPa)', rotation=90, fontsize=14,
                       fontweight='bold')

        # Plot colors for individual points so there is less ambiguity
        cmap = plt.cm.get_cmap("viridis")
        for loc_ref in range(0, len(results[:, column])):
            if results[loc_ref, column] != 0:
                # Get the color by the value's proximity to the upper and lower
                # limits used for levels (which set the colorbar limits)
                rgba = cmap((results[loc_ref, column] - np.min(levels))
                             / (np.max(levels) - np.min(levels)))
                plt.plot(x_loc[loc_ref] / 1000, y_loc[loc_ref] / 1000, marker='o',
                         markerfacecolor=rgba[0:3], markeredgecolor='k',
                         markeredgewidth=1.5, markersize=9)
    plt.xlim((np.min(x_loc) - dx) / 1000.0,
             (np.max(x_loc) + dx) / 1000.0)
    plt.ylim((np.min(y_loc) - dy) / 1000.0,
             (np.max(y_loc) + dy) / 1000.0)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.xlabel('Easting (km)', fontsize=14, fontweight='bold')
    plt.ylabel('Northing (km)', fontsize=14, fontweight='bold')
    plt.title('Maximum Increase in Pressure at each point (gray: 0 MPa)',
              fontsize=14, fontweight='bold')
    plt.tight_layout()
    to_save = True
    if to_save:
        plt.savefig(os.sep.join([output_directory,
                                 "Pressure_Inc{}.png".format(args.run)]), dpi=300)


    # CO2 Saturation Figure
    column = 2
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    ax.set_facecolor([0.67, 0.67, 0.67])
    plt.plot(x_loc / 1000, y_loc / 1000, linestyle='none', marker='o',
             color='k', markeredgewidth=1, markersize=6,
             markerfacecolor='none')
    if np.max(results[:, column]) != 0:
        levels = np.arange(np.min(results[:, column][results[:, column] > 0]),
                           np.max(results[:, column]) * 1.01,
                           ((np.max(results[:, column]) * 1.01)
                            - np.min(results[:, column][results[:, column] > 0])) / 100)
        plt.tricontourf(x_loc / 1000.0, y_loc / 1000.0, results[:, column], levels,
                        locator=ticker.MaxNLocator(nbins=100, prune='lower'),
                        cmap="plasma")
        cbar = plt.colorbar()
        cbar.set_label('CO$_2$ Saturation [-]', rotation=90,
                       fontsize=14, fontweight='bold')

        # Plot colors for individual points so there is less ambiguity
        cmap = plt.cm.get_cmap("plasma")
        for loc_ref in range(0, len(results[:, column])):
            if results[loc_ref, column] != 0:
                # Get the color by the value's proximity to the upper and lower
                # limits used for levels (which set the colorbar limits)
                rgba = cmap((results[loc_ref, column] - np.min(levels))
                             / (np.max(levels) - np.min(levels)))
                plt.plot(x_loc[loc_ref] / 1000, y_loc[loc_ref] / 1000, marker='o',
                         markerfacecolor=rgba[0:3], markeredgecolor='k',
                         markeredgewidth=1.5, markersize=9)
    plt.xlim((np.min(x_loc) - dx) / 1000.0,
             (np.max(x_loc) + dx) / 1000.0)
    plt.ylim((np.min(y_loc) - dy) / 1000.0,
             (np.max(y_loc) + dy) / 1000.0)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.xlabel('Easting (km)', fontsize=14, fontweight='bold')
    plt.ylabel('Northing (km)', fontsize=14, fontweight='bold')
    plt.title('Maximum CO$_2$ Saturation at each point (gray: 0)',
              fontsize=14, fontweight='bold')
    plt.tight_layout()
    to_save = True
    if to_save:
        plt.savefig(os.sep.join([
            output_directory, "CO2_Sat{}.png".format(args.run)]), dpi=300)


    # pH Volume Figure
    column = 3
    fig, ax = plt.subplots(figsize=(13,8), dpi=300)
    ax.set_facecolor([0.67, 0.67, 0.67])
    plt.plot(x_loc / 1000, y_loc / 1000, linestyle='none', marker='o',
             color='k', markeredgewidth=1, markersize=6,
             markerfacecolor='none')
    if np.max(results[:, column]) != 0:
        levels = np.arange(np.min(results[:, column][results[:, column] > 0]),
                           np.max(results[:, column]) * 1.01,
                           ((np.max(results[:, column]) * 1.01)
                            - np.min(results[:, column][results[:, column] > 0])) / 100)
        plt.tricontourf(x_loc / 1000.0, y_loc / 1000.0, results[:, column], levels,
                        locator=ticker.MaxNLocator(nbins=100, prune='lower'), cmap="GnBu")
        cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
        cbar.set_label('pH plume volume (m$^3$)', rotation=90, fontsize=14,
                       fontweight='bold')

        # Plot colors for individual points so there is less ambiguity
        cmap = plt.cm.get_cmap("GnBu")
        for loc_ref in range(0, len(results[:, column])):
            if results[loc_ref, column] != 0:
                # Get the color by the value's proximity to the upper and lower
                # limits used for levels (which set the colorbar limits)
                rgba = cmap((results[loc_ref, column] - np.min(levels))
                             / (np.max(levels) - np.min(levels)))
                plt.plot(x_loc[loc_ref] / 1000, y_loc[loc_ref] / 1000, marker='o',
                         markerfacecolor=rgba[0:3], markeredgecolor='k',
                         markeredgewidth=1.5, markersize=9)
    plt.xlim((np.min(x_loc) - dx) / 1000.0,
             (np.max(x_loc) + dx) / 1000.0)
    plt.ylim((np.min(y_loc) - dy) / 1000.0,
             (np.max(y_loc) + dy) / 1000.0)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.xlabel('Easting (km)', fontsize=14, fontweight='bold')
    plt.ylabel('Northing (km)', fontsize=14, fontweight='bold')
    plt.title('Maximum pH plume volume at each point (gray: 0 m$^3$)',
              fontsize=14, fontweight='bold')
    plt.tight_layout()
    to_save = True
    if to_save:
        plt.savefig(os.sep.join([
            output_directory, "AoR_pH{}.png".format(args.run)]), dpi=300)


    # TDS Volume Figure
    column = 4
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    ax.set_facecolor([0.67, 0.67, 0.67])
    plt.plot(x_loc / 1000, y_loc / 1000, linestyle='none', marker='o',
             color='k', markeredgewidth=1, markersize=6,
             markerfacecolor='none')
    if np.max(results[:, 4]) != 0:
        levels=np.arange(np.min(results[:, column][results[:, column] > 0]),
                           np.max(results[:, column]) * 1.01,
                           ((np.max(results[:, column]) * 1.01)
                            - np.min(results[:, column][results[:, column] > 0])) / 100)
        plt.tricontourf(x_loc / 1000.0, y_loc / 1000.0, results[:, column], levels,
                        locator=ticker.MaxNLocator(nbins=100, prune='lower'),
                        cmap="YlOrRd")
        cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
        cbar.set_label('TDS plume volume (m$^3$)', rotation=90,
                       fontsize=14, fontweight='bold')

        # Plot colors for individual points so there is less ambiguity
        cmap = plt.cm.get_cmap("YlOrRd")
        for loc_ref in range(0, len(results[:, column])):
            if results[loc_ref, column] != 0:
                # Get the color by the value's proximity to the upper and lower
                # limits used for levels (which set the colorbar limits)
                rgba = cmap((results[loc_ref, column] - np.min(levels))
                             / (np.max(levels) - np.min(levels)))
                plt.plot(x_loc[loc_ref] / 1000, y_loc[loc_ref] / 1000, marker='o',
                         markerfacecolor=rgba[0:3], markeredgecolor='k',
                         markeredgewidth=1.5, markersize=9)
    plt.xlim((np.min(x_loc) - dx) / 1000.0,
             (np.max(x_loc) + dx) / 1000.0)
    plt.ylim((np.min(y_loc) - dy) / 1000.0,
             (np.max(y_loc) + dy) / 1000.0)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.xlabel('Easting (km)', fontsize=14, fontweight='bold')
    plt.ylabel('Northing (km)', fontsize=14, fontweight='bold')
    plt.title('Maximum TDS plume volume at each point (gray: 0 m$^3$)',
              fontsize=14, fontweight='bold')
    plt.tight_layout()
    to_save = True
    if to_save:
        plt.savefig(os.sep.join([
            output_directory, "AoR_TDS{}.png".format(args.run)]), dpi=300)
