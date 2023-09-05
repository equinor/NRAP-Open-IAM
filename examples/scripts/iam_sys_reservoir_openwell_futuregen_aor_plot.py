'''
NRAP IAM AOR FutureGen Area of Review (AoR)
This example couples the simple reservoir, open wellbore, and FutureGen2 Aquifer
components. The saturation/pressure output produced by lookup table
reservoir model is used to drive leakage from a single open wellbore
model, which is passed to the input of an adapter that provides well
coordinates, CO2 and brine leakage rates and cumulative mass fluxes to the
FutureGen2 aquifer model.  The strata component defines the thickness of the aquifer
and shale layers above the reservoir, which are then used as input parameters
for the open wellbore and aquifer components. The risk-based Area of Review (AoR)
is defined based on the probability of whether there is an impact to groundwater
at each grid location.

This examples was created to serve as a script-based example of the analyses and
figure creation performed by control file example ControlFile_ex31.yaml.
That control file uses the area_of_review_plot() function within the area_of_review.py file.

This example was based on the file iam_sys_lutreservoir_openwell_futuregen_aor.py
created by Diana Bacon (diana.bacon@pnnl.gov).

Usage examples:
$ python iam_sys_reservoir_openwell_futuregen_aor_plot.py
'''
# @author: Nate Mitchell
# Nathaniel.Mitchell@netl.doe.gov

import os
import sys
import datetime
import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, SimpleReservoir, OpenWellbore,
                     FutureGen2Aquifer, RateToMassAdapter, Stratigraphy)

if __name__ == "__main__":
    # For multiprocessing in Spyder
    __spec__ = None

    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # These values are used to create the x and y locations.
    xmin = 5000
    xmax = 145000
    dx = 10000

    ymin = 5000
    ymax = 145000
    dy = 10000

    x_vals = np.arange(xmin, xmax + dx, dx)
    y_vals = np.arange(ymin, ymax + dy, dy)

    x_loc = np.zeros(shape=(len(x_vals) * len(y_vals)))
    y_loc = np.zeros(shape=(len(x_vals) * len(y_vals)))

    row_ref = -1
    for x_ref in range(0, len(x_vals)):
        for y_ref in range(0, len(y_vals)):
            row_ref += 1
            x_loc[row_ref] = x_vals[x_ref]
            y_loc[row_ref] = y_vals[y_ref]

    print('Number of locations: ', len(x_loc))

    output_directory = os.sep.join([
        '..', '..', 'Output', 'sres_ow_fgen_AoR_{date_time_stamp}'.format(
            date_time_stamp=datetime.datetime.now().strftime('%Y-%m-%d_%H'))])
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add stratigraphy component
    strata = sm.add_component_model_object(Stratigraphy(name='strata', parent=sm))

    # Add parameters of stratigraphy component model
    strata.add_par('numberOfShaleLayers', value=3, vary=False)
    strata.add_par('shale1Thickness', value=750, vary=False)
    strata.add_par('shale2Thickness', value=950, vary=False)
    strata.add_par('shale3Thickness', value=365, vary=False)
    strata.add_par('aquifer1Thickness', value=85, vary=False)
    strata.add_par('aquifer2Thickness', value=85, vary=False)
    strata.add_par('reservoirThickness', value=150, vary=False)
    strata.add_par('datumPressure', value=101325, vary=False)

    sress = []
    ows = []
    adapts = []
    fgas = []

    for loc_ref in range(0, len(x_loc)):
        # Add reservoir component
        sress.append(sm.add_component_model_object(SimpleReservoir(
            name='sres' + str(loc_ref), parent=sm,
            locX=x_loc[loc_ref], locY=y_loc[loc_ref])))

        # Add parameters of simple reservoir component
        sress[-1].add_par('injRate', value=1, vary=False)
        sress[-1].add_par('logResPerm', value=-10.15, vary=False)
        sress[-1].add_par('reservoirPorosity', value=0.3, vary=False)

        sress[-1].add_par_linked_to_par(
            'reservoirThickness', strata.deterministic_pars['reservoirThickness'])
        sress[-1].add_par_linked_to_par(
            'numberOfShaleLayers', strata.deterministic_pars['numberOfShaleLayers'])
        sress[-1].add_par_linked_to_par(
            'shale1Thickness', strata.deterministic_pars['shale1Thickness'])
        sress[-1].add_par_linked_to_par(
            'shale2Thickness', strata.deterministic_pars['shale2Thickness'])
        sress[-1].add_par_linked_to_par(
            'shale3Thickness', strata.deterministic_pars['shale3Thickness'])
        sress[-1].add_par_linked_to_par(
            'aquifer1Thickness', strata.deterministic_pars['aquifer1Thickness'])
        sress[-1].add_par_linked_to_par(
            'aquifer2Thickness', strata.deterministic_pars['aquifer2Thickness'])
        sress[-1].add_par_linked_to_par(
            'datumPressure', strata.deterministic_pars['datumPressure'])

        # Add observations of reservoir component model
        sress[-1].add_obs('pressure')
        sress[-1].add_obs('CO2saturation')
        sress[-1].add_obs_to_be_linked('pressure')
        sress[-1].add_obs_to_be_linked('CO2saturation')

        # Add open wellbore component
        ows.append(sm.add_component_model_object(
            OpenWellbore(name='ow' + str(loc_ref), parent=sm)))

        # Add parameters of open wellbore component
        ows[-1].add_par('wellRadius', value=0.05, vary=False)
        ows[-1].add_par('logReservoirTransmissivity', value=-10.97, vary=False)
        ows[-1].add_par('logAquiferTransmissivity', value=-10.0, vary=False)
        ows[-1].add_par('brineSalinity', value=0.2, vary=False)

        # Add keyword arguments of the open wellbore component model
        ows[-1].add_kwarg_linked_to_obs('pressure', sress[-1].linkobs['pressure'])
        ows[-1].add_kwarg_linked_to_obs('CO2saturation', sress[-1].linkobs['CO2saturation'])

        # Add composite parameter of open wellbore component
        ows[-1].add_composite_par(
            'reservoirDepth',
            expr='+'.join(['strata.shale1Thickness', 'strata.shale2Thickness',
                            'strata.shale3Thickness',
                            'strata.aquifer1Thickness', 'strata.aquifer2Thickness']))
        ows[-1].add_composite_par(
            'wellTop', expr='strata.shale3Thickness + strata.aquifer2Thickness')

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
        adapts.append(sm.add_component_model_object(RateToMassAdapter(
            name='adapt' + str(loc_ref), parent=sm)))
        adapts[-1].add_kwarg_linked_to_collection(
            'CO2_aquifer', [ows[-1].linkobs['CO2_aquifer'],
                            ows[-1].linkobs['CO2_atm']])
        adapts[-1].add_kwarg_linked_to_collection(
            'brine_aquifer', [ows[-1].linkobs['brine_aquifer'],
                              ows[-1].linkobs['brine_atm']])
        adapts[-1].add_obs_to_be_linked('mass_CO2_aquifer')
        adapts[-1].add_obs_to_be_linked('mass_brine_aquifer')
        adapts[-1].add_obs('mass_CO2_aquifer')
        adapts[-1].add_obs('mass_brine_aquifer')

        # Add futuregen aquifer model object and define parameters
        fgas.append(sm.add_component_model_object(FutureGen2Aquifer(
            name='fga' + str(loc_ref), parent=sm)))

        # St. Peter Sandstone
        fgas[-1].add_par('aqu_thick', value=85, vary=False)
        fgas[-1].add_par('depth', value=450, vary=False)
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
            sress[loc_ref], 'pressure') / 1.0e6

        results[loc_ref, 0] = max(pressure)
        results[loc_ref, 1] = max([p - pressure[0] for p in pressure])
        results[loc_ref, 2] = max(
            sm.collect_observations_as_time_series(sress[loc_ref], 'CO2saturation'))
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

    # This function is used to format the labels on colorbars
    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    # Pressure Figure
    column = 1
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    ax.set_facecolor([0.67, 0.67, 0.67])
    plt.plot(x_loc / 1000, y_loc / 1000, linestyle='none', marker='o',
             color='k', markeredgewidth=1, markersize=6, markerfacecolor='none')
    if np.max(results[:, column]) != 0:
        # I use results[:, column] * 1.01 here because having the maximum level
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
        cbar.set_label('Increase in Pressure (MPa)', rotation=90,
                       fontsize=14, fontweight='bold')

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
        plt.savefig(os.sep.join([output_directory, "Pressure_Inc.png"]), dpi=300)

    # CO2 Saturation Figure
    column = 2
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    ax.set_facecolor([0.67, 0.67, 0.67])
    plt.plot(x_loc / 1000, y_loc / 1000, linestyle='none', marker='o',
             color='k', markeredgewidth=1, markersize=6, markerfacecolor='none')
    if np.max(results[:, column]) != 0:
        levels = np.arange(np.min(results[:, column][results[:, column] > 0]),
                           np.max(results[:, column]) * 1.01,
                           ((np.max(results[:, column]) * 1.01)
                            - np.min(results[:, column][results[:, column] > 0])) / 100)
        plt.tricontourf(x_loc / 1000.0, y_loc / 1000.0, results[:, column], levels,
                        locator=ticker.MaxNLocator(nbins=100, prune='lower'), cmap="plasma")
        cbar = plt.colorbar()
        cbar.set_label('CO$_2$ Saturation [-]', rotation=90, fontsize=14, fontweight='bold')

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
        plt.savefig(os.sep.join([output_directory, "CO2_Sat.png"]), dpi=300)

    # pH Volume Figure
    column = 3
    fig, ax = plt.subplots(figsize=(13,8), dpi=300)
    ax.set_facecolor([0.67, 0.67, 0.67])
    plt.plot(x_loc / 1000, y_loc / 1000, linestyle='none', marker='o',
             color='k', markeredgewidth=1, markersize=6, markerfacecolor='none')
    if np.max(results[:, column]) != 0:
        levels = np.arange(np.min(results[:, column][results[:, column] > 0]),
                           np.max(results[:, column]) * 1.01,
                           ((np.max(results[:, column]) * 1.01)
                            - np.min(results[:, column][results[:, column] > 0])) / 100)
        plt.tricontourf(x_loc / 1000.0, y_loc / 1000.0, results[:, column], levels,
                        locator=ticker.MaxNLocator(nbins=100, prune='lower'), cmap="GnBu")
        cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
        cbar.set_label('pH plume volume (m$^3$)', rotation=90, fontsize=14, fontweight='bold')

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
        plt.savefig(os.sep.join([output_directory, "AoR_pH.png"]), dpi=300)

    # TDS Volume Figure
    column = 4
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    ax.set_facecolor([0.67, 0.67, 0.67])
    plt.plot(x_loc / 1000, y_loc / 1000, linestyle='none', marker='o',
             color='k', markeredgewidth=1, markersize=6, markerfacecolor='none')
    if np.max(results[:, column]) != 0:
        levels = np.arange(np.min(results[:, column][results[:, column] > 0]),
                           np.max(results[:, column]) * 1.01,
                           ((np.max(results[:, column]) * 1.01)
                            - np.min(results[:, column][results[:, column] > 0])) / 100)
        plt.tricontourf(x_loc / 1000.0, y_loc / 1000.0, results[:, column], levels,
                        locator=ticker.MaxNLocator(nbins=100, prune='lower'), cmap="YlOrRd")
        cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
        cbar.set_label('TDS plume volume (m$^3$)', rotation=90, fontsize=14, fontweight='bold')

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
        plt.savefig(os.sep.join([output_directory, "AoR_TDS.png"]), dpi=300)
