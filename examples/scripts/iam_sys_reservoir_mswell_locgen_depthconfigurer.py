"""
This example demonstrates the use of the LocationGenerator and WellDepthRiskConfigurer
components. The x, y, and z values (easting, northing, and depth) for wellbores
are randomly created by the LocationGenerator. The well depths are then provided
to the WellDepthRiskConfigurer - if the well does not reach the top of the
reservoir, then the WellDepthRiskConfigurer disables the well by setting the
run_frequency to 0 (so the wellbore component never runs).

This simulation uses a forward analysis.

Note that randomly assigning coordinates with the LocationGenerator is different
from first creating random x, y, and z values and then directly providing those
to a reservoir component. The difference lies in the fact that the LocationGenerator
creates the x, y, and z values during the actual simulation - if you run a
simulation with Latin Hypercube Sampling, each realization will have different
x, y, and z values because they are generated during the realization. Having
random well locations generated in each realization can be important if you know
that well records are incomplete - this approach allows one to assess the potential
significance of unknown wellbores. Some of the unknown wells may not penetrate
the storage reservoir and may therefore be excluded as potential leakage pathways.
"""
# Author: Nate Mitchell
# Nathaniel.Mitchell@netl.doe.gov
# Contributor: Veronika Vasylkivska
# Veronika.Vasylkivska@netl.doe.gov

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import ticker
from matplotlib.lines import Line2D

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, Stratigraphy, SimpleReservoir,
                     MultisegmentedWellbore, LocationGenerator,
                     WellDepthRiskConfigurer)


def get_distance_from_site_km(sm, comp_ref, injection_x_m, injection_y_m):
    """
    Returns the current component's distance from the injection site.
    """
    x_val_m = sm.obs['gen.locX{}_0'.format(comp_ref)].sim
    y_val_m = sm.obs['gen.locY{}_0'.format(comp_ref)].sim

    distance_from_site_km = (
        (((injection_x_m - x_val_m) ** 2) + ((injection_y_m - y_val_m) ** 2))
        ** 0.5) / 1000

    return distance_from_site_km


if __name__ == "__main__":

    # Simulation input
    num_wells = 50

    seed_val = 11

    num_years = 50

    shale1Thickness = 250
    shale2Thickness = 250
    shale3Thickness = 250
    aquifer1Thickness = 100
    aquifer2Thickness = 100
    reservoirThickness = 50

    # Define boundaries of random wells domains
    x_min = -5000.0
    x_max = 5000.0

    y_min = -5000.0
    y_max = 5000.0

    injection_x_m = 0
    injection_y_m = 0

    # Boundaries defining well depths
    # Here, the z values can range from the bottom of aquifer 1 (top of shale 2) to
    # the bottom of the reservoir.

    # Bottom of aquifer 1
    z_min = shale3Thickness + aquifer2Thickness + shale2Thickness + aquifer1Thickness

    # Bottom of the reservoir
    z_max = (shale3Thickness + aquifer2Thickness + shale2Thickness
              + aquifer1Thickness + shale1Thickness + reservoirThickness)

    # The reservoirDepth is the depth to the top of the reservoir
    reservoirDepth = (shale3Thickness + aquifer2Thickness + shale2Thickness
                      + aquifer1Thickness + shale1Thickness)

    # Define output directory for gridded observations
    output_dir = os.path.join('..', '..', 'output', 'script_well_depth_config')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Define keyword arguments of the system model
    time_array = 365.25 * np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add stratigraphy component
    strata = sm.add_component_model_object(Stratigraphy(name='strata', parent=sm))

    # Add parameters of stratigraphy component model
    strata.add_par('numberOfShaleLayers', value=3, vary=False)
    strata.add_par('shale1Thickness', value=shale1Thickness, vary=False)
    strata.add_par('shale2Thickness', value=shale2Thickness, vary=False)
    strata.add_par('shale3Thickness', value=shale3Thickness, vary=False)
    strata.add_par('aquifer1Thickness', value=aquifer1Thickness, vary=False)
    strata.add_par('aquifer2Thickness', value=aquifer2Thickness, vary=False)
    strata.add_par('reservoirThickness', value=reservoirThickness, vary=False)

    # Add generator component
    gen = sm.add_component_model_object(LocationGenerator(
        name='gen', parent=sm, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max,
        num_locations=num_wells, reproducible=True, z_min=z_min, z_max=z_max))

    # Add parameters of generator component
    gen.add_par('seed', value=seed_val, vary=False)

    for obs_nm in ['locX', 'locY', 'locZ']:
        gen.add_obs_to_be_linked(obs_nm, obs_type='grid')
        gen.add_grid_obs(obs_nm, 'array', output_dir=output_dir)

    # Locations can be obtained as observations of the generator component
    # Since the observations considered to be of gridded type the outputs
    # can only be obtained by adding them using add_obs method and
    # then reading the files with data after simulation is complete

    well_names = ['msw{}'.format(comp_ref) for comp_ref in range(num_wells)]

    # Add configurer component
    config = sm.add_component_model_object(WellDepthRiskConfigurer(
        name='cnfg', parent=sm, cmpnt_nms=well_names))

    config.add_par('reservoirDepth', value=reservoirDepth, vary=False)

    config.add_kwarg_linked_to_obs(
        'wellDepth', gen.linkobs['locZ'], obs_type='grid', constr_type='array')
    config.add_obs('num_cmpnts_on', index=[0])

    sres = []
    msw = []

    for comp_ref in range(num_wells):
        for obs_nm in ['locX', 'locY', 'locZ']:
            gen.add_local_obs(obs_nm+'{}'.format(comp_ref), grid_obs_name=obs_nm,
                              constr_type='array', loc_ind=comp_ref, index=[0])

        # Add simple reservoir component
        # We don't specify the location for the reservoir component, although
        # there is a default reservoir setup with locX=100, locY=100
        sres.append(sm.add_component_model_object(SimpleReservoir(
            name='sres{}'.format(comp_ref), parent=sm,
            injX=injection_x_m, injY=injection_y_m)))

        # Add parameters of reservoir component model
        # All parameters of the reservoir component are deterministic so all
        # uncertainty in the simulation comes from the uncertainty of the well location
        sres[-1].add_par('numberOfShaleLayers', value=3, vary=False)
        sres[-1].add_par('injRate', value=0.5, vary=False)
        sres[-1].add_par('shale1Thickness', value=shale1Thickness, vary=False)
        sres[-1].add_par('shale2Thickness', value=shale2Thickness, vary=False)
        sres[-1].add_par('shale3Thickness', value=shale3Thickness, vary=False)
        sres[-1].add_par('aquifer1Thickness', value=aquifer1Thickness, vary=False)
        sres[-1].add_par('aquifer2Thickness', value=aquifer2Thickness, vary=False)
        sres[-1].add_par('reservoirThickness', value=reservoirThickness, vary=False)

        # Simple reservoir component has keyword arguments of the model method:
        # locX and locY which we would link to the output of the generator component
        # Generator component outputs 5 random locations according to the setup above.
        # We can link reservoir component locX and locY to any of these produced locations
        # using arguments constr_type='array' and loc_ind=[comp_ref].
        for obs_nm in ['locX', 'locY']:
            sres[-1].add_kwarg_linked_to_obs(obs_nm, gen.linkobs[obs_nm],
                                             obs_type='grid', constr_type='array',
                                             loc_ind=[comp_ref])

        # Add observations of reservoir component model
        for obs_nm in ['pressure', 'CO2saturation']:
            sres[-1].add_obs(obs_nm)
            sres[-1].add_obs_to_be_linked(obs_nm)

        # Add multisegmented wellbore component
        msw.append(sm.add_component_model_object(MultisegmentedWellbore(
            name=well_names[comp_ref], parent=sm)))

        # Add parameters of multisegmented wellbore component
        msw[-1].add_par('wellRadius', value=0.015, vary=False)
        msw[-1].add_par('numberOfShaleLayers', value=3, vary=False)
        msw[-1].add_par('shale1Thickness', value=shale1Thickness, vary=False)
        msw[-1].add_par('shale2Thickness', value=shale2Thickness, vary=False)
        msw[-1].add_par('shale3Thickness', value=shale3Thickness, vary=False)
        msw[-1].add_par('aquifer1Thickness', value=aquifer1Thickness, vary=False)
        msw[-1].add_par('aquifer2Thickness', value=aquifer2Thickness, vary=False)
        msw[-1].add_par('reservoirThickness', value=reservoirThickness, vary=False)
        msw[-1].add_par('logWellPerm', value=-13.0, vary=False)
        msw[-1].add_par('logAquPerm', value=-11.0, vary=False)

        # Add keyword arguments linked to the output provided by reservoir model
        for obs_nm in ['pressure', 'CO2saturation']:
            msw[-1].add_kwarg_linked_to_obs(obs_nm, sres[-1].linkobs[obs_nm])

        # Add observations of multisegmented wellbore component model
        for obs_nm in ['CO2_aquifer1', 'CO2_aquifer2', 'brine_aquifer1', 'brine_aquifer2']:
            msw[-1].add_obs(obs_nm)

    # Run forward simulation
    sm.forward()

    # Print the observations
    num_turn_on_wells = sm.obs['cnfg.num_cmpnts_on_0'].sim
    print('Number of components on:', num_turn_on_wells)

    # Load produced depth data
    locZ = sm.collect_gridded_observations_as_time_series(
            gen, 'locZ', output_dir, indices=[0], rlzn_number=0)[0]

    genfontsize = 8
    font = {'family': 'Arial',
            'weight': 'normal',
            'size': genfontsize}
    plt.rc('font', **font)
    dpi_ref = 300
    plt.rcParams['figure.dpi'] = dpi_ref
    plt.rcParams['image.cmap'] = 'jet'

    # Plot depth data
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.plot(range(1, num_wells+1), locZ, 'or', label='well depths')
    ax1.plot([1, num_wells], [z_min, z_min], '--k', label='zmin, zmax')
    ax1.plot([1, num_wells], [z_max, z_max], '--k')
    ax1.plot([1, num_wells], [reservoirDepth, reservoirDepth], '--g',
             label='reservoir depth')
    ax1.legend()
    ax1.invert_yaxis()
    ax1.set_xlabel('Well index')
    ax1.set_ylabel('Depth, [m]')

    # Check run frequency of wellbore components
    turn_off_wells_inds = []
    turn_on_wells_inds = []
    for comp_ref in range(num_wells):
        if msw[comp_ref].run_frequency == 0:
            turn_off_wells_inds.append(comp_ref)
        else:
            turn_on_wells_inds.append(comp_ref)
    turn_off_wells_inds = np.array(turn_off_wells_inds)
    turn_on_wells_inds = np.array(turn_on_wells_inds)

    print('Expected number of turned off wells:', num_wells-num_turn_on_wells)
    print('Obtained number of turned off wells:', len(turn_off_wells_inds))

    z_less_res_depth_inds = np.where(locZ < reservoirDepth)[0]

    # Print indices of turned off wells
    print(turn_off_wells_inds)
    print(z_less_res_depth_inds)

    print('__________________Results of forward simulation___________________')
    print('The simulation used {} wells with randomly generated x, y, and z.'.format(num_wells))
    print('Wells penetrating the reservoir: ')
    for comp_ref in turn_on_wells_inds:
        if locZ[comp_ref] >= reservoirDepth:

            print('comp_ref: ', comp_ref)
            print('x = {:.2f} m, y = {:.2f} m'.format(
                sm.obs['gen.locX{}_0'.format(comp_ref)].sim,
                sm.obs['gen.locY{}_0'.format(comp_ref)].sim))
            print('z = {:.2f} m'.format(sm.obs['gen.locZ{}_0'.format(comp_ref)].sim))
            print('    Pressure [Pa]: ',
                  sm.collect_observations_as_time_series(
                      sres[comp_ref], 'pressure'), sep='\n')
            print('    Brine Leakage Rate [kg/s]: ',
                  sm.collect_observations_as_time_series(
                      msw[comp_ref], 'brine_aquifer1'), sep='\n')
            print(' ')

    distance_range_km = [0, 0]
    # First, get the range of x and y values. The range is used for the colorbar below.
    for comp_ref in range(num_wells):
        # Making them both positive just for conceptual clarity
        if sm.obs['gen.locZ{}_0'.format(comp_ref)].sim >= reservoirDepth:
            distance_from_site_km = get_distance_from_site_km(
                sm, comp_ref, injection_x_m, injection_y_m)

            # If it's the first time, set the value
            if distance_range_km[0] == 0:
                distance_range_km[0] = distance_from_site_km

            if distance_range_km[1] == 0:
                distance_range_km[1] = distance_from_site_km

            if distance_from_site_km < distance_range_km[0]:
                distance_range_km[0] = distance_from_site_km

            if distance_from_site_km > distance_range_km[1]:
                distance_range_km[1] = distance_from_site_km

    min_level = distance_range_km[0]
    max_level = distance_range_km[1]

    if min_level == 0 and max_level == 0:
        max_level = 1

    else:
        if min_level == max_level:
            max_level *= 1.1

    interval = (max_level - min_level) / 100

    distance_levels_km = np.arange(min_level, max_level + interval, interval)

    # Figure formatting input

    axislabelfontsize = 10
    titlefontsize = 10
    lgndfontsize = 6

    labelfontweight = 'bold'
    colormap = 'RdYlBu_r'

    lineWidth = 1.5

    edgeWidth = 1
    edgeWidthSite = 2

    unusedWellColor = [0.33, 0.33, 0.33]
    unusedWellEdgeColor = [unusedWellColor[0] * 0.5, unusedWellColor[1] * 0.5,
                           unusedWellColor[2] * 0.5]

    unusedWellMarkerSize = 4
    usedWellMarkerSize = 7
    injectionSiteMarkerSize = 5

    injectionSiteMarker = 's'
    wellMarker = 'o'

    cmap = plt.cm.get_cmap(colormap)

    for comp_ref in range(num_wells):
        # Making them both positive just for conceptual clarity
        if sm.obs['gen.locZ{}_0'.format(comp_ref)].sim < reservoirDepth:
            # Well does not penetrate the reservoir
            plt.figure(2)
            plt.plot(sm.obs['gen.locX{}_0'.format(comp_ref)].sim / 1000,
                     sm.obs['gen.locY{}_0'.format(comp_ref)].sim / 1000,
                     linestyle='none', markeredgewidth=edgeWidth, marker=wellMarker,
                     markerfacecolor=unusedWellColor, markersize=unusedWellMarkerSize,
                     color=unusedWellEdgeColor, zorder=1)

        else:
            # Well penetrates the reservoir
            distance_from_site_km = get_distance_from_site_km(
                sm, comp_ref, injection_x_m, injection_y_m)

            rgba = cmap(
                (distance_from_site_km  - np.min(distance_range_km))
                / (np.max(distance_range_km) - np.min(distance_range_km)))

            markerfacecolor = rgba[0:3]
            markeredgecolor = (np.array(rgba[0:3]) * 0.5).tolist()

            plt.figure(2)
            plt.plot(sm.obs['gen.locX{}_0'.format(comp_ref)].sim / 1000,
                     sm.obs['gen.locY{}_0'.format(comp_ref)].sim / 1000,
                     linestyle='none', markeredgewidth=edgeWidth, marker=wellMarker,
                     markerfacecolor=markerfacecolor, markersize=usedWellMarkerSize,
                     color=markeredgecolor, zorder=2)

            pressure = sm.collect_observations_as_time_series(
                sres[comp_ref], 'pressure')

            plt.figure(3)
            plt.plot(time_array / 365.25, pressure / 1.0e+6, color=markerfacecolor,
                      linewidth=lineWidth)

            brine_aquifer1 = sm.collect_observations_as_time_series(
                msw[comp_ref], 'brine_aquifer1')

            plt.figure(4)
            plt.plot(time_array / 365.25, brine_aquifer1, color=markerfacecolor,
                      linewidth=lineWidth)

    fig = plt.figure(2)

    plt.axis('equal')
    plt.plot(injection_x_m / 1000, injection_y_m / 1000, color='k',
              linestyle='none', marker=injectionSiteMarker, markerfacecolor='none',
              markersize=injectionSiteMarkerSize, markeredgewidth=edgeWidthSite, zorder=3)

    # Manually create the legend items for figure 1
    fig1_handle_list = []
    fig1_label_list = []

    injectionSite = Line2D(
        [0], [0], color='k', markeredgewidth=edgeWidthSite, linestyle='none',
        marker=injectionSiteMarker, markerfacecolor='none',
        markersize=injectionSiteMarkerSize)

    fig1_handle_list.append(injectionSite)
    fig1_label_list.append('Injection Site')

    unusedWells = Line2D(
        [0], [0], linestyle='none', marker=wellMarker, markerfacecolor=unusedWellColor,
        markersize=unusedWellMarkerSize, color=unusedWellEdgeColor,
        markeredgewidth=edgeWidth)

    fig1_handle_list.append(unusedWells)
    fig1_label_list.append('Well, No Connection to Reservoir')

    rgba = cmap(1)

    usedWells = Line2D(
        [0], [0], linestyle='none', marker=wellMarker, markerfacecolor=rgba[0:3],
        markersize=usedWellMarkerSize, color=(np.array(rgba[0:3]) * 0.5).tolist(),
        markeredgewidth=edgeWidth)

    fig1_handle_list.append(usedWells)
    fig1_label_list.append('Well Penetrating the Reservoir')

    # Format plots
    fig = plt.figure(2)

    plt.xlabel('x, [km]', fontsize=axislabelfontsize, fontweight=labelfontweight)
    plt.ylabel('y, [km]', fontsize=axislabelfontsize, fontweight=labelfontweight)
    plt.title('Map-View Image of the Site', fontsize=titlefontsize,
              fontweight=labelfontweight)

    # Make the colorbar
    cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=plt.gca(),
                        values=distance_levels_km)
    cbar.set_label('Distance From Injection Site, [km]', rotation=90,
                   fontsize=axislabelfontsize, fontweight=labelfontweight)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()

    plt.xlim(x_min / 1000, x_max / 1000)
    plt.xlim(y_min / 1000, y_max / 1000)

    plt.subplots_adjust(left=0.1, bottom=0.2, right=0.7,
                        top=0.9, wspace=0.1, hspace=0.1)

    # Add Legend
    plt.gca().legend(fig1_handle_list, fig1_label_list, fancybox=False,
                     fontsize=lgndfontsize, ncol=1,
                     edgecolor=[0, 0, 0], loc='upper center',
                     bbox_to_anchor=(1.55, 1), framealpha=0.67)

    fig = plt.figure(3)

    plt.xlabel('Time, [years]', fontsize=axislabelfontsize,
               fontweight=labelfontweight)
    plt.ylabel('Pressure, [MPa]', fontsize=axislabelfontsize,
               fontweight=labelfontweight)
    plt.title('Reservoir Pressure Over Time', fontsize=titlefontsize,
              fontweight=labelfontweight)

    # Make the colorbar
    cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=plt.gca(),
                        values=distance_levels_km)
    cbar.set_label('Distance From Injection Site, [km]', rotation=90,
                   fontsize=axislabelfontsize,
                    fontweight=labelfontweight)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()

    plt.subplots_adjust(left=0.2, bottom=0.2, right=0.9,
                        top=0.9, wspace=0.1, hspace=0.1)

    fig = plt.figure(4)

    plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0), useMathText = True)

    plt.xlabel('Time, [years]', fontsize=axislabelfontsize,
               fontweight=labelfontweight)
    plt.ylabel('Brine Leakage Rate, [kg/s]', fontsize=axislabelfontsize,
               fontweight=labelfontweight)
    plt.title('Brine Leakage Over Time', fontsize=titlefontsize,
              fontweight=labelfontweight)

    # Make the colorbar
    cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=plt.gca(),
                        values=distance_levels_km)
    cbar.set_label('Distance From Injection Site, [km]', rotation=90,
                   fontsize=axislabelfontsize, fontweight=labelfontweight)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()

    plt.subplots_adjust(left=0.2, bottom=0.2, right=0.9,
                        top=0.9, wspace=0.1, hspace=0.1)
