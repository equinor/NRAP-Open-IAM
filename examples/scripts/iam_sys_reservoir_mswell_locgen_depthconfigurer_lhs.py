"""
This example demonstrates the use of the LocationGenerator and WellDepthRiskConfigurer
components. The x, y, and z values (easting, northing, and depth) for wellbores
are randomly created by the LocationGenerator. The well depths are then provided
to the WellDepthRiskConfigurer. If the well does not reach the top of the
reservoir, then the WellDepthRiskConfigurer disables the well by setting the
run_frequency to 0 (so the wellbore component never runs).

This simulation uses the Latin Hypercube Sampling (lhs) analysis type.

Note that randomly assigning coordinates with the LocationGenerator is different
from first creating random x, y, and z values and then directly providing those
to a reservoir component. The difference lies in the fact that the LocationGenerator
creates the x, y, and z values during the actual simulation. If you run a
simulation with Latin Hypercube Sampling, each realization will have different
x, y, and z values because they are generated during the realization. Having
random well locations generated in each realization can be important if you know
that well records are incomplete: this approach allows one to assess the potential
significance of unknown wellbores. Some of the unknown wells may not penetrate
the storage reservoir and may, therefore, be excluded as potential leakage pathways.
"""
# Author: Nate Mitchell
# Nathaniel.Mitchell@netl.doe.gov
# Contributor: Veronika Vasylkivska
# Veronika.Vasylkivska@netl.doe.gov


import sys
import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import ticker
from matplotlib.lines import Line2D

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, Stratigraphy, SimpleReservoir,
                     MultisegmentedWellbore, LocationGenerator,
                     WellDepthRiskConfigurer)


def get_distance_from_site_km(injection_x_m, injection_y_m, locX, locY):
    """
    Returns the current component's distance from the injection site.
    """
    distance_from_site_km = (
        (((injection_x_m - locX) ** 2) + ((injection_y_m - locY) ** 2))
        ** 0.5) / 1000

    return distance_from_site_km

if __name__ == "__main__":
    # Simulation input
    num_wells = 25

    # These are the realizations that will have map-view images of the site.
    # The randomly generated wells will be different in each realization.
    selected_realizations = [0, 1, 2]

    # This is used as the seed input for the LocationGenerator
    loc_seed = 100
    min_loc_seed = 1
    max_loc_seed = 1000

    # This is the seed used when running the LHS simulation
    lhs_seed_val = 589

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

    # Figure formatting input
    figsize_unused_wells = (10, 6)
    figsize_used_wells = (10, 6)

    dpi_ref = 200

    genfontsize = 8
    axislabelfontsize = 10
    titlefontsize = 10
    lgndfontsize = 8

    labelfontweight = 'bold'
    colormap = 'RdYlBu_r'

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

    reals_color = 'darkslateblue'
    reals_alpha = 0.8
    reals_linewidth = 2.5

    output_folder = 'output_locgen_depthconfig_lhs_{date_time_stamp}'.format(
        date_time_stamp=datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))

    IAM_output_dir = os.path.join(os.getcwd(), '..', '..', 'output')
    output_dir = os.path.join(IAM_output_dir, output_folder)

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
    gen.add_par('seed', value=loc_seed, min=min_loc_seed, max=max_loc_seed)

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
        sres[-1].add_par('logResPerm', value=-12, vary=False) # min=-12.5, max=-11.5)

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
        msw[-1].add_par('logWellPerm', value=-13.0, vary=False) # min=-13.5, max=-12.5)
        msw[-1].add_par('logAquPerm', value=-11.0, vary=False)

        # Add keyword arguments linked to the output provided by reservoir model
        for obs_nm in ['pressure', 'CO2saturation']:
            msw[-1].add_kwarg_linked_to_obs(obs_nm, sres[-1].linkobs[obs_nm])

        # Add observations of multisegmented wellbore component model
        for obs_nm in ['CO2_aquifer1', 'CO2_aquifer2', 'brine_aquifer1', 'brine_aquifer2']:
            msw[-1].add_obs(obs_nm)

    # Setup parameters of LHS simulations
    num_samples = 30
    ncpus = 3
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=lhs_seed_val)

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)

    # Load produced well coordinates data
    coords = {}
    for obs_nm in ['locX', 'locY', 'locZ']:
        coords[obs_nm] = np.zeros((num_samples, num_wells))
        for real_num in range(1, num_samples+1):
            coords[obs_nm][real_num-1, :] = \
                sm.collect_gridded_observations_as_time_series(
                    gen, obs_nm, output_dir, indices=[0], rlzn_number=real_num)[0]

    distance_range_km = [0, 0]

    # Load pressure data from sample set
    pressure = {}
    for comp_ref in range(num_wells):
        pressure[comp_ref] = s.collect_observations_as_time_series(
            cmpnt=sres[comp_ref], obs_nm='pressure')[
                'sres{}.pressure'.format(comp_ref)]

    # Load brine leakage data from sample set
    brine_aquifer1 = {}
    for comp_ref in range(num_wells):
        brine_aquifer1[comp_ref] = s.collect_observations_as_time_series(
            cmpnt=msw[comp_ref], obs_nm='brine_aquifer1')[
                'msw{}.brine_aquifer1'.format(comp_ref)]

    # First, get the range of x and y values for the selected realization. The
    # range is used for the colorbar below.
    for comp_ref in range(num_wells):

        # Go through all realizations
        for real_num in range(num_samples):
            # Extract the x, y, and z values for the realization chosen for Figure 1
            locX_current = coords['locX'][real_num, comp_ref]
            locY_current = coords['locY'][real_num, comp_ref]
            locZ_current = coords['locZ'][real_num, comp_ref]

            # Making them both positive just for conceptual clarity
            if locZ_current >= reservoirDepth:
                distance_from_site_km = get_distance_from_site_km(
                    injection_x_m, injection_y_m, locX_current, locY_current)

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
    elif min_level == max_level:
        max_level *= 1.1

    interval = (max_level - min_level) / 100

    distance_levels_km = np.arange(min_level, max_level + interval, interval)

    # Setup color map
    cmap = plt.cm.get_cmap(colormap)

    font = {'family': 'Arial',
            'weight': 'normal',
            'size': genfontsize}
    plt.rc('font', **font)

    plt.rcParams['figure.dpi'] = dpi_ref
    plt.rcParams['image.cmap'] = 'jet'

    # First, plot the map-view images for each of the selected realizations
    for real_num in selected_realizations:
        plt.figure(100 + real_num, dpi=dpi_ref)
        ax = plt.gca()

        plt.plot(injection_x_m / 1000, injection_y_m / 1000, color='k',
                 linestyle='none', marker=injectionSiteMarker,
                 markerfacecolor='none', markersize=injectionSiteMarkerSize,
                 markeredgewidth=edgeWidthSite, zorder=3)

        for comp_ref in range(num_wells):
            locX_current = coords['locX'][real_num, comp_ref]
            locY_current = coords['locY'][real_num, comp_ref]
            locZ_current = coords['locZ'][real_num, comp_ref]

            if locZ_current < reservoirDepth:
                # Well does not penetrate the reservoir
                ax.plot(locX_current / 1000, locY_current / 1000,
                        linestyle='none', markeredgewidth=edgeWidth,
                        marker=wellMarker, markerfacecolor=unusedWellColor,
                        markersize=unusedWellMarkerSize,
                        color=unusedWellEdgeColor, zorder=1)
            else:
                # Well penetrates the reservoir
                distance_from_site_km = get_distance_from_site_km(
                    injection_x_m, injection_y_m, locX_current, locY_current)

                rgba = cmap(
                    (distance_from_site_km  - np.min(distance_range_km))
                    / (np.max(distance_range_km) - np.min(distance_range_km)))

                markerfacecolor = rgba[0:3]
                markeredgecolor = (np.array(rgba[0:3]) * 0.5).tolist()

                ax.plot(locX_current / 1000, locY_current / 1000,
                        linestyle='none', markeredgewidth=edgeWidth,
                        marker=wellMarker, markerfacecolor=markerfacecolor,
                        markersize=usedWellMarkerSize,
                        color=markeredgecolor, zorder=2)

    # Now plot the output for all realizations and all components
    for comp_ref in range(num_wells):

        for real_num in range(num_samples):
            locX_current = coords['locX'][real_num, comp_ref]
            locY_current = coords['locY'][real_num, comp_ref]
            locZ_current = coords['locZ'][real_num, comp_ref]

            if locZ_current < reservoirDepth:
                # Well does not penetrate the reservoir
                plt.figure(2, figsize=figsize_unused_wells, dpi=dpi_ref)

                color = unusedWellColor
            else:
                # Well penetrates the reservoir
                plt.figure(3, figsize=figsize_used_wells, dpi=dpi_ref)

                distance_from_site_km = get_distance_from_site_km(
                    injection_x_m, injection_y_m, locX_current, locY_current)

                rgba = cmap(
                    (distance_from_site_km  - np.min(distance_range_km))
                    / (np.max(distance_range_km) - np.min(distance_range_km)))

                color = rgba[0:3]
                color = (np.array(rgba[0:3]) * 0.5).tolist()

            ax = plt.subplot(1, 2, 1)
            ax.plot(time_array / 365.25, pressure[comp_ref][real_num, :] / 1.0e+6,
                    color=color, linestyle='-', linewidth=reals_linewidth,
                    alpha=reals_alpha, zorder=100)

            ax = plt.subplot(1, 2, 2)
            ax.plot(time_array / 365.25, brine_aquifer1[comp_ref][real_num, :],
                    color=color, linestyle='-', linewidth=reals_linewidth,
                    alpha=reals_alpha, zorder=100)


    # Manually create the legend items for the map-view figures
    mapview_fig_handle_list = []
    mapview_fig_label_list = []

    injectionSite = Line2D(
        [0], [0], color='k', markeredgewidth=edgeWidthSite, linestyle='none',
        marker=injectionSiteMarker, markerfacecolor='none',
        markersize=injectionSiteMarkerSize)

    mapview_fig_handle_list.append(injectionSite)
    mapview_fig_label_list.append('Injection Site')

    unusedWells = Line2D(
        [0], [0], linestyle='none', marker=wellMarker, markerfacecolor=unusedWellColor,
        markersize=unusedWellMarkerSize, color=unusedWellEdgeColor,
        markeredgewidth=edgeWidth)

    mapview_fig_handle_list.append(unusedWells)
    mapview_fig_label_list.append('Well, No Connection to Reservoir')

    rgba = cmap(1)

    usedWells = Line2D(
        [0], [0], linestyle='none', marker=wellMarker, markerfacecolor=rgba[0:3],
        markersize=usedWellMarkerSize, color=(np.array(rgba[0:3]) * 0.5).tolist(),
        markeredgewidth=edgeWidth)

    mapview_fig_handle_list.append(usedWells)
    mapview_fig_label_list.append('Well Penetrating the Reservoir')

    # Format and save plots
    for real_num in selected_realizations:
        fig = plt.figure(100 + real_num)
        plt.axis('equal')

        plt.xlabel('x, [km]', fontsize=axislabelfontsize, fontweight=labelfontweight)
        plt.ylabel('y, [km]', fontsize=axislabelfontsize, fontweight=labelfontweight)
        plt.title('Map-View Image of the Site,\nRealization {}'.format(real_num),
                  fontsize=titlefontsize, fontweight=labelfontweight)

        # Make the colorbar
        cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=plt.gca(),
                            values=distance_levels_km)
        cbar.set_label('Distance From Injection Site, [km]', rotation=90,
                       fontsize=axislabelfontsize,
                       fontweight=labelfontweight)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()

        plt.xlim(x_min / 1000, x_max / 1000)
        plt.ylim(y_min / 1000, y_max / 1000)

        plt.subplots_adjust(left=0.1, bottom=0.2, right=0.7,
                            top=0.9, wspace=0.1, hspace=0.1)

        # Add Legend
        plt.gca().legend(mapview_fig_handle_list, mapview_fig_label_list,
                         fancybox=False, fontsize=lgndfontsize-1, ncol=1,
                         edgecolor=[0, 0, 0], loc='upper center',
                         bbox_to_anchor=(1.55, 1), framealpha=0.67)

        plt.savefig(os.path.join(
            output_dir, 'Wells_Map_View_Realizaiton{}.png'.format(real_num)),
            dpi=dpi_ref)

    # All Realizations, Only Wells Without a Connection to the Reservoir
    fig = plt.figure(2)

    ax = plt.subplot(1, 2, 1)

    plt.xlabel('Time, [years]', fontsize=axislabelfontsize, fontweight=labelfontweight)
    plt.ylabel('Pressure, [MPa]', fontsize=axislabelfontsize, fontweight=labelfontweight)

    ax = plt.subplot(1, 2, 2)

    plt.xlabel('Time, [years]', fontsize=axislabelfontsize, fontweight=labelfontweight)
    plt.ylabel('Brine Leakage Rate, [kg/s]', fontsize=axislabelfontsize,
               fontweight=labelfontweight)

    plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0),
                         useMathText = True)

    plt.suptitle('Results Over Time, Only Wells Without a Connection to the Reservoir',
                 fontsize=titlefontsize, fontweight=labelfontweight)

    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.95,
                        top=0.9, wspace=0.3, hspace=0.1)

    plt.savefig(os.path.join(output_dir, 'Results_Over_Time_Wells_No_Connection.png'),
                dpi=dpi_ref)


    # All Realizations, Only Wells Penetrating the Reservoir
    fig = plt.figure(3)

    ax = plt.subplot(1, 2, 1)

    plt.xlabel('Time, [years]', fontsize=axislabelfontsize, fontweight=labelfontweight)
    plt.ylabel('Pressure, [MPa]', fontsize=axislabelfontsize, fontweight=labelfontweight)

    # Make the colorbar
    cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=plt.gca(),
                        values=distance_levels_km)
    cbar.set_label('Distance From Injection Site, [km]', rotation=90,
                   fontsize=axislabelfontsize, fontweight=labelfontweight)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()

    ax = plt.subplot(1, 2, 2)

    plt.xlabel('Time, [years]', fontsize=axislabelfontsize,
               fontweight=labelfontweight)
    plt.ylabel('Brine Leakage Rate, [kg/s]', fontsize=axislabelfontsize,
               fontweight=labelfontweight)

    plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0),
                         useMathText = True)

    # Make the colorbar
    cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=plt.gca(),
                        values=distance_levels_km)
    cbar.set_label('Distance From Injection Site, [km]', rotation=90,
                   fontsize=axislabelfontsize,
                   fontweight=labelfontweight)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()

    plt.suptitle('Results Over Time, Only Wells Penetrating the Reservoir',
                 fontsize=titlefontsize,
                 fontweight=labelfontweight)

    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.95,
                        top=0.9, wspace=0.3, hspace=0.1)

    plt.savefig(os.path.join(
        output_dir, 'Results_Over_Time_Wells_Penetrating_Reservoir.png'),
        dpi=dpi_ref)
