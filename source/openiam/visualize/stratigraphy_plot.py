"""
Code to create figures showing the stratigraphy created with
strata.py.

Examples illustrating applications or setup of stratigraphy_plot method:
    ControlFile_ex32b.yaml
    ControlFile_ex32c.yaml
    ControlFile_ex33a.yaml
    ControlFile_ex33b.yaml
    ControlFile_ex34.yaml
    ControlFile_ex35.yaml
    ControlFile_ex36.yaml
    ControlFile_ex37.yaml
    ControlFile_ex38.yaml

Created: September 15th, 2022
Last Modified: February 10th, 2023

@author: Nate Mitchell (Nathaniel.Mitchell@NETL.DOE.GOV)
@contributor: Veronika Vasylkivska (Veronika.Vasylkivska@NETL.DOE.GOV)
LRST (Battelle|Leidos) supporting NETL
"""

import os
import sys
import logging
import csv

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as clrs

SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(SOURCE_DIR)

import openiam.cfi.commons as iamcommons
import openiam.cfi.strata as iam_strata

reservoir_components = ['LookupTableReservoir',
                        'SimpleReservoir',
                        'AnalyticalReservoir']

wellbore_components = ['MultisegmentedWellbore',
                       'CementedWellbore',
                       'OpenWellbore',
                       'GeneralizedFlowRate']


def stratigraphy_plot(yaml_data, model_data, sm,
                      name='strata_Figure1', savefig=None, title=None,
                      figsize=(12, 10), figdpi=100, genfontsize=12,
                      axislabelfontsize=14, titlefontsize=14,
                      boldlabels=True, plot_wellbore_locations=True,
                      view_elev=None, view_azimuth=None,
                      LUTS_max_z=0, plot_well_labels=True,
                      plot_SandD_symbol=True, SandD_location=None,
                      plot_indiv_strat_comps=False, save_stratigraphy=False,
                      plot_injection_sites=False,
                      plot_injection_site_labels=False,
                      EnforceXandYLims=False, EnforceXandYGridLims=False,
                      Enforce_SandD_symbol_length=False, numberXTicks=5,
                      numberYTicks=5):
    """
    Makes a figure showing the domain's stratigraphy. This function works for
    flat-lying stratigraphy (no dip) or stratigraphy with a specified strike and
    dip.

    :param yaml_data: Dictionary of input values from the entire .yaml file
    :type yaml_data: dict

    :param model_data: Input from the 'ModelParams' section of the .yaml file
    :type model_data: dict

    :param sm: OpenIAM System model for which plots are created
    :type sm: openiam.SystemModel object

    :param name: Figure Name to be used/created.
    :type name: str

    :param savefig: Filename to save resulting figure to. No file saved if None.
    :type savefig: str

    :param title: Optional Title for figure
    :type title: str

    :param figsize: width and height of the figure (width, height) in inches.
        Default value is (12, 10).
    :type figsize: tuple

    :param figdpi: dpi (dots-per-inch) for the figure
    :type figdpi: float or int

    :param genfontsize: fontsize for tick labels, etc.
    :type genfontsize: float or int

    :param axislabelfontsize: fontsize for x and y axis labels
    :type axislabelfontsize: float or int

    :param titlefontsize: fontsize for the title
    :type titlefontsize: float or int

    :param colormap: string designation for a particular colormap (e.g., 'viridis')
    :type colormap: str

    :param boldlabels: option to use bold x and y labels and bold titles. Set to
        True for bold labels, False for normal labels.
    :type boldlabels: bool

    :param plot_wellbore_locations: option to plot each wellbore component. The
        component is shown as a vertical line extending from the well base to
        the surface. OpenWellbores have a wellTop parameter representing the
        depth above which the wellbore is no longer open. For example, the
        cement used to seal the well may be damaged beneath the wellTop depth,
        allowing leakage through. If an OpenWellbore has a wellTop value of
        zero, then the well is not effectively sealed between the well bottom
        (reservoirDepth parameter) and the surface. For an OpenWellbore with a
        nonzero wellTop value, a blue vertical line is plotted from the
        reservoirDepth to the wellTop. Then, a light blue vertical line is
        plotted from the wellTop depth to the surface.
    :type plot_wellbore_locations: bool

    :param view_elev: Value or values for the angle setting the 'elevation' in
        the 3D plot. Note that this value is in degrees, with higher angles
        moving the perspective futher above the area. Multiple values can be
        provided within a list, although the length should also match that of
        the view_azimuth input.
    :type view_elev: list, int, or float

    :param view_azimuth: Value or values for the angle setting the 'azimuth' in
        the 3D plot. Note that this value is in degrees. Multiple values can be
        provided within a list, although the length should also match that of
        the view_elev input.
    :type view_azimuth: list, int, or float

    :param LUTS_max_z: Upper limit used for the z axis when the
        LookupTableStratigrapy option is used. Default value is 100 m.
    :type LUTS_max_z: int or float

    :param plot_well_labels: option to the name of each wellbore
    :type plot_well_labels: bool

    :param plot_SandD_symbol: option to plot a strike and dip symbol
    :type plot_SandD_symbol: bool

    :param SandD_location: tuple containing the x and y coordinates at which the
        strike and dip symbol will be plotted. The first value is for x, the
        second is for y. The default value is None.
    :type SandD_location: tuple or None

    :param plot_indiv_strat_comps: option to plot the stratigraphy component
        created for each location within the control file interface. This
        option is provided to help verify that the functions within
        strata.py are behaving as expected. The tops of each unit
        are plotted as circles. Judging the match between the 3d planes and the
        points can be made easier by adjusting the view_elev and view_azimuth
        values (e.g., view_elev = 0 and view_azimuth is oriented along strike).
        This option is here to help check the data.
    :type plot_indiv_strat_comps: bool

    :param save_stratigraphy: option to save the thicknesses and depths for all
        units within .csv files. Values are saved for the grid locations used
        to make the 3D planes, x_loc and y_loc. The x and y values are also
        saved in separate files. This option is here to help check the data.
    :type save_stratigraphy: bool

    :param plot_injection_sites: option to plot the injection site location(s).
        Note that if using a Lookup Table Reservoir, InjectionCoordx and
        InjectionCoordy values must be provided in the .yaml file - see Control
        File example 35.
    :type plot_injection_sites: bool

    :param plot_injection_site_labels: option to display the text "Injection
        Site" by each plotted injection site (if plot_injection_sites is True).
    :type plot_injection_site_labels: bool

    :param EnforceXandYLims: option to enforce particular limits for the x and
        y axes. This option is enabled by having 'SpecifyXandYLims' in the
        .yaml file. The limits are set by xLims and yLims values listed under
        'SpecifyXandYLims' (i.e., on the next two lines and indented 4 spaces).
        The xLims and yLims should each be a list of length 2 (e.g., [0, 1000])
        with units of meters.
    :type EnforceXandYLims: bool

    :param EnforceXandYGridLims: option to enforce particular x and y limits
        for the 3D planes used to portray the tops of each unit. This option is
        enabled by having 'SpecifyXandYGridLims' in the .yaml file. The
        limits are set by gridXLims and gridYLims values listed under
        'SpecifyXandYGridLims'. The gridXLims and gridYLims should each
        be a list of length 2 (e.g., [100, 900]) with units of meters.
    :type EnforceXandYGridLims: bool

    :param Enforce_SandD_symbol_length: option to enforce particular length for
        the strike and dip symbol. This option is enabled by having 'length'
        under the 'StrikeAndDipSymbol' section within the Stratigraphy plot of
        a .yaml file. The length is given in meters.
    :type Enforce_SandD_symbol_length: bool

    :param numberXTicks: number of ticks along the x axis. Used only if the x
        axis limits are enforced.
    :type numberXTicks: int

    :param numberYTicks: number of ticks along the y axis. Used only if the y
        axis limits are enforced.
    :type numberYTicks: int

    :return: None
    """
    if view_elev is None:
        view_elev = [10, 30]

    if view_azimuth is None:
        view_azimuth = [300, 300]

    if boldlabels:
        selected_labelfontweight = 'bold'
    else:
        selected_labelfontweight = 'normal'

    # Get the stratigrapy information from the .yaml file
    strata_var_info = iam_strata.get_strata_var_info_from_yaml(yaml_data)

    var_type = strata_var_info['var_type']
    strike = strata_var_info['strike']
    dip = strata_var_info['dip']
    dipDirection = strata_var_info['dipDirection']
    coordxReferencePoint = strata_var_info['coordxReferencePoint']
    coordyReferencePoint = strata_var_info['coordyReferencePoint']

    if coordxReferencePoint is None:
        coordxReferencePoint = 0

    if coordyReferencePoint is None:
        coordyReferencePoint = 0

    if var_type == 'strikeAndDip':
        dipDirectionDegrees = iam_strata.obtain_dip_direction_degrees(
            strike, dipDirection)

    # Get the stratigraphy plot input from the .yaml file
    strat_plot_yaml_input = read_strata_plot_yaml_input(yaml_data, name)

    if strat_plot_yaml_input['dpi_input'] is not None:
        figdpi = float(strat_plot_yaml_input['dpi_input'])

    if strat_plot_yaml_input['EnforceXandYLims'] is not None:
        EnforceXandYLims = strat_plot_yaml_input['EnforceXandYLims']
        xLims = strat_plot_yaml_input['xLims']
        yLims = strat_plot_yaml_input['yLims']

    if strat_plot_yaml_input['EnforceXandYGridLims'] is not None:
        EnforceXandYGridLims = strat_plot_yaml_input['EnforceXandYGridLims']
        gridXLims = strat_plot_yaml_input['gridXLims']
        gridYLims = strat_plot_yaml_input['gridYLims']

    if strat_plot_yaml_input['xGridSpacing'] is not None:
        x_grid_spacing = strat_plot_yaml_input['xGridSpacing']
        set_xgrid_spacing = True
    else:
        set_xgrid_spacing = False

    if strat_plot_yaml_input['yGridSpacing'] is not None:
        y_grid_spacing = strat_plot_yaml_input['yGridSpacing']
        set_ygrid_spacing = True
    else:
        set_ygrid_spacing = False

    if strat_plot_yaml_input['plot_injection_sites'] is not None:
        plot_injection_sites = strat_plot_yaml_input['plot_injection_sites']

    if strat_plot_yaml_input['InjectionCoordx'] is not None and \
            strat_plot_yaml_input['InjectionCoordy'] is not None:
        InjectionCoordx = strat_plot_yaml_input['InjectionCoordx']
        InjectionCoordy = strat_plot_yaml_input['InjectionCoordy']

    if strat_plot_yaml_input['plot_injection_site_labels'] is not None:
        plot_injection_site_labels = strat_plot_yaml_input['plot_injection_site_labels']

    if strat_plot_yaml_input['plot_wellbore_locations'] is not None:
        plot_wellbore_locations = strat_plot_yaml_input['plot_wellbore_locations']

    if strat_plot_yaml_input['plot_well_labels'] is not None:
        plot_well_labels = strat_plot_yaml_input['plot_well_labels']

    if strat_plot_yaml_input['plot_indiv_strat_comps'] is not None:
        plot_indiv_strat_comps = strat_plot_yaml_input['plot_indiv_strat_comps']

    if strat_plot_yaml_input['save_stratigraphy'] is not None:
        save_stratigraphy = strat_plot_yaml_input['save_stratigraphy']

    if strat_plot_yaml_input['plot_SandD_symbol'] is not None:
        plot_SandD_symbol = strat_plot_yaml_input['plot_SandD_symbol']

    if strat_plot_yaml_input['SandD_location'] is not None:
        SandD_location = strat_plot_yaml_input['SandD_location']

    if strat_plot_yaml_input['Enforce_SandD_symbol_length'] is not None:
        Enforce_SandD_symbol_length = strat_plot_yaml_input['Enforce_SandD_symbol_length']
        SandD_symbol_length = strat_plot_yaml_input['SandD_symbol_length']

    if strat_plot_yaml_input['view_elev']  is not None and \
            strat_plot_yaml_input['view_azimuth']  is not None:
        view_elev = strat_plot_yaml_input['view_elev']
        view_azimuth = strat_plot_yaml_input['view_azimuth']

    # If LookupTableStratigrapy is being used, a component's location is needed
    # later in the code (any component). This location is used to read the table
    # - using an invalid location (too far away) would produce an error.
    components = list(sm.component_models.values())
    if var_type == 'LookupTable':
        plot_indiv_strat_comps = True
        LUTS_loc_check = True
        for comp in components:
            if LUTS_loc_check:
                if comp.class_type in reservoir_components:
                    x_vals = comp.locX
                    y_vals = comp.locY

                    if LUTS_loc_check:
                        x_val_LUTS = x_vals
                        y_val_LUTS = y_vals
                        LUTS_loc_check = False

    # Create the x and y values used for the 3D planes
    if EnforceXandYGridLims:
        min_x_val = gridXLims[0]
        max_x_val = gridXLims[1]

        min_y_val = gridYLims[0]
        max_y_val = gridYLims[1]

        x_buffer = 0
        y_buffer = 0

    else:
        # These are overwritten within the loop below
        min_x_val = 9e99
        max_x_val = 0
        min_y_val = 9e99
        max_y_val = 0

        for comp in components:
            if comp.class_type in reservoir_components:
                x_vals = comp.locX
                y_vals = comp.locY

                if np.min(x_vals) < min_x_val:
                    min_x_val = np.min(x_vals)

                if np.max(x_vals) > max_x_val:
                    max_x_val = np.max(x_vals)

                if np.min(y_vals) < min_y_val:
                    min_y_val = np.min(y_vals)

                if np.max(y_vals) > max_y_val:
                    max_y_val = np.max(y_vals)

            if comp.class_type in reservoir_components and plot_injection_sites:
                if comp.class_type == 'LookupTableReservoir':
                    x_vals = InjectionCoordx
                    y_vals = InjectionCoordy

                    try:
                        x_vals = float(x_vals)
                        y_vals = float(y_vals)

                    except TypeError:
                        x_vals = np.array(list(x_vals), dtype=float)
                        y_vals = np.array(list(y_vals), dtype=float)

                else:
                    x_vals = comp.injX
                    y_vals = comp.injY

                if np.min(x_vals) < min_x_val:
                    min_x_val = np.min(x_vals)

                if np.max(x_vals) > max_x_val:
                    max_x_val = np.max(x_vals)

                if np.min(y_vals) < min_y_val:
                    min_y_val = np.min(y_vals)

                if np.max(y_vals) > max_y_val:
                    max_y_val = np.max(y_vals)

        # If there's only one point, the plot will look weird. Here, I make
        # the surfaces extend a minimum of 200 m in the x direction and
        # 200 m in the y direction.
        if (max_x_val - min_x_val) < 200:
            max_x_val += 100
            min_x_val -= 100

        if (max_y_val - min_y_val) < 200:
            max_y_val += 100
            min_y_val -= 100

        x_buffer = (max_x_val - min_x_val) / 10
        y_buffer = (max_y_val - min_y_val) / 10

    if set_xgrid_spacing:
        num_x_points = int(np.ceil((max_x_val - min_x_val)
                                   / x_grid_spacing) + 1)
    else:
        num_x_points = 101

    if set_ygrid_spacing:
        num_y_points = int(np.ceil((max_y_val - min_y_val)
                                   / y_grid_spacing) + 1)
    else:
        num_y_points = 101

    uniq_x = np.linspace(min_x_val - x_buffer,
                         max_x_val + x_buffer,
                         num_x_points,
                         endpoint = True)
    uniq_y = np.linspace(min_y_val - y_buffer,
                         max_y_val + y_buffer,
                         num_y_points,
                         endpoint = True)

    x_loc, y_loc = np.meshgrid(uniq_x, uniq_y)

    if SandD_location:
        reset_SandD_location = check_SandD_location(SandD_location,
                                                    min_x_val, max_x_val,
                                                    min_y_val, max_y_val,
                                                    x_buffer, y_buffer)
        if reset_SandD_location:
            SandD_location = None

    if plot_SandD_symbol:
        if SandD_location is None:
            # I don't want the symbol to be right in the middle of the area b/c
            # there will often be a well there. Having the well right beneath
            # the symbol makes it seem like the symbol is representing something else.
            SandD_loc_x = np.mean([max_x_val, min_x_val]) - ((max_x_val - min_x_val) / 4)
            SandD_loc_y = np.mean([max_y_val, min_y_val]) + ((max_y_val - min_y_val) / 4)
            SandD_location = [SandD_loc_x, SandD_loc_y]

    if var_type == 'LookupTable':
        file_name = yaml_data['Stratigraphy']['spatiallyVariable'][
            'LookupTableStratigraphy']['FileName']
        file_directory = yaml_data['Stratigraphy']['spatiallyVariable'][
            'LookupTableStratigraphy']['FileDirectory']

        LUTStrat_dict = iam_strata.get_lut_stratigraphy_dict(
            file_name, file_directory, x_val_LUTS, y_val_LUTS)

        numShaleLayers = LUTStrat_dict['numberOfShaleLayers']

    else:
        if var_type == 'strikeAndDip':
            strataReferencePoint = sm.component_models['strataRefPoint']
        else:
            strataReferencePoint = sm.component_models['strata']

        strata_dict = iam_strata.get_strata_info_from_component(strataReferencePoint)

        numShaleLayers = strata_dict['numberOfShaleLayers']

    reservoirColor, reservoirAlpha, reservoirAlphaFill, reservoirLabel, \
        shaleColor, shaleAlpha, shaleAlphaFill, shaleLabel, \
            aquiferColor, aquiferAlpha, aquiferAlphaFill, aquiferLabel, \
                wellColor, wellAlpha, wellAlphaFill, wellLabel = \
                    iam_strata.get_plotting_setup_for_units(
                        strat_plot_yaml_input, numShaleLayers)

    # Make the figure
    font = {'family': 'Arial',
            'weight': 'normal',
            'size': genfontsize}
    plt.rc('font', **font)

    # 3D Stratigraphy Figure
    plt3d = plt.figure(figsize=figsize, dpi=figdpi)
    plt3d_ax = plt3d.add_subplot(projection='3d')

    if var_type == 'LookupTable':
        plot_SandD_symbol = False
        save_stratigraphy = False
        plot_indiv_strat_comps = True

        first_time_check = True
        for comp in components:
            if comp.class_type in reservoir_components:
                strat_comp_temp = sm.component_models['strata' + comp.name]

                maxDepth_temp = iam_strata.get_unit_depth_from_component(
                    numShaleLayers, strat_comp_temp,
                    unitType='reservoir', top_or_bottom='bottom')

                if first_time_check:
                    maxDepth = maxDepth_temp
                    first_time_check = False
                elif maxDepth_temp > maxDepth:
                    maxDepth = maxDepth_temp

        z0_plane = np.zeros((x_loc.shape[0], x_loc.shape[1]))
        plt3d_ax.plot_surface(x_loc / 1000, y_loc / 1000, z0_plane,
                              color=shaleColor[-1], alpha=shaleAlphaFill[-1], shade=False)

        plt3d_ax.text(x_loc[0, 0] / 1000, y_loc[0, 0] / 1000, z0_plane[0, 0],
                      'Shale ' + str(numShaleLayers), zdir='y', color=shaleColor[-1],
                      alpha=shaleAlpha[-1], fontsize=genfontsize - 2, fontweight='bold')

    else:
        shaleThicknessesReferencePoint = strata_dict['shaleThicknesses']

        aquiferThicknessesReferencePoint = strata_dict['aquiferThicknesses']

        reservoirThicknessReferencePoint = strata_dict['reservoirThickness']

        # Dictionary containing the arrays used for the 3D graph
        stratigraphy_by_loc = create_strata_planes_dict(
            x_loc, y_loc, numShaleLayers,
            shaleThicknessesReferencePoint,
            aquiferThicknessesReferencePoint,
            reservoirThicknessReferencePoint,
            var_type, strataReferencePoint,
            coordxReferencePoint=coordxReferencePoint,
            coordyReferencePoint=coordyReferencePoint,
            strike=strike, dip=dip,
            dipDirection=dipDirection)

        # Plot the top of each unit as a 3D surface
        surface = -stratigraphy_by_loc['resBottomDepth'][:, :]
        plt3d_ax.plot_surface(x_loc / 1000, y_loc / 1000, surface,
                              color=reservoirColor, alpha=reservoirAlphaFill,
                              shade=False)

        surface = -stratigraphy_by_loc['resTopDepth'][:, :]
        plt3d_ax.plot_surface(x_loc / 1000, y_loc / 1000, surface,
                              color=reservoirColor, alpha=(reservoirAlphaFill / 2),
                              shade=False)
        plt3d_ax.text(x_loc[0, 0] / 1000, y_loc[0, 0] / 1000, surface[0, 0],
                      reservoirLabel, zdir='y', color=reservoirColor,
                      alpha=reservoirAlpha, fontsize=genfontsize - 2,
                      fontweight='bold')

        for shaleRef in range(1, numShaleLayers + 1):
            surface = -stratigraphy_by_loc['shale{}TopDepth'.format(shaleRef)][:, :]
            plt3d_ax.plot_surface(x_loc / 1000, y_loc / 1000, surface,
                                  color=shaleColor[shaleRef - 1],
                                  alpha=shaleAlphaFill[shaleRef - 1], shade=False)

            plt3d_ax.text(x_loc[0, 0] / 1000, y_loc[0, 0] / 1000, surface[0, 0],
                          shaleLabel[shaleRef - 1], zdir='y', color=shaleColor[shaleRef - 1],
                          alpha=shaleAlpha[shaleRef - 1],
                          fontsize=genfontsize - 2, fontweight='bold')

            if shaleRef < numShaleLayers:
                surface = -stratigraphy_by_loc['aquifer{}TopDepth'.format(shaleRef)][:, :]
                plt3d_ax.plot_surface(x_loc / 1000, y_loc / 1000, surface,
                                      color=aquiferColor[shaleRef - 1],
                                      alpha=aquiferAlphaFill[shaleRef - 1],
                                      shade=False)

                plt3d_ax.text(x_loc[0, 0] / 1000, y_loc[0, 0] / 1000, surface[0, 0],
                              aquiferLabel[shaleRef - 1], zdir='y',
                              color=aquiferColor[shaleRef - 1],
                              alpha=aquiferAlpha[shaleRef - 1],
                              fontsize=genfontsize - 2, fontweight='bold')

    if var_type == 'LookupTable':
        zlim = [-maxDepth, LUTS_max_z]
    else:
        zlim = plt3d_ax.get_zlim()

    if plot_injection_sites:
        for comp in components:
            if comp.class_type in reservoir_components:
                if comp.class_type == 'LookupTableReservoir':
                    x_vals_res = InjectionCoordx
                    y_vals_res = InjectionCoordy

                    try:
                        x_vals_res = float(x_vals_res)
                        y_vals_res = float(y_vals_res)
                    except TypeError:
                        x_vals_res = np.array(list(x_vals_res), dtype=float)
                        y_vals_res = np.array(list(y_vals_res), dtype=float)

                else:
                    x_vals_res = float(comp.injX)
                    y_vals_res = float(comp.injY)

                if not isinstance(x_vals_res, float):
                    for x_res_ref, x_res_val in enumerate(x_vals_res):
                        if var_type == 'strikeAndDip':
                            # Get the reservoir top depth
                            updatedStratigraphy = \
                                iam_strata.update_stratigraphy_by_strike_and_dip(
                                    numberOfShaleLayers=numShaleLayers,
                                    shaleThicknessList=shaleThicknessesReferencePoint[:],
                                    aquiferThicknessList=aquiferThicknessesReferencePoint[:],
                                    reservoirThickness=reservoirThicknessReferencePoint,
                                    strike=strike, dip=dip, dipDirection=dipDirection,
                                    coordxRefPoint=coordxReferencePoint,
                                    coordyRefPoint=coordyReferencePoint,
                                    location_x=x_res_val,
                                    location_y=y_vals_res[x_res_ref],
                                    strataRefPoint=strataReferencePoint)

                            resTopDepth_temp = updatedStratigraphy['reservoirTopDepth']

                        elif var_type == 'noVariation':
                            resTopDepth_temp = sum(shaleThicknessesReferencePoint) \
                                + sum(aquiferThicknessesReferencePoint)
                        elif var_type == 'LookupTable':
                            # Don't try to plot a vertical line from the
                            # surface to the reservoir when using
                            # LookupTable Stratigraphy.
                            pass

                        # Plot a 'shadow' beneath the point
                        plt3d_ax.plot(x_res_val / 1000,
                                      y_vals_res[x_res_ref] / 1000,
                                      zlim[0], marker='s', markersize=3,
                                      color=[0.67, 0.67, 0.67], linewidth=1)

                        if var_type != 'LookupTable':
                            plt3d_ax.plot([x_res_val / 1000,
                                           x_res_val / 1000],
                                          [y_vals_res[x_res_ref] / 1000,
                                           y_vals_res[x_res_ref] / 1000],
                                          [-resTopDepth_temp, 0], marker='s',
                                          markersize=3,
                                          color=[0.25, 0.25, 0.25], linewidth=1)

                        plt3d_ax.plot(x_res_val / 1000,
                                      y_vals_res[x_res_ref] / 1000,
                                      0, marker='s', markersize=3,
                                      color=[0.25, 0.25, 0.25], linewidth=1,
                                      zorder=98)

                        if plot_injection_site_labels and \
                                x_res_ref == (len(x_vals_res) - 1):
                            plt3d_ax.text(x_res_val / 1000,
                                          y_vals_res[x_res_ref] / 1000, 0,
                                          'Injection\nSites', zdir='x', color='k',
                                          fontsize=genfontsize - 4, fontweight='bold',
                                          zorder=99)

                else:
                    if var_type == 'strikeAndDip':
                        # Get the reservoir top depth
                        updatedStratigraphy = \
                            iam_strata.update_stratigraphy_by_strike_and_dip(
                                numberOfShaleLayers=numShaleLayers,
                                shaleThicknessList=shaleThicknessesReferencePoint[:],
                                aquiferThicknessList=aquiferThicknessesReferencePoint[:],
                                reservoirThickness=reservoirThicknessReferencePoint,
                                strike=strike, dip=dip, dipDirection=dipDirection,
                                coordxRefPoint=coordxReferencePoint,
                                coordyRefPoint=coordyReferencePoint,
                                location_x=x_vals_res,
                                location_y=y_vals_res,
                                strataRefPoint=strataReferencePoint)

                        resTopDepth_temp = updatedStratigraphy['reservoirTopDepth']

                    elif var_type == 'noVariation':
                        resTopDepth_temp = sum(shaleThicknessesReferencePoint) \
                            + sum(aquiferThicknessesReferencePoint)
                    elif var_type == 'LookupTable':
                        # Don't try to plot a vertical line from the
                        # surface to the reservoir when using
                        # LookupTable Stratigraphy.
                        pass

                    # Plot a 'shadow' beneath the point
                    plt3d_ax.plot(x_vals_res / 1000, y_vals_res / 1000,
                                  zlim[0], marker='s', markersize=3,
                                  color=[0.67, 0.67, 0.67], linewidth=1)

                    if var_type != 'LookupTable':
                        plt3d_ax.plot([x_vals_res / 1000, x_vals_res / 1000],
                                      [y_vals_res / 1000, y_vals_res / 1000],
                                      [-resTopDepth_temp, 0], marker='s',
                                      markersize=3, color=[0.25, 0.25, 0.25],
                                      linewidth=1)

                    plt3d_ax.plot(x_vals_res / 1000, y_vals_res / 1000,
                                  0, marker='s', markersize=3,
                                  color=[0.25, 0.25, 0.25], linewidth=1,
                                  zorder=98)

                    if plot_injection_site_labels:
                        plt3d_ax.text(x_vals_res / 1000, y_vals_res / 1000, 0,
                                      'Injection\nSite', zdir='x', color='k',
                                      fontsize=genfontsize - 4, fontweight='bold',
                                      zorder=99)

    if plot_wellbore_locations:
        for comp in components:
            if comp.class_type in wellbore_components:
                res_comp = sm.component_models[yaml_data[comp.name]['Connection']]
                x_vals_well = res_comp.locX
                y_vals_well = res_comp.locY

                if not isinstance(x_vals_well, (int, float)):
                    x_vals_well = float(x_vals_well)

                if not isinstance(y_vals_well, (int, float)):
                    y_vals_well = float(y_vals_well)

                if comp.class_type == 'OpenWellbore':
                    number = int(comp.name[(comp.name.index('_') + 1):None])
                    compName = 'Open\nWellbore {}'.format(number)

                    z_min = iamcommons.get_parameter_val(comp, 'reservoirDepth',
                                                         sm=sm, yaml_data=yaml_data)
                    z_max = iamcommons.get_parameter_val(comp, 'wellTop',
                                                         sm=sm, yaml_data=yaml_data)

                elif comp.class_type == 'MultisegmentedWellbore':
                    number = int(comp.name[(comp.name.index('_') + 1):None])
                    compName = 'M.S.\nWellbore {}'.format(number)

                    if var_type in ['strikeAndDip', 'LookupTable']:
                        strata_temp = sm.component_models['strata' + comp.name]

                        z_min = iam_strata.get_unit_depth_from_component(
                            numShaleLayers, strata_temp,
                            unitType='reservoir', top_or_bottom='top')

                    elif var_type == 'noVariation':
                        z_min=iam_strata.get_unit_depth_from_component(
                            numShaleLayers, strataReferencePoint,
                            unitType='reservoir', top_or_bottom='top')

                    z_max = 0

                elif comp.class_type == 'CementedWellbore':
                    number = int(comp.name[(comp.name.index('_') + 1):None])
                    compName = 'Cemented\nWellbore {}'.format(number)

                    z_min = iamcommons.get_parameter_val(comp, 'wellDepth',
                                                         sm=sm, yaml_data=yaml_data)
                    z_max = 0

                # Plot a 'shadow' beneath the point
                plt3d_ax.plot(x_vals_well / 1000, y_vals_well / 1000,
                              zlim[0], marker='o', markersize=3,
                              color=[0.67, 0.67, 0.67], linewidth=1)

                plt3d_ax.plot([x_vals_well / 1000, x_vals_well / 1000],
                              [y_vals_well / 1000, y_vals_well / 1000],
                              [-z_min, -z_max], marker='o', markersize=3,
                              color=wellColor, alpha=wellAlphaFill, linewidth=1)

                if plot_well_labels:
                    # Have the label slightly darker, so it doesn't blend in
                    # too much with the line for the well
                    rgbWell = clrs.to_rgba(wellColor[:])
                    rgbWell = np.array(list(rgbWell[:]))
                    rgbWell *= 0.925
                    # alpha is one
                    rgbWell[-1] = 1

                    if not wellLabel is None:
                        if '{}' in wellLabel:
                            plt3d_ax.text(
                                x_vals_well / 1000, y_vals_well / 1000, -z_min,
                                str(wellLabel.format(number)), zdir='x', color=rgbWell,
                                alpha=wellAlpha, fontsize=genfontsize - 4,
                                fontweight='bold')
                        else:
                            plt3d_ax.text(
                                x_vals_well / 1000, y_vals_well / 1000, -z_min,
                                wellLabel, zdir='x', color=rgbWell,
                                alpha=wellAlpha, fontsize=genfontsize - 4,
                                fontweight='bold')
                    else:
                        plt3d_ax.text(x_vals_well / 1000, y_vals_well / 1000, -z_min,
                                compName, zdir='x', color=rgbWell, alpha=wellAlpha,
                                fontsize=genfontsize - 4, fontweight='bold')

                if comp.class_type == 'OpenWellbore':
                    if z_max != 0:
                        plt3d_ax.plot([x_vals_well / 1000, x_vals_well / 1000],
                                      [y_vals_well / 1000, y_vals_well / 1000],
                                      [-z_max, 0], marker='o', markersize=3,
                                      color=wellColor, alpha=0.5, linewidth=1)

                    if plot_well_labels:
                        if z_max != 0:
                            plt3d_ax.text(
                                x_vals_well / 1000, y_vals_well / 1000, -z_max,
                                'Well\nTop', zdir='x', color=rgbWell,
                                alpha=wellAlpha, fontsize=genfontsize - 4,
                                fontweight='bold')
                        else:
                            plt3d_ax.text(
                                x_vals_well / 1000, y_vals_well / 1000, -z_max,
                                'Well\nTop', zdir='x', color=rgbWell,
                                alpha=wellAlpha, fontsize=genfontsize - 4,
                                fontweight='bold', zorder=200)

    if plot_indiv_strat_comps:
        for comp in components:
            if comp.class_type in wellbore_components:
                res_comp = sm.component_models[yaml_data[comp.name]['Connection']]
                x_vals_well = res_comp.locX
                y_vals_well = res_comp.locY

                if var_type in ['strikeAndDip', 'LookupTable']:
                    strata_temp = sm.component_models['strata' + comp.name]
                elif var_type == 'noVariation':
                    strata_temp = sm.component_models['strata']

                z_temp = iam_strata.get_unit_depth_from_component(
                    numShaleLayers, strata_temp, unitType='reservoir',
                    top_or_bottom='bottom')

                plt3d_ax.plot(x_vals_well / 1000, y_vals_well / 1000,
                              -z_temp, marker='s', markersize=2,
                              color=reservoirColor, linewidth=1)

                z_temp = iam_strata.get_unit_depth_from_component(
                    numShaleLayers, strata_temp, unitType='reservoir',
                    top_or_bottom='top')

                plt3d_ax.plot(x_vals_well / 1000, y_vals_well / 1000,
                              -z_temp, marker='s', markersize=2,
                              color=reservoirColor, linewidth=1)

                z_temp = iam_strata.get_unit_depth_from_component(
                    numShaleLayers, strata_temp, unitType='shale',
                    unitNumber=numShaleLayers, top_or_bottom='top')

                plt3d_ax.plot(x_vals_well / 1000, y_vals_well / 1000,
                              -z_temp, marker='s', markersize=2,
                              color=shaleColor[-1], linewidth=1)

                for shaleRef in range(numShaleLayers - 1):
                    z_temp = iam_strata.get_unit_depth_from_component(
                        numShaleLayers, strata_temp, unitType='shale',
                        unitNumber=(shaleRef + 1), top_or_bottom='top')

                    plt3d_ax.plot(x_vals_well / 1000, y_vals_well / 1000,
                                  -z_temp, marker='s', markersize=2,
                                  color=shaleColor[shaleRef], linewidth=1)

                    if (shaleRef + 1) < numShaleLayers:
                        z_temp = iam_strata.get_unit_depth_from_component(
                            numShaleLayers, strata_temp, unitType='aquifer',
                            unitNumber=(shaleRef + 1), top_or_bottom='top')

                        plt3d_ax.plot(x_vals_well / 1000, y_vals_well / 1000,
                                      -z_temp, marker='s', markersize=2,
                                      color=aquiferColor[shaleRef], linewidth=1)

    # I think 3D graphs are more clear if you add shadows beneath the features
    # (i.e., at the same x and y, but at the minimum z value). The shadows help
    # you judge each feature's position in the space. I also add vertical lines
    # at the edges of the surfaces.
    plt3d_ax.plot([x_loc[0, 0] / 1000, x_loc[0, 0] / 1000],
                 [y_loc[0, 0] / 1000, y_loc[0, 0] / 1000],
                 [zlim[0], 0], color='k', linewidth=1)

    plt3d_ax.plot([x_loc[-1, -1] / 1000, x_loc[-1, -1] / 1000],
                 [y_loc[0, 0] / 1000, y_loc[0, 0] / 1000],
                 [zlim[0], 0], color='k', linewidth=1)

    plt3d_ax.plot([x_loc[-1, -1] / 1000, x_loc[-1, -1] / 1000],
                 [y_loc[-1, -1] / 1000, y_loc[-1, -1] / 1000],
                 [zlim[0], 0], color='k', linewidth=1)

    plt3d_ax.plot([x_loc[0, 0] / 1000, x_loc[0, 0] / 1000],
                 [y_loc[-1, -1] / 1000, y_loc[-1, -1] / 1000],
                 [zlim[0], 0], color='k', linewidth=1)

    # Plot a 'shadow' beneath the planes
    if var_type == 'LookupTable':
        zmin_plane = np.ones((x_loc.shape[0], x_loc.shape[1])) * zlim[0]
        plt3d_ax.plot_surface(x_loc / 1000, y_loc / 1000, zmin_plane,
                              color=[0.5, 0.5, 0.5], alpha=0.25, shade=False)
    else:
        surface[:, :] = zlim[0]
        plt3d_ax.plot_surface(x_loc / 1000, y_loc / 1000, surface,
                              color=[0.5, 0.5, 0.5], alpha=0.25, shade=False)

    if zlim[0] >= -500:
        z_interval = 100
    elif -1000 <= zlim[0] < -500:
        z_interval = 250
    elif -3000 <= zlim[0] < -1000:
        z_interval = 500
    elif zlim[0] < -3000:
        z_interval = 1000

    closest_mult_above = (np.ceil(zlim[0] / -z_interval) * -z_interval) + z_interval
    ticks_pt2 = np.arange(closest_mult_above, 0, z_interval)
    ticks_pt2 = ticks_pt2.tolist()
    ticks_pt2.append(zlim[0])
    ticks_pt2.append(zlim[1])
    ticks_pt2.append(0)
    zticks = np.unique(ticks_pt2)

    zticks = zticks.tolist()
    plt3d_ax.set_zticks(zticks)

    if EnforceXandYLims:
        xticks = np.linspace(xLims[0] / 1000, xLims[1] / 1000, numberXTicks,
                             endpoint=True)
        xticks = xticks.tolist()
        plt3d_ax.set_xticks(xticks)

        yticks = np.linspace(yLims[0] / 1000, yLims[1] / 1000, numberYTicks,
                             endpoint=True)
        yticks = yticks.tolist()
        plt3d_ax.set_yticks(yticks)

    axislabel_pad = 10
    plt.xlabel('Easting (km)', fontsize=axislabelfontsize,
               fontweight='bold', labelpad=axislabel_pad)
    plt.ylabel('Northing (km)', fontsize=axislabelfontsize,
               fontweight='bold', labelpad=axislabel_pad)
    plt3d_ax.set_zlabel('Depth (m)', fontsize=axislabelfontsize,
                        fontweight='bold', labelpad=axislabel_pad)

    if not title:
        if var_type == 'LookupTable':
            plt3d_ax.set_title(
                'Stratigraphy for the Study Area,\nRed Squares: '
                + 'Top of Shales, Blue Squares: Top of Aquifers,\nGray '
                + 'Squares: Top and Bottom of the Reservoir',
                fontsize=titlefontsize, fontweight=selected_labelfontweight,
                pad=0)
        else:
            plt3d_ax.set_title(
                'Stratigraphy for the Study Area,\nTop of Each Unit Shown',
                fontsize=titlefontsize, fontweight=selected_labelfontweight,
                pad=0)
    else:
        plt3d_ax.set_title(title, fontweight=selected_labelfontweight,
                           fontsize=titlefontsize, pad=0)

    # Make the x and y axes have equal aspects
    xlim = plt3d_ax.get_xlim()
    ylim = plt3d_ax.get_ylim()

    x_range = xlim[1] - xlim[0]
    y_range = ylim[1] - ylim[0]

    if x_range > y_range:
        y_center = (ylim[1] + ylim[0]) / 2
        new_min_y = y_center - (x_range / 2)
        new_max_y = y_center + (x_range / 2)
        plt3d_ax.set_ylim(new_min_y, new_max_y)

    elif y_range > x_range:
        x_center = (xlim[1] + xlim[0]) / 2
        new_min_x = x_center - (y_range / 2)
        new_max_x = x_center + (y_range / 2)
        plt3d_ax.set_xlim(new_min_x, new_max_x)

    if plot_SandD_symbol:
        # Make the strike and dip symbol by using 5 known points. The point
        # locations are calculated with dx and dy values, which are calculated
        # with a length scale L. L is made to scale with the area, but it can
        # also be specified in the .yaml file.
        if Enforce_SandD_symbol_length:
            L = SandD_symbol_length / 2
        else:
            L = 2 * (((max_x_val - min_x_val) * (max_y_val - min_y_val)) ** 0.25)

        if var_type == 'strikeAndDip':

            x_points, y_points, z_points = iam_strata.make_xyz_points_from_strike_and_dip(
                dip, dipDirectionDegrees, L1=L * 0.75, L2=L, L3=L, L4=L,
                point0_xyz=[SandD_location[0], SandD_location[1], 0])

            x_p0 = x_points[0]
            x_p1 = x_points[1]
            x_p2 = x_points[2]
            x_p3 = x_points[3]
            x_p4 = x_points[4]

            y_p0 = y_points[0]
            y_p1 = y_points[1]
            y_p2 = y_points[2]
            y_p3 = y_points[3]
            y_p4 = y_points[4]

            z_p0 = z_points[0]
            z_p2 = z_points[2]
            z_p3 = z_points[3]

            plt3d_ax.plot([x_p2 / 1000, x_p3 / 1000],
                          [y_p2 / 1000, y_p3 / 1000],
                          [z_p2, z_p3],
                          color='k', linewidth=2, zorder=100)

            plt3d_ax.plot([x_p0 / 1000, x_p1 / 1000],
                          [y_p0 / 1000, y_p1 / 1000],
                          [z_p0, z_p0],
                          color='k', linewidth=2, zorder=101)

            plt3d_ax.text(x_p4 / 1000, y_p4 / 1000, z_p0,
                          str(dip), zdir='x', color='w',
                          fontsize=genfontsize - 2, fontweight='bold',
                          zorder=102)

        elif var_type == 'noVariation':
            degrees_for_circle = np.arange(0, 360, 1)

            circleX = np.zeros(len(degrees_for_circle))
            circleY = np.zeros(len(degrees_for_circle))

            for ind, d_val in enumerate(degrees_for_circle):
                circleX[ind] = SandD_location[0] + np.cos(np.radians(d_val)) * L
                circleY[ind] = SandD_location[1] + np.sin(np.radians(d_val)) * L

            plt3d_ax.plot(circleX / 1000, circleY / 1000, 0, marker=None,
                          color='k', linewidth=1, zorder=100)

            plt3d_ax.plot([circleX[0] / 1000, circleX[179] / 1000],
                          [circleY[0] / 1000, circleY[179] / 1000],
                          [0, 0], marker=None, color='k', linewidth=1, zorder=101)

            plt3d_ax.plot([circleX[89] / 1000, circleX[269] / 1000],
                          [circleY[89] / 1000, circleY[269] / 1000],
                          [0, 0], marker=None, color='k', linewidth=1, zorder=102)

    plt3d_ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    plt3d_ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

    plt.tight_layout()

    if save_stratigraphy:
        save_stratigraphy_by_loc_to_csv(model_data['OutputDirectory'],
                                        x_loc, y_loc, numShaleLayers,
                                        stratigraphy_by_loc)

    if savefig:
        for saveRef, elev_val in enumerate(view_elev):
            plt3d_ax.view_init(elev_val, view_azimuth[saveRef])

            if '.' in name:
                name_pt1 = name[0:name.index('.')]
                name_pt2 = name[name.index('.'):None]
                fig_dir_and_name = os.path.join(
                    model_data['OutputDirectory'],
                    name_pt1 + '_View' + str(saveRef + 1) + name_pt2)
            else:
                fig_dir_and_name = os.path.join(model_data['OutputDirectory'],
                                                name + '_View' + str(saveRef + 1))

            try:
                plt3d.savefig(fig_dir_and_name, bbox_inches='tight', dpi=figdpi)
            except ValueError:
                # User has specified plot with a '.' in name but no extension.
                # Add .png as output format.
                fig_dir_and_name += '.png'
                plt3d.savefig(fig_dir_and_name, bbox_inches='tight', dpi=figdpi)
    else:
        plt.show()


def not_boolean_debug_message(input_name, name, default_value):
    """
    Returns string delivering debug message regarding a variable not being
    of boolean type and setting it to the default value (True or False).

    input_name: string
    default_value: True or False
    """
    msg = ''.join(['Input for {} within the Stratigraphy plot {} was not of ',
                   'boolean type. Using the default value of {}.']).format(
                       input_name, name, default_value)
    return msg


def not_of_length_2_message(input_name, name):
    """
    Returns string delivering debug message regarding a list not being
    of length 2.

    input_name: string. Name of input not satisfying the right format.
    name: string. Name of plot setup incorrectly.
    """
    msg = ''.join([
        'The {} provided for the Stratigraphy plot {} is not a list ',
        'of length 2. Check your inputs in the .yaml file.']).format(input_name,
                                                                     name)
    return msg


def first_present_second_missing_message(input_name1, input_name2, name):
    """
    Returns string delivering debug message regarding expectations for user input
    to provide two of the relevant inputs when one is missing.

    input_name1: string. Name of input provided in the setup.
    input_name2: string. Name of input missing in the setup.
    name: string. Name of plot setup incorrectly.
    """
    extra_message = ''
    if input_name1 in ['ViewAngleElevation', 'ViewAngleAzimuth']:
        extra_message = 'Default values will be used.'
    elif input_name1 in ['InjectionCoordx', 'InjectionCoordy']:
        extra_message = 'Check your input. Injection sites will not be displayed.'

    msg = ''.join(['{} was provided for the Stratigraphy plot {} in the .yaml file',
                   'but {} was not. {}']).format(input_name1, input_name2,
                                                 name, extra_message)
    return msg


def not_between_message(input_name, low_bound, upper_bound, name,
                        units='degrees', value_is_list=True):
    """
    Returns string delivering debug message when the input value or values
    is not within the required boundaries.

    input_name: string. Name of input provided in the setup.
    low_bound: int or float. Lower boundary for the input value.
    upper_bound: int or float. Upper boundary for the input value.
    name: string. Name of plot setup incorrectly.
    units: string. Name of units for input value(s).
    value_is_list: bool. Falg indicating whether the message refers to a single
        value or list of values
    """
    if value_is_list:
        msg = ''.join(['At least one of the {} values provided for ',
                       'the Stratigraphy plot {} does not fall between {} ',
                       'and {} {}. Check the input in the .yaml file.']).format(
                           input_name, name, low_bound, upper_bound, units)
    else:
        msg = ''.join(['The {} value provided for the Stratigraphy plot {} ',
                       'does not fall between {} and {} {}. Check the input ',
                       'in the .yaml file.']).format(
                           input_name, name, low_bound, upper_bound, units)
    return msg


def read_strata_plot_yaml_input(yaml_data, name):
    """
    Function that reads the Stratigraphy plot input provided in the .yaml file
    and returns a dictionary containing the input.

    :param yaml_data: Dictionary of input values from the entire .yaml file
    :type yaml_data: dict

    :param name: Name of the Stratigraphy plot provided in the Plots section of
        the .yaml file
    :type name: str

    :returns: yaml_input
    """
    # Default values
    yaml_input_keys = [
        'dpi_input', 'EnforceXandYLims', 'xLims', 'yLims', 'EnforceXandYGridLims',
        'gridXLims', 'gridYLims', 'xGridSpacing', 'yGridSpacing',
        'plot_injection_sites','InjectionCoordx', 'InjectionCoordy',
        'plot_injection_site_labels', 'plot_wellbore_locations',
        'plot_well_labels', 'plot_indiv_strat_comps', 'save_stratigraphy',
        'plot_SandD_symbol', 'SandD_location', 'Enforce_SandD_symbol_length',
        'SandD_symbol_length', 'view_elev', 'view_azimuth']

    # Initialize output
    yaml_input = {key: None for key in yaml_input_keys}

    # Get shortcut to data to be analyzed
    strata_plot_data = yaml_data['Plots'][name]['Stratigraphy']

    if strata_plot_data is not None:
        # Checks for and adds colors, alpha values, and labels provided for
        # units. Excludes any invalid input provided.
        yaml_input = iam_strata.check_color_alpha_label_yaml_input(
            yaml_input, strata_plot_data, name)

        if 'FigureDPI' in strata_plot_data:
            yaml_input['dpi_input'] = strata_plot_data['FigureDPI']

        if 'SpecifyXandYLims' in strata_plot_data:
            yaml_input['EnforceXandYLims'] = True
            yaml_input['xLims'] = strata_plot_data['SpecifyXandYLims']['xLims']
            yaml_input['yLims'] = strata_plot_data['SpecifyXandYLims']['yLims']

            for input_key in ['xLims', 'yLims']:
                if len(yaml_input[input_key]) != 2:
                    debug_msg = not_of_length_2_message(input_key, name)
                    logging.debug(debug_msg)
                    yaml_input['EnforceXandYLims'] = False

        if 'SpecifyXandYGridLims' in strata_plot_data:
            yaml_input['EnforceXandYGridLims'] = True
            yaml_input['gridXLims'] = \
                strata_plot_data['SpecifyXandYGridLims']['gridXLims']
            yaml_input['gridYLims'] = \
                strata_plot_data['SpecifyXandYGridLims']['gridYLims']

            for input_key in ['gridXLims', 'gridYLims']:
                if len(yaml_input[input_key]) != 2:
                    debug_msg = not_of_length_2_message(input_key, name)
                    logging.debug(debug_msg)
                    yaml_input['EnforceXandYGridLims'] = False

        if 'xGridSpacing' in strata_plot_data:
            yaml_input['xGridSpacing'] = strata_plot_data['xGridSpacing']

        if 'yGridSpacing' in strata_plot_data:
            yaml_input['yGridSpacing'] = strata_plot_data['yGridSpacing']

        if 'PlotInjectionSites' in strata_plot_data:
            yaml_input['plot_injection_sites'] = strata_plot_data['PlotInjectionSites']
            if not isinstance(yaml_input['plot_injection_sites'], bool):
                debug_msg = not_boolean_debug_message(
                    'PlotInjectionSites', name, False)
                logging.debug(debug_msg)
                yaml_input['plot_injection_sites'] = False

        if 'InjectionCoordx' in strata_plot_data:
            yaml_input['InjectionCoordx'] = strata_plot_data['InjectionCoordx']

        if 'InjectionCoordy' in strata_plot_data:
            yaml_input['InjectionCoordy'] = strata_plot_data['InjectionCoordy']
            if yaml_input['InjectionCoordx'] is None:
                debug_msg = first_present_second_missing_message(
                    'InjectionCoordy', 'InjectionCoordx', name)
                logging.debug(debug_msg)
                yaml_input['plot_injection_sites'] = False

            elif len(yaml_input['InjectionCoordx']) != len(yaml_input['InjectionCoordy']):
                debug_msg = ''.join([
                    'The InjectionCoordy provided for Stratigraphy plot ', name,
                    'was not of the same length as the InjectionCoordx provided. ',
                    'Check your input. Injection sites will not be displayed.'])
                logging.debug(debug_msg)
                yaml_input['plot_injection_sites'] = False

        if yaml_input['InjectionCoordx'] is not None and \
                yaml_input['InjectionCoordy'] is None:
            debug_msg = first_present_second_missing_message(
                'InjectionCoordx', 'InjectionCoordy', name)
            logging.debug(debug_msg)
            yaml_input['plot_injection_sites'] = False

        if 'PlotInjectionSiteLabels' in strata_plot_data:
            yaml_input['plot_injection_site_labels'] = \
                strata_plot_data['PlotInjectionSiteLabels']
            if not isinstance(yaml_input['plot_injection_site_labels'], bool):
                debug_msg = not_boolean_debug_message(
                    'PlotInjectionSiteLabels', name, False)
                logging.debug(debug_msg)
                yaml_input['plot_injection_site_labels'] = False

        if 'PlotWellbores' in strata_plot_data:
            yaml_input['plot_wellbore_locations'] = strata_plot_data['PlotWellbores']
            if not isinstance(yaml_input['plot_wellbore_locations'], bool):
                debug_msg = not_boolean_debug_message('PlotWellbores', name, True)
                logging.debug(debug_msg)
                yaml_input['plot_wellbore_locations'] = True

        if 'PlotWellLabels' in strata_plot_data:
            yaml_input['plot_well_labels'] = strata_plot_data['PlotWellLabels']
            if not isinstance(yaml_input['plot_well_labels'], bool):
                debug_msg = not_boolean_debug_message('PlotWellLabels', name, True)
                logging.debug(debug_msg)
                yaml_input['plot_well_labels'] = True

        if 'PlotStratComponents' in strata_plot_data:
            yaml_input['plot_indiv_strat_comps'] = \
                strata_plot_data['PlotStratComponents']
            if not isinstance(yaml_input['plot_indiv_strat_comps'], bool):
                debug_msg = not_boolean_debug_message(
                    'PlotStratComponents', name, False)
                logging.debug(debug_msg)
                yaml_input['plot_indiv_strat_comps'] = False

        if 'SaveCSVFiles' in strata_plot_data:
            yaml_input['save_stratigraphy'] = strata_plot_data['SaveCSVFiles']
            if not isinstance(yaml_input['save_stratigraphy'], bool):
                debug_msg = not_boolean_debug_message('SaveCSVFiles', name, False)
                logging.debug(debug_msg)
                yaml_input['save_stratigraphy'] = False

        if 'StrikeAndDipSymbol' in strata_plot_data:
            if 'PlotSymbol' in strata_plot_data['StrikeAndDipSymbol']:
                yaml_input['plot_SandD_symbol'] = \
                    strata_plot_data['StrikeAndDipSymbol']['PlotSymbol']

            if not isinstance(yaml_input['plot_SandD_symbol'], bool):
                debug_msg = not_boolean_debug_message(
                    'StrikeAndDipSymbol: PlotSymbol', name, True)
                logging.debug(debug_msg)
                yaml_input['plot_SandD_symbol'] = True

            if 'length' in strata_plot_data['StrikeAndDipSymbol']:
                yaml_input['SandD_symbol_length'] = \
                    strata_plot_data['StrikeAndDipSymbol']['length']
                yaml_input['Enforce_SandD_symbol_length'] = True

            if 'coordx' in strata_plot_data['StrikeAndDipSymbol'] and \
                    'coordy' in strata_plot_data['StrikeAndDipSymbol'] and \
                        yaml_input['plot_SandD_symbol']:
                yaml_input['SandD_location'] = [
                    strata_plot_data['StrikeAndDipSymbol']['coordx'],
                    strata_plot_data['StrikeAndDipSymbol']['coordy']]

        if 'View' in strata_plot_data:
            if 'ViewAngleElevation' in strata_plot_data['View'] and \
                    'ViewAngleAzimuth' not in strata_plot_data['View']:
                debug_msg = first_present_second_missing_message(
                    'ViewAngleElevation', 'ViewAngleAzimuth', name)
                logging.debug(debug_msg)

            if 'ViewAngleElevation' not in strata_plot_data['View'] and \
                    'ViewAngleAzimuth' in strata_plot_data['View']:
                debug_msg = first_present_second_missing_message(
                    'ViewAngleAzimuth', 'ViewAngleElevation', name)
                logging.debug(debug_msg)

            if 'ViewAngleElevation' not in strata_plot_data['View'] and \
                    'ViewAngleAzimuth' not in strata_plot_data['View']:
                debug_msg = ''.join(['View was provided for the Stratigraphy ',
                                     'plot ', name, ', but ViewAngleElevation ',
                                     'and ViewAngleElevation were not included ',
                                     'within the View option. Using default values.'])
                logging.debug(debug_msg)

            if 'ViewAngleElevation' in strata_plot_data['View'] and \
                    'ViewAngleAzimuth' in strata_plot_data['View']:
                yaml_input['view_elev'] = strata_plot_data['View']['ViewAngleElevation']
                yaml_input['view_azimuth'] = strata_plot_data['View']['ViewAngleAzimuth']
                if len(yaml_input['view_elev']) != len(yaml_input['view_azimuth']):
                    err_msg = ''.join([
                        'The number of ViewAngleElevation values provided for ',
                        'the Stratigraphy plot ', name, ' does not match the ',
                        'number of ViewAngleAzimuth values provided. ',
                        'Check the input in the .yaml file.'])
                    raise ValueError(err_msg)

                if isinstance(yaml_input['view_elev'], list):
                    for view_elevs in yaml_input['view_elev']:
                        if not -90 < view_elevs < 90:
                            err_msg = not_between_message(
                                'ViewAngleElevation', -90, 90, name,
                                units='degrees', value_is_list=True)
                            raise ValueError(err_msg)

                elif not -90 < yaml_input['view_elev'] < 90:
                    err_msg = not_between_message(
                        'ViewAngleElevation', -90, 90, name,
                        units='degrees', value_is_list=False)
                    raise ValueError(err_msg)

                if isinstance(yaml_input['view_azimuth'], list):
                    for view_azimuths in yaml_input['view_azimuth']:
                        if not 0 <= view_azimuths <= 360:
                            err_msg = not_between_message(
                                'ViewAngleAzimuth', 0, 360, name,
                                units='degrees', value_is_list=True)
                            raise ValueError(err_msg)

                elif not 0 <= yaml_input['view_azimuth'] <= 360:
                    err_msg = not_between_message(
                        'ViewAngleAzimuth', 0, 360, name,
                        units='degrees', value_is_list=False)
                    raise ValueError(err_msg)

    return yaml_input


def check_SandD_location(SandD_location, min_x_val, max_x_val, min_y_val,
                         max_y_val, x_buffer, y_buffer):
    """
    Function that checks if the location provided for the strike and dip symbol
    (SandD_location) is within the grid used for the 3D planes. If the location
    is not within that grid, this function returns reset_SandD_location = True.

    :param SandD_location: tuple containing the x and y coordinates at which the
        strike and dip symbol will be plotted. The first value is for x, the
        second is for y. The default value is None.
    :type SandD_location: tuple or None

    :param min_x_val: Minimum x value (m) for the grid used to create 3D planes
    :type min_x_val: float or int

    :param max_x_val: Maximum x value (m) for the grid used to create 3D planes
    :type max_x_val: float or int

    :param min_y_val: Minimum y value (m) for the grid used to create 3D planes
    :type min_y_val: float or int

    :param max_x_val: Maximum y value (m) for the grid used to create 3D planes
    :type max_x_val: float or int

    :param x_buffer: Distance (m) used to adjust the min_x_val and max_x_val
        values. It is subtracted from min_x_val and added to max_x_val.
    :type x_buffer: float or int

    :param y_buffer: Distance (m) used to adjust the min_y_val and max_y_val
        values. It is subtracted from min_y_val and added to max_y_val.
    :type y_buffer: float or int

    :returns: reset_SandD_location
    """
    reset_SandD_location = False

    if SandD_location[0] < (min_x_val - x_buffer) or \
            SandD_location[0] > (max_x_val + x_buffer):
        warning_msg = ''.join(['The x value provided for the strike and dip ',
                               'symbol was beyond the x values used within ',
                               'the domain. Using the default location.'])
        logging.debug(warning_msg)
        reset_SandD_location = True

    if SandD_location[1] < (min_y_val - y_buffer) or \
            SandD_location[1] > (max_y_val + y_buffer):
        warning_msg = ''.join(['The y value provided for the strike and dip ',
                               'symbol was beyond the y values used within ',
                               'the domain. Using the default location.'])
        logging.debug(warning_msg)
        reset_SandD_location = True

    return reset_SandD_location


def create_strata_planes_dict(x_loc, y_loc, numShaleLayers,
                              shaleThicknessesReferencePoint,
                              aquiferThicknessesReferencePoint,
                              reservoirThicknessReferencePoint,
                              var_type, strataReferencePoint,
                              coordxReferencePoint=0, coordyReferencePoint=0,
                              strike=None, dip=None, dipDirection=None):
    """
    Function that creates a dictionary containing gridded values of unit
    thicknesses, top depths, and bottom depths.

    :param x_loc: 2-dimensional grid of x values (m) used to create and display
        the 3D planes
    :type x_loc: numpy.ndarray

    :param y_loc: 2-dimensional grid of y values (m) used to create and display
        the 3D planes
    :type y_loc: numpy.ndarray

    :param numShaleLayers: number of shale layers
    :type numShaleLayers: int or float

    :param shaleThicknessesReferencePoint: List of shale thicknesses, where the
        first entry is the lowest shale and the last is the highest. If strike
        and dip values are being used, these thicknesses only apply at
        the reference point.
    :type shaleThicknessesReferencePoint: list

    :param aquiferThicknessesReferencePoint: List of aquifer thicknesses, where
        the first entry is the lowest aquifer and the last is the highest. If
        strike and dip values are being used, these thicknesses only apply at
        the reference point.
    :type aquiferThicknessesReferencePoint: list

    :param reservoirThicknessReferencePoint: Thickness of the reservoir at the
        reference point.
    :type reservoirThicknessReferencePoint: list

    :param var_type: Option used for the domain's stratigraphy. 'noVariation'
        means the strata are flat and unchanging over space while
        'strikeAndDip' means that unit thicknesses change according to strike,
        dip, and dipDirection values. 'LookupTable' means that the stratigraphy
        at each component location (e.g., wells) is specified in a table - that
        option is not handled by this function.
    :type var_type: str

    :param strike: The units' strike in degrees clockwise from north, so that
        90 is east, 180 is south, and 270 is west. The default value is None.
    :type strike: int, float, or None

    :param dip: The units' dip in degrees.
    :type dip: int, float, or None

    :param dipDirection: The units' dip direction in a cardinal direction: N,
        E, S, W, NE, SE, SW, or NW.
    :type dip: str or None

    :returns: stratigraphy_by_loc
    """

    stratigraphy_by_loc = dict()

    stratigraphy_by_loc['resTopDepth'] = np.zeros((x_loc.shape[0], x_loc.shape[1]))
    stratigraphy_by_loc['resBottomDepth'] = np.zeros((x_loc.shape[0], x_loc.shape[1]))

    for shaleRef in range(numShaleLayers - 1, -1, -1):
        nm_v1 = 'shale{}'.format(shaleRef + 1)
        stratigraphy_by_loc[nm_v1 + 'Thickness'] = np.zeros((x_loc.shape[0],
                                                             x_loc.shape[1]))
        stratigraphy_by_loc[nm_v1 + 'TopDepth'] = np.zeros((x_loc.shape[0],
                                                            x_loc.shape[1]))
        if (shaleRef + 1) < numShaleLayers:
            nm_v2 = 'aquifer{}'.format(shaleRef + 1)
            stratigraphy_by_loc[nm_v2 + 'Thickness'] = np.zeros((x_loc.shape[0],
                                                                 x_loc.shape[1]))
            stratigraphy_by_loc[nm_v2 + 'TopDepth'] = np.zeros((x_loc.shape[0],
                                                                x_loc.shape[1]))

    if var_type == 'strikeAndDip':
        for x_ref in range(x_loc.shape[0]):
            for y_ref in range(x_loc.shape[1]):
                x_loc_temp = float(x_loc[x_ref, y_ref])
                y_loc_temp = float(y_loc[x_ref, y_ref])

                updatedStratigraphy = iam_strata.update_stratigraphy_by_strike_and_dip(
                    numberOfShaleLayers=numShaleLayers,
                    shaleThicknessList=shaleThicknessesReferencePoint[:],
                    aquiferThicknessList=aquiferThicknessesReferencePoint[:],
                    reservoirThickness=reservoirThicknessReferencePoint,
                    strike=strike, dip=dip, dipDirection=dipDirection,
                    coordxRefPoint=coordxReferencePoint,
                    coordyRefPoint=coordyReferencePoint,
                    location_x=x_loc_temp, location_y=y_loc_temp,
                    strataRefPoint=strataReferencePoint)

                shaleThicknessListUpdated = updatedStratigraphy['shaleThicknessList']
                aquiferThicknessListUpdated = updatedStratigraphy['aquiferThicknessList']

                shaleTopDepthListUpdated = updatedStratigraphy['shaleTopDepthList']
                aquiferTopDepthListUpdated = updatedStratigraphy['aquiferTopDepthList']
                resTopDepthUpdated = updatedStratigraphy['reservoirTopDepth']

                resBottomDepthUpdated = updatedStratigraphy['reservoirBottomDepth']

                stratigraphy_by_loc['resTopDepth'][x_ref, y_ref] = resTopDepthUpdated

                stratigraphy_by_loc['resBottomDepth'][x_ref, y_ref] = resBottomDepthUpdated

                for shaleRef in range(numShaleLayers - 1, -1, -1):
                    nm_v1 = 'shale{}'.format(shaleRef + 1)
                    stratigraphy_by_loc[nm_v1+'Thickness'][x_ref, y_ref] = \
                        shaleThicknessListUpdated[shaleRef]

                    stratigraphy_by_loc[nm_v1+'TopDepth'][x_ref, y_ref] = \
                        shaleTopDepthListUpdated[shaleRef]

                    if (shaleRef + 1) < numShaleLayers:
                        nm_v2 = 'aquifer{}'.format(shaleRef + 1)
                        stratigraphy_by_loc[nm_v2 + 'Thickness'][x_ref, y_ref] = \
                            aquiferThicknessListUpdated[shaleRef]

                        stratigraphy_by_loc[nm_v2 + 'TopDepth'][x_ref, y_ref] = \
                            aquiferTopDepthListUpdated[shaleRef]

    elif var_type == 'noVariation':
        stratigraphy_by_loc['resThickness'] = np.ones(
            (x_loc.shape[0], x_loc.shape[1])) * reservoirThicknessReferencePoint

        for shaleRef in range(numShaleLayers - 1, -1, -1):
            stratigraphy_by_loc['shale{}Thickness'.format(shaleRef + 1)] = np.ones(
                (x_loc.shape[0], x_loc.shape[1])) * shaleThicknessesReferencePoint[shaleRef]

            if (shaleRef + 1) < numShaleLayers:
                stratigraphy_by_loc['aquifer{}Thickness'.format(shaleRef + 1)] = np.ones(
                    (x_loc.shape[0], x_loc.shape[1])) * aquiferThicknessesReferencePoint[shaleRef]

        stratigraphy_by_loc['shale{}TopDepth'.format(numShaleLayers)][:, :] = 0

        for shaleRef in range(numShaleLayers - 1, 0, -1):
            stratigraphy_by_loc['aquifer{}TopDepth'.format(shaleRef)][:, :] = \
                stratigraphy_by_loc['shale{}TopDepth'.format(shaleRef+1)][:, :] \
                                        + stratigraphy_by_loc['shale' + str(shaleRef + 1)
                                                        + 'Thickness'][:, :]

            stratigraphy_by_loc['shale{}TopDepth'.format(shaleRef)][:, :] = \
                stratigraphy_by_loc['aquifer{}TopDepth'.format(shaleRef)][:, :]  \
                    + stratigraphy_by_loc['aquifer{}Thickness'.format(shaleRef)][:, :]

        stratigraphy_by_loc['resTopDepth'][:, :] = \
            stratigraphy_by_loc['shale1TopDepth'][:, :] \
                + stratigraphy_by_loc['shale1Thickness'][:, :]

        stratigraphy_by_loc['resBottomDepth'][:, :] = \
            stratigraphy_by_loc['resTopDepth'][:, :] \
                + stratigraphy_by_loc['resThickness'][:, :]

    return stratigraphy_by_loc


def save_stratigraphy_by_loc_to_csv(output_dir, x_loc, y_loc, numShaleLayers,
                                    stratigraphy_by_loc):
    """
    Function that saves the stratigraphy information in stratigraphy_by_loc in
    .csv files.

    :param output_dir: Directory in which to save the .csv files
    :type output_dir: str

    :param x_loc: 2-dimensional grid of x values (m) used to create and display
        the 3D planes
    :type x_loc: numpy.ndarray

    :param y_loc: 2-dimensional grid of y values (m) used to create and display
        the 3D planes
    :type y_loc: numpy.ndarray

    :param numShaleLayers: number of shale layers
    :type numShaleLayers: int or float

    :param stratigraphy_by_loc: Dictionary containing the 2D- distribution of
        unit thicknesses, top depths, and bottom depths. This dictionary is
        produced by the function create_strata_planes_dict().
    :type stratigraphy_by_loc: dict

    :returns: None
    """
    if not os.path.exists(os.path.join(output_dir, 'csv_files')):
        os.mkdir(os.path.join(output_dir, 'csv_files'))

    strata_filename = os.path.join(output_dir, 'csv_files', 'Strata_x_values_m.csv')
    with open(strata_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        for row_ref in range(x_loc.shape[0]):
            writer.writerow(x_loc[row_ref, :])

    strata_filename = os.path.join(output_dir, 'csv_files', 'Strata_y_values_m.csv')
    with open(strata_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        for row_ref in range(y_loc.shape[0]):
            writer.writerow(y_loc[row_ref, :])

    for shaleRef in range(1, numShaleLayers + 1):
        thickness_par_name = 'shale{}Thickness'.format(shaleRef)
        shaleThicknesses_temp = stratigraphy_by_loc[thickness_par_name][:, :]

        strata_filename = os.path.join(
            output_dir, 'csv_files', thickness_par_name+'.csv')
        with open(strata_filename, 'w', newline='') as f:
            writer = csv.writer(f)
            for row_ref in range(x_loc.shape[0]):
                writer.writerow(shaleThicknesses_temp[row_ref, :])

        depth_par_name ='shale{}TopDepth'.format(shaleRef)
        shaleTopDepth_temp = stratigraphy_by_loc[depth_par_name][:, :]

        strata_filename = os.path.join(
            output_dir, 'csv_files', depth_par_name+'.csv')
        with open(strata_filename, 'w', newline='') as f:
            writer = csv.writer(f)
            for row_ref in range(x_loc.shape[0]):
                writer.writerow(shaleTopDepth_temp[row_ref, :])

        if shaleRef < numShaleLayers:
            thickness_par_name = 'aquifer{}Thickness'.format(shaleRef)
            aquiferThicknesses_temp = stratigraphy_by_loc[thickness_par_name][:, :]

            strata_filename = os.path.join(
                output_dir, 'csv_files', thickness_par_name+'.csv')
            with open(strata_filename, 'w', newline='') as f:
                writer = csv.writer(f)
                for row_ref in range(x_loc.shape[0]):
                    writer.writerow(aquiferThicknesses_temp[row_ref, :])

            depth_par_name ='aquifer{}TopDepth'.format(shaleRef)
            aquiferTopDepth_temp = stratigraphy_by_loc[depth_par_name][:, :]

            strata_filename = os.path.join(
                output_dir, 'csv_files', depth_par_name+'.csv')
            with open(strata_filename, 'w', newline='') as f:
                writer = csv.writer(f)
                for row_ref in range(x_loc.shape[0]):
                    writer.writerow(aquiferTopDepth_temp[row_ref, :])
