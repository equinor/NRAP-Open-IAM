'''
This example demonstrates the use of the time to first detection (TTFD) plot
type. Using ttfd_plot() from ttfd.py can actually create a large variety of
results. It creates maps showing the spatiotemporal evolution of contaminant
plumes, the TTFD for those plumes given the location of sensors in the monitoring
network, and a map showing the probability of a plume's presence at each location
within the domain. This last plot is only created if the simulation uses Latin
Hypercube Sampling (lhs). The probability is calculated as the number of
realizataions in which a plume occurred at the location divided by the total
number of realizations (set by the num_samples parameter here).

If the MonitoringLocations section is present, as shown here, then separate
graphs and .csv files will be made with the TTFD at the monitor locations.
The HorizontalWindow and VerticalWindow are the horizontal and vertical distances
(in meters) used to search for plume timings from the monitoring locations. These
search windows are used because TTFD results are assessed within grids, and the
x, y, and z values at grid point will not be exact matches of the monitoring
locations' x, y, and z values. If a monitoring location is close enough to a
grid point (i.e., within the HorizontalWindow and VerticalWindow), then the
results at that grid point will be included for the monitoring location. If not
entered, the HorizontalWindow and VerticalWindow will both resort to the default
value of 1 m. The VerticalWindow can be considered a length
along a monitoring well with sensor equipment, with the sensors covering the
distance from coordz - VerticalWindow to coordz + VerticalWindow.

Note that one can control the x- and y-spacings of grid points with the
xGridSpacing and yGridSpacing entries, which should be in meters. The minimum
and maximum x and y values can be set with the SpecifyXandYGridLims, gridXLims,
and gridLims entires. The z values are set by unit depths, but the z-spacing
can be altered with the NumZPointsWithinAquifers and NumZPointsWithinShales
entries. These parameters are the number of grid points extending from the
bottom to the top of each aquifer and shale, respectively. These values must
be entered as integers. The x, y, and z grid values used for the TTFD plot
results are saved in .csv files.

Example of run:
$ python iam_sys_reservoir_mswell_fgen_ttfd.py
'''

import sys
import os
import numpy as np
import random
import datetime
sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, Stratigraphy, SimpleReservoir,
                     MultisegmentedWellbore, FutureGen2AZMI, FutureGen2Aquifer,
                     RateToMassAdapter)

import openiam as iam
import openiam.openiam_cf_strata as iam_strata
import openiam.visualize as iam_vis


if __name__ == "__main__":
    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25 * np.arange(0.0, num_years + 1)
    sm_model_kwargs = {'time_array': time_array}   # time is given in days

    # These are used for setting up the LHS simulation
    num_samples = 30
    ncpus = 4

    # Stratigraphy information
    datumPressure = 101325
    numberOfShaleLayers = 3
    shale1Thickness = 650
    shale2Thickness = 450
    shale3Thickness = 150
    aquifer1Thickness = 50
    aquifer2Thickness = 50
    reservoirThickness = 75

    aquifer2Depth = shale3Thickness + aquifer2Thickness
    aquifer1Depth = aquifer2Depth + shale2Thickness + aquifer1Thickness

    # These lists are required by stratigraphy_plot()
    shaleThicknessList = [shale1Thickness, shale2Thickness, shale3Thickness]
    aquiferThicknessList = [aquifer1Thickness, aquifer2Thickness]

    injectionX = 1000
    injectionY = 1000

    max_x_value = 2000
    max_y_value = 2000

    # Set up the dictionaries required by stratigraphy_plot()
    model_data = dict()
    model_data['OutputDirectory'] = os.path.join(
        iam.IAM_DIR, 'output', 'ttfd_script_example_'
        + str(datetime.date.today()))
    model_data['Analysis'] = dict()
    model_data['Analysis']['type'] = 'lhs'
    model_data['Analysis']['siz'] = num_samples

    yaml_data = dict()

    yaml_data['Stratigraphy'] = dict()
    yaml_data['Stratigraphy']['shale1Thickness'] = shale1Thickness
    yaml_data['Stratigraphy']['shale2Thickness'] = shale2Thickness
    yaml_data['Stratigraphy']['shale2Thickness'] = shale3Thickness
    yaml_data['Stratigraphy']['aquifer1Thickness'] = aquifer1Thickness
    yaml_data['Stratigraphy']['aquifer2Thickness'] = aquifer2Thickness
    yaml_data['Stratigraphy']['reservoirThickness'] = reservoirThickness

    plotName = 'Plot1.png'
    yaml_data['Plots'] = dict()
    yaml_data['Plots'][plotName] = dict()
    yaml_data['Plots'][plotName]['TTFD'] = dict()

    # Required entries within the TTFD plot
    # This entry specifies the type of plume. The options are Pressure, pH, TDS,
    # Dissolved_CO2, Dissolved_salt, and CarbonateAquifer. These different
    # plume types correspond with the output names of plume dimension metrics.
    # For example, the FutureGen2Aquifer produces TDS_dx, TDS_dy, and TDS_dz
    # while the GenericAquifer produces Dissolved_salt_dr and Dissolved_salt_dz.
    # The CarbonateAquifer only produces dx and dy, so TTFD plots using Carbonate
    # Aquifer output need plume_type = 'CarbonateAquifer.' This example is set
    # up to handle 'Pressure', 'pH', or 'TDS' plumes from thr FutureGen2Aquifer
    # and FutureGen2AZMI components.
    yaml_data['Plots'][plotName]['TTFD']['plume_type'] = 'pH'
    plume_type = yaml_data['Plots'][plotName]['TTFD']['plume_type']
    # If you added more aquifer components and wanted to use their output in
    # the TTFD plot, you would have to add the component names here. In the
    # code below, I create entries in yaml_data with the component names -
    # these entries are required for the TTFD plot.
    yaml_data['Plots'][plotName]['TTFD']['aquifer_name_list'] = []

    use_specific_settings = True
    # I do this to demonstrate that the remaining options are not required. Try
    # setting use_specific_settings to False.
    if use_specific_settings:
        yaml_data['Plots'][plotName]['TTFD']['MonitoringLocations'] = dict()
        # There are 8 sensors located at 4 (x,y) locations. At each location,
        # there are two sensor: one in the middle of aquifer 1, and one in the
        # middle of aquifer 2.
        yaml_data['Plots'][plotName]['TTFD']['MonitoringLocations'][
            'coordx'] = [900, 900, 1100, 1100, 900, 900, 1100, 1100]
        yaml_data['Plots'][plotName]['TTFD']['MonitoringLocations'][
            'coordy'] = [900, 1100, 1100, 900, 900, 1100, 1100, 900]
        # Note that these coordz values must be negative.
        aquifer1_coordz = -aquifer1Depth + (aquifer1Thickness / 2)
        aquifer2_coordz = -aquifer2Depth + (aquifer2Thickness / 2)
        yaml_data['Plots'][plotName]['TTFD']['MonitoringLocations'][
            'coordz'] = [aquifer1_coordz, aquifer1_coordz,
                         aquifer1_coordz, aquifer1_coordz,
                         aquifer2_coordz, aquifer2_coordz,
                         aquifer2_coordz, aquifer2_coordz]

        # This controls how close the monitor must be to a grid location for it
        # to 'detect' a plume at that location
        yaml_data['Plots'][plotName]['TTFD']['MonitoringLocations'][
            'HorizontalWindow'] = 1
        # Because the sensors are in the middle of each aquifer and VerticaWindow
        # is (aquifer2Thickness / 4), it will detect any plumes in the central
        # 50% of the aquifers. Even so, there will likely be cases where a plume
        # is missed!
        yaml_data['Plots'][plotName]['TTFD']['MonitoringLocations'][
            'VerticalWindow'] = (aquifer2Thickness / 4)
        yaml_data['Plots'][plotName]['TTFD']['PlotInjectionSites'] = True
        yaml_data['Plots'][plotName]['TTFD']['xGridSpacing'] = 10
        yaml_data['Plots'][plotName]['TTFD']['yGridSpacing'] = 10
        yaml_data['Plots'][plotName]['TTFD']['NumZPointsWithinAquifers'] = 10
        yaml_data['Plots'][plotName]['TTFD']['NumZPointsWithinShales'] = 3
        yaml_data['Plots'][plotName]['TTFD']['WriteDreamOutput'] = False
        yaml_data['Plots'][plotName]['TTFD']['SaveCSVFiles'] = True

        yaml_data['Plots'][plotName]['TTFD']['SpecifyXandYLims'] = dict()
        yaml_data['Plots'][plotName]['TTFD']['SpecifyXandYLims'][
            'xLims'] = [0, max_x_value]
        yaml_data['Plots'][plotName]['TTFD']['SpecifyXandYLims'][
            'yLims'] = [0, max_y_value]

        yaml_data['Plots'][plotName]['TTFD']['SpecifyXandYGridLims'] = dict()
        yaml_data['Plots'][plotName]['TTFD']['SpecifyXandYGridLims'][
            'gridXLims'] = [100, max_x_value - 100]
        yaml_data['Plots'][plotName]['TTFD']['SpecifyXandYGridLims'][
            'gridYLims'] = [100, max_y_value - 100]

    # Create system model
    sm = SystemModel(model_kwargs = sm_model_kwargs)

    # Lists of components
    strata = []
    sres = []
    ms = []
    adapts = []
    aquifer1Comps = []
    aquifer2Comps = []

    # Add stratigraphy component
    strata.append(sm.add_component_model_object(Stratigraphy(
        name = 'strata', parent = sm)))

    # Add parameters for the reference point stratrigraphy component
    strata[-1].add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)
    strata[-1].add_par('shale1Thickness', value=shale1Thickness, vary=False)
    strata[-1].add_par('shale2Thickness', value=shale2Thickness, vary=False)
    strata[-1].add_par('shale3Thickness', value=shale3Thickness, vary=False)
    strata[-1].add_par('aquifer1Thickness', value=aquifer1Thickness, vary=False)
    strata[-1].add_par('aquifer2Thickness', value=aquifer2Thickness, vary=False)
    strata[-1].add_par('reservoirThickness', value=reservoirThickness, vary=False)
    strata[-1].add_par('datumPressure', value=datumPressure, vary=False)

    numberOfWells = 4
    well_x_values = [900, 900, 1100, 1100]
    well_y_values = [900, 1100, 1100, 900]

    for locRef in range(0, numberOfWells):
        # The names for the reservoir and wellbore components at this location.
        # This naming convention is required by ttfd_plot() because the
        # convention is used by the control file interface. If the locRef values
        # exceeded 9, you would have to remove one of the 0's so there would
        # still be a total of 3 digits.
        sresName = 'SimpleReservoir_00' + str(locRef)
        msName = 'MultisegmentedWellbore_00' + str(locRef)
        adapterName = 'Adapter_00' + str(locRef)
        # There is a name for aquifer 1 and 2 in case the unit thicknesses are
        # changed enough to require a different component type.
        FGen2Aq1Name = 'FutureGen2Aquifer1_00' + str(locRef)
        FGen2AZMI1Name = 'FutureGen2AZMI1_00' + str(locRef)
        FGen2Aq2Name = 'FutureGen2Aquifer2_00' + str(locRef)
        FGen2AZMI2Name = 'FutureGen2AZMI2_00' + str(locRef)

        # Add reservoir component for the current well
        sres.append(sm.add_component_model_object(SimpleReservoir(
            name = sresName, parent = sm, locX = well_x_values[locRef],
            locY = well_y_values[locRef], injX = injectionX, injY = injectionY)))

        # Add parameters of reservoir component model
        sres[-1].add_par('injRate', value=0.5, vary=False)
        sres[-1].add_par('brineDensity', value=1100, vary=False)
        sres[-1].add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)
        sres[-1].add_par('shale1Thickness', value=shaleThicknessList[0], vary=False)
        sres[-1].add_par('shale2Thickness', value=shaleThicknessList[1], vary=False)
        sres[-1].add_par('shale3Thickness', value=shaleThicknessList[2], vary=False)
        sres[-1].add_par('aquifer1Thickness', value=aquiferThicknessList[0], vary=False)
        sres[-1].add_par('aquifer2Thickness', value=aquiferThicknessList[1], vary=False)
        sres[-1].add_par('reservoirThickness', value=reservoirThickness, vary=False)
        sres[-1].add_par('datumPressure', value=datumPressure, vary=False)

        sres[-1].add_obs('pressure')
        sres[-1].add_obs('CO2saturation')
        sres[-1].add_obs_to_be_linked('pressure')
        sres[-1].add_obs_to_be_linked('CO2saturation')

        # Add multisegmented wellbore component
        ms.append(sm.add_component_model_object(
            MultisegmentedWellbore(name=msName, parent=sm)))

        # Add parameters of multisegmented wellbore component model
        ms[-1].add_par('wellRadius', value=0.1, vary=False)
        ms[-1].add_par('logWellPerm', value=-12.0, min=-12.25, max=-11.75)

        # Add linked parameters: common to reservoir and wellbore components
        ms[-1].add_par_linked_to_par(
            'numberOfShaleLayers', sres[-1].deterministic_pars['numberOfShaleLayers'])
        ms[-1].add_par_linked_to_par(
            'shale1Thickness', sres[-1].deterministic_pars['shale1Thickness'])
        ms[-1].add_par_linked_to_par(
            'shale2Thickness', sres[-1].deterministic_pars['shale2Thickness'])
        ms[-1].add_par_linked_to_par(
            'shale3Thickness', sres[-1].deterministic_pars['shale3Thickness'])
        ms[-1].add_par_linked_to_par(
            'aquifer1Thickness', sres[-1].deterministic_pars['aquifer1Thickness'])
        ms[-1].add_par_linked_to_par(
            'aquifer2Thickness', sres[-1].deterministic_pars['aquifer2Thickness'])
        ms[-1].add_par_linked_to_par(
            'reservoirThickness', sres[-1].deterministic_pars['reservoirThickness'])
        ms[-1].add_par_linked_to_par(
            'datumPressure', sres[-1].deterministic_pars['datumPressure'])

        # Add keyword arguments linked to the output provided by reservoir model
        ms[-1].add_kwarg_linked_to_obs('pressure', sres[-1].linkobs['pressure'])
        ms[-1].add_kwarg_linked_to_obs('CO2saturation', sres[-1].linkobs['CO2saturation'])

        # Add observations of multisegmented wellbore component model
        ms[-1].add_obs('brine_aquifer1')
        ms[-1].add_obs('CO2_aquifer1')
        ms[-1].add_obs('brine_aquifer2')
        ms[-1].add_obs('CO2_aquifer2')
        ms[-1].add_obs('brine_atm')
        ms[-1].add_obs('CO2_atm')
        ms[-1].add_obs_to_be_linked('brine_aquifer1')
        ms[-1].add_obs_to_be_linked('CO2_aquifer1')
        ms[-1].add_obs_to_be_linked('brine_aquifer2')
        ms[-1].add_obs_to_be_linked('CO2_aquifer2')
        ms[-1].add_obs_to_be_linked('brine_atm')
        ms[-1].add_obs_to_be_linked('CO2_atm')

        # Add adapter that transforms leakage rates to accumulated mass
        adapts.append(sm.add_component_model_object(
            RateToMassAdapter(name=adapterName, parent=sm)))
        adapts[-1].add_kwarg_linked_to_collection(
            'CO2_aquifer1', [ms[-1].linkobs['CO2_aquifer1'],
                             ms[-1].linkobs['CO2_aquifer2']])
        adapts[-1].add_kwarg_linked_to_collection(
            'CO2_aquifer2', [ms[-1].linkobs['CO2_aquifer2'],
                             ms[-1].linkobs['CO2_atm']])
        adapts[-1].add_kwarg_linked_to_collection(
            'brine_aquifer1', [ms[-1].linkobs['brine_aquifer1'],
                              ms[-1].linkobs['brine_aquifer2']])
        adapts[-1].add_kwarg_linked_to_collection(
            'brine_aquifer2', [ms[-1].linkobs['brine_aquifer2'],
                              ms[-1].linkobs['brine_atm']])
        adapts[-1].add_obs_to_be_linked('mass_CO2_aquifer1')
        adapts[-1].add_obs_to_be_linked('mass_brine_aquifer1')
        adapts[-1].add_obs_to_be_linked('mass_CO2_aquifer2')
        adapts[-1].add_obs_to_be_linked('mass_brine_aquifer2')
        adapts[-1].add_obs('mass_CO2_aquifer1')
        adapts[-1].add_obs('mass_brine_aquifer1')
        adapts[-1].add_obs('mass_CO2_aquifer2')
        adapts[-1].add_obs('mass_brine_aquifer2')

        # Add aquifer compopnent for aquifer 1. FutureGen2AZMI applies to depths
        # > 700 m. This is set up to change with unit thickness values.
        if aquifer1Depth < 700:
            aquifer1Comps.append(sm.add_component_model_object(
                FutureGen2AZMI(name=FGen2AZMI1Name, parent=sm)))

        else:
            aquifer1Comps.append(sm.add_component_model_object(
                FutureGen2Aquifer(name=FGen2Aq1Name, parent=sm)))

        # Add parameters of FutureGen2AZMI or FutureGen2Aquifer component model
        aquifer1Comps[-1].add_par('aqu_thick', value=aquifer1Thickness, vary=False)
        aquifer1Comps[-1].add_par('depth', value=-12.0, vary=False)
        aquifer1Comps[-1].add_par('por', value=0.12, min=0.09, max=0.15)
        aquifer1Comps[-1].add_par('log_permh', value=-13.39, min=-13.89,
                            max=-12.89)
        aquifer1Comps[-1].add_par('rel_vol_frac_calcite', value=0.025,
                                  min=0.01, max=0.05)

        # Add aquifer component's keyword argument co2_rate linked to the
        # collection created above
        aquifer1Comps[-1].add_kwarg_linked_to_obs('co2_rate',
                                         ms[-1].linkobs['CO2_aquifer1'])
        aquifer1Comps[-1].add_kwarg_linked_to_obs('brine_rate',
                                         ms[-1].linkobs['brine_aquifer1'])
        aquifer1Comps[-1].add_kwarg_linked_to_obs('co2_mass',
                                         adapts[-1].linkobs['mass_CO2_aquifer1'])
        aquifer1Comps[-1].add_kwarg_linked_to_obs('brine_mass',
                                         adapts[-1].linkobs['mass_brine_aquifer1'])

        # Add observations (output) from the futuregen aquifer model
        aquifer1Comps[-1].add_obs(plume_type + '_dx')
        aquifer1Comps[-1].add_obs(plume_type + '_dy')
        aquifer1Comps[-1].add_obs(plume_type + '_dz')

        # Add the component for aquifer 2
        if aquifer2Depth < 700:
            aquifer2Comps.append(sm.add_component_model_object(
                FutureGen2AZMI(name = FGen2AZMI2Name, parent = sm)))

        else:
            aquifer2Comps.append(sm.add_component_model_object(
                FutureGen2Aquifer(name = FGen2Aq2Name, parent = sm)))

        # Add parameters of FutureGen2AZMI or FutureGen2Aquifer component model
        aquifer2Comps[-1].add_par('aqu_thick', value=aquifer1Thickness, vary=False)
        aquifer2Comps[-1].add_par('depth', value=-12.0, vary=False)
        aquifer2Comps[-1].add_par('por', value=0.12, min=0.09, max=0.15)
        aquifer2Comps[-1].add_par('log_permh', value=-13.39, min=-13.89,
                            max=-12.89)
        aquifer2Comps[-1].add_par('rel_vol_frac_calcite', value=0.025,
                                  min=0.01, max=0.05)

        # Add aquifer component's keyword argument co2_rate linked to the
        # collection created above
        aquifer2Comps[-1].add_kwarg_linked_to_obs('co2_rate',
                                         ms[-1].linkobs['CO2_aquifer2'])
        aquifer2Comps[-1].add_kwarg_linked_to_obs('brine_rate',
                                         ms[-1].linkobs['brine_aquifer2'])
        aquifer2Comps[-1].add_kwarg_linked_to_obs('co2_mass',
                                         adapts[-1].linkobs['mass_CO2_aquifer2'])
        aquifer2Comps[-1].add_kwarg_linked_to_obs('brine_mass',
                                         adapts[-1].linkobs['mass_brine_aquifer2'])

        # Add observations (output) from the futuregen aquifer model
        aquifer2Comps[-1].add_obs(plume_type + '_dx')
        aquifer2Comps[-1].add_obs(plume_type + '_dy')
        aquifer2Comps[-1].add_obs(plume_type + '_dz')

        # These yaml_data entries are needed for the TTFD plot. To see how this
        # approach is similar to the control files, compare these entries with
        # those in ControlFile_ex39.yaml.
        yaml_data['strata'] = dict()
        yaml_data[sres[-1].name] = dict()
        yaml_data[sres[-1].name]['Type'] = 'SimpleReservoir'

        yaml_data[ms[-1].name] = dict()
        yaml_data[ms[-1].name]['Type'] = 'MultisegmentedWellbore'
        yaml_data[ms[-1].name]['Connection'] = sres[-1].name

        yaml_data[adapts[-1].name] = dict()

        yaml_data[aquifer1Comps[-1].name] = dict()
        if aquifer1Depth > 700:
            yaml_data[aquifer1Comps[-1].name]['Type'] = 'FutureGen2AZMI'
        else:
            yaml_data[aquifer1Comps[-1].name]['Type'] = 'FutureGen2Aquifer'
        yaml_data[aquifer1Comps[-1].name]['Connection'] = ms[-1].name
        yaml_data[aquifer1Comps[-1].name]['AquiferName'] = 'aquifer1'

        yaml_data[aquifer2Comps[-1].name] = dict()
        if aquifer2Depth > 700:
            yaml_data[aquifer2Comps[-1].name]['Type'] = 'FutureGen2AZMI'
        else:
            yaml_data[aquifer2Comps[-1].name]['Type'] = 'FutureGen2Aquifer'
        yaml_data[aquifer2Comps[-1].name]['Connection'] = ms[-1].name
        yaml_data[aquifer2Comps[-1].name]['AquiferName'] = 'aquifer2'

        yaml_data['Plots'][plotName]['TTFD']['aquifer_name_list'].append(
            aquifer1Comps[-1].name)
        yaml_data['Plots'][plotName]['TTFD']['aquifer_name_list'].append(
            aquifer2Comps[-1].name)

    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=random.randint(500, 1100))

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)

    if not os.path.exists(model_data['OutputDirectory']):
        os.mkdir(model_data['OutputDirectory'])

    savefig = os.path.join(model_data['OutputDirectory'], plotName)

    plots = yaml_data['Plots']

    iam_vis.ttfd_plot(
        yaml_data, model_data, sm, s, model_data['OutputDirectory'],
        name=plotName, analysis='lhs')
