'''
This example demonstrates the creation of dipping stratigraphy through the
use of the functions in cfi/strata.py and stratigraphy_plot.py. The
stratigraphy_plot function requires the input dictionaries yaml_data and
model_data. These dictionaries can be made with .yaml Control Files, but they
can also be made manually as demonstrated here. Note that the function requires
yaml_data['Plots'][plotName]['Stratigraphy'], where plotName is an arbitrary
user defined string (e.g., 'Plot1' or 'StratPlot'). The function does not
require any of the plotting options that can be contained within
yaml_data['Plots'][plotName]['Stratigraphy']. If any options are not present,
the default settings will be used.

This example creates a stratigraphy plot. Unlike the script example
iam_sys_reservoir_mswell_stratplot_dipping_strata.py, however, this example
does not use spatially variable stratigraphy (i.e., no dip).

Example of run:
$ python iam_sys_reservoir_mswell_stratplot_no_dip.py
'''

import sys
import os
import numpy as np
import random
import datetime
sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, Stratigraphy, SimpleReservoir,
                     MultisegmentedWellbore)

import openiam as iam
import openiam.cfi.strata as iam_strata
import openiam.visualize as iam_vis


if __name__ == "__main__":
    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25 * np.arange(0.0, num_years + 1)
    sm_model_kwargs = {'time_array': time_array}   # time is given in days

    # Stratigraphy information
    datumPressure = 101325
    numberOfShaleLayers = 3
    shale1Thickness = 650
    shale2Thickness = 450
    shale3Thickness = 150
    aquifer1Thickness = 200
    aquifer2Thickness = 100
    reservoirThickness = 75

    # These lists are required by stratigraphy_plot()
    shaleThicknessList = [shale1Thickness, shale2Thickness, shale3Thickness]
    aquiferThicknessList = [aquifer1Thickness, aquifer2Thickness]

    strike = 315
    dip = 5
    dipDirection = 'NE'
    coordxRefPoint = 0
    coordyRefPoint = 0

    dipDirectionDegrees = iam_strata.obtain_dip_direction_degrees(
        strike, dipDirection)

    injectionX = 2500
    injectionY = 2500

    max_x_value = 5000
    max_y_value = 5000

    # Set up the dictionaries required by stratigraphy_plot()
    model_data = dict()
    model_data['OutputDirectory'] = os.path.join(
        iam.IAM_DIR, 'output', 'stratplot_example_no_dip_'
        + str(datetime.date.today()))

    yaml_data = dict()
    yaml_data['Stratigraphy'] = {
        'shale1Thickness': shale1Thickness,
        'shale2Thickness': shale2Thickness,
        'shale3Thickness': shale3Thickness,
        'aquifer1Thickness': aquifer1Thickness,
        'aquifer2Thickness': aquifer2Thickness,
        'reservoirThickness': reservoirThickness}

    plotName = 'Plot1.png'
    yaml_data['Plots'] = dict()
    yaml_data['Plots'][plotName] = dict()
    yaml_data['Plots'][plotName]['Stratigraphy'] = dict()

    use_specific_settings = True
    # I do this to demonstrate that the remaining options are not required. Try
    # setting use_specific_settings to False.
    if use_specific_settings:
        for key in ['PlotWellbores', 'PlotWellLabels', 'PlotStratComponents',
                    'PlotInjectionSites', 'PlotInjectionSiteLabels']:
            yaml_data['Plots'][plotName]['Stratigraphy'][key] = True

        yaml_data['Plots'][plotName]['Stratigraphy']['FigureDPI'] = 100

        yaml_data['Plots'][plotName]['Stratigraphy']['SaveCSVFiles'] = False

        yaml_data['Plots'][plotName]['Stratigraphy']['SpecifyXandYLims'] = \
            {'xLims': [-100, max_x_value + 100],
             'yLims': [-100, max_y_value + 100]}

        yaml_data['Plots'][plotName]['Stratigraphy']['SpecifyXandYGridLims'] = \
            {'gridXLims': [0, max_x_value], 'gridYLims': [0, max_y_value]}

        yaml_data['Plots'][plotName]['Stratigraphy']['View'] = \
            {'ViewAngleElevation': [10, 10, 10],
             'ViewAngleAzimuth': [300, 315, 330]}

        yaml_data['Plots'][plotName]['Stratigraphy']['StrikeAndDipSymbol'] = \
            {'PlotSymbol': True, 'coordx': 1000, 'coordy': 1000, 'length': 750}

    # Create system model
    sm = SystemModel(model_kwargs = sm_model_kwargs)

    # Lists of components
    strata = []
    sres = []
    ms = []

    # List of all components except for stratigraphy components. This list is
    # needed for the stratigraphy_plot() function.
    components = []

    # Add stratigraphy component
    strata.append(sm.add_component_model_object(Stratigraphy(
        name='strata', parent=sm)))

    # Add parameters for the reference point stratrigraphy component
    strata[-1].add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)
    strata[-1].add_par('shale1Thickness', value=shale1Thickness, vary=False)
    strata[-1].add_par('shale2Thickness', value=shale2Thickness, vary=False)
    strata[-1].add_par('shale3Thickness', value=shale3Thickness, vary=False)
    strata[-1].add_par('aquifer1Thickness', value=aquifer1Thickness, vary=False)
    strata[-1].add_par('aquifer2Thickness', value=aquifer2Thickness, vary=False)
    strata[-1].add_par('reservoirThickness', value=reservoirThickness, vary=False)
    strata[-1].add_par('datumPressure', value=datumPressure, vary=False)
    # This is meant for the Control File interface, but it also calculates
    # parameters like unit depths. Those unit depths are used in stratigraphy_plot().
    strata[-1].connect_with_system()

    numberOfWells = 6
    np.random.seed(random.randint(0, 1.0e6))
    well_x_values = np.random.rand(1, numberOfWells) * max_x_value
    well_x_values = well_x_values.tolist()[0]
    well_y_values = np.random.rand(1, numberOfWells) * max_y_value
    well_y_values = well_y_values.tolist()[0]

    for locRef in range(0, numberOfWells):
        # The names for the reservoir and wellbore components at this location.
        # This naming convention is required by stratigraphy_plot() because the
        # convention is used by the control file interface.
        sresName = 'SimpleReservoir_{0:03}'.format(locRef)
        msName = 'MultisegmentedWellbore_{0:03}'.format(locRef)

        # Add reservoir component for the current well
        sres.append(sm.add_component_model_object(SimpleReservoir(
            name=sresName, parent=sm, locX=well_x_values[locRef],
            locY=well_y_values[locRef], injX=injectionX, injY=injectionY)))

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
        ms[-1].add_par('logWellPerm', value=-12.0, vary=False)

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

        components.append(sres[-1])
        components.append(ms[-1])

        yaml_data[sres[-1].name] = dict()
        yaml_data[sres[-1].name]['Type'] = 'SimpleReservoir'
        yaml_data[ms[-1].name] = dict()
        yaml_data[ms[-1].name]['Type'] = 'MultisegmentedWellbore'
        yaml_data[ms[-1].name]['Connection'] = sres[-1].name

    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    print('pressure, Well 1\n',
          sm.collect_observations_as_time_series(sres[0], 'pressure'))
    print('------------------------------------------------------------------')
    print('pressure, Well 2\n',
          sm.collect_observations_as_time_series(sres[1], 'pressure'))
    print('------------------------------------------------------------------')
    print('CO2_aquifer1, Well 1\n',
          sm.collect_observations_as_time_series(ms[0], 'CO2_aquifer1'))
    print('------------------------------------------------------------------')
    print('brine_aquifer1, Well 2\n',
          sm.collect_observations_as_time_series(ms[1], 'brine_aquifer1'))
    print('------------------------------------------------------------------')

    if not os.path.exists(model_data['OutputDirectory']):
        os.mkdir(model_data['OutputDirectory'])

    savefig = os.path.join(model_data['OutputDirectory'], plotName)

    iam_vis.stratigraphy_plot(yaml_data, model_data, sm,
                              name=plotName, savefig=savefig)
