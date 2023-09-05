'''
This example demonstrates the use of the time to the first detection (TTFD) plot
type through a script (instead of a control file). Additionally, this example
uses Latin Hypercube Sampling (lhs) to vary the 'logWellPerm' and 'logAqu1Perm'
parameters of the MultisegmentedWellbore components.

To use the TTFD plot type in a script, one has to set up the yaml_data dictionary
to match the format used by a control file. For an example, see Control File
example 39a.

The TTFD plot type has more requirements than most other plot types. For example,
the aquifer components corresponding with the 'ComponentNameList' entry in
yaml_data need to each produce the plume metrics required - here, 'pH_dx',
'pH_dy', and 'pH_dz'. The plume metrics required correspond with the 'PlumeType'
entry provided (here, 'pH').

Example of run:
$ python iam_sys_reservoir_mswell_futuregen_ttfdplot_no_dip_lhs.py
'''

import sys
import os
import random
import datetime
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, Stratigraphy, SimpleReservoir,
                     MultisegmentedWellbore, RateToMassAdapter, FutureGen2AZMI)

import openiam as iam
import openiam.visualize as iam_vis


if __name__ == "__main__":
    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25 * np.arange(0.0, num_years + 1)
    sm_model_kwargs = {'time_array': time_array}   # time is given in days

    # Stratigraphy information
    datumPressure = 101325
    numberOfShaleLayers = 3
    shale1Thickness = 750
    shale2Thickness = 950
    shale3Thickness = 365
    aquifer1Thickness = 85
    aquifer2Thickness = 85
    reservoirThickness = 150

    num_aquifers = 2

    # These lists are required by stratigraphy_plot()
    shaleThicknesses = [shale1Thickness, shale2Thickness, shale3Thickness]
    aquiferThicknesses = [aquifer1Thickness, aquifer2Thickness]

    aquiferTopDepths = [shaleThicknesses[1] + shaleThicknesses[2] + aquifer2Thickness,
                        shaleThicknesses[2]]

    injectionX = 0
    injectionY = 0

    # Set up the dictionaries required by stratigraphy_plot()
    model_data = dict()
    model_data['OutputDirectory'] = os.path.join(
        iam.IAM_DIR, 'output', 'ttfdplot_example_no_dip_lhs_'
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
    yaml_data['Plots'][plotName]['TTFD'] = dict()

    # These entries are required
    yaml_data['Plots'][plotName]['TTFD']['PlumeType'] = 'pH'
    yaml_data['Plots'][plotName]['TTFD']['ComponentNameList'] = ['FutureGen2AZMI1']

    use_specific_settings = True
    # I do this to demonstrate that the remaining options are not required. Try
    # setting use_specific_settings to False.
    if use_specific_settings:
        yaml_data['Plots'][plotName]['TTFD']['PlotInjectionSites'] = True
        yaml_data['Plots'][plotName]['TTFD']['WriteDreamOutput'] = False

        monitor_coordx = [50, 75, 150]
        monitor_coordy = [50, 75, 150]

        # Near the top of aquifer 1
        monitor_coordz = [-aquiferTopDepths[0] - 5, -aquiferTopDepths[0] - 5,
                          -aquiferTopDepths[0] - 5]

        yaml_data['Plots'][plotName]['TTFD']['MonitoringLocations'] = {
            'coordx': monitor_coordx, 'coordy': monitor_coordy, 'coordz': monitor_coordz,
            'HorizontalWindow': 1, 'VerticalWindow': 5}

        yaml_data['Plots'][plotName]['TTFD']['NumZPointsWithinAquifers'] = 10
        yaml_data['Plots'][plotName]['TTFD']['NumZPointsWithinShales'] = 3

        yaml_data['Plots'][plotName]['TTFD']['xGridSpacing'] = 10
        yaml_data['Plots'][plotName]['TTFD']['yGridSpacing'] = 10

        yaml_data['Plots'][plotName]['TTFD']['FigureDPI'] = 200

        yaml_data['Plots'][plotName]['TTFD']['SaveCSVFiles'] = True

        yaml_data['Plots'][plotName]['TTFD']['SpecifyXandYLims'] = \
            {'xLims': [-75, 675],
             'yLims': [-75, 675]}

        yaml_data['Plots'][plotName]['TTFD']['SpecifyXandYGridLims'] = \
            {'gridXLims': [0, 650], 'gridYLims': [0, 650]}

    # Create system model
    sm = SystemModel(model_kwargs = sm_model_kwargs)

    # Lists of components
    strata = []
    sres = []
    ms = []
    adapt = []
    fgaq = []

    # List of all components except for stratigraphy components. This list is
    # needed for the stratigraphy_plot() function.
    components = []

    # Add stratigraphy component for the reference point
    strata.append(sm.add_component_model_object(Stratigraphy(
        name='strata', parent=sm)))

    # Add parameters for the reference point stratrigraphy component
    strata[-1].add_par('numberOfShaleLayers',
                       value=numberOfShaleLayers, vary=False)
    strata[-1].add_par('shale1Thickness',
                       value=shale1Thickness, vary=False)
    strata[-1].add_par('shale2Thickness',
                       value=shale2Thickness, vary=False)
    strata[-1].add_par('shale3Thickness',
                       value=shale3Thickness, vary=False)
    strata[-1].add_par('aquifer1Thickness',
                       value=aquifer1Thickness, vary=False)
    strata[-1].add_par('aquifer2Thickness',
                       value=aquifer2Thickness, vary=False)
    strata[-1].add_par('reservoirThickness',
                       value=reservoirThickness, vary=False)
    strata[-1].add_par('datumPressure', value=datumPressure, vary=False)

    well_x_values = [100, 200, 300, 400, 500]
    well_y_values = [100, 200, 300, 400, 500]
    numberOfWells = len(well_x_values)

    for locRef in range(0, numberOfWells):
        # The names for the reservoir and wellbore components at this location.
        # This naming convention is required by stratigraphy_plot() because the
        # convention is used by the control file interface.
        sresName = 'SimpleReservoir1_{0:03}'.format(locRef)
        msName = 'MultisegmentedWellbore1_{0:03}'.format(locRef)
        fgaqName = 'FutureGen2AZMI1_{0:03}'.format(locRef)
        adaptName = 'Adapter_{0:03}'.format(locRef)

        # Add reservoir component for the current well
        sres.append(sm.add_component_model_object(
            SimpleReservoir(name=sresName, parent=sm,
                            locX=well_x_values[locRef], locY=well_y_values[locRef],
                            injX=injectionX, injY=injectionY)))

        # Add parameters of reservoir component model
        sres[-1].add_par('injRate', value=10, vary=False)
        sres[-1].add_par('logResPerm', value=-12, vary=False)
        sres[-1].add_par('brineDensity', value=1030.9, vary=False)
        sres[-1].add_par('CO2Density', value=775.0, vary=False)
        sres[-1].add_par('brineViscosity', value=7.5e-4, vary=False)
        sres[-1].add_par('CO2Viscosity', value=6.6e-5, vary=False)
        sres[-1].add_par('numberOfShaleLayers',
                         value=numberOfShaleLayers, vary=False)
        sres[-1].add_par('shale1Thickness',
                         value=shaleThicknesses[0], vary=False)
        sres[-1].add_par('shale2Thickness',
                         value=shaleThicknesses[1], vary=False)
        sres[-1].add_par('shale3Thickness',
                         value=shaleThicknesses[2], vary=False)
        sres[-1].add_par('aquifer1Thickness',
                         value=aquiferThicknesses[0], vary=False)
        sres[-1].add_par('aquifer2Thickness',
                         value=aquiferThicknesses[1], vary=False)
        sres[-1].add_par('reservoirThickness',
                         value=reservoirThickness, vary=False)
        sres[-1].add_par('datumPressure',
                         value=datumPressure, vary=False)

        sres[-1].add_obs('pressure')
        sres[-1].add_obs('CO2saturation')
        sres[-1].add_obs_to_be_linked('pressure')
        sres[-1].add_obs_to_be_linked('CO2saturation')

        # Add multisegmented wellbore component
        ms.append(sm.add_component_model_object(
            MultisegmentedWellbore(name=msName, parent=sm)))

        # Add parameters of multisegmented wellbore component model
        ms[-1].add_par('wellRadius', value=0.1, vary=False)
        ms[-1].add_par('logWellPerm', value=-13.0, min=-14.0, max=-11.0)
        ms[-1].add_par('logAqu1Perm', value=-15.39, min=-15.89, max=-14.89)
        ms[-1].add_par('logAqu2Perm', value=-13.05, vary=False)
        ms[-1].add_par('brineDensity', value=1030.9, vary=False)
        ms[-1].add_par('CO2Density', value=775, vary=False)
        ms[-1].add_par('brineViscosity', value=7.5e-4, vary=False)
        ms[-1].add_par('CO2Viscosity', value=6.6e-5, vary=False)

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

        adapt.append(sm.add_component_model_object(
            RateToMassAdapter(name=adaptName, parent=sm)))

        # CO2 leakage rates for mass adapter
        for num in range(1, num_aquifers):
            obs_nm1 = 'CO2_aquifer{}'.format(num)
            obs_nm2 = 'CO2_aquifer{}'.format(num+1)
            adapt[-1].add_kwarg_linked_to_collection(
                obs_nm1, [ms[-1].linkobs[obs_nm1], ms[-1].linkobs[obs_nm2]])

        obs_nm = 'CO2_aquifer{}'.format(num_aquifers)
        adapt[-1].add_kwarg_linked_to_collection(
            obs_nm, [ms[-1].linkobs[obs_nm], ms[-1].linkobs['CO2_atm']])

        # Brine leakage rates for mass adapter
        for num in range(1, num_aquifers):
            obs_nm1 = 'brine_aquifer{}'.format(num)
            obs_nm2 = 'brine_aquifer{}'.format(num+1)
            adapt[-1].add_kwarg_linked_to_collection(
                obs_nm1, [ms[-1].linkobs[obs_nm1], ms[-1].linkobs[obs_nm2]])

        obs_nm = 'brine_aquifer{}'.format(num_aquifers)
        adapt[-1].add_kwarg_linked_to_collection(
            obs_nm, [ms[-1].linkobs[obs_nm], ms[-1].linkobs['brine_atm']])

        for num in range(1, num_aquifers+1):
            adapt[-1].add_obs_to_be_linked('mass_CO2_aquifer{}'.format(num))
            adapt[-1].add_obs_to_be_linked('mass_brine_aquifer{}'.format(num))
            adapt[-1].add_obs('mass_CO2_aquifer{}'.format(num))
            adapt[-1].add_obs('mass_brine_aquifer{}'.format(num))

        # Add FutureGen2AZMI component
        fgaq.append(sm.add_component_model_object(
            FutureGen2AZMI(name=fgaqName, parent=sm)))

        fgaq[-1].add_par('aqu_thick', value=aquiferThicknesses[0], vary=False)
        fgaq[-1].add_par('depth', value=aquiferTopDepths[0], vary=False)
        fgaq[-1].add_par('por', value=0.18, vary=False)
        fgaq[-1].add_par('log_permh', value=-13.0, vary=False)
        fgaq[-1].add_par('log_aniso', value=0.3, vary=False)
        fgaq[-1].add_par('rel_vol_frac_calcite', value=0.1, vary=False)

        fgaq[-1].add_kwarg_linked_to_obs('co2_rate', ms[-1].linkobs['CO2_aquifer1'])
        fgaq[-1].add_kwarg_linked_to_obs('brine_rate', ms[-1].linkobs['brine_aquifer1'])

        fgaq[-1].add_kwarg_linked_to_obs('co2_mass', adapt[-1].linkobs['mass_CO2_aquifer1'])
        fgaq[-1].add_kwarg_linked_to_obs('brine_mass', adapt[-1].linkobs['mass_brine_aquifer1'])

        # If you change the PlumeType entry, change these outputs to match
        fgaq[-1].add_obs('pH_volume')
        fgaq[-1].add_obs('pH_dx')
        fgaq[-1].add_obs('pH_dy')
        fgaq[-1].add_obs('pH_dz')

        components.append(sres[-1])
        components.append(ms[-1])
        components.append(fgaq[-1])

        yaml_data[sres[-1].name] = dict()
        yaml_data[sres[-1].name]['Type'] = 'SimpleReservoir'

        yaml_data[ms[-1].name] = dict()
        yaml_data[ms[-1].name]['Type'] = 'MultisegmentedWellbore'
        yaml_data[ms[-1].name]['LeakTo'] = 'aquifer1'
        yaml_data[ms[-1].name]['Connection'] = sres[-1].name

        yaml_data[adapt[-1].name] = dict()
        yaml_data[adapt[-1].name]['Type'] = 'RateToMassAdapter'

        yaml_data[fgaq[-1].name] = dict()
        yaml_data[fgaq[-1].name]['Type'] = 'FutureGen2AZMI'
        yaml_data[fgaq[-1].name]['AquiferName'] = 'aquifer1'
        yaml_data[fgaq[-1].name]['Connection'] = ms[-1].name

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

    analysis = 'lhs'
    num_samples = 30
    ncpus = 1

    model_data['Analysis'] = {'type': analysis, 'siz': num_samples}

    analysis_dict = model_data['Analysis']
    analysis = analysis_dict.pop('type').lower()
    analysis_dict['seed'] = random.randint(500, 1100)

    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(**analysis_dict)

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose='progress', logfile=None)

    iam_vis.ttfd_plot(yaml_data, model_data, sm, s, model_data['OutputDirectory'],
                      name=plotName, analysis=analysis)
