'''
NRAP IAM AOR FutureGen Area of Review (AoR)
This examples couples the lookup table reservoir, strata, open wellbore and
FutureGen aquifer components.  The saturation/pressure output produced by lookup table
reservoir model is used to drive leakage from a single open wellbore
model, which is passed to the input of an adapter that provides well
coordinates, CO2 and brine leakage rates and cumulative mass fluxes to the
futuregen aquifer model.  The strata component defines the thickness of the aquifer
and shale layers above the reservoir, which are then used as input parameters
for the open wellbore and aquifer components. The risk-based Area of Review (AoR) is defined
based on the probability of whether there is an impact to groundwater at each grid location.

This example requires the additional FutureGen 2.0 data set.
FutureGen 2.0 data set can be downloaded from the following source:
https://edx.netl.doe.gov/dataset/futuregen-2-0-1008-simulation-reservoir-lookup-table

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/FutureGen2/1008_sims

Usage examples:
$ python iam_sys_lutreservoir_openwell_futuregen_aor.py --run 1
'''
# @author: Diana Bacon
# diana.bacon@pnnl.gov

import os
import sys
import argparse
import datetime
import numpy as np
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
            'Check this example description for more information.'])
        print(msg)

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
    aqu_thick = [33.2, 84.1, 31.1, 61.6]
    por = [0.118, 0.038, 0.132, 0.18]
    log_permh = [-13.39, -11.05, -12.48, -11.92]
    log_aniso = [0.30, 1.00, 0.30, 0.30]
    rel_vol_frac_calcite = [0.1, 0.5, 0.1, 0.1]
    num_shales = 5
    shale = [198.7, 74.4, 110.3, 118.9, 530.4]

    # Get grid locations from selected data file
    locations = np.genfromtxt(
        os.path.join(file_directory, 'fg{}.csv'.format(run)),
        delimiter=',', skip_header=1, dtype=None, usecols=(0, 1, 2))

    x_loc = locations[:, 0]
    y_loc = locations[:, 1]
    z_loc = locations[:, 2]

    print('Number of locations', len(x_loc))

    xmin = min(x_loc)
    xmax = max(x_loc)
    ymin = min(y_loc)
    ymax = max(y_loc)
    zmin = min(z_loc)
    zmax = max(z_loc)

    output_directory = os.sep.join(['..', '..', 'Output', 'AoR_{date_time_stamp}'.format(
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

    def max_impact(i):
        print('Starting {}'.format(i))
        x = x_loc[i]
        y = y_loc[i]
        z = z_loc[i]

        if x <= xmin or x >= xmax or y <= ymin or y >= ymax:
            return [0, 0, 0, 0, 0]

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
        ind = run-1
        signature = {par_names[j]: float(signature_data[ind+1, j+1]) for j in range(num_pars)}

        sm.add_interpolator(ReservoirDataInterpolator(
            name='int{}'.format(ind+1), parent=sm,
            header_file_dir=file_directory,
            time_file='time_points.csv',
            data_file='fg{}.csv'.format(ind+1),
            index=int(signature_data[ind+1, 0]),
            signature=signature), intr_family='reservoir')

        # Add reservoir component
        ltres = sm.add_component_model_object(LookupTableReservoir(
            name='ltres', parent=sm, intr_family='reservoir', locX=x, locY=y))

        # Add parameters of reservoir component model
        for j in range(num_pars):
            # add arbitrary line of values from signature_file
            ltres.add_par(par_names[j], value=float(signature_data[run, j+1]), vary=False)

        # Add observations of reservoir component model
        ltres.add_obs('pressure')
        ltres.add_obs('CO2saturation')
        ltres.add_obs_to_be_linked('pressure')
        ltres.add_obs_to_be_linked('CO2saturation')

        # Add stratigraphy component
        strata = sm.add_component_model_object(Stratigraphy(name='strata', parent=sm))

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

        # Add open wellbore component
        ow = sm.add_component_model_object(OpenWellbore(name='ow', parent=sm))

        # Add parameters of open wellbore component
        ow.add_par('wellRadius', value=0.005, vary=False)
        ow.add_par('logReservoirTransmissivity', min=-11.0, max=-9.0, value=-10.0)
        ow.add_par('logAquiferTransmissivity', min=-11.0, max=-9.0, value=-10.0)
        ow.add_par('brineSalinity', value=0.0475, vary=False)

        # Add keyword arguments of the open wellbore component model
        ow.add_kwarg_linked_to_obs('pressure', ltres.linkobs['pressure'])
        ow.add_kwarg_linked_to_obs('CO2saturation', ltres.linkobs['CO2saturation'])

        # Add composite parameter of open wellbore component
        ow.add_composite_par(
            'reservoirDepth',
            expr='+'.join(['strata.shale1Thickness', 'strata.shale2Thickness',
                           'strata.shale3Thickness', 'strata.shale4Thickness',
                           'strata.shale5Thickness', 'strata.aquifer1Thickness',
                           'strata.aquifer2Thickness', 'strata.aquifer3Thickness',
                           'strata.aquifer4Thickness']))
        ow.add_composite_par(
            'wellTop', expr='strata.shale5Thickness + strata.aquifer4Thickness')

        # Add observations of open wellbore component model
        ow.add_obs_to_be_linked('CO2_aquifer')
        ow.add_obs_to_be_linked('brine_aquifer')
        ow.add_obs_to_be_linked('brine_atm')
        ow.add_obs_to_be_linked('CO2_atm')
        ow.add_obs('CO2_aquifer')
        ow.add_obs('brine_aquifer')
        ow.add_obs('CO2_atm') # zero since well top is in aquifer
        ow.add_obs('brine_atm') # zero since well top is in aquifer

        # Add adapter that transforms leakage rates to accumulated mass
        adapt = sm.add_component_model_object(RateToMassAdapter(name='adapt', parent=sm))
        adapt.add_kwarg_linked_to_collection('CO2_aquifer',
                                             [ow.linkobs['CO2_aquifer'],
                                              ow.linkobs['CO2_atm']])
        adapt.add_kwarg_linked_to_collection('brine_aquifer',
                                             [ow.linkobs['brine_aquifer'],
                                              ow.linkobs['brine_atm']])
        adapt.add_obs_to_be_linked('mass_CO2_aquifer')
        adapt.add_obs_to_be_linked('mass_brine_aquifer')
        adapt.add_obs('mass_CO2_aquifer')
        adapt.add_obs('mass_brine_aquifer')

        # Add futuregen aquifer model object and define parameters
        if float(depth[3]) < 700:
            fga = sm.add_component_model_object(FutureGen2Aquifer(name='fga', parent=sm))
        else:
            fga = sm.add_component_model_object(FutureGen2AZMI(name='fga', parent=sm))

        # St. Peter Sandstone
        fga.add_par('aqu_thick', value=aqu_thick[3])
        fga.add_par('depth', value=depth[3])
        fga.add_par('por', value=0.18)
        fga.add_par('log_permh', value=-11.92)
        fga.add_par('log_aniso', value=0.30)
        fga.add_par('rel_vol_frac_calcite', value=0.01)

        # Add aquifer component's keyword argument co2_rate linked to the collection created above
        fga.add_kwarg_linked_to_obs('co2_rate', ow.linkobs['CO2_aquifer'])
        fga.add_kwarg_linked_to_obs('brine_rate', ow.linkobs['brine_aquifer'])
        fga.add_kwarg_linked_to_obs('co2_mass', adapt.linkobs['mass_CO2_aquifer'])
        fga.add_kwarg_linked_to_obs('brine_mass', adapt.linkobs['mass_brine_aquifer'])

        # Add observations (output) from the futuregen aquifer model
        fga.add_obs('pH_volume')
        fga.add_obs('TDS_volume')

        # Run system model using current values of its parameters
        sm.forward()  # system model is run deterministically

        # Get results
        pressure = sm.collect_observations_as_time_series(ltres, 'pressure')/1.0e6 # convert to MPa
        pdiff = [p - pressure[0] for p in pressure]
        sat = sm.collect_observations_as_time_series(ltres, 'CO2saturation')
        pH = sm.collect_observations_as_time_series(fga, 'pH_volume')
        TDS = sm.collect_observations_as_time_series(fga, 'TDS_volume')

        # Find maximum value for all time steps
        return [max(pressure), max(pdiff), max(sat), max(pH), max(TDS)]

    from sys import platform
    if platform == "win32":
        results = []
        # Loop replaces parallel execution of the simulations on Windows
        for ind in range(len(x_loc)):
            results.append(max_impact(ind))
    else:
        # The following code should work on Mac but not on Windows.
        from multiprocessing import Pool
        with Pool(processes=4) as pool:
            results = pool.map(max_impact, list(range(len(x_loc))))

    results = np.array(results)

    fig, ax = plt.subplots()
    plt.tricontourf(x_loc/1000.0, y_loc/1000.0, results[:, 3],
                    locator=ticker.MaxNLocator(nbins=100, prune='lower'), cmap=plt.cm.Greens)
    plt.colorbar(format='%.1e')
    levels = np.arange(0, .56, .08)
    CS = plt.tricontour(x_loc/1000.0, y_loc/1000.0, results[:, 2], levels=levels, colors='k')
    plt.xlim(233137/1000.0, 240379/1000.0)
    plt.ylim(4406487/1000.0, 4413729/1000.0)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
    ax = plt.gca()
    ax.set_aspect('equal')
    fmt = ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    fmt1 = ticker.FuncFormatter(lambda x, p: "%.2f" % x)
    clabels = plt.clabel(CS, inline=1, fontsize=8, fmt=fmt1)
    for txt in clabels:
        txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    plt.xlabel('Easting (km)')
    plt.ylabel('Northing (km)')
    plt.gca().get_xaxis().set_major_formatter(fmt)
    plt.gca().get_yaxis().set_major_formatter(fmt)
    plt.title('Maximum pH Volume, m$^3$\nGas Saturation (contours)')
    plt.tight_layout()
    to_save = True
    if to_save:
        plt.savefig(os.sep.join([output_directory, "AoR_pH{}.png".format(args.run)]), dpi=300)
        plt.savefig(os.sep.join([output_directory, "AoR_pH{}.pdf".format(args.run)]), dpi=300)

    fig, ax = plt.subplots()
    plt.tricontourf(x_loc/1000.0, y_loc/1000.0, results[:, 4],
                    locator=ticker.MaxNLocator(nbins=1000, prune='lower'), cmap=plt.cm.Blues)
    plt.colorbar(format='%.1e')
    levels = np.array([0.1, 0.2, 0.5, 1.0, 2.0])
    CS = plt.tricontour(x_loc/1000.0, y_loc/1000.0, results[:, 1], levels=levels, colors='k')
    plt.xlim(200548/1000.0, 272968/1000.0)
    plt.ylim(4373897/1000.0, 4446318/1000.0)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
    ax = plt.gca()
    ax.set_aspect('equal')
    fmt = ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    fmt1 = ticker.FuncFormatter(lambda x, p: "%.2f" % x)
    clabels = plt.clabel(CS, inline=1, fontsize=8, fmt=fmt1)
    for txt in clabels:
        txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0))
    plt.xlabel('Easting (km)')
    plt.ylabel('Northing (km)')
    plt.gca().get_xaxis().set_major_formatter(fmt)
    plt.gca().get_yaxis().set_major_formatter(fmt)
    plt.title('Maximum TDS Volume, m$^3$\nMaximum Pressure Increase, MPa (contours)')
    plt.tight_layout()
    to_save = True
    if to_save:
        plt.savefig(os.sep.join([output_directory, "AoR_TDS{}.png".format(args.run)]), dpi=300)
