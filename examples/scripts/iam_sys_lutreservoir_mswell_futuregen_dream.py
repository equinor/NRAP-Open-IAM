'''
NRAP IAM DREAM output example
This examples couples the FutureGen2 lookup table reservoir, multisegmented wellbore and
FutureGen2 aquifer and AZMI models. The saturation/pressure output produced by the
reservoir model is used to drive leakage from two multisegmented wellbore
models, which are passed to the input of an adapter that provides well
coordinates, CO2 and brine leakage rates and cumulative mass fluxes to three
FutureGen2 AZMI and one FutureGen2 aquifer models. A matrix of times to first detection
in the aquifer for each realization are output.

This example requires the additional FutureGen 2.0 data set.
FutureGen 2.0 data set can be downloaded from the following source:
https://edx.netl.doe.gov/dataset/futuregen-2-0-1008-simulation-reservoir-lookup-table

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/FutureGen2/1008_sims

Usage examples:
$ python iam_sys_lutreservoir_mswell_futuregen_dream.py --run 1 --metric Pressure
$ python iam_sys_lutreservoir_mswell_futuregen_dream.py --run 5 --metric Dissolved_CO2
'''
# @author: Diana Bacon
# diana.bacon@pnnl.gov

import sys
import os
import argparse
import logging
import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import ticker

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, ReservoirDataInterpolator, LookupTableReservoir,
                     MultisegmentedWellbore, FutureGen2AZMI, FutureGen2Aquifer,
                     RateToMassAdapter)


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
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
        logging.error(msg)

    # Input arguments
    parser = argparse.ArgumentParser(description='Calculate time to first detection for DREAM')
    parser.add_argument('--run', default='5', help='run number to process')
    parser.add_argument('--metric', default='Dissolved_CO2',
                        choices=['Pressure', 'Temperature', 'pH', 'TDS', 'Dissolved_CO2'])
    args = parser.parse_args()
    run = int(args.run)
    metric = args.metric

    output_directory = os.sep.join([
        '..', '..', 'Output', 'DREAM_{date_time_stamp}'.format(
            date_time_stamp=datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))])
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Define keyword arguments of the system model
    num_years = 70
    time_array = 365.25*np.arange(0.0, num_years+1, 5)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

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

    intpr = sm.add_interpolator(ReservoirDataInterpolator(
        name='int'+str(ind+1), parent=sm,
        header_file_dir=os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                     'lookuptables', 'FutureGen2', '1008_sims']),
        time_file='time_points.csv',
        data_file='fg{ind1}.csv'.format(ind1=ind+1),
        index=int(signature_data[ind+1, 0]),
        signature=signature), intr_family='reservoir')

    # Setup leak locations
    xloc = []
    yloc = []
    # Stratigraphic well
    xloc.append(782979. * 0.3048)
    yloc.append(14471060. * 0.3048)
    # Injection Well
    xloc.append(776765. * 0.3048)
    yloc.append(14468857. * 0.3048)

    ltres = []
    for well, value in enumerate(xloc):
        # Add reservoir component
        ltres.append(sm.add_component_model_object(LookupTableReservoir(
            name='ltres{}'.format(well+1), parent=sm,
            intr_family='reservoir', locX=xloc[well], locY=yloc[well])))

        # Add parameters of reservoir component model
        for j in range(num_pars):
            ltres[-1].add_par(par_names[j], value=float(signature_data[run, j+1]), vary=False)

        # Add observations of reservoir component model
        ltres[-1].add_obs('pressure')
        ltres[-1].add_obs('CO2saturation')
        ltres[-1].add_obs_to_be_linked('pressure')
        ltres[-1].add_obs_to_be_linked('CO2saturation')

    num_aquifers = 4
    # Define aquifer properties for Ironton-Galesville, Potosi, New Richmond, St Peter
    aquifer_name = ['Ironton-Galesville', 'Potosi', 'New Richmond', 'St Peter']
    aqu_thick = [33.2, 84.1, 31.1, 61.6]
    por = [0.118, 0.038, 0.132, 0.18]
    log_permh = [-13.39, -11.05, -12.48, -11.92]
    log_aniso = [0.30, 1.00, 0.30, 0.30]
    rel_vol_frac_calcite = [0.1, 0.5, 0.1, 0.1]

    num_shales = 5
    shale = [198.7, 74.4, 110.3, 118.9, 530.4]

    depth = [0, 0, 0, 0]  # list with four elements to be defined below
    depth[3] = shale[4] + aqu_thick[3]
    depth[2] = depth[3] + shale[3] + aqu_thick[2]
    depth[1] = depth[2] + shale[2] + aqu_thick[1]
    depth[0] = depth[1] + shale[1] + aqu_thick[0]
    # depth is [1044.0, 936.4, 742.0, 592.0]

    ms = []
    for well, value in enumerate(xloc):
        # Add multisegmented wellbore component
        ms.append(sm.add_component_model_object(
            MultisegmentedWellbore(name='ms{}'.format(well+1), parent=sm)))

        # Add parameters of multisegmented wellbore component
        ms[-1].add_par('wellRadius', value=0.05715, vary=False) # 4.5 inch
        ms[-1].add_par('numberOfShaleLayers', value=5, vary=False)
        ms[-1].add_par('shale1Thickness', value=shale[0], vary=False)
        ms[-1].add_par('shale2Thickness', value=shale[1], vary=False)
        ms[-1].add_par('shale3Thickness', value=shale[2], vary=False)
        ms[-1].add_par('shale4Thickness', value=shale[3], vary=False)
        ms[-1].add_par('shale5Thickness', value=shale[4], vary=False)
        ms[-1].add_par('aquifer1Thickness', value=aqu_thick[0], vary=False)
        ms[-1].add_par('aquifer2Thickness', value=aqu_thick[1], vary=False)
        ms[-1].add_par('aquifer3Thickness', value=aqu_thick[2], vary=False)
        ms[-1].add_par('aquifer4Thickness', value=aqu_thick[3], vary=False)
        ms[-1].add_par('reservoirThickness', value=7., vary=False)
        ms[-1].add_par('logAqu1Perm', min=log_permh[0]-0.5, max=log_permh[0]+0.5,
                       value=log_permh[0]) # Ironton-Galesville
        ms[-1].add_par('logAqu2Perm', min=log_permh[1]-0.5, max=log_permh[1]+0.5,
                       value=log_permh[1]) # Potosi
        ms[-1].add_par('logAqu3Perm', min=log_permh[2]-0.5, max=log_permh[2]+0.5,
                       value=log_permh[2]) # New Richmond
        ms[-1].add_par('logAqu4Perm', min=log_permh[3]-0.5, max=log_permh[3]+0.5,
                       value=log_permh[3]) # St Peter
        ms[-1].add_par('logWellPerm', min=-14.0, max=-11.0, value=-13.0)
        ms[-1].add_par('brineDensity', value=1030.9, vary=False)
        ms[-1].add_par('CO2Density', value=775.0, vary=False)
        ms[-1].add_par('brineViscosity', value=7.5e-4, vary=False)
        ms[-1].add_par('CO2Viscosity', value=6.6e-5, vary=False)

        # Add keyword arguments linked to the output provided by reservoir model
        ms[-1].add_kwarg_linked_to_obs('pressure', ltres[well].linkobs['pressure'])
        ms[-1].add_kwarg_linked_to_obs('CO2saturation', ltres[well].linkobs['CO2saturation'])

        # Add observations of multisegmented wellbore component model
        for num in range(1, num_aquifers+1):
            ms[-1].add_obs_to_be_linked('CO2_aquifer{}'.format(num))
            ms[-1].add_obs_to_be_linked('brine_aquifer{}'.format(num))
            ms[-1].add_obs_to_be_linked('mass_CO2_aquifer{}'.format(num))
        ms[-1].add_obs_to_be_linked('brine_atm')
        ms[-1].add_obs_to_be_linked('CO2_atm')

    adapt = []
    for well, value in enumerate(xloc):
        # Add adapter that transforms leakage rates to accumulated mass
        adapt.append(sm.add_component_model_object(
            RateToMassAdapter(name='adapt'+str(well+1), parent=sm)))

        for num in range(1, num_aquifers):
            obs_nm1 = 'CO2_aquifer{}'.format(num)
            obs_nm2 = 'CO2_aquifer{}'.format(num+1)
            adapt[-1].add_kwarg_linked_to_collection(
                obs_nm1, [ms[well].linkobs[obs_nm1], ms[well].linkobs[obs_nm2]])

        obs_nm = 'CO2_aquifer{}'.format(num_aquifers)
        adapt[-1].add_kwarg_linked_to_collection(
            obs_nm, [ms[well].linkobs[obs_nm], ms[well].linkobs['CO2_atm']])

        for num in range(1, num_aquifers):
            obs_nm1 = 'brine_aquifer{}'.format(num)
            obs_nm2 = 'brine_aquifer{}'.format(num+1)
            adapt[-1].add_kwarg_linked_to_collection(
                obs_nm1, [ms[well].linkobs[obs_nm1], ms[well].linkobs[obs_nm2]])

        obs_nm = 'brine_aquifer{}'.format(num_aquifers)
        adapt[-1].add_kwarg_linked_to_collection(
            obs_nm, [ms[well].linkobs[obs_nm], ms[well].linkobs['brine_atm']])

        for num in range(1, num_aquifers+1):
            adapt[-1].add_obs_to_be_linked('mass_CO2_aquifer{}'.format(num))
            adapt[-1].add_obs_to_be_linked('mass_brine_aquifer{}'.format(num))
            adapt[-1].add_obs('mass_CO2_aquifer{}'.format(num))
            adapt[-1].add_obs('mass_brine_aquifer{}'.format(num))

    fga = []
    CO2_rate_obs_list = []
    brine_rate_obs_list = []
    CO2_mass_obs_list = []
    brine_mass_obs_list = []
    for aqz in range(num_aquifers):
        fga.append([])
        for well in range(len(xloc)):
            # Add futuregen AZMI or aquifer model object and define parameters
            if depth[aqz] > 700:
                fga[-1].append(sm.add_component_model_object(
                    FutureGen2AZMI(name='fga{}{}'.format(aqz+1, well+1), parent=sm)))
            else:
                fga[-1].append(sm.add_component_model_object(
                    FutureGen2Aquifer(name='fga{}{}'.format(aqz+1, well+1), parent=sm)))
            fga[-1][-1].add_par('aqu_thick', value=aqu_thick[aqz], vary=False)
            fga[-1][-1].add_par('depth', value=depth[aqz], vary=False)
            fga[-1][-1].add_par('por', value=por[aqz], vary=False)
            fga[-1][-1].add_par('log_permh', value=log_permh[aqz], vary=False)
            fga[-1][-1].add_par('log_aniso', value=log_aniso[aqz], vary=False)
            fga[-1][-1].add_par('rel_vol_frac_calcite',
                                value=rel_vol_frac_calcite[aqz], vary=False)

            # Add aquifer component's keyword arguments
            # linked to the observations of wellbore and adapter components
            fga[-1][-1].add_kwarg_linked_to_obs(
                'co2_rate', ms[well].linkobs['CO2_aquifer{}'.format(aqz+1)])
            fga[-1][-1].add_kwarg_linked_to_obs(
                'brine_rate', ms[well].linkobs['brine_aquifer{}'.format(aqz+1)])
            fga[-1][-1].add_kwarg_linked_to_obs(
                'co2_mass', adapt[well].linkobs['mass_CO2_aquifer{}'.format(aqz+1)])
            fga[-1][-1].add_kwarg_linked_to_obs(
                'brine_mass', adapt[well].linkobs['mass_brine_aquifer{}'.format(aqz+1)])

            # Add observations (output) of the AZMI component
            fga[-1][-1].add_obs('pH_dx')
            fga[-1][-1].add_obs('pH_dy')
            fga[-1][-1].add_obs('pH_dz')
            fga[-1][-1].add_obs('Pressure_dx')
            fga[-1][-1].add_obs('Pressure_dy')
            fga[-1][-1].add_obs('Pressure_dz')
            fga[-1][-1].add_obs('TDS_dx')
            fga[-1][-1].add_obs('TDS_dy')
            fga[-1][-1].add_obs('TDS_dz')
            fga[-1][-1].add_obs('Dissolved_CO2_dx')
            fga[-1][-1].add_obs('Dissolved_CO2_dy')
            fga[-1][-1].add_obs('Dissolved_CO2_dz')
            if aqz != num_aquifers-1:
                fga[-1][-1].add_obs('Temperature_dx')
                fga[-1][-1].add_obs('Temperature_dy')
                fga[-1][-1].add_obs('Temperature_dz')

    print('------------------------------------------------------------------')
    print('                                UQ ')
    print('------------------------------------------------------------------')

    num_samples = 24
    realization = str((run-1)*num_samples+1)
    ncpus = 8
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=530242488)
    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)

    # print('------------------------------------------------------------------')
    # print('                               Plots ')
    # print('------------------------------------------------------------------')

#    def plot_time_series(ax, x, y, xlabel, ylabel, title, fontsize=12):
#        ax.plot(x, y)
#        ax.locator_params(nbins=3)
#        ax.set_xlabel(xlabel, fontsize=fontsize)
#        ax.set_ylabel(ylabel, fontsize=fontsize)
#        ax.set_title(title, fontsize=fontsize)
#
#    def plot_ms_results(ind_list, well, result):
#        fig, axs = plt.subplots(nrows=2, ncols=2, constrained_layout=True)
#        for aqz, ax in enumerate(axs.flatten()):
#            aquifer_rates = np.array([s.recarray['ms{}.{}_aquifer{}_{}'.format(
#                well+1, result, aqz+1, indd)] for indd in ind_list])
#            outfilename = os.sep.join([
#                output_directory, '{}_m{}_r{}_a{}_w{}.csv'.format(
#                    result.lower(), metric, realization, aqz+1, well)])
#            np.savetxt(outfilename, aquifer_rates, delimiter=",")
#            plot_time_series(ax, time_array/365.25, aquifer_rates,
#                             'Time, yr', '{} rate, kg/s'.format(result.capitalize()),
#                             aquifer_name[aqz])
#
#        plotfilename = os.sep.join([
#            output_directory, '{}rate_m{}_r{}_a{}_w{}.pdf'.format(
#                result.lower(), metric, realization, aqz+1, well)])
#
#        plt.savefig(plotfilename, dpi=300)
#        plt.close('all')
#
#    def plot_aquifer_results(ind_list, well, result):
#        fig, axs = plt.subplots(nrows=2, ncols=2, constrained_layout=True)
#        for aqz, ax in enumerate(axs.flatten()):
#            if (aqz == num_aquifers-1 and metric == 'Temperature'):
#                continue
#            dim = np.array([s.recarray['fga{}{}.{}_{}_{}'.format(
#                aqz+1, well+1, metric, result, indd)] for indd in ind_list])
#
#            outfilename = os.sep.join([
#                output_directory, '{}_m{}_r{}_a{}_w{}.csv'.format(
#                    result, metric, realization, aqz+1, well)])
#
#            np.savetxt(outfilename, dim, delimiter=",")
#            plot_time_series(ax, time_array/365.25, dim,
#                             'Time, yr', '{}, m'.format(result), aquifer_name[aqz])
#
#        filename = os.sep.join([output_directory, '{}_m{}_r{}_w{}.pdf'.format(
#            result, metric, realization, well)])
#
#        plt.savefig(filename, dpi=300)
#        plt.close('all')
#
#    inputs = []
#    ind_list = list(range(len(time_array)))
#    for well, value in enumerate(xloc):
#        for result in ['CO2', 'brine']:
#            inputs.append((ind_list, well, result))
#
#    with Pool(processes=4) as pool:
#        pool.starmap(plot_ms_results, inputs)
#
#    inputs = []
#    ind_list = list(range(len(time_array)))
#    for well, value in enumerate(xloc):
#        for result in ['dx', 'dy', 'dz']:
#            inputs.append((ind_list, well, result))
#
#    with Pool(processes=ncpus) as pool:
#        pool.starmap(plot_aquifer_results, inputs)

    print('------------------------------------------------------------------')
    print('                              Grid ')
    print('------------------------------------------------------------------')

    # Hardwiring for now, should add to ROM as a characteristic or parameter, DO NOT CHANGE
    threshold = {'pH': 0.2, 'Pressure': 0.00065, 'Temperature': 0.0003,
                 'TDS': 0.1, 'Dissolved_CO2': 0.2}
    indicator = {'pH': 'absolute', 'Pressure': 'relative',
                 'Temperature': 'relative', 'TDS': 'relative',
                 'Dissolved_CO2': 'relative'}

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Define grid coordinates
    xmin = 235758
    xmax = 239652
    ymin = 4409108
    ymax = 4411779
    xstep = 10
    ystep = 10
    x = np.arange(xmin, xmax, xstep)
    y = np.arange(ymin, ymax, ystep)[:, None]  # y values of interest, as a "column" array

    # Define layers based on system model
    z1 = np.linspace(-depth[0], -depth[0]+aqu_thick[0], 10)
    z2 = np.linspace(z1[-1], z1[-1]+shale[1], 3)
    z3 = np.linspace(-depth[1], -depth[1]+aqu_thick[1], 10)
    z4 = np.linspace(z3[-1], z3[-1]+shale[2], 3)
    z5 = np.linspace(-depth[2], -depth[2]+aqu_thick[2], 10)
    z6 = np.linspace(z5[-1], z5[-1]+shale[3], 3)
    z7 = np.linspace(-depth[3], -depth[3]+aqu_thick[3], 10)
    z = np.concatenate((z1[:-1], z2[:-1], z3[:-1], z4[:-1],
                        z5[:-1], z6[:-1], z7[:-1]), axis=0)[:, None, None]

    # Write grid in feet for DREAM
    filename = os.sep.join([output_directory, 'iam.grid'])
    if not os.path.isfile(filename):
        f = open(filename, "w+")
        ix = -1
        for xx in x:
            ix = ix + 1
            f.write('{},'.format(xx/.3048))
        f.write('\n')

        iy = -1
        for yy in y:
            iy = iy + 1
            f.write('{},'.format(yy[0]/.3048))
        f.write('\n')

        iz = -1
        for zz in z:
            iz = iz + 1
            f.write('{},'.format(zz[0, 0]/.3048))
        f.write('\n')
        f.close()

    print('------------------------------------------------------------------')
    print('                               Map ')
    print('------------------------------------------------------------------')

    def map_plumes(sample):

        realization = str((run-1)*num_samples+sample+1)
        print(metric, 'run', run, 'sample', sample, 'realization', realization)

        # Initialize time to first detection grid to large number
        ttfd = np.ones((len(z), len(y), len(x)))*1.0e30

        for aqz in range(num_aquifers):
            if (aqz == num_aquifers-1 and metric == 'Temperature'):
                continue
            z0 = -depth[aqz]+aqu_thick[aqz]
            for well in range(len(xloc)):
                x0 = ltres[well].locX
                y0 = ltres[well].locY

                # Get plume dimensions
                for indd, tt in enumerate(time_array):
                    a = s.recarray['fga{}{}.{}_dx_{}'.format(
                        aqz+1, well+1, metric, indd)][sample]
                    b = s.recarray['fga{}{}.{}_dy_{}'.format(
                        aqz+1, well+1, metric, indd)][sample]
                    c = s.recarray['fga{}{}.{}_dz_{}'.format(
                        aqz+1, well+1, metric, indd)][sample]

                    # Calculate time to first detection
                    if (a > 0 and b > 0 and c > 0):
                        ellipsoid = ((x-x0)/(a*0.5))**2 + ((y-y0)/(b*0.5))**2 + (
                            (z-z0)/(c*0.5))**2 <= 1
                        in_aquifer = z-z0 < 0
                        half_ellipsoid = np.logical_and(ellipsoid, in_aquifer)
                        mask = np.logical_not(half_ellipsoid)*1.0e30
                        plume = half_ellipsoid * tt
                        ttfd = np.minimum(ttfd, mask + plume)

        # Write output (in feet) for DREAM
        filename = os.sep.join([output_directory, 'ttfd_{}_{}.iam'.format(
            metric, realization)])
        f = open(filename, "w+")
        f.write('IAM,{},{},{},{},\n'.format(
            realization, metric, indicator[metric], threshold[metric]))
        ix = -1
        for xx in x:
            ix = ix + 1
            iy = -1
            for yy in y:
                iy = iy + 1
                iz = -1
                for zz in z:
                    iz = iz + 1
                    if ttfd[iz][iy][ix] != 1.0e+30:
                        f.write('{},{},{},{}\n'.format(
                            xx/.3048, yy[0]/.3048, zz[0, 0]/.3048, ttfd[iz][iy][ix]))
        f.close()

        # Make plots
        fig, ax = plt.subplots()
        X, Y = np.meshgrid(x, y)
        cp = plt.contourf(X/1000.0, Y/1000.0, ttfd[-1, :, :], levels=time_array)
        plt.colorbar(cp)
        ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
        fmt = ticker.FuncFormatter(lambda x, p: format(int(x), ','))
        plt.xlabel('Easting (km)')
        plt.ylabel('Northing (km)')
        plt.gca().get_xaxis().set_major_formatter(fmt)
        plt.gca().get_yaxis().set_major_formatter(fmt)
        plt.title('Realization {}\nTime to Detection (days)'.format(realization))
        plt.tight_layout()
        to_save = True
        if to_save:
            file_nm = 'ttfd_{}_{}{}'
            plt.savefig(os.sep.join([output_directory, file_nm.format(
                metric, realization, '.png')]), dpi=300)

        plt.close()

    from sys import platform
    if platform == "win32":
        # Loop replaces parallel execution of the simulations on Windows
        for s_ind in range(num_samples):
            map_plumes(s_ind)
    else:
        # The following code should work on Mac and Linux but not on Windows.
        from multiprocessing import Pool
        with Pool(processes=ncpus) as pool:
            pool.map(map_plumes, list(range(num_samples)))
