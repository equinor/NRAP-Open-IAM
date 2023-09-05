'''
This example couples the analytical reservoir, multisegmented wellbore and
generic aquifer models. The saturation/pressure output produced by analytical
reservoir model is used to drive leakage from a single multisegmented wellbore
model, which is passed to the input of an adapter that provides well
coordinates, |CO2| and brine leakage rates and cumulative mass fluxes to the
generic aquifer model. HDF5 files for input to DREAM are created.

Example of run:
$ python iam_sys_reservoir_mswell_generic_lhs.py
'''

import sys
import os
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import h5py
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, AnalyticalReservoir, MultisegmentedWellbore,
                     GenericAquifer, RateToMassAdapter)
import openiam.enmesh as en


if __name__ == "__main__":
    # For multiprocessing in Spyder
    __spec__ = None
    # Define keyword arguments of the system model
    num_years = 10
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}   # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # legacy well location
    xloc = 200
    yloc = 200

    # Add reservoir component
    ares = sm.add_component_model_object(AnalyticalReservoir(
        name='ares', parent=sm, injX=0., injY=0., locX=xloc, locY=yloc))

    # Add parameters of reservoir component model
    ares.add_par('numberOfShaleLayers', value=3, vary=False)
    ares.add_par('shale1Thickness', value=100.0, vary=False)
    ares.add_par('aquifer1Thickness', value=100.0, vary=False)
    ares.add_par('shale2Thickness', value=100.0, vary=False)
    ares.add_par('aquifer2Thickness', value=100.0, vary=False)
    ares.add_par('shale3Thickness', value=500.0, vary=False)
    ares.add_par('injRate', value=1.0, vary=False)

    # Add observations of reservoir component model
    ares.add_obs_to_be_linked('pressure')
    ares.add_obs_to_be_linked('CO2saturation')
    ares.add_obs('pressure')
    ares.add_obs('CO2saturation')
    ares.add_obs('mass_CO2_reservoir')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))
    ms.add_par('logWellPerm', min=-14.0, max=-12.0, value=-13.0)

    # Add linked parameters: common to reservoir and wellbore components
    ms.add_par_linked_to_par('numberOfShaleLayers',
                             ares.deterministic_pars['numberOfShaleLayers'])
    ms.add_par_linked_to_par('shale1Thickness',
                             ares.deterministic_pars['shale1Thickness'])
    ms.add_par_linked_to_par('shale2Thickness',
                             ares.deterministic_pars['shale2Thickness'])
    ms.add_par_linked_to_par('shale3Thickness',
                             ares.deterministic_pars['shale3Thickness'])
    ms.add_par_linked_to_par('aquifer1Thickness',
                             ares.deterministic_pars['aquifer1Thickness'])
    ms.add_par_linked_to_par('aquifer2Thickness',
                             ares.deterministic_pars['aquifer2Thickness'])
    ms.add_par_linked_to_par('reservoirThickness',
                             ares.default_pars['reservoirThickness'])
    ms.add_par_linked_to_par('datumPressure',
                             ares.default_pars['datumPressure'])

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', ares.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', ares.linkobs['CO2saturation'])

    # Add observations of multisegmented wellbore component model
    ms.add_obs_to_be_linked('CO2_aquifer1')
    ms.add_obs_to_be_linked('CO2_aquifer2')
    ms.add_obs_to_be_linked('brine_aquifer1')
    ms.add_obs_to_be_linked('brine_aquifer2')
    ms.add_obs_to_be_linked('mass_CO2_aquifer1')
    ms.add_obs_to_be_linked('mass_CO2_aquifer2')
    ms.add_obs_to_be_linked('brine_atm')
    ms.add_obs_to_be_linked('CO2_atm')
    ms.add_obs('brine_aquifer1')
    ms.add_obs('CO2_aquifer1')

    # Add adapter that transforms leakage rates to accumulated mass
    adapt = sm.add_component_model_object(RateToMassAdapter(name='adapt',
                                                            parent=sm))
    adapt.add_kwarg_linked_to_collection('CO2_aquifer1',
        [ms.linkobs['CO2_aquifer1'], ms.linkobs['CO2_aquifer2']])
    adapt.add_kwarg_linked_to_collection('CO2_aquifer2',
        [ms.linkobs['CO2_aquifer2'], ms.linkobs['CO2_atm']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer1',
        [ms.linkobs['brine_aquifer1'], ms.linkobs['brine_aquifer2']])
    adapt.add_kwarg_linked_to_collection('brine_aquifer2',
        [ms.linkobs['brine_aquifer2'], ms.linkobs['brine_atm']])
    adapt.add_obs_to_be_linked('mass_CO2_aquifer1')
    adapt.add_obs_to_be_linked('mass_CO2_aquifer2')
    adapt.add_obs_to_be_linked('mass_brine_aquifer1')
    adapt.add_obs_to_be_linked('mass_brine_aquifer2')
    adapt.add_obs('mass_CO2_aquifer1')
    adapt.add_obs('mass_brine_aquifer1')
    adapt.add_obs('mass_CO2_aquifer2')
    adapt.add_obs('mass_brine_aquifer2')

    # Add generic aquifer model object and define parameters
    ga = sm.add_component_model_object(GenericAquifer(name='ga', parent=sm))
    ga.add_par_linked_to_par('aqu_thick', ares.deterministic_pars['aquifer1Thickness'])
    ga.add_composite_par('top_depth',
        expr='ares.shale2Thickness + ares.shale3Thickness' +
        '+ ares.aquifer2Thickness')
    ga.add_par('por', value=1.965259282453879763e-01, vary=False)
    ga.add_par('log_permh', value=-1.191464515905165555e+01, vary=False)
    ga.add_par('log_aniso', value=8.046003470121247947e-01, vary=False)
    ga.add_par('aquifer_salinity', value=1.267995018132549341e-02, vary=False)
    ga.add_par('reservoir_salinity', value=4.159677791928499679e-02, vary=False)
    ga.add_par('dissolved_salt_threshold', value=0.015, vary=False)
    ga.add_par('dissolved_co2_threshold', value=0.001, vary=False)

    ga.add_kwarg_linked_to_obs('co2_mass', adapt.linkobs['mass_CO2_aquifer1'])
    ga.add_kwarg_linked_to_obs('brine_mass', adapt.linkobs['mass_brine_aquifer1'])

    # Add observations (output) from the generic aquifer model
    ga.add_obs('Dissolved_salt_volume')
    ga.add_obs('Dissolved_CO2_volume')

    # Define output folder to keep data files with gridded observations
    output_dir = 'dream_data'

    # Add gridded observations of the aquifer component
    ga.add_grid_obs('r_coordinate', constr_type='matrix', output_dir=output_dir)
    ga.add_grid_obs('z_coordinate', constr_type='matrix', output_dir=output_dir)
    ga.add_grid_obs('Dissolved_CO2_mass_fraction',
                    constr_type='matrix', output_dir=output_dir)
    ga.add_grid_obs('Dissolved_salt_mass_fraction',
                    constr_type='matrix', output_dir=output_dir)

    print('------------------------------------------------------------------')
    print('                      UQ illustration ')
    print('------------------------------------------------------------------')

    import random
    num_samples = 25
    ncpus = 1
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=random.randint(500, 1100))

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)

    print('------------------------------------------------------------------')
    print('                 Write DREAM File ')
    print('------------------------------------------------------------------')

    # set boundaries of site grid
    top = ares.deterministic_pars['aquifer2Thickness'].value + \
        ares.deterministic_pars['shale2Thickness'].value + \
            ares.deterministic_pars['shale3Thickness'].value
    bottom = top + ares.deterministic_pars['aquifer1Thickness'].value
    grid = en.Grid(xmin=0, xmax=1000, ymin=0, ymax=1000, zmin=top, zmax=bottom)

    # set constant grid spacing for vertical axis
    grid.get_axis('z').add_zone(min=grid.zmin, max=grid.zmax).add_ticks(min_spacing=10., mul=0)

    # add wells
    grid.add_well(x=xloc, y=yloc, name='Well 1')

    # refine grid around wells
    grid.refine_grid_around_wells(min_x_spacing=5, min_y_spacing=5, x_mul=1.5, y_mul=1.5)

    # get 3D mesh with radial distance from well
    xx, yy, zz = grid.get_vertices()
    cx, cy, cz = grid.get_centroids()
    cr = grid.get_radial_distance(xloc, yloc, cx, cy)
    positions = np.vstack(list(zip(cr.ravel(), cz.ravel())))

    # Read Gridded Observation Coordinates
    rcoord = np.load(os.path.join(
        output_dir, 'ga_r_coordinate_sim_1_time_0.npz'))['data']
    zcoord = np.load(os.path.join(
        output_dir, 'ga_z_coordinate_sim_1_time_0.npz'))['data']
    r_coord = rcoord[:, 0]
    z_coord = np.flip(zcoord[0, :])

    # start simulation loop
    for sim in np.arange(1, num_samples+1):
        print('Processing realization', sim)
        years = time_array/365.25
        steps = years.astype(int)

        # Open DREAM file
        hdf5 = h5py.File(os.path.join(output_dir, 'ga_sim_'+str(sim)+'.h5'), 'w')

        # write grid and geologic data
        g1 = hdf5.create_group('data')
        g1.create_dataset(
            'porosity',
            data=ga.deterministic_pars['por'].value*np.ones_like(cx),
            dtype='float32')
        g1.create_dataset('steps', data=steps, dtype='float32')
        g1.create_dataset('times', data=time_array, dtype='float32')
        g1.create_dataset('vertex-x', data=np.array(xx)[:, 0, 0], dtype='float32')
        g1.create_dataset('vertex-y', data=np.array(yy)[0, :, 0], dtype='float32')
        g1.create_dataset('vertex-z', data=np.array(zz)[0, 0, :], dtype='float32')
        g1.create_dataset('x', data=np.array(cx)[:, 0, 0], dtype='float32')
        g1.create_dataset('y', data=np.array(cy)[0, :, 0], dtype='float32')
        g1.create_dataset('z', data=np.array(cz)[0, 0, :], dtype='float32')
        g1['x'].attrs['units'] = 'm'
        g1['y'].attrs['units'] = 'm'
        g1['z'].attrs['units'] = 'm'
        g1['vertex-x'].attrs['units'] = 'm'
        g1['vertex-y'].attrs['units'] = 'm'
        g1['vertex-z'].attrs['units'] = 'm'
        g1['z'].attrs['positive'] = 'down'
        g1['vertex-z'].attrs['positive'] = 'down'

        # write gridded observations
        def write_obs(name, unit, positions, shape):

            def snake_to_camel(name):
                temp = name.split('_')
                res = ''.join(ele.title() for ele in temp)
                return res

            name2 = snake_to_camel(name)
            unit2 = snake_to_camel(unit)

            means = []
            mins = []
            maxs = []

            for step in steps:

                # read gridded observation files
                obs = np.load(os.path.join(
                    output_dir,
                    'ga_'+name+'_'+unit+'_sim_'+str(sim)+'_time_'+str(step)+'.npz'))['data']

                # Interpolate Gridded Observations to Site Grid
                interpolator = RegularGridInterpolator(
                    (r_coord, z_coord), obs, bounds_error=False)
                interpolated = interpolator(positions)
                reshaped = interpolated.reshape(shape)

                means.append(np.mean(reshaped))
                mins.append(np.min(reshaped))
                maxs.append(np.max(reshaped))

                g2 = hdf5.require_group('plot%i'%step)
                g2.create_dataset(name2, data=reshaped, dtype='float32')

                g2[name2].attrs['unit'] = unit2

            # calculate min,mean,max over all time steps
            g3 = hdf5.require_group('statistics')
            g3.create_dataset(name2,
                              data=np.array([
                                  np.min(np.array(mins)),
                                  np.mean(np.array(means)),
                                  np.max(np.array(maxs))]),
                              dtype='float32')

        write_obs('dissolved_co2', 'mass_fraction', positions, cx.shape)
        write_obs('dissolved_salt', 'mass_fraction', positions, cx.shape)

        hdf5.close()

        # end sim loop
