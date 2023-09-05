"""
Example illustrates simple linking between NRAP-Open-IAM and SOSAT.
For NRAP-Open-IAM we simply use the lookup table component to get the pressure field that
will be used in SOSAT.

This example requires the additional Kimberlina data set.
Kimberlina data set (pressure-in-pa-kimb-54-sims.zip or
Pressure_In_Pa_Kimb_54_sims.zip) can be downloaded from one of the following places:
1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam
2. https://gitlab.com/NRAP/Kimberlina_data

The downloaded (unzipped) data set should be placed here:
    source/components/reservoir/lookuptables/Kimb_54_sims

This example also requires that you have SOSAT package installed.
You can install SOSAT by running this command on the terminal:
pip install SOSAT

Example of run:
$ python iam_sys_lutreservoir_sosat.py
"""
import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import uniform

from SOSAT import StressState
# from SOSAT.constraints import StressMeasurement
from SOSAT.constraints import FaultConstraint
from SOSAT.constraints import FaultingRegimeConstraint, SU
from SOSAT.risk_analysis import CriticalFaultActivation

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, LookupTableReservoir

if __name__ == "__main__":

    # Setup location of lookup table data set
    file_directory = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                  'lookuptables', 'Kimb_54_sims'])

    if not os.path.exists(os.sep.join([file_directory, 'Reservoir_data_sim01.csv'])):
        msg = ''.join(['\nKimberlina data set can be downloaded ',
                       'from one of the following places:\n',
                       '1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam \n',
                       '2. https://gitlab.com/NRAP/Kimberlina_data  \n',
                       'Check this example description for more information.'])
        logging.error(msg)

    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25 * np.arange(0.0, num_years+1)
    seal_model_kwargs = {'time_array': time_array} # time is given in days

    # Define file type for saving gridded observations
    grid_save_type = 'npz'

    # Create system model
    sm = SystemModel(model_kwargs=seal_model_kwargs)

    # Read file with signatures of interpolators and names of files
    # with the corresponding data.
    signature_data = np.genfromtxt(
        os.path.join(file_directory, 'parameters_and_filenames.csv'),
        delimiter=",", dtype='str')

    import openiam.enmesh as en
    # Get xmin and xmax from a random lookup table
    # Get filename of lookuptable
    lookup_filename = 'Reservoir_data_sim05.csv'
    file = os.sep.join([file_directory, lookup_filename])
    x, y = np.loadtxt(file, skiprows=1, delimiter=',', unpack=True, usecols=[0, 1])
    f2, ax = plt.subplots()
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95, wspace=0.25)
    ax.plot(x, y, 'kx')
    ax.set_aspect('equal')
    plt.show()

    # Setting up a grid to plot results
    auxDepth = 0.
    nx = 21
    ny = 27
    minX = 13000
    maxX = 60000
    minY = 35000
    maxY = 65000
    min_x_spacing = (maxX - minX) / nx
    min_y_spacing = (maxY - minY) / ny
    grid = en.Grid(xmin=minX, xmax=maxX, ymin=minY, ymax=maxY, zmin=auxDepth, zmax=10.)
    grid.refine_grid_around_wells(min_x_spacing=min_x_spacing,
                                  min_y_spacing=min_y_spacing, x_mul=1.0, y_mul=1.0)
    grid.plot_xy()
    xx, yy, zz = grid.get_vertices()
    nVertices = xx.shape[0] * xx.shape[1] * xx.shape[2]
    nVerticesXY = xx.shape[0] * xx.shape[1]
    auxX = xx[:, :, 0]
    loc_X = np.reshape(auxX, nVerticesXY)
    auxY = yy[:, :, 1]
    loc_Y = np.reshape(auxY, nVerticesXY)

    # Add reservoir component
    ltres = sm.add_component_model_object(
        LookupTableReservoir(name='ltres', parent=sm, locX=loc_X, locY=loc_Y))

    ltres.build_and_link_interpolators(file_directory=file_directory,
                                       intr_family='reservoir',
                                       default_values={'salinity': 0.1,
                                                       'temperature': 50.0},
                                       recompute_triangulation=False,
                                       build_on_the_fly=False)

    # Add parameter of the reservoir component indicating what lookup table
    # will be used
    ltres.add_par('index', value=5, vary=False)

    # Setup output folder to keep data files with gridded observations
    output_dir = os.sep.join(['..', '..', 'output', 'sosat_sim_data'])
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Add gridded observations of the reservoir component
    ltres.add_grid_obs('pressure', constr_type='array',
                       output_dir=output_dir, save_type=grid_save_type)

    print('Starting NRAP-Open-IAM simulation.')
    # Run forward simulation
    sm.forward()
    print('NRAP-Open-IAM simulation is finished.')

    # Collect gridded observations
    pressure_grid = sm.collect_gridded_observations_as_time_series(
        ltres, 'pressure', output_dir, rlzn_number=0, save_type=grid_save_type)

    print('Starting SOSAT analysis.')

    # Pick a time (in years) to analyze results in SOSAT
    time_step = 50

    # Since stratigraphy not present, estimate depth based on initial pressure
    brine_density = 1080.0
    depth = (pressure_grid[0]-101325.0)/(9.8*brine_density)

    pore_pressure=pressure_grid[time_step]/1.0e6  # conversion from Pa to MPa

    # Estimate overburden density in kg/m^3
    avg_overburden_density = np.full_like(pore_pressure, 2580.0)
    Pfail = np.full_like(pore_pressure, -1)

    for i, dval in enumerate(depth):
        ss = StressState(depth=dval,
                         avg_overburden_density=avg_overburden_density[i],
                         pore_pressure=pore_pressure[i])

        fc = FaultConstraint()
        ss.add_constraint(fc)

        # Create a faulting regime constraint
        frc = FaultingRegimeConstraint(SU(w_NF=0.01, w_SS=0.24, w_TF=0.75,
                                          theta1=np.sqrt(2.0) * 0.5, k1=300.0,
                                          theta2=-np.sqrt(2.0) * 0.5, k2=100.0))
        ss.add_constraint(frc)

        # Add stress measurement constraint if there is any field data
        # smc = StressMeasurement(shmin_dist=uniform(loc=25.0,
        #                                            scale=5.0))
        # ss.add_constraint(smc)

        dPmax = 0.0
        gamma_dist = uniform(0.4, (0.6 - 0.4))
        cfa = CriticalFaultActivation(ss, dPmax, gamma_dist)
        pressures, Pfail[i] = cfa.EvaluatePfail(Npressures=1, Nsamples=1e5)

    print('SOSAT analysis is finished.')

    print('------------------------------------------------------------------')
    print('                  Start plotting ')
    print('------------------------------------------------------------------')

    f2, ax = plt.subplots()
    t = ax.tricontourf(loc_X, loc_Y, Pfail, vmin=0.2295, vmax=0.255)
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    plt.title('Probability of failure')
    t.set_clim(vmin=0.2295, vmax=0.255)
    plt.colorbar(t, ax=ax)
    plt.savefig(output_dir+'/probFailure_50years.pdf', dpi=300)
    plt.show()

    deltaPressure = (pressure_grid[time_step] - pressure_grid[0])/1.0e6
    f2, ax = plt.subplots()
    t = ax.tricontourf(loc_X, loc_Y, deltaPressure)
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    plt.title('Change in pore pressure [MPa]')
    f2.colorbar(t)
    plt.savefig(output_dir+'/pressure_50years.pdf', dpi=300)
    plt.show()
