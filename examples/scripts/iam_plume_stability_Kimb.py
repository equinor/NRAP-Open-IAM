"""
Example illustrates plume stability analysis functionality.

This example requires the additional Kimberlina data set.
Kimberlina data set (pressure-in-pa-kimb-54-sims.zip or
Pressure_In_Pa_Kimb_54_sims.zip) can be downloaded from one of the following places:
1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam
2. https://gitlab.com/NRAP/Kimberlina_data

The downloaded data set should be placed here:
    source/components/reservoir/lookuptables/Kimb_54_sims

Example of run:
$ python iam_plume_stability_Kimb.py
"""

import sys
import os
import matplotlib.pyplot as plt
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam.mesh2D import read_Mesh2D_data
from openiam.reservoir_data_interpolator import read_time_points


if __name__ == '__main__':
    ##############
    # User input #
    ##############
    datafolder = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                              'lookuptables', 'Kimb_54_sims'])
    if not os.path.exists(os.sep.join([datafolder, 'Reservoir_data_sim01.csv'])):
        msg = ''.join(['\nKimberlina data set can be downloaded ',
                       'from one of the following places:\n',
                       '1. https://edx.netl.doe.gov/dataset/nrap-open-source-iam \n',
                       '2. https://gitlab.com/NRAP/Kimberlina_data  \n',
                       'Check this example description for more information.'])
        print(msg)

    dP_thresh = 1.e6   # Overpressure threshold [Pa]
    S_thresh = 0.01    # CO2 saturation threshold [m^3/m^3]
    ##################
    # End user input #
    ##################

    ##################
    # PREPARING DATA #
    ##################

    # Read file with time points data
    _, time_points = read_time_points(False, datafolder, 'time_points.csv')
    # Read a single realization into a Mesh2D object
    M = read_Mesh2D_data(os.path.join(datafolder, 'Reservoir_data_sim01.csv'),
                         ['pressure', 'CO2saturation'], time_points, False)

    # Not necessary, but handy to create independent python variables
    # to store Mesh2D object variables for ease of use below
    P = M.variables['pressure']
    S = M.variables['CO2saturation']

    # Create overpressure variable by subtracting each pressure dataset from the first pressure dataset
    # This is a good example of how derived variables can be created from the dataset.
    O = M.add_variable('overpressure')
    for t, d in P.datasets.items():
        O.add_dataset(d.data - P.datasets[P.times[0]].data, t)

    ############
    # Plotting #
    ############
    # Now you can utilize the methods of the class to facilitate plotting plume stability metrics

    # Plot plume areas
    f, ax = plt.subplots(1)
    ax.plot(O.times, O.plume_areas(dP_thresh)/1000**2)
    ax.set_xlabel('Time [y]')
    ax.set_ylabel(r'Pressure plume area [km$^2$]')
    f.show()
    f, ax = plt.subplots(1)
    ax.plot(S.times, S.plume_areas(S_thresh)/1000**2)
    ax.set_xlabel('Time [y]')
    ax.set_ylabel(r'CO$_2$ area [km$^2$]')
    f.show()

    # Plot temporal derivative of plume area
    f, ax = plt.subplots(1)
    ax.plot(O.times, O.plume_areas_dt(dP_thresh)/1000**2, label='Pressure')
    ax.plot(S.times, S.plume_areas_dt(S_thresh)/1000**2, label=r'CO$_2$')
    ax.legend()
    ax.set_xlabel('Time [y]')
    ax.set_ylabel(r'dA/dt [km$^2$/y]')
    f.show()

    # Plot temporal derivative of plume centriod
    f, ax = plt.subplots(1)
    ax.plot(O.times, O.plume_centroids_dt(dP_thresh)[:, 0], label='x direction')
    ax.plot(O.times, O.plume_centroids_dt(dP_thresh)[:, 1], label='y direction')
    ax.legend()
    ax.set_xlabel('Time [y]')
    ax.set_ylabel(r'Pressure plume mobility [m/y]')
    f.show()
    f, ax = plt.subplots(1)
    ax.plot(S.times, S.plume_centroids_dt(S_thresh)[:, 0], label='x direction')
    ax.plot(S.times, S.plume_centroids_dt(S_thresh)[:, 1], label='y direction')
    ax.legend()
    ax.set_xlabel('Time [y]')
    ax.set_ylabel(r'Saturation plume mobility [m/y]')
    f.show()

    # Plot temporal derivative of plume spread
    f, ax = plt.subplots(1)
    ax.plot(O.times, O.plume_spreads_dt(dP_thresh)[:, 0]/1000**2, label=r'$\lambda$1')
    ax.plot(O.times, O.plume_spreads_dt(dP_thresh)[:, 1]/1000**2, label=r'$\lambda$2')
    ax.legend()
    ax.set_xlabel('Time [y]')
    ax.set_ylabel(r'Pressure plume spreading [km$^2$/y]')
    f.show()
    f, ax = plt.subplots(1)
    ax.plot(S.times, S.plume_spreads_dt(S_thresh)[:, 0]/1000**2, label=r'$\lambda$1')
    ax.plot(S.times, S.plume_spreads_dt(S_thresh)[:, 1]/1000**2, label=r'$\lambda$2')
    ax.legend()
    ax.set_xlabel('Time [y]')
    ax.set_ylabel(r'Saturation plume spreading [km$^2$/y]')
    f.show()

    #############
    # Animation #
    #############
    # This currently requires ffmpeg to be installed
    # These wrap up all the metrics into one animation
    # Note: filename='filename.mp4' can be added to these calls to save the animations.

    O.plume_stability_animation(dP_thresh)
    S.plume_stability_animation(S_thresh)

    #################
    # Interpolation #
    #################
    # Scipy efficiently interpolates on a numpy meshgrid
    # The following uses the scipy functionality to plot interpolated values over time

    # Choose some points to interpolate at
    pts = [[30000, 60000], [32000, 61000], [29000, 58000]]

    f, ax = plt.subplots(1)
    ax.plot(O.times, O.interpolate(pts)[:, 0], label='Point 1')
    ax.plot(O.times, O.interpolate(pts)[:, 1], label='Point 2')
    ax.plot(O.times, O.interpolate(pts)[:, 2], label='Point 3')
    ax.legend()
    ax.set_xlabel('Time [y]')
    ax.set_ylabel(r'Pressure [Pa]')
    f.show()

    f, ax = plt.subplots(1)
    ax.plot(S.times, S.interpolate(pts)[:, 0], label='Point 1')
    ax.plot(S.times, S.interpolate(pts)[:, 1], label='Point 2')
    ax.plot(S.times, S.interpolate(pts)[:, 2], label='Point 3')
    ax.legend()
    ax.set_xlabel('Time [y]')
    ax.set_ylabel(r'CO$_2$ saturation [-]')
    f.show()
