"""
Example illustrates plume stability analysis functionality.
While this example requires a dataset not distributed with the OpenIAM
it shows how a set of data can be read in from any format and used
to populate a Mesh2D object for plume analysis.

Example of run:
$ python iam_plume_stability.py
"""

import sys
import os
from glob import glob
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam.mesh2D import Mesh2D

##############
# User input #
##############
# Where the mesh and data files are located
datafolder = '/project/nrap/dharp/plume_stability/example'
dP_thresh = 1. # Overpressure threshold [MPa]
S_thresh = 0.01 # CO2 saturation threshold [m^3/m^3]
grid_name = 'rsu_grid.dat' # Name of grid in David Dempsey format
data_file_glob = '_*yr.dat' # Glob pattern for name of data files in David Dempsey format
output_file_name = 'RSU_plume_stability.dat' # Name of file to put results in
##################
# End user input #
##################

##################
# PREPARING DATA #
##################
# This is an example using Rock Springs Uplift files created by David Dempsey.
# The goal is to create a Mesh2D object and populate its variables with your data.
# THIS WILL NEED TO BE MODIFIED FOR DATA FILES IN OTHER FORMATS!!!

# Collect data files and make sure they are sorted by time
fls = sorted(glob(os.path.join(datafolder, data_file_glob)),
             key=lambda x: int(x.split('_')[-1].split('y')[0]))
# Limit number of files for testing
#fls = fls[0:2]

# Read in mesh
print('\nReading mesh...')
g = np.genfromtxt(os.path.join(datafolder, grid_name), delimiter=',', names=True)

# Get unique x and y coordinates, need to round to nearest 10 meters due to stacked mesh offsets
xu = np.unique(np.round(g['x_m']/10)*10)
yu = np.unique(np.round(g['y_m']/10)*10)
xyu = np.unique(g[['x_m', 'y_m']]).tolist()
# Create mesh grid
xv,yv = np.meshgrid(xu, yu, sparse=True)
# Get cell horizontal area assuming that they are the same for all cells
# (in this case, they are all the same (200 m)).
# If they are all different, a numpy meshgrid of non-identical cell areas can be constructed
# and used to create the Mesh2D object below.
dxc = np.mean(np.diff(xu))
dyc = np.mean(np.diff(yu))
#dAc = dxc*dyc # horizontal area of cells

# Calculate cell areas
# Find midpoints between nodes and add nodes on each end
xv_mid = np.column_stack(
    [xv[:, 0] - np.diff(xv)[:, 0]/2.,
     xv[:, 0:-1] + np.diff(xv)/2.,
     xv[:, -1] + np.diff(xv)[:, -1]/2.])
yv_mid = np.row_stack(
    [yv[0, :] - np.diff(yv, axis=0)[0, :]/2.,
     yv[0:-1, :] + np.diff(yv, axis=0)/2.,
     yv[-1, :] + np.diff(yv, axis=0)[-1, :]/2.])
# Multiply difference between x and y midpoint coordinates to get cell areas
dAc = np.diff(xv_mid) * np.diff(yv_mid, axis=0)

# Collect indices along x,y (i,j) columns (it only has to be done once)
print('\nCalculating col_inds. it may take a while...')
col_inds = {}
# Generalized for irregular mesh
for i, x in enumerate(xu):
    for j, y in enumerate(yu):
        if (x, y) in xyu:
            col_inds[(i, j)] = np.intersect1d(
                np.where(np.abs(g['x_m']-x)<2.0)[0],
                np.where(np.abs(g['y_m']-y)<2.0)[0])

# Create Mesh2D object and add pressure co2 saturation variable types
M = Mesh2D(xu, yu, As=dAc)
O = M.add_variable('overpressure')
S = M.add_variable('co2sat')
# Mesh2D variables can be accessed by Python variables ("O" and "S")
# or by the Mesh2D variables dictionary
# ("M.variables['overpressure'] and M.variables['co2sat']), whichever you prefer.

# Read in data at time 0 and add overpressure and saturation variables
# Assume saturations are all zero at time zero (overpressures are zero at time 0 by definition)
d0 = np.genfromtxt(fls[0], delimiter=',', names=True)
t = float(fls[0].split('_')[-1].split('yr')[0])
O.add_dataset(M.zeros_like_data(), t)
S.add_dataset(M.zeros_like_data(), t)

# Add data at each timesteps greater than zero
print('\nProcessing data...')
for fl in fls[1:]:
    print(fl)
    # Collect time from filename
    t = float(fl.split('_')[-1].split('yr')[0])
    # Read in current file
    dn = np.genfromtxt(fl, delimiter=',', names=True)
    # Collect max overpressure and saturation along each x,y column
    Pv = M.zeros_like_data()
    Sv = M.zeros_like_data()
    for xs,nds in col_inds.items():
        Pv[xs[0], xs[1]] = np.max(dn[nds]['P_MPa'] - d0[nds]['P_MPa'])
        Sv[xs[0], xs[1]] = np.max(dn[nds]['S_co2'])
	# Add datasets to Mesh2D variable
    O.add_dataset(Pv,t)
    S.add_dataset(Sv,t)

############
# Plotting #
############
# Now you can utilize the methods of the class to facilitate plotting plume stability metrics

## Plot plume areas
#f, ax = plt.subplots(1)
PA = O.plume_areas(dP_thresh)/1000**2
#ax.plot(O.times, PA)
#ax.set_xlabel('Time [y]')
#ax.set_ylabel(r'Pressure plume area [km$^2$]')
#f.show()
#f, ax = plt.subplots(1)
SA = S.plume_areas(S_thresh)/1000**2
#ax.plot(S.times, SA)
#ax.set_xlabel('Time [y]')
#ax.set_ylabel(r'CO$_2$ area [km$^2$]')
#f.show()
#
## Plot temporal derivative of plume area
#f,ax = plt.subplots(1)
dPAdt = O.plume_areas_dt(dP_thresh)/1000**2
#ax.plot(O.times, dPAdt, label='Pressure')
dSAdt = S.plume_areas_dt(S_thresh)/1000**2
#ax.plot(S.times, dSAdt, label=r'CO$_2$')
#ax.legend()
#ax.set_xlabel('Time [y]')
#ax.set_ylabel(r'dA/dt [km$^2$/y]')
#f.show()
#
## Plot temporal derivative of plume centriod
#f,ax = plt.subplots(1)
dPxcsdt = O.plume_centroids_dt(dP_thresh)[:, 0]
#ax.plot(O.times, dPxcsdt, label='x direction')
dPycsdt = O.plume_centroids_dt(dP_thresh)[:, 1]
#ax.plot(O.times, dPycsdt, label='y direction')
#ax.legend()
#ax.set_xlabel('Time [y]')
#ax.set_ylabel(r'Pressure plume mobility [m/y]')
#f.show()
#f, ax = plt.subplots(1)
dSxcsdt = S.plume_centroids_dt(S_thresh)[:, 0]
#ax.plot(S.times, dSxcsdt, label='x direction')
dSycsdt = S.plume_centroids_dt(S_thresh)[:, 1]
#ax.plot(S.times, dSycsdt, label='y direction')
#ax.legend()
#ax.set_xlabel('Time [y]')
#ax.set_ylabel(r'Saturation plume mobility [m/y]')
#f.show()
#
## Plot temporal derivative of plume spread
#f,ax = plt.subplots(1)
dPlm1dt = O.plume_spreads_dt(dP_thresh)[:, 0]/1000**2
#ax.plot(O.times,dPlm1dt,label=r'$\lambda$1')
dPlm2dt = O.plume_spreads_dt(dP_thresh)[:, 1]/1000**2
#ax.plot(O.times, dPlm2dt, label=r'$\lambda$2')
#ax.legend()
#ax.set_xlabel('Time [y]')
#ax.set_ylabel(r'Pressure plume spreading [km$^2$/y]')
#f.show()
#f, ax = plt.subplots(1)
dSlm1dt = S.plume_spreads_dt(S_thresh)[:, 0]/1000**2
#ax.plot(S.times, dSlm1dt, label=r'$\lambda$1')
dSlm2dt = S.plume_spreads_dt(S_thresh)[:, 1]/1000**2
#ax.plot(S.times, dSlm2dt, label=r'$\lambda$2')
#ax.legend()
#ax.set_xlabel('Time [y]')
#ax.set_ylabel(r'Saturation plume spreading [km$^2$/y]')
#f.show()

# Print Run Summary (RSM) file
print('\nPrinting run summary (RSM) file')
O_c = O.plume_centroids(dP_thresh)
S_c = S.plume_centroids(S_thresh)
dum = np.vstack([np.array(O.times), np.array(PA), np.array(SA), np.array(dPAdt),
                 np.array(dSAdt), np.array(O_c[:, 0]), np.array(O_c[:, 1]),
                 np.array(S_c[:, 0]), np.array(S_c[:, 1]) ,np.array(dPxcsdt),
                 np.array(dPycsdt), np.array(dSxcsdt), np.array(dSycsdt),
                 np.array(dPlm1dt), np.array(dPlm2dt), np.array(dSlm1dt),
                 np.array(dSlm2dt), np.transpose(O.plume_angles(dP_thresh)),
                 np.transpose(S.plume_angles(S_thresh))])
dum = dum.transpose()
RSM = np.matrix(dum)
del dum
with open(output_file_name, 'w') as f:
    f.write(''.join([
        'Time[y], PA[km^2], SA[km^2], dPA/dt[km^2/y], dSA/dt[km^2/y], Pxc[m], ',
        'Pyc[m], Sxc[m], Syc[m], dPxcs/dt[m/y], dPycs/dt[m/y], dSxcs/dt[m/y], ',
        'dSycs/dt[m/y], dP1m1/dt[km^2/y], dP1m2/dt[km^2/y], dS1m1/dt[km^2/y], ',
        'dS1m2/dt[km^2/y], Pprim_angle[deg], Psec_angle[deg], Sprim_angle[deg], ',
        'Ssec_angle[deg]\n']))
    for line in RSM:
        np.savetxt(f, line, delimiter=',', fmt='%10.6e')

#############
# Animation #
#############
# This currently requires ffmpeg to be installed
# These wrap up all the metrics into one animation

O.plume_stability_animation(dP_thresh, filename='overpressure.mp4')
S.plume_stability_animation(S_thresh, filename='saturation.mp4')
