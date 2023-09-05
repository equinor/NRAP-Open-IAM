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
from re import sub
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.sep.join(['..',' ..', 'source']))
from openiam.mesh2D import Mesh2D


##############
# User input #
##############
datafolder = '/scratch/nts/dharp/projects/best/co2runs/igda06/run.1/prod1_rate'
dP_thresh = 1. # Overpressure threshold [MPa]
S_thresh = 0.01 # CO2 saturation threshold [m^3/m^3]
##################
# End user input #
##################

##################
# PREPARING DATA #
##################
# This is an example using FEHM TecPlot Rock Springs Uplift.
# The goal is to create a Mesh2D object and populate its variables with your data.
# THIS WILL NEED TO BE MODIFIED FOR DATA FILES IN OTHER FORMATS!!!

# Collect data files and make sure they are sorted by time
fls = sorted(glob(os.path.join(datafolder, 'run.*_sca_node.dat')),
             key=lambda x: int(x.split('.')[-2].split('_')[0]))
# Limit number of files for testing
#fls = fls[0:2]

# Read in mesh
with open('/scratch/fwo/tam/wy_2015/dharp_mesh_v2/tet.fehmn', 'r') as fh:
    fh.readline()
    ncrds = int(fh.readline().strip())
    g = np.genfromtxt(fh, names=['node', 'x_m', 'y_m', 'z_m'], max_rows=ncrds)

# Get unique x and y coordinates, need to round to nearest 10 meters due to stacked mesh offsets
xu = np.unique(np.round(g['x_m']/10)*10)
yu = np.unique(np.round(g['y_m']/10)*10)
# Create mesh grid
# xv,yv = np.meshgrid(xu,yu,sparse=True)
# Get cell horizontal area assuming that they are the same for all cells
# (in this case, they are all the same (200 m)).
# If they are all different, a numpy meshgrid of non-identical cell areas
# can be constructed and used to create the Mesh2D object below.
dxc = np.mean(np.diff(xu))
dyc = np.mean(np.diff(yu))
dAc = dxc*dyc # horizontal area of cells

# In this case, we read in 3D data and search for the max along each vertical column of cells
# at each timestep.
# Collect indices along x,y (i,j) columns, this way, it only has to be done once
col_inds = {}
for i, x in enumerate(xu):
    for j, y in enumerate(yu):
        col_inds[(i, j)] = np.intersect1d(
            np.where(np.abs(g['x_m']-x)<10)[0],
            np.where(np.abs(g['y_m']-y)<10)[0])

# Create Mesh2D object and add pressure co2 saturation variable types
M = Mesh2D(xu, yu, As=dAc)
O = M.add_variable('overpressure')
S = M.add_variable('co2sat')
# Mesh2D variables can be accessed by Python variables ("O" and "S")
# or by the Mesh2D variables dictionary ("M.variables['overpressure']
# and M.variables['co2sat']), whichever you prefer.

# Read in data at time 0 and add overpressure and saturation variables
# Assume saturations are all zero at time zero (overpressures are zero at time 0 by definition)
with open(fls[0], 'r') as fh:
    fh.readline()
    names = [sub('"', '', v) for v in fh.readline().strip().split('= "')[-1].split('" "')]
    t = float(fh.readline().split('time')[-1].split('days')[0])/365.25 # Convert to years
    d0 = np.genfromtxt(fh, names=names, max_rows=ncrds)
    #t = float(fls[0].split('_')[-1].split('yr')[0])
O.add_dataset(M.zeros_like_data(), t)
S.add_dataset(M.zeros_like_data(), t)
if 'Zone' in names:
    names.remove('Zone')

# Add data at each timesteps greater than zero
for fl in fls[1:]:
    with open(fl, 'r') as fh:
        t = float(fh.readline().split('time')[-1].split('days')[0])/365.25 # Convert to years
        dn = np.genfromtxt(fh, names=names[3:], max_rows=ncrds)
    # Collect max overpressure and saturation along each x,y column
    Pv = M.zeros_like_data()
    Sv = M.zeros_like_data()
    for xs, nds in col_inds.items():
        Pv[xs[0], xs[1]] = np.max(dn[nds]['Liquid_Pressure_MPa'] - d0[nds]['Liquid_Pressure_MPa'])
        Sv[xs[0], xs[1]] = np.max(dn[nds]['SuperCriticalLiquid_CO2_Saturation'])
    # Add datasets to Mesh2D variable
    O.add_dataset(Pv, t)
    S.add_dataset(Sv, t)

############
# Plotting #
############
# Now you can utilize the methods of the class to facilitate plotting plume stability metrics

# Plot plume areas
f,ax = plt.subplots(1)
PA = O.plume_areas(dP_thresh)/1000**2
ax.plot(O.times, PA)
ax.set_xlabel('Time [y]')
ax.set_ylabel(r'Pressure plume area [km$^2$]')
f.show()
f,ax = plt.subplots(1)
SA = S.plume_areas(S_thresh)/1000**2
ax.plot(S.times,SA)
ax.set_xlabel('Time [y]')
ax.set_ylabel(r'CO$_2$ area [km$^2$]')
f.show()

# Plot temporal derivative of plume area
f,ax = plt.subplots(1)
dPAdt = O.plume_areas_dt(dP_thresh)/1000**2
ax.plot(O.times, dPAdt, label='Pressure')
dSAdt = S.plume_areas_dt(S_thresh)/1000**2
ax.plot(S.times, dSAdt, label=r'CO$_2$')
ax.legend()
ax.set_xlabel('Time [y]')
ax.set_ylabel(r'dA/dt [km$^2$/y]')
f.show()

# Plot temporal derivative of plume centriod
f,ax = plt.subplots(1)
dPxcsdt = O.plume_centroids_dt(dP_thresh)[:, 0]
ax.plot(O.times, dPxcsdt, label='x direction')
dPycsdt = O.plume_centroids_dt(dP_thresh)[:, 1]
ax.plot(O.times, dPycsdt, label='y direction')
ax.legend()
ax.set_xlabel('Time [y]')
ax.set_ylabel(r'Pressure plume mobility [m/y]')
f.show()
f,ax = plt.subplots(1)
dSxcsdt = S.plume_centroids_dt(S_thresh)[:, 0]
ax.plot(S.times, dSxcsdt, label='x direction')
dSycsdt = S.plume_centroids_dt(S_thresh)[:, 1]
ax.plot(S.times, dSycsdt, label='y direction')
ax.legend()
ax.set_xlabel('Time [y]')
ax.set_ylabel(r'Saturation plume mobility [m/y]')
f.show()

# Plot temporal derivative of plume spread
f,ax = plt.subplots(1)
dPlm1dt = O.plume_spreads_dt(dP_thresh)[:, 0]/1000**2
ax.plot(O.times, dPlm1dt, label=r'$\lambda$1')
dPlm2dt = O.plume_spreads_dt(dP_thresh)[:, 1]/1000**2
ax.plot(O.times, dPlm2dt, label=r'$\lambda$2')
ax.legend()
ax.set_xlabel('Time [y]')
ax.set_ylabel(r'Pressure plume spreading [km$^2$/y]')
f.show()
f,ax = plt.subplots(1)
dSlm1dt = S.plume_spreads_dt(S_thresh)[:, 0]/1000**2
ax.plot(S.times, dSlm1dt, label=r'$\lambda$1')
dSlm2dt = S.plume_spreads_dt(S_thresh)[:, 1]/1000**2
ax.plot(S.times, dSlm2dt, label=r'$\lambda$2')
ax.legend()
ax.set_xlabel('Time [y]')
ax.set_ylabel(r'Saturation plume spreading [km$^2$/y]')
f.show()

# Print Run Summary (RSM) file
print('\nPrinting run summary (RSM) file')
dum = np.vstack([np.array(O.times), np.array(PA), np.array(SA), np.array(dPAdt),
               np.array(dSAdt), np.array(dPxcsdt), np.array(dPycsdt),
               np.array(dSxcsdt), np.array(dSycsdt), np.array(dPlm1dt),
               np.array(dPlm2dt), np.array(dSlm1dt), np.array(dSlm2dt)])
dum = dum.transpose()
RSM = np.matrix(dum)
del dum
with open('iam_plume_stability.RSM','w') as f:
    f.write(''.join([
        'Time[y], PA[km^2], SA[km^2], dPA/dt[km^2/y], dSA/dt[km^2/y], ',
        'dPxcs/dt[m/y], dPycs/dt[m/y], dSxcs/dt[m/y], dSycs/dt[m/y], ',
        'dP1m1/dt[km^2/y], dP1m2/dt[km^2/y], dS1m1/dt[km^2/y], dS1m2/dt[km^2/y]\n']))
    for line in RSM:
        np.savetxt(f, line, delimiter=',', fmt='%10.6e')

#############
# Animation #
#############
# This currently requires ffmpeg to be installed
# These wrap up all the metrics into one animation

O.plume_stability_animation(dP_thresh)
S.plume_stability_animation(S_thresh)

#################
# Interpolation #
#################
# Scipy efficiently interpolates on a numpy meshgrid
# The following uses the scipy functionality to plot interpolated values over time

# Choose some points to interpolate at
pts = [[9267, 6694],
       [8005, 7999],
       [9003, 8343]]

f,ax = plt.subplots(1)
ax.plot(O.times, O.interpolate(pts)[:, 0], label='Point 1')
ax.plot(O.times, O.interpolate(pts)[:, 1], label='Point 2')
ax.plot(O.times, O.interpolate(pts)[:, 2], label='Point 3')
ax.legend()
ax.set_xlabel('Time [y]')
ax.set_ylabel(r'Pressure [MPa]')
f.show()

f,ax = plt.subplots(1)
ax.plot(S.times, S.interpolate(pts)[:, 0], label='Point 1')
ax.plot(S.times, S.interpolate(pts)[:, 1], label='Point 2')
ax.plot(S.times, S.interpolate(pts)[:, 2], label='Point 3')
ax.legend()
ax.set_xlabel('Time [y]')
ax.set_ylabel(r'CO$_2$ saturation [-]')
f.show()
