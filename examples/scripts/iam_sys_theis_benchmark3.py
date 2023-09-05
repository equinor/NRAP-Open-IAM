'''
Example illustrates system model containing four Theis reservoir models,
each with a different time-varying injection rate.  System model is run
for an array of time points and an array of grid locations. Pressure changes
predicted by each reservoir model are added to calculate resulting pressure
at the observation location.

Examples of run:
$ python iam_sys_theis_benchmark3.py
'''

import sys
import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
try:
    from tqdm import tqdm
except ModuleNotFoundError:
    err_msg = "Your environment is missing tqdm library. Please update and try again."
    raise ModuleNotFoundError(err_msg)

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, TheisReservoir
from openiam.enmesh import Grid

if __name__=='__main__':

    # Benchmark data
    time_array = np.array([
        0.,   6.,   12.,  18.,  24.,  30.,  36.,  42.,  48.,  54.,  60.,
        66.,  72.,  78.,  84.,  90.,  96.,  102., 108., 114., 120., 126.,
        132., 138., 144., 150., 156., 162., 168., 174., 180., 186., 192.,
        198., 204., 210., 216., 222., 228., 234., 240.])/24.
    stomp_pressures = np.array([
        600000., 599327., 599311., 599310., 599311., 599176., 599173.,
        599173., 599173., 599039., 599035., 599035., 599035., 598901.,
        598897., 598897., 598897., 599974., 600001., 600001., 600000.,
        600000., 600000., 600000., 600000., 600000., 600000., 600000.,
        600000., 600000., 600000., 600000., 600000., 600000., 600000.,
        600000., 600000., 600000., 600000., 600000., 600000.])
    rate = np.array([
        0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.006, 0.006, 0.007,
        0.007, 0.007, 0.007, 0.008, 0.008, 0.008, 0.008, 0.   , 0.   ,
        0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
        0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
        0.   , 0.   , 0.   , 0.   , 0.   ])

    # Specify x, y coordinates of four injection wells
    wx = [739, 584, 287, 519]
    wy = [423, 822, 333, 856]

    # Specify x, y coordinates of observation well
    obsx = 532
    obsy = 608

    num_wells = len(wx)

    # create site grid refined around injection wells
    grid = Grid(xmin=-2000, xmax=3000, ymin=-2000, ymax=3000, zmin=0, zmax=50)
    grid.add_well(obsx, obsy, name='Obs')
    for i in range(num_wells):
        grid.add_well(wx[i], wy[i], name=str(i))
    grid.refine_grid_around_wells(min_x_spacing=5, min_y_spacing=5, x_mul=1.5, y_mul=1.5)

    # print out grid and well locations for STOMP
    print('Grid card for STOMP')
    grid.print_stomp()
    print('Well grid indices for STOMP')
    locs = grid.get_well_grid_locations()
    for loc in locs:
        print(loc)

    X, Y, Z = grid.get_centroids()

    # Specify injection rates over time for four wells
    rates = []
    rates.append(rate)
    rates.append(-rate)
    rates.append(2*rate)
    rates.append(-2*rate)

    initialPressure = 600000.

    def run_model_at_loc(x, y):

        # Create system model
        sm_model_kwargs = {'time_array': time_array}  # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Create four Theis well component models, each predicting pressure at the same locations
        tress = []
        for i in range(num_wells):

            # Add reservoir component
            tress.append(sm.add_component_model_object(
                TheisReservoir(name='tres'+str(i), parent=sm,
                               injX=wx[i], injY=wy[i],
                               locX=x, locY=y,
                               injTimes=time_array, injRates=rates[i])))

            # Set input parameters
            tress[i].add_par('initialPressure', value=initialPressure) # Pa
            tress[i].add_par('reservoirThickness', value=50)           # m
            tress[i].add_par('logResPerm', value=-10.6271524099)       # m^2
            tress[i].add_par('reservoirPorosity', value=.35)
            tress[i].add_par('brineDensity', value=9.98664E+02)        # kg/m^3
            tress[i].add_par('brineViscosity', value=1.01759E-03)      # Pa*s
            # Assume we're injecting water, not CO2
            tress[i].add_par('CO2Density', value=9.98664E+02)          # kg/m^3
            tress[i].add_par('compressibility', value=4.3531E-11)      # 1/Pa

            # Add observations of reservoir component model
            tress[i].add_obs('pressure')

        # Run system model using current values of its parameters
        sm.forward()

        # add pressure changes from all four wells
        pressures = sum(sm.collect_observations_as_time_series(tres, 'pressure')\
                        - initialPressure for i, tres in enumerate(tress)) + initialPressure

        return pressures

    # Get grid coordinates
    xx = X.squeeze()  # remove z-axis
    yy = Y.squeeze()  # remove z-axis
    xy_coordinates = zip(np.ravel(xx), np.ravel(yy))

    # run model and collect results
    pressures = np.array([run_model_at_loc(x, y) for x, y in tqdm(
        xy_coordinates, total=xx.size)]).reshape(
            xx.shape[0], yy.shape[1], len(time_array))

    # Pressure change in MPa
    dP = (pressures - initialPressure)/1.0e6

    # Animated contour plot of pressure change over time
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlabel('E-W Distance, m')
    ax.set_ylabel('N-S Distance, m')
    CS = ax.contourf(xx[:, :], yy[:, :], dP[:, :, 0],
                     levels=np.linspace(dP.min(), dP.max(), 7))
    cbar = fig.colorbar(CS)
    cbar.set_label('Pressure Change, MPa')
    ax.scatter(wx, wy, c='w', marker='.')
    for i, well in enumerate(grid.wells):
        ax.annotate(well.name, (well.x, well.y), c='w')

    def contour_plot(j):
        ax.clear()
        _ = ax.contourf(xx[:, :], yy[:, :], dP[:, :, j],
                         levels=np.linspace(dP.min(), dP.max(), 10))
        ax.scatter(wx, wy, c='w', marker='.')
        for well in grid.wells:
            ax.annotate(well.name, (well.x, well.y), c='w')
        ax.set_title(str(time_array[j])+' days')

    ani = FuncAnimation(fig, contour_plot, len(time_array), interval=1000, blit=False)

    plt.show()

    # Compare pressures predicted by Theis solution and STOMP-W at observation well
    pressures = run_model_at_loc(obsx, obsy)
    fig, ax = plt.subplots()
    ax.plot(time_array, pressures, color='red', label='Theis')
    ax.scatter(time_array, stomp_pressures, color='blue', label='STOMP-W')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Pressure (Pa)')
    ax.legend()
    plt.tight_layout()
    plt.show()
