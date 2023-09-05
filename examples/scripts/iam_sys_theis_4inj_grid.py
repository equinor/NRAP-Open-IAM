'''
Example illustrates system model containing four Theis reservoir models,
each with a different time-varying injection rate. System model is run
for an array of time points and an array of grid locations. Pressure changes
predicted by each reservoir model are added to calculate resulting pressure
at the observation location.

Examples of run:
$ python iam_sys_theis_4inj_grid.py
'''

import sys
import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, TheisReservoir
from openiam.enmesh import Grid

if __name__=='__main__':

    # Specify x, y coordinates of four wells
    wx = [739, 584, 287, 519]
    wy = [423, 822, 333, 856]
    num_wells = len(wx)

    # create site grid refined around injection wells
    grid = Grid(xmin=-2000, xmax=3000, ymin=-2000, ymax=3000, zmin=970, zmax=1000)
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
    rates.append(np.array([5.0e-3, 6.0e-3, 7.0e-3, 8.0e-3, 0, 0, 0, 0, 0, 0, 0, 0,]))
    rates.append(-rates[0])
    rates.append(2*rates[0])
    rates.append(-2*rates[0])

    time_array = np.array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.])

    def run_model_at_loc(x, y):

        # print('Solving at location', x, y)

        # Create system model
        sm_model_kwargs = {'time_array': time_array}  # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Create four Theis well component models, each predicting pressure at the same location
        tress = []
        for i in range(num_wells):

            # Add reservoir component
            tress.append(sm.add_component_model_object(
                TheisReservoir(name='tres'+str(i), parent=sm,
                               injX=wx[i], injY=wy[i], locX=x, locY=y,
                               injTimes=time_array, injRates=rates[i])))

            # Set input parameters
            tress[i].add_par('initialPressure', value=1.0e6)     # Pa
            tress[i].add_par('reservoirThickness', value=30)     # m
            tress[i].add_par('logResPerm', value=-10.69897)      # m^2
            tress[i].add_par('reservoirPorosity', value=.2)
            tress[i].add_par('brineDensity', value=1000)         # kg/m^3
            tress[i].add_par('brineViscosity', value=2.535e-3)   # Pa*s
            tress[i].add_par('CO2Density', value=800)            # kg/m^3
            tress[i].add_par('compressibility', value=2.46e-9)   # 1/Pa

            # Add observations of reservoir component model
            tress[i].add_obs('pressure')

        # Run system model using current values of its parameters
        sm.forward()

        # add pressure changes from all four wells
        pressures = sum(sm.collect_observations_as_time_series(tres, 'pressure') \
                        - 1.0e6 for i, tres in enumerate(tress)) + 1.0e6

        return pressures

    # Grid coordinates with z-axis index removed
    xx = X.squeeze()
    yy = Y.squeeze()

    # run model and collect results at all grid locations
    pressures = np.array([run_model_at_loc(x, y) for x, y in zip(
        np.ravel(xx), np.ravel(yy))]).reshape(xx.shape[0], xx.shape[1], len(time_array))

    # Check results near center of grid
    print('Results at location', xx[20, 22], yy[20, 22])
    for j, time in enumerate(time_array):
        print(rates[0][j], rates[1][j], rates[2][j], rates[3][j], pressures[20, 22, j])

    # Pressure change in MPa
    dP = (pressures - 1.0e6)/1.0e6

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
