# @author: Jaisree Iyer
# iyer5@llnl.gov

"""
System model contains two linked component models: reservoir and chemical
well sealing models. The location of the wellbores are random. Some paremeters
of the second component are obtained from the observations of the first component.

Examples of run:
$ python iam_sys_reservoir_chem_well_seal.py
"""

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import (SystemModel, SimpleReservoir, ChemicalWellSealing)
from matk import pyDOE


if __name__ == '__main__':
    # For multiprocessing in Spyder
    __spec__ = None

    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}   # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Define keyword arguments of the second system model and create it
    time_array = 365.25*np.arange(0.0, 2.0)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    sm2 = SystemModel(model_kwargs=sm_model_kwargs)

    # Create randomly located leaky well locations within a box
    # The entries in wellPropsMins and wellPropMaxs are the minimum and maximum
    # limits on x- and y-coordinates of the location, fracture aperture, and
    # fracture length (all in m)
    wellPropMins = np.array([-100., -100., 0.00008, 10.])
    wellPropMaxs = np.array([100., 100., 0.0003, 100.])
    # Sample the locations and properties of wells
    num_wells = 200
    well_props = wellPropMins + pyDOE.lhs(4, samples=num_wells)*(
        wellPropMaxs-wellPropMins)

    # Initialize the simple reservoir and the chemical well sealing models
    sress = []
    cwss = []
    for i, crds in enumerate(well_props):

        # Add reservoir components with an injector at 0,0 and desired well locations
        sress.append(sm.add_component_model_object(
            SimpleReservoir(name='sres'+str(i), parent=sm,
                            injX=0., injY=0., locX=crds[0], locY=crds[1])))

        # Add parameters of reservoir component model
        sress[-1].add_par('numberOfShaleLayers', value=3, vary=False)
        sress[-1].add_par('injRate', value=0.04, vary=False)
        sress[-1].add_par('shale1Thickness', min=30.0, max=50., value=40.0)
        sress[-1].add_par('aquifer1Thickness', min=20.0, max=60., value=50.0)
        sress[-1].add_par('aquifer2Thickness', min=30.0, max=50., value=45.0)
        sress[-1].add_par('logResPerm', min=-13., max=-11., value=-12.5)

        # Add observations of reservoir component model
        sress[-1].add_obs('pressure')

    # Run the reservoir model
    sm.forward()
    # Initialize the maximum overpressure input to the chemical sealing wellbore ROM
    maxOverpressure_values = []

    # Create the chemical well sealing component for each well
    for i, sres in enumerate(sress):
        pressure = sm.collect_observations_as_time_series(sres, 'pressure')
        # maxOverpressure = maximum increase in pressure experienced by the well
        maxOverpressure_value = np.amax(pressure)-pressure[0]
        maxOverpressure_values.append((np.amax(pressure)-pressure[0])/1000000)
        # Add chemicalWellSealing ROM components
        cwss.append(sm2.add_component_model_object(ChemicalWellSealing(
            name='cws'+str(i), parent=sm2)))

        # Add parameters of chemical well sealing component
        cwss[-1].add_par('fractureAperture', value=well_props[i][2])
        cwss[-1].add_par('fractureLength', value=well_props[i][3])
        cwss[-1].add_par('maxOverpressure', value=maxOverpressure_value)
        cwss[-1].add_obs('seal_flag', index=[0])
        cwss[-1].add_obs('seal_time', index=[0])

    # Run the chemical well sealing component
    sm2.forward()

    # Initialize variables needed for plotting and printing
    seal_flags = []
    seal_times = []
    minimumBrineResidenceTimes = []
    viscosity = 5.154*0.0001 # viscosity of 1m NaCl solution at 60 C
    for i, cws in enumerate(cwss):
        # Print inputs for each well
        print('Well location in m (x):', well_props[i][0])
        print('Well location in m (y):', well_props[i][1])
        print('fractureLength in m:', well_props[i][3])
        print('fractureAperture in micron:', (well_props[i][2])*1000000)
        print('maxOverpressure in MPa:', maxOverpressure_values[i])

        # Set aperture to 10 microns if it is lower than that
        aperture = well_props[i][2]*1000000
        if well_props[i][2]*1000000 < 10:
            aperture = 10
        # Compute brine residence times
        maximumBrineVelocity = (aperture*10**(-6))**2*maxOverpressure_values[i]*\
            10**6/12/viscosity/well_props[i][3]
        minimumBrineResidenceTimes.append(well_props[i][3]/maximumBrineVelocity/60)
        # Print outputs for each well
        seal_flags.append(sm2.obs['cws'+str(i)+'.seal_flag_0'].sim)
        seal_times.append(sm2.obs['cws'+str(i)+'.seal_time_0'].sim/86400)
        print('seal_flag:', sm2.obs['cws'+str(i)+'.seal_flag_0'].sim)
        print('seal_time:', sm2.obs['cws'+str(i)+'.seal_time_0'].sim/86400)

    # Plot results
    # Setup plot parameters
    line_width = 1
    xy_label_size = 14
    title_size = 18
    fig = plt.figure(figsize=(10, 8))

    # plot aperture distribution
    ax = fig.add_subplot(221)
    sc = plt.scatter(well_props[:, 0], well_props[:, 1], c=well_props[:, 2]*1000000)
    plt.colorbar(sc)
    plt.xlabel('Location along x (m)', fontsize=xy_label_size)
    plt.ylabel('Location along y (m)', fontsize=xy_label_size)
    plt.title('Aperture distribution of wells (microns)', fontsize=title_size)
    plt.tight_layout()
    plt.tick_params(labelsize=12)

    # plot fracture length distribution
    ax = fig.add_subplot(222)
    sc = plt.scatter(well_props[:, 0], well_props[:, 1], c=well_props[:, 3])
    plt.colorbar(sc)
    plt.xlabel('Location along x (m)', fontsize=xy_label_size)
    plt.ylabel('Location along y (m)', fontsize=xy_label_size)
    plt.title('Cement coverage of wells (m)', fontsize=title_size)
    plt.tight_layout()
    plt.tick_params(labelsize=12)

    # plot maximum overpressure distribution
    ax = fig.add_subplot(223)
    sc = plt.scatter(well_props[:, 0], well_props[:, 1], c=maxOverpressure_values)
    plt.colorbar(sc)
    plt.xlabel('Location along x (m)', fontsize=xy_label_size)
    plt.ylabel('Location along y (m)', fontsize=xy_label_size)
    plt.title('Maximium overpressure experienced by wells (MPa)', fontsize=title_size)
    plt.tight_layout()
    plt.tick_params(labelsize=12)

    # plot which wells are predicted to seal
    ax = fig.add_subplot(224)
    sc = plt.scatter(well_props[:, 0], well_props[:, 1], c=seal_flags)
    plt.colorbar(sc)
    plt.xlabel('Location along x (m)', fontsize=xy_label_size)
    plt.ylabel('Location along y (m)', fontsize=xy_label_size)
    plt.title('Sealing wells', fontsize=title_size)
    plt.tight_layout()
    plt.tick_params(labelsize=12)

    to_save = True
    if to_save:
        plt.savefig('sealingWells.png', dpi=300)

    # plot sealing map for QA-QC
    fig2 = plt.figure(figsize=(4, 4))
    sc = plt.scatter(well_props[:, 2]*1000000, minimumBrineResidenceTimes, c=seal_flags)
    plt.plot([10, 50, 100, 250, 500, 1000, 2000],
             [0.04, 0.8, 4, 16, 80, 200, 1000],
             color="black", linewidth=line_width)
    plt.colorbar(sc)
    plt.xlabel('Fracture aperture (micron)', fontsize=xy_label_size)
    plt.ylabel('Residence time (minutes)', fontsize=xy_label_size)
    plt.title('Sealing map', fontsize=title_size)
    plt.tight_layout()
    plt.tick_params(labelsize=12)
    plt.xlim([80, 300])
    plt.ylim([0.01, 100.0])
    plt.yscale('log')
    plt.xscale('log')

    to_save = True
    if to_save:
        plt.savefig('sealingMap.png', dpi=300)

    # plot sealing time for QA-QC
    fig3 = plt.figure(figsize=(4, 4))
    sc = plt.scatter(well_props[:, 2]*1000000, seal_times, c=seal_flags)
    plt.plot(
        [0, 50, 100, 150],
        [9.54*10**(-3)*0**1.79, 9.54*10**(-3)*50**1.79,
         9.54*10**(-3)*100**1.79, 9.54*10**(-3)*150**1.79],
        color="black", linewidth=line_width)
    plt.xlabel('Fracture aperture (micron)', fontsize=xy_label_size)
    plt.ylabel('Time (days)', fontsize=xy_label_size)
    plt.title('Sealing time', fontsize=title_size)
    plt.tight_layout()
    plt.tick_params(labelsize=12)
    plt.xlim([0, 150])
    plt.ylim([10, 60.0])

    to_save = True
    if to_save:
        plt.savefig('sealingTime.png', dpi=300)
