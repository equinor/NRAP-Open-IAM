'''
NRAP-Open-IAM Time To First Detection example
This examples couples the simple reservoir, multisegmented wellbore and
carbonate aquifer models.  The saturation/pressure output produced by simple
reservoir model is used to drive leakage from a single open wellbore
model, which is passed to the input of an adapter that provides well
coordinates, CO2 and brine leakage rates and cumulative mass fluxes to the
carbonate aquifer model.  A matrix of time to first detection in the aquifer for
each realization are output.

DREAM accepts text files with an .iam extension as input.
The format should be:

IAM, Scenario#, monitoring parameter, threshold indicator, threshold value
x1, y1 , z1 , time to detection 1
.
.
.
xn, yn , zn, time to detection n

Here, the threshold indicator should be “above”, “below”, “relative”, “absolute”.
For example, if we want to identify TDS > 1200 mg/L as a leak, it would be:
    IAM, Scenario#, TDS, above, 1200

If we want to identify pH < 7 as a leak, it would be:
    IAM, Scenario#, pH, below, 7

If we want to identify a pressure change > 10,000 psi as a leak, it would be:
    IAM, Scenario#, aqueous pressure, absolute, 10000

you can specify a negative 10% change as:
    IAM, Scenario#, pH, relative, -0.1
or a positive 10% change as
    IAM, Scenario#, pH, relative, +0.1

Not including a sign will specify that a change in either direction triggers a leak.

Usage examples:
$ python iam_ttfd_example.py

'''
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, SimpleReservoir, MultisegmentedWellbore,
                     CarbonateAquifer, RateToMassAdapter)


if __name__ == "__main__":
    # For multiprocessing in Spyder
    __spec__ = None
    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    sres = sm.add_component_model_object(SimpleReservoir(
        name='sres', parent=sm, injX=0., injY=0., locX=500., locY=500.))

    # Add parameters of reservoir component model
    sres.add_par('numberOfShaleLayers', value=3, vary=False)
    sres.add_par('shale1Thickness', min=30.0, max=150., value=45.0)
    sres.add_par('injRate', min=0.1, max=10.0, value=1.0)

    # Add observations of reservoir component model
    sres.add_obs_to_be_linked('pressure')
    sres.add_obs_to_be_linked('CO2saturation')
    sres.add_obs('pressure')
    sres.add_obs('CO2saturation')
    sres.add_obs('mass_CO2_reservoir')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))
    ms.add_par('wellRadius', min=0.025, max=0.075, value=0.05)
    ms.add_par('logWellPerm', min=-14.0, max=-13.0, value=-13.5)

    # Add linked parameters: common to both components
    ms.add_par_linked_to_par('numberOfShaleLayers',
                             sres.deterministic_pars['numberOfShaleLayers'])
    ms.add_par_linked_to_par('shale1Thickness', sres.pars['shale1Thickness'])
    ms.add_par_linked_to_par('shale2Thickness',
                             sres.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('shale3Thickness',
                             sres.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('aquifer1Thickness',
                             sres.default_pars['aquiferThickness'])
    ms.add_par_linked_to_par('aquifer2Thickness',
                             sres.default_pars['aquiferThickness'])
    ms.add_par_linked_to_par('reservoirThickness',
                             sres.default_pars['reservoirThickness'])
    ms.add_par_linked_to_par('datumPressure',
                             sres.default_pars['datumPressure'])

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

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
    adapt = sm.add_component_model_object(RateToMassAdapter(name='adapt', parent=sm))
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

    # add carbonate aquifer model object and define parameters
    ca = sm.add_component_model_object(CarbonateAquifer(name='ca', parent=sm))
    ca.add_par('perm_var', min=0.017, max=1.89, value=0.9535)
    ca.add_par('corr_len', min=1.0, max=3.95, value=2.475)
    ca.add_par('aniso', min=1.1, max=49.1, value=25.1)
    ca.add_par('mean_perm', min=-13.8, max=-10.3, value=-12.05)
    ca.add_par('aqu_thick', value=300., vary=False)
    ca.add_par('hyd_grad', min=2.88e-4, max=1.89e-2, value=9.59e-3)
    ca.add_par('cl', min=0.1, max=6.025, value=0.776)
    ca.model_kwargs['x'] = [0.]
    ca.model_kwargs['y'] = [0.]

    CO2_rate_obs_list = []
    brine_rate_obs_list = []
    CO2_mass_obs_list = []
    brine_mass_obs_list = []
    CO2_rate_obs_list.append(ms.linkobs['CO2_aquifer1'])
    brine_rate_obs_list.append(ms.linkobs['brine_aquifer1'])
    CO2_mass_obs_list.append(adapt.linkobs['mass_CO2_aquifer1'])
    brine_mass_obs_list.append(adapt.linkobs['mass_brine_aquifer1'])

    # Add aquifer component's keyword argument co2_rate linked to the collection created above
    ca.add_kwarg_linked_to_collection('co2_rate', CO2_rate_obs_list)
    ca.add_kwarg_linked_to_collection('brine_rate', brine_rate_obs_list)
    ca.add_kwarg_linked_to_collection('co2_mass', CO2_mass_obs_list)
    ca.add_kwarg_linked_to_collection('brine_mass', brine_mass_obs_list)

    # Add observations (output) from the carbonate aquifer model
    ca.add_obs('dx')
    ca.add_obs('dy')

    # no-impact or MCL threshold value
    if ca.default_pars['ithresh'].value == 1:
        threshold = 6.5
    else:
        threshold = 6.6

    output_directory = os.sep.join(['..', '..', 'Output', 'dream_output'])
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Define grid
    nx = 200
    ny = 200
    nz = 10

    xmin = 0
    xmax = 1000
    ymin = 0
    ymax = 1000
    zmin = 0
    zmax = -ca.deterministic_pars['aqu_thick'].value

    x = np.linspace(xmin, xmax, nx)  # x values of interest
    y = np.linspace(ymin, ymax, ny)[:, None]  # y values of interest, as a "column" array
    z = np.linspace(zmin, zmax, nz)

    print('------------------------------------------------------------------')
    print('                          UQ illustration ')
    print('------------------------------------------------------------------')

    import random
    num_samples = 25
    ncpus = 2
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=random.randint(500, 1100))
    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)

    # get plumes
    x0 = sres.locX
    y0 = sres.locY
    ind_list = list(range(len(time_array)))
    dxs = np.array([s.recarray['ca.dx_'+str(indd)] for indd in ind_list])
    dys = np.array([s.recarray['ca.dy_'+str(indd)] for indd in ind_list])

    # for each sample
    for j in range(num_samples):
        print('sample', j+1)
        dx = dxs[:, j]
        dy = dys[:, j]

        # Calculate time to first detection
        ttfd = np.ones((len(x), len(y)))*1.e30
        for tt, a, b in zip(time_array, dx, dy):
            if (a > 0 and b > 0):
                # True for points inside the ellipse
                ellipse = ((x-x0)/(a*0.5))**2 + ((y-y0)/(b*0.5))**2 <= 1
                mask = np.logical_not(ellipse)*1.0e30
                plume = ellipse * tt
                ttfd = np.minimum(ttfd, mask + plume)

        # Write output for DREAM
        filename = os.sep.join([output_directory, 'ttfd_' + str(j+1) + '.iam'])
        f = open(filename, "w+")
        metric = 'pH'
        indicator = 'below'
        f.write('IAM,' + str(j+1) + ',' + metric + ',' + indicator + ',' + \
                str(threshold) + ',\n')
        ix = -1
        for xx in x:
            ix = ix + 1
            iy = -1
            for yy in y:
                iy = iy + 1
                for iz in z:
                    f.write(str(xx) + ',' + str(yy[0]) + ',' + str(iz) + ',' + \
                            str(ttfd[ix][iy]) + '\n')
        f.close()

        # Make plots
        X, Y = np.meshgrid(x, y)
        cp = plt.contourf(X, Y, ttfd, levels=time_array)
        plt.colorbar(cp)
        plt.title('Realization ' + str(j+1) + '\nTime to Detection (days)')
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        # plt.rcParams.update({'font.size': 18})
        filename = os.sep.join([output_directory, 'ttfd_' + str(j+1) + '.png'])
        plt.savefig(filename, dpi=300)
        plt.gcf().clear()
