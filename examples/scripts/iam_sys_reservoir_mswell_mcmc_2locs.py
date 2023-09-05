'''
This example performs MCMC on the simple reservoir model to produce
posterior distributions for the reservoir permeability (logResPerm),
reservoir porosity (reservoirPorosity) and reservoir thickness
(reservoirThickness) - parameters constrained by measured values of pressure
and |CO2| saturation at monitoring well. The posterior parameter values are
then run through the multisegmented wellbore model in order to propagate
the uncertainty to the brine and |CO2| leakage at leaking well.
Permeability of the leaking well is defined as a stochastic parameter.

Example of run:
$ python iam_sys_reservoir_mswell_mcmc_2locs.py
'''
import sys
import os
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import SystemModel, SimpleReservoir, MultisegmentedWellbore
from matk import logposteriorwithvariance


if __name__=='__main__':
    # Define keyword arguments of the system model
    param_dist = 1  # 1 for Gaussian distribution; 2 for Uniform distribution
    nTypeData = 2  # The number of measurement type
    num_years = 5  # Total time for acquiring data
    DataFreq = 12  # The number of measurements per year
    locX_monitor = 50.0  # X-coordinate for monitoring well
    locY_monitor = 50.0  # Y-coordinate for monitoring well
    # X-coordinate for the target well that the uncertainty of wellbore leakage will be quantified
    locX_leakwell = 200.0
    # Y-coordinate for the target well that the uncertainty of wellbore leakage will be quantified
    locY_leakwell = 200.0
    # Measurement or data assimilation error tolerance for pressure data (137895 pa = 20psi)
    Pres_err = 137895
    # Measurement or data assimilation error tolerance for saturation data
    Sat_err = 0.02

    num_timeSeries = num_years*DataFreq
    time_array = 365.25/DataFreq*np.arange(0.0, num_timeSeries+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    nwalkers = 20  # Number of MCMC chains
    nsamples_per_walker = 200  # Number of samples (exclusive burnin samples) per chain (walker)
    burnin = 50  # Number of samples to throw away at the beginning of chains

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component model
    sres = sm.add_component_model_object(
        SimpleReservoir(name='sres', parent=sm, locX=locX_monitor, locY=locY_monitor))

    # Add parameters of reservoir component model
    sres.add_par('logResPerm', min=-13., max=-11., value=-12.)
    sres.add_par('reservoirPorosity', min=0.05, max=0.35, value=0.2)
    sres.add_par('reservoirThickness', min=20., max=40., value=30)

    # Add observations of reservoir component model to be used by the next component
    sres.add_obs('pressure')  # by default observations at all time points are added
    sres.add_obs('CO2saturation')

    # Create synthetic observations by running the model forward
    # with specified parameter 'values' defined above in add_par calls
    # Observation data will be read in from files in the future
    sm.forward()
    pres_list = []
    sat_list = []
    for ind in range(len(time_array)):
        pres_list.append(sres.obs['pressure'+'_'+str(ind)].sim)
        sat_list.append(sres.obs['CO2saturation'+'_'+str(ind)].sim)

    # Add noise (measurement error) to data
    nData = nTypeData*len(time_array)
    normNum = np.random.normal(0, 1, nData)
    presNoise = Pres_err*normNum[0:len(time_array)]
    satNoise = Sat_err*normNum[len(time_array):nData]
    pres_list = pres_list + presNoise
    sat_list = sat_list + satNoise

    sres.add_obs('pressure', index=list(range(num_timeSeries+1)), value=pres_list)
    sres.add_obs('CO2saturation', index=list(range(num_timeSeries+1)), value=sat_list)
    sres.add_obs_to_be_linked('pressure')
    sres.add_obs_to_be_linked('CO2saturation')

    # Run MCMC to produce parameter posterior distributions based on consistency
    # with sres model observations 'pressure' and 'CO2saturation'
    # 137895 Pa and 0.02, respectively, are the history matching
    # error tolerances for pressure and saturation data. 137895 Pa = 20 psi
    mylnprob = logposteriorwithvariance(
        sm, var=[Pres_err]*len(time_array) + [Sat_err]*len(time_array))
    ss = sm.emcee(nwalkers=nwalkers, lnprob=mylnprob,
                  nsamples=nsamples_per_walker, burnin=burnin)
    print('MCMC sampling posterior: Done!')

    # Check the convergence of MCMC chain
    # Calculate the mean and variance of the MCMC samples
    nMCMCSamples = nsamples_per_walker*nwalkers
    step_size = nwalkers*2
    step_num = int(nMCMCSamples/step_size)
    Mean_chain = np.zeros((3, step_num))
    var_chain = np.zeros((3, step_num))

    xx = []
    for i in range(0, step_num):
        Mean_samples0 = np.mean(ss.chain[:, 0:step_size*(i+1), 0])
        var_samples0 = np.var(ss.chain[:, 0:step_size*(i+1), 0])
        Mean_samples1 = np.mean(ss.chain[:, 0:step_size*(i+1), 1])
        var_samples1 = np.var(ss.chain[:, 0:step_size*(i+1), 1])
        Mean_samples2 = np.mean(ss.chain[:, 0:step_size*(i+1), 2])
        var_samples2 = np.var(ss.chain[:, 0:step_size*(i+1), 2])
        val = (i+1)*step_size
        xx.append(val)
        Mean_chain[0][i] = Mean_samples0
        var_chain[0][i] = var_samples0
        Mean_chain[1][i] = Mean_samples1
        var_chain[1][i] = var_samples1
        Mean_chain[2][i] = Mean_samples2
        var_chain[2][i] = var_samples2

    f1, ax = plt.subplots(3, 1, sharex=True, figsize=(10, 10))
    ax[0].plot(xx, Mean_chain[0], color='k', alpha=0.4)
    ax[0].set_ylabel('logResPerm')
    ax[0].set_ylim(-13, -11)
    ax[1].plot(xx, Mean_chain[1], color='k', alpha=0.4)
    ax[1].set_ylabel('Porosity')
    ax[1].set_ylim(0.05, 0.35)
    ax[2].plot(xx, Mean_chain[2], color='k', alpha=0.4)
    ax[2].set_xlabel('Number of MCMC Samples')
    ax[2].set_ylabel('reservoirThickness')
    ax[2].set_ylim(20., 40.)
    f1.savefig('chain_mean.png')
    f1.show()

    f2,ax = plt.subplots(3, 1, sharex=True, figsize=(10, 10))
    ax[0].plot(xx, var_chain[0], color='k', alpha=0.4)
    ax[0].set_ylabel('logResPerm')
    ax[1].plot(xx, var_chain[1], color='k', alpha=0.4)
    ax[1].set_ylabel('Porosity')
    ax[2].plot(xx, var_chain[2], color='k', alpha=0.4)
    ax[2].set_xlabel('Number of MCMC Samples')
    ax[2].set_ylabel('reservoirThickness')
    f2.savefig('chain_var.png')
    f2.show()

    # Create trace plots. This is done somewhat manually here, but could
    # be easily encapsulated in a SystemModel method
    # Perm plots
    f3,ax = plt.subplots(3, 1, sharex=True, figsize=(10, 10))
    ax[0].plot(ss.chain[:, :, 0].T, color='k', alpha=0.4)
    ax[0].axhline(-12, color='gray', lw=2)
    ax[0].set_ylabel('logResPerm')
    ax[0].set_ylim(-13, -11)
    # Porosity plots
    ax[1].plot(ss.chain[:, :, 1].T, color='k', alpha=0.4)
    ax[1].axhline(0.2, color='gray', lw=2)
    ax[1].set_ylabel('Porosity')
    ax[1].set_ylim(0.05, 0.35)
    # Thickness plots
    ax[2].plot(ss.chain[:, :, 2].T, color='k', alpha=0.4)
    ax[2].axhline(30, color='gray', lw=2)
    ax[2].set_xlabel('Step number')
    ax[2].set_ylabel('reservoirThickness')
    ax[2].set_ylim(20., 40.)
    f3.show()
    f3.savefig('mcmc_traces.png')

    # Create panels plot to explore the histograms and correlations of
    # posterior parameter distributions
    sc = sm.create_sampleset(ss.chain[:, :, :].reshape((-1, len(sm.pars))))
    sc.panels(bins=50, figsize=(10, 10), fontsize=16, ms=1,
              filename='mcmc_posteriors_panels.png')

    # Create sample set including wellbore permeability
    mcmc_pars = ss.chain[:, :, :].reshape((-1, len(sm.pars)))
    if param_dist == 1:  # Guassian distribution for wellbore perm
        wbperms = np.random.normal(-13, 0.15, len(mcmc_pars))
    elif param_dist == 2:  # Uniform distribution for wellbore perm
        wbperm_min = -14
        wbperm_max = -12
        wbperms = wbperm_min+(wbperm_max-wbperm_min)*np.random.random(len(mcmc_pars))
    else:
        print('There is no such option for param_dist, please reset!')
    sc_pars = np.column_stack([mcmc_pars, wbperms])

    num_years = 10
    time_array1 = 365.25*np.arange(0.0, num_years+1, step=1./4.)
    sm2_model_kwargs = {'time_array': time_array1} # time is given in days

    # Create system model
    sm2 = SystemModel(model_kwargs=sm2_model_kwargs)

    # Add reservoir component to system model
    sres = sm2.add_component_model_object(
        SimpleReservoir(name='sres', parent=sm2, locX=locX_leakwell, locY=locY_leakwell))

    # Add parameters of reservoir component model
    if param_dist == 1:  # Guassian distribution for wellbore perm
        sres.add_par('logResPerm', dist='norm', mean=-12, std=0.33)
        sres.add_par('reservoirPorosity', dist='norm', mean=0.2, std=0.033)
        sres.add_par('reservoirThickness', dist='norm', mean=30, std=3.3)
    elif param_dist==2:  # Uniform distribution for wellbore perm
        sres.add_par('logResPerm', min=-13., max=-11.)
        sres.add_par('reservoirPorosity', min=0.05, max=0.35)
        sres.add_par('reservoirThickness', min=20., max=40.)
    else:
        print('There is no such option for param_dist, please reset!')

    # Add observations of reservoir component model to be used by the next component
    sres.add_obs_to_be_linked('pressure')
    sres.add_obs_to_be_linked('CO2saturation')

    # Add multisegmented wellbore component so that reservoir uncertainty
    # can be propagated through a "default" wellbore to estimate the uncertainty
    # in CO2 and brine leakage.
    # Note that the uncertainty in wellbore permeability is also considered.
    ms = sm2.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm2))

    if param_dist == 1:  # Guassian distribution for wellbore perm
        ms.add_par('logWellPerm', dist='norm', mean=-13, std=0.33)
    elif param_dist==2:  # Uniform distribution for wellbore perm
        ms.add_par('logWellPerm', min=-14, max=-12)
    else:
        print('There is no such option for param_dist, please reset!')

    # Add parameters of wellbore component linked to the parameters
    # of the same name from reservoir model. In this example most of wellbore
    # component parameters are default parameters, i.e., linked to the
    # corresponding parameters of the reservoir model
    ms.add_par_linked_to_par('numberOfShaleLayers', sres.default_pars['numberOfShaleLayers'])
    ms.add_par_linked_to_par('shale1Thickness', sres.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('shale2Thickness', sres.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('shale3Thickness', sres.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('aquifer1Thickness', sres.default_pars['aquiferThickness'])
    ms.add_par_linked_to_par('aquifer2Thickness', sres.default_pars['aquiferThickness'])
    ms.add_par_linked_to_par('reservoirThickness', sres.pars['reservoirThickness'])
    ms.add_par_linked_to_par('datumPressure', sres.default_pars['datumPressure'])

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])
    ms.add_obs('brine_aquifer1')
    ms.add_obs('CO2_aquifer1')

    # Run posterior parameter distributions with multisegmented wellbore model
    # added to the system model
    sc = sm2.create_sampleset(sc_pars)
    sc.run(cpus=4, verbose=False)

    # Run prior (uniform) distribution through coupled model to determine
    # prior uncertainty in leakage. The number samples equals the number of
    # mcmc samples minus the burnin for each chain
    sprior = sm2.lhs(siz=nwalkers*nsamples_per_walker)
    sprior.run(cpus=2,verbose=False)

    # Plot resulting CO2 and brine leakage histograms to explore
    # the resulting leakage uncertainty due to reservoir permeability,
    # porosity, and thickness uncertainty
    obs_ind = '_20'   # specify the index of time point of interest; 20 corresponds to 5 years
    obs_names = ['ms.brine_aquifer1'+obs_ind,'ms.CO2_aquifer1'+obs_ind]
    xlabels = [r"Log$_{10}$ brine leakage rate (kg/s)", r"Log$_{10}$ CO$_2$ leakage rate (kg/s)"]

    f4, ax = plt.subplots(1, 2, figsize=(13.5, 6))

    data1 = sprior.recarray[obs_names[0]]
    data2 = sc.recarray[obs_names[0]]
    range_min = min(np.min(np.log10(data1)), np.min(np.log10(data2)))
    range_max = max(np.max(np.log10(data1)), np.max(np.log10(data2)))
    ax[0].hist(np.log10(data1), bins=60, label='Prior',
               range=(range_min, range_max), alpha=0.7, density=True)
    ax[0].hist(np.log10(data2), bins=60, label='Posterior',
               range=(range_min, range_max), alpha=0.7, density=True)

    data3 = sprior.recarray[obs_names[1]]
    data4 = sc.recarray[obs_names[1]]
    range_min = min(np.min(np.log10(data3)), np.min(np.log10(data4)))
    range_max = max(np.max(np.log10(data3)), np.max(np.log10(data4)))
    ax[1].hist(np.log10(data3), bins=60, label='Prior',
               range=(range_min, range_max), alpha=0.7, density=True)
    ax[1].hist(np.log10(data4), bins=60, label='Posterior',
               range=(range_min,range_max), alpha=0.7, density=True)

    ax[0].set_ylim(0.0, 3.0)
    ax[1].set_ylim(0.0, 3.0)
    ax[0].legend(fontsize=16)
    ax[1].legend(fontsize=16)
    ax[0].set_ylabel("Density", fontsize=18)
    ax[0].set_xlabel(xlabels[0], fontsize=18)
    ax[1].set_xlabel(xlabels[1], fontsize=18)
    xticks_brine = np.linspace(-7., -3., num=5)
    xticks_CO2 = np.linspace(-5., -2., num=4)
    yticks = np.linspace(0., 3., num=7)
    ax[0].set_xticks(xticks_brine)
    ax[1].set_xticks(xticks_CO2)
    ax[0].set_xticklabels(xticks_brine, fontsize=16)
    ax[1].set_xticklabels(xticks_CO2, fontsize=16)
    ax[0].set_yticks(yticks)
    ax[1].set_yticks(yticks)
    ax[0].set_yticklabels(yticks, fontsize=16)
    ax[1].set_yticklabels(yticks, fontsize=16)

    f4.show()
    f4.savefig('mcmc_leakage_panels.png')

    # plot CDF for prior and posterior distribution of observations
    f5, ax = plt.subplots(1, 2, figsize=(12, 6))
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    num_bins = 100
    hist, bin_edges = np.histogram(np.log10(data1), bins=num_bins)
    cdf = np.cumsum(hist)
    cdf1 = np.zeros(num_bins)
    for i in range(0, num_bins):
        cdf1[i] = float(cdf[i])/nMCMCSamples
    hist_post, bin_edges_post = np.histogram(np.log10(data2), bins=num_bins)
    cdf_post = np.cumsum(hist_post)
    cdf1_post = np.zeros(num_bins)
    for i in range(0, num_bins):
        cdf1_post[i] = float(cdf_post[i])/nMCMCSamples

    ax[0].plot(bin_edges[1:], cdf1, label="Prior", alpha=0.8, linewidth=4)

    ax[0].plot(bin_edges_post[1:], cdf1_post, label="Posterior", alpha=0.8,
               linewidth=4, color='orange')
    ax[0].plot([bin_edges_post[-1], bin_edges[-1]], [cdf1_post[-1], 1], alpha=0.8,
               linewidth=4, color='orange')
    ax[0].plot([bin_edges[0], bin_edges_post[0]], [0.0, 0.0], alpha=0.8,
               linewidth=4, color='orange')
    prob = 0.95  # 0.956
    ax[0].plot([bin_edges[0], bin_edges[-1]], [prob, prob], '--k', linewidth=2.5)
    ax[0].plot([bin_edges[1:][cdf1>=prob][0], bin_edges[1:][cdf1>=prob][0]],
               [0.0, prob], '--k', linewidth=2.5)
    ax[0].plot([bin_edges_post[1:][cdf1_post>=prob][0],
                bin_edges_post[1:][cdf1_post>=prob][0]],
               [0.0, prob], '--k', linewidth=2.5)
    ax[0].plot([bin_edges[0], bin_edges[-1]], [0.0, 0.0], '--k', linewidth=2.5)
    ax[0].set_xlabel(xlabels[0], fontsize=18)
    ax[0].set_ylabel('CDF', fontsize=18)
    ax[0].legend(loc='center left', fontsize=16)

    hist, bin_edges = np.histogram(np.log10(data3), bins=num_bins)
    cdf = np.cumsum(hist)
    cdf1 = np.zeros(num_bins)
    for i in range(0, num_bins):
        cdf1[i] = float(cdf[i])/nMCMCSamples

    hist_post, bin_edges_post = np.histogram(np.log10(data4), bins=num_bins)
    cdf_post = np.cumsum(hist_post)
    cdf1_post = np.zeros(num_bins)
    for i in range(0, num_bins):
        cdf1_post[i] = float(cdf_post[i])/nMCMCSamples
    ax[1].plot(bin_edges[1:], cdf1, label="Prior", alpha=0.8, linewidth=4)

    ax[1].plot(bin_edges_post[1:], cdf1_post, label="Posterior", alpha=0.8,
               linewidth=4, color='orange')
    ax[1].plot([bin_edges_post[-1], bin_edges[-1]], [cdf1_post[-1], 1], alpha=0.8,
               linewidth=4, color='orange')
    ax[1].plot([bin_edges[0], bin_edges_post[0]], [0.0, 0.0], alpha=0.8,
               linewidth=4, color='orange')
    prob = 0.95
    ax[1].plot([bin_edges[0], bin_edges[-1]], [prob, prob], '--k', linewidth=2.5)
    ax[1].plot([bin_edges[1:][cdf1>=prob][0], bin_edges[1:][cdf1>=prob][0]],
               [0.0, prob], '--k', linewidth=2.5)
    ax[1].plot([bin_edges_post[1:][cdf1_post>=prob][0],
                bin_edges_post[1:][cdf1_post>=prob][0]],
               [0.0, prob], '--k', linewidth=2.5)
    ax[1].plot([bin_edges[0], bin_edges[-1]], [0.0, 0.0], '--k', linewidth=2.5)
    ax[1].set_xlabel(xlabels[1], fontsize=18)
    ax[1].legend(loc='center left', fontsize=16)
    f5.savefig("CDF.png")
    f5.show()
    print('Uncertainty Quantification: Done!')
