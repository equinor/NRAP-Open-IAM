'''
Example illustrates system model containing only simple reservoir model.

Examples of run:
$ python iam_concordance.py
'''
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, SimpleReservoir

def robustness(Pc, Ps, kus):
    '''
    Calculate robustness against future measurements exceeding a critical pressure

    :param Pc: Critical pressure threshold
    :type Pc: float
    :param Ps: List of arrays containing the maximum pressure prediction associated
    :type Ps: list(array(float)) or array(array(float))
    :param kus: List of arrays of deviation from nominal in units of standard errors
    :type kus: list(array(float))
    '''
    hhats = []
    for Pn1, ku in zip(Ps, kus):
        # If the minimum max pressure for the ensemble is greater than Pc,
        # set the robustness to zero
        if Pc < np.min(Pn1):
            hhats.append(0)
        else:
            # Create function to interpolate the horizon of uncertainty
            # as a function of max pressure prediction
            fi = interp1d(Pn1, ku)
            hhat = fi(Pc)
            # If horizon of uncertainty is less than zero, set to zero (zero robustness).
            if hhat < 0:
                hhats.append(0)
            else: # else, set robustness to hhat
                hhats.append(hhat)
    return hhats

if __name__ == '__main__':
    __spec__ = None
    #######################################################
    # EVALUATE ROBUSTNESS FOR ALTERNATIVE INJECTION RATES #
    #######################################################

    # Assume that you start checking conformance after the second measurement
    # is obtained (time_array[2]).
    # Perform the following analysis looping through remaining times:
    #  1. Fit the model to the available measurements, producing an estimate
    #     of reservoir permeability and standard error of the estimate.
    #  2. Using the estimated permeability, predict future pressures.
    #  3. Predict pressures for an ensemble of permeabilities.
    #  4. Calculate the number of standard errors of each k in the sample
    #     from the current nominal value (pf['sres.logResPerm']).
    #  5. Collect the max pressure for predictions for each k in the sample
    arr_stderrs = [] # least square standard error of perm estimate
    # deviation of res perm from least square estimate divided by least square standard error
    arr_kus = []
    arr_kfits = []
    arr_Pnps = [] # max pressure predictions for permeability samples
    arr_Ps = [] # Simulated pressures for all time step
    #injRates = [0.075,0.1,0.125] # Injection rates, m^3/s
    injRates_Mty = np.array([0.75, 1.0, 1.25]) # Injection rates, Mt/y
    # Assuming 20C at the surface and a 0.03C/m geothermal gradient and
    # 100000 Pa at the surface and 1000. kg/m3 water density, the density of CO2
    # in the reservoir would be around 328 kg/m3. We can then calculate
    # the injection rate in kg/s as
    injRates = injRates_Mty * 1e9 / (365 * 86400 * 328) # Injection rates, m^3/s
    noisy_pressures = []
    for i, injRate in enumerate(injRates):
        #####################
        # CREATE NOISY DATA #
        #####################

        # Create system model up to present (180 days)
        time_array = np.linspace(0, 480, 480//15+1)   # in days
        sm_model_kwargs = {'time_array': time_array}  # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)
        # Add reservoir component
        sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))
        # Add parameters of reservoir component model
        sres.add_par('logResPerm', min=-14, max=-9, value=-11)
        # Set injection rate parameter, but do not vary
        sres.add_par('injRate', min=0.05, max=0.15, value=injRate, vary=False)
        # Add observations of reservoir component model
        sres.add_obs('pressure')
        # Run system model using current values of its parameters
        sm.forward()
        # Assign observations of the model to pressure and CO2saturation variables
        pressure = sm.collect_observations_as_time_series(sres, 'pressure')
        # Add noise to data and add as observations to system model
        np.random.seed(3000)
        noisy_pressure = pressure + np.random.normal(loc=0, scale=5e4, size=len(pressure))
        noisy_pressures.append(noisy_pressure)
        for o, n in zip(sm.obs.values(), noisy_pressure):
            o.value = n

        stderrs = [] # least square standard error of perm estimate
        # deviation of res perm from least square estimate divided by least square standard error
        kus = []
        kfits = []
        Pnps = [] # max pressure predictions for permeability samples
        Ps = [] # Simulated pressures for all time step
        for t in time_array[2:-1]:
            # Create system model with only three observations
            # time is given in days
            sm1 = SystemModel(model_kwargs={'time_array': np.arange(0, t+1, 15)})
            # Add reservoir component
            sres1 = sm1.add_component_model_object(SimpleReservoir(name='sres', parent=sm1))
            # Add parameters of reservoir component model
            sres1.add_par('logResPerm', min=-14, max=-9, value=-10.5)
            # Set injection rate parameter, but do not vary
            sres1.add_par('injRate', min=0.05, max=0.15, value=injRate, vary=False)
            # Add observations of reservoir component model
            # Note: Since sm1.obs.values() is shorter then noisy_pressure,
            # the 'zip' will truncate to sm1.obs.values()
            sres1.add_obs('pressure')
            for o, n in zip(sm1.obs.values(), noisy_pressure):
                o.value = n

            # Perform least squares fit and collect stderr
            lf, pf = sm1.lmfit(verbose=False, difference_type='central')
            kfits.append(pf['sres.logResPerm'].value)
            stderrs.append(pf['sres.logResPerm'].stderr)
            # Run system model using current values of its parameters
            # to ensure that optimal pressures are in observations
            sm1.forward()

            # Create system model with all times to predict future pressures
            # to perform robustness analysis
            # Add remaining noisy pressures to observations
            sm2 = SystemModel(model_kwargs={'time_array': time_array})
            sres2 = sm2.add_component_model_object(SimpleReservoir(name='sres', parent=sm2))
            # Add parameters of reservoir component model
            sres2.add_par('logResPerm', min=-14, max=-9, value=kfits[-1])
            # Set injection rate parameter, but do not vary
            sres2.add_par('injRate', min=0.05, max=0.15, value=injRate, vary=False)
            # Add observations of reservoir component model
            # Note: Since sm1.obs.values() is shorter then noisy_pressure,
            # the 'zip' will truncate to sm1.obs.values()
            sres2.add_obs('pressure')
            for o, n in zip(sm2.obs.values(), noisy_pressure):
                o.value = n
            # Run model forward to collect predicted pressures
            Ps.append(np.array(list(sm2.forward().values())))
            # Run sampleset of reservoir perms
            ks = np.linspace(-10., -13., 21)
            ss = sm2.create_sampleset(ks.reshape(21, 1))
            ss.run(verbose=False)
            # Calculate robustnesses
            kus.append((pf['sres.logResPerm'].value - ks)/stderrs[-1])
            # Collect pressure predictions
            Pnames = [v for v in ss.responses.recarray.dtype.names if v not in sm1.obsnames]
            # Collect max pressure predictions for each ks and convert to MPa
            Pnps.append(np.max(ss.responses.recarray[Pnames].tolist(), axis=1)/1.0e+6)

            del sm1, sm2, ss, lf, pf
        arr_stderrs.append(stderrs)
        arr_kus.append(kus)
        arr_kfits.append(kfits)
        arr_Pnps.append(Pnps)
        arr_Ps.append(Ps)

    # Plot robustness over time for different injection rates
    # Assuming 20C at the surface and a 0.03C/m geothermal gradient and 100000 Pa
    # at the surface and 1000. kg/m3 water density, the density of CO2
    # in the reservoir would be around 328 kg/m3. We can then approximate
    # the injection rate in Mt/y as
    # injRates_Mty = np.round(365 * 24*3600 * np.array(injRates) * 328 / 1.e9, decimals=2)
    f,ax = plt.subplots(1)
    for aPnps, akus, injRate, injRate_Mty in zip(arr_Pnps, arr_kus, injRates, injRates_Mty):
        ax.plot(time_array[2:-1], robustness(10.2, aPnps, akus),
                label=r'$Q$='+str(injRate_Mty)+r' Mt/y')
    ax.fill_between(time_array[2:-1], 0, 1, color='gray', alpha=0.35)
    ax.fill_between(time_array[2:-1], 1, 2, color='gray', alpha=0.25)
    ax.fill_between(time_array[2:-1], 2, 3, color='gray', alpha=0.15)
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Conformance Robustness')
    ax.axhline(1, ls='--', lw=0.5, c='grey')
    ax.text(460, 0.4,r'P($\mu-\sigma \leq k \leq \mu$) = 34.1 %', ha='right')
    ax.axhline(2, ls='--', lw=0.5, c='grey')
    ax.text(460, 1.4,r'P($\mu-2\sigma\leq k \leq \mu-\sigma$) = 13.6 %', ha='right')
    ax.axhline(3, ls='--', lw=0.5, c='grey')
    ax.text(460, 2.4, r'P($\mu-3\sigma\leq k \leq \mu-2\sigma$) = 2.1 %', ha='right')
    ax.text(460, 3.4, r'P($k \leq \mu-3\sigma$ ) = 0.1 %', ha='right')
    ax.set_xlim(30,time_array[-2])
    ax.set_ylim(bottom=0)
    f.show()
    f.savefig('rob_Q_ts.pdf')

    ###################################
    # Make plots for Q=0.1 m^2/s case #
    ###################################

    # Plot measured pressures, pressure prediction curves, and critical pressure thresholds
    f,ax = plt.subplots(1)
    ax.scatter(time_array, noisy_pressures[1]/1.0e+6, label='noisy pressure',
               edgecolors='k', facecolors="None")
    ax.axhline(10.0, ls='--', c='b', lw=1)
    ax.axhline(10.1, ls='--', c='r', lw=1)
    ax.axhline(10.2, ls='--', c='g', lw=1)
    ax.set_xlabel('Time, [days]')
    ax.set_ylabel('Pressure, (MPa)')
    ax.text(0., 10.005, r'$P_c$ = 10.0 MPa', color='b')
    ax.text(0., 10.105, r'$P_c$ = 10.1 MPa', color='r')
    ax.text(0., 10.205, r'$P_c$ = 10.2 MPa', color='g')
    # Only plot predicted pressures, skipping pressure up to the time the fit is performed
    for i,P in enumerate(arr_Ps[1]):
        ax.plot(time_array[2+i:], P[2+i:]/1.0e+6, lw=0.7, c='k')
    #for i,P in enumerate(arr_Ps[1]): ax.plot(time_array[2+i:],P[2+i:]/1.0e+6,lw=0.7,alpha=0.80)
    f.show()
    f.savefig('data_Pc.pdf')

    # Plot the standard error of the estimate over time as more measurements become available
    f,ax = plt.subplots(1)
    ax.scatter(time_array[2:-1], arr_stderrs[1], facecolors="None", edgecolors='k')
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Log Permeability Standard Error')
    f.show()
    f.savefig('stderr.pdf')

    # Plot permeability estimates with standard error bars over time
    # as more measurements are obtained
    f,ax = plt.subplots(1)
    ax.errorbar(time_array[2:-1], arr_kfits[1], yerr=arr_stderrs[1], fmt='ko')
    ax.set_xlabel('Time [days]')
    ax.set_ylabel(r'Log Reservoir Permeability [log$_{10}]$')
    ax.axhline(-11, ls='--', c='k', lw=1)
    ax.text(100, -10.80, 'Circles are least-squares permeability estimates.',
            style='italic')
    ax.text(100,-10.82, 'Vertical lines represent the least-squares standard error.',
            style='italic')
    ax.text(100,-10.84, "The dashed line is the 'true' log reservoir permeability.",
            style='italic')
    f.tight_layout()
    f.show()
    f.savefig('uncert_k.pdf')

    # Plot robustness over time for different critical pressures
    f,ax = plt.subplots(1)
    for Pc in [10.0,10.1,10.2]:
        ax.plot(time_array[2:-1], robustness(Pc, arr_Pnps[1], arr_kus[1]),
                label=r'$P_c$='+str(Pc)+' MPa')
    ax.fill_between(time_array[2:-1], 0, 1, color='gray', alpha=0.35)
    ax.fill_between(time_array[2:-1], 1, 2, color='gray', alpha=0.25)
    ax.fill_between(time_array[2:-1], 2, 3, color='gray', alpha=0.15)
    ax.legend()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Conformance Robustness')
    ax.axhline(1, ls='--', lw=0.5, c='grey')
    ax.text(460, 0.4 ,r'P($\mu-\sigma \leq k \leq \mu$) = 34.1 %', ha='right')
    ax.axhline(2, ls='--', lw=0.5, c='grey')
    ax.text(460, 1.4, r'P($\mu-2\sigma\leq k \leq \mu-\sigma$) = 13.6 %', ha='right')
    ax.axhline(3, ls='--', lw=0.5, c='grey')
    ax.text(460, 2.4, r'P($\mu-3\sigma\leq k \leq \mu-2\sigma$) = 2.1 %', ha='right')
    ax.text(460, 3.4, r'P($k \leq \mu-3\sigma$ ) = 0.1 %', ha='right')
    ax.set_xlim(30, time_array[-2])
    ax.set_ylim(bottom=0)
    f.show()
    f.savefig('rob_Pc_ts.pdf')

    # Plot pressure prediction bands as more monitoring data is collected
    sm = SystemModel(model_kwargs={'time_array': time_array})
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))
    # Add parameters of reservoir component model
    sres.add_par('logResPerm', min=-14, max=-9, value=kfits[-1])
    # Set injection rate parameter, but do not vary
    sres.add_par('injRate', min=0.05, max=0.15, value=injRates[1], vary=False)
    # Add observations of reservoir component model
    sres.add_obs('pressure')
    # Collect pressures from model with best fit k and best fit k plus
    # and minus the standard error
    f, ax = plt.subplots(3, sharey=True, sharex=True)
    for axi, i in enumerate([0, 8, 16]):
        pu = sm.forward(pardict={'sres.logResPerm':arr_kfits[1][i]+arr_stderrs[1][i]})
        p = sm.forward(pardict={'sres.logResPerm':arr_kfits[1][i]})
        pl = sm.forward(pardict={'sres.logResPerm':arr_kfits[1][i]-arr_stderrs[1][i]})
        ax[axi].scatter(time_array[:2+i], noisy_pressures[1][:2+i]/1.0e+6, 5,
                        label='noisy pressure', edgecolors='k', facecolors="None")
        ax[axi].plot(time_array[1+i:], np.array(list(p.values())[1+i:])/1.0e+6,
                     color='blue')
        ax[axi].fill_between(time_array[1+i:], np.array(list(pu.values())[1+i:])/1.0e+6,
                             np.array(list(pl.values())[1+i:])/1.0e+6,
                             color='blue', alpha=0.20)
        ax[axi].plot(time_array[:2+i], np.array(list(p.values())[:2+i])/1.0e6,
                     color='gray')
        ax[axi].fill_between(time_array[:2+i], np.array(list(pu.values())[:2+i])/1.0e+6,
                             np.array(list(pl.values())[:2+i])/1.0e+6, color='gray', alpha=0.20)
        ax[axi].text(50, 9.45, r'$\sigma$ = '+str(np.round(arr_stderrs[1][i], decimals=3)))
    ax[0].text(475, 9.45, r'2 measurements available', ha='right')
    ax[0].text(475, 10.05, r'P($\mu-\sigma \leq k \leq \mu+\sigma$) = 68.2 %', ha='right')
    ax[1].text(475, 9.45, r'10 measurements available', ha='right')
    ax[0].text(150, 10.4, r'Gray - best fit; Blue - prediction')
    ax[2].text(475, 9.45, r'18 measurements available', ha='right')
    ax[2].set_xlabel('Time, [days]')
    ax[1].set_ylabel('Pressure, (MPa)')
    f.show()
    f.savefig('P_unc.pdf')
