'''
This example illustrates linking of simple reservoir, open wellbore and
atmosphere models. The saturation/pressure output produced by several simple
reservoir components is used to drive leakage from open wellbores.
|CO2| leakage rates are passed to the atmosphere model.
Plots of relevant observations are created.
Example illustrates setup for 20000 realizations Latin Hypercube Sampling study.

Example of run:
$ python iam_sys_reservoir_openwell_atmosphere_3locs_lhs.py
'''

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import SystemModel, SimpleReservoir, OpenWellbore, AtmosphericROM
import openiam.visualize as viz


if __name__ == '__main__':
    # For multiprocessing in Spyder
    __spec__ = None

    # Define keyword arguments of the system model
    num_periods = 4
    num_years = 1
    num_points = num_years*num_periods
    time_array = 7*np.arange(0.0, num_points+1) # in days
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Create randomly located leaky well locations
    # within box defined by xmin,xmax,ymin,ymax
    num_wells = 3
    well_xys = np.array([[600, 630], [490, 450], [270, 340]])

    sress = []
    mss = []
    ow = []
    co2_leakrates_collector = []
    for i, crds in enumerate(well_xys):

        # Add reservoir components
        sress.append(sm.add_component_model_object(
            SimpleReservoir(name='sres'+str(i), parent=sm,
                            injX=350., injY=400., locX=crds[0], locY=crds[1])))

        # Add parameters of reservoir component model
        sress[-1].add_par('numberOfShaleLayers', value=3, vary=False)
        # sress[-1].add_par('injRate', min=0.04, max=0.2, value=0.05)
        sress[-1].add_par('injRate', value=0.1, vary=False)
        if i == 0:
            sress[-1].add_par('reservoirThickness', value=20.0, vary=False)
            sress[-1].add_par('shale1Thickness', min=300.0, max=500.,
                              value=400.0, vary=False)
            sress[-1].add_par('aquifer1Thickness', min=20.0, max=60.,
                              value=50.0, vary=False)
            sress[-1].add_par('shale2Thickness', min=700.0, max=900.,
                              value=800.0, vary=False)
            sress[-1].add_par('aquifer2Thickness', min=60.0, max=80.,
                              value=75.0, vary=False)
            sress[-1].add_par('shale3Thickness', min=150.0, max=250.,
                              value=200.0, vary=False)
            sress[-1].add_par('logResPerm', min=-12.5, max=-11.0, value=-12.)
        else:
            # sress[-1].add_par_linked_to_par(
            #     'shale1Thickness', sress[0].pars['shale1Thickness'])
            # sress[-1].add_par_linked_to_par(
            #     'aquifer1Thickness', sress[0].pars['aquifer1Thickness'])
            # sress[-1].add_par_linked_to_par(
            #     'shale2Thickness', sress[0].pars['shale2Thickness'])
            # sress[-1].add_par_linked_to_par(
            #     'aquifer2Thickness', sress[0].pars['aquifer2Thickness'])
            # sress[-1].add_par_linked_to_par(
            #     'shale3Thickness', sress[0].pars['shale3Thickness'])
            sress[-1].add_par_linked_to_par(
                'reservoirThickness', sress[0].deterministic_pars['reservoirThickness'])
            sress[-1].add_par_linked_to_par(
                'shale1Thickness', sress[0].deterministic_pars['shale1Thickness'])
            sress[-1].add_par_linked_to_par(
                'aquifer1Thickness', sress[0].deterministic_pars['aquifer1Thickness'])
            sress[-1].add_par_linked_to_par(
                'shale2Thickness', sress[0].deterministic_pars['shale2Thickness'])
            sress[-1].add_par_linked_to_par(
                'aquifer2Thickness', sress[0].deterministic_pars['aquifer2Thickness'])
            sress[-1].add_par_linked_to_par(
                'shale3Thickness', sress[0].deterministic_pars['shale3Thickness'])

            sress[-1].add_par_linked_to_par(
                'logResPerm', sress[0].pars['logResPerm'])

        # Add observations of reservoir component model to be used by the next component
        sress[-1].add_obs_to_be_linked('pressure')
        sress[-1].add_obs_to_be_linked('CO2saturation')

        # Add observations of reservoir component model
        sress[-1].add_obs('pressure')
        sress[-1].add_obs('CO2saturation')

        # Add open wellbore component
        ow.append(sm.add_component_model_object(
            OpenWellbore(name='ow'+str(i), parent=sm)))

        # Add parameters of open wellbore component
        ow[-1].add_par('wellRadius', min=0.01, max=0.2, value=0.05, vary=False)
        ow[-1].add_par('wellTop', value=0.0, vary=False)
        if i == 0:
            ow[-1].add_composite_par('logReservoirTransmissivity',
                                     expr='13010299/10000000 + sres0.logResPerm')
            # ow[-1].add_par('logReservoirTransmissivity', min=-11.2, max=-8.8, value=-10.0)
            ow[-1].add_par('logAquiferTransmissivity', min=-11.2, max=-8.8,
                           value=-10.0, vary=False)
            ow[-1].add_par('brineSalinity', min=0.02, max=0.18,
                           value=0.1, vary=False)
        else:
            ow[-1].add_composite_par('logReservoirTransmissivity',
                                     expr='13010299/10000000 + sres0.logResPerm')
            # ow[-1].add_par_linked_to_par('logReservoirTransmissivity',
            #                              ow[0].pars['logReservoirTransmissivity'])
            ow[-1].add_par_linked_to_par('logAquiferTransmissivity',
                                         ow[0].deterministic_pars['logAquiferTransmissivity'])
            ow[-1].add_par_linked_to_par('brineSalinity',
                                         ow[0].deterministic_pars['brineSalinity'])

        # Add keyword arguments of the open wellbore component model
        ow[-1].add_kwarg_linked_to_obs('pressure', sress[-1].linkobs['pressure'])
        ow[-1].add_kwarg_linked_to_obs('CO2saturation', sress[-1].linkobs['CO2saturation'])

        # Add composite parameter of open wellbore component
        ow[-1].add_composite_par(
            'reservoirDepth', expr='+'.join(['sres0.shale1Thickness',
                                             'sres0.shale2Thickness',
                                             'sres0.shale3Thickness',
                                             'sres0.aquifer1Thickness',
                                             'sres0.aquifer2Thickness']))

        # Add observations of open wellbore component model
        ow[-1].add_obs_to_be_linked('CO2_atm')
        # ow[-1].add_obs('CO2_aquifer')  # zero since well top is in atm
        # ow[-1].add_obs('brine_aquifer')  # zero since well top is in atm
        ow[-1].add_obs('CO2_atm')
        ow[-1].add_obs('brine_atm')

        # Create collector for leakrates
        co2_leakrates_collector.append(ow[-1].linkobs['CO2_atm'])

    # Add Atm ROM component
    satm = sm.add_component_model_object(AtmosphericROM(name='satm', parent=sm))
    satm.add_par('T_amb', value=15.0, vary=False)
    satm.add_par('P_amb', min=0.7, max=1.02, value=1.01)
    satm.add_par('V_wind', min=2.0, max=18.0, value=5.0)
    satm.add_par('C0_critical', value=0.01, vary=False)

    satm.model_kwargs['x_receptor'] = [250, 300, 350, 400, 580]
    satm.model_kwargs['y_receptor'] = [580, 660, 480, 525, 580]

    satm.model_kwargs['x_coor'] = well_xys[:, 0]
    satm.model_kwargs['y_coor'] = well_xys[:, 1]

    satm.add_kwarg_linked_to_collection('co2_leakrate', co2_leakrates_collector)

    # Add observations for receptors
    for i, r in enumerate(satm.model_kwargs['x_receptor']):
        satm.add_obs('outflag_r{0:03}'.format(i))

    n_sources = len(satm.model_kwargs['x_coor'])
    for i in range(n_sources):
        satm.add_obs('x_new_s{0:03}'.format(i))
        satm.add_obs('y_new_s{0:03}'.format(i))
        satm.add_obs('critical_distance_s{0:03}'.format(i))

    satm.add_obs('num_sources')

    # There is no forward method call as the corresponding results are not of interest
    print('Latin Hypercube sampling study started...')
    import random
    num_samples = 20000
    ncpus = 4
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=random.randint(123, 1456))
    # s = sm.parstudy(nvals=5)
    # num_samples = 5**4

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)
    print(results)
    print('------------------------------------------------------------------')
    print('                  End ')
    print('------------------------------------------------------------------')

    if not os.path.exists('atm_plumes'):
        os.mkdir('atm_plumes')
    # Plot single realization of plume plot
    to_plot_plumes = 0
    if to_plot_plumes:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        t1 = 1
        t2 = 2
        for time_index in range(t1, t2):
            for realzn in range(120, 170):
                fig, ax = plt.subplots(1, 2, sharey=True, figsize=(13.5, 6))
                fig.subplots_adjust(left=0.1, right=0.95, wspace=0.1)

                extent = ([200, 700], [300, 800])
                save_fig_name = os.path.join(
                    'atm_plumes', 'atm_forward_plume_sample_{}_t'.format(realzn)+str(time_index))

                plumes, ex = viz.get_plumes(s, time_index=time_index,
                                            sample=realzn, atm_name=satm.name)
                # print(plumes[-1].x, plumes[-1].y, ex)

                receptors = list(zip(satm.model_kwargs['x_receptor'],
                                     satm.model_kwargs['y_receptor']))

                for r in receptors:
                    ax[0].plot(r[0], r[1], 'k.', label='Marker', markersize=15)

                for p in plumes:
                    ax[0].plot(p.x, p.y, 'r.', label='Source', markersize=15)
                    e1 = mpatches.Ellipse((p.x, p.y), p.dx*2, p.dy*2, alpha=0.3)
                    ax[0].add_artist(e1)

                    ax[0].annotate("{0:.1f}m".format(p.dx),
                                    xy=(p.x, p.y),
                                    xytext=(p.x, p.y+1), size=16)

                ax[0].set_aspect(1)

                ax[0].set(xlim=extent[0], ylim=extent[1])
                ax[0].set_title('a.', fontsize=20)

                ax[0].set_xlabel('x, [m]', fontsize=18)
                ax[0].set_ylabel('y, [m]', fontsize=18)
                ax[0].set_xticklabels([int(val) for val in ax[0].get_xticks()], fontsize=16)
                ax[0].set_yticklabels([int(val) for val in ax[0].get_yticks()], fontsize=16)

                l, h = ax[0].get_legend_handles_labels()
                ptch = mpatches.Patch(alpha=0.3, label="Critical\n region")
                ax[0].legend(handles=[l[0], l[-1], ptch], loc='lower right', fontsize=16)

                # Second subplot: calculate once and comment out to process first subplots
                second_plot = 1
                if second_plot:
                    time_index = 1
                    plumes_grid, xy = viz.prob_plume_grid(s, satm, time_index+1,
                                                          extent=extent)
                    im = ax[1].tricontourf(xy[:, 0], xy[:, 1], plumes_grid,
                                           levels=np.linspace(0, 1.001, num=10),
                                           cmap='jet')

                    divider = make_axes_locatable(ax[1])
                    cax = divider.append_axes('right', size='5%', pad=0.05)

                    cb = fig.colorbar(im, cax=cax, orientation='vertical', format='%.1f')
                    ax[1].set_aspect(1)
                    ax[1].set(xlim=extent[0], ylim=extent[1])
                    cb.ax.tick_params(labelsize=16)
                    ax[1].set_xlabel('x, [m]', fontsize=18)
                    ax[1].set_xticklabels([int(val) for val in ax[1].get_xticks()],
                                          fontsize=16)

                    ax[1].set_title('b.', fontsize=20)
                    for r in receptors:
                        ax[1].plot(r[0], r[1], 'k.', label='Marker', markersize=15)

                fig.savefig(save_fig_name+'.png', dpi=300)
                plt.close(fig)

    to_print = 0
    if to_print:

        outputs = s.collect_observations_as_time_series()
        pressure = {}
        CO2saturation = {}

        for ind in range(3):
            pressure[ind] = outputs['sres{}.pressure'.format(ind)]
            CO2saturation[ind] = outputs['sres{}.CO2saturation'.format(ind)]

        line_width = 1.5
        alpha = 0.8
        fig, ax = plt.subplots(1, 2, figsize=(13.5,6))
        fig.subplots_adjust(left=0.1, right=0.95, wspace=0.1)
        colors = ['black', 'red', 'blue']
        for ind in range(3):
            for sim_ind in range(num_samples):
                ax[0].plot(time_array/365.25, pressure[ind][sim_ind, :],
                           color=colors[ind], linewidth=line_width, alpha=alpha)

                ax[1].plot(time_array/365.25, CO2saturation[ind][sim_ind, :],
                           color=colors[ind], linewidth=line_width, alpha=alpha)

        fig, ax = plt.subplots(1, 3, sharey=True, figsize=(13.5,6))
        fig.subplots_adjust(left=0.1, right=0.95, wspace=0.1)
        for ind in range(3):
            for sim_ind in range(num_samples):
                ax[ind].plot(time_array/365.25, pressure[ind][sim_ind, :],
                             color=colors[ind], linewidth=line_width, alpha=alpha)

        # Add capture point so sensitivity is calculated at time other than ending
        t_ind = 1
        title = ['a.', 'b.', 'c.']
        for ind in range(3):
            sens_obs = 'satm.critical_distance_s00{}_{}'.format(ind, t_ind)
            # Run sensitivity analysis at capture point
            atm_sensitivity = s.rbd_fast(obsname=sens_obs)

            viz.simple_sensitivities_barplot(
                atm_sensitivity, sm,
                title=title[ind],
                # ylabel='First Order Sensitivity\n(Calculated using RBD-Fast Method)',
                ylabel='First order sensitivity',
                savefig=os.path.join('atm_plumes', 'sens_crit_dist{}.png'.format(ind)),
                outfile='cr_d_sensitivities{}.txt'.format(ind))

        # Another set of plots
        for ind in range(3):
            viz.time_series_sensitivities(
                'satm.critical_distance_s00{}'.format(ind), sm, s, time_array,
                title='a.',#title[ind],
                savefig=os.path.join('atm_plumes',
                                     'critical_distance_s00{}.png'.format(ind)),
                outfile='crt_dist{}.csv'.format(ind))


        # Plot for many observations
        time_point = 1
        viz.multi_sensitivities_barplot([
            'ow0.CO2_atm_{}'.format(time_point),
            'ow1.CO2_atm_{}'.format(time_point),
            'ow2.CO2_atm_{}'.format(time_point),
            'satm.critical_distance_s000_{}'.format(time_point),
            'satm.critical_distance_s001_{}'.format(time_point),
            'satm.critical_distance_s002_{}'.format(time_point)],
            sm, s,
            savefig=os.path.join('atm_plumes',
                                 'multi-bar-sensitivities_b_{}.png'.format(time_point)),
            title='b.',
            outfile='multi_sens_b.csv')

        # Correlation coefficients
        excludes_obs = ['satm.x_new_s00{}'.format(ind) for ind in range(3)]+[
            'satm.y_new_s00{}'.format(ind) for ind in range(3)]+[
                'satm.outflag_r00{}'.format(ind) for ind in range(5)]+[
                    'satm.num_sources']+[
                        'ow{}.brine_atm'.format(ind) for ind in range(3)]+[
                            'sres{}.pressure'.format(ind) for ind in range(3)]+[
                                'sres{}.CO2saturation'.format(ind) for ind in range(3)]

        corr_coeffs = viz.correlations_at_time(
            s, capture_point=2,
            excludes=excludes_obs,
            ctype='pearson',
            plot=True, figsize=(6.75, 6),
            printout=False, xrotation=90,
            title='',
            savefig=os.path.join('atm_plumes',
                                  'Pearson_Corr_coeff_at_time_2.png'),
            outfile='Pearson_Corr_coeffs.csv')

        corr_coeffs = viz.correlations_at_time(
            s, capture_point=2,
            excludes=excludes_obs,
            ctype='spearman',
            plot=True, figsize=(6.75, 6),
            printout=False, xrotation=90,
            title='b.',
            savefig=os.path.join('atm_plumes', 'spearman_Corr_coeff_at_time_2.png'),
            outfile='spearman_Corr_coeffs.csv')


        s.panels(figsize=(15, 15))
        plt.savefig(os.path.join('atm_plumes', 'panels.png'), dpi=300)
