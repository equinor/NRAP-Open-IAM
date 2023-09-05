import os
import sys
import unittest
import warnings
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import h5py
from shutil import copyfile
import matplotlib.pyplot as plt

try:
    from openiam import (SystemModel, Stratigraphy,
                         SimpleReservoir, AnalyticalReservoir,
                         GenericReservoir, TheisReservoir,
                         MultisegmentedWellbore, CementedWellbore,
                         CementedWellboreWR, OpenWellbore, KimberlinaWellbore,
                         GeneralizedFlowRate, HydrocarbonLeakage, RateToMassAdapter,
                         CarbonateAquifer, GenericAquifer, AtmosphericROM,
                         ReservoirDataInterpolator, LookupTableReservoir,
                         AlluviumAquifer, AlluviumAquiferLF,
                         DeepAlluviumAquifer, DeepAlluviumAquiferML,
                         FutureGen2Aquifer, FutureGen2AZMI,
                         PlumeStability, SealHorizon, FaultFlow, FaultLeakage,
                         ChemicalWellSealing, SHPermeabilitySampler,
                         SHThicknessSampler)
    import openiam.visualize as iam_vis
except ModuleNotFoundError:
    try:
        sys.path.append(os.sep.join(['..', 'source']))
        from openiam import (SystemModel, Stratigraphy,
                             SimpleReservoir, AnalyticalReservoir,
                             GenericReservoir, TheisReservoir,
                             MultisegmentedWellbore, CementedWellbore,
                             CementedWellboreWR, OpenWellbore, KimberlinaWellbore,
                             GeneralizedFlowRate, HydrocarbonLeakage, RateToMassAdapter,
                             CarbonateAquifer, GenericAquifer, AtmosphericROM,
                             ReservoirDataInterpolator, LookupTableReservoir,
                             AlluviumAquifer, AlluviumAquiferLF,
                             DeepAlluviumAquifer, DeepAlluviumAquiferML,
                             FutureGen2Aquifer, FutureGen2AZMI,
                             PlumeStability, SealHorizon, FaultFlow, FaultLeakage,
                             ChemicalWellSealing, SHPermeabilitySampler,
                             SHThicknessSampler)
        import openiam.visualize as iam_vis
    except ImportError as err:
        print('Unable to load IAM class module: {}'.format(err))

CURRENT_WORK_DIR = os.getcwd()


class Tests(unittest.TestCase):

    def setUp(self):
        """Defines the actions performed before each test."""
        # return to original directory
        os.chdir(CURRENT_WORK_DIR)

    def shortDescription(self):
        """Defines information that is printed about each of the tests.

        Method redefines method of TestCase class. In the original method
        only the first line of docstring is printed. In this method
        the full docstring of each test in Tests class is printed.
        """
        doc = self._testMethodDoc
        return doc

    def test_alluvium_aquifer(self):
        """Tests basic alluvium aquifer.

        Tests an alluvium aquifer component in a
        forward model against expected output for 2 years of data.
        """

        time_array = 365.25*np.arange(0.0, 2.0)
        sm_model_kwargs = {'time_array': time_array} # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add alluvium aquifer model object and define parameters
        aa = sm.add_component_model_object(AlluviumAquifer(name='aa', parent=sm))

        aa.add_par('sandFraction', value=0.422)
        aa.add_par('correlationLengthX', value=361.350)
        aa.add_par('correlationLengthZ', value=17.382)
        aa.add_par('permeabilityClay', value=-16.340)
        aa.add_par('NaMolality', value=0.121)
        aa.add_par('PbMolality', value=-5.910)
        aa.add_par('benzeneMolality', value=0.109)
        aa.add_par('tMitigation', value=87.914)
        aa.add_par('CEC', value=32.073)
        aa.add_par('Asbrine', value=-5.397)
        aa.add_par('Babrine', value=-3.397)
        aa.add_par('Cdbrine', value=-8.574)
        aa.add_par('Pbbrine', value=-7.719)
        aa.add_par('Benzene_brine', value=-8.610)
        aa.add_par('Benzene_kd', value=-3.571)
        aa.add_par('Benzene_decay', value=-2.732)
        aa.add_par('PAH_brine', value=-7.118)
        aa.add_par('PAH_kd', value=-0.985)
        aa.add_par('PAH_decay', value=-3.371)
        aa.add_par('phenol_brine', value=-6.666)
        aa.add_par('phenol_kd', value=-1.342)
        aa.add_par('phenol_decay', value=-3.546)
        aa.add_par('porositySand', value=0.468)
        aa.add_par('densitySand', value=2165.953)
        aa.add_par('VG_mSand', value=0.627)
        aa.add_par('VG_alphaSand', value=-4.341)
        aa.add_par('permeabilitySand', value=-12.430)
        aa.add_par('Clbrine', value=-0.339)
        aa.add_par('calcitevol', value=0.165)
        aa.add_par('V_goethite', value=0.004)
        aa.add_par('V_illite', value=0.006)
        aa.add_par('V_kaolinite', value=0.004)
        aa.add_par('V_smectite', value=0.010)

        aa.model_kwargs['co2_rate'] = 1.90e-02 # kg/s
        aa.model_kwargs['co2_mass'] = 6.00e+05 # kg
        aa.model_kwargs['brine_rate'] = 4.62e-03 # kg/s
        aa.model_kwargs['brine_mass'] = 1.46e+05 # kg

        # Add observations (output) from the alluvium aquifer model
        aa.add_obs('TDS_volume')
        aa.add_obs('pH_volume')

        # Run the system model
        sm.forward()

        # Assign observations of the model to variables
        tds = sm.collect_observations_as_time_series(aa, 'TDS_volume')
        ph = sm.collect_observations_as_time_series(aa, 'pH_volume')

        # True values
        true_tds = [0.0, 6878026.69860782]
        true_ph = [0.0, 5742563782.01]

        # Test with helpful message!
        for ts, s, t in zip(true_tds, tds, time_array):
            self.assertTrue(abs((ts-s)) < 1.0,
                            'TDS volume at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))
        for ts, s, t in zip(true_ph, ph, time_array):
            self.assertTrue(abs((ts-s)) < 1.0,
                            'pH volume at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))

    def test_alluvium_aquifer_lf(self):
        """Tests an alluvim aquifer for low input fluxes.

        Tests an alluvium aquifer component in a forward simulation
        against expected output available for 6 years of data.
        """
        # Create system model
        time_array = 365.25*np.arange(0.0, 6.0)
        sm_model_kwargs = {'time_array': time_array} # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add alluvium aquifer model object and define parameters
        aalf = sm.add_component_model_object(AlluviumAquiferLF(name='aalf', parent=sm))
        aalf.add_par('sandFraction', value=0.7087)
        aalf.add_par('correlationLengthX', value=473.9740)
        aalf.add_par('correlationLengthZ', value=22.5966)
        aalf.add_par('logK_sand', value=-13.3634)
        aalf.add_par('logK_clay', value=-15.5075)
        aalf.add_par('NaMolality', value=0.9279)
        aalf.add_par('PbMolality', value=-6.2021)
        aalf.add_par('benzeneMolality', value=-7.7138)
        aalf.add_par('tMitigation', value=7.7387)

        aalf.model_kwargs['co2_rate'] = 0.000006883      # kg/s
        aalf.model_kwargs['co2_mass'] = 10**2.73481      # kg
        aalf.model_kwargs['brine_rate'] = 0.000145878    # kg/s
        aalf.model_kwargs['brine_mass'] = 10**4.36206    # kg

        # Add observations (output) of the alluvium aquifer model
        aalf.add_obs('TDS_volume')
        aalf.add_obs('pH_volume')
        aalf.add_obs('As_volume')
        aalf.add_obs('Pb_volume')
        aalf.add_obs('Cd_volume')
        aalf.add_obs('Ba_volume')
        aalf.add_obs('Benzene_volume')
        aalf.add_obs('Naphthalene_volume')
        aalf.add_obs('Phenol_volume')

        # Run the system model
        sm.forward()

        # True values
        true_vals = {}
        true_vals['TDS_volume'] = [0., 243.60472506, 315.42085149, 681.51000667,
                                   1464.27143244, 2817.58361957]
        true_vals['pH_volume'] = [0., 135.55234468, 547.96051906, 605.7159738,
                                  606.48105472, 678.17416717]
        true_vals['As_volume'] = [0., 10., 17.32817682, 25.83795442,
                                  53.71518848, 97.3487519]
        true_vals['Pb_volume'] = [0., 179.47389199, 1554.30216956, 5972.94243509,
                                  13346.66092414, 24999.17432054]
        true_vals['Cd_volume'] = [0., 105.02705854, 2362.73319012, 9383.10921071,
                                  11061.49479283, 8258.10081888]
        true_vals['Ba_volume'] = [0., 3973.21653329, 8082.43147709, 10288.42580438,
                                  12179.40618369, 13643.12399265]
        true_vals['Benzene_volume'] = [0., 62.15296578, 353.59654512, 714.91883593,
                                       1223.97332442, 2131.52468324]
        true_vals['Naphthalene_volume'] = [0., 38.76136875, 133.63125123, 538.19287788,
                                           1639.28247224, 3103.13930929]
        true_vals['Phenol_volume'] = [0., 4436.94112399, 4309.9468996, 4803.11442366,
                                      7642.13761712, 12014.17764713]

        indices = [5]
        for key, vals in true_vals.items():
            for ind in indices:
                true_val = vals[ind]
                sim_val = aalf.obs[key+'_{}'.format(ind)].sim
                self.assertTrue(
                    abs(true_val-sim_val)/true_val < 0.01,
                    'Affected {} volume at {} days is {} but should be {}.'.format(
                        key, time_array[ind], sim_val, true_val))

    def test_analytical_reservoir_forward(self):
        """Tests analytical reservoir component.

        Tests the system model with a analytical reservoir component in a
        forward model against expected output for 4 years of data.
        """
        # Create system model
        num_years = 4
        time_array = 365.25*np.arange(0, num_years+1)

        model_kwargs = {'time_array': time_array}  # time is given in days
        sm = SystemModel(model_kwargs=model_kwargs)

        # Add analytical reservoir component
        res = sm.add_component_model_object(AnalyticalReservoir(
            name='res', parent=sm, injX=100., injY=0., locX=0., locY=0.))

        res.add_par('reservoirThickness', value=30.0)
        res.add_par('shaleThickness', value=970.0)
        res.add_par('aquiferThickness', value=30.0)
        res.add_par('logResPerm', value=-13.69897)
        res.add_par('reservoirPorosity', value=0.15)
        res.add_par('datumPressure', value=101325.0)
        res.add_par('brineDensity', value=1045.0)
        res.add_par('CO2Density', value=479.0)
        res.add_par('brineViscosity', value=2.535e-4)
        res.add_par('CO2Viscosity', value=3.95e-5)
        res.add_par('brineResSaturation', value=0.0)
        res.add_par('injRate', value=8.87/479.0)
        res.add_par('reservoirRadius', value=500)
        res.add_par('brineCompressibility', value=1.e-11)

        # Add observations of analytical reservoir component
        res.add_obs('pressure')
        res.add_obs('CO2saturation')
        res.add_obs('pressureAve')

        # Run system model using current values of its parameters
        sm.forward()

        # Assign observations of the model to pressure and CO2saturation variables
        # Obtain pressure and CO2 saturation as lists
        pressure = sm.collect_observations_as_time_series(res, 'pressure')
        CO2saturation = sm.collect_observations_as_time_series(res, 'CO2saturation')
        pressureAve = sm.collect_observations_as_time_series(res, 'pressureAve')

        # True values
        true_pressure = [30517095., 31336671., 30942776., 30729670., 30600297.]
        true_saturation = [0.0, 0.76610985, 1.0, 1.0, 1.0]
        true_pressureAve = [30670710., 31426544., 31013189., 30800083., 30670710.]

        # Test with helpful message!
        for tp, p, t in zip(true_pressure, pressure, time_array):
            self.assertTrue(abs((tp-p)/tp) < 0.01,
                            'Pressure at time t={} days is {} Pa but should be {} Pa'
                            .format(str(t), str(p), str(tp)))
        for ts, s, t in zip(true_saturation[1:], CO2saturation[1:], time_array[1:]):
            self.assertTrue(abs((ts-s)/ts) < 0.01,
                            'Saturation at time t={} days is {} but should be {}'
                            .format(str(t), str(s), str(ts)))
        for tx, q, t in zip(true_pressureAve, pressureAve, time_array):
            self.assertTrue(abs((tx-q)/tx) < 0.01,
                            'Pressure at time t={} days is {} Pa but should be {} Pa'
                            .format(str(t), str(q), str(tx)))

    @unittest.skipIf(AtmosphericROM.model_data_check() == False,
                     "Atmospheric ROM component library file is not present.")
    def test_atm(self):
        """Tests atmospheric ROM component.

        Tests the forward system model with only an atmospheric model component driven
        by known fixed flow rates against known data.
        """
        time_array = 365.25*np.arange(0, 7.0)
        sm_model_kwargs = {'time_array': time_array}

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add atmopsheric rom component
        atm = sm.add_component_model_object(AtmosphericROM(name='atm', parent=sm))
        # Add parameters
        atm.add_par('T_amb', value=25.0)
        atm.add_par('P_amb', value=0.9869233)
        atm.add_par('V_wind', value=5.0)
        atm.add_par('C0_critical', value=0.01)
        atm.add_par('T_source', value=25.0)

        # Setup model keyword arguments
        atm.model_kwargs['x_receptor'] = [820.0]
        atm.model_kwargs['y_receptor'] = [1010.0]
        atm.model_kwargs['x_coor'] = [900.0, 1000.0]
        atm.model_kwargs['y_coor'] = [1000.0, 1050.0]

        # Add dynamic input parameters
        atm.add_dynamic_kwarg('co2_leakrate', [[0.0, 0.0],
                                               [0.5, 1.0],
                                               [2.0, 3.0],
                                               [4.5, 5.0],
                                               [6.0, 7.0],
                                               [8.0, 9.0],
                                               [10.0, 10.0]])
        # Add observations
        atm.add_obs('num_sources')
        atm.add_obs('x_new_s000')
        atm.add_obs('y_new_s000')
        atm.add_obs('critical_distance_s000')
        atm.add_obs('x_new_s001')
        atm.add_obs('y_new_s001')
        atm.add_obs('critical_distance_s001')
        atm.add_obs('outflag_r000')

        # Run model
        sm.forward()

        num_sources = sm.collect_observations_as_time_series(atm, 'num_sources')
        true_sources = [2, 2, 1, 1, 1, 1, 1]
        for tv, sv, time in zip(true_sources, num_sources, time_array):
            self.assertTrue(abs((tv-sv)) < 0.01,
                            'Number of sources at time {} is {} but should be {}'
                            .format(str(time), str(sv), str(tv)))
        x_new_0 = sm.collect_observations_as_time_series(atm, 'x_new_s000')
        y_new_0 = sm.collect_observations_as_time_series(atm, 'y_new_s000')
        crit_dist_0 = sm.collect_observations_as_time_series(atm, 'critical_distance_s000')
        x_true = [900., 900., 956.03217926, 951.57461885, 952.30292214,
                  951.76012878, 950.]
        y_true = [1000., 1000., 1028.01608963, 1025.78730943, 1026.15146107,
                  1025.88006439, 1025.]
        crit_dist_true = [0., 50.01595939, 198.20226211, 291.59765037, 341.10956199,
                          390.0737076, 423.09433987]
        test_now = [(x_true, x_new_0), (y_true, y_new_0), (crit_dist_true, crit_dist_0)]
        for true_vals, sim_vals in test_now:
            for tv, sv, time in zip(true_vals, sim_vals, time_array):
                self.assertTrue(abs((tv-sv)) < 0.01,
                                'Atmospheric ROM output at time {} is {} but should be {}'
                                .format(str(time), str(sv), str(tv)))
        out_flag = sm.collect_observations_as_time_series(atm, 'outflag_r000')
        true_flag = [0, 0, 1, 2, 2, 2, 2]
        for tv, sv, time in zip(true_flag, out_flag, time_array):
            self.assertTrue(abs((tv-sv)) < 0.01,
                            'Atmospheric ROM out_flag at time {} is {} but should be {}'
                            .format(str(time), str(sv), str(tv)))

    @unittest.skipIf(CarbonateAquifer.model_data_check() == False,
                     "Carbonate Aquifer component library file is not present.")
    def test_carb_aquifer(self):
        """Tests the carbonate aquifer component.

        Tests the forward system model with only an aquifer component driven
        by known fixed flow rates against known data.
        """
        # Create system model
        time_array = 365.25*np.arange(0.0, 2.0)
        sm_model_kwargs = {'time_array': time_array}  # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add carbonate aquifer model object and define parameters
        ca = sm.add_component_model_object(CarbonateAquifer(name='ca', parent=sm))
        ca.add_par('perm_var', value=1.24e+00)
        ca.add_par('corr_len', value=2.01E+00)
        ca.add_par('aniso', value=4.36e+01)
        ca.add_par('mean_perm', value=-1.14E+01)
        ca.add_par('aqu_thick', value=4.71E+02)
        ca.add_par('hyd_grad', value=1.63e-02)
        ca.add_par('calcite_ssa', value=4.13E-03)
        ca.add_par('cl', value=4.59E+00)
        ca.add_par('organic_carbon', value=5.15e-03)
        ca.add_par('benzene_kd', value=1.53e+00)
        ca.add_par('nap_kd', value=2.78e+00)
        ca.add_par('phenol_kd', value=1.29e+00)
        ca.add_par('benzene_decay', value=1.58e+00)
        ca.add_par('nap_decay', value=3.65e-01)
        ca.add_par('phenol_decay', value=2.61e-01)
        ca.model_kwargs['x'] = [0.]
        ca.model_kwargs['y'] = [0.]
        ca.model_kwargs['co2_rate'] = [1.90e-02]  # kg/s
        ca.model_kwargs['co2_mass'] = [6.00e+05]  # kg
        ca.model_kwargs['brine_rate'] = [4.62e-03]  # kg/s
        ca.model_kwargs['brine_mass'] = [1.46e+05]  # kg

        # Add observations (output) from the carbonate aquifer model
        ca.add_obs('TDS_volume')

        # Run system model using current values of its parameters
        sm.forward()

        # Assign observations of the model to variables
        tds = sm.collect_observations_as_time_series(ca, 'TDS_volume')

        # True values of observations at time 0
        true_tds = [18940957.051663853, 18925088.43222359]

        # Test with helpful message!
        for ts, s, t in zip(true_tds, tds, time_array):
            self.assertTrue(abs((ts-s)/ts) < 0.01,
                            'TDS at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))

    def test_cemented_wellbore_forward(self):
        """Tests cemented wellbore component.

        Tests the system model with a cemented wellbore component provided
        with dynamic input of pressure and CO2 saturation against expected output
        for 10 years of data.
        """
        # Create system model
        # Define keyword arguments of the system model
        num_years = 10
        time_array = 365.25*np.arange(0.0, num_years+1)
        sm_model_kwargs = {'time_array': time_array}   # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Define dynamic input
        pressure = np.array([10516765., 15017301.38843215, 15290401.37235173,
                             15450154.92651543, 15563501.96157327, 15651420.7877774,
                             15723255.7175046, 15783991.29577689, 15836602.85344631,
                             15883009.54069186, 15924521.74018078])
        saturation = np.array([0., 0.21384356, 0.30898839, 0.38199558, 0.44354349,
                               0.49776829, 0.54679126, 0.59187253, 0.63383315,
                               0.67324343, 0.71051859])

        # Add cemented wellbore component
        cw = sm.add_component_model_object(CementedWellbore(name='cw', parent=sm))

        # Add parameters of cemented wellbore component
        cw.add_par('logWellPerm', value=-12.5, vary=False)
        cw.add_par('logThiefPerm', value=-13.0, vary=False)
        cw.add_par('wellDepth', value=1062.0, vary=False)
        cw.add_par('depthRatio', value=0.48607, vary=False)

        # Add dynamic input of cemented wellbore component
        cw.add_dynamic_kwarg('pressure', pressure)
        cw.add_dynamic_kwarg('CO2saturation', saturation)

        # Add observations of the cemented wellbore components
        obs_names = ['CO2_aquifer1', 'CO2_aquifer2', 'brine_aquifer1', 'brine_aquifer2']
        for key in obs_names:
            cw.add_obs(key)

        # Run system model
        sm.forward()

        true_vals = {}
        true_vals['CO2_aquifer1'] = np.array([
            0., 0.00013478, 0.00017323, 0.00020062, 0.00022284, 0.00024191,
            0.00025884, 0.00027419, 0.00028831, 0.00030145, 0.0003138])
        true_vals['CO2_aquifer2'] = np.array([
            0., 2.32835622e-05, 2.70753038e-05, 2.99848093e-05, 3.24376359e-05,
            3.45986197e-05, 3.65522984e-05, 3.83488914e-05, 4.00211192e-05,
            4.15917100e-05, 4.30772114e-05])
        true_vals['brine_aquifer1'] = np.array([
            0., 3.12597455e-05, 3.10693469e-05, 3.05445631e-05,
            2.99302301e-05, 2.92885762e-05, 2.86418491e-05, 2.79993213e-05,
            2.73651082e-05, 2.67410110e-05, 2.61277104e-05])
        true_vals['brine_aquifer2'] = np.array([
            0., 3.94622270e-06, 4.00488449e-06, 4.03921021e-06, 4.06356573e-06,
            4.08245772e-06, 4.09789380e-06, 4.11094493e-06, 4.12225041e-06,
            4.13222260e-06, 4.14114306e-06])

        collected_vals = {}
        for key in obs_names:
            collected_vals[key] = sm.collect_observations_as_time_series(cw, key)

        for key in obs_names:
            for true_v, col_v, t in zip(true_vals[key][1:],
                                        collected_vals[key][1:],
                                        time_array[1:]):
                self.assertTrue(
                    abs((col_v-true_v)/true_v) < 0.001,
                    'Observation {} at t={} days is {} kg/s but should be {} kg/s'\
                                .format(key, t, col_v, true_v))

    def test_cemented_wellbore_wr_forward(self):
        """Tests cemented wellbore (wr) component.

        Tests the system model with a cemented wellbore (wr) component provided
        with dynamic input of pressure and CO2 saturation against expected output
        for 10 years of data.
        """
        # Create system model
        # Define keyword arguments of the system model
        num_years = 10
        time_array = 365.25*np.arange(0.0, num_years+1)
        sm_model_kwargs = {'time_array': time_array}   # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Define dynamic input
        pressure = np.array([10516765., 15017301.38843215, 15290401.37235173,
                             15450154.92651543, 15563501.96157327, 15651420.7877774,
                             15723255.7175046, 15783991.29577689, 15836602.85344631,
                             15883009.54069186, 15924521.74018078])
        saturation = np.array([0., 0.21384356, 0.30898839, 0.38199558, 0.44354349,
                               0.49776829, 0.54679126, 0.59187253, 0.63383315,
                               0.67324343, 0.71051859])

        # Add cemented wellbore component
        cw = sm.add_component_model_object(CementedWellboreWR(name='cw', parent=sm))

        # Add parameters of cemented wellbore component
        cw.add_par('logWellPerm', value=-12.5, vary=False)
        cw.add_par('logThiefPerm', value=-13.0, vary=False)
        cw.add_par('wellDepth', value=1062.0, vary=False)
        cw.add_par('depthRatio', value=0.48607, vary=False)
        cw.add_par('aquiferThickness', value=22.4, vary=False)
        cw.add_par('thiefZoneThickness', value=19.2, vary=False)
        cw.add_par('reservoirThickness', value=51.2, vary=False)

        # Add dynamic input of cemented wellbore component
        cw.add_dynamic_kwarg('pressure', pressure)
        cw.add_dynamic_kwarg('CO2saturation', saturation)

        # Add observations of the cemented wellbore components
        obs_names = ['CO2_aquifer1', 'CO2_aquifer2', 'brine_aquifer1', 'brine_aquifer2']
        for key in obs_names:
            cw.add_obs(key)

        # Run system model
        sm.forward()

        true_vals = {}
        true_vals['CO2_aquifer1'] = np.array([
            0., 1.54549001e-06, 2.65823820e-06, 4.03646479e-06, 5.74332264e-06,
            7.83837685e-06, 1.03851095e-05, 1.34532151e-05, 1.71198220e-05,
            2.14704367e-05, 2.65998210e-05])
        true_vals['CO2_aquifer2'] = np.array([
            0., 8.12962701e-09, 6.87638942e-09, 6.28456749e-09, 5.89702945e-09,
            5.61164683e-09, 5.38726758e-09, 5.20324338e-09, 5.04778528e-09,
            4.91355077e-09, 4.79566908e-09])
        true_vals['brine_aquifer1'] = np.array([
            0., 8.67330716e-07, 9.75758547e-07, 1.09362038e-06, 1.20339979e-06,
            1.30886292e-06, 1.41190025e-06, 1.51361684e-06, 1.61472437e-06,
            1.71571347e-06, 1.81694008e-06])
        true_vals['brine_aquifer2'] = np.array([
            0., 1.87002640e-07, 1.35664923e-07, 1.06362201e-07, 8.66162656e-08,
            7.24389716e-08, 6.18977679e-08, 5.37264563e-08, 4.71992092e-08,
            4.18652791e-08, 3.74278144e-08])

        collected_vals = {}
        for key in obs_names:
            collected_vals[key] = sm.collect_observations_as_time_series(cw, key)

        for key in obs_names:
            for true_v, col_v, t in zip(true_vals[key][1:],
                                        collected_vals[key][1:],
                                        time_array[1:]):
                self.assertTrue(
                    abs((col_v-true_v)/true_v) < 0.001,
                    'Observation {} at t={} days is {} kg/s but should be {} kg/s'\
                                .format(key, t, col_v, true_v))

    def test_chemical_sealing_coupled_simple_reservoir(self):
        """Tests the chemical sealing component.

        Tests the system model with a chemical well sealing component coupled to a
        simple reservoir component.
        """
        # Create system model
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

        # Add reservoir component
        sres = sm.add_component_model_object(
            SimpleReservoir(name='sres', parent=sm, injX=0., injY=0.,
                            locX=50., locY=100.))

        # Add parameters of reservoir component model
        sres.add_par('numberOfShaleLayers', value=3, vary=False)
        sres.add_par('injRate', value=0.04, vary=False)
        sres.add_par('shale1Thickness', min=30.0, max=50., value=40.0)
        sres.add_par('aquifer1Thickness', min=20.0, max=60., value=50.0)
        sres.add_par('aquifer2Thickness', min=30.0, max=50., value=45.0)
        sres.add_par('logResPerm', min=-13., max=-11., value=-12.5)

        # Add pressure observation of reservoir component model
        sres.add_obs('pressure')

        # Run the reservoir model
        sm.forward()

        # Collect pressure data from the reservoir component
        pressure = sm.collect_observations_as_time_series(sres, 'pressure')
        # maxOverpressure = maximum increase in pressure experienced by the well
        maxOverpressure_value = np.amax(pressure)-pressure[0]

        # Add chemical well sealing component
        cws = sm2.add_component_model_object(ChemicalWellSealing(name='cws', parent=sm2))

        # Add parameters of chemical well sealing component
        fracAperture = 0.000150
        fracLength = 20.0
        cws.add_par('fractureAperture', value=fracAperture)
        cws.add_par('fractureLength', value=fracLength)
        cws.add_par('maxOverpressure', value=maxOverpressure_value)

        # Add observations (output) of the chemical well sealing component
        cws.add_obs('seal_flag', index=[0])
        cws.add_obs('seal_time', index=[0])

        # Run the system model
        sm2.forward()

        # True values: simulation results
        true_results = [0, 0]

        # Test with helpful message
        self.assertTrue(abs((true_results[0]-sm2.obs['cws.seal_flag_0'].sim)) < 0.01,
                        'Sealing flag is {} but should be {}'.format(
                            true_results[0], sm2.obs['cws.seal_flag_0'].sim))
        self.assertTrue(abs((true_results[1]-sm2.obs['cws.seal_time_0'].sim)) < 0.01,
                        'Sealing time is {} but should be {}'.format(
                            true_results[1], sm2.obs['cws.seal_time_0'].sim))

    def test_chemical_sealing_not_seal_forward(self):
        """Tests the chemical sealing component.

        Tests the system model with a chemical well sealing component for the case
        when a fracture pathway is not expected to seal.
        """
        # Create system model
        time_array = 365.25*np.arange(0.0, 2.0)
        sm_model_kwargs = {'time_array': time_array} # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add chemical well sealing component and define parameters
        cws = sm.add_component_model_object(ChemicalWellSealing(name='cws', parent=sm))
        cws.add_par('fractureAperture', value=0.00066)
        cws.add_par('fractureLength', value=25)
        cws.add_par('maxOverpressure', value=6800000)

        # Add observations (output) of the chemical well sealing component
        cws.add_obs('seal_flag', index=[0])
        cws.add_obs('seal_time', index=[0])

        # Run the system model
        sm.forward()

        # True values: simulation results
        true_results = [0, 0]

        # Test with helpful message
        self.assertTrue(abs((true_results[0]-sm.obs['cws.seal_flag_0'].sim)) < 0.01,
                        'Sealing flag is {} but should be {}'.format(
                            true_results[0], sm.obs['cws.seal_flag_0'].sim))
        self.assertTrue(abs((true_results[1]-sm.obs['cws.seal_time_0'].sim)) < 0.01,
                        'Sealing time is {} but should be {}'.format(
                            true_results[1], sm.obs['cws.seal_time_0'].sim))

    def test_chemical_sealing_seal_forward(self):
        """Tests the chemical sealing component.

        Tests the system model with a chemical well sealing component for the case
        when a fracture pathway is expected to seal.
        """
        # Create system model
        time_array = 365.25*np.arange(0.0, 2.0)
        sm_model_kwargs = {'time_array': time_array} # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add chemical well sealing component and define parameters
        cws = sm.add_component_model_object(ChemicalWellSealing(name='cws', parent=sm))
        cws.add_par('fractureAperture', value=0.000043)
        cws.add_par('fractureLength', value=160)
        cws.add_par('maxOverpressure', value=14000000)

        # Add observations (output) of the chemical well sealing component
        cws.add_obs('seal_flag', index=[0])
        cws.add_obs('seal_time', index=[0])

        # Run the system model
        sm.forward()

        # True values: simulation results
        true_results = [1, 691784.2222]

        # Test with helpful message
        self.assertTrue(
            abs((true_results[0]-sm.obs['cws.seal_flag_0'].sim)/true_results[0]) < 0.01,
            'Sealing flag is {} but should be {}'.format(
                true_results[0], sm.obs['cws.seal_flag_0'].sim))
        self.assertTrue(
            abs((true_results[1]-sm.obs['cws.seal_time_0'].sim)/true_results[1]) < 0.01,
            'Sealing time is {} but should be {}'.format(
                true_results[1], sm.obs['cws.seal_time_0'].sim))

    def test_coupled_analytical_reservoir_ms_wellbore_forward(self):
        """Tests coupling of reservoir and multisegmented wellbore components.

        Tests the system model with an analytical reservoir component
        coupled with a multi-segmented wellbore component against 5 years of data.
        """
        num_years = 5
        time_array = 365.25*np.arange(0.0, num_years+1)
        sm_model_kwargs = {'time_array': time_array}  # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add reservoir component
        res = sm.add_component_model_object(AnalyticalReservoir(
            name='res', parent=sm, injX=100., injY=0., locX=0., locY=0.))

        # Add parameters of reservoir component model
        res.add_par('injRate', value=0.0185)
        res.add_par('reservoirRadius', value=500.0)
        res.add_par('reservoirThickness', value=50.0)
        res.add_par('logResPerm', value=-13.69897)
        res.add_par('reservoirPorosity', value=0.15)
        res.add_par('aquiferThickness', value=30.0)
        res.add_par('shaleThickness', value=970.0)

        res.add_par('brineDensity', value=1045.0)
        res.add_par('CO2Density', value=479.0)
        res.add_par('brineViscosity', value=2.535e-4)
        res.add_par('CO2Viscosity', value=3.95e-5)
        res.add_par('brineResSaturation', value=0.0)
        res.add_par('brineCompressibility', value=1.e-11)

        res.add_par('datumPressure', value=101325.0)

        # Add observations of reservoir component model
        res.add_obs_to_be_linked('pressure')
        res.add_obs_to_be_linked('CO2saturation')

        res.add_obs('pressure')
        res.add_obs('CO2saturation')
        res.add_obs('mass_CO2_reservoir')

        # Add multisegmented wellbore component
        ms = sm.add_component_model_object(MultisegmentedWellbore(
            name='ms', parent=sm))
        ms.add_par('wellRadius', value=0.15)

        ms.add_par('numberOfShaleLayers', value=5)
        ms.add_par('shale1Thickness', value=100.0)
        ms.add_par('shale2Thickness', value=100.0)
        ms.add_par('shale3Thickness', value=100.0)
        ms.add_par('shale4Thickness', value=100.0)
        ms.add_par('shale5Thickness', value=2450.0)

        ms.add_par('logWell1Perm', value=-12)
        ms.add_par('logWell2Perm', value=-12)
        ms.add_par('logWell3Perm', value=-12)
        ms.add_par('logWell4Perm', value=-12)
        ms.add_par('logWell5Perm', value=-100)

        ms.add_par('aquifer1Thickness', value=30.0)
        ms.add_par('aquifer2Thickness', value=30.0)
        ms.add_par('aquifer3Thickness', value=30.0)
        ms.add_par('aquifer4Thickness', value=30.0)

        ms.add_par('logAqu1Perm', value=-12)
        ms.add_par('logAqu2Perm', value=-12)
        ms.add_par('logAqu3Perm', value=-12)
        ms.add_par('logAqu4Perm', value=-12)

        ms.add_par('brineResSatAquifer1', value=0.18)
        ms.add_par('brineResSatAquifer2', value=0.4)
        ms.add_par('brineResSatAquifer3', value=0.5)
        ms.add_par('brineResSatAquifer4', value=0.0)

        ms.add_par_linked_to_par('reservoirThickness',
                                 res.default_pars['reservoirThickness'])
        ms.add_par_linked_to_par('aqu0BrineResSaturation',
                                 res.default_pars['brineResSaturation']) # reservoir
        ms.add_par_linked_to_par('datumPressure',
                                 res.default_pars['datumPressure'])

        ms.add_kwarg_linked_to_obs('pressure', res.linkobs['pressure'])
        ms.add_kwarg_linked_to_obs('CO2saturation', res.linkobs['CO2saturation'])

        ms.add_obs('CO2_aquifer1')
        ms.add_obs('CO2_aquifer4')
        ms.add_obs('brine_aquifer1')
        ms.add_obs('brine_aquifer4')

        sm.forward()  # system model is run deterministically

        test_CO2_aquifer1 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer1')
        test_CO2_aquifer4 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer4')
        test_brine_aquifer1 = sm.collect_observations_as_time_series(ms, 'brine_aquifer1')
        test_brine_aquifer4 = sm.collect_observations_as_time_series(ms, 'brine_aquifer4')

        # True values
        true_CO2_aquifer1 = [0., 0.01390931, 0.01474191, 0.0162474, 0.01625035]
        true_CO2_aquifer4 = [0., 0., 0., 0., 0.00254786]
        true_brine_aquifer1 = [0., 3.04154360e-04, 5.99924958e-05, 0., 0.]
        true_brine_aquifer4 = [0., -1.02938960e-07, -1.09542806e-07,
                               -4.55770797e-06, 0.0e+00]

        # Test with helpful message!
        for tp, p, t in zip(true_CO2_aquifer1, test_CO2_aquifer1, time_array):
            self.assertTrue(abs((tp-p)) < 0.01,
                            'CO2_aquifer1 at time t={} days is {} kg/s but should be {} kg/s'
                            .format(str(t), str(p), str(tp)))

        for tp, p, t in zip(true_CO2_aquifer4, test_CO2_aquifer4, time_array):
            self.assertTrue(abs((tp-p)) < 0.01,
                            'CO2_aquifer4 at time t={} days is {} kg/s but should be {} kg/s'
                            .format(str(t), str(p), str(tp)))

        for tp, p, t in zip(true_brine_aquifer1, test_brine_aquifer1, time_array):
            self.assertTrue(abs((tp-p)) < 0.01,
                            'brine_aquifer1 at time t={} days is {} kg/s but should be {} kg/s'
                            .format(str(t), str(p), str(tp)))

        for tp, p, t in zip(true_brine_aquifer4, test_brine_aquifer4, time_array):
            self.assertTrue(abs((tp-p)) < 0.01,
                            'brine_aquifer4 at time t={} days is {} kg/s but should be {} kg/s'
                            .format(str(t), str(p), str(tp)))

    def test_coupled_components_with_grid_obs(self):
        """ Tests coupling of two components utilizing gridded observations.

        Tests the system model with two coupled components one of which has
        gridded observation against known data. Tests the content of the file created
        for the storage of gridded observation simulated values at one
        of the time points.
        """
        def comp_model1(p, x=None, y=None, time_point=0.0):
            if (x is not None) and (y is not None):
                xx = x
                yy = y
            else:
                x = np.linspace(0., 10., 6)
                y = np.linspace(0., 10., 6)

            if 'var1' in p:
                var1 = p['var1']
            else:
                var1 = 2.

            out = {'plane': np.zeros((len(x)*len(y),))}
            ind = 0
            for xv in xx:
                for yv in yy:
                    out['plane'][ind] = (
                        np.log10(time_point+1)*(xv+yv)/(10*var1) +
                        4*np.log10(time_point+1))
                    ind = ind + 1
            out['time_squared'] = 4*np.sin(time_point)**2/1.0e+3 + var1*time_point/10.
            return out

        def comp_model2(p, time_point=0.0, kw1=1.0, kw2=None, kw3=1.0, kw4=None):
            # Check which parameters are provided in the input dictionary
            if 'var1' in p:
                var1 = p['var1']
            else:
                var1 = 1.0  # default values if the desired parameter value was not passed

            if 'var2' in p:
                var2 = p['var2']
            else:
                var2 = 1.0

            if 'var3' in p:
                var3 = p['var3']
            else:
                var3 = 1.0

            # Check whether keyword arguments kw2 or kw4 are provided
            if kw2 is None:
                kw2 = [1.0]
            if kw4 is None:
                kw4 = [1.3, 1.5]

            # Define output of the function (output is also a dictionary)
            output = dict()
            output['output1'] = np.abs((kw1*var1**2 + kw3*var2*var3)*np.sin(time_point))
            output['output2'] = np.sum(kw2) + np.sum(kw4)

            # Component model should return a dictionary of outputs
            return output

        # Define keyword arguments of the system model
        num_years = 4
        time_array = 365.25*np.arange(1.0, num_years+1)
        sm_model_kwargs = {'time_array': time_array}  # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add component 1
        N = 100
        x = np.linspace(1., 10., N+1)

        cm1 = sm.add_component_model(
            name='cm1', model=comp_model1,
            model_kwargs={'x': x, 'y': x, 'time_point': 365.25})
        cm1.grid_obs_keys = ['plane']
        cm1.add_par('var1', min=2., max=5., value=3.)
        cm1.add_grid_obs(
            'plane', constr_type='array', output_dir='test_output', index=[3])
        cm1.add_local_obs('plane_loc1', 'plane', constr_type='array', loc_ind=65)
        cm1.add_obs('time_squared')
        cm1.add_obs_to_be_linked('time_squared')
        cm1.add_obs_to_be_linked('plane', obs_type='grid')

        # Add component 2
        cm2 = sm.add_component_model(name='cm2', model=comp_model2,
                                     model_kwargs={'time_point': 365.25})
        cm2.add_par('var1', value=1., vary=False)
        cm2.add_kwarg_linked_to_obs('kw1', cm1.linkobs['time_squared'])
        cm2.add_kwarg_linked_to_obs('kw2', cm1.linkobs['plane'], obs_type='grid')
        cm2.add_kwarg_linked_to_obs('kw3', cm1.linkobs['plane'], obs_type='grid',
                                    constr_type='array', loc_ind=[45])
        cm2.add_kwarg_linked_to_obs('kw4', cm1.linkobs['plane'], obs_type='grid',
                                    constr_type='array', loc_ind=[23, 33, 43, 53])
        cm2.add_obs('output1')
        sm.forward()

        # Check saved files
        filename = '_'.join(['cm1', 'plane', 'sim_0', 'time_3'])+'.npz'
        d = np.load(os.sep.join(['test_output', filename]))
        file_data = d['data'][10:15]
        d.close()
        true_file_data = [12.9657344, 12.97522925, 12.98472409, 12.99421893, 13.00371377]

        plane_loc1 = sm.collect_observations_as_time_series(cm1, 'plane_loc1')
        time_squared = sm.collect_observations_as_time_series(cm1, 'time_squared')
        true_plane_loc1 = [10.92596568, 12.20632674, 12.95592541, 13.48795072]
        true_time_squared = [109.57715925, 219.15397464, 328.7265263, 438.3001008]

        output1 = sm.collect_observations_as_time_series(cm2, 'output1')
        true_output1 = [88.42290161, 230.4544242, 210.95033575, 71.68746713]

        for t, v, tv in zip(time_array, plane_loc1, true_plane_loc1):
            self.assertTrue(abs((v-tv)/tv) < 0.01,
                            'Observation plane_loc1 at time {} is {} but should be {}'
                            .format(str(t), str(v), str(tv)))

        for t, v, tv in zip(time_array, time_squared, true_time_squared):
            self.assertTrue(abs((v-tv)/tv) < 0.01,
                            'Observation time_squared at time {} is {} but should be {}'
                            .format(str(t), str(v), str(tv)))

        for t, v, tv in zip(time_array, output1, true_output1):
            self.assertTrue(abs((v-tv)/tv) < 0.01,
                            'Observation output1 at time {} is {} but should be {}'
                            .format(str(t), str(v), str(tv)))

        for v, tv in zip(file_data, true_file_data):
            self.assertTrue(abs((v-tv)/tv) < 0.01,
                            'Observation written to file is {} but should be {}'
                            .format(str(v), str(tv)))

    def test_coupled_reservoir_cemented_wellbore_forward(self):
        """Tests coupling of reservoir and cemented wellbore components.

        Verifies flow rates for a system model with a simple reservoir
        component coupled to a cemented wellbore component against 5 years
        of data.

        .. image:: test_output/test_coupled_reservoir_cemented_wellbore_forward.png
            :align: center
            :width: 400

        """
        # Create system model
        num_years = 5
        time_array = 365.25*np.arange(0.0, num_years+1)
        sm_model_kwargs = {'time_array': time_array}  # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add reservoir component
        sres = sm.add_component_model_object(SimpleReservoir(name='sres',
                                                             parent=sm,
                                                             locX=550.,
                                                             locY=100.))

        # Add parameters of reservoir component model
        sres.add_par('numberOfShaleLayers', value=3, vary=False)
        sres.add_par('shale1Thickness', min=500.0, max=550., value=525.0)
        sres.add_par('shale2Thickness', min=450.0, max=500., value=475.0)
        # Shale 3 has a fixed thickness of 11.2 m
        sres.add_par('shale3Thickness', value=11.2, vary=False)
        # Aquifer 1 (thief zone has a fixed thickness of 22.4)
        sres.add_par('aquifer1Thickness', value=22.4, vary=False)
        # Aquifer 2 (shallow aquifer) has a fixed thickness of 19.2
        sres.add_par('aquifer2Thickness', value=19.2, vary=False)
        # Reservoir has a fixed thickness of 51.2
        sres.add_par('reservoirThickness', value=51.2, vary=False)

        # Add observations of reservoir component model
        sres.add_obs('pressure')
        sres.add_obs('CO2saturation')
        sres.add_obs_to_be_linked('pressure')
        sres.add_obs_to_be_linked('CO2saturation')

        # Add cemented wellbore component
        cw = sm.add_component_model_object(CementedWellbore(name='cw', parent=sm))

         # Add parameters of cemented wellbore component
        cw.add_par('logWellPerm', min=-13.9, max=-10.0, value=-10.1)
        cw.add_par('logThiefPerm', value=-12.00035, vary=False)

        # Add keyword arguments of the cemented wellbore component model
        cw.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
        cw.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

        cw.add_composite_par('wellDepth', expr='sres.shale1Thickness' +
                             '+sres.shale2Thickness + sres.shale3Thickness' +
                             '+sres.aquifer1Thickness+ sres.aquifer2Thickness')
        cw.add_composite_par('depthRatio',
                             expr='(sres.shale2Thickness+sres.shale3Thickness' +
                             '+ sres.aquifer2Thickness + sres.aquifer1Thickness/2)/cw.wellDepth')
        cw.add_composite_par('initPressure',
                             expr='sres.datumPressure + cw.wellDepth*cw.g*sres.brineDensity')

        # Add observations of the cemented wellbore component
        cw.add_obs('CO2_aquifer1')
        cw.add_obs('CO2_aquifer2')
        cw.add_obs('CO2_atm')
        cw.add_obs('brine_aquifer1')
        cw.add_obs('brine_aquifer2')
        sm.forward()

        # True values: simulation results for the last time point
        true_flowrates = [0.0007111896900558063, 0.0006156910515266846,
                          0.0001732564889916223, 0.011944097940789211,
                          2.6686475096907638e-05]

        # Test with helpful message!
        labels = ['CO2_aquifer1', 'CO2_aquifer2', 'CO2_atm', 'brine_aquifer1', 'brine_aquifer2']
        # Collect results as dictionary of lists for every observation in labels list
        results = dict()
        for label in labels:
            results[label] = sm.collect_observations_as_time_series(cw, label)

        # Create plot for documentation
        sim_flowrates = [res[num_years-1] for res in results.values()]
        f, ax = plt.subplots(1, figsize=(8, 6))
        ax.plot(true_flowrates, sim_flowrates, "o")
        ax.plot([0, max(max(sim_flowrates, true_flowrates))],
                [0, max(max(sim_flowrates, true_flowrates))], "--", lw=0.5, c="k")
        ax.set_xlabel("Reference Flowrates [kg/s]")
        ax.set_ylabel("ROM Flowrates [kg/s]")
        f.tight_layout()
        fname = os.path.join(
            "test_output", "test_coupled_reservoir_cemented_wellbore_forward.png")
        f.savefig(fname)
        copyfile(fname, os.path.join(
            os.path.join("..", "documentation", "qaqc", "source"), fname))

        for tf, l in zip(true_flowrates, labels):
            self.assertTrue(abs((tf-results[l][num_years-1])/tf) < 0.01,
                            'Flow rate {} is {} but should be {}'
                            .format(l, str(results[l][num_years-1]), str(tf)))

    def test_coupled_reservoir_ms_wellbore_forward(self):
        """Tests coupling of reservoir and multisegmented wellbore components.

        Tests the system model with a simple reservoir component
        coupled with a multi-segmented wellbore component against 5 years of data.
        """
        num_years = 5
        time_array = 365.25*np.arange(0.0, num_years+1)
        sm_model_kwargs = {'time_array': time_array}  # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add reservoir component
        sres = sm.add_component_model_object(SimpleReservoir(name='sres',
                                                             parent=sm,
                                                             locX=550.,
                                                             locY=100.))

        # Add parameters of reservoir component model
        sres.add_par('numberOfShaleLayers', value=3, vary=False)
        sres.add_par('shale1Thickness', min=30.0, max=150., value=45.0)

        # Add observations of reservoir component model
        sres.add_obs_to_be_linked('pressure')
        sres.add_obs_to_be_linked('CO2saturation')
        sres.add_obs('pressure')
        sres.add_obs('CO2saturation')

        # Add multisegmented wellbore component
        ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))
        ms.add_par('wellRadius', min=0.01, max=0.02, value=0.015)

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
        ms.add_obs('CO2_aquifer1')
        ms.add_obs('CO2_aquifer2')
        ms.add_obs('CO2_atm')
        ms.add_obs('brine_aquifer1')
        ms.add_obs('brine_aquifer2')

        sm.forward()

        # True values
        true_flowrates = [8.32098502e-06, 1.90750350e-07,
                          1.91397466e-07, 2.26169428e-06,
                          1.86168527e-09]

        # Test with helpful message!
        labels = ['CO2_aquifer1', 'CO2_aquifer2', 'CO2_atm', 'brine_aquifer1', 'brine_aquifer2']
        # Collect results as dictionary of lists for every observation in labels list
        results = dict()
        for label in labels:
            results[label] = sm.collect_observations_as_time_series(ms, label)

        for tf, l in zip(true_flowrates, labels):
            self.assertTrue((tf-results[l][num_years-1] == 0.0) or
                            abs((tf-results[l][num_years-1])/tf) < 0.01,
                            'Flow rate {} is {} but should be {}'
                            .format(l, str(results[l][num_years-1]), str(tf)))

    @unittest.skipIf(CarbonateAquifer.model_data_check() == False,
                     "Carbonate Aquifer component library file is not present.")
    def test_coupled_reservoir_open_carbonate_forward(self):
        """Tests coupling of reservoir, wellbore and aquifer components.

        Tests the system model with coupled simple reservoir, open wellbore,
        adapter and carbonate aquifer components against the expected input.
        """
        # Define keyword arguments of the system model
        num_years = 5
        time_array = 365.25*np.arange(0.0, num_years+1)
        sm_model_kwargs = {'time_array': time_array}  # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add reservoir component
        sres = sm.add_component_model_object(SimpleReservoir(name='sres',
                                                             parent=sm,
                                                             locX=225.,
                                                             locY=100.))

        # Add parameters of reservoir component model
        sres.add_par('numberOfShaleLayers', value=3, vary=False)
        sres.add_par('shale1Thickness', value=620.0, vary=False)
        sres.add_par('shale2Thickness', value=510.0, vary=False)
        sres.add_par('shale3Thickness', value=400., vary=False)
        sres.add_par('aquifer1Thickness', value=10., vary=False)
        sres.add_par('aquifer2Thickness', value=80, vary=False)
        sres.add_par('reservoirThickness', value=55.0, vary=False)
        sres.add_par('injRate', value=0.0045, vary=False)
        sres.add_par('CO2Density', value=550, vary=False)

        # Add observations of reservoir component model
        sres.add_obs('pressure')
        sres.add_obs('CO2saturation')
        sres.add_obs_to_be_linked('pressure')
        sres.add_obs_to_be_linked('CO2saturation')

        # Add open wellbore component
        ow = sm.add_component_model_object(OpenWellbore(name='ow', parent=sm))

        # Add parameters of open wellbore component
        ow.add_par('wellRadius', min=0.025, max=0.05, value=0.025)
        ow.add_par('logReservoirTransmissivity', min=-11.2, max=-9.0, value=-11.2)
        ow.add_par('logAquiferTransmissivity', min=-11.2, max=-9.0, value=-11.2)
        ow.add_par('brineSalinity', value=0.004, vary=False)

        # Add keyword arguments of the open wellbore component model
        ow.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
        ow.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

        # Add composite parameter of open wellbore component
        ow.add_composite_par('reservoirDepth',
                             expr='+'.join(['sres.shale1Thickness',
                                            'sres.shale2Thickness',
                                            'sres.shale3Thickness',
                                            'sres.aquifer1Thickness',
                                            'sres.aquifer2Thickness']))

        ow.add_composite_par(
            'wellTop', expr='+'.join(['sres.shale3Thickness',
                                      'sres.aquifer2Thickness']))

        # Add observations of multisegmented wellbore component model
        ow.add_obs_to_be_linked('CO2_aquifer')
        ow.add_obs_to_be_linked('brine_aquifer')
        ow.add_obs_to_be_linked('brine_atm')
        ow.add_obs_to_be_linked('CO2_atm')
        ow.add_obs('brine_aquifer')
        ow.add_obs('CO2_atm')  # zero since well top is in aquifer
        ow.add_obs('CO2_aquifer')
        ow.add_obs('brine_atm')  # zero since well top is in aquifer

        # Add adapter that transforms leakage rates to accumulated mass
        adapt = sm.add_component_model_object(RateToMassAdapter(name='adapt', parent=sm))
        adapt.add_kwarg_linked_to_collection('CO2_aquifer',
                                             [ow.linkobs['CO2_aquifer'],
                                              ow.linkobs['CO2_atm']])
        adapt.add_kwarg_linked_to_collection('brine_aquifer',
                                             [ow.linkobs['brine_aquifer'],
                                              ow.linkobs['brine_atm']])
        adapt.add_obs_to_be_linked('mass_CO2_aquifer')
        adapt.add_obs_to_be_linked('mass_brine_aquifer')
        adapt.add_obs('mass_CO2_aquifer')
        adapt.add_obs('mass_brine_aquifer')

        # Add carbonate aquifer model object and define parameters
        ca = sm.add_component_model_object(CarbonateAquifer(name='ca', parent=sm))
        ca.add_par('perm_var', min=0.017, max=1.89, value=0.9535)
        ca.add_par('corr_len', min=1.0, max=3.95, value=2.475)
        ca.add_par('aniso', min=1.1, max=49.1, value=25.1)
        ca.add_par('mean_perm', min=-13.8, max=-10.3, value=-12.05)
        ca.add_par('aqu_thick', min=100., max=500., value=300.)
        ca.add_par('hyd_grad', min=2.88e-4, max=1.89e-2, value=9.59e-3)
        ca.add_par('calcite_ssa', min=0, max=1.e-2, value=5.5e-03)
        ca.add_par('organic_carbon', min=0, max=1.e-2, value=5.5e-03)
        ca.add_par('benzene_kd', min=1.49, max=1.73, value=1.61)
        ca.add_par('benzene_decay', min=0.15, max=2.84, value=1.5)
        ca.add_par('nap_kd', min=2.78, max=3.18, value=2.98)
        ca.add_par('nap_decay', min=-0.85, max=2.04, value=0.595)
        ca.add_par('phenol_kd', min=1.21, max=1.48, value=1.35)
        ca.add_par('phenol_decay', min=-1.22, max=2.06, value=0.42)
        ca.add_par('cl', min=0.1, max=6.025, value=0.776)
        ca.model_kwargs['x'] = [100.]
        ca.model_kwargs['y'] = [100.]

        CO2_rate_obs_list = []
        brine_rate_obs_list = []
        CO2_mass_obs_list = []
        brine_mass_obs_list = []
        CO2_rate_obs_list.append(ow.linkobs['CO2_aquifer'])
        brine_rate_obs_list.append(ow.linkobs['brine_aquifer'])
        CO2_mass_obs_list.append(adapt.linkobs['mass_CO2_aquifer'])
        brine_mass_obs_list.append(adapt.linkobs['mass_brine_aquifer'])

        # Add aquifer component's keyword argument co2_rate linked to the collection created above
        ca.add_kwarg_linked_to_collection('co2_rate', CO2_rate_obs_list)
        ca.add_kwarg_linked_to_collection('brine_rate', brine_rate_obs_list)
        ca.add_kwarg_linked_to_collection('co2_mass', CO2_mass_obs_list)
        ca.add_kwarg_linked_to_collection('brine_mass', brine_mass_obs_list)

        # Add observations (output) from the carbonate aquifer model
        ca.add_obs('pH_volume')
        ca.add_obs('TDS_volume')

        # Run system model using current values of its parameters
        sm.forward()  # system model is run deterministically

        tds = sm.collect_observations_as_time_series(ca, 'TDS_volume')
        ph = sm.collect_observations_as_time_series(ca, 'pH_volume')

        # True values
        true_tds = [0., 0., 0., 0., 284241.39903277, 2278197.61547477]
        true_ph = [0., 0., 0., 0., 0., 946552.25898049]

        # Suppress Numpy FutureWarning
        warnings.simplefilter(action='ignore', category=FutureWarning)

        # Test with helpful message
        for ts, s, t in zip(true_tds, tds, time_array):
            self.assertTrue((ts-s == 0.0) or abs((ts-s)/ts) < 0.01,
                            'TDS at time {} is {} but should be {}'.format(
                                str(t), str(s), str(ts)))

        for tp, p, t in zip(true_ph, ph, time_array):
            self.assertTrue((tp-p == 0.0) or abs((tp-p)/tp) < 0.01,
                            'pH at time {} is {} but should be {}'.format(
                                str(t), str(p), str(tp)))

    def test_coupled_reservoir_open_wellbore_forward(self):
        """Tests coupling of reservoir and open wellbore components.

        Tests the system model with a simple reservoir component coupled
        to an open wellbore component against 5 years of data.
        """
        # Create system model
        num_years = 5
        time_array = 365.25*np.arange(0.0, num_years+1)
        sm_model_kwargs = {'time_array': time_array}  # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add reservoir component
        sres = sm.add_component_model_object(SimpleReservoir(name='sres',
                                                             parent=sm,
                                                             locX=550.,
                                                             locY=100.))

        # Add parameters of reservoir component model
        sres.add_par('numberOfShaleLayers', value=3, vary=False)
        sres.add_par('shale1Thickness', min=500.0, max=550., value=525.0)
        sres.add_par('shale2Thickness', min=510.0, max=550., value=525.0)
        sres.add_par('shale3Thickness', min=10.0, max=2.0, value=1.5)
        sres.add_par('aquifer1Thickness', value=25.0, vary=False)
        sres.add_par('aquifer2Thickness', value=200.0, vary=False)
        sres.add_par('reservoirThickness', value=60.0, vary=False)

        # Add observations of reservoir component model
        sres.add_obs('pressure')
        sres.add_obs('CO2saturation')
        sres.add_obs_to_be_linked('pressure')
        sres.add_obs_to_be_linked('CO2saturation')

        # Add open wellbore component
        ow = sm.add_component_model_object(OpenWellbore(name='ow', parent=sm))

        # Add parameters of open wellbore component
        ow.add_par('logReservoirTransmissivity', value=-11.0, vary=False)
        ow.add_par('logAquiferTransmissivity', value=-11.0, vary=False)
        ow.add_par('brineSalinity', value=0.0, vary=False)

        # Add keyword arguments of the open wellbore component model
        ow.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
        ow.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

        # Add composite parameter of open wellbore component
        ow.add_composite_par(
            'reservoirDepth', expr='+'.join(['sres.shale1Thickness',
                                             'sres.shale2Thickness',
                                             'sres.shale3Thickness',
                                             'sres.aquifer1Thickness',
                                             'sres.aquifer2Thickness']))
        ow.add_composite_par(
            'wellTop', expr='sres.shale3Thickness + sres.aquifer2Thickness')

        # Add observations of the open wellbore component
        ow.add_obs('CO2_aquifer')
        ow.add_obs('brine_aquifer')

        sm.forward()

        # True values: simulation results for the last time point
        true_flowrates = [0.47971159914465, 3.0733673108455997]

        # Test with helpful message!
        labels = ['CO2_aquifer', 'brine_aquifer']
        for tf, l in zip(true_flowrates, labels):
            # The observations at the last time point can be accessed directly,
            # e.g. as ow.obs['CO2_aquifer_5'].sim (ow.obs[l+'_'+str(num_years)].sim) or
            # as sm.collect_observations_as_time_series(ow,'CO2_aquifer')[num_years]
            self.assertTrue(abs((tf-ow.obs[l+'_'+str(num_years)].sim)/tf) < 0.01,
                            'Flow rate {} is {} but should be {}'
                            .format(l, str(ow.obs[l+'_'+str(num_years)].sim), str(tf)))

    def test_deep_alluvium_aquifer(self):
        """Tests basic deep alluvium aquifer.

        Tests a simple deep alluvium aquifer component in a
        forward model against expected output for 2 years of data.
        """

        time_array = 365.25*np.arange(0.0, 2.0)
        sm_model_kwargs = {'time_array': time_array} # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add deep alluvium aquifer model object and define parameters
        daa = sm.add_component_model_object(
            DeepAlluviumAquifer(name='daa', parent=sm))

        daa.add_par('logK_sand1', value=-11.92098495)
        daa.add_par('logK_sand2', value=-11.7198002)
        daa.add_par('logK_sand3', value=-11.70137252)
        daa.add_par('logK_caprock', value=-15.69758676)
        daa.add_par('correlationLengthX', value=1098.994284)
        daa.add_par('correlationLengthZ', value=79.8062668)
        daa.add_par('sandFraction', value=0.800121364)
        daa.add_par('groundwater_gradient', value=0.001333374)
        daa.add_par('leak_depth', value=885.5060281)

        daa.model_kwargs['brine_rate'] = 3.21903E-05 # kg/s
        daa.model_kwargs['brine_mass'] = 10**4.71081307 # kg
        daa.model_kwargs['co2_rate'] = 0.060985038 # kg/s
        daa.model_kwargs['co2_mass'] = 10**6.737803184 # kg

        # Add observations (output) for the deep alluvium aquifer model
        daa.add_obs('TDS_volume')
        daa.add_obs('Pressure_dx')
        daa.add_obs('pH_dy')

        # Run the system model
        sm.forward()

        # Assign observations of the model to variables
        tds = sm.collect_observations_as_time_series(daa, 'TDS_volume')
        pressdx = sm.collect_observations_as_time_series(daa, 'Pressure_dx')
        phdy = sm.collect_observations_as_time_series(daa, 'pH_dy')

        # True values
        true_tds = [0.0, 2436214.13519926]
        true_press = [0.0, 176.71076959]
        true_pH = [0.0, 88.08772333]

        # Test with helpful message!
        for ts, s, t in zip(true_tds, tds, time_array):
            self.assertTrue(abs((ts-s)) < 1.0,
                            'TDS volume at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))

        for ts, s, t in zip(true_press, pressdx, time_array):
            self.assertTrue(abs((ts-s)) < 0.01,
                            'Pressure dx at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))
        for ts, s, t in zip(true_pH, phdy, time_array):
            self.assertTrue(abs((ts-s)) < 0.01,
                            'pH dy at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))

    def test_deep_alluvium_aquifer_ml(self):
        """Tests a deep alluvium aquifer component utilizing ML model.

        Tests a deep alluvium aquifer component in a forward simulation
        against expected output for 4 years of data.
        """
        # Create system model
        time_array = 365.25*np.arange(0.0, 4.0)
        sm_model_kwargs = {'time_array': time_array} # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add deep alluvium aquifer model object and define parameters
        daaml = sm.add_component_model_object(DeepAlluviumAquiferML(name='daaml', parent=sm))
        daaml.add_par('logK_sand1', value=-11.3845)
        daaml.add_par('logK_sand2', value=-11.9252)
        daaml.add_par('logK_sand3', value=-10.8862)
        daaml.add_par('logK_caprock', value=-14.731)
        daaml.add_par('correlationLengthX', value=520.721)
        daaml.add_par('correlationLengthZ', value=112.442)
        daaml.add_par('sandFraction', value=0.743)
        daaml.add_par('groundwater_gradient', value=0.00105408)
        daaml.add_par('leak_depth', value=715.99)

        # Add model keyword arguments
        daaml.model_kwargs['brine_rate'] = 6.24E-05    # kg/s
        daaml.model_kwargs['brine_mass'] = 2646900     # 10**6.422737534  kg
        daaml.model_kwargs['co2_rate'] = 0.150771      # kg/s
        daaml.model_kwargs['co2_mass'] = 171604000.1   # 10**8.234527407 kg

        # Add observations (output) of the deep alluvium aquifer model
        daaml.add_obs('TDS_volume')
        daaml.add_obs('TDS_dx')
        daaml.add_obs('TDS_dy')
        daaml.add_obs('TDS_dz')

        daaml.add_obs('Pressure_volume')
        daaml.add_obs('Pressure_dx')
        daaml.add_obs('Pressure_dy')
        daaml.add_obs('Pressure_dz')

        daaml.add_obs('pH_volume')
        daaml.add_obs('pH_dx')
        daaml.add_obs('pH_dy')
        daaml.add_obs('pH_dz')

        # Run the system model
        sm.forward()

        # True values
        true_vals = {}
        true_vals['TDS_volume'] = [0., 1602563.55016136,
                                   5438865.18640179, 8059917.15962449]
        true_vals['TDS_dx'] = [0., 65.69592945, 122.73450292, 213.45443131]
        true_vals['TDS_dy'] = [0., 182.430035, 225.89221453, 252.38214608]
        true_vals['TDS_dz'] = [0., 46.41604274, 104.95408265, 161.20099329]
        true_vals['Pressure_volume'] = [0., 14855869.2546651,
                                        11719586.45504402, 10645598.30293803]
        true_vals['Pressure_dx'] = [0., 278.70072389, 335.24981865, 379.97843937]
        true_vals['Pressure_dy'] = [0., 269.11794614, 302.68000236, 314.19433156]
        true_vals['Pressure_dz'] = [0., 348.06364686, 370.36887966, 389.19799003]
        true_vals['pH_volume'] = [0., 26076001.86850136,
                                  54044350.52776746, 67349593.38404109]
        true_vals['pH_dx'] = [0., 519.86851253, 656.91480486, 664.48539853]
        true_vals['pH_dy'] = [0., 362.15100201, 523.36799963, 582.61274082]
        true_vals['pH_dz'] = [0., 842.09423191, 1037.53341233, 1082.90360577]

        indices = [1, 2, 3]
        for key, vals in true_vals.items():
            for ind in indices:
                true_val = vals[ind]
                sim_val = daaml.obs[key+'_{}'.format(ind)].sim
                self.assertTrue(
                    abs(true_val-sim_val)/true_val < 0.01,
                    'Observation {} at {} days is {} but should be {}.'.format(
                        key, time_array[ind], sim_val, true_val))

    def test_fault_flow(self):
        """Tests FaultFlow component.

        Tests FaultFlow component in a forward model against expected output
        for 4 years of data.
        """
        # Define keyword arguments of the system model.
        # 4 time points
        time_array = 365.25 * np.array([
            0.00E+00, 1.00E+00, 2.00E+00, 4.00E+00])
        sm_model_kwargs = {'time_array': time_array}     # time is given in days

        # Create system model.
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        pressure = np.array([
            3.035851266631769761e+07, 3.035851266631769761e+07,
            3.035858161389059946e+07, 3.036002951292150095e+07]).reshape(4, 1)

        CO2saturation = np.array([
            2.000000000000000042e-03, 2.000000000000000042e-03,
            2.000000000000000042e-03, 1.999920000000000083e-03]).reshape(4, 1)

        ffc = sm.add_component_model_object(FaultFlow(name='ffc', parent=sm))
        ffc.add_dynamic_kwarg('pressure', pressure)
        ffc.add_dynamic_kwarg('CO2saturation', CO2saturation)

        # Add parameters of the FaultFlow component
        ffc.add_par('strike', value=95.4, vary=False)
        ffc.add_par('dip', value=90.0, vary=False)
        ffc.add_par('length', value=100.0, vary=False)
        ffc.add_par('xStart', value=388618.7, vary=False)
        ffc.add_par('yStart', value=3436341.1, vary=False)
        ffc.add_par('nSegments', value=1, vary=False)
        ffc.add_par('SGR', value=0.0, vary=False)
        ffc.add_par('stateVariable', value=1.0, vary=False)

        ffc.model_kwargs['aperture'] = [1.50E-05]

        ffc.add_par('aquiferDepth', value=500.0, vary=False)
        ffc.add_par('injectDepth', value=2884.31, vary=False)
        ffc.add_par('aquiferPressure', value=4.9e+06, vary=False)
        ffc.add_par('fieldPressure', value=2.85e+7, vary=False)
        ffc.add_par('injectPressure', value=3.0778e+7, vary=False)
        ffc.add_par('finalPressure', value=3.90e+7, vary=False)
        ffc.add_par('aquiferTemperature', value=30.0, vary=False)
        ffc.add_par('injectTemperature', value=98.9, vary=False)
        ffc.add_par('injectX', value=388505.9, vary=False)
        ffc.add_par('injectY', value=3434629.9, vary=False)

        ffc.add_par('salinity', value=50000.0, vary=False)
        ffc.add_par('CO2Density', value=430.0, vary=False)
        ffc.add_par('CO2Viscosity', value=3.72e-5, vary=False)
        ffc.add_par('brineDensity', value=988.00, vary=False)
        ffc.add_par('brineViscosity', value=4.36e-4, vary=False)
        ffc.add_par('CO2Solubility', value=2.0e-3, vary=False)

        ffc.add_par('brineResSaturation', value=0.15, vary=False)
        ffc.add_par('CO2ResSaturation', value=0.0, vary=False)
        ffc.model_kwargs['relativeModel'] = 'BC'
        ffc.add_par('permRatio', value=0.6, vary=False)
        ffc.add_par('entryPressure', value=1.0e+05, vary=False)
        ffc.add_par('lambda', value=2.5, vary=False)

        ffc.add_par('maxHorizontal', value=3.0e+07, vary=False)
        ffc.add_par('minHorizontal', value=2.0e+07, vary=False)
        ffc.add_par('maxTrend', value=55.0, vary=False)

        # Add scalar observations
        ffc.add_obs('CO2_aquifer_total')
        ffc.add_obs('brine_aquifer_total')

        # Run model
        sm.forward()

        # Collect scalar observations
        CO2_aquifer_total = sm.collect_observations_as_time_series(
            ffc, 'CO2_aquifer_total')
        brine_aquifer_total = sm.collect_observations_as_time_series(
            ffc, 'brine_aquifer_total')

        true_CO2_aquifer_total = [3.25610295e-11, 3.25610295e-11,
                                  3.25611982e-11, 3.25608339e-11]
        true_brine_aquifer_total = [0.00012755, 0.00012755,
                                    0.00012756, 0.0001276]

        for tv, v, tp in zip(true_CO2_aquifer_total, CO2_aquifer_total, time_array):
            self.assertTrue(
                abs((tv-v)) < 1.0,
                'Total CO2 leakage rate at time {} is {} but should be {}'
                .format(str(tp), str(v), str(tv)))
        for tv, v, tp in zip(true_brine_aquifer_total, brine_aquifer_total, time_array):
            self.assertTrue(
                abs((tv-v)) < 1.0,
                'Total brine leakage rate at time {} is {} but should be {}'
                .format(str(tp), str(v), str(tv)))

    def test_fault_leakage(self):
        """Tests FaultLeakage component.

        Tests FaultLeakage component in a forward model against expected output
        for 20 data points over 30 years.
        """
        # Define keyword arguments of the system model.
        times = np.linspace(0.0, 30.0, num=20)*365.25 # time in days
        sm_model_kwargs = {'time_array': times} # time must be given in ***days***

        # Create the system model.
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # well_index is the well location index (not to be confused with Peaceman's well index)
        params = {'damage_zone_perm': -13.0, 'damage_zone_por': 0.05,
                  'shallow_aquifer_perm': -13.0, 'deep_aquifer_perm': -12.0,
                  'shallow_aquifer_por': 0.30, 'deep_aquifer_por': 0.20,
                  'well_index': 0.0, 'well_rate': 20.0, 'dip_angle': 40.0,
                  'injection_time': 20.0,
                  'geothermal_gradient': 35.0}

        ffc = sm.add_component_model_object(FaultLeakage(name='ffc', parent=sm))

        # Add parameters
        for key in params:
            ffc.add_par(key, value=params[key], vary=False)

        # Add scalar observations
        ffc.add_obs('brine_aquifer')
        ffc.add_obs('mass_brine_aquifer')
        ffc.add_obs('CO2_aquifer')
        ffc.add_obs('mass_CO2_aquifer')

        sm.forward()

        # Collect observations and convert from kg/sec to tonne/year
        kgs_to_ty = 1.0e-3*365.25*24*3600
        kg_to_tonne = 1.0e-3
        brine_rate_obs = kgs_to_ty * sm.collect_observations_as_time_series(
            ffc, 'brine_aquifer')
        cum_brine_obs = kg_to_tonne * sm.collect_observations_as_time_series(
            ffc, 'mass_brine_aquifer')
        CO2_rate_obs = kgs_to_ty * sm.collect_observations_as_time_series(
            ffc, 'CO2_aquifer')
        cum_CO2_obs = kg_to_tonne * sm.collect_observations_as_time_series(
            ffc, 'mass_CO2_aquifer')

        # # Print the results.
        # print('Brine leakage rate (t/yr)', brine_rate_obs,sep='\n')
        # print('Cumulative brine leakage (t)', cum_brine_obs,sep='\n')
        # print('CO2 leakage rate (t/yr)', CO2_rate_obs,sep='\n')
        # print('Cumulative CO2 leakage (t)', cum_CO2_obs,sep='\n')

        # Perform an accuracy test with a known case.
        true_brine_rate = np.array([0., 175979.9, 145872.62, 129209.84, 125567.67,
               127366.84, 134765.28, 140630.7, 146373.23, 149061.16,
               150104.33, 148167.31, 146717.6, 124379.67, 84806.8,
                68673.25, 58678.89, 48936.89, 40774.6, 33600.992])  # t/yr

        true_CO2_rate = np.array([0., 239869.08, 445098.28, 467009.75, 476176.1,
               478896.06, 479526.47, 482197.6, 489237., 504251.44,
               529424.25, 557741.75, 576606.5, 414577.22, 214732.17,
               110940.52, 69257.9, 51961.707, 46253.23, 42923.3])   # t/yr

        true_cum_brine = np.array([0., 138931.50493421, 393025.60855263,
                610195.97861842, 811336.12253289, 1011021.26644737,
               1217967.68092105, 1435385.57565789, 1661967.63157895,
               1895205.29605263, 2131388.58552632, 2366866.18421053,
               2599670.05756579, 2813694.20230263, 2978841.41447368,
               3100009.87253289, 3200551.03618421, 3285510.86348684,
               3356335.72574013, 3415053.29975329])  # t/yr

        true_cum_CO2 = np.array([0., 189370.32483553, 730134.04194079,
                1450219.30509868, 2194839.73273026, 2948844.04194079,
                3705493.38404605, 4464749.22286184, 5231671.29523026,
                6016004.27220395, 6832064.02549342, 7690352.97286184,
                8585891.06496711, 9368404.55180921, 9865227.74259869,
               10122337.75904606, 10264599.67105263, 10360299.36266448,
               10437837.47121711, 10508239.99588816]) # t/yr

        tol = 1
        self.assertTrue(np.max(abs(brine_rate_obs[1:]-true_brine_rate[1:])/\
                               true_brine_rate[1:]) < tol,
            'Brine leakage rates failed the test for the fault leakage component.')
        self.assertTrue(np.max(abs(CO2_rate_obs[1:]-true_CO2_rate[1:])/true_CO2_rate[1:]) < tol,
            'CO2 leakage rates failed the test for the fault leakage component.')
        self.assertTrue(np.max(abs(cum_brine_obs[1:]-true_cum_brine[1:])/true_cum_brine[1:]) < tol,
            'Total brine leakage failed the test for the fault leakage component.')
        self.assertTrue(np.max(abs(cum_CO2_obs[1:]-true_cum_CO2[1:])/true_cum_CO2[1:]) < tol,
            'Total CO2 leakage failed the test for the fault leakage component.')

        # Test the zeros separately
        self.assertTrue(brine_rate_obs[0] < tol,
            'Brine leakage rates failed the test for the fault leakage component.')
        self.assertTrue(CO2_rate_obs[0] < tol,
            'CO2 leakage rates failed the test for the fault leakage component.')
        self.assertTrue(cum_brine_obs[0] < tol,
            'Total brine leakage failed the test for the fault leakage component.')
        self.assertTrue(cum_CO2_obs[0] < tol,
            'Total CO2 leakage failed the test for the fault leakage component.')

    def test_futuregen_aquifer(self):
        """Tests FutureGen2 aquifer component.

        Tests FutureGen2 aquifer component in a
        forward model against expected output for 2 years of data.
        """
        # Create system model
        time_array = 365.25*np.arange(0, 2)
        sm_model_kwargs = {'time_array': time_array} # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add FutureGen aquifer model object and define parameters
        fga = sm.add_component_model_object(FutureGen2Aquifer(name='fga', parent=sm))

        # St Peter Sandstone
        fga.add_par('aqu_thick', value=61.6)
        fga.add_par('depth', value=591.9)
        fga.add_par('por', value=0.18)
        fga.add_par('log_permh', value=-11.92)
        fga.add_par('log_aniso', value=0.30)
        fga.add_par('rel_vol_frac_calcite', value=0.01)

        fga.model_kwargs['brine_rate'] = 1.0e-3 # kg/s
        fga.model_kwargs['co2_rate'] = 1.0e-2   # kg/s
        fga.model_kwargs['brine_mass'] = 1.0e-3 * 1 *86400*365.25 # kg
        fga.model_kwargs['co2_mass'] = 1.0e-2 * 1 *86400*365.25   # kg

        # Add observations (output) from the aquifer model
        fga.add_obs('TDS_volume')
        fga.add_obs('Pressure_dx')
        fga.add_obs('pH_dz')

        # Run the system model
        sm.forward()

        # Assign observations of the model to variables
        tdsvol = sm.collect_observations_as_time_series(fga, 'TDS_volume')
        pressdx = sm.collect_observations_as_time_series(fga, 'Pressure_dx')
        phdz = sm.collect_observations_as_time_series(fga, 'pH_dz')

        # True values
        true_tds = [0.0, 221.65165841]
        true_press = [0.0, 2.36196175]
        true_pH = [0.0, 40.99599331338357]

        # Test with helpful message!
        for ts, s, t in zip(true_tds, tdsvol, time_array):
            self.assertTrue(abs((ts-s)) < 1.0,
                            'TDS volume at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))

        for ts, s, t in zip(true_press, pressdx, time_array):
            self.assertTrue(abs((ts-s)) < 0.01,
                            'Pressure dx at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))
        for ts, s, t in zip(true_pH, phdz, time_array):
            self.assertTrue(abs((ts-s)) < 0.01,
                            'pH dz at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))

    def test_futuregen_azmi(self):
        """Tests FutureGen2 AZMI component.

        Tests FutureGen2 AZMI component in a
        forward model against expected output for 2 years of data.
        """
        # Create system model
        time_array = 365.25*np.arange(0, 2)
        sm_model_kwargs = {'time_array': time_array} # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add FutureGen aquifer model object and define parameters
        fga = sm.add_component_model_object(FutureGen2AZMI(name='fga', parent=sm))

        # Ironton-Galesville
        fga.add_par('aqu_thick', value=33.2)
        fga.add_par('depth', value=1043.9)
        fga.add_par('por', value=0.118)
        fga.add_par('log_permh', value=-13.39)
        fga.add_par('log_aniso', value=0.30)
        fga.add_par('rel_vol_frac_calcite', value=0.01)

        fga.model_kwargs['brine_rate'] = 1.0e-3 # kg/s
        fga.model_kwargs['co2_rate'] = 1.0e-2 # kg/s
        fga.model_kwargs['brine_mass'] = 1.0e-3 * 1 *86400*365.25 # kg
        fga.model_kwargs['co2_mass'] = 1.0e-2 * 1 *86400*365.25 # kg

        # Add observations (output) from the aquifer model
        fga.add_obs('TDS_volume')
        fga.add_obs('Pressure_dx')
        fga.add_obs('pH_dz')

        # Run the system model
        sm.forward()

        # Assign observations of the model to variables
        tdsvol = sm.collect_observations_as_time_series(fga, 'TDS_volume')
        pressdx = sm.collect_observations_as_time_series(fga, 'Pressure_dx')
        phdz = sm.collect_observations_as_time_series(fga, 'pH_dz')

        # True values
        true_tds = [0.0, 192.45109653]
        true_press = [0.0, 7.61720814]
        true_pH = [0.0, 30.950531971071516]

        # Test with helpful message!
        for ts, s, t in zip(true_tds, tdsvol, time_array):
            self.assertTrue(abs((ts-s)) < 1.0,
                            'TDS volume at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))

        for ts, s, t in zip(true_press, pressdx, time_array):
            self.assertTrue(abs((ts-s)) < 0.01,
                            'Pressure dx at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))
        for ts, s, t in zip(true_pH, phdz, time_array):
            self.assertTrue(abs((ts-s)) < 0.01,
                            'pH dz at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))

    def test_generalized_flow_rate_cmpnt(self):
        """Tests generalized flow rate component.

        Tests a generalized flow rate component in a
        forward simulation against expected output for 5 years of data.
        """
        # Define keyword arguments of the system model
        num_years = 5
        time_array = 365.25*np.arange(0.0, num_years+1)
        sm_model_kwargs = {'time_array': time_array} # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add generalized flow rate component
        gfr = sm.add_component_model_object(GeneralizedFlowRate(name='gfr', parent=sm))

        # Add parameters of generalized flow rate component model
        gfr.add_par('logPeakCO2Rate', value=-4.0, vary=False)
        gfr.add_par('timePeakCO2Rate', value=10.0, vary=False)
        gfr.add_par('durationPeakCO2Rate', value=15.0, vary=False)
        gfr.add_par('durationPeakZeroCO2Rate', value=100.0, vary=False)
        gfr.add_par('logInitBrineRate', value=-10.0, vary=False)
        gfr.add_par('logFinalBrineRate', value=-11.5, vary=False)
        gfr.add_par('durationInitBrineRate', value=2.0, vary=False)
        gfr.add_par('durationInitFinalBrineRate', value=10.0, vary=False)
        gfr.add_par('mitigationTime', value=45.0, vary=False)

        # Add observations of generalized flow rate component model
        gfr.add_obs('CO2_aquifer1')
        gfr.add_obs('brine_aquifer1')
        gfr.add_obs('mass_CO2_aquifer1')
        gfr.add_obs('mass_brine_aquifer1')

        sm.forward()  # system model is run deterministicall

        # To avoid comparing with zero we start with the third element
        CO2_aquifer = sm.collect_observations_as_time_series(
            gfr, 'CO2_aquifer1')[2:]
        brine_aquifer = sm.collect_observations_as_time_series(
            gfr, 'brine_aquifer1')[2:]
        mass_CO2_aquifer = sm.collect_observations_as_time_series(
            gfr, 'mass_CO2_aquifer1')[2:]
        mass_brine_aquifer = sm.collect_observations_as_time_series(
            gfr, 'mass_brine_aquifer1')[2:]

        true_CO2_aquifer = [0.e+00, 1.e-05, 2.e-05, 3.e-05, 4.e-05, 5.e-05]
        true_brine_aquifer = [
            1.0e-10, 1.0e-10, 1.0e-10, 9.03162278e-11, 8.06324555e-11, 7.09486833e-11]
        true_mass_CO2_aquifer = [0., 0., 315.576, 946.728, 1893.456, 3155.76]
        true_mass_brine_aquifer = [
            0., 0.00315576, 0.00631152, 0.00946728, 0.01231744, 0.01486201]

        # Test the produced values against expected
        for tr, r, t in zip(true_CO2_aquifer[2:], CO2_aquifer, time_array[2:]):
            self.assertTrue(
                abs((tr-r)/tr) < 0.01,
                'At time t {} CO2 flow rate is {} kg/s but should be {} kg/s.'
                .format(str(t/365.25), str(r), str(tr)))

        for tr, r, t in zip(true_brine_aquifer[2:], brine_aquifer, time_array[2:]):
            self.assertTrue(
                abs((tr-r)/tr) < 0.01,
                'At time t {} brine flow rate is {} kg/s but should be {} kg/s.'
                .format(str(t/365.25), str(r), str(tr)))

        for tm, m, t in zip(true_mass_CO2_aquifer[2:], mass_CO2_aquifer, time_array[2:]):
            self.assertTrue(
                abs((tm-m)/tm) < 0.01,
                'At time t {} CO2 mass is {} kg but should be {} kg.'
                .format(str(t/365.25), str(m), str(tm)))

        for tm, m, t in zip(true_mass_brine_aquifer[2:], mass_brine_aquifer, time_array[2:]):
            self.assertTrue(
                abs((tm-m)/tm) < 0.01,
                'At time t {} brine mass is {} kg but should be {} kg.'
                .format(str(t/365.25), str(m), str(tm)))

    @unittest.skipIf(GenericAquifer.model_data_check() == False,
                     "Generic Aquifer component model files are not present.")
    def test_generic_aquifer(self):
        """Tests generic aquifer component.

        Tests generic aquifer component in a
        forward model against expected output for 2 years of data.
        """
        # Create system model
        time = 1
        time_array = 365.25*np.arange(0, time+1)
        sm_model_kwargs = {'time_array': time_array} # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add Generic aquifer model object and define parameters
        ga = sm.add_component_model_object(GenericAquifer(name='ga', parent=sm))

        # parameters from run 3263 of 6250 (2100-2499m)
        par = [1.0260E+02, -2.3066E+03, 7.9001E-02, -1.2434E+01, 2.2254E+00,
               5.7080E-03, 4.0938E-01, -5.0149E+00, 4.0617E-02]

        ga.add_par('aqu_thick', value=par[0])
        ga.add_par('top_depth', value=-par[1])
        ga.add_par('por', value=par[2])
        ga.add_par('log_permh', value=par[3])
        ga.add_par('log_aniso', value=par[4])
        ga.add_par('aquifer_salinity', value=par[5])
        ga.add_par('reservoir_salinity', value=par[8])
        ga.add_par('dissolved_salt_threshold', value=0.02)
        ga.add_par('dissolved_co2_threshold', value=0.01)

        log_co2_rate = par[6] # kg/s
        log_brine_rate = par[7] # kg/s
        ga.model_kwargs['brine_mass'] = 10**log_brine_rate * time  * 365.25 * 86400   # kg
        ga.model_kwargs['co2_mass'] = 10**log_co2_rate * time * 365.25 * 86400   # kg

        # Add observations (output) from the aquifer model
        ga.add_obs('Dissolved_salt_volume')
        ga.add_obs('Dissolved_salt_dr')
        ga.add_obs('Dissolved_salt_dz')
        ga.add_obs('Dissolved_CO2_volume')
        ga.add_obs('Dissolved_CO2_dr')
        ga.add_obs('Dissolved_CO2_dz')

        # Run the system model
        sm.forward()

        # Assign observations of the model to variables
        dissolved_salt_volume = sm.collect_observations_as_time_series(
            ga, 'Dissolved_salt_volume')
        dissolved_co2_volume = sm.collect_observations_as_time_series(
            ga, 'Dissolved_CO2_volume')

        # True values
        true_dissolved_salt_volume = [0., 2706.93]
        true_dissolved_co2_volume = [0., 5373325.513648777]

        # Test with helpful message!
        for ts, s, t in zip(true_dissolved_salt_volume, dissolved_salt_volume, time_array):
            self.assertTrue(abs((ts-s)) < 1.0,
                            'Dissolved salt volume at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))

        for ts, s, t in zip(true_dissolved_co2_volume, dissolved_co2_volume, time_array):
            self.assertTrue(abs((ts-s)) < 1.0,
                            'Dissolved CO2 volume at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))

    @unittest.skipIf(GenericReservoir.model_data_check() == False,
                     "Generic Reservoir component model files are not present.")
    def test_generic_reservoir_forward(self):
        """Tests generic reservoir component.

        Tests the system model with a generic reservoir component in a
        forward model against expected output for 75 years of data.
        """
        # Create system model
        distance = 2000 # meter

        num_years = 75
        time_array = 365.0*np.arange(0, num_years+1, 5)

        sm_model_kwargs = {'time_array': time_array} # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add reservoir component
        res = sm.add_component_model_object(GenericReservoir(
            name='res', parent=sm, injX=0., injY=0., locX=distance, locY=0.))

        # Add parameters of reservoir component model
        res.add_par('reservoirDepth', value=2000.0, vary=False)
        res.add_par('reservoirThickness', value=50.0, vary=False)
        res.add_par('logResPerm', value=-14.0, vary=False)
        res.add_par('resTempGradient', value=30.0, vary=False)
        res.add_par('injRate', value=100, vary=False)
        res.add_par('initialSalinity', value=0.05, vary=False)
        res.add_par('wellRadius', value=0.05, vary=False)
        res.add_par('reservoirPorosity', value=0.1, vary=False)

        # Add observations
        res.add_obs('pressure')
        res.add_obs('CO2saturation')

        # Run system model using current values of its parameters
        sm.forward()

        # Assign observations of the model to pressure and CO2saturation variables
        # Obtain pressure and CO2 saturation as lists
        pressure = sm.collect_observations_as_time_series(res, 'pressure')

        # saturation(defined time) at defined location
        saturation = sm.collect_observations_as_time_series(res, 'CO2saturation')

        # True values
        true_pressure = [
            1.91833e+07, 4.92151e+07, 5.48526e+07, 5.77966e+07, 5.96888e+07,
            6.10506e+07, 4.28154e+07, 3.70276e+07, 3.37106e+07, 3.16047e+07,
            3.00346e+07, 2.88335e+07, 2.78929e+07, 2.71535e+07, 2.653e+07,
            2.60014e+07]

        true_saturation = [
           0.00000e+00, 8.20977e-01, 9.80491e-01, 9.92484e-01, 9.94269e-01,
           9.95659e-01, 9.97513e-01, 9.97742e-01, 9.97767e-01, 9.97757e-01,
           9.97733e-01, 9.97702e-01, 9.97664e-01, 9.97621e-01, 9.97576e-01,
           9.97530e-01]

        # Test with helpful message!
        for tp, p, t in zip(true_pressure, pressure, time_array):
            self.assertTrue(abs((tp-p)/tp) < 0.01,
                            'Pressure at time t={} days is {} Pa but should be {} Pa'
                            .format(str(t), str(p), str(tp)))
        for ts, s, t in zip(true_saturation[1:], saturation[1:], time_array[1:]):
            self.assertTrue(abs((ts-s)/ts) < 0.01,
                            'Saturation at time t={} days is {} but should be {}'
                            .format(str(t), str(s), str(ts)))

    def test_hydrocarbon_leakage_forward(self):
        """
        Tests a hydrocarbon leakage component.

        Tests a hydrocarbon leakage component in a forward simulation
        against expected output for 6 time points of data.
        """
        # Define keyword arguments of the system model
        num_years = 50
        t0 = 200.0  # initial time point
        time_array = 365.25*np.arange(t0, t0+num_years+1, 10)

        sm_model_kwargs = {'time_array': time_array}  # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add hydrocarbon leakage component
        hcl_comp = sm.add_component_model_object(HydrocarbonLeakage(
            name='hcl', parent=sm))

        # Add parameters of hydrocarbon component model
        hcl_comp.add_par('reservoirDepth', value=1046.0)
        hcl_comp.add_par('NTG', value=0.444)
        hcl_comp.add_par('logResPerm', value=-13.23)
        hcl_comp.add_par('reservoirPressureMult', value=1.13)
        hcl_comp.add_par('logWellPerm', value=-12.44)
        hcl_comp.add_par('avgWaterSaturation', value=0.719)
        hcl_comp.add_par('FCO2', value=0.647)
        hcl_comp.add_par('FC1', value=0.062)
        hcl_comp.add_par('FC4', value=0.082)
        hcl_comp.add_par('FC7Plus', value=0.209)

        # Add observations (output) of the hydrocarbon leakage component
        hcl_observations = ['mass_CO2_aquifer', 'mass_CO2_gas_aquifer',
                            'mass_methane_oil_aquifer', 'mass_methane_gas_aquifer',
                            'mass_oil_aquifer', 'mass_gas_aquifer']
        for obs_nm in hcl_observations:
            hcl_comp.add_obs(obs_nm)

        # Run forward simulation
        sm.forward()

        # Collect ouputs into a dictionary
        outputs = {obs_nm: sm.collect_observations_as_time_series(hcl_comp, obs_nm)
                   for obs_nm in hcl_observations}

        # True values
        true_vals = {
            'mass_CO2_aquifer': [3298355.25, 4538579., 5796666.,
                                 7057549., 8318430.5, 9579311.],
            'mass_CO2_gas_aquifer': [47086264., 55086008., 62576168.,
                                     69997032., 77417912., 82863688.],
            'mass_methane_oil_aquifer': [287800.0625, 484279.71875, 680759.4375,
                                         876004.8125, 966723.1875, 1057442.125],
            'mass_methane_gas_aquifer': [ 8337887.5, 9125749., 9945482.,
                                         10863910., 12104090., 13668065.],
            'mass_oil_aquifer': [6.33236960e+07, 7.29319360e+07, 8.25402000e+07,
                                 9.21484480e+07, 1.01756696e+08, 1.11364960e+08],
            'mass_gas_aquifer': [6.59736560e+07, 7.41845440e+07, 8.23954240e+07,
                                 9.03425440e+07, 9.76368080e+07, 1.04827896e+08]}

        for key, vals in true_vals.items():
            for ind in range(6):
                true_val = vals[ind]
                sim_val = outputs[key][ind]
                self.assertTrue(
                    abs(true_val-sim_val)/true_val < 0.01,
                    'Observation {} at {} days is {} but should be {}.'.format(
                        key, time_array[ind], sim_val, true_val))

    def test_kimberlina_wellbore(self):
        """
        Tests a Kimberlina wellbore component.

        Tests a Kimberlina wellbore component in a forward simulation
        against expected output for 4 years of data.
        """
        # Create system model
        time_array = 365.25*np.arange(0.0, 4.0)
        sm_model_kwargs = {'time_array': time_array} # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add Kimberlina wellbore model object and define parameters
        kwb = sm.add_component_model_object(KimberlinaWellbore(name='kwb', parent=sm))

        kwb.add_par('logK_well_permeability', value=-10.58833)

        kwb.model_kwargs['pressure'] = 28228195.8 # Pa
        kwb.model_kwargs['CO2saturation'] = 0.45741
        kwb.model_kwargs['previous_co2_mass'] = 348953.5334519825 # kg
        kwb.model_kwargs['previous_brine_mass'] = 2754.645894381841 # kg

        kwb.add_obs('CO2_mass')
        kwb.add_obs('brine_mass')
        kwb.add_obs('CO2_rate')
        kwb.add_obs('brine_rate')

        # Run the system model
        sm.forward()

        # True values
        true_vals = {}
        true_vals['CO2_mass'] = [0.0, 6.13886995e-01, 1.36926020e+01, 1.52343591e+02]
        true_vals['brine_mass'] = [0.0, 39.96582597, 52.93537418, 68.05336453]
        true_vals['CO2_rate'] = [0.0, 1.91360241e-08, 4.14439469e-07, 4.39358472e-06]
        true_vals['brine_rate'] = [0.0, 1.26612372e-06, 4.10980183e-07, 4.79060206e-07]

        indices = [1, 2, 3]
        for key, vals in true_vals.items():
            for ind in indices:
                true_val = vals[ind]
                sim_val = kwb.obs[key+'_{}'.format(ind)].sim
                self.assertTrue(
                    abs(true_val-sim_val)/true_val < 0.01,
                    'Observation {} at {} days is {} but should be {}.'.format(
                        key, time_array[ind], sim_val, true_val))

    def test_lhs(self):
        """Tests Latin Hypercube sampling method of system model class.

        Tests the system model latin hypercube sampling with a simple reservoir
        component coupled to a cemented wellbore component. Tests
        the simulated mean and standard deviation against expected values.

        .. image:: test_output/sens_test_1.png
            :align: center
            :width: 400

        .. image:: test_output/sens_test_2.png
            :align: center
            :width: 400

        """
        # Create system model
        sm = SystemModel()

        # Add reservoir component
        sres = sm.add_component_model_object(SimpleReservoir(name='sres',
                                                             parent=sm,
                                                             locX=550.,
                                                             locY=100.))

        # Add parameters of reservoir component model
        sres.add_par('numberOfShaleLayers', value=3, vary=False)
        sres.add_par('shale1Thickness', value=525.0, vary=False)
        sres.add_par('shale2Thickness', value=475.0, vary=False)
        # Shale 3 has a fixed thickness of 11.2 m
        sres.add_par('shale3Thickness', value=11.2, vary=False)
        # Aquifer 1 (thief zone has a fixed thickness of 22.4)
        sres.add_par('aquifer1Thickness', value=22.4, vary=False)
        # Aquifer 2 (shallow aquifer) has a fixed thickness of 19.2
        sres.add_par('aquifer2Thickness', value=19.2, vary=False)
        # Reservoir has a fixed thickness of 51.2
        sres.add_par('reservoirThickness', value=51.2, vary=False)

        # Add observations of reservoir component model
        sres.add_obs_to_be_linked('pressure')
        sres.add_obs_to_be_linked('CO2saturation')

        # Add cemented wellbore component
        cw = sm.add_component_model_object(CementedWellbore(name='cw', parent=sm))

        # Add parameters of cemented wellbore component
        cw.add_par('logWellPerm', min=-13.9, max=-11., value=-12.8)
        cw.add_par('logThiefPerm', min=-13.9, max=-12.1, value=-12.2)

        # Add keyword arguments of the cemented wellbore component model
        cw.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
        cw.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

        cw.add_composite_par('wellDepth', expr='sres.shale1Thickness' +
                             '+sres.shale2Thickness + sres.shale3Thickness' +
                             '+sres.aquifer1Thickness+ sres.aquifer2Thickness')
        cw.add_composite_par('depthRatio',
                             expr='(sres.shale2Thickness+sres.shale3Thickness' +
                             '+ sres.aquifer2Thickness + sres.aquifer1Thickness/2)/cw.wellDepth')
        cw.add_composite_par('initPressure',
                             expr='sres.datumPressure + cw.wellDepth*cw.g*sres.brineDensity')

        cw.add_obs('CO2_aquifer1')
        cw.add_obs('CO2_aquifer2')
        cw.add_obs('CO2_atm')
        cw.add_obs('brine_aquifer1')
        cw.add_obs('brine_aquifer2')

        # Decided to fix seed: just to make sampling more determined
        s = sm.lhs(siz=1000, seed=465)
        s.run(verbose=False)

        # True values of means and standard deviations obtained by averaging
        # over 100000 simulations
        true_means = [2.849026600358214e-05, 1.4710877772422628e-05,
                      9.722353566396697e-06, 0.00018406411834396143,
                      2.6032908047605337e-06]
        true_stds = [6.523818950580592e-05, 1.6530682199456962e-05,
                     1.6664940620365873e-05, 0.00029734016947531854,
                     8.116163380809236e-07]

        # Add tolerance appropriate for each observation
        # Tolerance calculations are based on probability theory:
        # true_st_dev/sampleSize^0.5/true_means
        tols = [0.03651502, 0.38714553, 0.29447871, 0.07299162, 0.00643125]

        # Test with helpful message!
        labels = s.recarray.dtype.names
        sim_means = np.mean(s.recarray[list(labels[2:])].tolist(), axis=0)
        sim_stds = np.std(s.recarray[list(labels[2:])].tolist(), axis=0)
        for l, tm, sim, tol in zip(labels[2:], true_means, sim_means, tols):
            self.assertTrue(
                abs((tm-sim)/tm) < tol,
                'Mean for {} is {} but should be {}'.format(l, str(sim), str(tm)))
        for l, ts, ss in zip(labels[2:], true_stds, sim_stds):
            self.assertTrue(
                abs((ts-ss)/ts) < 0.05,
                'Standard deviation for {} is {} but should be {}'.format(
                    l, str(ss), str(ts)))

        if not os.path.exists('test_output'):
            os.mkdir('test_output')
        base_outfile = os.path.join('test_output', 'sens_test_{n}.{fmt}')
        # Run RBD_fast sensitivity analysis on CO2 leaked into aquifer
        CO2_aq1 = s.rbd_fast(obsname='cw.CO2_aquifer1', print_to_console=False)
        CO2_aq1_ssens = iam_vis.simple_sensitivities_barplot(
            CO2_aq1, sm,
            title='RBD-Fast Sensitivity\nfor Leakage Rates of $CO_2$ into the Aquifer 1',
            savefig=base_outfile.format(n=1, fmt='png'),
            outfile=base_outfile.format(n=1, fmt='csv'))
        true_CO2aq1_logWell = 0.930413351316639
        self.assertTrue(os.path.exists(base_outfile.format(n=1, fmt='png')),
                        'Output image not created')
        self.assertTrue(os.path.exists(base_outfile.format(n=1, fmt='csv')),
                        'Output data file not created')
        self.assertTrue(abs(true_CO2aq1_logWell - CO2_aq1_ssens[0]) < 0.01,
                        'CO2_aq1-logWellPerm is {ov} but should be {tv}'
                        .format(ov=CO2_aq1_ssens[0], tv=true_CO2aq1_logWell))

        ms = iam_vis.multi_sensitivities_barplot(
            ['cw.CO2_aquifer1', 'cw.brine_aquifer1'], sm, s,
            title='Sensitivity Coefficients',
            savefig=base_outfile.format(n=2, fmt='png'),
            outfile=base_outfile.format(n=2, fmt='csv'))
        true_brine_aq1_logWell = 0.99

        # Copy files to documentation folder
        copyfile(base_outfile.format(n=1, fmt='png'),
                 os.path.join(os.path.join("..", "documentation", "qaqc", "source"),
                              base_outfile.format(n=1, fmt='png')))
        copyfile(base_outfile.format(n=2, fmt='png'),
                 os.path.join(os.path.join("..", "documentation", "qaqc", "source"),
                              base_outfile.format(n=2, fmt='png')))

        self.assertTrue(os.path.exists(base_outfile.format(n=2, fmt='png')),
                        'Output image not created')
        self.assertTrue(os.path.exists(base_outfile.format(n=2, fmt='csv')),
                        'Output data file not created')
        ov = ms['cw.CO2_aquifer1'][0]
        self.assertTrue(abs(true_CO2aq1_logWell - ov) < 0.01,
                        'CO2_aq1-logWellPerm is {ov} but should be {tv}'
                        .format(ov=ov, tv=true_CO2aq1_logWell))
        ov = ms['cw.brine_aquifer1'][0]
        self.assertTrue(abs(true_brine_aq1_logWell - ov) < 0.01,
                        'brine_aq1-logWellPerm is {ov} but should be {tv}'
                        .format(ov=ov, tv=true_brine_aq1_logWell))

    def test_lookup_table_reservoir(self):
        """
        Tests lookup table reservoir component.

        Tests the system model with a lookup table reservoir component
        linked to two interpolators in a forward model against the expected
        output for 6 years of data.
        """

        # Define keyword arguments of the system model
        time_array = 365.25*np.linspace(0.0, 50.0, 6)
        sm_model_kwargs = {'time_array': time_array}  # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Define directory containing lookup tables
        file_dir = os.path.join('..', 'source', 'components', 'reservoir',
                                'lookuptables', 'Test_2d')

        # Read file with signatures of interpolators and names of files with the corresponding data
        signature_data = np.genfromtxt(
            os.path.join(file_dir, 'parameters_and_filenames.csv'),
            delimiter=",", dtype='str')

        # The first row (except the last element) of the file contains names of the parameters
        par_names = signature_data[0, 1:-1]
        num_pars = len(par_names)

        # Create and add two interpolators to the system model
        for ind in range(2):
            signature = {
                par_names[j]: float(signature_data[ind+1, j+1]) for j in range(num_pars)}

            sm.add_interpolator(
                ReservoirDataInterpolator(
                    name='int'+str(ind+1), parent=sm,
                    header_file_dir=file_dir,
                    time_file='time_points.csv',
                    data_file='Reservoir_data_sim0{}.csv'.format(ind+1),
                    index=int(signature_data[ind+1, 0]),
                    signature=signature),
                intr_family='reservoir')

        # Add reservoir component
        ltres = sm.add_component_model_object(
            LookupTableReservoir(
                name='ltres', parent=sm, intr_family='reservoir',
                locX=37478.0, locY=48233.0,
                parameter_filename='parameters_and_filenames.csv'))

        # Add parameters of reservoir component model
        for j in range(num_pars):
            # add first line with values from signature_file
            ltres.add_par(par_names[j], value=float(signature_data[1, j+1]), vary=False)

        # Add observations of reservoir component model
        ltres.add_obs('pressure')
        ltres.add_obs('CO2saturation')

        # Run system model using current values of its parameters
        sm.forward()

        # Collect pressure and saturation observations
        pressure = sm.collect_observations_as_time_series(ltres, 'pressure')[:] / 1.0e+6
        saturation = sm.collect_observations_as_time_series(ltres, 'CO2saturation')[:]

        true_pressure = [27.11992857, 33.16685714, 34.10442857, 34.16092857,
                         34.08864286, 33.98235714]
        true_saturation = [0., 0.34310857, 0.45430571, 0.47404286,
                           0.48760786, 0.50069071]

        # Test with helpful message!
        for tp, p, t in zip(true_pressure, pressure, time_array):
            self.assertTrue(abs((tp-p)/tp) < 0.01,
                            'Pressure at time {} is {} MPa but should be {} MPa'
                            .format(str(t), str(p), str(tp)))
        for ts, s, t in zip(true_saturation[1:], saturation[1:], time_array[1:]):
            self.assertTrue(abs((ts-s)/ts) < 0.01,
                            'Saturation at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))


    def test_lookup_table_reservoir_h5(self):
        """
        Tests lookup table reservoir component with .h5 files instead of
        .csv files.

        Tests the system model with a lookup table reservoir component
        linked to two interpolators in a forward model against the expected
        output for 6 years of data.
        """

        # Define keyword arguments of the system model
        time_array = 365.25*np.linspace(0.0, 50.0, 6)
        sm_model_kwargs = {'time_array': time_array}  # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Define directory containing lookup tables
        file_dir = os.path.join('..', 'source', 'components', 'reservoir',
                                'lookuptables', 'Test_2d')

        # Read file with signatures of interpolators and names of files with the corresponding data
        signature_data = np.genfromtxt(
            os.path.join(file_dir, 'parameters_and_filenames.csv'),
            delimiter=",", dtype='str')

        # The first row (except the last element) of the file contains names of the parameters
        par_names = signature_data[0, 1:-1]
        num_pars = len(par_names)

        # Create and add two interpolators to the system model
        for ind in range(2):
            signature = {
                par_names[j]: float(signature_data[ind+1, j+1]) for j in range(num_pars)}

            sm.add_interpolator(
                ReservoirDataInterpolator(
                    name='int'+str(ind+1), parent=sm,
                    header_file_dir=file_dir,
                    time_file='time_points.csv',
                    data_file='Reservoir_data_sim0{}.h5'.format(ind+1),
                    index=int(signature_data[ind+1, 0]),
                    signature=signature),
                intr_family='reservoir')

        # Add reservoir component
        ltres = sm.add_component_model_object(
            LookupTableReservoir(
                name='ltres', parent=sm, intr_family='reservoir',
                locX=37478.0, locY=48233.0, parameter_filename=''))

        # Add parameters of reservoir component model
        for j in range(num_pars):
            # add first line with values from signature_file
            ltres.add_par(par_names[j], value=float(signature_data[1, j+1]), vary=False)

        # Add observations of reservoir component model
        ltres.add_obs('pressure')
        ltres.add_obs('CO2saturation')

        # Run system model using current values of its parameters
        sm.forward()

        # Collect pressure and saturation observations
        pressure = sm.collect_observations_as_time_series(ltres, 'pressure')[:] / 1.0e+6
        saturation = sm.collect_observations_as_time_series(ltres, 'CO2saturation')[:]

        true_pressure = [27.11992857, 33.16685714, 34.10442857, 34.16092857,
                         34.08864286, 33.98235714]
        true_saturation = [0., 0.34310857, 0.45430571, 0.47404286,
                           0.48760786, 0.50069071]

        # Test with helpful message!
        for tp, p, t in zip(true_pressure, pressure, time_array):
            self.assertTrue(abs((tp-p)/tp) < 0.01,
                            'Pressure at time {} is {} MPa but should be {} MPa'
                            .format(str(t), str(p), str(tp)))
        for ts, s, t in zip(true_saturation[1:], saturation[1:], time_array[1:]):
            self.assertTrue(abs((ts-s)/ts) < 0.01,
                            'Saturation at time {} is {} but should be {}'
                            .format(str(t), str(s), str(ts)))


    def test_lookup_table_reservoir_2d(self):
        """
        Tests lookup table reservoir component 2d interpolation.

        Tests the system model with several lookup table reservoir components
        linked to one interpolator in a forward model against the expected
        output for all years of data provided in lookup table.
        """
        # Define directory containing lookup tables
        file_dir = os.path.join('..', 'source', 'components', 'reservoir',
                                'lookuptables', 'Test_2d')

        # Read file with time points
        time_data = np.genfromtxt(
            os.path.join(file_dir, 'time_points_QAQC.csv'),
            delimiter=",")

        # Define keyword arguments of the system model
        time_array = 365.25*time_data
        sm_model_kwargs = {'time_array': time_array}  # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Create and add interpolator to the system model
        sm.add_interpolator(
            ReservoirDataInterpolator(
                name='int1', parent=sm,
                header_file_dir=file_dir,
                time_file='time_points_QAQC.csv',
                data_file='Reservoir_data_QAQC.csv',
                index=1,
                signature={}),
            intr_family='reservoir')

        # Read file with locations and data
        res_data = np.genfromtxt(
            os.path.join(file_dir, 'Reservoir_data_QAQC.csv'),
            delimiter=",", skip_header=1)

        # Generate randomly indices of points at which observations
        # from reservoir component will be compared to the values provided
        # in the tables
        # np.random.seed(5) # one can use a seed
        num_points = 5
        loc_indices = np.random.randint(0, res_data.shape[0], size=num_points)

        # Initialize list of reservoir components
        ltress = []

        # Initialize list of true observations
        true_pressure = []
        true_saturation = []

        # Add reservoir component for each selected location
        for comp_ind in range(num_points):
            loc_ind = loc_indices[comp_ind]
            ltress.append(sm.add_component_model_object(
                LookupTableReservoir(
                    name='ltres'+str(comp_ind+1), parent=sm,
                    intr_family='reservoir',
                    locX=res_data[loc_ind, 0],
                    locY=res_data[loc_ind, 1],
                    parameter_filename='')))

            ltress[-1].add_par('index', value=1, vary=False)

            # Add observations of reservoir component model
            ltress[-1].add_obs('pressure')
            ltress[-1].add_obs('CO2saturation')

            # Extract pressure/saturation data from lookup table
            true_pressure.append(res_data[loc_ind, 2:2+len(time_array)])
            true_saturation.append(res_data[loc_ind, 2+len(time_array):])

        # Run system model
        sm.forward()

        # Initialize list of computed observations
        pressure = []
        saturation = []
        for comp_ind in range(num_points):
            pressure.append(sm.collect_observations_as_time_series(
                ltress[comp_ind], 'pressure'))
            saturation.append(sm.collect_observations_as_time_series(
                ltress[comp_ind], 'CO2saturation'))

        for comp_ind in range(num_points):
            for tp, p, t in zip(true_pressure[comp_ind], pressure[comp_ind], time_array):
                self.assertTrue(abs((tp-p)/tp) < 0.001,
                                'Pressure at time t={} days is {} Pa but should be {} Pa'
                                .format(str(t), str(p), str(tp)))

            for ts, s, t in zip(true_saturation[comp_ind], saturation[comp_ind], time_array):
                # Small value at the bottom to avoid division by zero for zero saturation
                self.assertTrue(abs((ts-s)/(ts+0.00001)) < 0.001,
                                'CO2 saturation at t={} days is {} but should be {}'
                                .format(t, s, ts))


    def test_lookup_table_reservoir_2d_h5(self):
        """
        Tests lookup table reservoir component 2d interpolation with an .h5 file
        instead of a .csv file.

        Tests the system model with several lookup table reservoir components
        linked to one interpolator in a forward model against the expected
        output for all years of data provided in lookup table.
        """
        # Define directory containing lookup tables
        file_dir = os.path.join('..', 'source', 'components', 'reservoir',
                                'lookuptables', 'Test_2d')

        # Read file with time points
        time_data = np.genfromtxt(
            os.path.join(file_dir, 'time_points_QAQC.csv'),
            delimiter=",")

        # Define keyword arguments of the system model
        time_array = 365.25*time_data
        sm_model_kwargs = {'time_array': time_array}  # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Create and add interpolator to the system model
        sm.add_interpolator(
            ReservoirDataInterpolator(
                name='int1', parent=sm,
                header_file_dir=file_dir,
                time_file='time_points_QAQC.csv',
                data_file='Reservoir_data_QAQC.h5',
                index=1,
                signature={}),
            intr_family='reservoir')

        # Read file with locations and data
        f = h5py.File(os.path.join(file_dir, 'Reservoir_data_QAQC.h5'),'r')

        x = f['x'][()]
        y = f['y'][()]

        # Generate randomly indices of points at which observations
        # from reservoir component will be compared to the values provided
        # in the tables
        # np.random.seed(5) # one can use a seed
        num_points = 5
        loc_indices = np.random.randint(0, len(x), size=num_points)

        # Initialize list of reservoir components
        ltress = []

        # Initialize list of true observations
        true_pressure = [None] * num_points
        true_saturation = [None] * num_points

        # Add reservoir component for each selected location
        for comp_ind in range(num_points):
            loc_ind = loc_indices[comp_ind]
            ltress.append(sm.add_component_model_object(
                LookupTableReservoir(
                    name='ltres'+str(comp_ind+1), parent=sm,
                    intr_family='reservoir',
                    locX=x[loc_ind], locY=y[loc_ind],
                    parameter_filename='')))

            ltress[-1].add_par('index', value=1, vary=False)

            # Add observations of reservoir component model
            ltress[-1].add_obs('pressure')
            ltress[-1].add_obs('CO2saturation')

            # Extract pressure/saturation data from lookup table
            pressure_compiled = np.zeros(len(time_array))
            sat_compiled = np.zeros(len(time_array))
            for tRef in range(0, len(time_array)):
                # Get the pressure and saturation for the current location over time
                pressure_temp = f['pressure_{:.0f}'.format(tRef + 1)][()]
                pressure_compiled[tRef] = pressure_temp[loc_ind]

                sat_temp = f['CO2saturation_{:.0f}'.format(tRef + 1)][()]
                sat_compiled[tRef] = sat_temp[loc_ind]

            true_pressure[comp_ind] = pressure_compiled
            true_saturation[comp_ind] = sat_compiled

        # Run system model
        sm.forward()

        # Initialize list of computed observations
        pressure = []
        saturation = []
        for comp_ind in range(num_points):
            pressure.append(sm.collect_observations_as_time_series(
                ltress[comp_ind], 'pressure'))
            saturation.append(sm.collect_observations_as_time_series(
                ltress[comp_ind], 'CO2saturation'))

        for comp_ind in range(num_points):
            for tp, p, t in zip(true_pressure[comp_ind], pressure[comp_ind], time_array):
                self.assertTrue(abs((tp-p)/tp) < 0.001,
                                'Pressure at time t={} days is {} Pa but should be {} Pa'
                                .format(str(t), str(p), str(tp)))

            for ts, s, t in zip(true_saturation[comp_ind], saturation[comp_ind], time_array):
                # Small value at the bottom to avoid division by zero for zero saturation
                self.assertTrue(abs((ts-s)/(ts+0.00001)) < 0.001,
                                'CO2 saturation at t={} days is {} but should be {}'
                                .format(t, s, ts))


    def test_lookup_table_reservoir_3d(self):
        """
        Tests lookup table reservoir component 3d interpolation.

        Tests the system model with several lookup table reservoir components
        linked to one interpolator in a forward model against the expected
        output for all years of data provided in lookup table.
        """
        # Define directory containing lookup tables
        file_dir = os.path.join('..', 'source', 'components', 'reservoir',
                                'lookuptables', 'Test_3d')

        # Read file with time points
        time_data = np.genfromtxt(
            os.path.join(file_dir, 'time_points.csv'),
            delimiter=",")

        # Define keyword arguments of the system model
        time_array = 365.25*time_data
        sm_model_kwargs = {'time_array': time_array}  # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Create and add interpolator to the system model
        sm.add_interpolator(
            ReservoirDataInterpolator(
                name='int1', parent=sm,
                header_file_dir=file_dir,
                time_file='time_points.csv',
                data_file='Reservoir_data_sim01.csv',
                index=1,
                signature={},
                interp_2d=False),
            intr_family='reservoir')

        # Read file with locations and data
        res_data = np.genfromtxt(
            os.path.join(file_dir, 'Reservoir_data_sim01.csv'),
            delimiter=",", skip_header=1)

        # Generate randomly indices of points at which observations
        # from reservoir component will be compared to the values provided
        # in the tables
        # np.random.seed(5) # one can use a seed
        num_points = 5
        loc_indices = np.random.randint(0, res_data.shape[0], size=num_points)

        # Initialize list of reservoir components
        ltress = []

        # Initialize list of true observations
        true_pressure = []
        true_saturation = []

        # Add reservoir component for each selected location
        for comp_ind in range(num_points):
            loc_ind = loc_indices[comp_ind]
            ltress.append(sm.add_component_model_object(
                LookupTableReservoir(
                    name='ltres'+str(comp_ind+1), parent=sm,
                    intr_family='reservoir',
                    locX=res_data[loc_ind, 0],
                    locY=res_data[loc_ind, 1],
                    locZ=res_data[loc_ind, 2],
                    parameter_filename='')))

            ltress[-1].add_par('index', value=1, vary=False)

            # Add observations of reservoir component model
            ltress[-1].add_obs('pressure')
            ltress[-1].add_obs('CO2saturation')

            # Extract pressure/saturation data from lookup table
            true_pressure.append(res_data[loc_ind, 3:3+len(time_array)])
            true_saturation.append(res_data[loc_ind, 3+len(time_array):])

        # Run system model
        sm.forward()

        # Initialize list of computed observations
        pressure = []
        saturation = []
        for comp_ind in range(num_points):
            pressure.append(sm.collect_observations_as_time_series(
                ltress[comp_ind], 'pressure'))
            saturation.append(sm.collect_observations_as_time_series(
                ltress[comp_ind], 'CO2saturation'))

        for comp_ind in range(num_points):
            for tp, p, t in zip(true_pressure[comp_ind], pressure[comp_ind], time_array):
                self.assertTrue(abs((tp-p)/tp) < 0.001,
                                'Pressure at time t={} days is {} Pa but should be {} Pa'
                                .format(str(t), str(p), str(tp)))

            for ts, s, t in zip(true_saturation[comp_ind], saturation[comp_ind], time_array):
                # Small value at the bottom to avoid division by zero for zero saturation
                self.assertTrue(abs((ts-s)/(ts+0.00001)) < 0.001,
                                'CO2 saturation at t={} days is {} but should be {}'
                                .format(t, s, ts))


    def test_lookup_table_reservoir_3d_h5(self):
        """
        Tests lookup table reservoir component 3d interpolation with an .h5 file
        instead of a .csv file.

        Tests the system model with several lookup table reservoir components
        linked to one interpolator in a forward model against the expected
        output for all years of data provided in lookup table.
        """
        # Define directory containing lookup tables
        file_dir = os.path.join('..', 'source', 'components', 'reservoir',
                                'lookuptables', 'Test_3d')

        # Read file with time points
        time_data = np.genfromtxt(
            os.path.join(file_dir, 'time_points.csv'),
            delimiter=",")

        # Define keyword arguments of the system model
        time_array = 365.25*time_data
        sm_model_kwargs = {'time_array': time_array}  # time is given in days

        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Create and add interpolator to the system model
        sm.add_interpolator(
            ReservoirDataInterpolator(
                name='int1', parent=sm,
                header_file_dir=file_dir,
                time_file='time_points.csv',
                data_file='Reservoir_data_sim01.h5',
                index=1,
                signature={},
                interp_2d=False),
            intr_family='reservoir')

        # Read file with locations and data
        # Read file with locations and data
        f = h5py.File(os.path.join(file_dir, 'Reservoir_data_sim01.h5'),'r')

        x = f['x'][()]
        y = f['y'][()]
        z = f['z'][()]

        # Generate randomly indices of points at which observations
        # from reservoir component will be compared to the values provided
        # in the tables
        # np.random.seed(5) # one can use a seed
        num_points = 5
        loc_indices = np.random.randint(0, len(x), size=num_points)

        # Initialize list of reservoir components
        ltress = []

        # Initialize list of true observations
        true_pressure = [None] * num_points
        true_saturation = [None] * num_points

        # Add reservoir component for each selected location
        for comp_ind in range(num_points):
            loc_ind = loc_indices[comp_ind]
            ltress.append(sm.add_component_model_object(
                LookupTableReservoir(
                    name='ltres'+str(comp_ind+1), parent=sm,
                    intr_family='reservoir',
                    locX=x[loc_ind], locY=y[loc_ind], locZ=z[loc_ind],
                    parameter_filename='')))

            ltress[-1].add_par('index', value=1, vary=False)

            # Add observations of reservoir component model
            ltress[-1].add_obs('pressure')
            ltress[-1].add_obs('CO2saturation')

            # Extract pressure/saturation data from lookup table

            # Extract pressure/saturation data from lookup table
            pressure_compiled = np.zeros(len(time_array))
            sat_compiled = np.zeros(len(time_array))
            for tRef in range(0, len(time_array)):
                # Get the pressure and saturation for the current location over time
                pressure_temp = f['pressure_{:.0f}'.format(tRef + 1)][()]
                pressure_compiled[tRef] = pressure_temp[loc_ind]

                sat_temp = f['CO2saturation_{:.0f}'.format(tRef + 1)][()]
                sat_compiled[tRef] = sat_temp[loc_ind]

            true_pressure[comp_ind] = pressure_compiled
            true_saturation[comp_ind] = sat_compiled

        # Run system model
        sm.forward()

        # Initialize list of computed observations
        pressure = []
        saturation = []
        for comp_ind in range(num_points):
            pressure.append(sm.collect_observations_as_time_series(
                ltress[comp_ind], 'pressure'))
            saturation.append(sm.collect_observations_as_time_series(
                ltress[comp_ind], 'CO2saturation'))

        for comp_ind in range(num_points):
            for tp, p, t in zip(true_pressure[comp_ind], pressure[comp_ind], time_array):
                self.assertTrue(abs((tp-p)/tp) < 0.001,
                                'Pressure at time t={} days is {} Pa but should be {} Pa'
                                .format(str(t), str(p), str(tp)))

            for ts, s, t in zip(true_saturation[comp_ind], saturation[comp_ind], time_array):
                # Small value at the bottom to avoid division by zero for zero saturation
                self.assertTrue(abs((ts-s)/(ts+0.00001)) < 0.001,
                                'CO2 saturation at t={} days is {} but should be {}'
                                .format(t, s, ts))

    @unittest.skipIf(CarbonateAquifer.model_data_check() == False,
                     "Carbonate Aquifer component library file is not present.")
    @unittest.skipIf(AtmosphericROM.model_data_check() == False,
                     "Atmospheric ROM component library file is not present.")
    def test_openiam_cf(self):
        """Tests Control File functionality of the NRAP-Open-IAM.

        Tests the NRAP-Open-IAM run from a control file with a simple reservoir
        component in a forward model.
        """
        import pickle
        from openiam.openiam_cf import main as openiam_cf

        # Change directories to find input and output files
        current_work_dir = os.getcwd()
        test_dir = os.path.dirname(os.path.abspath(__file__))
        os.chdir(test_dir)

        outfiles = ['SimpleReservoir1_000.pressure.txt',
                    'SimpleReservoir1_000.CO2saturation.txt']

        for ofile in outfiles:
            odfile = os.path.join('test_output', ofile)
            if os.path.exists(odfile):
                os.remove(odfile)

        # Run Test Control File
        test_control_file = 'test_control_file.yaml'
        self.assertTrue(openiam_cf(test_control_file),
                        'NRAP-Open-IAM failed to run properly')

        # Test output directory was created
        self.assertTrue(os.path.exists('test_output'),
                        'Output directory not created')
        # Read yaml_dump file
        with open(os.sep.join(['test_output', 'combined_output.pkl']), 'rb') as yaml_out:
            output = pickle.load(yaml_out)

        # Read output pressure and saturation files
        out_vals = {}
        out_vals['pressure'] = np.loadtxt(os.sep.join(
            ['test_output', 'csv_files', 'SimpleReservoir1_000.pressure.csv']),
            skiprows=1, delimiter=',')[:, 1]
        out_vals['saturation'] = np.loadtxt(os.sep.join(
            ['test_output', 'csv_files', 'SimpleReservoir1_000.CO2saturation.csv']),
            skiprows=1, delimiter=',')[:, 1]
        out_vals['co2aq1'] = np.loadtxt(os.sep.join(
            ['test_output', 'csv_files', 'CementedWellbore1_000.CO2_aquifer1.csv']),
            skiprows=1, delimiter=',')[:, 1]
        out_vals['co2aq2'] = np.loadtxt(os.sep.join(
            ['test_output', 'csv_files', 'CementedWellbore1_000.CO2_aquifer2.csv']),
            skiprows=1, delimiter=',')[:, 1]
        out_vals['co2atm'] = np.loadtxt(os.sep.join(
            ['test_output', 'csv_files', 'CementedWellbore1_000.CO2_atm.csv']),
            skiprows=1, delimiter=',')[:, 1]
        out_vals['brine_aq1'] = np.loadtxt(os.sep.join(
            ['test_output', 'csv_files', 'CementedWellbore1_000.brine_aquifer1.csv']),
            skiprows=1, delimiter=',')[:, 1]
        out_vals['brine_aq2'] = np.loadtxt(os.sep.join(
            ['test_output', 'csv_files', 'CementedWellbore1_000.brine_aquifer2.csv']),
            skiprows=1, delimiter=',')[:, 1]
        out_vals['TDS_volume'] = np.loadtxt(os.sep.join(
            ['test_output', 'csv_files', 'CarbonateAquifer1.TDS_volume.csv']),
            skiprows=1, delimiter=',')[:, 1]
        out_vals['msw_co2'] = np.loadtxt(os.sep.join(
            ['test_output', 'csv_files', 'msw1_000.CO2_aquifer1.csv']),
            skiprows=1, delimiter=',')[:, 1]
        out_vals['msw_brine'] = np.loadtxt(os.sep.join(
            ['test_output', 'csv_files', 'msw1_000.brine_aquifer1.csv']),
            skiprows=1, delimiter=',')[:, 1]
        out_vals['crit_dist'] = np.loadtxt(os.sep.join(
            ['test_output', 'csv_files', 'atmRom1.critical_distance_s000.csv']),
            skiprows=1, delimiter=',')[:, 1]

        # Change back to original working directory
        os.chdir(current_work_dir)

        # Test model parameters set correctly
        modelparams = output['ModelParams']
        self.assertTrue(modelparams['EndTime'] == 5,
                        'EndTime not loaded correctly')
        self.assertTrue(modelparams['Analysis_type'] == 'forward',
                        'Analysis_type not loaded correctly')
        self.assertTrue(modelparams['TimeStep'] == 1.0,
                        'TimeStep not loaded correctly')

        # Test Stratigraphy parameters set correctly
        true_sr_params = {
            'aquifer1Thickness': {'value': 142.4, 'vary': False},
            'aquifer2Thickness': {'value': 19.2, 'vary': False},
            'numberOfShaleLayers': {'value': 3, 'vary': False},
            'reservoirThickness': {'value': 51.2, 'vary': False},
            'shale1Thickness': {'max': 550.0, 'min': 400.0, 'value': 405.0},
            'shale2Thickness': {'max': 500.0, 'min': 450.0, 'value': 475.0},
            'shale3Thickness': {'value': 11.2, 'vary': False},
            'datumPressure': {'value': 101325.0, 'vary': False}}
        stratparams = output['Stratigraphy']
        for key, value in true_sr_params.items():
            self.assertTrue(stratparams[key] == value,
                            '{key} not loaded into Stratigraphy properly'.format(key=key))

        # Test pressure and saturation values
        true_vals = {}
        true_vals['pressure'] = [10418765., 14919301.38843215, 15192401.37235173,
                                 15352154.92651543, 15465501.96157327, 15553420.7877774]
        true_vals['saturation'] = [
            0., 0.21384356, 0.30898839, 0.38199558, 0.44354349, 0.49776829]
        true_vals['co2aq1'] = [
            0.000000000000000000e+00, 4.009259450441665808e-05, 3.328369138093752655e-05,
            2.693611103739417147e-05, 2.107207114594151191e-05, 1.560442836007955292e-05]
        true_vals['co2aq2'] = [
            0.000000000000000000e+00, 1.858034134482971605e-06, 7.995150506733141689e-06,
            1.277841678447248961e-05, 1.684263101312697962e-05, 2.044104130297448487e-05]
        true_vals['co2atm'] = [
            0.000000000000000000e+00, -1.030710967084584508e-08, -2.143717012310900791e-07,
            1.353328298553025458e-06, 2.674960534125898739e-06, 3.268333388857705667e-06]
        true_vals['brine_aq1'] = [
            0.000000000000000000e+00, 1.174587248555037712e-05, 1.181246277069288336e-05,
            1.185186206817051792e-05, 1.187983716205768980e-05, 1.190154380165251380e-05]
        true_vals['brine_aq2'] = [
            0.000000000000000000e+00, 1.029048012603127806e-06, 1.076541008069456788e-06,
            1.104337101947985265e-06, 1.124060054107538356e-06, 1.139358847736474137e-06]
        true_vals['TDS_volume'] = [
            0.00000000000e+00, 4953997.530131221, 4959028.588699669,
            4965434.020670503, 4972596.912364215, 4980559.7430431545]
        true_vals['msw_co2'] = [
            0.000000000000000000e+00, 1.685779384170747642e-06, 2.332231253976781943e-06,
            2.439835615310001397e-06, 2.501344676892450799e-06, 2.555921287545392011e-06]
        true_vals['msw_brine'] = [
            0.000000000000000000e+00, 7.625411632876481455e-06, 7.970275967325878469e-06,
            7.988670590201842547e-06, 7.958216246107291283e-06, 7.914684507949380389e-06]
        true_vals['crit_dist'] = [
            1.290765821949e+02, 2.067820300959e+02, 2.273035260713e+02,
            2.273035260713e+02, 2.350109081766e+02, 2.387713214789e+02]

        for key, vals in true_vals.items():
            for tv, ov in zip(vals, out_vals[key]):
                self.assertTrue(
                    abs(tv - ov) < 0.01,
                    '{key} is {ov} but should be {tv}'.format(
                        key=key, ov=ov, tv=tv))

    def test_parstudy(self):
        """Test parameter study method of system model class.

        Tests the system model parameter study method with a simple reservoir
        component coupled to a cemented wellbore component.
        """
        # Create system model
        sm = SystemModel()

        # Add reservoir component
        sres = sm.add_component_model_object(
            SimpleReservoir(name='sres', parent=sm, locX=550., locY=100.))

        # Add parameters of reservoir component model
        sres.add_par('numberOfShaleLayers', value=3, vary=False)
        sres.add_par('shale1Thickness', value=525.0, vary=False)
        sres.add_par('shale2Thickness', value=475.0, vary=False)
        # Shale 3 has a fixed thickness of 11.2 m
        sres.add_par('shale3Thickness', value=11.2, vary=False)
        # Aquifer 1 (thief zone has a fixed thickness of 22.4)
        sres.add_par('aquifer1Thickness', value=22.4, vary=False)
        # Aquifer 2 (shallow aquifer) has a fixed thickness of 19.2
        sres.add_par('aquifer2Thickness', value=19.2, vary=False)
        # Reservoir has a fixed thickness of 51.2
        sres.add_par('reservoirThickness', value=51.2, vary=False)

        # Add observations of reservoir component model
        sres.add_obs_to_be_linked('pressure')
        sres.add_obs_to_be_linked('CO2saturation')

        # Add cemented wellbore component
        cw = sm.add_component_model_object(CementedWellbore(name='cw', parent=sm))

        # Add parameters of cemented wellbore component
        cw.add_par('logWellPerm', min=-13.9, max=-11., value=-12.8)
        cw.add_par('logThiefPerm', min=-13.9, max=-12.1, value=-12.2)

        # Add keyword arguments of the cemented wellbore component model
        cw.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
        cw.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

        cw.add_composite_par('wellDepth', expr='sres.shale1Thickness' +
                             '+sres.shale2Thickness + sres.shale3Thickness' +
                             '+sres.aquifer1Thickness+ sres.aquifer2Thickness')
        cw.add_composite_par(
            'depthRatio',
            expr='(sres.shale2Thickness+sres.shale3Thickness' +
            '+ sres.aquifer2Thickness + sres.aquifer1Thickness/2)/cw.wellDepth')
        cw.add_composite_par(
            'initPressure',
            expr='sres.datumPressure + cw.wellDepth*cw.g*sres.brineDensity')

        cw.add_obs('CO2_aquifer1')
        cw.add_obs('CO2_aquifer2')
        cw.add_obs('CO2_atm')
        cw.add_obs('brine_aquifer1')
        cw.add_obs('brine_aquifer2')
        s = sm.parstudy(nvals=3)
        s.run(verbose=False)

        # True values
        truth = np.array([
            [-1.39000000e+01, -1.39000000e+01, -0.00010903861900342573,
             -8.779162134229658e-06, -3.6713236636101406e-07,
             6.8195568862251325e-06, 1.0755547742394922e-06],
            [-1.39000000e+01, -1.30000000e+01, -0.00010903861900342573,
             -8.779162134229658e-06, -3.6713236636101406e-07,
             6.8195568862251325e-06, 1.0755547742394922e-06],
            [-1.39000000e+01, -1.21000000e+01, -0.0001090386190034257,
             -8.779162134229658e-06, -3.6713236636101406e-07,
             6.8195568862251325e-06, 1.0755547742394922e-06],
            [-1.24500000e+01, -1.39000000e+01, 3.779689757666963e-05,
             1.7609250171761864e-05, -1.9731563542096083e-07,
             2.884248876918178e-05, 3.503127221225778e-06],
            [-1.24500000e+01, -1.30000000e+01, 3.779689757666963e-05,
             1.7609250171761864e-05, -1.9731563542096083e-07,
             2.884248876918178e-05, 3.503127221225778e-06],
            [-1.24500000e+01, -1.21000000e+01, 3.779689757666963e-05,
             1.7609250171761864e-05, -1.9731563542096083e-07,
             2.884248876918178e-05, 3.503127221225778e-06],
            [-1.10000000e+01, -1.39000000e+01, 0.00013106710530757463,
             5.909689485163282e-05, 6.061493326793466e-05,
             0.0013189921053484584, 4.008159997299794e-06],
            [-1.10000000e+01, -1.30000000e+01, 0.00013106710530757463,
             5.909689485163282e-05, 6.061493326793466e-05,
             0.0013189921053484584, 4.008159997299794e-06],
            [-1.10000000e+01, -1.21000000e+01, 0.00013106710530757463,
             5.909689485163282e-05, 6.061493326793466e-05,
             0.0013189921053484584, 4.008159997299794e-06]])

        # Test with helpful message!
        labels = s.recarray.dtype.names
       # Test that parameter values are exactly correct
        for i, l in enumerate(labels[:2]):
            for j, t, sim in zip(list(range(truth.shape[0])), truth[:, i], s.recarray[l]):
                self.assertEqual(
                    t, sim,
                    'Parameter {} for sample {} is {} but should be {}'.format(
                        l, j+1, str(sim), str(t)))
        # Test simulated values are within tolerance
        for i, l in enumerate(labels[2:]):
            for j, t, sim in zip(list(range(truth.shape[0])), truth[:, i+2], s.recarray[l]):
                self.assertTrue(abs((t-sim)/(t+0.000001)) < 0.01,
                                '{} for sample {} is {} but should be {}'
                                .format(l, j+1, str(sim), str(t)))

    def test_plume_stability_cmpnt(self):
        """Tests work of plume stability component.

        Tests system model containing only a plume stability component present
        with system model time array equal to the time points of the linked
        dataset.
        """
        file_directory = os.sep.join(['..', 'source', 'components', 'reservoir',
                                      'lookuptables', 'Test_2d'])
        # Time array is the same as in the data set but converted to days
        time_array = np.genfromtxt(
            os.sep.join([file_directory, 'time_points_PSA.csv']), delimiter=',')*365.25

        sm_model_kwargs = {'time_array': time_array}   # time is given in days
        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add plume stability analysis component for pressure
        sps = sm.add_component_model_object(
            PlumeStability(
                name='sps', parent=sm, file_directory=file_directory,
                variable_names=['pressure'], thresholds={'pressure': 3.8e7},
                parameter_filename='parameters_and_filenames_PSA.csv',
                time_file='time_points_PSA.csv'))
        sps.add_par('index', value=1, vary=False)

        obs_names = ['pressure_{}'.format(nm) for nm in ['areas', 'areas_dt',
                                                         'mobility', 'spreading']]
        for nm in obs_names:
            sps.add_obs(nm)

        sm.forward()

        # True values: simulation results for the 6th time point (counting from 0):
        # area, its derivative, mobility, spreading
        true_results = [69090000., 784000., 6.82994397e-01, -3.32444377e+02]

        time_index = 6
        sim_data = [sps.obs['{}_{}'.format(nm, time_index)].sim for nm in obs_names]

        for tr, sd, obs in zip(true_results, sim_data, obs_names):
            self.assertTrue(
                abs((tr-sd)/tr) < 0.01,
                'At time t={} years {} is {} but should be {}'.format(
                    time_array[time_index]/365.25, obs, sd, tr))

    def test_plume_stability_cmpnt_h5(self):
        """Tests work of plume stability component with .h5 files instead of
        .csv files.

        Tests system model containing only a plume stability component present
        with system model time array equal to the time points of the linked
        dataset.
        """
        file_directory = os.sep.join(['..', 'source', 'components', 'reservoir',
                                      'lookuptables', 'Test_2d'])
        # Time array is the same as in the data set but converted to days
        time_array = np.genfromtxt(
            os.sep.join([file_directory, 'time_points_PSA.csv']), delimiter=',')*365.25

        sm_model_kwargs = {'time_array': time_array}   # time is given in days
        # Create system model
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add plume stability analysis component for pressure
        sps = sm.add_component_model_object(
            PlumeStability(
                name='sps', parent=sm, file_directory=file_directory,
                variable_names=['pressure'], thresholds={'pressure': 3.8e7},
                parameter_filename='parameters_and_filenames_PSA_h5.csv',
                time_file='time_points_PSA.csv'))
        sps.add_par('index', value=1, vary=False)

        obs_names = ['pressure_{}'.format(nm) for nm in ['areas', 'areas_dt',
                                                         'mobility', 'spreading']]
        for nm in obs_names:
            sps.add_obs(nm)

        sm.forward()

        # True values: simulation results for the 6th time point (counting from 0):
        # area, its derivative, mobility, spreading
        true_results = [69090000., 784000., 6.82994397e-01, -3.32444377e+02]

        time_index = 6
        sim_data = [sps.obs['{}_{}'.format(nm, time_index)].sim for nm in obs_names]

        for tr, sd, obs in zip(true_results, sim_data, obs_names):
            self.assertTrue(
                abs((tr-sd)/tr) < 0.01,
                'At time t={} years {} is {} but should be {}'.format(
                    time_array[time_index]/365.25, obs, sd, tr))

    def test_rate_to_mass_adapter(self):
        """Tests adapter component.

        Tests system model containing only an adapter component present
        with fixed in time input parameters against expected output.
        """
        # Create system model
        num_years = 5
        time_array = 365.25*np.arange(0.0, num_years+1)
        sm_model_kwargs = {'time_array': time_array}  # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add adapter component
        adapt = sm.add_component_model_object(RateToMassAdapter(name='adapt', parent=sm))
        adapt.model_kwargs['CO2_aquifer1'] = [2.0e-6, 4.0e-7]
        adapt.model_kwargs['CO2_aquifer2'] = [4.0e-7, 1.0e-8]
        adapt.add_obs('mass_CO2_aquifer1')
        adapt.add_obs('mass_CO2_aquifer2')
        sm.forward()
        true_masses = [302.95296, 73.844784]

        # Test with helpful message!
        labels = ['mass_CO2_aquifer1', 'mass_CO2_aquifer2']
        for tf, l in zip(true_masses, labels):
            # The observations at the last time point can be accessed directly,
            # e.g. as adapt.obs['mass_CO2_aquifer1_5'].sim;
            # general form: adapt.obs[l+'_'+str(num_years)].sim.
            self.assertTrue(abs((tf-adapt.obs[l+'_'+str(num_years)].sim)/tf) < 0.01,
                            'Flowrate at {} is {} but should be {}'
                            .format(l, str(adapt.obs[l+'_'+str(num_years)].sim), str(tf)))

    def test_reservoir_data_interpolator(self):
        """Tests reservoir data interpolator.

        Tests reservoir data interpolator for a sample data set
        against expected output for 6 years of data.
        """
        # Create system model
        sm = SystemModel()

        # Create interpolator
        int1 = sm.add_interpolator(
            ReservoirDataInterpolator(
                name='int1', parent=sm,
                header_file_dir=os.path.join(
                    '..', 'source', 'components', 'reservoir',
                    'lookuptables', 'Test_2d'),
                time_file='time_points.csv',
                data_file='Reservoir_data_sim01.csv',
                index=1,
                signature={'logResPerm': -13.2773,
                           'reservoirPorosity': 0.2145,
                           'logShalePerm': -18.7}),
            intr_family='reservoir')

        # Setup location of interest
        locX, locY = [37478.0, 48233.0]

        from openiam.iam_gridded_observation import interp_weights
        # Calculate weights of the location of interest
        vertices, weights = interp_weights(int1.triangulation, np.array([[locX, locY]]))

        # Calculate pressure and saturation for several time points
        time_array = 365.25*np.linspace(0, 50, 6)
        pressure = []
        saturation = []

        for t in time_array:
            output = int1(t, vertices, weights)
            pressure.append(output['pressure'][0]/1.0e+6)
            saturation.append(output['CO2saturation'][0])

        # True values
        true_pressure = [
            27.119928571428566, 33.16685714285714, 34.104428571428571,
            34.160928571428578, 34.088642857142858, 33.98235714285714]
        true_saturation = [
            0.0, 0.34310857142857143, 0.45430571428571431,
            0.4740428571428571, 0.4876078571428572, 0.50069071428571432]

        # Test with helpful message!
        for tp, p, t in zip(true_pressure, pressure, time_array):
            self.assertTrue(abs((tp-p)/tp) < 0.01,
                            'Pressure at time {} is {} MPa but should be {} MPa'
                            .format(str(t/365.25), str(p), str(tp)))
        # To avoid comparing zero saturations we start with the second element
        for ts, s, t in zip(true_saturation[1:], saturation[1:], time_array[1:]):
            self.assertTrue(abs((ts-s)/ts) < 0.01,
                            'Saturation at time {} is {} but should be {}'
                            .format(str(t/365.25), str(s), str(ts)))

    def test_seal_horizon(self):
        """ Tests SealHorizon component.

        Test the system model with a seal horizon component in a forward
        simulation against expected output for 2 years of data.
        """
        import components.seal.seal_units as sunit
        # Time constants.
        num_years = 1
        time_array = 365.25 * np.arange(0.0, num_years+1)
        seal_model_kwargs = {'time_array': time_array} # time is given in days

        # Create system model.
        sm = SystemModel(model_kwargs=seal_model_kwargs)
        shc = sm.add_component_model_object(
            SealHorizon(name='shc', parent=sm,
                        locX=[50.0, 1150.0], locY=[550.0, 50.0], area=[100.0, 100.0]))

        # Parameters to vary. Arrays size of keyword arguments coincide with the size
        # of arrays containing coordinates of the cell coordinates
        shc.model_kwargs['thickness'] = [75.0, 109.38422609226659]
        shc.model_kwargs['permeability'] = 2*[2.5e-2 * sunit.microd_to_metersq()]
        shc.model_kwargs['baseDepth'] = [1130.0, 1130.0]
        shc.add_par('aveBaseDepth', value=1130.0, vary=False)
        shc.add_par('aveBasePressure', value=32000000.0, vary=False)

        shc.add_par('influenceModel', value=0, vary=False)
        shc.add_par('aveTemperature', value=50.0, vary=False)
        shc.add_par('salinity', value=15000.0, vary=False)
        shc.add_par('staticDepth', value=1000.0, vary=False)
        shc.add_par('staticPressure', value=10000000.0, vary=False)

        # Fluid parameters
        shc.add_par('brineDensity', value=1004.0, vary=False)
        shc.add_par('CO2Density', value=597.8, vary=False)
        shc.add_par('brineViscosity', value=5.634e-4, vary=False)
        shc.add_par('CO2Viscosity', value=4.52e-5, vary=False)
        shc.add_par('CO2Solubility', value=0.0, vary=False)

        # Two-phase model parameters for L-E-T model
        shc.add_par('wetting1', value=1.0, vary=False)
        shc.add_par('wetting2', value=10.0, vary=False)
        shc.add_par('wetting3', value=1.25, vary=False)
        shc.add_par('nonwet1', value=1.05, vary=False)
        shc.add_par('nonwet2', value=3.5, vary=False)
        shc.add_par('nonwet3', value=1.25, vary=False)
        shc.add_par('capillary1', value=0.2, vary=False)
        shc.add_par('capillary2', value=2.8, vary=False)
        shc.add_par('capillary3', value=0.43, vary=False)
        shc.add_par('maxCapillary', value=10000000.0, vary=False)

        # Additional parameters for two-phase
        shc.add_par('permRatio', value=0.6, vary=False)
        shc.add_par('entryPressure', value=5.0e+3, vary=False)
        shc.add_par('brineResSaturation', value=0.15)
        shc.add_par('CO2ResSaturation', value=0.0)

        # Add time varying keyword arguments
        shc.add_dynamic_kwarg('CO2saturation', [2*[0.5], 2*[0.5]])
        shc.add_dynamic_kwarg('pressure', [2*[3.20E+07], 2*[3.20E+07]])

        # Add gridded observations of the SealHorizon component
        output_dir = 'test_output'
        shc.add_grid_obs('CO2_aquifer', constr_type='array', output_dir=output_dir)
        shc.add_grid_obs('brine_aquifer', constr_type='array', output_dir=output_dir)

        # Run forward simulation
        sm.forward()

        # Collect gridded observations
        grid_outputs = {}
        for nm in ['CO2_aquifer', 'brine_aquifer']:
            filename = 'shc_{}_sim_0_time_1.npz'.format(nm)
            d = np.load(os.sep.join([output_dir, filename]))
            grid_outputs[nm] = d['data']
            d.close()
            grid_outputs[nm] = np.array(grid_outputs[nm])

        # Outputs at time t = 1 year
        true_outputs = {
            'CO2_aquifer': 597.8*np.array(
                [2.580212167001341e-09, 1.7827398785245965e-09]),
            'brine_aquifer': 1004.0*np.array(
                [8.954687353788802e-11, 6.13983912970832e-11])}

        for nm in ['CO2_aquifer', 'brine_aquifer']:
            for tv, v in zip(true_outputs[nm], grid_outputs[nm]):
                self.assertTrue(
                    abs((tv-v)/tv) < 0.01,
                    'Observation {} at time t=1 year is {} but should be {}.'.format(
                        nm, v, tv))

    def test_sh_thickness_sampler(self):
        """ Tests SH Thickness Sampler.

        Test the system model with a thickness sampler in a forward
        simulation against expected output.
        """
        model_kwargs = {'time_point': 0.0} # time is given in days

        current_work_dir = os.getcwd()
        output_dir = os.path.join(current_work_dir, 'test_output')
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # Create system model.
        sm = SystemModel(model_kwargs=model_kwargs)

        # Define parameter values
        ave = 500.0
        std_dev = 50.0
        t_min = 450
        t_max = 575

        # Add sampler
        t_sampler = sm.add_component_model_object(
            SHThicknessSampler(name='ts', parent=sm, grid_shape=(1000, )))
        # Add sampler parameters
        t_sampler.add_par('seed', value=234, vary=False)
        t_sampler.add_par('aveThickness', value=ave, vary=False)
        t_sampler.add_par('stdDevThickness', value=std_dev, vary=False)
        t_sampler.add_par('minThickness', value=t_min, vary=False)
        t_sampler.add_par('maxThickness', value=t_max, vary=False)
        t_sampler.add_grid_obs('thickness', constr_type='array', index=[0],
                               output_dir=output_dir)

        sm.forward()

        # Collect saved observation
        data = sm.collect_gridded_observations_as_time_series(
            t_sampler, 'thickness', output_dir,
            indices=[0], rlzn_number=0)

        diff = ave-np.mean(data)
        min_val = np.min(data)
        max_val = np.max(data)

        self.assertTrue(
            diff < std_dev,
            ''.join(['The difference of means is {} but should be ',
                     'less than {}']).format(diff, std_dev))

        self.assertTrue(
            min_val > t_min,
            ''.join(['The minimum value of thickness array is {} but ',
                     'should be greater than {}']).format(min_val, t_min))

        self.assertTrue(
            max_val < t_max,
            ''.join(['The maximum value of thickness array is {} but ',
                     'should be less than {}']).format(max_val, t_max))

    def test_sh_permeability_sampler(self):
        """ Tests SH Permeability Sampler.

        Test the system model with a permeability sampler in a forward
        simulation against expected output.
        """
        model_kwargs = {'time_point': 0.0} # time is given in days

        current_work_dir = os.getcwd()
        output_dir = os.path.join(current_work_dir, 'test_output')
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # Create system model.
        sm = SystemModel(model_kwargs=model_kwargs)

        # Define parameter values
        ave = 1.0e-18
        std_dev = 1.0e-19
        t_min = 1.0e-19
        t_max = 1.0e-17

        # Add permeability sampler
        perm_sampler = sm.add_component_model_object(
            SHPermeabilitySampler(name='ps', parent=sm, grid_shape=(1000, )))
        # Add sampler parameters
        perm_sampler.add_par('seed', value=234, vary=False)
        perm_sampler.add_par('avePermeability', value=ave, vary=False)
        perm_sampler.add_par('stdDevPermeability', value=std_dev, vary=False)
        perm_sampler.add_par('minPermeability', value=t_min, vary=False)
        perm_sampler.add_par('maxPermeability', value=t_max, vary=False)
        perm_sampler.add_grid_obs('permeability', constr_type='array', index=[0],
                                  output_dir=output_dir)

        sm.forward()

        # Collect saved observation
        data = sm.collect_gridded_observations_as_time_series(
            perm_sampler, 'permeability', output_dir,
            indices=[0], rlzn_number=0)

        diff = ave-np.mean(data)
        min_val = np.min(data)
        max_val = np.max(data)

        self.assertTrue(
            diff < std_dev,
            ''.join(['The difference of means is {} but should be ',
                     'less than {}']).format(diff, std_dev))

        self.assertTrue(
            min_val > t_min,
            ''.join(['The minimum value of permeability array is {} but ',
                     'should be greater than {}']).format(min_val, t_min))

        self.assertTrue(
            max_val < t_max,
            ''.join(['The maximum value of permeability array is {} but ',
                     'should be less than {}']).format(max_val, t_max))

    def test_simple_reservoir_forward(self):
        """Tests simple reservoir component.

        Tests the system model with a simple reservoir component in a
        forward model against expected output for 4 years of data.
        """
        # Create system model
        time_array = 365.25*np.arange(0, 5)
        model_kwargs = {'time_array': time_array}  # time is given in days
        sm = SystemModel(model_kwargs=model_kwargs)

        # Add simple reservoir component
        sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))
        # Add observations of simple reservoir component model
        sres.add_obs('pressure')
        sres.add_obs('CO2saturation')

        # Run system model using current values of its parameters
        sm.forward()

        # Assign observations of the model to pressure and CO2saturation variables
        # Obtain pressure and CO2 saturation as lists
        pressure = sm.collect_observations_as_time_series(sres, 'pressure')
        CO2saturation = sm.collect_observations_as_time_series(sres, 'CO2saturation')

        # True values
        true_pressure = [9411325., 14019874.26175452, 14299528.64528817,
                         14463116.2847518, 14579183.64865103]
        true_saturation = [0., 0.21658361, 0.31286341, 0.3867415, 0.4490236]

        # Test with helpful message!
        for tp, p, t in zip(true_pressure, pressure, time_array):
            self.assertTrue(abs((tp-p)/tp) < 0.01,
                            'Pressure at time t={} days is {} Pa but should be {} Pa'
                            .format(str(t), str(p), str(tp)))
        for ts, s, t in zip(true_saturation[1:], CO2saturation[1:], time_array[1:]):
            self.assertTrue(abs((ts-s)/ts) < 0.01,
                            'Saturation at time t={} days is {} but should be {}'
                            .format(str(t), str(s), str(ts)))

    def test_stratigraphy(self):
        """Tests stratigraphy component.

        Tests whether parameters set in the stratigraphy component passed correctly
        to the subsequent components.
        """
        # Create system model 1 with both stratigraphy and reservoir component
        time_array = np.arange(0, 5)
        model_kwargs = {'time_array': time_array}  # time is given in days
        sm1 = SystemModel(model_kwargs=model_kwargs)

        # Add stratigraphy component
        strata = sm1.add_component_model_object(Stratigraphy(name='strata', parent=sm1))

        # Define parameters values
        num_shale_layers = 3
        true_shale_thickness = [100.0, 120.0, 140.0]
        true_aq_thickness = [50.0, 75.0]

        # Add parameters of stratigraphy component model: some of the parameters
        # are stochastic, some are deterministic
        strata.add_par('numberOfShaleLayers', value=num_shale_layers, vary=False)
        for ind in range(num_shale_layers):
            strata.add_par('shale{}Thickness'.format(ind+1),
                           value=true_shale_thickness[ind], vary=False)

        for ind in range(num_shale_layers-1):
            strata.add_par('aquifer{}Thickness'.format(ind+1),
                           # min and max differ from value by 10
                           min=true_aq_thickness[ind]-10.0,
                           max=true_aq_thickness[ind]+10.0,
                           value=true_aq_thickness[ind])

        # Add reservoir component
        sres1 = sm1.add_component_model_object(
            SimpleReservoir(name='sres1', parent=sm1))

        # Add deterministic parameters of reservoir component model
        # linked to the stratigraphy parameters
        sres1.add_par_linked_to_par(
            'numberOfShaleLayers',
            strata.deterministic_pars['numberOfShaleLayers'])
        for ind in range(num_shale_layers):
            par_name = 'shale{}Thickness'.format(ind+1)
            sres1.add_par_linked_to_par(par_name,
                                        strata.deterministic_pars[par_name])

        # Add stochastic parameters of reservoir component model
        # linked to the stratigraphy parameters
        for ind in range(num_shale_layers-1):
            par_name = 'aquifer{}Thickness'.format(ind+1)
            sres1.add_par_linked_to_par(par_name,
                                        strata.pars[par_name])

        # Add observations of reservoir component model
        sres1.add_obs('pressure')
        sres1.add_obs('CO2saturation')

        # Run system model 1 deterministically
        sm1.forward()

        # Collect pressure and saturation from system model 1
        pressure1 = sm1.collect_observations_as_time_series(
            sres1, 'pressure')
        CO2saturation1 = sm1.collect_observations_as_time_series(
            sres1, 'CO2saturation')

        # Create system model 2 with reservoir component only
        sm2 = SystemModel(model_kwargs=model_kwargs)

        # Add reservoir component
        sres2 = sm2.add_component_model_object(SimpleReservoir(name='sres2', parent=sm2))

        # Add deterministic parameters of reservoir component model
        # linked to the stratigraphy parameters
        sres2.add_par('numberOfShaleLayers', value=num_shale_layers, vary=False)
        for ind in range(num_shale_layers):
            sres2.add_par('shale{}Thickness'.format(ind+1),
                          value=true_shale_thickness[ind], vary=False)

        # Add stochastic parameters of reservoir component model
        # linked to the stratigraphy parameters
        for ind in range(num_shale_layers-1):
            sres2.add_par('aquifer{}Thickness'.format(ind+1),
                          # min and max differ from value by 10
                          min=true_aq_thickness[ind]-10.0,
                          max=true_aq_thickness[ind]+10.0,
                          value=true_aq_thickness[ind])

        # Add observations of reservoir component model
        sres2.add_obs('pressure')
        sres2.add_obs('CO2saturation')

        # Run system model 2 deterministically
        sm2.forward()

        # Collect pressure and saturation from system model 1
        pressure2 = sm2.collect_observations_as_time_series(
            sres2, 'pressure')
        CO2saturation2 = sm2.collect_observations_as_time_series(
            sres2, 'CO2saturation')

        true_pressure = [4854325., 7083141.10892644, 7362344.32621009,
                         7525781.52051495, 7641773.6512542]
        true_saturation = [0., 0., 0.00134371, 0.00520934, 0.00846822]

        # Compare observations of the second system model 2 to the true values
        for p2, tp, t in zip(pressure2, true_pressure, time_array):
            self.assertTrue(abs((p2-tp)/tp) < 0.001,
                            'Pressure at t={} days is {} Pa but should be {} Pa'
                            .format(t, p2, tp))

        for s2, ts, t in zip(CO2saturation2, true_saturation, time_array):
            # Small value at the bottom to avoid division by zero for zero saturation
            self.assertTrue(abs((s2-ts)/(ts+0.00001)) < 0.001,
                            'CO2 saturation at t={} days is {} but should be {}'
                            .format(t, s2, ts))

        # Compare pressure values between two system models
        for p1, p2, t in zip(pressure1, pressure2, time_array):
            self.assertTrue(abs((p1-p2)/p2) < 0.001,
                            'Pressure at t={} days is {} Pa but should be {} Pa'
                            .format(t, p1, p2))

        # Compare saturation values between two system models
        for s1, s2, t in zip(CO2saturation1, CO2saturation2, time_array):
            self.assertTrue(abs((s1-s2)/(s2+0.00001)) < 0.001,
                            'CO2 saturation at t={} days is {} but should be {}'
                            .format(t, s1, s2))

    def test_theis_reservoir_forward(self):
        """Tests Theis reservoir component.

        Tests the system model with a Theis reservoir component in a
        forward model against expected output for 11 years of data.
        """
        # Create system model
        time_array = np.array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.])
        model_kwargs = {'time_array': time_array}  # time is given in days
        sm = SystemModel(model_kwargs=model_kwargs)

        # Add Theis reservoir component
        rates = [5.0e-3, 6.0e-3, 7.0e-3, 8.0e-3, 0, 0, 0, 0, 0, 0, 0, 0,]
        tres = sm.add_component_model_object(
            TheisReservoir(name='tres', parent=sm,
                           injX=100, injY=100, locX=100, locY=45,
                           injTimes=time_array, injRates=rates))

        # Add observations of Theis reservoir component model
        tres.add_obs('pressure')

        # Set input parameters
        tres.add_par('initialPressure', value=1.0e6)       # Pa
        tres.add_par('reservoirThickness', value=30)       # m
        tres.add_par('logResPerm', value=-10.69897)        # log_(10) m^2
        tres.add_par('reservoirPorosity', value=0.2)
        tres.add_par('brineDensity', value=1000)           # kg/m^3
        tres.add_par('brineViscosity', value=2.535e-3)     # Pa*s
        tres.add_par('CO2Density', value=800)              # kg/m^3
        tres.add_par('compressibility', value=2.46e-9)     # 1/Pa

        # Run system model using current values of its parameters
        sm.forward()

        # Assign observations of the model to pressure and CO2saturation variables
        # Obtain pressure and CO2 saturation as lists
        pressure = sm.collect_observations_as_time_series(tres, 'pressure')

        # True values
        true_pressure = [
            1000000., 1007119.96734507, 1009474.2473876, 1011628.95843403,
            1013734.5164918, 1003014.77979786, 1002017.61984966, 1001539.69430367,
            1001251.12645661, 1001055.98536937, 1000914.54243831, 1000807.03684403]

        # Test with helpful message!
        for tp, p, t in zip(true_pressure, pressure, time_array):
            self.assertTrue(abs((tp-p)/tp) < 0.01,
                            'Pressure at time t={} days is {} Pa but should be {} Pa'
                            .format(str(t), str(p), str(tp)))


BASE_TESTS = [
    'test_openiam_cf',
    # Do not move test_openiam_cf test anywhere down the list:
    # one of the tests below causes it to fail
    'test_alluvium_aquifer',
    'test_analytical_reservoir_forward',
    'test_atm',
    'test_carb_aquifer',
    'test_cemented_wellbore_forward',
    'test_cemented_wellbore_wr_forward',
    'test_chemical_sealing_coupled_simple_reservoir',
    'test_chemical_sealing_not_seal_forward',
    'test_chemical_sealing_seal_forward',
    'test_coupled_analytical_reservoir_ms_wellbore_forward',
    'test_coupled_components_with_grid_obs',
    'test_coupled_reservoir_cemented_wellbore_forward',
    'test_coupled_reservoir_ms_wellbore_forward',
    'test_coupled_reservoir_open_carbonate_forward',
    'test_coupled_reservoir_open_wellbore_forward',
    'test_deep_alluvium_aquifer',
    'test_fault_flow',
    'test_fault_leakage',
    'test_futuregen_aquifer',
    'test_futuregen_azmi',
    'test_generalized_flow_rate_cmpnt',
    'test_generic_aquifer',
    'test_generic_reservoir_forward',
    'test_hydrocarbon_leakage_forward',
    'test_lhs',
    'test_lookup_table_reservoir',
    'test_lookup_table_reservoir_h5',
    'test_lookup_table_reservoir_2d',
    'test_lookup_table_reservoir_2d_h5',
    'test_lookup_table_reservoir_3d',
    'test_lookup_table_reservoir_3d_h5',
    'test_parstudy',
    'test_plume_stability_cmpnt',
    'test_plume_stability_cmpnt_h5',
    'test_rate_to_mass_adapter',
    'test_reservoir_data_interpolator',
    'test_seal_horizon',
    'test_sh_thickness_sampler',
    'test_sh_permeability_sampler',
    'test_simple_reservoir_forward',
    'test_stratigraphy',
    'test_theis_reservoir_forward',
    ]

KERAS_TESTS = [
    'test_alluvium_aquifer_lf',
    'test_deep_alluvium_aquifer_ml',
    'test_fault_leakage',
    'test_generic_aquifer',
    'test_generic_reservoir_forward',
    'test_hydrocarbon_leakage_forward',
    'test_kimberlina_wellbore'
    ]


def suite(case):
    """Determines set of tests to be performed."""

    suite = unittest.TestSuite()

    if case == 'all':
        for test_nm in BASE_TESTS:
            suite.addTest(Tests(test_nm))
        for test_nm in KERAS_TESTS:
            suite.addTest(Tests(test_nm))
        return suite
    if case == 'base':
        for test_nm in BASE_TESTS:
            suite.addTest(Tests(test_nm))
        return suite
    if case == 'keras':
        for test_nm in KERAS_TESTS:
            suite.addTest(Tests(test_nm))
        return suite

    return suite

if __name__ == '__main__':
    # For multiprocessing in Spyder
    __spec__ = None

    if len(sys.argv) > 1:
        case = sys.argv[1]
    else:
        case = 'base'

    # Second keyword argument makes testing messages printed
    # when running test suite in IPython console
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stderr)
    test_suite = suite(case)
    result = runner.run(test_suite)
    if not result.wasSuccessful():
        sys.exit(1)
