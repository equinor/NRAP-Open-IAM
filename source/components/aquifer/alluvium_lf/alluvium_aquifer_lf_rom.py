# -*- coding: utf-8 -*-
"""
The module contains the solution class for the alluvium aquifer (low flux) ROM
based on the High Plains Aquifer (rev 2019) Synthetic Model dataset

Authors: class Kayyum Mansoor

Date: 07/03/2019
Last modified: 07/03/2019

outputs:
TDS: affected volume of aquifer above TDS threshold in mg/L (TDS > 1300 mg/L).
pH: affected volume of aquifer below pH threshold (pH < 7.0 ).
As: affected volume of aquifer above arsenic threshold in μg/L (arsenic > 9.3 μg/L).
Ba: affected volume of aquifer above barium threshold in μg/L (barium > 140 μg L/L).
Cd: affected volume of aquifer above cadmium threshold in μg/L (cadmium > 0.25 μg/L).
Pb: affected volume of aquifer above lead threshold in μg/L (lead  > 0.63 μg/L).
Benzene: affected volume of aquifer above benzene threshold (benzene  > 0.03 μg/L).
Naphthalene: affected volume of aquifer above naphthalene threshold (napthalene  > 0.2 μg/L).
Phenol: affected volume of aquifer above phenol threshold (phenol  > 0.003 μg/L).


Component model input definitions:

* **simtime** [years] (1 to 200) - Simulation time (default 5.0 years)

* **brn_rate** [kg/s] (0 to 0.00075) - Brine flux (default 0.000145878 kg/s)

* **brn_mass** [|log10| (kg)] (3.20 to 6.23) - Cumulative brine mass
  (default 4.36206 |log10| (kg))

* **co2_rate** [kg/s] (0 to 0.005001) - CO2 flux (default 0.000006883 kg/s)

* **co2_mass** [|log10| (kg)] (0.53 to 6.75) - Cumulative CO2 mass
  (default 2.73481 |log10| (kg))

* **sandFraction** [-] (0.60 to 0.90) - Sand volume fraction (default 0.7087)

* **correlationLengthX** [m] (200 to 2500) - Correlation length in  X-direction
  (default 473.9740 m)

* **correlationLengthZ** [m] (0.5 to 25) - Correlation length in  Z-direction
  (default 22.5966 m)

* **permeabilitySand** [|m^2|] (-14 to -10) - Permeability [|m^2|]
  (default -13.3634  |m^2|)

* **permeabilityClay** [|log10| (|m^2|)] (-18 to -15) - Permeability in  clay
  (default -15.5075 |log10| (|m^2|))

* **NaMolality** [|log10| (Molality)] (-3 to 1) - Sodium molality
  (default 0.9279  |log10| (Molality))

* **PbMolality** [|log10| (Molality)] (-9.1836 to -5.6836) - Lead molality
  (default  -6.2021 |log10| (Molality))

* **benzeneMolality** [|log10| (Molality)] (-10 to -6) - Benzene molality
  (default  -7.7138 |log10| (Molality))

* **tMitigation** [years] (1 to 100) - Mitigation time(default 7.7387 years)

Observations from the Alluvium Aquifer (low flux) component are:

* **TDS** [|m^3|] - volume of aquifer above TDS threshold in mg/L (TDS > 1300 mg/L).
* **pH** [|m^3|] - volume of aquifer below pH threshold (pH < 7.0 ).
* **As** [|m^3|] - volume of aquifer above arsenic threshold in μg/L (arsenic > 9.3 μg/L).
* **Ba** [|m^3|] - volume of aquifer above barium threshold in μg/L (barium > 140 μg L/L).
* **Cd** [|m^3|] - volume of aquifer above cadmium threshold in μg/L (cadmium > 0.25 μg/L).
* **Pb** [|m^3|] - volume of aquifer above lead threshold in μg/L (lead  > 0.63 μg/L).
* **Benzene** [|m^3|] - volume of aquifer above benzene threshold (benzene  > 0.03 μg/L).
* **Naphthalene** [|m^3|] - volume of aquifer above naphthalene threshold (napthalene  > 0.2 μg/L).
* **Phenol** [|m^3|] - volume of aquifer above phenol threshold (phenol  > 0.003 μg/L).

"""
# Created on: July 3, 2019
# @author: Kayyum Mansoor
# email: mansoor1@llnl.gov
import numpy as np
import os

AA_LF_FILE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)))

from tensorflow.keras.models import load_model

AA_MODELS = {
    'tds': load_model(os.path.join(AA_LF_FILE_PATH, 'rom_tds.h5')),
    'ph': load_model(os.path.join(AA_LF_FILE_PATH, 'rom_ph.h5')),
    'as': load_model(os.path.join(AA_LF_FILE_PATH, 'rom_as.h5')),
    'ba': load_model(os.path.join(AA_LF_FILE_PATH, 'rom_ba.h5')),
    'cd': load_model(os.path.join(AA_LF_FILE_PATH, 'rom_cd.h5')),
    'pb': load_model(os.path.join(AA_LF_FILE_PATH, 'rom_pb.h5')),
    'obz': load_model(os.path.join(AA_LF_FILE_PATH, 'rom_obz.h5')),
    'onp': load_model(os.path.join(AA_LF_FILE_PATH, 'rom_onp.h5')),
    'opl': load_model(os.path.join(AA_LF_FILE_PATH, 'rom_opl.h5'))}

def aa_model(name):
    return AA_MODELS[name]


class Solution(object):

    def __init__(self):
        """ Create an instance of the Solution class. """
        self.Outputs = None

        self.pars_bounds = dict()
        self.pars_bounds['log_time'] = [np.log10(1.0), np.log10(200.0)]
        # 0.0000002272 to 0.0007500000 kg/s
        self.pars_bounds['log_brn_flux'] = [-6.643679611, -3.124938737]
        self.pars_bounds['log_brn_mass'] = [3.20, 6.23]
        # 0.000000002 to 0.005000000 kg/s
        self.pars_bounds['log_co2_flux'] = [-8.704490, -2.301030]
        self.pars_bounds['log_co2_mass'] = [0.53, 6.75]
        self.pars_bounds['sandFraction'] = [0.60, 0.90]
        self.pars_bounds['log_correlationLengthX'] = [np.log10(200.0), np.log10(2500.0)]
        self.pars_bounds['log_correlationLengthZ'] = [np.log10(0.5), np.log10(25.0)]
        self.pars_bounds['permeabilitySand'] = [-14.0, -10.0]
        self.pars_bounds['permeabilityClay'] = [-18.0, -15.0]
        self.pars_bounds['NaMolality'] = [-3.0, 1.0]
        self.pars_bounds['PbMolality'] = [-9.1836, -5.6836]
        self.pars_bounds['benzeneMolality'] = [-10, -6]
        self.pars_bounds['log_tMitigation'] = [np.log10(1.0), np.log10(100.0)]

        self.pars_bounds['tds_log_vol'] = [1.000000, 6.075595e+00]
        self.pars_bounds['ph_log_vol'] = [1.000000, 8.840981e+00]
        self.pars_bounds['as_log_vol'] = [1.000000, 4.763311e+00]
        self.pars_bounds['ba_log_vol'] = [1.000000, 6.707547e+00]
        self.pars_bounds['cd_log_vol'] = [1.000000, 5.808999e+00]
        self.pars_bounds['pb_log_vol'] = [1.000000, 6.352865e+00]
        self.pars_bounds['obz_log_vol'] = [1.000000, 7.011281e+00]
        self.pars_bounds['onp_log_vol'] = [1.000000, 6.898353e+00]
        self.pars_bounds['opl_log_vol'] = [1.000000, 7.256172e+00]

    def find(self, inputArray1):
        """ Find output of the ROM for the provided input parameters."""
        self.Outputs = np.ones(9)*1e-10

        self.simtime = inputArray1[0]

        if inputArray1[0] > 0:
            self.brine_rate = inputArray1[1]
            self.brine_mass = inputArray1[2]
            self.co2_rate = inputArray1[3]
            self.co2_mass = inputArray1[4]

            self.sandFraction = inputArray1[5]
            self.correlationLengthX = inputArray1[6]
            self.correlationLengthZ = inputArray1[7]
            self.permeabilitySand = inputArray1[8]
            self.permeabilityClay = inputArray1[9]
            self.NaMolality = inputArray1[10]
            self.PbMolality = inputArray1[11]
            self.benzeneMolality = inputArray1[12]
            self.tMitigation = inputArray1[13]

            self.norm_time = np.interp(np.log10(self.simtime),
                                       (self.pars_bounds['log_time']), (0.0, 1.0))

            self.norm_log_brine_rate = np.interp(np.log10(self.brine_rate),
                                                 (self.pars_bounds['log_brn_flux']),
                                                 (0.0, 1.0))
            self.norm_log_brine_mass = np.interp(np.log10(self.brine_mass),
                                                 (self.pars_bounds['log_brn_mass']),
                                                 (0.0, 1.0))
            self.norm_log_co2_rate = np.interp(np.log10(self.co2_rate),
                                               (self.pars_bounds['log_co2_flux']),
                                               (0.0, 1.0))
            self.norm_log_co2_mass = np.interp(np.log10(self.co2_mass),
                                               (self.pars_bounds['log_co2_mass']),
                                               (0.0, 1.0))

            self.norm_sandFraction = np.interp(self.sandFraction,
                                               (self.pars_bounds['sandFraction']),
                                               (0.0, 1.0))
            self.norm_log_correlationLengthX = np.interp(
                np.log10(self.correlationLengthX),
                (self.pars_bounds['log_correlationLengthX']), (0.0, 1.0))
            self.norm_log_correlationLengthZ = np.interp(
                np.log10(self.correlationLengthZ),
                (self.pars_bounds['log_correlationLengthZ']), (0.0, 1.0))

            self.norm_permeabilitySand = np.interp(self.permeabilitySand,
                                                   (self.pars_bounds['permeabilitySand']),
                                                   (0.0, 1.0))
            self.norm_permeabilityClay = np.interp(self.permeabilityClay,
                                                   (self.pars_bounds['permeabilityClay']),
                                                   (0.0, 1.0))
            self.norm_NaMolality = np.interp(self.NaMolality,
                                             (self.pars_bounds['NaMolality']),
                                             (0.0, 1.0))
            self.norm_PbMolality = np.interp(self.PbMolality,
                                             (self.pars_bounds['PbMolality']),
                                             (0.0, 1.0))
            self.norm_BenzeneMolality = np.interp(self.benzeneMolality,
                                                  (self.pars_bounds['benzeneMolality']),
                                                  (0.0, 1.0))

            self.norm_log_tMitigation = np.interp(np.log10(self.tMitigation),
                                                  (self.pars_bounds['log_tMitigation']),
                                                  (0.0, 1.0))

            # Setup input array for the model
            self.input_array = [
                self.norm_time, self.norm_log_brine_rate, self.norm_log_brine_mass,
                self.norm_log_co2_rate, self.norm_log_co2_mass, self.norm_sandFraction,
                self.norm_log_correlationLengthX, self.norm_log_correlationLengthZ,
                self.norm_permeabilitySand, self.norm_permeabilityClay,
                self.norm_NaMolality, self.norm_PbMolality, self.norm_BenzeneMolality,
                self.norm_log_tMitigation]

            import tensorflow as tf
            self.input_array = tf.convert_to_tensor(np.reshape(self.input_array, (1, 14)))
            self.model_tds = (aa_model('tds')(self.input_array, training=False).numpy()[0])
            self.model_ph = (aa_model('ph')(self.input_array, training=False).numpy()[0])
            self.model_as = (aa_model('as')(self.input_array, training=False).numpy()[0])
            self.model_ba = (aa_model('ba')(self.input_array, training=False).numpy()[0])
            self.model_cd = (aa_model('cd')(self.input_array, training=False).numpy()[0])
            self.model_pb = (aa_model('pb')(self.input_array, training=False).numpy()[0])
            self.model_obz = (aa_model('obz')(self.input_array, training=False).numpy()[0])
            self.model_onp = (aa_model('onp')(self.input_array, training=False).numpy()[0])
            self.model_opl = (aa_model('opl')(self.input_array, training=False).numpy()[0])

            self.Outputs[0] = (np.interp(
                self.model_tds[0], (0.0, 1.0), (self.pars_bounds['tds_log_vol'])))
            self.Outputs[1] = (np.interp(
                self.model_ph[0], (0.0, 1.0), (self.pars_bounds['ph_log_vol'])))
            self.Outputs[2] = (np.interp(
                self.model_as[0], (0.0, 1.0), (self.pars_bounds['as_log_vol'])))
            self.Outputs[3] = (np.interp(
                self.model_ba[0], (0.0, 1.0), (self.pars_bounds['ba_log_vol'])))
            self.Outputs[4] = (np.interp(
                self.model_cd[0], (0.0, 1.0), (self.pars_bounds['cd_log_vol'])))
            self.Outputs[5] = (np.interp(
                self.model_pb[0], (0.0, 1.0), (self.pars_bounds['pb_log_vol'])))
            self.Outputs[6] = (np.interp(
                self.model_obz[0], (0.0, 1.0), (self.pars_bounds['obz_log_vol'])))
            self.Outputs[7] = (np.interp(
                self.model_onp[0], (0.0, 1.0), (self.pars_bounds['onp_log_vol'])))
            self.Outputs[8] = (np.interp(
                self.model_opl[0], (0.0, 1.0), (self.pars_bounds['opl_log_vol'])))
            self.Outputs = 10**(self.Outputs)


if __name__ == "__main__":
    simtime = 5.0
    brine_rate = 0.000145878 # kg/s
    brine_mass = 10**4.36206 # kg
    co2_rate = 0.000006883 # kg/s
    co2_mass = 10**2.73481 # kg
    sandFraction = 0.7087
    correlationLengthX = 473.9740
    correlationLengthZ = 22.5966
    permeabilitySand = -13.3634
    permeabilityClay = -15.5075
    NaMolality = 0.9279
    PbMolality = -6.2021
    benzeneMolality = -7.7138
    tMitigation = 7.7387

    inputArray1 = [simtime, brine_rate, brine_mass, co2_rate, co2_mass,
                   sandFraction, correlationLengthX, correlationLengthZ,
                   permeabilitySand, permeabilityClay,
                   NaMolality, PbMolality, benzeneMolality, tMitigation]

    sol = Solution()
    sol.find(inputArray1)

    labels = ['TDS volume', 'pH volume', 'As volume', 'Ba volume',
              'Cd volume', 'Pb volume', 'Benzene volume', 'Naphthalene volume',
              'Phenol volume']

    results = sol.Outputs

    print('results:')
    for i, v in enumerate(labels):
        print('%-20s  %14.2f'%(v, results[i]))

    # The answer should be
    # results:
    # TDS volume                   2817.58
    # pH volume                     678.17
    # As volume                      97.35
    # Ba volume                   24999.17
    # Cd volume                    8258.10
    # Pb volume                   13643.12
    # Benzene volume               2131.52
    # Naphthalene volume           3103.14
    # Phenol volume               12014.18
