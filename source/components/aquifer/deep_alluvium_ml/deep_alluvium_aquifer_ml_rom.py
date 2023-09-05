# -*- coding: utf-8 -*-
"""
The module contains the solution class for the deep alluvium aquifer ROM
based on the Kimberlina Synthetic Model dataset

Authors: class Kayyum Mansoor

Date: 07/17/2018
Last modified: 07/27/2018

outputs:
TDSOutputs: affected volume, length, width and height of TDS plume with change
greater than 100 mg/L from baseline
PressureOutputs: affected volume, length, width and height of pressure plume
with change greater than 500 Pascals from baseline
pHOutputs: affected volume, length, width and height of plume with pH = 6.75


Component model input definitions:

* **simtime** * [years] (0 to 200) - Simulation time (default 20.0 years)

* **brine_rate** * [kg/s] (0 to 0.017) - Brine flux (default 0.00030 kg/s)

* **brine_mass** * [|log10| (kg)] (2.337 to 6.939) - Cumulative brine mass
(default 4.928 |log10| (kg))

* **co2_rate** * [kg/s] (0 to 0.385) - |CO2| flux (default 0.045 kg/s)

* **co2_mass** * [|log10| (kg)] (0.001 to 9.210) - Cumulative |CO2| mass
  (default 7.214 |log10| (kg))

* **logK_sand1** * [|log10| (|m^2|)] (-12.92 to -10.92) - Permeability
  in layer 1 10-546 m (default -11.91 |log10| (|m^2|))

* **logK_sand2** * [|log10| (|m^2|)] (-12.72 to -10.72) - Permeability
  in layer 2 546-1225 m (default -11.71 |log10| (|m^2|))

* **logK_sand3** * [|log10| (|m^2|)] (-12.70 to -10.70) - Permeability
  in layer 3 1225-1411 m (default -11.69 |log10| (|m^2|))

* **logK_caprock** * [|log10| (|m^2|)] (-16.70 to -14.70) - Caprock permeability
  0-5 m (default -15.70 |log10| (|m^2|))

* **correlationLengthX** * [m] (200 to 2000 m) - Correlation length
  in x-direction (default 1098.235 m)

* **correlationLengthZ** * [m] (10 to 150 m) - Correlation length
  in z-direction (default 79.827 m)

* **sandFraction** * [-] (0.70 to 0.90) - Sand volume fraction (default 0.800)

* **groundwater_gradient** * [-] (0.0010 to 0.0017) - Regional groundwater
  gradient (dh/dx=change in hydraulic head/distance) (default 0.0013)

* **leak_depth** * [-] (424.4 to 1341.5 m) - Depth of leakage interval (default 883.3 m)

Observations from the Deep Alluvium Aquifer component are:

* **TDS_volume** [|m^3|] - volume of plume above baseline TDS change in mg/L
  (change in TDS > 100 mg/L).
* **TDS_dx** [|m|] - length of plume above baseline TDS change in mg/L
  (change in TDS > 100 mg/L).
* **TDS_dy** [|m|] - width of plume above baseline TDS change in mg/L
  (change in TDS > 100 mg/L).
* **TDS_dz** [|m|] - height of plume above baseline TDS change in mg/L
  (change in TDS > 100 mg/L).

* **Pressure_volume** [|m^3|] - volume of plume above baseline pressure change
  in Pa (change in pressure > 500 Pa).
* **Pressure_dx** [|m|] - length of plume above baseline pressure change in Pa
  (change in pressure > 500 Pa).
* **Pressure_dy** [|m|] - width of plume above baseline pressure change in Pa
  (change in pressure > 500 Pa).
* **Pressure_dz** [|m|] - height of plume above baseline pressure change in Pa
  (change in pressure > 500 Pa).

* **pH_volume** [|m^3|] - volume of plume below pH threshold (pH < 6.75).
* **pH_dx** [|m|] - length of plume below pH threshold (pH < 6.75).
* **pH_dy** [|m|] - width of plume below pH threshold (pH < 6.75).
* **pH_dz** [|m|] - height of plume below pH threshold (pH < 6.75).

"""
# Created on: June 10, 2019
# @author: Kayyum Mansoor
# email: mansoor1@llnl.gov
import numpy as np
import os

DAA_ML_FILE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)))

from tensorflow.keras.models import load_model

DAA_MODELS = {
    'tds': load_model(os.path.join(DAA_ML_FILE_PATH, 'keras.TDS100_voldxdydz.h5')),
    'ph': load_model(os.path.join(DAA_ML_FILE_PATH, 'keras.pH675_voldxdydz.h5')),
    'p': load_model(os.path.join(DAA_ML_FILE_PATH, 'keras.P500_voldxdydz.h5'))}

def daa_model(name):
    return DAA_MODELS[name]


class Solution(object):
    """ Solution class for Deep Alluvium Aquifer ML ROM. """
    def __init__(self):
        """ Create an instance of the Solution class. """
        self.Outputs = None

        self.pars_bounds = dict()
        self.pars_bounds['log_time'] = [np.log10(1), np.log10(200)]
        self.pars_bounds['brine_rate'] = [3.4951E-07, 0.017002313]
        self.pars_bounds['log_brine_rate'] = [np.log10(3.4951E-07), np.log10(0.017002313)]
        self.pars_bounds['brine_mass'] = [10**1.37525, 10**6.93934]
        self.pars_bounds['log_brine_mass'] = [1.37525, 6.93934]

        self.pars_bounds['co2_rate'] = [5.000E-11, 0.384984283]
        self.pars_bounds['log_co2_rate'] = [np.log10(5.000E-11), np.log10(0.384984283)]
        self.pars_bounds['co2_mass'] = [10**-2.80193, 10**9.20977]
        self.pars_bounds['log_co2_mass'] = [-2.80193, 9.20977]

        self.pars_bounds['logK_sand1'] = [-12.92, -10.92]
        self.pars_bounds['logK_sand2'] = [-12.72, -10.72]
        self.pars_bounds['logK_sand3'] = [-12.7, -10.7]
        self.pars_bounds['logK_caprock'] = [-16.699, -14.699]

        self.pars_bounds['log_correlationLengthX'] = [np.log10(200), np.log10(2000)]
        self.pars_bounds['log_correlationLengthZ'] = [np.log10(10), np.log10(150)]
        self.pars_bounds['correlationLengthX'] = [200, 2000]
        self.pars_bounds['correlationLengthZ'] = [10, 150]

        self.pars_bounds['sandFraction'] = [0.7, 0.9]

        self.pars_bounds['groundwater_gradient'] = [0.001000, 0.001667]
        self.pars_bounds['log_leak_depth'] = [np.log10(424.36), np.log10(1341.48)]
        self.pars_bounds['leak_depth'] = [424.36, 1341.48]

        self.pars_bounds['tds_log_vol'] = [1.022833, 9.418434]
        self.pars_bounds['tds_log_dx'] = [1.003280, 3.938895]
        self.pars_bounds['tds_log_dy'] = [1.003280, 3.694605]
        self.pars_bounds['tds_log_dz'] = [1.004485, 3.148661]

        self.pars_bounds['ph_log_vol'] = [1.025828, 10.775760]
        self.pars_bounds['ph_log_dx'] = [1.000000, 4.003784]
        self.pars_bounds['ph_log_dy'] = [1.000000, 3.706504]
        self.pars_bounds['ph_log_dz'] = [1.000000, 3.148081]

        self.pars_bounds['p_log_vol'] = [1.014235, 10.608166]
        self.pars_bounds['p_log_dx'] = [1.000000, 3.972550]
        self.pars_bounds['p_log_dy'] = [1.000000, 3.693507]
        self.pars_bounds['p_log_dz'] = [1.000000, 3.148852]

    def find(self, inputArray1):
        """ Find output of the ROM for the provided input parameters."""

        self.Outputs = np.ones(12)*1e-10
        self.simtime = inputArray1[0]

        if inputArray1[0] > 0:

            self.brine_rate = inputArray1[1]
            self.log_brine_mass = inputArray1[2]
            self.co2_rate = inputArray1[3]
            self.log_co2_mass = inputArray1[4]
            self.logK_sand1 = inputArray1[5]
            self.logK_sand2 = inputArray1[6]
            self.logK_sand3 = inputArray1[7]
            self.logK_caprock = inputArray1[8]
            self.correlationLengthX = inputArray1[9]
            self.correlationLengthZ = inputArray1[10]
            self.sandFraction = inputArray1[11]
            self.groundwater_gradient = inputArray1[12]
            self.leak_depth = inputArray1[13]

            self.norm_time = np.interp(np.log10(self.simtime),
                                       (self.pars_bounds['log_time']),
                                       (0.0, 1.0))
            self.norm_log_brine_rate = np.interp(np.log10(self.brine_rate),
                                                 (self.pars_bounds['log_brine_rate']),
                                                 (0.0, 1.0))
            self.norm_log_brine_mass = np.interp(np.log10(self.log_brine_mass),
                                                 (self.pars_bounds['log_brine_mass']),
                                                 (0.0, 1.0))

            self.norm_log_co2_rate = np.interp(np.log10(self.co2_rate),
                                               (self.pars_bounds['log_co2_rate']),
                                               (0.0, 1.0))
            self.norm_log_co2_mass = np.interp(np.log10(self.log_co2_mass),
                                               (self.pars_bounds['log_co2_mass']),
                                               (0.0, 1.0))

            self.norm_logK_sand1 = np.interp(self.logK_sand1,
                                             (self.pars_bounds['logK_sand1']),
                                             (0.0, 1.0))
            self.norm_logK_sand2 = np.interp(self.logK_sand2,
                                             (self.pars_bounds['logK_sand2']),
                                             (0.0, 1.0))
            self.norm_logK_sand3 = np.interp(self.logK_sand3,
                                             (self.pars_bounds['logK_sand3']),
                                             (0.0, 1.0))
            self.norm_logK_caprock = np.interp(self.logK_caprock,
                                               (self.pars_bounds['logK_caprock']),
                                               (0.0, 1.0))

            self.norm_correlationLengthX = np.interp(
                np.log10(self.correlationLengthX),
                (self.pars_bounds['log_correlationLengthX']),
                (0.0, 1.0))
            self.norm_correlationLengthZ = np.interp(
                np.log10(self.correlationLengthZ),
                (self.pars_bounds['log_correlationLengthZ']),
                (0.0, 1.0))

            self.norm_sandFraction = np.interp(self.sandFraction,
                                               (self.pars_bounds['sandFraction']),
                                               (0.0, 1.0))

            self.norm_groundwater_gradient = np.interp(
                self.groundwater_gradient,
                (self.pars_bounds['groundwater_gradient']),
                (0.0, 1.0))

            self.norm_leak_depth = np.interp(np.log10(self.leak_depth),
                                             (self.pars_bounds['log_leak_depth']),
                                             (0.0, 1.0))

            self.input_array = [
                self.norm_time, self.norm_log_brine_rate, self.norm_log_brine_mass,
                self.norm_log_co2_rate, self.norm_log_co2_mass,
                self.norm_logK_sand1, self.norm_logK_sand2, self.norm_logK_sand3,
                self.norm_logK_caprock,
                self.norm_correlationLengthX, self.norm_correlationLengthZ,
                self.norm_sandFraction, self.norm_groundwater_gradient,
                self.norm_leak_depth]

            import tensorflow as tf
            self.input_array = tf.convert_to_tensor(np.reshape(self.input_array, (1, 14)))
            self.model_tds = (daa_model('tds')(self.input_array, training=False).numpy()[0])
            self.model_ph = (daa_model('ph')(self.input_array, training=False).numpy()[0])
            self.model_p = (daa_model('p')(self.input_array, training=False).numpy()[0])

            self.Outputs[0] = (np.interp(self.model_tds[0], (0.0, 1.0),
                                         (self.pars_bounds['tds_log_vol'])))
            self.Outputs[1] = (np.interp(self.model_tds[1], (0.0, 1.0),
                                         (self.pars_bounds['tds_log_dx'])))
            self.Outputs[2] = (np.interp(self.model_tds[2], (0.0, 1.0),
                                         (self.pars_bounds['tds_log_dy'])))
            self.Outputs[3] = (np.interp(self.model_tds[3], (0.0, 1.0),
                                         (self.pars_bounds['tds_log_dz'])))

            self.Outputs[4] = (np.interp(self.model_ph[0], (0.0, 1.0),
                                         (self.pars_bounds['ph_log_vol'])))
            self.Outputs[5] = (np.interp(self.model_ph[1], (0.0, 1.0),
                                         (self.pars_bounds['ph_log_dx'])))
            self.Outputs[6] = (np.interp(self.model_ph[2], (0.0, 1.0),
                                         (self.pars_bounds['ph_log_dy'])))
            self.Outputs[7] = (np.interp(self.model_ph[3], (0.0, 1.0),
                                         (self.pars_bounds['ph_log_dz'])))

            self.Outputs[8] = (np.interp(self.model_p[0], (0.0, 1.0),
                                         (self.pars_bounds['p_log_vol'])))
            self.Outputs[9] = (np.interp(self.model_p[1], (0.0, 1.0),
                                         (self.pars_bounds['p_log_dx'])))
            self.Outputs[10] = (np.interp(self.model_p[2], (0.0, 1.0),
                                          (self.pars_bounds['p_log_dy'])))
            self.Outputs[11] = (np.interp(self.model_p[3], (0.0, 1.0),
                                          (self.pars_bounds['p_log_dz'])))
            self.Outputs = 10**(self.Outputs)


if __name__ == "__main__":

    simtime = 50
    dt = 365.25
    brine_rate = 6.24E-05
    brine_mass = 2646900   # 10**6.422737534
    co2_rate = 0.150771
    co2_mass = 171604000   # 10**8.234527407
    logK_sand1 = -11.3845
    logK_sand2 = -11.9252
    logK_sand3 = -10.8862
    logK_caprock = -14.731
    correlationLengthX = 520.721
    correlationLengthZ = 112.442
    sandFraction = 0.743
    groundwater_gradient = 0.00105408
    leak_depth = 715.99

    inputArray1 = [simtime, brine_rate, brine_mass, co2_rate, co2_mass,
                   logK_sand1, logK_sand2, logK_sand3, logK_caprock,
                   correlationLengthX, correlationLengthZ, sandFraction,
                   groundwater_gradient, leak_depth]

    sol = Solution()
    sol.find(inputArray1)

    labels = ['TDS_volume', 'TDS_dx', 'TDS_dy', 'TDS_dz',
              'Pressure_volume', 'Pressure_dx', 'Pressure_dy', 'Pressure_dz',
              'pH_volume', 'pH_dx', 'pH_dy', 'pH_dz']

    results = sol.Outputs

    print('results:')
    for i, v in enumerate(labels):
        print('%-15s  %14.2f'%(v, results[i]))

    # The answer should be
    # results:
    # TDS_volume         100467876.61
    # TDS_dx                  1061.26
    # TDS_dy                   719.13
    # TDS_dz                   807.42
    # Pressure_volume     41944085.19
    # Pressure_dx             1130.33
    # Pressure_dy              592.68
    # Pressure_dz              738.88
    # pH_volume           81631821.69
    # pH_dx                   1152.09
    # pH_dy                    968.55
    # pH_dz                    824.11
