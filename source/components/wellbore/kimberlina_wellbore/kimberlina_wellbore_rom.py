# -*- coding: utf-8 -*-
"""
The module contains the solution class for the deep Kimberlina reservoir
based on the Kimberlina wellbore model

Authors: class Kayyum Mansoor

Date: 11/12/2018
Last modified: 11/12/2018

outputs:
BrineMass
CO2Mass

    norm_pressure = 0.08825
    norm_co2_saturation = 0.02199
    norm_logK_well_permeability = 0.13772
    norm_previous_log10_co2_mass = 0.09542
    norm_previous_log10_brine_mass = 0.09542

Component model input definitions:

* **norm_simtime** * [years] (0.0 - 1.0 | full range: 0 to 200) - simulation time
  (default 20.0 years)

* **norm_pressure** * [Pa] (0.0 - 1.0 | full range: 17130510 to 37483630)
  - bottom hole pressure (default 28228195.8 (Pa))

* **norm_co2_saturation** * [ ] (0.0 - 1.0 | full range: 0.0 to 0.634324) -
  bottom home CO2 saturation (default 0.45741 (-))

* **norm_logK_well_permeability** * [|log10| (|m^2|)] (0.0 - 1.0 | full range:
  -11.8468 to -10.001) - well permeability 0-5 m (default -11.5926 |log10| (|m^2|))

* **norm_previous_log10_co2_mass** * [|log10| (kg)] (0.0 - 1.0 | full range:
  -2.80 to 9.21) - cumulative brine mass from previous timestep
  (default 5.54 |log10| (kg))

* **norm_previous_log10_brine_mass** * [|log10| (kg)] (0.0 - 1.0 | full range:
  1.37 to 6.94) - cumulative CO2 mass from previous timestep
  (default 3.44 |log10| (kg))

Observations from the Kimberlina Wellbore ROM component are:

* **norm_log10_co2_mass** [|log10| (kg)] - Cumulative brine mass).
* **norm_log10_brine_mass** [|log10| (kg)] - Cumulative CO2 mass).
"""
# Created on: Nov 13, 2018
# @author: Kayyum Mansoor
# email: mansoor1@llnl.gov
import numpy as np
import os

KW_FILE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)))

from tensorflow.keras.models import load_model

KW_MODELS = {
    'co2': load_model(os.path.join(KW_FILE_PATH, 'keras.SYN_co2Mass.h5')),
    'brine': load_model(os.path.join(KW_FILE_PATH, 'keras.SYN_brnMass.h5'))}

def kw_model(name):
    return KW_MODELS[name]


class Solution(object):
    """ Solution class for Kimberlina wellbore ROM. """
    def __init__(self):
        """ Create an instance of the Solution class. """
        self.Outputs = None

        self.pars_bounds = dict()
        self.pars_bounds['pressure'] = [17130510, 37483630]
        self.pars_bounds['co2_saturation'] = [0, 0.634324]
        self.pars_bounds['logK_well_permeability'] = [-11.8468, -10.001]

        self.pars_bounds['co2_mass'] = [10**(-2.80192602858), 10**(9.20977325986)]
        self.pars_bounds['brn_mass'] = [10**(1.37524597484), 10**(6.93933611696)]

        self.pars_bounds['log10_co2_mass'] = [-2.80192602858, 9.20977325986]
        self.pars_bounds['log10_brn_mass'] = [1.37524597484, 6.93933611696]

    def find(self, inputArray1, inputArray2):
        """ Find output of the ROM for the provided input parameters."""

        self.Outputs = np.ones(4)*1e-10
        self.simtime, self.dt = inputArray1[4], inputArray1[5]

        if self.simtime > 0:

            self.pressure, self.co2_saturation, self.logK_well_permeability = (
                inputArray1[0], inputArray1[1], inputArray1[2])
            self._previous_co2_mass, self._previous_brine_mass = (
                inputArray1[3], inputArray2[3])

            if self._previous_co2_mass <= 0:
                self._previous_co2_mass = 0.01
            if self._previous_brine_mass <= 0:
                self._previous_brine_mass = 0.01

            self.norm_pressure = np.interp(
                self.pressure, (self.pars_bounds['pressure']), (0, 1))
            self.norm_co2_saturation = np.interp(
                self.co2_saturation, (self.pars_bounds['co2_saturation']), (0, 1))
            self.norm_logK_well_permeability = np.interp(
                self.logK_well_permeability,
                (self.pars_bounds['logK_well_permeability']), (0, 1))

            self.norm_previous_log10_co2_mass = np.interp(
                np.log10(self._previous_co2_mass),
                (self.pars_bounds['log10_co2_mass']), (0, 1))
            self.norm_previous_log10_brine_mass = np.interp(
                np.log10(self._previous_brine_mass),
                (self.pars_bounds['log10_brn_mass']), (0, 1))

            self.co2_array = np.reshape([self.norm_pressure,
                                         self.norm_co2_saturation,
                                         self.norm_logK_well_permeability,
                                         self.norm_previous_log10_co2_mass],
                                        (1, 1, 4))
            self.brn_array = np.reshape([self.norm_pressure,
                                         self.norm_co2_saturation,
                                         self.norm_logK_well_permeability,
                                         self.norm_previous_log10_brine_mass],
                                        (1, 1, 4))

            self.Outputs[0] = 10**(np.interp(
                kw_model('co2').predict(self.co2_array)[0, 0, 3],
                (0, 1), (self.pars_bounds['log10_co2_mass'])))
            self.Outputs[1] = 10**(np.interp(
                kw_model('brine').predict(self.brn_array)[0, 0, 3],
                (0, 1), (self.pars_bounds['log10_brn_mass'])))

            self.Outputs[2] = (self.Outputs[0] - self._previous_co2_mass)/self.dt/86400
            self.Outputs[3] = (self.Outputs[1] - self._previous_brine_mass)/self.dt/86400

if __name__ == "__main__":
    # Hydrology ROM
    simtime = 50
    dt = 365.25
    pressure = 28228195.8
    co2_saturation = 0.45741
    logK_well_permeability = -11.5926
    previous_co2_mass = 348953.5334519825
    previous_brine_mass = 2754.645894381841

    input_array1 = [pressure, co2_saturation, logK_well_permeability,
                    previous_co2_mass, simtime, dt]
    input_array2 = [pressure, co2_saturation, logK_well_permeability,
                    previous_brine_mass, simtime, dt]

    # Input parameters are not normalized
    # print(input_array1)
    # print(input_array2)
    # [28228195.8, 0.45741, -11.5926, 348953.5334519825, 50, 365.25]
    # [28228195.8, 0.45741, -11.5926, 2754.645894381841, 50, 365.25]

    # The answer should be
    # results:
    # co2_mass           349538.492290 kg
    # brine_mass           2811.847020 kg
    # co2_rate                0.000019 kg/s
    # brine_rate              0.000002 kg/s

    sol = Solution()
    sol.find(input_array1, input_array2)

    labels = ['CO2_mass', 'brine_mass', 'CO2_rate', 'brine_rate']
    units = ['kg', 'kg', 'kg/s', 'kg/s']

    results = sol.Outputs

    print('    results:')
    for i, v in enumerate(labels):
        print('   %-15s  %15.6f %s'%(v, results[i], units[i]))

    time_array = 365.25*np.arange(0.0, 50.0)

    pressures = [24569431.4, 24904877.8, 25240324.3, 25329010.4, 25417696.6,
                 25506382.8, 25573948.0, 25641513.3, 25709078.6, 25776643.8,
                 25844209.1, 25924208.6, 26004208.2, 26084207.7, 26164207.3,
                 26244206.8, 26339739.1, 26435271.3, 26530803.5, 26626335.8,
                 26721868.0, 26796140.6, 26870413.1, 26944685.6, 27018958.1,
                 27093230.6, 27152614.3, 27211997.9, 27271381.6, 27330765.3,
                 27390148.9, 27439616.0, 27489083.2, 27538550.3, 27588017.4,
                 27637484.5, 27679685.5, 27721886.4, 27764087.4, 27806288.4,
                 27848489.3, 27888139.9, 27927790.4, 27967440.9, 28007091.5,
                 28046742.0, 28083032.8, 28119323.5, 28155614.3, 28191905.1,
                 28228195.8]

    co2_saturations = [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
                       0.00465, 0.00930, 0.01395, 0.01860, 0.02325, 0.05142,
                       0.07959, 0.10776, 0.13593, 0.16410, 0.20317, 0.24224,
                       0.28131, 0.32038, 0.35946, 0.36743, 0.37541, 0.38339,
                       0.39137, 0.39935, 0.40331, 0.40727, 0.41123, 0.41518,
                       0.41914, 0.42182, 0.42450, 0.42718, 0.42986, 0.43254,
                       0.43456, 0.43658, 0.43860, 0.44062, 0.44264, 0.44426,
                       0.44587, 0.44748, 0.44910, 0.45071, 0.45205, 0.45339,
                       0.45473, 0.45607, 0.45741]

    print('\n'*2)
    header = '     Time    CO2_Mass(kg)     Brine_Mass(kg)    CO2_Rate(kg/s)  Brine_Rate(kg/s)'
    print(header)

    for i, simtime in enumerate(time_array):
        dt = simtime - time_array[i-1]

        input_array1 = [pressure, co2_saturation, logK_well_permeability,
                        previous_co2_mass, simtime, dt]
        input_array2 = [pressure, co2_saturation, logK_well_permeability,
                        previous_brine_mass, simtime, dt]

        sol.find(input_array1, input_array2)
        results = sol.Outputs
        previous_co2_mass = results[0]
        previous_brine_mass = results[1]

        print('  %5.1d  %16.2f   %16.2f   %14.8f  %14.8f'%(
            simtime/365.25, results[0], results[1], results[2], results[3]))

    # The answer should be

    #     Time    CO2_Mass(kg)     Brine_Mass(kg)    CO2_Rate(kg/s)  Brine_Rate(kg/s)
    #      0              0.00               0.00       0.00000000      0.00000000
    #      1              0.07              50.22       0.00000000      0.00000159
    #      2              1.09              73.46       0.00000003      0.00000074
    #      3              9.68              97.60       0.00000027      0.00000077
    #      4             54.62             125.66       0.00000142      0.00000089
    #      5            221.75             157.56       0.00000530      0.00000101
    #      6            698.51             193.15       0.00001511      0.00000113
    #      7           1803.03             232.27       0.00003500      0.00000124
    #      8           3970.78             274.71       0.00006869      0.00000134
    #      9           7692.33             320.23       0.00011793      0.00000144
    #     10          13421.45             368.60       0.00018154      0.00000153
    #     11          21484.96             419.59       0.00025552      0.00000162
    #     12          32022.25             472.96       0.00033391      0.00000169
    #     13          44966.54             528.46       0.00041018      0.00000176
    #     14          60064.91             585.89       0.00047844      0.00000182
    #     15          76924.75             645.02       0.00053426      0.00000187
    #     16          95071.02             705.66       0.00057502      0.00000192
    #     17         114002.21             767.61       0.00059989      0.00000196
    #     18         133234.04             830.70       0.00060942      0.00000200
    #     19         152335.86             894.75       0.00060530      0.00000203
    #     20         170946.20             959.61       0.00058973      0.00000206
    #     21         188780.73            1025.15       0.00056514      0.00000208
    #     22         205630.92            1091.22       0.00053395      0.00000209
    #     23         221359.40            1157.72       0.00049841      0.00000211
    #     24         235887.99            1224.52       0.00046038      0.00000212
    #     25         249190.00            1291.53       0.00042152      0.00000212
    #     26         261273.67            1358.65       0.00038291      0.00000213
    #     27         272179.60            1425.80       0.00034559      0.00000213
    #     28         281965.25            1492.90       0.00031009      0.00000213
    #     29         290704.30            1559.88       0.00027692      0.00000212
    #     30         298472.17            1626.69       0.00024615      0.00000212
    #     31         305353.29            1693.25       0.00021805      0.00000211
    #     32         311429.96            1759.52       0.00019256      0.00000210
    #     33         316781.47            1825.47       0.00016958      0.00000209
    #     34         321482.64            1891.03       0.00014897      0.00000208
    #     35         325604.51            1956.18       0.00013061      0.00000206
    #     36         329213.23            2020.89       0.00011435      0.00000205
    #     37         332364.60            2085.12       0.00009986      0.00000204
    #     38         335116.07            2148.84       0.00008719      0.00000202
    #     39         337513.96            2212.04       0.00007598      0.00000200
    #     40         339603.59            2274.70       0.00006622      0.00000199
    #     41         341420.67            2336.79       0.00005758      0.00000197
    #     42         343000.85            2398.30       0.00005007      0.00000195
    #     43         344374.82            2459.23       0.00004354      0.00000193
    #     44         345566.81            2519.55       0.00003777      0.00000191
    #     45         346602.90            2579.26       0.00003283      0.00000189
    #     46         347500.57            2638.35       0.00002845      0.00000187
    #     47         348277.68            2696.81       0.00002463      0.00000185
    #     48         348953.53            2754.65       0.00002142      0.00000183
    #     49         349538.49            2811.85       0.00001854      0.00000181
