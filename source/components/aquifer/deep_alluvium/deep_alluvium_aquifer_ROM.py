"""
The module contains the solution class for the deep alluvium aquifer ROM
based on the Kimberlina Synthetic Model dataset

Authors: class Kayyum Mansoor

Date: 07/17/2018
Last modified: 07/27/2018

outputs:
TDSOutputs: affected volume, length, width and height of TDS plume with
change greater than 100 mg/L from baseline
PressureOutputs: affected volume, length, width and height of pressure plume
with change greater than 500 Pascals from baseline
pHOutputs: affected volume, length, width and height of plume with pH = 6.75


Component model input definitions:

* **simtime** * [years] (0 to 200) - Simulation Time (default 20.0 years)

* **brine_rate** * [kg/s] (0 to 0.017) - Brine rate (default 0.00030 kg/s)

* **brine_mass** * [|log10| (kg)] (2.337 to 6.939) - Cumulative brine mass
  (default 4.928 |log10| (kg))

* **co2_rate** * [kg/s] (0 to 0.385) - CO2 rate (default 0.045 kg/s)

* **co2_mass** * [|log10| (kg)] (0.001 to 9.210) - Cumulative CO2 mass
  (default 7.214 |log10| (kg))

* **logK_sand1** * [|log10| (|m^2|)] (-12.92 to -10.92) - Permeability in layer 1
  10-546 m (default -11.92 |log10| (|m^2|))

* **logK_sand2** * [|log10| (|m^2|)] (-12.72 to -10.72) - Permeability in layer 2
  546-1225 m (default -11.72 |log10| (|m^2|))

* **logK_sand3** * [|log10| (|m^2|)] (-12.70 to -10.70) - Permeability in layer 3
  1225-1411 m (default -11.70 |log10| (|m^2|))

* **logK_caprock** * [|log10| (|m^2|)] (-16.70 to -14.70) - Caprock Permeability 0-5 m
  (default -15.70 |log10| (|m^2|))

* **correlationLengthX** * [m] (200 to 2000 m) - Correlation length in X-direction
  (default 1098.99 m)

* **correlationLengthZ** * [m] (10 to 150 m) - Correlation length in Z-direction
  (default 79.81 m)

* **sandFraction** * [-] (0.70 to 0.90) - Sand volume fraction (default 0.800)

* **groundwater_gradient** * [-] (0.0010 to 0.0017) - Regional groundwater
  gradient (dh/dx=change in hydraulic head/distance) (default 0.0013)

* **leak_depth** * [-] (424.4 to 1341.5 m) - Depth of leakage interval (default 885.51 m)

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
* **Pressure_dx** [|m|] - length of plume above baseline pressure change
  in Pa (change in pressure > 500 Pa).
* **Pressure_dy** [|m|] - width of plume above baseline pressure change in Pa
  (change in pressure > 500 Pa).
* **Pressure_dz** [|m|] - height of plume above baseline pressure change in Pa
  (change in pressure > 500 Pa).

* **pH_volume** [|m^3|] - volume of plume below pH threshold (pH < 6.75).
* **pH_dx** [|m|] - length of plume below pH threshold (pH < 6.75).
* **pH_dy** [|m|] - width of plume below pH threshold (pH < 6.75).
* **pH_dz** [|m|] - height of plume below pH threshold (pH < 6.75).

"""
import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from deep_alluvium import (tds_log_dx, tds_log_dy, tds_log_dz, tds_log_vol,
                           pressure_log_dx, pressure_log_dy,
                           pressure_log_dz, pressure_log_vol,
                           ph_log_dx, ph_log_dy, ph_log_dz, ph_log_vol)


class Solution(object):
    """ Solution class for Deep Alluvium Aquifer ROM."""
    def __init__(self):
        """ Constructor method. """

        self.Outputs = None

        self.maxdx, self.maxdy, self.maxdz = 5000.0, 10000.0, 1411.0
        self.maxvol = self.maxdx*self.maxdy*self.maxdz

        self.mindx, self.mindy, self.mindz = 50, 50, 40
        self.minvol = self.mindx*self.mindy*self.mindz

    def find(self, inputArray):
        """ Find solution for the provided parameters. """

        self.Outputs = np.zeros(12)
        self.Outputs[0] = 10**([v for v in tds_log_vol.model([inputArray])][0])
        self.Outputs[1] = 10**([v for v in tds_log_dx.model([inputArray])][0])
        self.Outputs[2] = 10**([v for v in tds_log_dy.model([inputArray])][0])
        self.Outputs[3] = 10**([v for v in tds_log_dz.model([inputArray])][0])

        self.Outputs[4] = 10**([v for v in pressure_log_vol.model([inputArray])][0])
        self.Outputs[5] = 10**([v for v in pressure_log_dx.model([inputArray])][0])
        self.Outputs[6] = 10**([v for v in pressure_log_dy.model([inputArray])][0])
        self.Outputs[7] = 10**([v for v in pressure_log_dz.model([inputArray])][0])

        self.Outputs[8] = 10**([v for v in ph_log_vol.model([inputArray])][0])
        self.Outputs[9] = 10**([v for v in ph_log_dx.model([inputArray])][0])
        self.Outputs[10] = 10**([v for v in ph_log_dy.model([inputArray])][0])
        self.Outputs[11] = 10**([v for v in ph_log_dz.model([inputArray])][0])

        # Constrain outputs
        self.Outputs[0] = min(self.maxvol, self.Outputs[0])
        self.Outputs[4] = min(self.maxvol, self.Outputs[4])
        self.Outputs[8] = min(self.maxvol, self.Outputs[8])

        self.Outputs[1] = min(self.maxdx, self.Outputs[1])
        self.Outputs[5] = min(self.maxdx, self.Outputs[5])
        self.Outputs[9] = min(self.maxdx, self.Outputs[9])

        self.Outputs[2] = min(self.maxdy, self.Outputs[2])
        self.Outputs[6] = min(self.maxdy, self.Outputs[6])
        self.Outputs[10] = min(self.maxdy, self.Outputs[10])

        self.Outputs[3] = min(self.maxdz, self.Outputs[3])
        self.Outputs[7] = min(self.maxdz, self.Outputs[7])
        self.Outputs[11] = min(self.maxdz, self.Outputs[11])


if __name__ == "__main__":

    # Hydrology ROM
    time = 20
    brineRate = 3.21903E-05
    log_brineMass = 4.71081307
    co2Rate = 0.060985038
    log_co2Mass = 6.737803184
    logK_sand1 = -11.92098495
    logK_sand2 = -11.7198002
    logK_sand3 = -11.70137252
    logK_caprock = -15.69758676
    corrlationLengthX = 1098.994284
    corrlationLengthZ = 79.8062668
    sandFraction = 0.800121364
    groundwater_gradient = 0.001333374
    leak_depth = 885.5060281

    input_array = np.array([time, brineRate, log_brineMass, co2Rate, log_co2Mass,
                            logK_sand1, logK_sand2, logK_sand3, logK_caprock,
                            corrlationLengthX, corrlationLengthZ,
                            sandFraction, groundwater_gradient, leak_depth])

    sol = Solution()
    sol.find(input_array)

    labels = ['TDS_volume', 'TDS_dx', 'TDS_dy', 'TDS_dz',
              'Pressure_volume', 'Pressure_dx', 'Pressure_dy', 'Pressure_dz',
              'pH_volume', 'pH_dx', 'pH_dy', 'pH_dz']

    results = sol.Outputs

    print('results:')
    for i, v in enumerate(labels):
        print('%-15s  %14.2f'%(v, results[i]))

    # The answer should be
    # results:
    # TDS_volume           2436214.14
    # TDS_dx                   121.08
    # TDS_dy                   139.62
    # TDS_dz                   257.08
    # Pressure_volume      2945213.33
    # Pressure_dx              176.71
    # Pressure_dy              152.40
    # Pressure_dz              260.96
    # pH_volume            2656675.38
    # pH_dx                     95.07
    # pH_dy                     88.09
    # pH_dz                    249.01
