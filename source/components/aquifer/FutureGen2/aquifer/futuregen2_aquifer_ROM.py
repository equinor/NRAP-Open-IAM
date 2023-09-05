"""
The module contains the solution class for the aquifer ROM
based on the FutureGen 2.0 dataset for the following unit:
St. Peter Sandstone
  aqu_thick = 61.6 m
  bottom depth = -591.9 m
  porosity = 0.18
  horizontal permeability = 10^-11.92
  anisotropy = 10^.30

Authors: Diana Bacon

Date: 08/11/2018
Last modified: 10/16/2018

Component model input definitions:

* **aqu_thick** * [m] (30 to 90 m) - Thickness of unit (default 33.2 m)

* **depth** * [m] (100 to 700 m) - Depth to bottom of unit (default 590.1 m)

* **por** * [-] (0.02 to 0.2) - Porosity of unit (default 0.118)

* **log_permh** * [|log10| (|m^2|)] (-14 to -11) - Horizontal permeability
  (default -13.39 |log10| (|m^2|))

* **log_aniso** * [|log10|] (0 to 3) - Anisotropy ratio (default 0.30 |log10|)

* **rel_vol_frac_calcite** * [ ] (0 to 1) - Relative volume fraction of calcite
  in solid phase (default 0.01)

* **log_co2_mass** * [|log10| (kg)] (0.5 to 10.84) - Cumulative CO2 mass
  (default 6.34 |log10| (kg))

* **log_brine_mass** * [|log10| (kg)] (0.5 to 10.84) - Cumulative brine mass
  (default 6.34 |log10| (kg))

Observations from the FutureGen2 Aquifer component are:

* **Pressure_volume** [|m^3|] - volume of plume where relative change in pressure > 0.065%

* **Pressure_dx** [|m|] - length of plume where relative change in pressure > 0.065%

* **Pressure_dy** [|m|] - width of plume where relative change in pressure > 0.065%

* **Pressure_dz** [|m|] - height of plume where relative change in pressure > 0.065%

* **pH_volume** [|m^3|] - volume of plume where absolute change in pH > 0.2

* **pH_dx** [|m|] - length of plume where absolute change in pH > 0.2

* **pH_dy** [|m|] - width of plume where absolute change in pH  > 0.2

* **pH_dz** [|m|] - height of plume where absolute change in pH  > 0.2

* **TDS_volume** [|m^3|] - volume of plume where relative change in TDS > 10%

* **TDS_dx** [|m|] - length of plume where relative change in TDS > 10%

* **TDS_dy** [|m|] - width of plume where relative change in TDS > 10%

* **TDS_dz** [|m|] - height of plume where relative change in TDS > 10%

* **Dissolved_CO2_volume** [|m^3|] - volume of plume where relative change
in dissolved |CO2| concentration > 20%

* **Dissolved_CO2_dx** [|m|] - length of plume where relative change
in dissolved |CO2| concentration > 20%

* **Dissolved_CO2_dy** [|m|] - width of plume where relative change
in dissolved |CO2| concentration > 20%

* **Dissolved_CO2_dz** [|m|] - height of plume where relative change
in dissolved |CO2| concentration > 20%

"""
import os
import sys
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from aquifer import (ph_log_dx, ph_log_dy, ph_log_dz, ph_log_vol,
                     tds_log_dx, tds_log_dy, tds_log_dz, tds_log_vol,
                     pres_log_dx, pres_log_dy, pres_log_dz, pres_log_vol,
                     dis_log_dx, dis_log_dy, dis_log_dz, dis_log_vol)


class Solution(object):
    """ NRAP-Open-IAM Solution class for FutureGen2 aquifer ROM. """
    def __init__(self):
        """ Create an instance of the Solution class. """
        self.Outputs = None

        # Maximum domain size
        self.maxdx, self.maxdy, self.maxdz = 80467.20*2, 80467.20*2, 90.
        self.maxvol = self.maxdx*self.maxdy*self.maxdz

        # Minimum grid cell size
        self.mindx, self.mindy, self.mindz = 0., 0., 0.
        self.minvol = self.mindx*self.mindy*self.mindz

    def find(self, inputArray):
        """ Find output of the ROM for the provided input parameters."""
        inputArray[1] = -inputArray[1]   # Change depth to negative value

        self.Outputs = np.zeros(16)
        self.Outputs[0] = np.expm1(next(tds_log_vol.model([inputArray])))
        self.Outputs[1] = np.expm1(next(tds_log_dx.model([inputArray])))
        self.Outputs[2] = np.expm1(next(tds_log_dy.model([inputArray])))
        self.Outputs[3] = np.expm1(next(tds_log_dz.model([inputArray])))

        self.Outputs[4] = np.expm1(next(ph_log_vol.model([inputArray])))
        self.Outputs[5] = np.expm1(next(ph_log_dx.model([inputArray])))
        self.Outputs[6] = np.expm1(next(ph_log_dy.model([inputArray])))
        self.Outputs[7] = np.expm1(next(ph_log_dz.model([inputArray])))

        self.Outputs[8] = np.expm1(next(pres_log_vol.model([inputArray])))
        self.Outputs[9] = np.expm1(next(pres_log_dx.model([inputArray])))
        self.Outputs[10] = np.expm1(next(pres_log_dy.model([inputArray])))
        self.Outputs[11] = np.expm1(next(pres_log_dz.model([inputArray])))

        self.Outputs[12] = np.expm1(next(dis_log_vol.model([inputArray])))
        self.Outputs[13] = np.expm1(next(dis_log_dx.model([inputArray])))
        self.Outputs[14] = np.expm1(next(dis_log_dy.model([inputArray])))
        self.Outputs[15] = np.expm1(next(dis_log_dz.model([inputArray])))

        # Constrain outputs
        for i in range(16):
            self.Outputs[i] = max(0, self.Outputs[i])

        # If no CO2 or brine has accumulated then outputs should be zero
        if (inputArray[9] == 0 and inputArray[10] == 0):
            for i in range(16):
                try:
                    self.Outputs[i] = 0
                except:
                    self.Outputs[i] = [0 for item in self.Outputs[i]]

        # If no CO2 has accumulated then pH and dissolved CO2 outputs should be zero
        if inputArray[9] == 0:
            for i in [4, 5, 6, 7, 12, 13, 14, 15]:
                try:
                    self.Outputs[i] = 0
                except:
                    self.Outputs[i] = [0 for item in self.Outputs[i]]

        # dz should not be greater than aquifer thickness

        for i in [3, 7, 11, 15]:
            self.Outputs[i] = min(inputArray[0], self.Outputs[i])

if __name__ == "__main__":

    # St. Peter
    aqu_thick = 61.6
    depth = 591.9
    por = 0.18
    log_permh = -11.92
    log_aniso = 0.3
    rel_vol_frac_calcite = 0.01

    co2_rate = 1.0e-2
    brine_rate = 1.0e-3
    time = 1.
    co2_mass = co2_rate*time*86400*365.25
    brine_mass = brine_rate*time*86400*365.25
    log_co2_rate = np.log1p(co2_rate)
    log_brine_rate = np.log1p(brine_rate)
    log_co2_mass = np.log1p(co2_mass)
    log_brine_mass = np.log1p(brine_mass)

    inputArray = np.array([aqu_thick, depth, por, log_permh, log_aniso,
                           rel_vol_frac_calcite, log_co2_rate, log_brine_rate,
                           time, log_co2_mass, log_brine_mass])

    sol = Solution()
    sol.find(inputArray)

    labels = ['TDS_volume', 'TDS_dx', 'TDS_dy', 'TDS_dz',
              'pH_volume', 'pH_dx', 'pH_dy', 'pH_dz',
              'Pressure_volume', 'Pressure_dx', 'Pressure_dy', 'Pressure_dz',
              'Dissolved_CO2_volume', 'Dissolved_CO2_dx',
              'Dissolved_CO2_dy', 'Dissolved_CO2_dz']

    results = sol.Outputs

    print('Results:')
    for i, v in enumerate(labels):
        print("{:20s}  {:12.2f}".format(v, results[i]))

    # The output should be
    # Results:
    # TDS_volume                  221.65
    # TDS_dx                        7.01
    # TDS_dy                        7.01
    # TDS_dz                       10.70
    # pH_volume                283401.12
    # pH_dx                       100.29
    # pH_dy                       100.29
    # pH_dz                        23.99
    # Pressure_volume              24.97
    # Pressure_dx                   2.36
    # Pressure_dy                   2.36
    # Pressure_dz                   1.76
    # Dissolved_CO2_volume     118108.91
    # Dissolved_CO2_dx             59.56
    # Dissolved_CO2_dy             59.56
    # Dissolved_CO2_dz             34.14
