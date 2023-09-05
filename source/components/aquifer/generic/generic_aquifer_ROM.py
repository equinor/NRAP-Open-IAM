"""
The module contains the solution class for the Generic aquifer ROM

Authors: Diana Bacon

Date: 08/20/2021
Last modified: 11/15/202021

Component model input definitions:

* **aqu_thick** [|m|] (25 to 250) - thickness of unit (default: 33.2);
  *linked to Stratigraphy*

* **top_depth** [|m|] (100 to 4100) - depth to the top of aquifer (default: 590.1);
  *linked to Stratigraphy*

* **por** [-] (0.02 to 0.25) - porosity of unit (default: 0.118)

* **log_permh** [|log10| |m^2|] (-14 to -10) - horizontal permeability
  (default: -13.39)

* **log_aniso** [|log10|] (0 to 3) - anisotropy ratio Kh/Kv (default: 0.3)

* **aquifer_salinity** [-] (0.0 to 0.015) - salt mass fraction in aquifer water
  (default: 0.005)

* **reservoir_salinity** [-] (0.015 to 0.05) - salt mass fraction in leak water
  (default: 0.03)

* **dissolved_salt_threshold** [-] (0.0 to 1.0) - threshold for salt mass fraction

* **dissolved_co2_threshold** [-] (0.0 to 1.0) - threshold for CO2 mass fraction

* **log_co2_mass** * [|log10| |kg|] (0.5 to 10.84) - cumulative CO2 mass leaked
  (default 6.34 |log10| (kg))

* **log_brine_mass** * [|log10| |kg|] (0.5 to 10.84) - cumulative brine mass leaked
  (default 6.34 |log10| (kg))

Gridded observations from the Generic Aquifer component are:

* **Dissolved_CO2_mass_fraction** [-] - mass fraction of |CO2| in aquifer
  pore water on a 100x10 radial grid surrounding the leaky well

* **Dissolved_salt_mass_fraction** [-] - mass fraction of salt in aquifer
  pore water on a 100x10 radial grid surrounding the leaky well

Observations from the Generic Aquifer component are:

* **Dissolved_salt_volume** [|m^3|] - volume of plume where relative change
  in salt mass fraction > dissolved_salt_threshold

* **Dissolved_salt_dr** [|m|] - radius of plume where relative change
  in salt mass fraction > dissolved_salt_threshold

* **Dissolved_salt_dz** [|m|] - height of plume where relative change
  in salt mass fraction > dissolved_salt_threshold

* **Dissolved_CO2_volume** [|m^3|] - volume of plume where dissolved
  |CO2| mass fraction > dissolved_co2_threshold

* **Dissolved_CO2_dr** [|m|] - radius of plume where dissolved
  |CO2| mass fraction > dissolved_co2_threshold

* **Dissolved_CO2_dz** [|m|] - height of plume where dissolved
  |CO2| mass fraction > dissolved_co2_threshold
"""

import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from model_input import generate_input, get_vertices, get_dimensions
from salt_predict import model as salt_model
from co2_predict import model as co2_model

class Solution():
    def __init__(self):
        """ Create an instance of the Solution class. """
        self.Outputs = np.zeros(6)
        self.GriddedOutputs = np.zeros((4, 100, 10))

    def find(self, inputArray):
        """ Find output of the ROM for the provided input parameters."""

        # grid centroids
        inp = generate_input(inputArray)
        r_coord = inp[0, :, 0, :, 0]
        z_coord = inp[0, :, 0, :, 1]

        time = inputArray[9]
        if time == 0.:
            self.Outputs = np.zeros(6)
            self.GriddedOutputs = np.zeros((4, 100, 10))
            self.GriddedOutputs[0, :, :] = r_coord
            self.GriddedOutputs[1, :, :] = -z_coord
            return

        self.GriddedOutputs[0, :, :] = r_coord
        self.GriddedOutputs[1, :, :] = -z_coord

        # gridded concentrations
        salt_mass_frac = salt_model(inputArray)[0, :, 0, :, 0]
        co2_mass_frac = co2_model(inputArray)[0, :, 0, :, 0]

        self.GriddedOutputs[2, :, :] = np.clip(salt_mass_frac, 0, None)
        self.GriddedOutputs[3, :, :] = np.clip(co2_mass_frac, 0, None)

        thick = inputArray[0]
        depth = inputArray[1]
        rv, zv = get_vertices(depth, thick, 10)
        dr, dz, vol = get_dimensions(rv, zv)

        salt_thresh = inputArray[10]
        co2_thresh = inputArray[11]

        # calculate size of impact
        def impact(dr, dz, vol, conc, thresh):
            impact_vol = np.ma.masked_where(conc < thresh, vol).filled(0).sum()
            impact_dr = np.ma.masked_where(conc < thresh, dr).filled(0).sum(axis=1).mean()
            impact_dz = np.ma.masked_where(conc < thresh, dz).filled(0).sum(axis=0).mean()
            return impact_dr, impact_dz, impact_vol

        salt_dr, salt_dz, salt_vol = impact(dr, dz, vol, salt_mass_frac, salt_thresh)
        co2_dr, co2_dz, co2_vol = impact(dr, dz, vol, co2_mass_frac, co2_thresh)

        self.Outputs[0] = salt_vol
        self.Outputs[1] = salt_dr
        # Limiting height of plume from above by aquifer thickness
        self.Outputs[2] = min(salt_dz, thick)
        if inputArray[5] == 0.0:  # if leaked CO2 mass is zero
            self.Outputs[3] = 0.0
            self.Outputs[4] = 0.0
            self.Outputs[5] = 0.0
        else:
            self.Outputs[3] = co2_vol
            self.Outputs[4] = co2_dr
            # Limiting height of plume from above by aquifer thickness
            self.Outputs[5] = min(co2_dz, thick)


if __name__ == "__main__":

    time = 1 # years
    dissolved_salt_threshold = 0.02
    dissolved_co2_threshold = 0.01

    # parameters from test runs for 2100-2400m depth interval
    params = {}

    #[run]=[thick, elevation, por, log_permh,
    #       log_aniso, aq.salinity, log co2_rate, log brine_rate, res. salinity]
    params['1006'] = [1.3822E+02, -2.2416E+03, 2.4010E-01, -1.0961E+01,
                      1.4663E+00, 2.4632E-03, -8.3395E+00, 1.2203E+00, 4.5073E-02]
    params['3263'] = [1.0260E+02, -2.3066E+03, 7.9001E-02, -1.2434E+01,
                      2.2254E+00, 5.7080E-03, 4.0938E-01, -5.0149E+00, 4.0617E-02]

    def convert_params(par, time, dissolved_salt_threshold, dissolved_co2_threshold):
        # convert list of input parameters for STOMP simulation into input array for ML model
        inp = []
        inp.append(par[0])  # thick
        inp.append(-par[1]) # elevation to top depth
        inp.append(par[2])  # por
        inp.append(par[3])  # log_permh
        inp.append(par[4])  # log_aniso
        inp.append(10**par[6] * time * 365.25 * 86400) # log_co2_rate to co2_mass_leaked in kg
        inp.append(10**par[7] * time * 365.25 * 86400) # log_brine_rate to brine_mass_leaked in kg
        inp.append(par[5]) # aquifer salinity
        inp.append(par[8]) # reservoir salinity
        inp.append(time)
        inp.append(dissolved_salt_threshold)
        inp.append(dissolved_co2_threshold)
        return np.array(inp)

    expected = {}
    expected['1006'] = [2269111.67, 6.86, 172.77, 0.00, 0.00, 0.00]
    expected['3263'] = [2706.93, 0.26, 8.21, 5373325.51, 12.31, 180.58]

    labels = ['Dissolved_salt_volume', 'Dissolved_salt_dr', 'Dissolved_salt_dz',
              'Dissolved_CO2_volume', 'Dissolved_CO2_dr', 'Dissolved_CO2_dz']

    sol = Solution()
    for run in ['1006','3263']:
        input_array = convert_params(
            params[run], time, dissolved_salt_threshold, dissolved_co2_threshold)
        # print(input_array)
        sol.find(input_array)

        results = sol.Outputs

        print('Run '+run+' Prediction         Result       Expected')
        for i, v in enumerate(labels):
            print("{:20s}  {:12.2f}  {:12.2f}".format(
                v, results[i], expected[run][i]))
