"""
The module contains the solution class for the open wellbore model.

Authors: Diana Bacon, Veronika Vasylkivska
Date: 01/09/2020
"""
import os
import logging

import numpy
from scipy.interpolate import RegularGridInterpolator


class Solution():
    """ Solution class for open wellbore ROM."""
    def __init__(self, header_file_directory):
        """ Create an instance of the Solution class."""
        # Read data from file
        file_path = os.path.join(header_file_directory, 'Open_Wellbore_2020.txt')
        with open(file_path) as f:
            self.top = [int(x) for x in next(f).split(' ')]
            self.depth = [float(x) for x in next(f).split(' ')]
            self.trans = [float(x) for x in next(f).split(' ')]
            self.aquif = [float(x) for x in next(f).split(' ')]
            self.dpres = [float(x) for x in next(f).split(' ')]
            self.sat = [float(x) for x in next(f).split(' ')]
            self.brine = [float(x) for x in next(f).split(' ')]

            ntop = len(self.top)
            ndepth = len(self.depth)
            ntrans = len(self.trans)
            naquif = len(self.aquif)
            npres = len(self.dpres)
            nsat = len(self.sat)
            nbrine = len(self.brine)

            co2Rate = numpy.zeros(shape=(
                ntop, ndepth, ntrans, naquif, npres, nsat, nbrine))
            brineRate = numpy.zeros(shape=(
                ntop, ndepth, ntrans, naquif, npres, nsat, nbrine))
            co2Surf = numpy.zeros(shape=(
                ndepth, ntrans, npres, nsat, nbrine))
            brineSurf = numpy.zeros(shape=(
                ndepth, ntrans, npres, nsat, nbrine))

            for d in range(ndepth):
                for r in range(ntrans):
                    for p in range(npres):
                        for s in range(nsat):
                            for b in range(nbrine):
                                co2Surf[d][r][p][s][b], brineSurf[d][r][p][s][b] = [
                                    float(x) for x in next(f).split(' ')]

            for t in range(ntop):
                for d in range(ndepth):
                    for r in range(ntrans):
                        for a in range(naquif):
                            for p in range(npres):
                                for s in range(nsat):
                                    for b in range(nbrine):
                                        co2Rate[t][d][r][a][p][s][b], brineRate[t][d][r][a][p][s][b] = [
                                            float(x) for x in next(f).split(' ')]

        # Create interpolating functions
        self.co2AtmInterpolator = RegularGridInterpolator(
            (self.depth, self.trans, self.dpres, self.sat, self.brine), co2Surf)
        self.brineAtmInterpolator = RegularGridInterpolator(
            (self.depth, self.trans, self.dpres, self.sat, self.brine), brineSurf)
        self.co2Interpolator = RegularGridInterpolator(
            (self.top, self.depth, self.trans, self.aquif,
             self.dpres, self.sat, self.brine),
            co2Rate)
        self.brineInterpolator = RegularGridInterpolator(
            (self.top, self.depth, self.trans, self.aquif,
             self.dpres, self.sat, self.brine),
            brineRate)

        # Placeholder for the future results
        self.CO2LeakageRates = None
        self.brineLeakageRates = None

    def find(self, input_array):
        """ Find solution of the ROM corresponding to the provided parameters."""
        self.CO2LeakageRates = numpy.zeros(2)
        self.brineLeakageRates = numpy.zeros(2)
        # Range checking
        inputPres = input_array[4]
        inputSat = input_array[5]
        if input_array[1] < self.depth[0]:
            input_array[1] = self.depth[0]
            logging.warning("Clipped to minimum depth of %d", self.depth[0])
        if input_array[1] > self.depth[-1]:
            input_array[1] = self.depth[-1]
            logging.warning("Clipped to maximum depth of %d", self.depth[-1])
        if input_array[2] < self.trans[0]:
            input_array[2] = self.trans[0]
            logging.warning("Clipped to minimum log reservoir transmissivity of %.2f",
                            self.trans[0])
        if input_array[2] > self.trans[-1]:
            input_array[2] = self.trans[-1]
            logging.warning("Clipped to maximum log reservoir transmissivity of %.2f",
                            self.trans[-1])
        if input_array[3] < self.aquif[0]:
            input_array[3] = self.aquif[0]
            logging.warning("Clipped to minimum log aquifer transmissivity of %.2f",
                            self.aquif[0])
        if input_array[3] > self.aquif[-1]:
            input_array[3] = self.aquif[-1]
            logging.warning("Clipped to maximum log aquifer transmissivity of %.2f",
                            self.aquif[-1])
        if input_array[4] < self.dpres[0]:
            input_array[4] = self.dpres[0]
            logging.warning("Scaling below minimum pressure differential of %d",
                            self.dpres[0])
        if input_array[4] > self.dpres[-1]:
            input_array[4] = self.dpres[-1]
            logging.warning("Clipped to maximum pressure differential of %d",
                            self.dpres[-1])
        if input_array[5] < self.sat[0]:
            input_array[5] = self.sat[0]
            logging.warning("Scaling below minimum saturation of %.2f",
                            self.sat[0])
        if input_array[5] > self.sat[-1]:
            input_array[5] = self.sat[-1]
            logging.warning("Clipped to maximum saturation of %.2f",
                            self.sat[-1])
        if input_array[6] < self.brine[0]:
            input_array[6] = self.brine[0]
            logging.warning("Clipped to minimum brine concentration of %.1f",
                            self.brine[0])
        if input_array[6] > self.brine[-1]:
            input_array[6] = self.brine[-1]
            logging.warning("Clipped to maximum brine concentration of %.1f",
                            self.brine[-1])

        # Save requested wellbore top depth
        requestedTop = input_array[0]

        if requestedTop == 0:
            surfaceArray = numpy.array([
                input_array[1], input_array[2], input_array[4],
                input_array[5], input_array[6]])
            # Return leakage rate to atmosphere
            self.CO2LeakageRates[0] = self.co2AtmInterpolator(surfaceArray)
            self.brineLeakageRates[0] = self.brineAtmInterpolator(surfaceArray)
            # Set leakage rate to aquifer to zero
            self.CO2LeakageRates[1] = 0.
            self.brineLeakageRates[1] = 0.

            # If CO2 saturation is less than minimum, then linearly interpolate to zero
            if inputSat < self.sat[0]:
                self.CO2LeakageRates[0] = self.CO2LeakageRates[0]*max(
                    inputSat, 0)/self.sat[0]

            # If pressure differential is less than minimum, then linearly interpolate to zero
            if inputPres < self.dpres[0]:
                self.CO2LeakageRates[0] = self.CO2LeakageRates[0]*max(
                    inputPres, 0)/self.dpres[0]
                self.brineLeakageRates[0] = self.brineLeakageRates[0]*max(
                    inputPres, 0)/self.dpres[0]
        else:
            # Set leakage rate to atmosphere to zero
            self.CO2LeakageRates[0] = 0.0
            self.brineLeakageRates[0] = 0.0

            # Leakage rate with wellbore top at 500 m depth
            input_array[0] = 500.
            CO2LeakageRate500 = self.co2Interpolator(input_array)
            brineLeakageRate500 = self.brineInterpolator(input_array)

            # Leakage rate with wellbore top at 100 m depth
            input_array[0] = 100.
            CO2LeakageRate100 = self.co2Interpolator(input_array)
            brineLeakageRate100 = self.brineInterpolator(input_array)

            # Interpolate or extrapolate to requested wellbore top depth
            self.CO2LeakageRates[1] = CO2LeakageRate100 + (
                CO2LeakageRate500-CO2LeakageRate100)/(500.-100.)*(requestedTop-100.)
            self.brineLeakageRates[1] = brineLeakageRate100 + (
                brineLeakageRate500-brineLeakageRate100)/(500.-100.)*(requestedTop-100.)

            # If CO2 saturation is less than minimum, then linearly interpolate to zero
            if inputSat < self.sat[0]:
                self.CO2LeakageRates[1] = self.CO2LeakageRates[1] * max(
                    inputSat, 0)/self.sat[0]

            # If pressure differential is less than minimum, then linearly interpolate to zero
            if inputPres < self.dpres[0]:
                self.CO2LeakageRates[1] = self.CO2LeakageRates[1] * max(
                    inputPres, 0)/self.dpres[0]
                self.brineLeakageRates[1] = self.brineLeakageRates[1] * max(
                    inputPres, 0)/self.dpres[0]

# The following if statement needed only if one wants to run
# the script as well as to use it as an importable module
# The code below is mainly used for testing purposes.
if __name__ == "__main__":
    test = 1
    if test == 1:
        # Atmospheric release
        top = 0
        depth = 1000
        trans = -11.27
        aquif = 0 # not used for atmospheric release
        brine = 0
        dpres = 0.10E+06
        sat = 0.01
        # The answer should be
        # WARNING:root:Clipped to minimum log reservoir transmissivity of -11.27
        # WARNING:root:Clipped to maximum log aquifer transmissivity of -8.40
        # CO2:  [0.0030171 0.       ]
        # Brine:  [0.13063954 0.        ]

    elif test == 2:
        # Aquifer release
        top = 600
        depth = 2000
        trans = -10
        aquif = -10
        brine = 0.1
        dpres = 1.05E+06
        sat = 0.45
        # The answer should be
        # CO2:  [0.         2.21264335]
        # Brine:  [0.         0.92816901]

    elif test == 3:
        # Deep reservoir
        top = 500
        depth = 4000
        trans = -11
        aquif = -11
        brine = 0.1
        dpres = 20E+06
        sat = 0.12
        # The answer should be
        # CO2:  [0.         8.41153151]
        # Brine:  [0.         4.79529928]

    elif test == 4:
        # Bad parameters to test logging
        top = 500
        depth = 5000
        trans = -12
        aquif = -10
        brine = 0.1
        dpres = 20E+06
        sat = 0.12
        # The answer should be
        # WARNING:root:Clipped maximum depth to 4000
        # WARNING:root:Clipped minimum log reservoir transmissivity to -11.27
        # CO2:  [0.         7.75155573]
        # Brine:  [0.         4.25500822]

    elif test == 5:
        # Low saturation scaling
        top = 700
        depth = 1500
        trans = -11.27
        aquif = -10.13
        brine = 0.0475
        dpres = 0.068E+06
        sat = 0
        # The answer should be
        # WARNING:root:Clipped to minimum log reservoir transmissivity of -11.27
        # WARNING:root:Scaling below minimum saturation of 0.01
        # CO2:  [0. 0.]
        # Brine:  [0.         0.00916077]

    sol = Solution('.')
    input_array = numpy.array([top, depth, trans, aquif, dpres, sat, brine])
    sol.find(input_array)

    print('CO2: ', sol.CO2LeakageRates)
    print('Brine: ', sol.brineLeakageRates)
