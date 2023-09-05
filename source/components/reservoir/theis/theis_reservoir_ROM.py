"""
The module contains the solution class for the Theis reservoir ROM.

The solution is based on Theis, C.V., 1935. The relation between the lowering
of the piezometric surface and the rate and duration of discharge of a well
using groundwater storage, Am. Geophys. Union Trans., vol. 16, pp. 519-524.

Author: Diana Bacon and Alex Hanna
Date: 10/18/2022
"""
import os
import sys
import numpy as np
from scipy.special import exp1
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
np.set_printoptions(linewidth=np.inf)

class Parameters():
    """ Parameters class for Theis reservoir ROM."""
    def __init__(self):
        """ Create an instance of the Parameters class. """

        # Initial Pressure
        self.initialPressure = None

        # Reservoir Thickness
        self.reservoirThickness = None

        # Permeability
        self.reservoirPermeability = None

        # Porosity
        self.reservoirPorosity = None

        # Brine Density
        self.brineDensity = None

        # Brine Viscosity
        self.brineViscosity = None

        # CO2 Density
        self.CO2Density = None

        # Compressibility
        self.compressibility = None

        # Distance between injection and leaking well
        self.distance = None

        # Injection rate at the injection well (+ injection, - pumping)
        self.injRates = None

        # Times corresponding to each injection rate
        self.injTimes = None

class Solution():
    """ Solution class for simple reservoir ROM."""
    def __init__(self, parameters):
        """ Create an instance of the Solution class. """

        # Setup the parameters that define a solution
        self.parameters = Parameters()
        self.parameters = parameters

        # Reserve a space for the solution pieces
        self.pressure = None

        self.gravity = 9.80665  # acceleration due to gravity, m/s^2
        self.T = None  # Transmissivity (m^2/day)
        self.Q = None  # Scaled pumping rate (m^3/day)
        self.s = None  # Drawdown (m)
        self.S = None  # Storage Coefficient
        self.u = None  # dimensionless time
        self.W = None  # Well Function

    def dimensionless_time(self, r, t):
        """
        Calculate and return the dimensionless time parameter, u.
        """
        with np.errstate(divide='ignore'):
            return r**2 * self.S / 4 / self.T / t

    def drawdown(self, t, Q, r):
        """
        Calculate and return the drawdown in meters

        This version uses the Theis equation, s(t,Q,r) = Q * W(u) / (4.pi.T),
        where W(u) is the Well function for u = Sr^2 / (4Tt).
        t is the time (days),
        S is the Storage Coefficient,
        T is the Transmissivity (m2/day),
        r is the distance from the well (m), and
        Q is the pumping rate (m3/day).

        Theis, C.V., 1935. The relation between the lowering of the piezometric surface
        and the rate and duration of discharge of a well using groundwater storage,
        Am. Geophys. Union Trans., vol. 16, pp. 519-524.
        """
        self.u = self.dimensionless_time(r, t)
        self.W = exp1(self.u)
        return Q * self.W / (4*np.pi*self.T)

    def transient_drawdown(self, t, Q, r):
        """
        t is the time (days),
        Q is the pumping rate (m3/day)
        r is the distance from the well (m)

        requirements:
            t and Q must be 1-dimensional lists or arrays of the same length
            r can be a scalar, list or array of multiple dimensions

        Birsoy, Yuksel K., and W. K. Summers. 1980. “Determination of Aquifer Parameters
        from Step Tests and Intermittent Pumping Data.” Ground Water 18 (2): 137–46.
        https://doi.org/10.1111/j.1745-6584.1980.tb03382.x.
        """

        dd = self.drawdown(t, Q[0], r)
        nt = len(t)
        for jt in range(1, nt):
            dQ = Q[jt] - Q[jt-1]
            if dQ != 0:
                for it in range(jt+1, nt):
                    dd[..., it] += self.drawdown(t[it] - t[jt], dQ, r).squeeze()
        return dd

    def pressure_change(self, s):
        """Calculate change in pressure from drawdown (change in head)

        s = drawdown (m)
        density = fluid mass density (kg/m3)
        g = acceleration due to gravity (m/s2)
        dp = change in pressure (Pa)

        Freeze, R. Allan, and John A. Cherry. 1979. Groundwater.
        Englewood Cliffs, New Jersey: Prentice-Hall, Inc.

        """
        dp = s * self.parameters.brineDensity * self.gravity
        return dp

    def update_distance(self, value):
        """
        Update distance between injection well and location of interest with new
        value.

        Parameters
        ----------
        value : float
            Distance between injection well and location of interest.

        Returns
        -------
        None.
        """
        self.parameters.distance = value

    def update_inj_rates(self, new_inj_rates):
        """
        Update injection rates with values from new array.

        Parameters
        ----------
        new_inj_rates : array-like
            list or array of new injection rates

        Returns
        -------
        None.
        """
        self.parameters.injRates = new_inj_rates

    def update_inj_times(self, new_inj_times):
        """
        Update injection times with values from new array.

        Parameters
        ----------
        new_inj_times : array-like
            list or array of new injection times

        Returns
        -------
        None.
        """
        self.parameters.injTimes = new_inj_times

    def find(self):
        """ Find solution of the ROM corresponding to the provided parameters."""

        # Scale injection rate to account for CO2 density
        self.Q = np.asarray(self.parameters.injRates) * self.parameters.CO2Density \
            / self.parameters.brineDensity * 86400

        waterCompressibility = 4.4e-10
        # Common term for storativity and transmissivity
        common_term = self.parameters.reservoirThickness \
            * self.parameters.brineDensity * self.gravity

        # Storativity
        self.S = common_term * (self.parameters.compressibility \
            + self.parameters.reservoirPorosity * waterCompressibility)

        # Transmissivity
        self.T = common_term * self.parameters.reservoirPermeability \
            / self.parameters.brineViscosity
        # Convert to the right time units
        self.T = self.T * 86400

        # Drawdown (change in hydraulic head)
        t = np.asarray(self.parameters.injTimes)
        self.s = self.transient_drawdown(t, self.Q, self.parameters.distance)

        # Pressure
        dp = self.pressure_change(self.s)
        self.pressure = self.parameters.initialPressure + dp
        self.dp = dp


if __name__ == "__main__":

    parameters = Parameters()
    parameters.initialPressure = 1.0e6                # Pa
    parameters.reservoirThickness = 30                # m
    parameters.reservoirPermeability = 10**-10.69897  # m^2
    parameters.reservoirPorosity = .2
    parameters.brineDensity = 1000                    # kg/m^3
    parameters.brineViscosity = 2.535e-3              # Pa*s
    parameters.CO2Density = 800                       # kg/m^3
    parameters.compressibility = 2.46e-9              # 1/Pa
    parameters.distance = 55                          # m
    parameters.injRates = [-5.0e-3]                   # m^3/s
    parameters.injTimes = [1]

    sol = Solution(parameters)
    sol.find()

    print('Q = ', sol.Q[0]/86400, 'm^3/s')
    print('r = ', parameters.distance, 'm')
    print('T = ', sol.T/86400, 'm^2/s')
    print('S = ', sol.S)
    print('s = ', sol.s[0], 'm')
    print('pressure =', sol.pressure, 'Pa')

    # Expected output:
    # Q =  -0.004 m^3/s
    # r =  55 m
    # T =  0.00232110061488997 m^2/s
    # S =  0.0007496203260000001
    # s =  -0.7260346137646179 m

    # Compare results to Freeze & Cherry, Figure 8.5b
    print('Figure 8.5 (b) Calculated curve of ho-h versus t.')
    print('t(s) h0-h')
    times = [1.0e2, 3.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 1.0e5]  # time in seconds
    parameters.injTimes[:] = np.array(times)/86400  # time in days
    # Constant injection rate for all (7) time points
    parameters.injRates = 7*[-5.0e-3]  # m^3/s
    sol = Solution(parameters)
    sol.find()
    for t, d in zip(times, -sol.s):
        print(t, d)

    # Expected output:
    # t(s) h0-h
    # 100.0 0.003687007860452865
    # 300.0 0.04152344851511029
    # 1000.0 0.14570801924046678
    # 3000.0 0.27575534374439825
    # 10000.0 0.4332526844525569
    # 30000.0 0.5816987331871969
    # 100000.0 0.7460290617634939

    # Compare results to Freeze & Cherry, Table 8.1
    print('Table 8.1 Values of W(u) for Various Values of u')
    print('u', end=' ')
    for coefficient in range(1, 10):
        print(coefficient, end=' ')
    print()
    for exponent in range(0, 16):
        print(10**-exponent, end=' ')
        for coefficient in range(1, 10):
            u = coefficient * 10**-exponent
            t = parameters.distance**2 * sol.S / 4 / sol.T / u
            parameters.injTimes = [t]
            parameters.injRates = [-5.0e-3]
            sol = Solution(parameters)
            sol.find()
            if exponent == 0:
                if coefficient < 4:
                    print('{:.3f}'.format(sol.W[0]), end=' ')
                elif coefficient < 6:
                    print('{:.4f}'.format(sol.W[0]), end=' ')
                elif coefficient < 8:
                    print('{:.5f}'.format(sol.W[0]), end=' ')
                else:
                    print('{:.6f}'.format(sol.W[0]), end=' ')
            else:
                print('{:.2f}'.format(sol.W[0]), end=' ')
        print()

    # # Expected output:
    # # u 1 2 3 4 5 6 7 8 9
    # # 1 0.219 0.049 0.013 0.0038 0.0011 0.00036 0.00012 0.000038 0.000012
    # # 0.1 1.82 1.22 0.91 0.70 0.56 0.45 0.37 0.31 0.26
    # # 0.01 4.04 3.35 2.96 2.68 2.47 2.30 2.15 2.03 1.92
    # # 0.001 6.33 5.64 5.23 4.95 4.73 4.54 4.39 4.26 4.14
    # # 0.0001 8.63 7.94 7.53 7.25 7.02 6.84 6.69 6.55 6.44
    # # 1e-05 10.94 10.24 9.84 9.55 9.33 9.14 8.99 8.86 8.74
    # # 1e-06 13.24 12.55 12.14 11.85 11.63 11.45 11.29 11.16 11.04
    # # 1e-07 15.54 14.85 14.44 14.15 13.93 13.75 13.59 13.46 13.34
    # # 1e-08 17.84 17.15 16.74 16.46 16.23 16.05 15.90 15.76 15.65
    # # 1e-09 20.15 19.45 19.05 18.76 18.54 18.35 18.20 18.07 17.95
    # # 1e-10 22.45 21.76 21.35 21.06 20.84 20.66 20.50 20.37 20.25
    # # 1e-11 24.75 24.06 23.65 23.36 23.14 22.96 22.81 22.67 22.55
    # # 1e-12 27.05 26.36 25.96 25.67 25.44 25.26 25.11 24.97 24.86
    # # 1e-13 29.36 28.66 28.26 27.97 27.75 27.56 27.41 27.28 27.16
    # # 1e-14 31.66 30.97 30.56 30.27 30.05 29.87 29.71 29.58 29.46
    # # 1e-15 33.96 33.27 32.86 32.58 32.35 32.17 32.02 31.88 31.76
