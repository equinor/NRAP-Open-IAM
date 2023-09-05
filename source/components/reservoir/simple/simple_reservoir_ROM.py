"""
The module contains the solution class for the simple reservoir ROM.

The solution is based on the Celia et al., 2011 paper
"Field-scale application of a semi-analytical model
for estimation of CO2 and brine leakage along old wells".

Author: Veronika S. Vasylkivska
Date: 06/25/2015
"""
import os
import sys
import numpy as np
import scipy.special as scm
from scipy.interpolate import interp1d
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from simple import units

class Parameters():
    """ Parameters class for simple reservoir ROM."""
    def __init__(self):
        """ Create an instance of the Parameters class. """
        self.numberOfShaleLayers = None

        self.shaleThickness = None
        self.aquiferThickness = None
        self.reservoirThickness = None

        # Permeability
        self.reservoirPermeability = None

        # Porosity
        self.reservoirPorosity = None

        # Land surface pressure
        self.datumPressure = None

        # Density
        self.brineDensity = None
        self.brineDensityForm = "Constant"     # constant or functional
        self.CO2Density = None
        self.CO2DensityForm = "Constant"       # constant or functional

        # Viscosity
        self.brineViscosity = None
        self.brineViscosityForm = "Constant"     # constant or functional
        self.CO2Viscosity = None
        self.CO2ViscosityForm = "Constant"       # constant or functional

        # Residual saturation and compressibility
        self.brineResidualSaturation = None
        self.compressibility = None

        # Input time point
        self.timePoint = None

        # Amount of CO2 present in reservoir before simulation
        self.prevCO2Volume = None

        # Time step in days
        self.timeStep = None

        # Distance between injection and leaking well
        self.distance = None

        # Rate at the injection well
        self.rate = None


class Solution():
    """ Solution class for simple reservoir ROM."""
    def __init__(self, parameters):
        """ Create an instance of the Solution class. """

        # Setup the parameters that define a solution
        self.parameters = Parameters()
        self.parameters = parameters

        # Reserve a space for the solution pieces
        self.pressure = None
        self.saturation = None
        self.g = 9.8  # acceleration due to gravity

    def get_depth_for_pressure(self):
        """ Calculate depth of all layers for the pressure calculations. """
        # The deepest shale has the first thickness in the array
        shaleThickness = self.parameters.shaleThickness
        aquiferThickness = self.parameters.aquiferThickness
        reservoirThickness = self.parameters.reservoirThickness

        # Get the sum of shale layers thicknesses
        shaleThicknessSum = sum(shaleThickness)

        # Number of aquifers includes reservoir so
        # number of aquifers = (number of shales - 1) + 1
        aquiferNum = self.parameters.numberOfShaleLayers

        # Get the sum of aquifer thicknesses
        if aquiferNum == 2:
            aquiferThicknessSum = aquiferThickness
        else:
            aquiferThicknessSum = sum(aquiferThickness)

        # Get the depth of the bottom of the reservoir
        self.depth = shaleThicknessSum + aquiferThicknessSum + reservoirThickness

    def setup_initial_conditions(self):
        """ Setup initial conditions at the abandoned well"""

        # Determine depths of the aquifers
        self.get_depth_for_pressure()

        # Setup initial hydrostatic pressure
        # datumPressure is atmospheric pressure or some reference pressure
        # at the top of the upper shale layer
        self.initialPressure = self.parameters.datumPressure + (
            self.parameters.brineDensity*self.g*self.depth)

        self.nSL = self.parameters.numberOfShaleLayers

        # Thickness of reservoir
        self.thicknessH = self.parameters.reservoirThickness

        # Pressure at the top of the reservoir
        self.initialTopPressure = self.parameters.datumPressure + (
            self.parameters.brineDensity*self.g*(self.depth-self.thicknessH))

        # Setup initial masses of CO2 in reservoir
        self.CO2Volume = self.parameters.prevCO2Volume

        # Setup initial interface height h at the leaking well
        self.interface = 0.0

        # Setup initial saturations of CO2 and brine in reservoir
        self.brineSaturationAq = 1.0
        self.CO2SaturationAq = 0.0

    def find(self):
        """ Find solution of the ROM corresponding to the provided parameters."""
        # Setup initial pressure, saturations, masses and interface height
        self.setup_initial_conditions()

        self.find_lambda()
        self.get_mobility_for_reservoir()

        # Setup time and time step
        timeStep = self.parameters.timeStep*units.day()
        t = self.parameters.timePoint*units.day()

        # Parameters for injection well function
        compressibility = self.parameters.compressibility
        r0 = self.parameters.distance

        # u0 - argument of injection well function
        u0Scalar = r0**2*compressibility*self.parameters.reservoirPorosity/(
            4*self.parameters.reservoirPermeability)
        QInj = self.parameters.rate
        tScalar = 1
        if t > 0.0:
            uInj = u0Scalar/(self.effectiveMobility*tScalar*t)
            GInj = well_function(uInj)
        else:
            uInj = 0.0
            GInj = 0.0

        fourPi = 4*np.pi
        permxThickness = self.thicknessH*self.parameters.reservoirPermeability

        # Coefficient next to the well function in the pressure equation
        denomS = fourPi*permxThickness*self.effectiveMobility
        pressure = self.initialPressure+QInj*GInj/denomS

        self.pressure = pressure-(
            self.thicknessH-self.interface)*self.parameters.brineDensity*self.g-(
                self.interface*self.parameters.CO2Density*self.g)

        # Update CO2 volume
        self.CO2Volume = self.CO2Volume+QInj*timeStep

        self.get_interface()

        self.saturation = self.interface/self.thicknessH

    def get_mobility_for_reservoir(self):
        """ Calculate mobility parameter for brine and CO2."""
        CO2Permeability = min(1-self.parameters.brineResidualSaturation,
                              self.interface/self.thicknessH)
        brinePermeability = 1/(1-self.parameters.brineResidualSaturation)*(
            1-self.parameters.brineResidualSaturation-CO2Permeability)

        self.CO2MobilityAq = CO2Permeability/self.parameters.CO2Viscosity
        self.brineMobilityAq = brinePermeability/self.parameters.brineViscosity
        self.effectiveMobility = (self.interface*self.CO2MobilityAq+(
            self.thicknessH-self.interface)*self.brineMobilityAq)/(
                self.thicknessH)

        # For use in well functions
        if self.interface == 0:
            self.effectiveMobility = 1/self.parameters.brineViscosity

    def find_lambda(self):
        """ Calculate lambda constant. """
        # Lambda is a constant in the formula for the outer egde of the plume
        # It is constant for each aquifer in the system unless
        # brine and CO2 have different relative permeabilities, residual
        # saturations or viscosities
        # For simple calculations, the mobility ratio Lambda is calculated at
        # the endpoint relative permeability values for the two phases.

        # Maximum value for CO2 relative permeability
        CO2Permeability = 1 - self.parameters.brineResidualSaturation
        # Maximum value for brine relative permeability
        brinePermeability = 1.

        brineViscosity = self.parameters.brineViscosity
        CO2Viscosity = self.parameters.CO2Viscosity

        self.Lambda = CO2Permeability*brineViscosity/(
            brinePermeability*CO2Viscosity)

    def update_lambda(self):
        """ Update lambda constant. """
        Sres = self.parameters.brineResidualSaturation
        CO2Permeability = min((1-Sres), self.interface/self.thicknessH)
        brinePermeability = 1-CO2Permeability

        brineViscosity = self.parameters.brineViscosity
        CO2Viscosity = self.parameters.CO2Viscosity

        self.Lambda = CO2Permeability*brineViscosity/(
            brinePermeability*CO2Viscosity)

    def get_interface(self):
        """Calculate height of the interface between brine and CO2. """
        H = self.parameters.reservoirThickness
        phi = self.parameters.reservoirPorosity
        Sres = self.parameters.brineResidualSaturation
        r = self.parameters.distance
        xLeak = 2*np.pi*H*phi*(1-Sres)*(r**2)/self.CO2Volume

        if xLeak < 2/self.Lambda:
            self.interface = (1-Sres)*H
        elif xLeak >= 2*self.Lambda:
            self.interface = 0.
        else:
            self.interface = 1/(self.Lambda-1)*(
                np.sqrt(2*self.Lambda/xLeak)-1)*(1-Sres)*H

def well_function(x):
    """
    Return the approximation of the well function with even number
    of terms in expansion. Expansions with an odd number of terms
    diverge to plus infinity without crossing zero.
    """
    u = x
    if u <= 1.0:
        W = (-0.577216-np.log(u)+u-u**2/(2*scm.factorial(2))+
             u**3/(3*scm.factorial(3))-u**4/(4*scm.factorial(4))+
             u**5/(5*scm.factorial(5))-u**6/(6*scm.factorial(6))+
             u**7/(7*scm.factorial(7))-u**8/(8*scm.factorial(8)))
    elif u <= 9:
        uu = np.linspace(1.0, 9.0, num=9)
        Wu = np.array([0.219, 0.049, 0.013, 0.0038, 0.0011,
                       0.00036, 0.00012, 0.000038, 0.000012])
        Wfun = interp1d(uu, Wu, kind='linear')
        W = Wfun(u)
    else:
        W = 0.000001
    return W

if __name__ == "__main__":
    x = 5.0
    W = well_function(x)
    print(W)
