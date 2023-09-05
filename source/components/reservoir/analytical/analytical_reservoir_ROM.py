"""
The module contains the solution class for the analytical reservoir ROM.

The solution is based on the Celia et al., 2011 paper
"Field-scale application of a semi-analytical model
for estimation of CO2 and brine leakage along old wells".

Further reading can be found in : Baek et al., NRAP-Open-IAM Analytical
Reservoir Model - Development and Testing, PNNL-31418, 2021.

Author: Seunghwan Baek, PNNL
Date: 05/19/2021
Last update: 09/03/2021

"""
import os
import sys
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from simple import units

class Parameters():
    """ Parameters class for analytical reservoir ROM."""
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
        self.CO2Density = None

        # Viscosity
        self.brineViscosity = None
        self.CO2Viscosity = None

        # Residual saturation
        self.brineResidualSaturation = None

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

        # CO2 Plume Profile : Plume height + Reservoir (Top) Depth
        self.plumeProfile = None

        # Brine compressibilty
        self.brineCompressibility = None

class Solution():
    """ Solution class for analytical reservoir ROM."""
    def __init__(self, parameters):
        """ Create an instance of the Solution class. """

        # Setup the parameters that define a solution
        self.parameters = Parameters()
        self.parameters = parameters

        # Reserve a space for the solution pieces
        self.pressure = None
        self.saturation = None
        self.pressureTop = None
        self.g = 9.8  # acceleration due to gravity

        # Reserve space for additional attributes
        self.depth = None
        self.initialPressure = None
        self.initialTopPressure = None
        self.thicknessH = None
        self.CO2Volume = None
        self.interface = None
        self.brineSaturationAq = None
        self.brineMobilityAq = None
        self.CO2SaturationAq = None
        self.CO2MobilityAq = None
        self.plumeProfile = None
        self.effectiveMobility = None
        self.Lambda = None

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
        # Pressure at bottom of the reservoir
        self.initialPressure = self.parameters.datumPressure + (
            self.parameters.brineDensity*self.g*self.depth)

        # Thickness of reservoir
        self.thicknessH = self.parameters.reservoirThickness

        # Pressure at top of the reservoir
        self.initialTopPressure = self.parameters.datumPressure + (
            self.parameters.brineDensity*self.g*(self.depth-self.thicknessH))

        # Setup initial masses of CO2 in reservoir
        self.CO2Volume = self.parameters.prevCO2Volume

        # Setup initial interface height h at the leaking well
        self.interface = 0.0

        # Setup initial saturations of CO2 and brine in reservoir
        self.brineSaturationAq = 1.0
        self.CO2SaturationAq = 0.0

        # Setup initial CO2 profile in reservoir
        self.plumeProfile = 0.0

    def find(self):
        """ Find solution of the ROM corresponding to the provided parameters."""
        # Setup time and time step
        timeStep = self.parameters.timeStep*units.day()
        QInj = self.parameters.rate

        # Setup initial pressure, saturations, masses and interface height
        self.setup_initial_conditions()

        self.find_lambda()

        # Mobility ratio constant
        self.get_mobility_for_reservoir()

        # Update CO2 Mass
        self.CO2Volume = self.CO2Volume+QInj*timeStep # it is volume

        # Calculate plume height and saturation at the leaky well
        self.get_interface()

        self.saturation = self.interface/self.thicknessH
        self.saturation = min([self.saturation,
                               1-self.parameters.brineResidualSaturation])

        # Calculate global pressure (vertically averaged) at the leaky well
        self.get_pressure()

        # Calculate top pressure at the leaky well
        averagedDensity = (1-self.saturation)*self.parameters.brineDensity + \
            self.saturation*self.parameters.CO2Density
        self.pressureTop = self.pressure - 0.5*averagedDensity*self.g*self.thicknessH

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

    def get_interface(self):
        """Calculate height of the interface between brine and CO2. """
        H = self.parameters.reservoirThickness
        phi = self.parameters.reservoirPorosity
        Sres = self.parameters.brineResidualSaturation
        r = self.parameters.distance

        xLeak = 2*np.pi*H*phi*(1-Sres)*(r**2)/(self.CO2Volume)

        if xLeak < 2/self.Lambda:
            self.interface = (1-Sres)*H
        elif xLeak >= 2*self.Lambda:
            self.interface = 0.
        else:
            self.interface = 1/(self.Lambda-1)*(np.sqrt(2*self.Lambda/xLeak)-1)*(1-Sres)*H

    def get_pressure(self):
        """Calculate global pressure based on Eqs.4-5 in Celia et al., 2011 """
        phi = self.parameters.reservoirPorosity
        Sres = self.parameters.brineResidualSaturation
        r = self.parameters.distance
        resOuterRadius = self.parameters.resOuterRadius

        # Function of Inj CO2 Volume
        xLeak = 2*np.pi*self.thicknessH*phi*(1-Sres)*(r**2)/(self.CO2Volume)

        delta_Density = (self.parameters.brineDensity - self.parameters.CO2Density)

        # Function of Inj CO2 Volume
        Psi_inf = 4.5*np.pi*self.thicknessH*phi*self.parameters.reservoirPermeability*(1-Sres)/(
            self.parameters.brineViscosity*self.parameters.rate* \
                self.parameters.brineCompressibility)
        # Function of Inj CO2 Volume
        Psi_fin = 2*np.pi*self.thicknessH*phi*(1-Sres)*resOuterRadius**2/self.CO2Volume
        Psi_OuterBD = min([Psi_inf, Psi_fin])

        # Constant
        Gamma = 2*np.pi*delta_Density*self.g*self.parameters.reservoirPermeability* \
            self.thicknessH**2/(self.parameters.brineViscosity*self.parameters.rate)

        # Function of saturation
        F_h = -1*self.Lambda/(self.Lambda-1)*(self.saturation+np.log(
            (self.Lambda-1)*self.saturation + 1)/(self.Lambda-1))

        if xLeak >= Psi_OuterBD:
            delta_pressure = 0.
            CaseCheck = 1

        elif Psi_OuterBD > xLeak >= 2*self.Lambda:
            delta_pressure = -1/2/Gamma*np.log(xLeak/Psi_OuterBD) + F_h
            CaseCheck = 2

        elif 2*self.Lambda > xLeak >= 2/self.Lambda:
            delta_pressure = 1/Gamma - np.sqrt(xLeak)/(Gamma*np.sqrt(2*self.Lambda))- \
                1/2/Gamma*np.log(2*self.Lambda/Psi_OuterBD)+F_h
            CaseCheck = 3

        elif 2/self.Lambda > xLeak:
            delta_pressure = -1/2/self.Lambda/Gamma*np.log(xLeak*self.Lambda/2)+ \
                1/Gamma-1/Gamma*1/self.Lambda-1/2/Gamma*np.log(2*self.Lambda/Psi_OuterBD)+F_h
            CaseCheck = 4

        else:
            warn_msg = 'Value of xLeak does not satisfy any of the conditions.'
            print(warn_msg)
            CaseCheck = 99
            delta_pressure = None

        # Global pressure, Pa
        initialAvePressure = (self.initialPressure+self.initialTopPressure)/2
        self.pressure = initialAvePressure + max(
            [delta_pressure*delta_Density*self.g*self.thicknessH, 0])

        # Plume Evolution with respect to actual depth
        # comparable to the result from Celia's web simulator.
        ResTopDepth = self.depth - self.thicknessH
        self.plumeProfile = -1.*(ResTopDepth+self.saturation*self.thicknessH)

def user_msg():
    msg1 = ''.join([
        'This is analytical_reservoir_ROM for estimation of pressure ',
        'propagation and plume evolution in reservoirs caused ',
        'by CO2 injection. Developed by PNNL in 2021.'])
    msg2 = 'Check "analytical_reservoir_component.py".'
    print(msg1)
    print(msg2)

if __name__ == "__main__":
    user_msg()
