"""
The module contains the solution class for the multisegmented_well_model.

The solution is based on the Celia et al., 2011 paper
"Field-scale application of a semi-analytical model
for estimation of CO2 and brine leakage along old wells".

Author: Veronika S. Vasylkivska
Date: 06/25/2015
"""
import logging
import os
import sys
import numpy as np
import scipy.special as scm
from scipy.interpolate import interp1d
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from multisegmented import units


class Parameters():
    """ Parameters class for multisegmented wellbore ROM."""
    def __init__(self):
        """ Create an instance of the Parameters class. """
        self.numberOfShaleLayers = None

        self.shaleThickness = None
        self.aquiferThickness = None
        self.reservoirThickness = None

        # Permeability
        self.shalePermeability = None
        self.aquiferPermeability = None

        # Land surface pressure
        self.datumPressure = None

        # Input pressure, CO2 saturation and time points
        self.pressure = None
        self.CO2saturation = None
        self.timePoint = None

        # Amount of CO2 present in aquifers before simulation
        self.prevCO2Volume = None

        # Time step in days
        self.timeStep = None

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

        # Well radius
        self.wellRadius = None
        self.flowArea = None


class Solution():
    """ Solution class for multisegmented wellbore ROM."""
    def __init__(self, parameters):
        """ Create an instance of the Solution class."""

        # Setup the parameters that define a solution
        self.parameters = Parameters()
        self.parameters = parameters

        # Reserve a space for the solution pieces
        self.CO2LeakageRates = None
        self.brineLeakageRates = None
        self.g = 9.8  # acceleration due to gravity

        # Reserve space for instance attributes
        self.nSL = None
        self.reservoirBottom = None
        self.depth = None

        self.interface = None
        self.thicknessH = None
        self.initialPressure = None

        self.CO2Volume = None
        self.CO2SaturationAq = None

        self.CO2MobilityAq = None
        self.brineMobilityAq = None

        self.CO2MobilityShale = None
        self.brineMobilityShale = None
        self.effectiveMobility = None

        self.CO2Q = None
        self.brineQ = None
        self.Lambda = None

    def get_depth_for_pressure(self):
        """ Calculate depth of all layers for the pressure calculations."""
        # The deepest shale has the first thickness in the array
        shaleThickness = self.parameters.shaleThickness
        aquiferThickness = self.parameters.aquiferThickness
        reservoirThickness = self.parameters.reservoirThickness

        # Get the sum of shale layers thicknesses
        shaleThicknessSum = sum(shaleThickness)

        # Number of aquifer includes reservoir so
        # number of aquifers = (number of shales - 1) + 1 = number of shales
        aquiferNum = self.parameters.numberOfShaleLayers

        # Get the sum of aquifer (not including reservoir) thicknesses
        if aquiferNum == 2:
            aquiferThicknessSum = aquiferThickness
        else:
            aquiferThicknessSum = sum(aquiferThickness)

        # Get the depth of the bottom of the reservoir
        self.reservoirBottom = (shaleThicknessSum + aquiferThicknessSum +
                                reservoirThickness)

        # We need to know the depths to the bottom of each aquifer
        # Initialize depth array for depths of the bottoms of the aquifers
        self.depth = np.zeros((aquiferNum,))

        self.depth[0] = self.reservoirBottom
        self.depth[1] = self.depth[0] - (
            shaleThickness[0] + reservoirThickness)

        if aquiferNum > 2:
            for ind in range(2, aquiferNum):
                self.depth[ind] = self.depth[ind-1] - (
                    shaleThickness[ind-1] + aquiferThickness[ind-2])

    def setup_initial_conditions(self):
        """ Setup initial conditions at the abandoned well."""
        # Determine depths of the aquifers
        self.get_depth_for_pressure()

        # Setup initial hydrostatic pressure
        # datumPressure is atmospheric pressure or some reference pressure
        # at the top of the upper shale layer
        self.initialPressure = self.parameters.datumPressure + (
            self.parameters.brineDensity*self.g*self.depth)

        pressureTopReservoir = self.initialPressure[0] - (
            self.parameters.brineDensity*self.g*self.parameters.reservoirThickness)

        if pressureTopReservoir > self.parameters.pressure:
            warn_msg = ''.join([
                'Hydrostatic pressure at the top of reservoir exceeds ',
                'supplied pressure.\n This could explain negative leakage rates.'])
            logging.warning(warn_msg)

        self.nSL = self.parameters.numberOfShaleLayers

        # Setup initial masses of CO2 in each aquifer and reservoir
        self.CO2Volume = self.parameters.prevCO2Volume

        # Setup array for the thickness of aquifers and reservoir
        self.thicknessH = np.ones(self.nSL)
        self.thicknessH[0] = self.parameters.reservoirThickness
        self.thicknessH[1:self.nSL] = self.parameters.aquiferThickness

        # Setup initial interface height h at the leaking well
        self.interface = np.zeros(self.nSL)
        for j in range(self.nSL):
            self.interface[j] = self.parameters.CO2saturation[j]*self.thicknessH[j]

        # Setup saturation in aquifers
        self.CO2SaturationAq = np.zeros(self.nSL-1)

    def get_mobility_for_aquifers(self):
        """ Calculate mobility parameter for brine and CO2 in the aquifers."""

        CO2RelPermAq = np.zeros(self.nSL)
        brineRelPermAq = np.ones(self.nSL)

        Sres = self.parameters.aquBrineResSaturation

        for ind in range(self.nSL):
            CO2RelPermAq[ind] = np.minimum((
                1-Sres[ind]), self.interface[ind]/self.thicknessH[ind])
            brineRelPermAq[ind] = 1/(1-Sres[ind])*(
                1-Sres[ind]-CO2RelPermAq[ind])

        self.CO2MobilityAq = CO2RelPermAq/self.parameters.CO2Viscosity
        self.brineMobilityAq = brineRelPermAq/self.parameters.brineViscosity

        self.effectiveMobility = (self.interface*self.CO2MobilityAq+(
            self.thicknessH-self.interface)*self.brineMobilityAq)/self.thicknessH

        # For use in well functions
        for ind in range(0, self.nSL):
            if self.interface[ind] == 0:
                self.effectiveMobility[ind] = 1/self.parameters.brineViscosity

    def find(self):
        """ Find solution of the ROM corresponding to the provided parameters."""
        timeStep = self.parameters.timeStep*units.day()

        # Setup initial pressure, saturations, masses and interface height
        self.setup_initial_conditions()

        # Setup lambda and mobility
        self.find_lambda()

        self.get_mobility_for_aquifers()

        self.get_mobility_for_shales()

        # Initialize matrix for system of equations. AA*Pressure(Aquifers) = BB
        # LUT coupled: the pressure at the bottom of reservoir (or aquifer 1)
        # is known so it's not a variable
        StartingIdx = 1

        # tridiagonal (N-1)*(N-1) (N: # of aquifers including the bottom reservoir =  self.nSL)
        AA = np.zeros((self.nSL-StartingIdx, self.nSL-StartingIdx))
        BB = np.zeros(self.nSL-StartingIdx) # vector (N-1)*1

        # Calculate coefficients of Darcy equation for each aquifer for the matrix build-up
        # Refer to the document for the detailed information
        CC = self.CC_shale()
        GG = self.GG_shale()
        FF = self.FF_aquifer()
        WV = self.WV_aquifer()

        for ind in range(StartingIdx, self.nSL):

            F_below = FF[ind-1]
            if ind == StartingIdx:
                F_below = 0     # nothing below 1st shale layer

            F_zero = FF[ind]

            if ind < self.nSL-1:
                F_above = FF[ind+1]
            else:
                F_above = 0 # nothing above top shale layer

            C_below = CC[ind-1]
            C_zero = CC[ind]

            G_below = GG[ind-1]
            G_zero = GG[ind]

            # aquifer2 (bottom). Do not calculate the bottom reservoir due to LUT coupled case
            if ind == StartingIdx:
                # Top reservoir, not averaged, pressure from LUT
                ReservoirTopPressure = self.parameters.pressure

                AA[ind-StartingIdx, ind-StartingIdx] = 1+WV[ind]*C_below-WV[ind]*C_zero
                AA[ind-StartingIdx, ind+1-StartingIdx] = WV[ind]*C_zero

                BB[ind-StartingIdx] = WV[ind]*(-C_below*F_zero-G_below+C_zero*F_zero+\
                    C_zero*F_above-G_zero-C_below*ReservoirTopPressure) # + delta_pressure[ind-1]

            elif ind == (self.nSL-1): # aquifer N (top)
                DatumPressure = self.parameters.datumPressure

                AA[ind-StartingIdx, ind-1-StartingIdx] = -WV[ind]*C_below
                AA[ind-StartingIdx, ind-StartingIdx] = 1+WV[ind]*C_below-WV[ind-1]*C_zero

                BB[ind-StartingIdx] = WV[ind-1]*(-C_below*F_below-C_below*F_zero-\
                    G_below+C_zero*F_zero-G_zero-C_zero*DatumPressure) #+ delta_pressure[ind-1]

            else: # aquifer3 to aquifer N-1 (intermediate)

                AA[ind-StartingIdx, ind-1-StartingIdx] = -WV[ind]*C_below
                AA[ind-StartingIdx, ind-StartingIdx] = 1+WV[ind]*C_below-WV[ind]*C_zero
                AA[ind-StartingIdx, ind+1-StartingIdx] = WV[ind]*C_zero

                BB[ind-StartingIdx] = WV[ind]*(-C_below*F_below-C_below*F_zero-\
                    G_below+C_zero*F_zero+C_zero*F_above-G_zero) #+ delta_pressure[ind-1]

        # Calculate coupled pressure of each aquifer as solution of linear system
        # Vertically averaged delta pressure (compared to the initial values) at each aquifer
        deltaPave = np.linalg.solve(AA, BB)

        if StartingIdx == 1: # LUT coupled
            # Pressure in aquifers
            Pave = self.initialPressure[StartingIdx:] + deltaPave - FF[StartingIdx:]

            Pbot = Pave + FF[StartingIdx:]
            # Inclusion of the bottom reservoir pressure
            Pbot = np.append(self.initialPressure[0], Pbot)

            Ptop = Pave - FF[StartingIdx:]
            # Inclusion of the top reservoir pressure from LUT
            Ptop = np.append(ReservoirTopPressure, Ptop)

        # Find change in pressure to calculate flow rate
        deltaP = np.zeros(self.nSL)
        deltaP[0:self.nSL-1] = Ptop[0:self.nSL-1] - Pbot[1:self.nSL]
        deltaP[self.nSL-1] = Ptop[self.nSL-1] - self.parameters.datumPressure

        # Coefficients in Darcy equation
        DarcyCoefCO2 = self.parameters.flowArea*self.parameters.shalePermeability*\
            self.CO2MobilityShale
        DarcyCoefCO2Shale = DarcyCoefCO2*(1/self.parameters.shaleThickness)
        GravityCO2Shale = DarcyCoefCO2*self.parameters.CO2Density*self.g

        DarcyCoefBrine = self.parameters.flowArea*self.parameters.shalePermeability*\
            self.brineMobilityShale
        DarcyCoefBrineShale = DarcyCoefBrine*(1/self.parameters.shaleThickness)
        GravityBrineShale = DarcyCoefBrine*self.parameters.brineDensity*self.g

        # Inflow volumetric rate along shale layers = inflow volumetric rate
        # into the aquifer above the shale layer
        self.CO2Q = np.zeros(self.nSL)
        self.brineQ = np.zeros(self.nSL)

        self.CO2Q[:] = (DarcyCoefCO2Shale*deltaP - GravityCO2Shale)
        self.brineQ[:] = (DarcyCoefBrineShale*deltaP - GravityBrineShale)

        # Inflow mass rate along shale layers = inflow mass rate
        # into the aquifer above the shale layer
        self.CO2LeakageRates = np.zeros(self.nSL)
        self.brineLeakageRates = np.zeros(self.nSL)

        self.CO2LeakageRates = self.CO2Q*self.parameters.CO2Density
        self.brineLeakageRates = self.brineQ*self.parameters.brineDensity

        # Update cumulative CO2 "volume" in each aquifer (not including the bottom reservoir)
        self.CO2Volume = self.CO2Volume + (
            self.CO2Q[0:self.nSL-1]-self.CO2Q[1:self.nSL])*timeStep

        self.CO2Volume = np.maximum(self.CO2Volume, np.zeros(self.nSL-1))

        self.get_interface()
        for j in range(self.nSL-1):
            self.CO2SaturationAq[j] = self.interface[j+1]/self.thicknessH[j+1]

    def get_mobility_for_shales(self):
        """ Calculate mobility parameter for brine and CO2 in the shale layers."""
        CO2RelPermShale = np.zeros(self.nSL)
        brineRelPermShale = np.ones(self.nSL)

        Sres = self.parameters.aquBrineResSaturation

        for ind in range(0, self.nSL):
            if self.interface[ind] > 0:
                CO2RelPermShale[ind] = min(1-Sres[ind],
                                           self.interface[ind]/self.thicknessH[ind])
                brineRelPermShale[ind] = 1/(1-Sres[ind])*(1-Sres[ind]-CO2RelPermShale[ind])
            else:
                CO2RelPermShale[ind] = 0.
                brineRelPermShale[ind] = 1.

        self.CO2MobilityShale = CO2RelPermShale/self.parameters.CO2Viscosity
        self.brineMobilityShale = brineRelPermShale/self.parameters.brineViscosity

    def find_lambda(self):
        """ Calculate lambda constant."""
        # Lambda is a constant in the formula for the outer egde of the plume
        # It is constant for each aquifer in the system unless
        # brine and CO2 have different relative permeabilities, residual
        # saturations or viscosities
        # For simple calculations, the mobility ratio Lambda is calculated at
        # the endpoint relative permeability values for the two phases.

        Sres = self.parameters.aquBrineResSaturation

        # Maximum value for brine relative permeability
        brinePermeability = np.ones(self.nSL)

        # Maximum value for CO2 relative permeability
        CO2Permeability = np.zeros(self.nSL)
        for ind in range(self.nSL):
            CO2Permeability[ind] = 1 - Sres[ind]

        brineViscosity = self.parameters.brineViscosity
        CO2Viscosity = self.parameters.CO2Viscosity

        self.Lambda = np.ones(self.nSL)*CO2Permeability*brineViscosity/(
            brinePermeability*CO2Viscosity)

    def get_interface(self):
        """ Calculate height of the interface between brine and CO2. """
        Sres = self.parameters.aquBrineResSaturation

        # Set only for aquifers (not including the reservoir layer)
        for j in range(1, self.nSL):

            if self.CO2Volume[j-1] > 0:
                x = 2*np.pi*(self.thicknessH[j])*(1-Sres[j])*0.15*(
                    self.parameters.wellRadius**2)/self.CO2Volume[j-1]

                if x < 2/self.Lambda[j]:
                    tempInterface = (1-Sres[j])*self.thicknessH[j]

                elif x >= 2*self.Lambda[j]:
                    tempInterface = 0.
                else:
                    tempInterface = 1/(self.Lambda[j]-1)*(
                        np.sqrt(2*self.Lambda[j]/x)-1)*self.thicknessH[j]

                self.interface[j] = min([tempInterface,
                                         (1-Sres[j])*self.thicknessH[j],
                                         self.interface[j-1]])

    def CC_shale(self):

        PermWellShale = self.parameters.shalePermeability
        ShaleThickness = self.parameters.shaleThickness

        CO2MobilityShale = self.CO2MobilityShale
        brineMobilityShale = self.brineMobilityShale

        CC = PermWellShale/ShaleThickness*(CO2MobilityShale+brineMobilityShale)

        return CC

    def GG_shale(self):

        PermWellShale = self.parameters.shalePermeability
        gg = self.g

        CO2Density = self.parameters.CO2Density
        brineDensity = self.parameters.brineDensity

        CO2MobilityShale = self.CO2MobilityShale
        brineMobilityShale = self.brineMobilityShale

        GG = PermWellShale*gg*(CO2MobilityShale*CO2Density+brineMobilityShale*brineDensity)

        return GG

    def FF_aquifer(self):

        SCO2 = self.interface*(1/self.thicknessH)
        gg = self.g

        CO2Density = self.parameters.CO2Density
        brineDensity = self.parameters.brineDensity

        AquiferThickness = self.thicknessH

        FF = AquiferThickness/2*gg*(SCO2*CO2Density+(1-SCO2)*brineDensity)

        return FF

    def WV_aquifer(self):

        SCO2 = self.interface*(1/self.thicknessH)

        CO2MobilityAq = self.CO2MobilityAq
        brineMobilityAq = self.brineMobilityAq

        EffMobility = SCO2*CO2MobilityAq+(1-SCO2)*brineMobilityAq

        cf_ave = self.parameters.compressibility

        tt = self.parameters.timePoint*units.day()
        EffTime = 0.92*tt

        rw = self.parameters.wellRadius
        # Reservoir perm not defined separately
        PermAquiferHor = np.append(self.parameters.aquiferPermeability[0],
                                   self.parameters.aquiferPermeability)

        uLeak = rw**2*cf_ave/(4*EffMobility*PermAquiferHor*EffTime)
        GLeak = well_function(uLeak)

        AquiferThickness = self.thicknessH
        WV = GLeak/(4*np.pi*EffMobility*AquiferThickness*PermAquiferHor)

        return WV

def read_data(filename):
    """ Routine used for reading pressure and saturation data files."""
    # Check whether the file with given name exists
    if os.path.isfile(filename):
        data = np.genfromtxt(filename)
        return data
    return None


def well_function(x):
    """
    Return the approximation of the well function with even number
    of terms in expansion. Expansions with an odd number of terms
    diverge to plus infinity without crossing zero.
    """
    W = np.zeros(len(x))
    for i, u in list(enumerate(x)):
        if u <= 1.0:
            W[i] = (-0.577216-np.log(u)+u-u**2/(2*scm.factorial(2))+
                    u**3/(3*scm.factorial(3))-u**4/(4*scm.factorial(4))+
                    u**5/(5*scm.factorial(5))-u**6/(6*scm.factorial(6))+
                    u**7/(7*scm.factorial(7))-u**8/(8*scm.factorial(8)))
        elif u <= 9:
            uu = np.linspace(1.0, 9.0, num=9)
            Wu = np.array([0.219, 0.049, 0.013, 0.0038, 0.0011,
                           0.00036, 0.00012, 0.000038, 0.000012])
            Wfun = interp1d(uu, Wu, kind='linear')
            W[i] = Wfun(u)
        else:
            W[i] = 0.000001
    return W

if __name__ == "__main__":
    xx = [5.0, 6, 4]
    w_val = well_function(xx)
    print(w_val)
