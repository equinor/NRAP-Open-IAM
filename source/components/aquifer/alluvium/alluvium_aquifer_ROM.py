"""
The module contains the solution class for the alluvium aquifer ROM

Authors: class Kayyum Mansoor

Date: 06/12/2018
Last modified: 06/12/2018

Component model input definitions:

* **co2_flux** * [kg/s] (0 to 0.5001) - CO2 flux (default 0.25005 kg/s)

* **co2_mass** * [|log10| (kg)] (2.23 to 9.058) - Cumulative CO2 mass (default
5.644 |log10| (kg))

* **sandFraction** * [-] (0.35 to 0.65) - Sand volume fraction (default 0.422)

* **correlationLengthX** * [m] (200 to 2500) - Correlation length in
X-direction (default 361.350 m)

* **correlationLengthZ** * [m] (0.5 to 25) - Correlation length in
Z-direction (default 17.382 m)

* **permeabilityClay** * [|log10| (|m^2|)] (-18 to -15) - Permeability in
clay (default -16.340 |log10| (|m^2|))

* **NaMolality** * [|log10| (Molality)] (-3 to 1) - Sodium molality (default
0.121 |log10| (Molality))

* **PbMolality** * [|log10| (Molality)] (-8.5 to -5) - Lead molality (default
-5.910 |log10| (Molality))

* **benzeneMolality** * [|log10| (Molality)] (0 to 1) - Benzene molality (default
0.109 |log10| (Molality))

* **tMitigation** * [years] (1 to 200) - Mitigation time (default 87.914 years)

* **brn_flux** * [kg/s] (0 to 0.075) - Brine flux (default 0.0375 kg/s)

* **brn_mass** * [|log10| (kg)] (4.5 to 8.124) - Cumulative brine mass (default
6.312 |log10| (kg))

* **simtime** * [years] (1 to 200) - Simulation Time (default 5.0 years)

* **CEC** * [meq/100g] (0.1 to 40) - The cation exchange capacity
(meq/100) (default 32.073 meq/100g)

* **Asbrine** * [|log10| (Molality)] (-9 to -5) - Arsenic concentration in the
leaking brine (mol/L) (default -5.397 |log10| (Molality))

* **Babrine** * [|log10| (Molality)] (-5.1 to -2.3) - Barium concentration in
the leaking brine (mol/L) (default -3.397 |log10| (Molality))

* **Cdbrine** * [|log10| (Molality)] (-9 to -6) - Cadmium concentration in the
leaking brine (mol/L) (default -8.574 |log10| (Molality))

* **Pbbrine** * [|log10| (Molality)] (-8.5 to -5) - Lead concentration in the
leaking brine (mol/L) (default -7.719 |log10| (Molality))

* **Benzene_brine** * [|log10| (Molality)] (-8.8927 to -4.8927) - Benzene
concentration in the leaking brine (mol/L) (default -8.610 |log10| (Molality))

* **Benzene_kd** * [L/kg] (-4.5 to 0.69) - Benzene distribution coefficient
(L/kg) (default -3.571 L/kg)

* **Benzene_decay** * [1/s] (-6.1 to 0) - Benzene degradation constant
(1/s) (default -2.732 1/s)

* **PAH_brine** * [mol/L] (-10 to -4.1) - Naphthalene concentration in the
leaking brine (mol/L) (default -7.118 mol/L)

* **PAH_kd** * [L/kg] (-3.1 to 1.98) - Naphthalene distribution coefficient
(L/kg) (default -0.985 L/kg)

* **PAH_decay** * [1/s] (-6.45 to 0) - Naphthalene degradation constant
(1/s) (default -3.371 1/s)

* **phenol_brine** * [mol/L] (-10 to -3.7) - Phenol concentration in the leaking
brine (mol/L) (default -6.666 mol/L)

* **phenol_kd** * [L/kg] (-6 to 0.15) - Phenol distribution coefficient
(L/kg) (default -1.342 L/kg)

* **phenol_decay** * [1/s] (-5.62999999999999 to 0) - Phenol degradation
constant (1/s) (default -3.546 1/s)

* **porositySand** * [-] (0.25 to 0.5) - Porosity [-] (default 0.375)

* **densitySand** * [kg/|m^3|] (1500 to 2500) - Density [kg/|m^3|] (default
2165.953 kg/|m^3|)

* **VG_mSand** * [-] (0.52 to 0.79) - Van Genuchten m [-](default 0.627)

* **VG_alphaSand** * [1/m] (-4.69 to -3.81) - Van Genuchten alpha [1/m] (default
-4.341 1/m)

* **permeabilitySand** * [|m^2|] (-14 to -10) - Permeability [|m^2|] (default
-12.430 |m^2|)

* **Clbrine** * [mol/L] (-2 to 0.73) - Chlorine concentration in the leaking
brine (mol/L) (default -0.339 mol/L)

* **calcitevol** * [-] (0.035 to 0.2) - Volume fraction of calcite [-] (default
0.165)

* **V_goethite** * [-] (0 to 0.2) - Volume fraction of goethite [-] (default
0.004)

* **V_illite** * [-] (0 to 0.3) - Volume fraction of illite [-] (default 0.006)

* **V_kaolinite** * [-] (0 to 0.2) - Volume fraction of kaolinite [-] (default
0.004)

* **V_smectite** * [-] (0 to 0.5) - Volume fraction of smectite [-] (default
0.010)


"""
import os
import numpy as np


class alluviumaqROMs(object):
    """ Class containing set of ROMs for AlluviumAquifer component. """
    def __init__(self, filename1, filename2):
        """ Constructor method. """
        self.filename1 = filename1
        self.filename2 = filename2

        self.loadroms(self.filename1, self.filename2)

    def loadroms(self, filename1, filename2):
        """ Load ROMs needed for the proper work of the component. """
        f1rommat = np.loadtxt('%s'%(filename1), delimiter=' ', comments="#", skiprows=1)
        self.f1const = f1rommat[:, -1]
        self.f1mat = f1rommat[:, 0:np.shape(f1rommat)[1]-1]

        self.f1pindices = []
        for i in range(np.int(self.f1mat.max())):
            self.f1pindices.append((self.f1mat == i+1.0).nonzero())
        self.f1zeroind = (self.f1mat == 0).nonzero()

        g1rommat = np.loadtxt('%s'%(filename2), delimiter=' ', comments="#", skiprows=1)
        self.g1const = g1rommat[:, -1]
        self.g1mat = g1rommat[:, 0:np.shape(g1rommat)[1]-1]

        self.g1pindices = []
        for i in range(np.int(self.g1mat.max())):
            self.g1pindices.append((self.g1mat == i+1.0).nonzero())
        self.g1zeroind = (self.g1mat == 0).nonzero()

    def romsolve(self, hydat, cfdat):
        """ Find solution. """

        hy = (self.f1const, self.f1mat, self.f1pindices, self.f1zeroind)
        cf = (self.g1const, self.g1mat, self.g1pindices, self.g1zeroind)

        parsmat = hy[1]*1.0
        adjparsmat = cf[1]*1.0

        for j, hydat_j in enumerate(hydat):
            parsmat[hy[2][j][0], hy[2][j][1]] = hydat_j

        parsmat[hy[3][0], hy[3][1]] = 1

        hyval = np.sum(np.prod(parsmat, axis=1)*hy[0])

        for j, cfdat_j in enumerate(cfdat):
            adjparsmat[cf[2][j][0], cf[2][j][1]] = cfdat_j

        adjparsmat[cf[3][0], cf[3][1]] = 1

        cfval = np.sum(np.prod(adjparsmat, axis=1)*cf[0])

        outval = (10**hyval) * cfval

        # Constrain values to output between model cell (50x50x4=10,000m3)
        # or total model volume (10,000x5,000x240=1.2e10m3)
        if (outval < 1.00e5) or (outval > 1.2e10):
            outval = np.nan

        return outval


class Solution(object):
    """ Solution class for alluvium aquifer ROM."""
    def __init__(self):
        """ Create Solution class object. """
        global cocs
        cocs = ['tds', 'ph', 'as', 'pb', 'cd', 'ba', 'benzene', 'pah', 'phenol']

        self.roms = {}
        self.Outputs = np.zeros(9)

        file_path = os.path.abspath(__file__)

        self.roms['tds'] = alluviumaqROMs(
            os.path.join(os.path.dirname(file_path), 'rom_nuft_em_tds.rom'),
            os.path.join(os.path.dirname(file_path), 'rom_cfact_tds.rom'))
        self.roms['ph'] = alluviumaqROMs(
            os.path.join(os.path.dirname(file_path), 'rom_nuft_em_ph.rom'),
            os.path.join(os.path.dirname(file_path), 'rom_cfact_ph.rom'))
        self.roms['as'] = alluviumaqROMs(
            os.path.join(os.path.dirname(file_path), 'rom_nuft_em_ph.rom'),
            os.path.join(os.path.dirname(file_path), 'rom_cfact_as.rom'))
        self.roms['pb'] = alluviumaqROMs(
            os.path.join(os.path.dirname(file_path), 'rom_nuft_em_ph.rom'),
            os.path.join(os.path.dirname(file_path), 'rom_cfact_pb.rom'))
        self.roms['cd'] = alluviumaqROMs(
            os.path.join(os.path.dirname(file_path), 'rom_nuft_em_ph.rom'),
            os.path.join(os.path.dirname(file_path), 'rom_cfact_cd.rom'))
        self.roms['ba'] = alluviumaqROMs(
            os.path.join(os.path.dirname(file_path), 'rom_nuft_em_ph.rom'),
            os.path.join(os.path.dirname(file_path), 'rom_cfact_ba.rom'))
        self.roms['benzene'] = alluviumaqROMs(
            os.path.join(os.path.dirname(file_path), 'rom_nuft_em_benzene.rom'),
            os.path.join(os.path.dirname(file_path), 'rom_cfact_benzene.rom'))
        self.roms['pah'] = alluviumaqROMs(
            os.path.join(os.path.dirname(file_path), 'rom_nuft_em_pah.rom'),
            os.path.join(os.path.dirname(file_path), 'rom_cfact_pah.rom'))
        self.roms['phenol'] = alluviumaqROMs(
            os.path.join(os.path.dirname(file_path), 'rom_nuft_em_phenol.rom'),
            os.path.join(os.path.dirname(file_path), 'rom_cfact_phenol.rom'))

    def find(self, xp):
        """ Find solution. """
        hy_ph = np.array([xp[2], xp[3], xp[4], xp[31], xp[5], xp[6], xp[7],
                          xp[8], xp[9], xp[0], xp[1], xp[12]])
        hy_tds = np.array([xp[2], xp[3], xp[4], xp[31], xp[5], xp[6], xp[7],
                           xp[8], xp[9], xp[10], xp[11], xp[12]])
        cf_ph = np.array([xp[27], xp[28], xp[29], xp[30], xp[31], xp[33],
                          xp[34], xp[35], xp[36], xp[37], xp[13], xp[12]])
        cf_tds = np.array([xp[27], xp[28], xp[29], xp[30], xp[31], xp[32],
                           xp[33], xp[34], xp[35], xp[36], xp[37], xp[12]])
        cf_as = np.array([xp[27], xp[28], xp[29], xp[30], xp[31], xp[14],
                          xp[33], xp[34], xp[35], xp[36], xp[37], xp[12]])
        cf_ba = np.array([xp[27], xp[28], xp[29], xp[30], xp[31], xp[15],
                          xp[33], xp[34], xp[35], xp[36], xp[37], xp[12]])
        cf_cd = np.array([xp[27], xp[28], xp[29], xp[30], xp[31], xp[16],
                          xp[33], xp[34], xp[35], xp[36], xp[37], xp[12]])
        cf_pb = np.array([xp[27], xp[28], xp[29], xp[30], xp[31], xp[17],
                          xp[33], xp[34], xp[35], xp[36], xp[37], xp[12]])
        cf_benzene = np.array([xp[27], xp[28], xp[29], xp[30], xp[31],
                               xp[18], xp[19], xp[20], xp[12]])
        cf_pah = np.array([xp[27], xp[28], xp[29], xp[30], xp[31],
                           xp[21], xp[22], xp[23], xp[12]])
        cf_phenol = np.array([xp[27], xp[28], xp[29], xp[30], xp[31],
                              xp[24], xp[25], xp[26], xp[12]])

        self.Outputs[0] = self.roms['tds'].romsolve(hy_tds, cf_tds)
        self.Outputs[1] = self.roms['ph'].romsolve(hy_ph, cf_ph)
        self.Outputs[2] = self.roms['as'].romsolve(hy_ph, cf_as)
        self.Outputs[3] = self.roms['pb'].romsolve(hy_ph, cf_pb)
        self.Outputs[4] = self.roms['cd'].romsolve(hy_ph, cf_cd)
        self.Outputs[5] = self.roms['ba'].romsolve(hy_ph, cf_ba)
        self.Outputs[6] = self.roms['benzene'].romsolve(hy_tds, cf_benzene)
        self.Outputs[7] = self.roms['pah'].romsolve(hy_tds, cf_pah)
        self.Outputs[8] = self.roms['phenol'].romsolve(hy_tds, cf_phenol)


if __name__ == "__main__":

    # Hydrology ROM pH - Normalized
    co2_flux = 0.002048
    co2_mass = 0.315539
    sandFraction = 0.238806
    correlationLengthX = 0.070152
    correlationLengthZ = 0.689042
    permeabilityClay = 0.553339
    NaMolality = 0.780295
    PbMolality = 0.739868
    benzeneMolality = 0.109203
    tMitigation = 0.436755
    brn_flux = 0.293401
    brn_mass = 0.394501
    simtime = 0.005
    CEC = 0.801316
    Asbrine = 0.900818
    Babrine = 0.608352
    Cdbrine = 0.141883
    Pbbrine = 0.223084
    Benzene_brine = 0.0706014
    Benzene_kd = 0.179084
    Benzene_decay = 0.552091
    PAH_brine = 0.48845
    PAH_kd = 0.416288
    PAH_decay = 0.47743
    phenol_brine = 0.529151
    phenol_kd = 0.757335
    phenol_decay = 0.370127
    porositySand = 0.873237
    densitySand = 0.665953
    VG_mSand = 0.394548
    VG_alphaSand = 0.396608
    permeabilitySand = 0.39261
    Clbrine = 0.608252
    calcitevol = 0.788496
    V_goethite = 0.02
    V_illite = 0.02
    V_kaolinite = 0.02
    V_smectite = 0.02

    oparam = np.array([
        co2_flux, co2_mass, sandFraction, correlationLengthX, correlationLengthZ,
        permeabilityClay, NaMolality, PbMolality, benzeneMolality, tMitigation,
        brn_flux, brn_mass, simtime, CEC, Asbrine, Babrine, Cdbrine, Pbbrine,
        Benzene_brine, Benzene_kd, Benzene_decay, PAH_brine, PAH_kd, PAH_decay,
        phenol_brine, phenol_kd, phenol_decay, porositySand, densitySand, VG_mSand,
        VG_alphaSand, permeabilitySand, Clbrine, calcitevol, V_goethite, V_illite,
        V_kaolinite, V_smectite])

    sol = Solution()
    sol.find(oparam)

    print('results:')
    for v in sol.Outputs:
        print(v)

    # The answer should be
    # results:
    # 226684.460033
    # 510826242.172
    # 32379693.3343
    # 556614383.306
    # 191831077.986
    # 540048193.753
    # nan
    # nan
    # nan
