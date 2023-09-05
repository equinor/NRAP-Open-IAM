"""
The module contains the solution class for the Cemented Wellbore WR model.

Authors: class WRCementedWellboreRS Dylan Harp, Jaileen Del Valle Maldonado
         class Solution Veronika Vasylkivska

Last modified: 12/05/2022
"""
import os
import numpy as np


class WRCementedWellboreRS():
    """ Class for wellbore response surface model. """
    def __init__(self, filename):
        """ Constructor method for WRCementedWellboreRS class. """
        self.filename = filename
        self.read(filename)

    def __repr__(self):
        """Return a string containing a printable representation of an object."""
        rs_str = '\nfile: '+self.filename +'\n'
        rs_str += 'Param: [Min, Max]\n'
        for n, mn, mx in zip(self.param, self._mins, self._maxs):
            rs_str += ''.join([n, ': [', str(mn), ', ', str(mx), ']'])
            rs_str += '\n'

        return rs_str

    def read(self, filename):
        """ Read model input data file. """
        f = open(filename, 'r')
        f.readline()
        self.param = f.readline().split(':')[1].split()
        f.readline()
        self._nt = int(f.readline().split()[1])
        f.readline()
        self._np = int(f.readline().split()[1])
        f.readline()
        self._ns = int(f.readline().split()[1])
        f.readline()

        n = self.readfile(f)
        self._cuts = np.zeros([self._nt, self._np])
        for i in range(self._np):
            next(n)
            self._scl = float(next(n))
            self._cuts[:, i] = np.array(
                [float(next(n)) for j in range(self._nt)])
            if self._scl != 0.0:
                self._cuts[:, i] += np.array(
                    [float(next(n)) for j in range(self._nt)])
                self._cuts[:, i] /= self._scl

        for i in range(4):
            next(n)

        self._dirs = np.zeros([self._nt, self._np])
        for i in range(self._np):
            self._dirs[:, i] = [int(next(n)) for j in range(self._nt)]

        for i in range(2):
            next(n)

        self._c = np.zeros(self._ns)
        for i in range(self._ns):
            self._c[i] = float(next(n))

        for i in range(3):
            next(n)

        self._terms = np.zeros(self._ns)
        for i in range(self._ns):
            self._terms[i] = int(next(n))

        for i in range(2):
            next(n)

        self._mins = np.zeros(self._np)
        for i in range(self._np):
            self._mins[i] = float(next(n))

        for i in range(2):
            next(n)

        self._maxs = np.zeros(self._np)
        for i in range(self._np):
            self._maxs[i] = float(next(n))

        for i in range(3):
            next(n)

        mars_inputs = np.zeros(self._np)
        for i in range(self._np):
            mars_inputs[i] = int(next(n))

        for i in range(4):
            next(n)

        self._trans = np.zeros(3)
        for i in range(3):
            self._trans[i] = float(next(n))

        f.close()

    def change_boundaries(self, ind, min_val, max_val):
        """
        Change the lower and upper boundaries of the model parameter ind.

        The method is used to modify the limits of the ROM. It can appear that
        after testing the model is found to be suitable for a larger or smaller
        (more probable) range of a particular parameter in comparison with
        the original range.

        :param ind: index of the model parameter. Possible values of ind:
            0 - wellDepth,
            1 - depth ratio,
            2 - log well permeability,
            3 - log thief zone permeability,
            4 - thickness of thief zone
            5 - thickness of aquifer
            6 - thickness of reservoir
            7 - pressure change,
            8 - first derivative of pressure,
            9 - second derivative of pressure,
            10 - saturation,
            11 - first derivative of saturation,
            12 - second derivative of saturation
        :type ind: int

        :param min_val: new minimum value that the parameter can not be less than
        :type min_val: float

        :param max_val: new maximum value that the parameter can not exceed
        :type max_val: float
        """
        if ind < self._np:  # ind should not exceed number of parameters of the given ROM
            self._mins[ind] = min_val
            self._maxs[ind] = max_val
        else:
            raise ValueError("Argument ind exceeds the number of the ROM parameters.")

    def mars_rsm(self, inputs):
        '''
        Response surface calculator

        :param inputs: input parameter for response surface
        :type inputs: lst(float)

        :returns: float
        '''
        y = 0.0
        bx = np.ones(self._ns)
        temp = 0.0

        x = np.array(inputs)
        for i in range(self._np):
            if x[i] < self._mins[i]:
                x[i] = self._mins[i]
            if x[i] > self._maxs[i]:
                x[i] = self._maxs[i]

        for i in range(self._ns):
            for j in range(self._np):
                if self._dirs[int(self._terms[i]-1), j] == 2:
                    bx[i] = bx[i] * x[j]
                elif (self._dirs[int(self._terms[i]-1), j] == -1) or (
                        self._dirs[int(self._terms[i]-1), j] == 1):
                    temp = self._dirs[int(self._terms[i]-1), j]*(
                        x[j] - self._cuts[int(self._terms[i]-1), j])
                    if temp > 0:
                        bx[i] = bx[i] * temp
                    else:
                        bx[i] = 0
            y = y + bx[i] * self._c[i]
        y = y + self._trans[0]

        if self._trans[1] != 0:
            y = 10**y
        y = y + self._trans[2]
        return y

    def readfile(self, filehandle):
        """ Read file. """
        for line in filehandle:
            vs = line.split()
            for v in vs:
                yield v


class Solution():
    """ Solution class for cemented wellbore ROM."""
    def __init__(self, header_file_directory):
        """ Create an instance of the Solution class."""
        self.roms = {}
        self.roms['CO2ThiefZone1'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2thf_k1_3.h'))
        self.roms['CO2Aquifer1'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2aqu_k1_3.h'))
        self.roms['CO2Atmosphere1'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2atm_k1_3.h'))
        self.roms['brineThiefZone1'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2thf_k1_3.h'))
        self.roms['brineAquifer1'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2aqu_k1_3.h'))

        self.roms['CO2ThiefZone2'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2thf_k2_3.h'))
        self.roms['CO2Aquifer2'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2aqu_k2_3.h'))
        self.roms['CO2Atmosphere2'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2atm_k2_3.h'))
        self.roms['brineThiefZone2'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2thf_k2_3.h'))
        self.roms['brineAquifer2'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2aqu_k2_3.h'))

        self.roms['CO2ThiefZone3'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2thf_k3_3.h'))
        self.roms['CO2Aquifer3'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2aqu_k3_3.h'))
        self.roms['CO2Atmosphere3'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2atm_k3_3.h'))
        self.roms['brineThiefZone3'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2thf_k3_3.h'))
        self.roms['brineAquifer3'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2aqu_k3_3.h'))

        self.roms['CO2ThiefZone4'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2thf_k4_3.h'))
        self.roms['CO2Aquifer4'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2aqu_k4_3.h'))
        self.roms['CO2Atmosphere4'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2atm_k4_3.h'))
        self.roms['brineThiefZone4'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2thf_k4_3.h'))
        self.roms['brineAquifer4'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2aqu_k4_3.h'))

        self.roms['CO2ThiefZone'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2thf_3.h'))
        self.roms['CO2Aquifer'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2aqu_3.h'))
        self.roms['CO2Atmosphere'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2atm_3.h'))
        self.roms['brineThiefZone'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2thf_3.h'))
        self.roms['brineAquifer'] = WRCementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2aqu_3.h'))

        self.roms_names = ['CO2ThiefZone', 'CO2Aquifer', 'CO2Atmosphere',
                           'brineThiefZone', 'brineAquifer']
        self.num_subroms = 4

        # Placeholder for the future results
        self.CO2LeakageRates = None
        self.brineLeakageRates = None

    def __repr__(self):
        str_repr = ''
        for key in self.roms:
            str_repr = str_repr+self.roms[key].__repr__()
        return str_repr

    def get_mins_and_maxs(self):
        mins = []
        maxs = []
        comb_mod_mins = []
        comb_mod_maxs = []

        # Get mins and maxs of individual models
        for nm in [self.roms_names[0]]:
            for ind in range(1, self.num_subroms+1):
                rom = self.roms[nm+str(ind)]
                for i in range(rom._np):
                    mins.append(rom._mins[i])
                    maxs.append(rom._maxs[i])

            for i in range(self.roms[nm]._np):
                comb_mod_mins.append(self.roms[nm]._mins[i])
                comb_mod_maxs.append(self.roms[nm]._maxs[i])

        self._mins = np.array(mins).reshape(self.num_subroms, -1)
        self._maxs = np.array(maxs).reshape(self.num_subroms, -1)
        self._comb_mod_mins = np.array(comb_mod_mins)
        self._comb_mod_maxs = np.array(comb_mod_maxs)

        # Choose the maximum minimum and minimum maximum across all models
        self._true_mins = np.amax(self._mins, axis=0)
        self._true_maxs = np.amin(self._maxs, axis=0)

    def find(self, inputArray):
        """ Find solution of the ROM corresponding to the provided parameters."""
        self.CO2LeakageRates = np.zeros(3)
        self.brineLeakageRates = np.zeros(3)
        if -11 <= inputArray[2] < -10:
        # Index corresponds to thief zone (1), shallow aquifer (2), atmosphere (3)
            self.CO2LeakageRates[0] = self.roms['CO2ThiefZone1'].mars_rsm(inputArray)
            self.brineLeakageRates[0] = self.roms['brineThiefZone1'].mars_rsm(inputArray)
            self.CO2LeakageRates[1] = self.roms['CO2Aquifer1'].mars_rsm(inputArray)
            self.brineLeakageRates[1] = self.roms['brineAquifer1'].mars_rsm(inputArray)
            self.CO2LeakageRates[2] = self.roms['CO2Atmosphere1'].mars_rsm(inputArray)

        if -12 <= inputArray[2] < -11:
            self.CO2LeakageRates[0] = self.roms['CO2ThiefZone2'].mars_rsm(inputArray)
            self.brineLeakageRates[0] = self.roms['brineThiefZone2'].mars_rsm(inputArray)
            self.CO2LeakageRates[1] = self.roms['CO2Aquifer2'].mars_rsm(inputArray)
            self.brineLeakageRates[1] = self.roms['brineAquifer2'].mars_rsm(inputArray)
            self.CO2LeakageRates[2] = self.roms['CO2Atmosphere2'].mars_rsm(inputArray)

        if -13 <= inputArray[2] < -12:
            self.CO2LeakageRates[0] = self.roms['CO2ThiefZone3'].mars_rsm(inputArray)
            self.brineLeakageRates[0] = self.roms['brineThiefZone3'].mars_rsm(inputArray)
            self.CO2LeakageRates[1] = self.roms['CO2Aquifer3'].mars_rsm(inputArray)
            self.brineLeakageRates[1] = self.roms['brineAquifer3'].mars_rsm(inputArray)
            self.CO2LeakageRates[2] = self.roms['CO2Atmosphere3'].mars_rsm(inputArray)

        if -16 <= inputArray[2] < -13:
            self.CO2LeakageRates[0] = self.roms['CO2ThiefZone4'].mars_rsm(inputArray)
            self.brineLeakageRates[0] = self.roms['brineThiefZone4'].mars_rsm(inputArray)
            self.CO2LeakageRates[1] = self.roms['CO2Aquifer4'].mars_rsm(inputArray)
            self.brineLeakageRates[1] = self.roms['brineAquifer4'].mars_rsm(inputArray)
            self.CO2LeakageRates[2] = self.roms['CO2Atmosphere4'].mars_rsm(inputArray)

        if inputArray[4] == 0.0: # zero thief zone thickness
            self.CO2LeakageRates[0] = 0.0
            self.brineLeakageRates[0] = 0.0

# The following if statement needed only if one wants to run
# the script as well as to use it as an importable module
# The code below is mainly used for testing purposes.
if __name__ == "__main__":
    # Create solution to get access to the backend ROMs
    sol = Solution('header_files_wr')
    test = 0
    if test == 0:
        # Get mins and maxs of each individual ROMS as well as the joint ROM
        for ind in range(1, 5):
            for nm in sol.roms_names:
                # Use built-in method __repr__ to print information about each ROM
                print(sol.roms[nm+str(ind)])

        for nm in sol.roms_names:
            # Use built-in method __repr__ to print information about joint ROM
            print(sol.roms[nm])

        # Compare mins and maxs of all parameters except the one used to split
        # into subroms: well permeability has index 2

        sol.get_mins_and_maxs()
        # print(sol._mins, '\n')
        print(sol._comb_mod_mins, '\n')
        # print(sol._true_mins, '\n')
        # print(sol._maxs, '\n')
        print(sol._comb_mod_maxs, '\n')
        # print(sol._true_maxs, '\n')

    else:
        if test == 1:
            wellDepth = 1420.8
            depthRatio = 0.367117117117117
            logWellPerm = -13.018235
            logThiefZonePermeability = -12.988711
            thiefZoneThickness = 50.0
            aquiferThickness = 45.0
            reservoirThickness = 70.0
            deltaP = 4.4173754
            pressurePrime = 0.000162739219712479
            pressureDPrime = -4.29398445836781e-07
            saturation = 0.00249977071
            saturationPrime = 1.44342907597535e-05
            saturationDPrime = 6.87653445433427e-08

            # The answer:
            # Current model
            # CO2: [1.21922931e-08 7.11796428e-10 5.11781142e-10]
            # Brine: [7.41216880e-06 6.56564176e-08 0.00000000e+00]
            # Previous model
            # CO2: [6.23635636e-08 9.02336945e-09 9.62457009e-09]
            # Brine: [9.74088449e-06 6.46375218e-08 0.00000000e+00]

        else:
            wellDepth = 1420.8
            depthRatio = 0.367117117117117
            logWellPerm = -13.018235
            logThiefZonePermeability = -11.988711
            deltaP = 4.4533206
            thiefZoneThickness = 50.0
            aquiferThickness = 45.0
            reservoirThickness = 70.0
            pressurePrime = 8.31745379876666e-05
            pressureDPrime = -1.68318793773187e-07
            saturation = 0.0125745982
            saturationPrime = 5.11821232032857e-05
            saturationDPrime = 1.58576542465501e-07

            # The answer:
            # Current model:
            # CO2: [3.61964363e-25 3.10462210e-08 1.31998535e-08]
            # Brine: [7.39264924e-06 2.97396152e-07 0.00000000e+00]

            # Previous model
            # CO2: [4.96374746e-07 5.28094267e-08 4.99247008e-07]
            # Brine: [1.34401770e-04 6.46375218e-08 0.00000000e+00]

        input_array = np.array([wellDepth, depthRatio,
                                logWellPerm, logThiefZonePermeability,
                                thiefZoneThickness, aquiferThickness,
                                reservoirThickness, deltaP,
                                pressurePrime, pressureDPrime,
                                saturation, saturationPrime, saturationDPrime])

        sol.find(input_array)
        print('CO2:', (sol.CO2LeakageRates))
        print('Brine:', (sol.brineLeakageRates))
