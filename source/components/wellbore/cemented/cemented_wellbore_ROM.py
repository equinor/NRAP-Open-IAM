"""
The module contains the solution class for the Cemented Wellbore model.

Authors: class CementedWellboreRS Dylan Harp, Jaileen Del Valle Maldonado
         class Solution Veronika Vasylkivska

Last modified: 12/05/2022
"""
import os
import numpy as np


class CementedWellboreRS():
    """ Class for wellbore response surface model. """
    def __init__(self, filename):
        """ Constructor method for CementedWellboreRS class. """
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
        for i in  range(3):
            self._trans[i] = int(next(n))

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
            4 - pressure change,
            5 - first derivative of pressure,
            6 - second derivative of pressure,
            7 - saturation,
            8 - first derivative of saturation,
            9 - second derivative of saturation
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
        self.roms['CO2ThiefZone1'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2thf_10.h'))
        self.roms['CO2Aquifer1'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2aqu_10.h'))
        self.roms['CO2Atmosphere1'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2atm_10.h'))
        self.roms['brineThiefZone1'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2thf_10.h'))
        self.roms['brineAquifer1'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2aqu_10.h'))

        self.roms['CO2ThiefZone2'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2thf_11.h'))
        self.roms['CO2Aquifer2'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2aqu_11.h'))
        self.roms['CO2Atmosphere2'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2atm_11.h'))
        self.roms['brineThiefZone2'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2thf_11.h'))
        self.roms['brineAquifer2'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2aqu_11.h'))

        self.roms['CO2ThiefZone3'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2thf_12.h'))
        self.roms['CO2Aquifer3'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2aqu_12.h'))
        self.roms['CO2Atmosphere3'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2atm_12.h'))
        self.roms['brineThiefZone3'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2thf_12.h'))
        self.roms['brineAquifer3'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2aqu_12.h'))

        self.roms['CO2ThiefZone4'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2thf_13.h'))
        self.roms['CO2Aquifer4'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2aqu_13.h'))
        self.roms['CO2Atmosphere4'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2atm_13.h'))
       	self.roms['brineThiefZone4'] = CementedWellboreRS(
               os.path.join(header_file_directory, 'rsm_fw_wb2thf_13.h'))
        self.roms['brineAquifer4'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2aqu_13.h'))

        self.roms['CO2ThiefZone'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2thf_05.h'))
        self.roms['CO2Aquifer'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2aqu_05.h'))
        self.roms['CO2Atmosphere'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fc_wb2atm_05.h'))
        self.roms['brineThiefZone'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2thf_05.h'))
        self.roms['brineAquifer'] = CementedWellboreRS(
            os.path.join(header_file_directory, 'rsm_fw_wb2aqu_05.h'))

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
        if inputArray[2] > -11:
        # Index corresponds to thief zone (1), shallow aquifer (2), atmosphere (3)
            self.CO2LeakageRates[0] = self.roms['CO2ThiefZone1'].mars_rsm(inputArray)
            self.brineLeakageRates[0] = self.roms['brineThiefZone1'].mars_rsm(inputArray)
            self.CO2LeakageRates[1] = self.roms['CO2Aquifer1'].mars_rsm(inputArray)
            self.brineLeakageRates[1] = self.roms['brineAquifer1'].mars_rsm(inputArray)
            self.CO2LeakageRates[2] = self.roms['CO2Atmosphere1'].mars_rsm(inputArray)

        if inputArray[2] <= -11 and inputArray[2] > -12:
            self.CO2LeakageRates[0] = self.roms['CO2ThiefZone2'].mars_rsm(inputArray)
            self.brineLeakageRates[0] = self.roms['brineThiefZone2'].mars_rsm(inputArray)
            self.CO2LeakageRates[1] = self.roms['CO2Aquifer2'].mars_rsm(inputArray)
            self.brineLeakageRates[1] = self.roms['brineAquifer2'].mars_rsm(inputArray)
            self.CO2LeakageRates[2] = self.roms['CO2Atmosphere2'].mars_rsm(inputArray)

        if inputArray[2] <= -12 and inputArray[2] > -13:
            self.CO2LeakageRates[0] = self.roms['CO2ThiefZone3'].mars_rsm(inputArray)
            self.brineLeakageRates[0] = self.roms['brineThiefZone3'].mars_rsm(inputArray)
            self.CO2LeakageRates[1] = self.roms['CO2Aquifer3'].mars_rsm(inputArray)
            self.brineLeakageRates[1] = self.roms['brineAquifer3'].mars_rsm(inputArray)
            self.CO2LeakageRates[2] = self.roms['CO2Atmosphere3'].mars_rsm(inputArray)

        if inputArray[2] <= -13:
            self.CO2LeakageRates[0] = self.roms['CO2ThiefZone4'].mars_rsm(inputArray)
            self.brineLeakageRates[0] = self.roms['brineThiefZone4'].mars_rsm(inputArray)
            self.CO2LeakageRates[1] = self.roms['CO2Aquifer4'].mars_rsm(inputArray)
            self.brineLeakageRates[1] = self.roms['brineAquifer4'].mars_rsm(inputArray)
            self.CO2LeakageRates[2] = self.roms['CO2Atmosphere4'].mars_rsm(inputArray)
# The following if statement needed only if one wants to run
# the script as well as to use it as an importable module
# The code below is mainly used for testing purposes.
if __name__ == "__main__":
    # Create solution to get access to the backend ROMs
    sol = Solution('header_files')
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
        print(sol._mins, '\n')
        print(sol._comb_mod_mins, '\n')
        print(sol._true_mins, '\n')
        print(sol._maxs, '\n')
        print(sol._comb_mod_maxs, '\n')
        print(sol._true_maxs, '\n')

    else:
        if test == 1:
            wellDepth = 1420.8
            depthRatio = 0.367117117117117
            logWellPerm = -13.018235
            logThiefZonePermeability = -12.988711
            deltaP = 4.4173754
            pressurePrime = 0.000162739219712479
            pressureDPrime = -4.29398445836781e-07
            saturation = 0.00249977071
            saturationPrime = 1.44342907597535e-05
            saturationDPrime = 6.87653445433427e-08

            # The answer should be
            # CO2: [-1.75549609e-05 -8.56839668e-07  1.82368337e-06]
            # Brine: [5.79721927e-05 3.50508137e-06 0.00000000e+00]

        else:
            wellDepth = 1420.8
            depthRatio = 0.0 # 0.367117117117117
            logWellPerm = -13.018235
            logThiefZonePermeability = -11.988711
            deltaP = 4.4533206
            pressurePrime = 8.31745379876666e-05
            pressureDPrime = -1.68318793773187e-07
            saturation = 0.0125745982
            saturationPrime = 5.11821232032857e-05
            saturationDPrime = 1.58576542465501e-07

            # The answer should be
            # CO2: [1.71905454e-05 -1.36546316e-06  1.82368337e-06]
            # Brine: [4.94859038e-05 3.53217598e-06 0.00000000e+00]

        input_array = np.array([wellDepth, depthRatio,
                                logWellPerm, logThiefZonePermeability,
                                deltaP, pressurePrime, pressureDPrime,
                                saturation, saturationPrime, saturationDPrime])

        sol.find(input_array)
        print((sol.CO2LeakageRates))
        print((sol.brineLeakageRates))
