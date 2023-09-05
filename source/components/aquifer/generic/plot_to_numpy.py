import os
import sys
import argparse
from math import ceil
import numpy as np
import glob

# STOMP run directory
class STOMP_Run:

    def __init__(self, dir_name):

        self.files = glob.glob(os.path.join(dir_name, 'plot.*'))
        self.files.sort()
        pf = Plot_File(self.files[0])
        self.nx = pf.nx
        self.ny = pf.ny
        self.nz = pf.nz
        self.nfield = pf.nfield
        self.nactive = pf.nactive
        self.nvert = pf.nvert
        self.input_file = Input_File(glob.glob(os.path.join(dir_name, 'input'))[0])

    # get input parameter
    def get_source_rate(self, name):
        return self.input_file.get_source_rate(name)

    # get variable values from all plot files in directory
    def get_variable(self, name, unit_conversion=None, func=None):
        variables = {}
        # read variables and times from each plot file
        for file_path in self.files:
            # create plot file object
            pf = Plot_File(file_path)
            time = pf.read_time('yr')
            var = pf.read_variable(
                name=name, unit_conversion=unit_conversion, func=func)
            variables.update({time : var})
        return variables

    def get_variable_change(self, name, base_time, unit_conversion=None, func=None):
        variables = self.get_variable(name, unit_conversion, func)
        times = self.get_times()
        changes = {}
        for key, var_val in variables.items():
            changes[key] = var_val - variables[base_time]
        return changes

    # get node centroids
    def get_centroids(self, direction):
        centroids = {}
        for file_path in self.files:
            pf = Plot_File(file_path)
            time = pf.read_time('yr')
            cent = pf.read_centroids(direction)
            centroids.update({time : cent})
        return centroids

    # get all plot file times
    def get_times(self):
        times = []
        for file_path in self.files:
            pf = Plot_File(file_path)
            time = pf.read_time('yr')
            times.append([time])
        return times


# STOMP plot file
class Plot_File:

    convert_length = {'m':1.0, 'ft':0.3048, 'km':1000., 'mi':1609.34,
                      'mile':1609.34, 'nm':1e-9, 'cm':0.01, 'mm':0.001,
                      'yd':0.9144, 'in':0.0254}
    convert_pressure = {'pa':1.0, 'psi':6894.76, 'mpa':1000000.,
                        'atm':101325., 'bar':100000.}

    def __init__(self, file_name):

        self.file_name = file_name
        self.input_list = self.read_file(self.file_name)
        self.nx = self.read_parameter('Number of X or R-Direction Nodes =')
        self.ny = self.read_parameter('Number of Y or Theta-Direction Nodes =')
        self.nz = self.read_parameter('Number of Z-Direction Nodes =')
        self.nfield = self.read_parameter('Number of Field Nodes =')
        self.nactive = self.read_parameter('Number of Active Nodes =')
        self.nvert = self.read_parameter('Number of Vertices =')

    # read file into list
    def read_file(self, file_name):
        input_list = []
        with open(file_name) as f:
            input_list = f.readlines()
            # remove whitespace characters
            input_list = [line.strip() for line in input_list]
        return input_list

    # read x, y, or z dimension
    def read_parameter(self, search_string):
        for line in self.input_list:
            if search_string in line:
                # remove whitespace characters
                line_list = [item.strip() for item in line.split('=')]
        return int(line_list[-1])

    # get line number with variable name
    def find_variable(self, search_string):
        m = 0
        for line in self.input_list:
            if search_string in line:
                break
            m = m + 1
        return m

    # read or calculate node centroids
    def read_centroids(self, direction):
        zlist = []
        if direction == 'Z':
            m = self.find_variable(direction + ' Node-Centroid Position')
            if m < len(self.input_list) - 1:
                start = m+1
                end = start+int(ceil(self.nfield/10.))
                unit = self.input_list[m].split()[-1]
                for line in self.input_list[start:end]:
                    line_list = [
                        float(item.strip())*self.convert_length[unit] for item in line.split()]
                    zlist.extend(line_list)
        elif direction == 'X':
            m = self.find_variable('X-Direction Nodal Vertices')
            start = m+1
            end = start+self.nfield
            unit = self.input_list[m].split()[-1]
            for line in self.input_list[start:end]:
                line_list = [float(item.strip()) for item in line.split()]
                centroid = self.convert_length[unit] * (line_list[0] + line_list[1]) / 2.
                zlist.append(centroid)
        elif direction == 'Y':
            m = self.find_variable('Y-Direction Nodal Vertices')
            start = m+1
            end = start+self.nfield
            unit = self.input_list[m].split()[-1]
            for line in self.input_list[start:end]:
                line_list = [float(item.strip()) for item in line.split()]
                centroid = self.convert_length[unit] * (line_list[0] + line_list[2]) / 2.
                zlist.append(centroid)
        else:
            sys.exit('Plot file must contain Z Node-Centroid Position')
        zarray = np.array(zlist).reshape((self.nx, self.ny, self.nz), order='F')
        return zarray

    def read_variable(self, name, unit_conversion=None, func=None):
        vlist = []
        m = self.find_variable(name)
        start = m+1
        end = start+int(ceil(self.nfield/10.))
        unit = self.input_list[m].split()[-1]
        if unit_conversion is not None:
            conversion_factor = unit_conversion[unit]
        else:
            conversion_factor = 1.0
        for line in self.input_list[start:end]:
            line_list = [float(item.strip())*conversion_factor for item in line.split()]
            vlist.extend(line_list)
        varray = np.array(vlist).reshape((self.nx, self.ny, self.nz), order='F')
        if func is not None:
            varray = func(varray)
        return varray


    def read_time(self,unit):
        m = self.find_variable('Time = ')
        line_list = [item.strip() for item in self.input_list[m].split()]
        time_list = [i.split(',') for i in line_list[2:]]
        time_dict = {x[1]: x[0] for x in time_list}
        return float(time_dict[unit])

class Input_File:
    def __init__(self, file_name):

        self.file_name = file_name
        self.input_list = self.read_file(self.file_name)

    # read file into list
    def read_file(self, file_name):
        input_list = []
        with open(file_name) as f:
            input_list = f.readlines()
            # remove whitespace characters
            input_list = [line.strip() for line in input_list]
        return input_list

    # get line number with text string
    def find_text(self, search_string):
        m = 0
        for line in self.input_list:
            if search_string in line:
                break
            m = m + 1
        return m

    def get_source_rate(self, search_string):
        m = self.find_text(search_string)
        line = self.input_list[m+1]
        line_list = [item.strip() for item in line.split(',')]
        return float(line_list[4])

if __name__ == '__main__':

    # input arguments
    parser = argparse.ArgumentParser(description='Convert STOMP plot files into numpy arrays')
    parser.add_argument('--dir', default='.', help='directory with results of STOMP simulation')
    args = parser.parse_args()

    run = STOMP_Run(args.dir)

#  INPUT
    # x-permeability
    kx = run.get_variable('X-Direction Intrinsic Permeability')
    # z-permeability
    kz = run.get_variable('Z-Direction Intrinsic Permeability')
    # porosity
    por = run.get_variable('Diffusive Porosity')
    # depth
    z = run.get_centroids('Z')
    # thickness
    # calcite volume fraction
    cvf = run.get_variable('calcite Volume Fraction')
    # co2 leakage rate
    co2 = run.get_variable(name='CO2 Mass Source Rate', func=lambda x: x * 360)
    # brine leakage rate
    brine = run.get_variable(name='Water Mass Source Rate', func=lambda x: x * 360)
    # cumulative mass
    brinemass = run.get_variable_change(name='Water Mass Source Rate',
                                        base_time=100.0, func=lambda x: x * 360)
    co2mass = run.get_variable_change(name='CO2 Mass Source Rate',
                                      base_time=100.0, func=lambda x: x * 360)

    input_co2_rate = run.get_source_rate('Gas Mass Source')
    input_brine_rate = run.get_source_rate('Aqueous Mass Source')
    print(input_co2_rate, input_brine_rate)

# OUTPUT
    # pH
    ph = run.get_variable(name='Aqueous h+ Concentration',
                          func=lambda x: -np.log10(x))
    # TDS
    salt = run.get_variable('Aqueous Salt Mass Fraction')
    dens = run.get_variable(name='Aqueous Density', func=lambda x: x * 1000)
    tds = {key : value * dens[key] for key, value in salt.items() if key in dens}

    # saturation
    sat = run.get_variable('Gas Saturation')
    # pressure
    pres = run.get_variable('Gas Pressure')
