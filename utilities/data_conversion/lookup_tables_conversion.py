"""
Module contains methods batch-converting csv lookup tables into h5 lookup tables.

Example at the end of the file illustrates conversion of lookup tables
from csv format into h5 for 3 data sets available for testing and application
in NRAP-Open-IAM: Kimberlina, Kimberlina closed and FutureGen2.
"""
import csv
import os
import h5py
import numpy as np
import pandas as pd

def from_csv_to_h5(file_directory, time_file='time_points.csv',
                   parameter_filename='parameters_and_filenames.csv'):
    # Read file with signatures of interpolators and names of files with the corresponding data
    signature_data = np.genfromtxt(
        os.path.join(file_directory, parameter_filename),
        delimiter=",", dtype='str')

    new_signature_data = np.copy(signature_data)
    new_signature_data[1:, -1] = np.array([nm[:-3]+'h5' for nm in signature_data[1:, -1]])

    time_points = np.loadtxt(os.sep.join([file_directory, time_file]),
                             delimiter=",", dtype='f8') # in years

    for ind, file_nm in enumerate(signature_data[1:, -1]):
        print(file_nm)
        res_data_ind = 2
        csv_file_nm = os.sep.join([file_directory, file_nm])
        # Read headers fro csv data file
        with open(csv_file_nm, "r") as f:
            reader = csv.reader(f)
            data_headers = next(reader)

        data = pd.read_csv(csv_file_nm)

        h5_file_nm = os.sep.join([file_directory, new_signature_data[1+ind, -1]])
        with h5py.File(h5_file_nm, 'w') as hf:
            hf.create_dataset('t', data=time_points)
            if 'x' in data:
                hf.create_dataset('x', data=data['x'])
            elif '# x' in data:
                hf.create_dataset('x', data=data['# x'])
            hf.create_dataset('y', data=data['y'])
            if 'z' in data:
                hf.create_dataset('z', data=data['z'])
                res_data_ind = 3
            for col_nm in data_headers[res_data_ind:]:
                hf.create_dataset(col_nm, data=data[col_nm])

    with h5py.File(h5_file_nm, 'r') as hf:
        h5_keys = list(hf.keys())
    print('Keys of data within the new h5 files:')
    print(h5_keys)

    # Save new signature file
    np.savetxt(os.sep.join([file_directory, 'parameters_and_filenames_h5.csv']),
               new_signature_data, delimiter=",", fmt="%s")


if __name__ == '__main__':
    test_case = 1

    if test_case == 1:

        file_directory = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                      'lookuptables', 'Kimb_54_sims'])
        time_file = 'time_points.csv'
        parameter_filename = 'parameters_and_filenames.csv'
    elif test_case == 2:
        file_directory = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                      'lookuptables', 'FutureGen2', '1008_sims'])
        time_file = 'time_points.csv'
        parameter_filename = 'parameters_and_filenames.csv'
    elif test_case == 3:
        file_directory = os.sep.join(['..', '..', 'source', 'components', 'reservoir',
                                      'lookuptables', 'Kimb_closed_200_sims'])
        time_file = 'time_points.csv'
        parameter_filename = 'parameters_and_filenames.csv'

    from_csv_to_h5(file_directory, time_file=time_file,
                   parameter_filename=parameter_filename)
