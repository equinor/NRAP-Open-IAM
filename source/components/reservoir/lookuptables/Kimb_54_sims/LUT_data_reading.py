# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 10:24:31 2018

@author: kings
"""
import numpy as np

signature_data = np.genfromtxt('parameters_and_filenames.csv',
                               delimiter=",", dtype='str')

# The first row (except the last element) of the file contains names of the parameters
par_names = signature_data[0,0:-1]
num_pars = len(par_names)
par_dict = {}
for i, par in enumerate(par_names):
    par_dict[par] = np.unique(signature_data[1:, i])
    
print('Number of parameters: {0}'.format(num_pars))
for par in par_dict:
    print('Parameter {0} values: {1}'.format(par, par_dict[par]))
