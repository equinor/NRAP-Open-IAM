# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 17:46:36 2018

@author: Seth King, Kayyum Mansoor
AECOM supporting NETL
Seth.King@NETL.DOE.GOV
"""

import os
import sys
import numpy as np
sys.path.insert(0, os.sep.join(['..', '..', 'source']))

try:
    from openiam import SystemModel, DeepAlluviumAquifer
except ImportError as err:
    print('Unable to load IAM class module: '+str(err))



if __name__ == "__main__":
    # Create system model
    time_array = 365.25*np.arange(0, 6)*10

    sm_model_kwargs = {'time_array': time_array} # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add carbonate aquifer model object and define parameters
    daa = sm.add_component_model_object(
        DeepAlluviumAquifer(name='daa', parent=sm))

    #self.pars_bounds['perm_var'] = [0.017, 1.89]
    daa.add_par('logK_sand1', value=-11.92098495)
    daa.add_par('logK_sand2', value=-11.7198002)
    daa.add_par('logK_sand3', value=-11.70137252)
    daa.add_par('logK_caprock', value=-15.69758676)
    daa.add_par('correlationLengthX', value=1098.994284)
    daa.add_par('correlationLengthZ', value=79.8062668)
    daa.add_par('sandFraction', value=0.800121364)
    daa.add_par('groundwater_gradient', value=0.001333374)
    daa.add_par('leak_depth', value=885.5060281)

    #daa.model_kwargs['brine_rate'] = 3.21903E-05 # kg/s
    #daa.model_kwargs['brine_mass'] = 10**4.71081307 # kg
    #daa.model_kwargs['co2_rate'] = 0.060985038 # kg/s
    #daa.model_kwargs['co2_mass'] = 10**6.737803184 # kg

    daa.add_dynamic_kwarg(
        'co2_rate', [0, 5.21E-06, 2.31E-02, 2.40E-02, 2.55E-02, 2.50E-02]) # kg/s
    daa.add_dynamic_kwarg(
        'co2_mass', [0, 744.5, 5282919.9, 12773616.8, 20674017.0, 28667754.3]) # kg
    daa.add_dynamic_kwarg(
        'brine_rate', [2.83E-04, 2.13E-04, 1.53E-05, 1.14E-05, 1.13E-05, 1.08E-05]) # kg/s
    daa.add_dynamic_kwarg(
        'brine_mass', [0, 95423.2, 120307.6, 124444.0, 128072.6, 131571.8]) # kg

    # Add observations (output) from the deep alluvium aquifer model
    # When add_obs method is called, system model automatically
    # (if no other options of the method is used) adds a list of observations
    # with names daa.obsnm_0, daa.obsnm_1,.... The final index is determined by the
    # number of time points in system model time_array attribute.
    # For more information, see the docstring of add_obs of ComponentModel class.

    daa.add_obs('TDS_volume')
    daa.add_obs('TDS_dx')
    daa.add_obs('TDS_dy')
    daa.add_obs('TDS_dz')

    daa.add_obs('Pressure_volume')
    daa.add_obs('Pressure_dx')
    daa.add_obs('Pressure_dy')
    daa.add_obs('Pressure_dz')

    daa.add_obs('pH_volume')
    daa.add_obs('pH_dx')
    daa.add_obs('pH_dy')
    daa.add_obs('pH_dz')

    # Run the system model
    sm.forward()

    # Print the observations
    # Since the observations at the particular time points are different variables,
    # method collect_observations_as_time_series creates lists of
    # values of observations belonging to a given component (e.g. daa) and having the same
    # common name (e.g. 'pH_volume', 'TDS_volume', etc) but differing in indices.
    # More details are given in the docstring and documentation to the method
    # collect_observations_as_time_series of SystemModel class.
    print('TDS_volume:',
          sm.collect_observations_as_time_series(daa, 'TDS_volume'))
    print('TDS_dx:',
          sm.collect_observations_as_time_series(daa, 'TDS_dx'))
    print('TDS_dy:',
          sm.collect_observations_as_time_series(daa, 'TDS_dy'))
    print('TDS_dz:',
          sm.collect_observations_as_time_series(daa, 'TDS_dz'))

    print('Pressure_volume:',
          sm.collect_observations_as_time_series(daa, 'Pressure_volume'))
    print('Pressure_dx:',
          sm.collect_observations_as_time_series(daa, 'Pressure_dx'))
    print('Pressure_dy:',
          sm.collect_observations_as_time_series(daa, 'Pressure_dy'))
    print('Pressure_dz:',
          sm.collect_observations_as_time_series(daa, 'Pressure_dz'))

    print('pH_volume:',
          sm.collect_observations_as_time_series(daa, 'pH_volume'))
    print('pH_dx:',
          sm.collect_observations_as_time_series(daa, 'pH_dx'))
    print('pH_dy:',
          sm.collect_observations_as_time_series(daa, 'pH_dy'))
    print('pH_dz:',
          sm.collect_observations_as_time_series(daa, 'pH_dz'))

#    # Expected output
#	   ('TDS_volume', array([       0.        ,  1009973.75523013,  2242342.09707877,
#	           3297433.38960866,  4176475.93617435,  5841693.86612737]))
#	   ('TDS_dx', array([   0.        ,  159.46251178,  134.13364332,  213.04373878,
#	           273.18094045,  325.11268209]))
#	   ('TDS_dy', array([   0.        ,  321.87577129,  122.82622137,  166.29934698,
#	           198.06515247,  222.26184085]))
#	   ('TDS_dz', array([   0.        ,  140.61472767,  253.51582479,  279.97725954,
#	           308.04509969,  328.33396187]))
#	   ('Pressure_volume', array([       0.        ,   993067.54933028,  2317554.79709008,
#	           3488075.9550932 ,  5163117.4456041 ,  6582598.49910734]))
#	   ('Pressure_dx', array([   0.        ,  158.32149693,  158.77047715,  276.8057456 ,
#	           373.68868184,  458.7054945 ]))
#	   ('Pressure_dy', array([   0.        ,   96.93891216,  125.22432141,  199.34890935,
#	           243.33834316,  279.16736562]))
#	   ('Pressure_dz', array([   0.        ,   70.81306891,  227.40256112,  295.04750013,
#	           342.19643013,  375.34425888]))
#	   ('pH_volume', array([  0.00000000e+00,   1.13135063e+09,   3.03848935e+06,
#	            5.77699635e+06,   8.12671714e+06,   1.03111644e+07]))
#	   ('pH_dx', array([    0.        ,  4683.46320217,   122.53057074,   222.11159   ,
#	            294.28220977,   367.50629984]))
#	   ('pH_dy', array([    0.        ,  5026.60802585,    91.661079  ,   139.13506168,
#	            174.64203653,   204.04482001]))
#	   ('pH_dz', array([   0.        ,  156.97613507,  228.04737627,  252.98494287,
#	           267.82421979,  287.18483688]))
