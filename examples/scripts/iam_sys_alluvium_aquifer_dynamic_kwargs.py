# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 09:23:01 2018

@author: Seth King
AECOM supporting NETL
Seth.King@NETL.DOE.GOV
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel, AlluviumAquifer

if __name__ == '__main__':
    time_array = 365.25*np.arange(0.0, 3.0)
    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add alluvium aquifer model object and define parameters
    aa = sm.add_component_model_object(AlluviumAquifer(name='aa', parent=sm))

    aa.add_par('sandFraction', value=0.538806)
    aa.add_par('correlationLengthX', value=300.0)
    aa.add_par('correlationLengthZ', value=10.689042)
    aa.add_par('permeabilityClay', value=-17.0)
    aa.add_par('NaMolality', value=0.100, vary=False)
    aa.add_par('PbMolality', value=-6.000, vary=False)
    aa.add_par('benzeneMolality', value=0.200, vary=False)
    aa.add_par('tMitigation', value=87.914, vary=False)
    aa.add_par('CEC', value=32.073, vary=False)
    aa.add_par('Asbrine', value=-5.397, vary=False)
    aa.add_par('Babrine', value=-3.397, vary=False)
    aa.add_par('Cdbrine', value=-8.574, vary=False)
    aa.add_par('Pbbrine', value=-7.719, vary=False)
    aa.add_par('Benzene_brine', value=-8.010, vary=False)
    aa.add_par('Benzene_kd', value=-3.571, vary=False)
    aa.add_par('Benzene_decay', value=-2.732, vary=False)
    aa.add_par('PAH_brine', value=-7.118, vary=False)
    aa.add_par('PAH_kd', value=-0.985, vary=False)
    aa.add_par('PAH_decay', value=-3.371, vary=False)
    aa.add_par('phenol_brine', value=-6.666, vary=False)
    aa.add_par('phenol_kd', value=-1.342, vary=False)
    aa.add_par('phenol_decay', value=-3.546, vary=False)
    aa.add_par('porositySand', value=0.468, vary=False)
    aa.add_par('densitySand', value=2165.953, vary=False)
    aa.add_par('VG_mSand', value=0.627, vary=False)
    aa.add_par('VG_alphaSand', value=-4.341, vary=False)
    aa.add_par('permeabilitySand', value=-12.430, vary=False)
    aa.add_par('Clbrine', value=-0.339, vary=False)
    aa.add_par('calcitevol', value=0.165, vary=False)
    aa.add_par('V_goethite', value=0.004, vary=False)
    aa.add_par('V_illite', value=0.006, vary=False)
    aa.add_par('V_kaolinite', value=0.004, vary=False)
    aa.add_par('V_smectite', value=0.010, vary=False)

    aa.add_dynamic_kwarg('co2_rate', [0.0, 1.90e-02, 1.90e-02])
    aa.add_dynamic_kwarg('co2_mass', [0.0, 6.00e+05, 6.00e+05])
    aa.add_dynamic_kwarg('brine_rate', [0.0, 4.62e-03, 4.62e-03])
    aa.add_dynamic_kwarg('brine_mass', [0.0, 1.46e+05, 2.46e+05])

    # Add observations (output) from the alluvium aquifer model
    aa.add_obs('TDS_volume')
    aa.add_obs('pH_volume')
    aa.add_obs('As_volume')
    aa.add_obs('Pb_volume')
    aa.add_obs('Cd_volume')
    aa.add_obs('Ba_volume')

    # Run the system model
    sm.forward()

    # Print the observations
    print('TDS:',
          sm.collect_observations_as_time_series(aa, 'TDS_volume'))
    print('pH:',
          sm.collect_observations_as_time_series(aa, 'pH_volume'))
    print('As:',
          sm.collect_observations_as_time_series(aa, 'As_volume'))
    print('Pb:',
          sm.collect_observations_as_time_series(aa, 'Pb_volume'))
    print('Cd:',
          sm.collect_observations_as_time_series(aa, 'Cd_volume'))
    print('Ba:',
          sm.collect_observations_as_time_series(aa, 'Ba_volume'))
