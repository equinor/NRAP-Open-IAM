"""
The module contains the solution class for the generic reservoir ROM.

The Generic Reservoir component model is a reduced-order-model to predict
pressure and CO2 saturation of the top part of the reservoir during CO2 injection
(25 years) and post-injection period (50 years). The model is a machine learning
regression model fitted to the results of STOMP-CO2E multiphase flow transport
simulations using Random Forest. The model predicts only for top part of the
storage reservoir. Homogeneous reservoir model with radius 150 km was used.

Author: Seunghwan Baek, PNNL
Date: 08/10/2022
Last update: 08/10/2022

"""
import os
import sys
import numpy as np
import pandas as pd

componentPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(componentPath)


class Solution():
    """ Solution class for generic reservoir ROM."""
    def __init__(self):
        """ Create an instance of the Solution class. """
        self.Outputs = {}

    def find(self, inputArray, x_scaler_P, y_scaler_P, model_P,
             x_scaler_S, y_scaler_S, model_S):
        """ Find solution of the ROM corresponding to the provided parameters."""

        Press_surface = inputArray[9] # Pa
        Press_grad = inputArray[10] #

        IniBotPress = (Press_surface + Press_grad*1e3*inputArray[0])*1e-6 # MPa
        IniTopPress = (Press_surface + Press_grad*1e3*(inputArray[0]-inputArray[1]))*1e-6 # MPa

        Temp_surface = inputArray[8] # C
        Temp_grad = inputArray[3] # C/km

        T_bot = (Temp_surface + Temp_grad*1e-3*inputArray[0]) # C
        T_top = (Temp_surface + Temp_grad*1e-3*(inputArray[0]-inputArray[1])) # C

        inputArray[2] = (10**(inputArray[2]))*1e+15 # log10 m^2 to mD
        inputArray[4] = inputArray[4]*3.154e+7*1e-9 # kg/s to MMT/yr
        inputArray[6] = inputArray[6]/np.sqrt(np.pi) # checked ok
        distance = inputArray[-1]
        InputData = inputArray[:-4]

        InputData0 = np.transpose(
            np.append(InputData[:4], [T_bot, T_top, IniBotPress, IniTopPress]))
        InputData0 = np.append(InputData0, InputData[4:])

        InputData1 = np.zeros([1, 13])
        InputData1[:, :-1] = InputData0
        InputData1[:, -1] = distance

        Input_df = pd.DataFrame(data=InputData1, columns=[
                    'BottomDepth,m', 'ResThickness,m', 'ResPerm,mD', 'TempGrad,C/km',
                    'T_bot,C', 'T_top,C', 'IniBotPress,MPa', 'IniTopPress,MPa',
                    'InjRate,MMT/yr', 'IniSalinity', 'WellRadius,m', 'ResPoro', 'rCoord'])

        TargetVars = ['Pressure', 'Saturation']

        for TargetVar in TargetVars:
            if TargetVar == 'Pressure':
                Input_reduced = x_scaler_P.transform(Input_df)
                Output_reduced = model_P.predict(Input_reduced)
                Output = y_scaler_P.inverse_transform(Output_reduced)

                self.Outputs['{}'.format(TargetVar.lower())] = (10**Output + IniTopPress)*1e6 # Pa

            elif TargetVar == 'Saturation':
                Input_reduced = x_scaler_S.transform(Input_df)
                Output_reduced = model_S.predict(Input_reduced)
                Output = y_scaler_S.inverse_transform(Output_reduced)

                Output = np.where(Output < 1e-3, 0.0, Output)
                self.Outputs['{}'.format(TargetVar.lower())] = np.where(
                    Output < 1e-3, 0.0, Output)

def user_msg():
    msg1 = ''.join([
        'This is generic_reservoir_ROM for estimation of pressure ',
        'propagation and plume evolution in reservoirs caused ',
        'by CO2 injection. Developed by Seunghwan Baek, PNNL in 2022.'])
    msg2 = 'Check "generic_reservoir_component.py".'
    print(msg1)
    print(msg2)

if __name__ == "__main__":
    user_msg()
