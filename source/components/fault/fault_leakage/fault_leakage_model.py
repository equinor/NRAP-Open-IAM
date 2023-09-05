# -*- coding: utf-8 -*-
"""

"""
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
from tensorflow.keras.models import load_model
FL_MODELS = {
    'brine': load_model(os.sep.join([
        dir_path, 'brine_model'])),
    'CO2': load_model(os.sep.join([
        dir_path, 'CO2_model']))}

def fl_model(name):

    return FL_MODELS[name]
import joblib

FL_SCALERS = {
    'brine': joblib.load(os.sep.join([
        dir_path, 'brine_scaler'])),
    'CO2': joblib.load(os.sep.join([
        dir_path, 'CO2_scaler']))}
def fl_scaler(name):

    return FL_SCALERS[name]
