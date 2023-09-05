#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module provides global file names & values across other modules.

Author: Ernest N. Lindner
Date: 08/18/2022

Module Name
    config.py

Contents
    config variables

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
# Directories for I/O
OUTPUT_DIRECTORY = "seal_output"        # Output directory - default
INPUT_DIRECTORY = "seal_input"          # Seal data input directory
RESERVOIR_DIR = "reservoir_input"       # Reservoir input directory

# Names of seal and fracture control YAML files
S_CONTROL = "ControlFile_seal.yaml"     # Name of seal control file
F_CONTROL = "ControlFile_frac.yaml"     # Name of fracture control file

# Reservoir Input File Names
SATURATION_NAME = "lookup_reservoir_co2_saturation.txt"  # Base saturation
PRESSURE_NAME = "lookup_reservoir_co2_pressure.txt"      # Base pressure
TIME_STEP_NAME = "time_points.txt"                       # Time steps

# Seal Input Files Names
ACTIVE_NAME = "seal_active.txt"         # Status input file
AREA_NAME = "seal_area.txt"             # Area input file
BASE_NAME = "seal_base_depth.txt"       # Depth input file
ENTRY_NAME = "seal_entry_pressure.txt"  # Threshold pressure input file
INFLUENCE_NAME = "seal_influence.txt"   # Influence input file
PERM_NAME = "seal_permeability.txt"     # Permeability input file
PRESS_NAME = "seal_top_pressure.txt"    # Top boundary pressure input
THICKNESS_NAME = "seal_thickness.txt"   # Thickness input file
X_COORD_NAME = "seal_X_coord.txt"       # X-coordinates file
Y_COORD_NAME = "seal_Y_coord.txt"       # Y-coordinates file
USER_NAME = "user_fracture_input"       # User-defined fractures - default

# Seal Output Names
SEAL_SUM_FILE = "_seal_summary"                     # Seal summary file name
FRAC_SUM_FILE = "_frac_summary"                     # Frac summary file name
RESULTS_FILE = "seal_results-simulation_No._"       # Time history name
BRINE_INTRO = "brine_grid_data-sim_No._"            # NPY/csv file for brine
CO2_INTRO = "CO2_grid_data-sim_No._"                # NPY/csv file for CO2
FRAC_NAME = "frac_results-sim_No._"                 # File name - fracs
TOTO_NAME = "total_permeability-sim_No._"           # Fle name - perm
DEBUG_FILE = "debug_results"                        # Debug file name

# Interpolation Tables
CO2_DENSITY_FILE = "data_co2_density.npy"           # CO2 density table
CO2_VISCOSITY_FILE = "data_co2_viscosity.npy"       # CO2 viscosity table
BRINE_DENSITY_FILE = "data_brine_density.npy"       # Brine density table
BRINE_VISCOSITY_FILE = "data_brine_viscosity.npy"   # Brine viscosity table
CO2_SOLUBILITY_FILE = "data_co2_solubility.npy"     # CO2 solubility table

# File Type Extensions
EXTENSION_TXT = ".txt"                              # Text file
EXTENSION_CSV = ".csv"                              # Comma-delimited file
EXTENSION_NPY = ".npy"                              # Data file
EXTENSION_PNG = ".png"                              # Figure file

# Global Variables
PLOT_NO = 1                                         # Default figure number
SEEDX = 10                                          # Default random seed
# Set SEEDX to >None< for random/different seed generation for each simulation


#
# -----------------------------------------------------------------------------
# - End of module
