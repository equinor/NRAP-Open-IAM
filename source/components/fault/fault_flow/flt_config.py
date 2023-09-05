#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module provides global values across other modules.

Author: Ernest N. Lindner
Date: 08/23/2022

Module Name
    flt_config

Contents (0)
    Global Variables

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
# Directories
OUTPUT_DIR = "fault_output"                        # Output directory
INPUT_DIR = "fault_input"                          # Fault data input directory
RESERVOIR_DIR = "reservoir_input"                  # Reservoir input directory

# Reservoir Input File Names
SATURATION_NAME = "lookup_reservoir_co2_saturation.txt"  # Base saturation
PRESSURE_NAME = "lookup_reservoir_co2_pressure.txt"      # Base pressure
TIME_STEP_NAME = "time_points.txt"                       # Time steps file name

# I/O File Names
PERM_NAME = "permeability_for_simulation_"         # Detailed file
PLATE_NAME = "plate-sim_"                          # Plate data file name
RESULTS_NAME = "fault_results-simulation_"         # Time history name prefix
APERTURE_NAME = "fault_apertures.txt"              # Aperture input file
SUMMARY_FILE = "_fault_ROM_summary"                # Run summary file
FAULT_EXIST_FILE = "_fault_existence_list"         # Run fault present file

# Debugging File Names
BRINE_INTRO = "Brine_Data_for_Sim_"                # Start - file of brine
CO2_INTRO = "CO2_Data_for_Sim_"                    # Start - file of CO2

# Interpolation Tables
CO2_DENSITY_FILE = "data_co2_density.npy"          # CO2 density table
CO2_VISCOSITY_FILE = "data_co2_viscosity.npy"      # CO2 viscosity table
BRINE_DENSITY_FILE = "data_brine_density.npy"      # Brine density table
BRINE_VISCOSITY_FILE = "data_brine_viscosity.npy"  # Brine viscosity table
CO2_SOLUBILITY_FILE = "data_co2_solubility.npy"    # CO2 solubility table

# File Type Extensions
EXTENSION_TXT = ".txt"                             # Text file
EXTENSION_CSV = ".csv"                             # Comma-delimited File
EXTENSION_NPY = ".npy"                             # Data file
EXTENSION_PNG = ".png"                             # Figure

# Tracking & Random Seed
PLOT_NO = 1
SEEDX = 12345    # Set to "None" for variable seed!

# -----------------------------------------------------------------------------
# - End of module
