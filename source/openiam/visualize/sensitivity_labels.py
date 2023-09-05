# -*- coding: utf-8 -*-
"""
Module containing dictionaries defining default x-axis and legend labels
for the figures produced by sensitivity_analysis.py.
"""

# Dictionary used for the labels of parameters in sensitivity_analysis.py.
# Labels that are particularly long include '\n' so they are on the split lines.
PAR_NAME_DICT = {
    # Stratigraphy component parameters
    'numberOfShaleLayers': 'Number of Shale Layers',
    'shaleThickness': 'Shale Thickness',
    'aquiferThickness': 'Aquifer Thickness',
    'reservoirThickness': 'Reservoir Thickness',
    'datumPressure': 'Datum Pressure',
    'depth': 'Depth', # I am intentionally keeping this vague, as depth is
    # reservoir top depth for Stratigraphy components and aquifer bottom depth
    # for FutureGen2Aqufier components
    # Reservoir component parameters
    # Simple Reservoir component
    'logResPerm': 'Log. Reservoir Permeability',
    'reservoirPorosity': 'Reservoir Porosity',
    'brineDensity': 'Brine Density',
    'CO2Density': 'CO$_2$ Density',
    'brineViscosity': 'Brine Viscosity',
    'CO2Viscosity': 'CO$_2$ Viscosity',
    'brineResSaturation': 'Residual Brine Saturation',
    'compressibility': 'Compressibility',
    'injRate': 'Injection Rate',
    # Analytical Reservoir component
    'logResPerm': 'Log. Reservoir Permeability',
    'reservoirPorosity': 'Reservoir Porosity',
    'reservoirRadius': 'Reservoir Radius',
    'brineDensity': 'Brine Density',
    'CO2Density': 'CO$_2$ Density',
    'brineViscosity': 'Brine Viscosity',
    'CO2Viscosity': 'CO$_2$ Viscosity',
    'brineResSaturation': 'Residual Brine Saturation',
    'brineCompressibility': 'Brine Compressibility',
    'injRate': 'Injection Rate',
    'numberOfShaleLayers': 'Number of Shale Layers',
    'shaleThickness': 'Shale Thickness', # Does not include variations like
    # shale2Thickness - they will default to the actual parameter name
    'aquiferThickness': 'Aquifer Thickness',
    'reservoirThickness': 'Reservoir Thickness',
    'datumPressure': 'Datum Pressure',
    # Generic Reservoir component
    'reservoirDepth': 'Reservoir Depth', # Note: Generic reservoir uses it for base
    # of reservoir, OpenWellbore uses it for top (wellbore base)
    'reservoirThickness': 'Reservoir Thickness',
    'resTempGradient': 'Reservoir Temp. Gradient',
    'initialSalinity': 'Initial Salinity',
    # Wellbore component parameters
    # Multisegmented Wellbore component
    'logWellPerm': 'Log. Wellbore Permeability',
    'logAquPerm': 'Log. Aquifer Permeability',
    'brineDensity': 'Brine Density',
    'CO2Density': 'CO$_2$ Density',
    'brineViscosity': 'Brine Viscosity',
    'CO2Viscosity': 'CO$_2$ Viscosity',
    'aquBrineResSaturation': 'Residual Brine Saturation\nin Aquifer',
    'compressibility': 'Compressibility',
    'wellRadius': 'Well Radius',
    # Cemented Wellbore component
    'logThiefPerm': 'Log. Thief\nZone Permeability',
    'initPressure': 'Initial Pressure\nof Wellbore Base',
    'wellDepth': 'Depth to Well Base',
    'depthRatio': 'Depth Ratio',
    # Open Wellbore component
    'logReservoirTransmissivity': 'Log. Reservoir Transmissivity',
    'logAquiferTransmissivity': 'Log. Aquifer Transmissivity',
    'brineSalinity': 'Brine Salinity',
    'wellTop': 'Well Top Depth',
    'criticalPressure': 'Critical Pressure',
    'reservoirDepth': 'Reservoir Depth',
     # Generalized Flow Rate component
    'logPeakCO2Rate': 'Log. Peak\nCO$_2$ Flow Rate',
    'timePeakCO2Rate': 'Time of Peak\nCO$_2$ Flow Rate',
    'durationPeakCO2Rate': 'Duration of Peak\nCO$_2$ Flow Rate',
    'durationPeakZeroCO2Rate': 'Duration from Peak CO$_2$\nFlow Rate to Zero',
    'logInitBrineRate': 'Log. Initial\nBrine Flow Rate',
    'logFinalBrineRate': 'Log. Final\nBrine Flow Rate',
    'durationInitBrineRate': 'Duration of Initial\nBrine Flow Rate',
    'durationInitFinalBrineRate': 'Duration of Decrease from Initial\nto Final Brine Flow Rate',
    'mitigationTime': 'Time of\nLeakage Remediation',
    # Seal Horizon component
    'area': 'Cell Area',
    'thickness': 'Cell Thickness',
    'baseDepth': 'Depth to Base Seal',
    'permeability': 'Cell Initial Permeability',
    'entryPressure': 'Entry Pressure',
    'aveThickness': 'Seal Layer Average Thickness',
    'stdThickness': 'Seal Layer Thickness\nStandard Deviation',
    'minThickness': 'Seal Layer Minimum Thickness',
    'maxThickness': 'Seal Layer Maximum Thickness',
    'avePermeability': 'Seal Layer Average Permeability',
    'stdPermeability': 'Seal Layer Permeability\nStandard Deviation',
    'minPermeability': 'Seal Layer Minimum Permeability',
    'maxPermeability': 'Seal Layer Maximum Permeability',
    'heterFactor': 'Permeability Heterogeneity Factor',
    # Didn't finish Seal Horizon, there are a TON of parameters that are unlikely to be explored.
    # Left out Fault Flow and Fault Leakage components
    # Aquifer component parameters
    # Carbonate Aquifer component
    'ithresh': 'Impact Threshold',
    'rmin': 'Minimum Distance\nBetween Separate Leaks',
    'perm_var': 'Log. Permeability Variance',
    'corr_len': 'Correlation Length',
    'aniso': 'Anisotropy Factor',
    'mean_perm': 'Log. Mean Permeability',
    'hyd_grad': 'Horizontal Hydraulic Gradient',
    'calcite_ssa': 'Calcite Surface Area',
    'organic_carbon': 'Organic Carbon\nVolume Fraction',
    'benzene_kd': 'Benzene Distribution Coefficient',
    'benzene_decay': 'Benzene Decay Constant',
    'nap_kd': 'Napthalene Distribution Coefficient',
    'nap_decay': 'Napthalene Decay Constant',
    'phenol_kd': 'Phenol Distribution Coefficient',
    'phenol_decay': 'Phenol Decay Constant',
    'cl': 'Brine Salinity',
    'logf': 'Log. (1) or Linear (0)\nTransform of Plume Volume',
    'aqu_thick': 'Aquifer Thickness',
    # Deep Alluvium Aquifer component
    'logK_sand1': 'Log. Permeability\nof Layer 1',
    'logK_sand2': 'Log. Permeability\nof Layer 2',
    'logK_sand3': 'Log. Permeability\nof Layer 3',
    'logK_caprock': 'Log. Permeability\nof Caprock',
    'correlationLengthX': 'Correlation Length\nin x-Direction',
    'correlationLengthZ': 'Correlation Length\nin z-Direction',
    'sandFraction': 'Sand Volume Fraction',
    'groundwater_gradient': 'Regional Groundwater Gradient',
    'leak_depth': 'Depth of Leakage Interval',
    # FutureGen2Aquifer and FutureGen2AZMI components
    'por': 'Porosity',
    'log_permh': 'Log. Horizontal Permeability',
    'log_aniso': 'Anisotropy Ratio',
    'rel_vol_frac_calcite': 'Calcite Relative\nVolume Fraction',
    'aquifer_salinity': 'Aquifer Salinity',
    'reservoir_salinity': 'Reservoir Salinity',
    'dissolved_salt_threshold': 'Dissolved Salt Threshold',
    'dissolved_co2_threshold': 'Dissolved CO$_2$ Threshold',
    # AtmosphericROM component
    'T_amb': 'Ambient Temperature',
    'P_amb': 'Ambient Pressure',
    'V_wind': 'Wind Velocity',
    'C0_critical': 'Critical Concentration',
    'T_source': 'Released CO$_2$ Temperature',
    'x_receptor': 'x-Coordinate of Receptor',
    'y_receptor': 'y-Coordinate of Receptor',
    # Chemical Well Sealing component
    'fractureAperture': 'Fracture Aperture',
    'fractureLength': 'Fracture Length',
    'maxOverpressure': 'Maximum Overpressure',
    # Fault Leakage
    'damage_zone_perm': 'Fault Permeability',
    'damage_zone_por': 'Fault Porosity',
    'shallow_aquifer_perm': 'Shallow Aquifer Permeability',
    'shallow_aquifer_por': 'Shallow Aquifer Porosity',
    'deep_aquifer_perm': 'Deep Aquifer Permeability',
    'deep_aquifer_por': 'Deep Aquifer Porosity',
    'well_index': 'Well Index',
    'well_rate': 'Injectio Rate',
    'dip_angle': 'Fault dip angle',
    'injection_time': 'Duration of Injection',
    'geothermal_gradient': 'Geothermal Gradient',
    }

# This is used for legends in sensitivity_analysis.py
OUTPUT_DICT = {
    # Reservoir components
    'pressure': 'Pressure',
    'CO2saturation': 'CO$_2$ saturation',
    'mass_CO2_reservoir': 'CO$_2$ mass',
    # Wellbore and adapter components
    'CO2_aquifer': 'CO$_2$ leakage rate to aquifer',
    'CO2_atm': 'CO$_2$ leakage rate to atmosphere',
    'brine_aquifer': 'Brine leakage rate to aquifer',
    'brine_atm': 'Brine leakage rate to atmosphere',
    'mass_CO2_aquifer': 'CO$_2$ mass leaked to aquifer',
    'mass_brine_aquifer': 'Brine mass leaked to aquifer',
    # Fault flow and Seal horizon components
    'CO2_aquifer_total': 'Cumulative CO$_2$ leakage rate to aquifer',
    'brine_aquifer_total': 'Cumulative brine leakage rate to aquifer',
    'mass_CO2_aquifer_total': 'Cumulative CO$_2$ mass leaked to aquifer',
    'mass_brine_aquifer_total': 'Cumulative brine mass leaked to aquifer',
    # Aquifer components
    'Benzene_volume': 'Benzene plume volume',
    'Naphthalene_volume': 'Naphthalene plume volume',
    'Phenol_volume': 'Phenol plume volume',
    'As_volume': 'Arsenic plume volume',
    'Pb_volume': 'Lead plume volume',
    'Cd_volume': 'Cadmium plume volume',
    'Ba_volume': 'Barium plume volume',
    'Flux': 'CO$_2$ leakage rate to atmosphere',
    # FutureGen2 aquifer components
    'TDS_volume': 'TDS plume volume',
    'TDS_dx': 'TDS plume length',
    'TDS_dy': 'TDS plume width',
    'TDS_dz': 'TDS plume height',
    'pH_volume': 'pH plume volume',
    'pH_dx': 'pH plume length',
    'pH_dy': 'pH plume width',
    'pH_dz': 'pH plume height',
    'Pressure_volume': 'Pressure plume volume',
    'Pressure_dx': 'Pressure plume length',
    'Pressure_dy': 'Pressure plume width',
    'Pressure_dz': 'Pressure plume height',
    'Dissolved_CO2_volume': 'Dissolved CO$_2$ volume',
    'Dissolved_CO2_dx': 'Dissolved CO$_2$ length',
    'Dissolved_CO2_dy': 'Dissolved CO$_2$ width',
    'Dissolved_CO2_dz': 'Dissolved CO$_2$ height',
    'Temperature_volume': 'Temperature plume volume',
    'Temperature_dx': 'Temperature plume length',
    'Temperature_dy': 'Temperature plume width',
    'Temperature_dz': 'Temperature plume height',
    # Generic aquifer component
    'Dissolved_salt_volume': 'Dissolved salt plume volume',
    'Dissolved_salt_dr': 'Dissolved salt plume radius',
    'Dissolved_salt_dz': 'Dissolved salt plume height',
    'Dissolved_salt_mass_fraction': 'Salt mass fraction in aquifer',
    'Dissolved_CO2_dr': 'Dissolved CO$_2$ plume radius',
    'Dissolved_CO2_mass_fraction': 'CO$_2$ mass fraction in aquifer',
    # Plume stability component
    'pressure_areas': 'Pressure plume area',
    'pressure_areas_dt': 'Change in pressure plume area',
    'pressure_mobility': 'Velocity of pressure plume centroid',
    'pressure_mobility_angles': 'Direction of pressure plume centroid',
    'pressure_spreading': 'Dispersion of pressure plume',
    'pressure_spreading_angles': 'Direction of pressure plume dispersion',
    'CO2saturation_areas': 'CO$_2$ saturation plume area',
    'CO2saturation_areas_dt': 'Change in CO$_2$ saturation plume area',
    'CO2saturation_mobility': 'Velocity of CO$_2$ saturation plume centroid',
    'CO2saturation_mobility_angles': 'Direction of CO$_2$ saturation plume centroid',
    'CO2saturation_spreading': 'Dispersion of CO$_2$ saturation plume',
    'CO2saturation_spreading_angles': 'Direction of CO$_2$ saturation plume dispersion',
    }
