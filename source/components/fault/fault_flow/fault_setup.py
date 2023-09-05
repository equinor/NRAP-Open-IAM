# Define default setup of the Fault Flow Component
FAULT_SETUP_DICT = {
    'ModelParams': {
        'runTitle': 'Example_setup',
        'startTime': 0.0,                # years
        'endTime': 50.0,                 # years
        'injectEndTime': 50.0,           # years; Time when injection stops
        'realizations': 5,               # -; Simulation cycles
        'timePoints': 5,                 # -; Including start point
        'timeInput': False,              # True/False; True = Time steps from File
        'useSI': 'Metric'},
    'Controls': {
        'apertureApproach': False,       # True/False; True = Get Apertures from file
        'considerNearApproach': False,   # True/False; Use near-surface model
        'profileType': 0,                # 0,1,2: >0 = alternative model
        'strikeApproach': False,         # True/False; True = Vary fault strike
        'dipApproach': False,            # True/False; True = Vary fault dip
        'pressureApproach': False,       # True/False; True = Vary aperture with pressure
        'interpolateApproach': False},   # True/False; True = Interpolate field conditions (fluid constants)
    'FaultCore': {
        'faultProbability': 100.0 ,              # %; ave. probability of existence
        'aveStrike': 90.0,                       # deg; strike of plane from North
        'spreadStrike': 0.0,                     # deg; spread in orientation
        'aveDip': 90.0,                          # deg; average dip
        'stdDip': 0.0,                           # deg; std deviation angle of dip
        'length': {'valu': 100.0, 'unit': 'm'},  # length at surface from start point
        'xStart': {'valu': 500.0, 'unit': 'm'},  # x-coordinate, fault start point
        'yStart': {'valu': 500.0, 'unit': 'm'},  # y-coordinate, left start
        'nSegments': 1},                         # number of fault divisions / plates
    'Field': {
        'aquiferDepth': {'valu': 240.0, 'unit': 'm'},         # Depth to base of deepest aquifer
        'aquiferTemperature': {'valu': 22.0, 'unit': 'oC'},   # Base of deepest aquifer
        'aquiferPressure': {'valu': 1.42E+07, 'unit': 'Pa'},  # At base of aquifer
        'injectDepth': {'valu': 1880.0, 'unit': 'm'},         # At top of injection layer
        'injectTemperature': {'valu': 95.0, 'unit': 'oC'},    # Temperature at inject. depth
        'fieldPressure': {'valu': 1.9140e+07, 'unit': 'Pa'},  # At injection depth
        'injectPressure': {'valu': 2.9290E+07, 'unit': 'Pa'}, # At depth during injection
        'finalPressure': {'valu': 1.9140e+07, 'unit': 'Pa'},  # At depth after injection
        'injectX': {'valu': 0.00, 'unit': 'm'},               # x-coord. of injection well
        'injectY': {'valu': 0.00, 'unit': 'm'}},              # y-coord. of injection well
    'Aperture': {
        'aveAperture': {'valu': 10.0, 'unit': 'mm'},    # Ave. fault aperture
        'stdDevAperture': {'valu': 0.0, 'unit': 'mm'},  # Std. deviation of aperture
        'minAperture': {'valu': 1.0e-4, 'unit': 'mm'},  # Minimum aperture
        'maxAperture': {'valu': 20.0, 'unit': 'mm'},    # Maximum aperture
        'SGR': 0.0,                                     # %; SGR (%) = shale ratio
        'stateVariable': 1.0},                          # Non-Isothermal correction
    'InjectConditions': {
        'salinity': 0.0,                                          # ppm; at injection horizon
        'CO2Density': {'valu': 673.84, 'unit': 'kg/m^3'} ,        # CO2 density
        'CO2Viscosity': {'valu': 5.5173e-05, 'unit': 'Pa*s'},     # CO2 viscosity
        'brineDensity': {'valu': 974.895, 'unit': 'kg/m^3'},      # Brine density
        'brineViscosity': {'valu': 3.0491e-04, 'unit': 'Pa*s'},   # Brine viscosity
        'CO2Solubility': {'valu': 0.035, 'unit': 'mol/kg'}},      # CO2 Solubility
    'AquiferConditions': {
        'CO2Density': {'valu': 886.44, 'unit': 'kg/m^3'},         # CO2 density
        'CO2Viscosity': {'valu': 8.8010e-05, 'unit': 'Pa*s'},     # CO2 viscosity
        'brineDensity': {'valu': 1004.10, 'unit': 'kg/m^3'},      # Brine density
        'brineViscosity': {'valu': 3.0221e-04, 'unit': 'Pa*s'}},  # Brine viscosity
    'RelativeFlowLimits': {
        'relativeModel': 'BC',                             # --; Model type; options: BC, LET
        'brineResSaturation': 0.15,                        # --; Residual wetting
        'CO2ResSaturation': 0.0,                           # --; Residual nonwetting; typically = 0.0
        'permRatio': 0.6,                                  # --; Ratio of nonwetting/wetting
        'entryPressure': {'valu': 5000.0, 'unit': 'Pa'}},  # Default threshold value
    'BrooksCoreyModel': {
        'lambda': 2.5},                  # --; Brooks Corey Model lamda parameter; used only if relativeModel = BC
    'LETModel': {
        'wetting1': 1.0,                 # -; Used only if relativeModel = LET
        'wetting2': 10.0,                # -; Used only if relativeModel = LET
        'wetting3': 1.25,                # -; Used only if relativeModel = LET
        'nonwet1': 1.05,                 # -; Used only if relativeModel = LET
        'nonwet2': 3.50,                 # -; Used only if relativeModel = LET
        'nonwet3': 1.25},                # -; Used only if relativeModel = LET
    'LETCapillaryModel': {
        'capillary1': 0.2,               # -; Used only if relativeModel = LET
        'capillary2': 2.8,               # -; Used only if relativeModel = LET
        'capillary3': 0.43,              # -; Used only if relativeModel = LET
        'maxCapillary': {'valu': 1.0e+7, 'unit': 'Pa'}}, # Pa; Used only if relativeModel = LET
    'Stress': {
        'maxHorizontal': {'valu': 3.0e+07, 'unit': 'Pa'},        # Pa; secondary at depth of top of injection
        'minHorizontal': {'valu': 2.0e+07, 'unit': 'Pa'},        # Pa; secondary at depth of top of injection
        'maxTrend': 55.0},               # deg; from North
    'FlowPlot': {
        'timeSeriesPlot': False,         # True/False; True = plot time history
        'skipOutput': True}}             # True/False; True = skip sim. print
