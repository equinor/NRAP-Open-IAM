# Define default setup of the Seal Horizon component
# Do not delete!!!
# We need this file separate from ControlFile_seal.yaml since the default setup
# for Seal Flux model is different from the default setup for Seal Horizon

SEAL_SETUP_DICT = {
    'ModelParams': {
        'runTitle': 'Example_setup',
        'useSI': 'Metric',                          # Metric/English for output
        'startTime': {'valu': 0.0, 'unit': 'yrs'},  # years
        'endTime': {'valu': 50.0, 'unit': 'yrs'},   # years
        'timePoints': 5,                            # including start point
        'realizations': 5,                          # simulation cycles
        'timeInput': False,                         # True/False; True = time steps from file
        'totalInject': {'valu': 5.0e+05, 'unit': 'tonne'}},  # total amount of CO2 injected
    'OutputDirectory': 'seal_output',
    'Description': {
        'numCells': 72,
        'area': {'valu': 1.0e+4, 'unit': 'm^2'},             # m2
        'salinity': {'valu': 15000.0, 'unit': 'ppm'},        # ppm
        'aveTemperature': {'valu': 50.0, 'unit': 'oC'},      # oC, temperature at depth
        'aveBaseDepth': {'valu': 1100.0, 'unit': 'm' },      # Seal base depth
        'aveBasePressure': {'valu': 3.30e+7, 'unit': 'Pa'},  # During injection
        'staticDepth': {'valu': 1000.0, 'unit': 'm'},        # m; reference depth, positive below grade
        'staticPressure': {'valu': 1.0e+7, 'unit': 'Pa'}},   # Pa; reference pressure
    'Grid': {
        'gridApproach': False,           # True/False; True = cells use grid
        'gridRows': 6,                   # y-axis divisions for grid
        'gridColumns': 12,               # x-axis divisions for grid
        'cellHeight': {'valu': 100.0, 'unit': 'm'},  # m; ave. ea. cell dimension N-S (y-axis)
        'cellWidth': {'valu': 100.0, 'unit': 'm'}},  # m; ave. ea. cell dimension E-W (x-axis)
    'Conditions': {
        'CO2Density': {'valu': 597.8, 'unit': 'kg/m^3'},        # kg/m^3; ignored if interpolateApproach = True
        'CO2Viscosity': {'valu': 4.452e-5, 'unit': 'Pa*s'},     # Pa*100.0s; ignored if interpolateApproach = True
        'brineDensity': {'valu': 1004.00, 'unit': 'kg/m^3'},    # kg/m^3; ignored if interpolateApproach = True
        'brineViscosity': {'valu': 5.634e-4, 'unit': 'Pa*s'},   # Pa*s; ignored if interpolateApproach = True
        'CO2Solubility': {'valu': 0.035, 'unit': 'mol/kg'}},    # mol/kg; ignored if interpolateApproach = True
    'FileInput': {
        'depthApproach': False,          # True/False; True = depths to cell base from file
        'layoutApproach': False,         # True/False; True = coordinates from file
        'areaApproach': False,           # True/False: True = area from file
        'upperBoundApproach': False,     # True/False; True = upper-bound pressure from file
        'activeCellsApproach': False,    # True/False; True = status from file
        'entryApproach': False,          # True/False; True = entry pressure from file
        'permApproach': False,           # True/False; True = permeability from file
        'thicknessApproach': False},     # True/False; True = input from File
    'Controls': {
        'fractureApproach': False,       # True/False; True = Use fractures
        'correlateApproach': False,      # True/False; True = correlate entry press. w. perm.,
        'interpolateApproach': False,    # True/False; True = interpolate fluid constants
        'initializeApproach': False},    # True/False; True = initialize model from file
    'Thickness': {
        'aveThickness': {'valu': 100.0 , 'unit': 'm'},      # m
        'stdDevThickness': {'valu': 0.0, 'unit': 'm'},      # m
        'minThickness': {'valu': 75.0, 'unit': 'm'},        # m
        'maxThickness': {'valu': 125.0, 'unit': 'm'}},      # m
    'Permeability': {
        'avePermeability': {'valu': 2.5e-4, 'unit': 'mD'},  # microDarcy
        'stdDevPermeability': {'valu': 0.0, 'unit': 'mD'},  # microDarcy
        'minPermeability': {'valu': 1.0e-6, 'unit': 'mD'},  # microDarcy
        'maxPermeability': {'valu': 1.0e-3, 'unit': 'mD'},  # microDarcy
        'heterogeneityApproach': False,  # True/False; True = include heterogeneity
        'heterFactor': 0.5},
    'RelativeFlowLimits': {
        'brineResSaturation': 0.15,      # residual wetting
        'CO2ResSaturation': 0.0,         # residual nonwetting; typically = 0.0
        'relativeModel': 'LET',          # model type; options: BC, LET
        'permRatio': 0.6},               # ratio of nonwetting/wetting
    'CapillaryPressure': {
        'entryPressure': {'valu': 5000.0, 'unit': 'Pa'}}, # Pa; default threshold value
    'BrooksCoreyModel': {
        'lambda':  2.5},                 # used only if relativeModel = BC
    'LETModel': {
        'wetting1': 1.0,                 # used only if relativeModel = LET
        'wetting2': 10.0,                # used only if relativeModel = LET
        'wetting3': 1.25,                # used only if relativeModel = LET
        'nonwet1': 1.05,                 # used only if relativeModel = LET
        'nonwet2': 3.50,                 # used only if relativeModel = LET
        'nonwet3': 1.25},                # used only if relativeModel = LET
    'LETCapillaryModel': {
        'capillary1': 0.2,               # used only if relativeModel = LET
        'capillary2': 2.8,               # used only if relativeModel = LET
        'capillary3': 0.43,              # used only if relativeModel = LET
        'maxCapillary': {'valu': 1.0e+07, 'unit': 'Pa'}}, # Pa; used only if relativeModel = LET
    'TimeModel': {
        'influenceModel': 0,             # options: 0/1/2
        'influence': 1.0,                # initial influence factor 0 to 1
        'rateEffect': 0.1,               # time model -> rate of effect
        'totalEffect': 0.1,              # time model -> total effect as ratio
        'reactivity': 8.0,               # value 1 to 10
        'clayContent': 60.0,             # percent
        'carbonateContent': 8.0,         # percent
        'clayType': 'smectite'},         # clay type; types: smectite/illite/chlorite
    'SealPlots': {
        'permeabilityPlot': False,       # True/False; True = plot permeabilities
        'timeSeriesPlot': False,         # True/False; True = plot time history
        'CO2ContourPlot': False,         # True/False; True = contour CO2 release
        'maxDrawTime': {'valu': 50.0, 'unit': 'yrs'}}}
