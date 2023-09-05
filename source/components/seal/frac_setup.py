#------------------------------------------------------------------------------
FRAC_SETUP_DICT = {
    'Connectivity': {
        'geometric': 1.000},                # Connectivity and roughness factor for site
    'Controls': {
        'randomFracApproach': True,        # True/False; True = Generate random fractures
        'userFracApproach': False,          # True/False; True = Input fractures from file
        'apertureLengthApproach': False,    # True/False; True = Aperture-length correl.
        'pressureApproach': False,          # True/False; True = Pressure/aperture concept
        'plotFracApproach': False},         # True/False; True = Generate fracture plots
    'RandomFracs': {
        'Density': {
            'aveDensity': {'valu': 2.0e-04, 'unit': 'per_m^2'},    # per_mˆ2
            'minDensity': {'valu': 1.0e-04, 'unit': 'per_m^2'},    # per_mˆ2
            'maxDensity': {'valu': 3.0e-04, 'unit': 'per_m^2'}},   # per_mˆ2
        'Orientation': {
            'aveTrend': 30.0,          # degrees from North;
            'spreadTrend': 20.0},      # approximately 2-sigma spread
        'Length': {
            'function': 'lognormal',                      # lognormal or powerlaw
            'aveLength': {'valu': 50.0, 'unit': 'm'},     # m; for lognormal
            'stdDevLength': {'valu': 20.0, 'unit': 'm'},  # m; for lognormal
            'expLength': {'valu': -0.8 , 'unit': 'm'},    # m; for powerlaw where x is in m
            'minLength': {'valu': 1.0e-1, 'unit': 'm'},   # m
            'maxLength': {'valu':  500.0 , 'unit': 'm'}}, # m
        'VaryAperture': {
            'aveAperture': {'valu': 5.0e-4, 'unit': 'mm'},     # mm; also reference aperture
            'stdDevAperture': {'valu': 1.0e-4, 'unit': 'mm'},  # mm
            'minAperture': {'valu': 1.0e-6, 'unit': 'mm'},     # mm
            'maxAperture': {'valu': 5.0e-3, 'unit': 'mm'}},    # mm
        'CorrelateAperture': {
            'alpha': 0.6,                               # exponent for correlation
            'beta': {'valu': 0.00010, 'unit': 'mm'}},   # mm; aperture of 1 meter fracture
        'Threshold': {
            'refEntry': {'valu': 5.0e+03, 'unit': 'Pa'}},  # Pa; reference to mean aperture
        'PressureCorrection': {
            'resAperture': {'valu': 0.005, 'unit': 'mm'},     # mm; residual hydraulic aperture
            'wideAperture': {'valu': 0.05 , 'unit': 'mm'},    # mm; maximum hydraulic aperture
            'stressLimit': {'valu': 2.5e+07, 'unit': 'Pa'},   # Pa; stress limit - nonlinear response
            'theta': 2.5}},                                   # stiffness history factor
    'InputFracs': {
        'Filename': 'user_fracture_input',  # user input directory
        'Threshold': {
            'refAperture': {'valu': 1.0e-2, 'unit': 'mm'},    # mm
            'refPressure': {'valu': 1.0e+06, 'unit': 'Pa'}}}, # Pa
    'RockMatrix': {
        'Permeability': {
            'aveMatrix': {'valu': 1.0e-08, 'unit': 'mD'},     # mD (microdarcys)
            'stdDevMatrix': {'valu': 1.0e-10, 'unit': 'mD'},  # mD (microdarcys)
            'minMatrix': {'valu': 1.0e-09, 'unit': 'mD'},     # mD (microdarcys)
            'maxMatrix': {'valu': 3.0e-06, 'unit': 'mD'}},    # mD (microdarcys)
        'Threshold': {
            'matrixPerm': {'valu': 1.0e-08, 'unit': 'mD'},     # mD (microdarcys)
            'matrixPress': {'valu': 1.0e+06, 'unit': 'Pa'}}}}  # Pa
