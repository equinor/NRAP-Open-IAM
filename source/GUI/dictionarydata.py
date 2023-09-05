"""
This file consists of all the variables that are globally accessible.
Some of the variables are constants.
"""
d = {}
componentVars = {}


# componentChoices contains the names of the components already added to the simulation and
# to which connections can be created
componentChoices = []
# componentTypeDictionary contains types of the components added to the simulation in the order they were added
componentTypeDictionary = []
# connectionsDictionary contains names of the components to which the connections were created
# by the components at the given index position
connectionsDictionary = []
# connections is updated with the names of new components upon addition of the latter to the simulation
# Contains unique names of the possible connections + Dynamic Parameters option
connections = ['Dynamic Parameters']
# connectionTypes contains names of the possible connections for the components
# of the system starting with the second added component and ending with the possible
# connections for the next component (to be added yet)
connectionTypes = []

savedDictionary = {}

# Width of some widgets
# Width of labels for distribution parameters, e.g. Value:, Mean:, etc.
DISTRIBUTION_ARG_LABEL_WIDTH = 8
# Width of text fields (entries) for distribution parameters
DISTRIBUTION_ARG_TEXTFIELD_WIDTH = 8
# Width of menu containing distribution options
DISTRIBUTION_MENU_WIDTH = 16
# Width of buttons on most of the pages
BUTTON_WIDTH = 18
# Width of entries for file path input
FILE_ENTRY_WIDTH = 40
# Width of labels containing names of parameters for all components except
# generalized flow rate component and stratigraphy
PARAMETER_LABEL_WIDTH = 31
# Width of short labels containing names of component outputs
OUTPUT_LABEL_WIDTH1 = 25
# Width of long labels containing names of component outputs
OUTPUT_LABEL_WIDTH2 = 40
# Width of longer labels containing names of component outputs
OUTPUT_LABEL_WIDTH3 = 45
# Padding added on the left and right of the output checkboxes
CB_PADX = (0, 1)
# Padding added on the left and right of the parameter frames
PARAMETER_FRAME_PADX = 5
# Width of labels containing names of parameters of generalized flow rate component
GFR_PARAMETER_LABEL_WIDTH = 37
# Width of labels containing names of parameters of stratigraphy
STRATA_PARAMETER_LABEL_WIDTH = 25
# Width of labels containing names of parameters of fault leakage component
FL_PARAMETER_LABEL_WIDTH = 35
# Width of labels containing names of parameters of hydrocarbon leakage component
HCL_PARAMETER_LABEL_WIDTH = 35


# Width of widgets on Model tab
MODEL_TAB_LARGE_LABEL_WIDTH = 30
MODEL_TAB_LABEL_WIDTH1 = 25
MODEL_TAB_LABEL_WIDTH2 = 20
MODEL_TAB_LABEL_WIDTH3 = 8
MODEL_TAB_ENTRY_WIDTH = 8
MODEL_TAB_MENU_WIDTH = 10

# Width of widgets on Add Components tab
SETUP_LABEL_WIDTH = 20
SETUP_ENTRY_WIDTH = 26
SETUP_MENU_WIDTH = 20

# Setup width of widgets on the PostProcessor_Page
POSTPROCESSOR_LABEL_WIDTH = 20
POSTPROCESSOR_SPINBOX_WIDTH = 25
POSTPROCESSOR_ENTRY_WIDTH = 40
POSTPROCESSOR_FILE_ENTRY_WIDTH = 40
POSTPROCESSOR_CHECKBOX_WIDTH = 20
POSTPROCESSOR_MENU_WIDTH = 25

# Default aquifers to which leakage is possible unless number of shale layers
# is greater than 3
aquifers = ['aquifer1', 'aquifer2']

# Types of analysis available in NRAP-Open-IAM
ANALYSIS_TYPES = ["Forward", "LHS", "Parstudy"]

# Levels of logging available in NRAP-Open-IAM
LOGGING_TYPES = ["Debug", "Info", "Warning", "Error"]

# Types of distribution available for most of the parameters
# Remove Normal and lognormal for now since they cannot be bounded by lower and upper boundaries
# and there is no equivalent of truncated lognormal distribution
#DISTRIBUTION_OPTIONS = ['Fixed Value', 'Uniform', 'Normal',
#                       'Lognormal', 'Truncated', 'Triangular', 'Discrete']
DISTRIBUTION_OPTIONS = ['Fixed Value', 'Uniform', 'Truncated',
                        'Triangular', 'Discrete']

SPEC_DISTRIBUTION_OPTIONS = ['Fixed Value', 'List', 'File Input']

# Labels for possible parameters characterizing a particular distribution
DISTRIBUTION_PARS_LABELS = {
    'Fixed Value': ['Value:'],
    'List':        ['Value(s):'],
    'Uniform':     ['Min:', 'Max:'],
    'Normal':      ['Mean:', 'SD:'],
    'Lognormal':   ['Mean:', 'SD:'],
    'Truncated':   ['Mean:', 'SD:', 'Min:', 'Max:'],
    'Triangular':  ['Min:', 'Max:', 'Mode:'],
    'Discrete':    ['Values:', 'Weights:']}

# Texts for tooltips corresponding to different distribution parameters
DISTRIBUTION_PARS_SETUPS = {
    'Value:':   ['value', 'Set value of {0}.'],
    'Value(s):': ['ordvalues', 'Set list of ordered values for {0}.'],
    'File Input': ['filename', 'Select file containing values for {0}.'],
    'Min:': ['min', 'Set minimum of {0}.'],
    'Max:': ['max', 'Set maximum of {0}.'],
    'Mean:':    ['mean', 'Set mean of {0}.'],
    'SD:': ['std', 'Set standard deviation of {0}.'],
    'Mode:':    ['mode', 'Set the most common value of {0}.'],
    'Values:':  ['values', 'Set list of values for {0}.'],
    'Weights:': ['weights', 'Set weights for each value entered in Values field.']}

# Components available in the current version of GUI for NRAP-Open-IAM
COMPONENT_TYPES = [
    'Analytical Reservoir',
    'Lookup Table Reservoir',
    'Generic Reservoir',
    'Theis Reservoir',
    # 'Simple Reservoir',
    'Multisegmented Wellbore',
    'Cemented Wellbore',
    'Cemented Wellbore (WR)',
    'Open Wellbore',
    'Generalized Flow Rate',
    'Fault Flow',
    'Fault Leakage',
    'Hydrocarbon Leakage',
    'Seal Horizon',
    'Carbonate Aquifer',
    # 'Alluvium Aquifer (LF)',
    'Deep Alluvium Aquifer',
    #'Deep Alluvium Aquifer (ML)',
    'FutureGen2 Aquifer',
    'FutureGen2 AZMI',
    'Generic Aquifer',
    'Atmospheric ROM',
    'Plume Stability',
    'Chemical Well Sealing']

# Types of processing/plotting analysis and/or capabilities available for a given simulation
# Plotting is default and available for all simulations
processingTypes = ['Plotting']
APP_SIZE = [1170, 750]
TAB_SIZE = [700, 550]

# Types of plots available for a given simulation
plotTypes = ['Time Series']

# Setup of fonts for different widgets in NRAP-Open-IAM
LABEL_FONT = ('Helvetica', 14, 'bold')
LABEL_FONT2 = ('Helvetica', 12, 'bold')
INSTRUCTIONS_FONT = ('Helvetica', 12, 'bold')
DASHBOARD_FONT = ('Helvetica', 12)
