# Changelog
All notable changes to this project are documented in this file.

**Note (for all current and future maintainers of the code):** Whenever the version
of the tool is updated and changes are recorded in the current changelog,
please perform the following actions:
  - update version of the tool in the following files:
    - source/openiam/\__init\__.py
    - documentation/user/source/conf.py
    - documentation/developer/source/conf.py
    - documentation/qaqc/source/conf.py
    - documentation/website/source/conf.py
    - source/GUI/Dashboard.py
    - replace image of front page of GUI with the one stating new version
    of the tool in the documentation.
  - recompile user and developer guides, and QAQC documents;
  - add tag to the master branch corresponding to the current change log record;
  - replace old version of the user and developer guides with the new recompiled
  files in the current folder.
  - do not forget to commit the changes.


## [a2.7.2 - 2023-08-25]
### Added
  - CFI: Added new section for workflows; examples illustrating AoR workflow
  setup
  - Auxiliary component: Well Depth Risk Configurer
  - Jupyter notebooks illustrating setup and run of control file examples
  including workflows, and TTFD and AoR types of analysis and visualization
  - Output of the current NRAP-Open-IAM version and runtime information
  to a text file for control and GUI files.
### Updated
  - Auxiliary Location Generator component: generation of z-coordinates,
  in addition to x- and y-coordinates
  - Logging handlers bug fix
  - Moved CFI related files to a separate folder
  - Documentation updates for webpage

## [a2.7.1 - 2023-06-30]
### Added
  - Hydrocarbon Leakage component: class and CFI
  - Cemented Wellbore (WR) component: GUI and an example illustrating setup
  - Fault Leakage component: CFI and GUI; CFI and GUI examples illustrating setup
  - CFI: Stratigraphy column type of plot; examples illlustrating setup
  - CFI: Gridded radial metric type of plot for Generic Aquifer component;
  examples illustrating setup
  - CFI: Gridded metric type of plot for Seal Horizon and Fault Flow components;
  examples illustrating setup
  - Lookup Table Reservoir component: option for the use of lookup tables
  in *.h5 format; related tests to the test suite
  - Plume Stability component: option for the component metrics to be
  calculated for 3d data; option for the use of lookup tables in *.h5
  format; related tests to the test suite
  - Auxiliary functionality: well data classifier, leak scorer from well integrity
  study
  - Auxiliary components: Parameter Setup, Pressure Based Risk Configurer,
  Data Based Risk Configurer, Monitoring Scheduler, Monitoring Tool
  - Functionality to save gridded observations in formats different from *.npy
  - Documentation:
   - NRAP-Open-IAM workflows illustration figures
   - Section on conceptual model overview
### Updated
  - CFI and GUI for Open Wellbore component: option to consider critical
  pressure for leakage estimates
  - Plotting functionality related to sensitivity analysis in CFI and GUI;
  new CFI examples illustrating setup
  - Stratigraphy, time series, Total time to the first detection (TTFD), Area of Review
  workflows plots
  - Documentation:
    - Description of existing and added visualization capabilities

## [a2.7.0 - 2023-02-15]
### Added
  - Generic Reservoir component: class, CFI and GUI; script and CFI examples illustrating
  setup
  - Theis Reservoir component: class and CFI
  - Cemented Wellbore (WR) component: class and CFI (updated version of the previous
  cemented wellbore component); script and CFI examples illustrating setup
  - Fault Leakage component: class
  - Script example illustrating coupling of NRAP-Open-IAM and SOSAT tool
  - Additional options for setup of time series and AoR plots in CFI
  - Total time to the first detection (TTFD) workflow plots in CFI; examples
  illustrating setup
  - CFI: Stratigraphy setup through input file and its visualization through
  new stratigraphy plot
  - Script and CFI examples illustrating new visualization capabilities
  - GUI: options to save simulation results files for individual or combined outputs
### Updated
  - Setup of AmtPlumeSingle and AtmPlumeEnsemble plots
  - Documentation:
    - Description of recently added components
    - Description of new CFI examples
    - Description of existing and added visualization capabilities

## [a2.6.1 - 2022-09-30]
### Added
  - GUI: time points setup through manual and file inputs
  - CFI: new options for setting up locations for all placeable components
  (reservoir, leakage pathways, impact)
  - Area of Review (AoR) calculation script and CFI examples; AoR plot option
  - CFI: options to save simulation results files for individual or combined outputs
### Updated
  - Seal Horizon component class: added option to sample thickness and permeability
  parameters inside component; for CFI added option of different controls inherited
  from standalone Seal Flux code
  - Fault Flow component class: added option to sample strike, dip and aperture
  parameters; for CFI added option of different controls inherited
  from standalone FaultFlo code
  - Generic Aquifer component class: fixed issue with results for zero
  CO2 leaked mass
  - Restructured code related to CFI: the methods are split among several files
  according to the functionality they cover

## [a2.6.0 - 2022-08-12]
### Added
  - Generic Aquifer: CFI and GUI; examples illustrating applications
  - Permeability Sampler (for Seal Horizon component)
  - Thickness Sampler (for Seal Horizon component)
  - Added instructions for non-Anaconda installation
  - GUI: dynamic input through relative path
  - CFI: 3d interpolation setup for LUTR and wellbore components connected to it
  - CFI: setup of additional metric outputs for LUTR component
### Updated
  - Seal Horizon component (due to changes in the standalone Seal Flux code):
  class, CFI and GUI
  - Fault Flow Component (due to changes in the standalone Fault Flo code):
  class, CFI and GUI
  - Requirements for dynamic input file formatting: for several locations
  the spatial data changes row wise, temporal data changes column wise.
  For a single location, temporal data can change both row and column wise.

## [a2.5.0 - 2022-03-10]
### Added
  - Generic Aquifer component: class
  - Added new method collect_gridded_observations_as_time_series to SystemModel
  class to simplify export of gridded observations both for deterministic
  and stochastic simulations
### Changed
  - Moved the code to support Python 3.9
### Updated
  - Installation instructions for Python 3.9
  - Documentation for Fault Flow component
  - Lookup Table Reservoir component:
    - added check that the indexing of reservoir data start with 1 not with 0;
    - added check that the required data for pressure and saturation are present;
    - added check that the names of required observations pressure and
    CO2saturation are spelled correctly;
    - added check that the indices are integers and return an error message
    if they are not

## [a2.4.0 - 2021-11-11]
### Added
  - Chemical Well Sealing component: class, CFI, and GUI
### Updated
  - Saving GUI setup file: removed saving incorrect information not pertinent
  to the setup
  - GUI Postprocessing tab: fixed clearing of checked boxes next to the output
  names left after load of previous results
  - Reservoir and multisegmented wellbore components: update of the code; fixed
  calculation of total mass of CO2 injected into reservoir; added extra accumulator
  for CO2 volume
  - GUI setup example files: replaced simple reservoir with analytical
  reservoir
  - Updated SampleSet's method collect_observations_as_time_series to allow
  to collect observations corresponding to the user provided indices

## [a2.3.2 - 2021-08-05]
### Added
 - GUI for seal horizon and fault flow components
### Updated
 - GUI for simple and analytical reservoir components: added interface
 for injection well location
 - GUI for all reservoir components: added option to enter locations
 of points of observations which can differ from wellbore locations

## [a2.3.1 - 2021-06-15]
### Updated
 - ROMs for FutureGen 2 aquifer components to avoid predicted dz larger
 than aquifer thickness
 - Changed observation names of several aquifer components for consistency
 among all aquifer components
 - Updated script example used for the manuscript
 - Updated interface of discrete distribution for consistency between parameters
 of LUTR and other components. Three options are available: (Values, Weights),
 (values, weights), and (discrete_vals)

## [a2.3.0 - 2021-05-24]
### Added
 - Analytical reservoir component: class, CFI, and GUI
 - Seal Horizon component: script and CFI
 - Fault Flow component: script and CFI
 - QA tests for stratigraphy and Lookup Table Reservoir components
### Updated
 - Save button interface in GUI: allows user to choose location to save GUI setup file
 - Parstudy setup: fixed bug in MATK and made updates in GUI
 - CFI: Added an option to enter time points as array (for varying time step option)
 or input file

## [a2.2.0 - 2021-02-12]
### Added
 - Index parameter is now required to be present in the first column in the file
 parameters_and_filenames.csv for Lookup Table Reservoir component setup.
 - Added 3d interpolation capability for Lookup Table Reservoir; added illustrative
 example
 - Started QAQC documentation for the tool
### Updated
 - CFI and GUI interface of Lookup Table Reservoir and Plume Stability Components
 - Documentation for FutureGen2 Aquifer and AZMI components
 - Kimberlina data set is updated to adhere to new requirements to lookup table
 data sets

## [a2.1.0 - 2020-05-22]
### Added
 - CFI and GUI for Plume Stability component
 - CFI and GUI: wellbore locations interface for all wellbore components
 replaced single setup on the Model tab
 - Parameters check is added for all components available in GUI
 - Method add_par added to the Lookup Table Reservoir component class
### Updated
 - GUI: Stratigraphy tab is updated
 - GUI: Lookup Table Reservoir tab is updated
 - Open wellbore component model; updated examples and CFI, GUI for the components
 - Cemented wellbore component model
 - Test suite is updated with new tests

## [a2.0.0 - 2020-02-02]
### Added
 - FutureGen2 Aquifer and AZMI components are added
 - Several script examples are added: (a) script that generates output for DREAM;
 (b) script for Risk-Based Area of Review
 - CFI and GUI for FutureGen2 Aquifer and AZMI components; control file and
 GUI example files are added
 - DeepAlluviumAquifer and AlluviumAquifer components utilizing ML techniques
 are added
 - Tests and scripts are added to illustrate work of new aquifer (Kimberlina
 based and FutureGen2) components
 - Added description of dynamic (temporal) inputs for some wellbore and aquifer
 components to the tool documentation
 - Preparation work for CFI development of Plume Stability Analysis
 component: updated component module documentation and renamed constructor arguments
### Updated
 - Stratigraphy component is updated to provide depth to the bottom
 of each shale and aquifer layer
 - Installation instructions for Windows both for standard and Anaconda
 Python distribution
 - Prep work for external alpha 2.0 release.
 - Removed lookup table data folders not used in the examples from the repository

## [a1.2.0 - 2019-11-08]
### Added
 - Location Generator component class; added example to the scripts folder
 - Added Jupyter Notebook examples covering simple and more complex examples
 of OpenIAM application
 - Plume Stability Analysis component updates: added example of use and test
### Updated
 - For proper work with Location Generator component, LUT Reservoir component
 was updated to allow providing locations as kwargs of the LUT reservoir component
 - Plume Stability Analysis component updates: added interpolation and possibility
 to analyze more than one input variable in a single component
 - Added missing description of control file examples and saved GUI setup files
 to the user manual
 - Updated postprocessing capabilities of GUI: fixed some bugs and added
 keywords for title and file names

## [a1.1.2 - 2019-08-23]
### Updated
 - GUI update for LUT reservoir component: handling of deterministic
   and discrete parameters values including signature index sampling

## [a1.1.1 - 2019-08-08]
### Added
 - Generalized Flow Rate Component: graphical user interface (GUI)
 - Plume Stability Analysis Component class
### Updated
 - Correct handling of more than 3 shale layers in CFI and GUI
 - GUI code refactoring: code for each component tab has been moved
 to a separate module
 - GUI fixes: handling varying number of dynamic parameters for different components

## [a1.1.0 - 2019-07-22]
### Added
 - Generalized Flow Rate Component: script and control file interface (CFI)
### Updated
 - Dynamic parameters can be provided through direct input and data files in CFI
 - Multiple groups of wellbore can be handled through CFI
 - Multiple aquifers can be handled through CFI
 - LUT Reservoir Component: added possibility of signature index sampling

## [b1.0.2 - 2019-05-24]
### Updated
 - Updated main page, initial tabs (Model, Stratigraphy and Add components),
 and all components tabs
 - Fixed issue with loading the saved simulation on top of previously loaded
 simulation
 - Fix inconsistency with triangular and lognormal parameter distribution;
 removed truncated distribution

## [b1.0.1 - 2019-04-03]
### Updated
 - Prep work for external beta release
 - Fix issue with saving different GUI *.OpenIAM files

## [a0.6.3 - 2019-03-11]
### Updated
 - Updated GUI for scrolling for low-res and large text screens

## [a0.6.2 - 2019-02-13]
### Updated
 - LUT GUI fix
 - Documentation updating

## [a0.6.1 - 2019-01-25]
### Updated
 - GUI fixes
 - Documentation Updates

## [a0.6.0 - 2019-01-04]
### Added
 - GUI
 - Post Processing capability for GUI
### Updated
 - Moved to python 3

## [a0.5.0 - 2018-07-30]
### Added
 - Deep Alluvium Aquifer ROM
 - Plume Stability Calculation
### Fixed/Updated
 - Modified Reservoir LUTs for optional input data

## [a0.4.0] - 2018-04-20
### Added
 - Atmospheric Dispersion ROM
 - Progress bar for ensemble simulations

### Fixed/Updated
 - Reservoir LUT gridded observations
 - Cemented Wellbore parameter limits
 - Multi-segmented wellbore pressure checking
 - Renamed example files for clarity/consistency
 - Updated logging messages
 - Updated MCMC model update example
 - Added temporal input checking to Carbonate Aquifer component

## [a0.3.0] - 2018-03-01
### Added
 - Sensitivity Analysis UI and visualization
 - Dynamic kwarg
 - Gridded Parameters (Needs UI)
 - Gridded Observations
 - Plotting updates including
   - Statistics plotting for large realization runs
   - User Specified Titles
   - Subplotting Options
 - Advanced MCMC example
 - Matk Updates
   - Added sample statistics
   - Added statistics file output
   - Added RBD-Fast sensitivity analysis

### Fixed/Updated
 - Mass of CO2 in aquifers reported from wellbore models
 - Statistics output from Matk functionality
 - Clarify points on install process and lookup table files in User's Guide
 - Removed dependency on make utility for mac and linux install
 - Modified examples based on User Feedback

## [a0.2.0] - 2018-01-25
### Added
 - Reservoir Lookup Table component added
 - Stratigraphy Component added
   - This moves stratigraphy specification to new section in control files.
 - Run time output from control file runs
 - Analysis log files created for restarting options
 - Parameters and Observations written to CSV file
 - Parameter and Observation statistics (min, max, mean, etc.) written to CSV file

### Updated
 - Updated MATK
   - Handles discrete parameters (LHS)
   - Sobol Sensitivity Analysis
   - Matplotlib display check fix
 - Component model method add_obs restructured

## [a0.1.0] - 2017-12-31
### Added
 - Initial OpenIAM testing version
 - Probabilistic Framework for integrated assessment modeling, including
   - forward modeling
   - latin hypercube sampling
   - parameter studys
 - Component models included:
   - Simple Reservoir Model (semi-analytical model)
   - Multi-segmented Wellbore Model (semi-analytical model)
   - Cemented Wellbore Model
   - Open Wellbore Model
   - Carbonate Aquifer Model
 - Control-file based User Interface
 - User's Guide
 - Setup Script
 - Test Suite
 - Examples
