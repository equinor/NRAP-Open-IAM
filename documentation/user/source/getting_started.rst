.. toctree::

***************
Getting Started
***************

The NRAP-Open-IAM has several ways for a user to build and
run simulations, including graphical user interface (GUI), text based control files
and python scripts. The simplest way to build and run simulations for NRAP-Open-IAM
is the GUI. This guide will primarily focus on using the GUI to interact
with the NRAP-Open-IAM. To launch the GUI open a command prompt in the *source/GUI*
directory and type::

    python NRAP_OPENIAM.py

.. include:: gui_operation.rst

..
   # TODO Reference and cite use of Matk here
..
   # TODO Include control file text somewhere?

.. include:: control_file.rst

Output
======

Output is written to the folder specified in the model definition with the
``Output directory``. If the path to the output directory is not absolute
(i.e., does not containt the drive letter) it is assumed
to start from the NRAP-Open-IAM root folder containing the tool distribution.
For each component model of the system ``Outputs`` can be specified. When an
output is specified for a forward model the values of that output are written
to a file (*output_name.txt*) in the ``Output Directory``. For a stochastic model
(LHS or parametric study analysis) the outputs are included in the results file
(*LHS_results.txt* or *parstudy_results.txt*) as well as a file for a statistical summary
of the input parameters and output observations (*LHS_statistics.txt* or
*parstudy_statistics.txt*). A copy of the input control file is also written
to the ``Output Directory`` folder. Through the GUI, this input file can be loaded
back in to rerun the simulation.

After a simulation is run, the post processing section of the GUI can be used
to generate plots of the results and run sensitivity analysis. When the post processing
page is first opened it asks for a folder specification. This folder is the output folder
of the results you want to analyze. Navigate to that output folder using the **Browse** button.

Plotting
--------
User can access post-processing capabilities of GUI by clicking on the
**Post Processing** button on the main page. The **Post Processor** window will
appear. After a folder containing the results of the simulation is selected
different options for the simulation results appear that are already
set for plotting. There are several types of plots that can be created
depending on what type of simulation was run and what components were specified.
The simplest plot is a ``Time Series`` plot where the output is plotted against
time, multiple realizations will be plotted as separate lines. Specify a title
to give the plot and a file name along with making a selection of what output to plot.
Pressing the plot button will generate the plot, it will be saved with
the filename given to the output directory for the results. If a simulation
was run with multiple realizations a ``Time Series Stats`` plot or a ``Time Series
and Stats`` plot will be options in the ``Plot Type`` menu. A ``Time Series Stats``
plot will shade the quadrants of the results along with plotting the mean and
median results, but will not plot the individual realizations. A ``Time Series
and Stats`` plot will overlay the realizations on the stats plots of the
shaded quadrants. If AtmosphericROM component was included in the simulation,
map-view plots of the plume for a single realization or probabilistic ensemble
can be generated.

Sensitivity Analysis
--------------------
If the simulation results are from a LHS simulation the ``Processing`` menu
will have options for several types of sensitivity analysis. Note that while
a sensitivity analysis can be run on simulations with a small number of
realizations, the results will most likely be inaccurate. If the sensitivity
coefficients do not sum to one, or if they vary largely through time,
the number of realizations might need to be increased. Generally, 500 to 1000
realizations are needed for a sensitivity analysis. However, this might
change depending on the complexity of the simulation. Each type of
sensitivity analysis will produce plots and/or text file output in the output directory.

``Correlation Coefficients`` option produces a plot matrix of either ``Pearson`` or
``Spearman`` correlation coefficients. Any system model observation can be
excluded from the analysis if needed, although no exclusions need to be made.

``Sensitivity Coefficients`` option calculates the sensitivity coefficients
for each selected output to all inputs. Selecting multiple outputs will run
the sensitivity coefficient calculation multiple times. The capture point
is the index for point in time at which the sensitivity coefficient are
to be calculated. The analysis produces a bar chart.

``Multiple Sensitivity Coefficients`` option calculates the impact of input parameters
on multiple outputs. Multiple outputs should be selected here. The capture point
is the index for point in time at which the sensitivity coefficient are
to be calculated. The analysis will produce a bar chart.

``Time Series Sensitivity`` option will produce a line graph illustrating how the impact
from input parameters changes over time with respect to an output value.
Selecting multiple output values will run the analysis multiple times.
The capture point determines the time at which the sensitivity coefficients
are compared and then ordered based on the comparison.

Analysis Options in Control File
================================

The NRAP-Open-IAM uses the Model Analysis ToolKit (MATK) :cite:`MATK` for
the basis of its probabilistic framework. More information about MATK can be found here:
`http://dharp.github.io/matk/ <http://dharp.github.io/matk/>`_. The MATK code repository can
be found here: `https://github.com/dharp/matk <https://github.com/dharp/matk>`_.

Parameter input and output interactions can be explored using the *Analysis* section
of the control file. Correlation coefficients can be calculated using the
``CorrelationCoeff`` keyword. Parameter sensitivity coefficients for any output
simulation value can be calculated using a Random-Balanced-Design Fourier
Amplitude Sensitivity Test (RBD-Fast) technique. The control file keywords
``SensitivityCoeff``, ``MultiSensitivities``, ``TimeSeriesSensitivity`` can be
used to access different sensitivity coefficient outputs. See *ControlFile_ex8.yaml*
for details on using the analysis section.

The Sensitivity Analysis is done with the ``SALib`` package :cite:`Herman2017`.
For more information on the RBD-Fast technique see :cite:`TISSOT2012205` and
:cite:`PLISCHKE2010354`. While not accessible through the control files, a
scripting interface is provided to a Sobol sensitivity analysis :cite:`SOBOL20093009`.

This section will be expanded in the future.

.. include:: cfi_visualization.rst

.. include:: units.rst
