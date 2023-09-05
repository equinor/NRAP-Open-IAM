.. highlight:: python
   :linenothreshold: 3

.. _control_file:

Control Files
=============

Control files are a method of getting user input into the NRAP-Open-IAM for setting
up a simulation scenario. Control files use a YAML format (extension *.yaml*).
Any line in the control file starting with a pound sign (#) is a comment and
is ignored by the program. The basic format of the control file is a parameter
name followed by a colon, space, and a value. For objects with several parameters
the object name is followed by a colon and the underlying parameters are listed
on the consecutive lines tabbed in. For example, consider this partial file

.. code-block:: python
   :lineno-start: 1

    #-------------------------------------------------
    # NRAP-Open-IAM Control File example
    #-------------------------------------------------
    ModelParams:
        EndTime: 50
        TimeStep: 1.0
        Analysis: forward
        Components: [SimpleReservoir1,
                     CementedWellbore1]
        OutputDirectory: ../../output/output_ex1a_{datetime}
        Logging: Debug
    Stratigraphy:
        numberOfShaleLayers:
            vary: False
            value: 3
        # Thickness is in meters
        shale1Thickness:
            min: 500.0
            max: 550.0
            value: 525.0
        shale2Thickness:
            min: 450.0
            max: 500.0
            value: 475.0
        shale3Thickness:
            vary: False
            value: 11.2
        aquifer1Thickness:
            vary: False
            value: 22.4
        aquifer2Thickness:
            vary: False
            value: 19.2
        reservoirThickness:
            vary: False
            value: 51.2

Here, the first three lines are comments that are not read by the code.
The fourth line defines the keyword ``ModelParams`` which describes parameters of
the system model. The subsequent lines contain parameters of ``ModelParams``.
A ``ModelParams`` section is required in all NRAP-Open-IAM control files. The ``EndTime``
keyword defines the ending time for the simulation in years (50 years in
the example). The ``TimeStep`` parameter defines the length of a time step
(1 year in the example). The type of analysis being run is a ``forward`` (deterministic) simulation.
Other possible options for ``Analysis`` parameter are ``lhs`` for Latin
Hypercube Sampling analysis and ``parstudy`` for a parameter study. The ``Components`` parameter
is a required entry that contains a list of component model names that are
defined later in the file. The component list will always begin with a square
bracket '[' followed by each of the component names that make up the system
separated by a comma ',' and closed by a square bracket ']'. The names
of the components are listed in the order they are supposed to be run. The next keyword
``OutputDirectory`` defines a directory for the output to be written into.
The output directory can be appended with a keyword ``{datetime}``.
When a simulation is run, the ``{datetime}`` keyword will be replaced with the date
and time of the simulation run. Note that the ``{datetime}`` keyword is optional: if it is
omitted subsequent runs of the simulation will overwrite past results.
That is, if there is a need to keep all results from re-running an NRAP-Open-IAM case,
the ``{datetime}`` keyword will easily facilitate this; if re-running an
NRAP-Open-IAM case should overwrite previous results, the ``{datetime}`` keyword should be
omitted. In the output folder the NRAP-Open-IAM places a copy of the input file, all outputs
from the component models written to text files, and *.png* images for all
graphics. The last keyword ``Logging`` defines what level of logging
information is written out to the logging files. Options for ``Logging`` levels
are ``Debug``, ``Info``, ``Warning``, and ``Error``. ``Info`` is the default level (if no
``Logging`` keyword is given the logging will be set to ``Info``) and will give
you a valuable information about when parameters go outside of permitted
ranges and when there are problems in the program. Debug is a good option if
you have problems with your IAM model and want more information to explore the
causes.

The next keyword section of the file is the required ``Stratigraphy`` section. In
this section any model parameters related to the stratigraphy of the |CO2| storage site is
defined. Any parameters for the stratigraphy are defined here with either
a deterministic value or a range to vary over. A fixed value of any
given parameter can be specified with the ``vary: False`` and ``value: ###``
specification shown here or simply ``parameterName: ###``.
The ``min`` and ``max`` specification gives a range for the parameter to vary over
if an analysis is done over multiple realizations. See the :ref:`stratigraphy_component`
section of this document for a list of all available parameters.

The next sections of the input file defines every component model
in the component model list specified earlier in the control file.
The first component listed in the example is SimpleReservoir1 defined as follows

.. code-block:: python
   :lineno-start: 37

    #-------------------------------------------------
    # SimpleReservoir1 is a user defined name for component;
    # the type SimpleReservoir is the ROM model name
    #-------------------------------------------------
    SimpleReservoir1:
        Type: SimpleReservoir
        Parameters:
            injRate: 0.1
        Outputs: [pressure,
                  CO2saturation]

This section of the file defines a ``SimpleReservoir`` Component model named
*SimpleReservoir1* to be part of the system model. The name *SimpleReservoir1* can
be replaced with any other name defined by user, but will not be a part of the system
model unless it is an element of the components list described in the previous section
``ModelParams``. The ``Type`` is a keyword that defines the component model
to be used and must match up with one of the component models currently
available in the NRAP-Open-IAM.
The ``Parameters`` section defines parameters of the component model.
Description of parameters available for the user to specify can be found in
the :ref:`components_description` chapter of the current documentation. The component model
parameters are specified in the same fashion as the ``Stratigraphy`` parameters.
The ``Outputs`` specifies the observations of the component model
that will be output from the simulation.
Please refer to the :ref:`components_description` chapter of this document to see which
parameters and outputs are available for user specification in the control file.

Generally, dynamic (time-varying) input to component models comes from the
output of other connected component models (e.g., the pressure and saturation
as an input to a wellbore leakage model comes from a reservoir model). In some instances
there may be a need to study a component model without the other attached component
models feeding the input. In this case dynamic input can be specified with
the ``DynamicParameters`` keyword. Under the ``DynamicParameters`` section each
input name is specified followed by a list of values (enclosed in square brackets [])
of the same length as the number of time points (a value for each time point, including
an initial value). See files *ControlFile_ex7a.yaml* and *ControlFile_ex7b.yaml*
for example of control files utilizing dynamic input for some components.

The next section of the input file is similar to the previous section and defines
the next component model *CementedWellbore1*.

.. code-block:: python
   :lineno-start: 47

    #-------------------------------------------------
    CementedWellbore1:
        Type: CementedWellbore
        Connection: SimpleReservoir1
        Number: 4
        Locations:
            coordx: [100, 540]
            coordy: [100, 630]
        RandomLocDomain:
            xmin: 150
            xmax: 250
            ymin: 200
            ymax: 300
        Parameters:
            logWellPerm:
                min: -14.0
                max: -12.0
                value: -13.0
        Outputs: [CO2_aquifer1,
                  CO2_aquifer2,
                  CO2_atm,
                  brine_aquifer1,
                  brine_aquifer2]

In this part of the example, ``CementedWellbore`` type component model is
specified. There are four wellbores of this type being added with ``Number: 4``:
two of the locations are given in the ``Locations`` part and other two are
generated randomly within the domain specified in ``RandomLocDomain`` part.

Unknown wellbore locations can be generated by specifying more wellbores
(with ``Number: ``) than the number of known wellbore locations.
To control the location of the random well placement, a ``RandomLocDomain``
section need to be used as

.. code-block:: python
   :lineno-start: 55

        RandomLocDomain:
            xmin: 150
            xmax: 250
            ymin: 200
            ymax: 300

This specification will limit the x-coordinate of random wells to be between 150
and 250, and the y-coordinate to be between 200 and 300. Sampling will be
from a uniform distribution on the domain defined by ``xmin`` (``ymin``)
and ``xmax`` (``ymax``).

All coordinate systems are assumed to have units of meters and are defined
by the reservoir component used. Known wells will be placed first;
after all known well coordinates are used wells will be
placed within the random wells domain.

Known wellbore coordinates are entered as a comma separated list. There must be
a comma between each coordinate. Random wellbores generated in an
area when more wells are specified than number of known coordinates. After
completing the model parameters proceed to the Stratigraphy tab.

For the SimpleReservoir component model the default injection
location is at [0, 0]. For the Lookup Table based reservoir component
the wellbore locations should fall within the domain of the reservoir simulations.

See *ControlFile_ex3.yaml* for additional example using random well
placement and *ControlFile_ex4a.yaml* for example using only known well locations.

The last section of the input file is used to specify a graphical output

.. code-block:: python
   :lineno-start: 70

    #-------------------------------------------------
    # Plot setup part of the control file
    #-------------------------------------------------
    Plots:
        CO2_Leakage1:
            TimeSeries: [CO2_aquifer1]
            subplot:
                ncols: 2
                use: True
        CO2_Leakage2:
            TimeSeries: [CO2_aquifer2]
            subplot:
                ncols: 2
                use: True
        Pressure_plot:
            TimeSeries: [pressure]
            subplot:
                ncols: 2
                use: True
                SimpleReservoir1_000.pressure: 'Pressure at well #1'
                SimpleReservoir1_001.pressure: 'Pressure at well #2'
                SimpleReservoir1_002.pressure: 'Pressure at well #3'
                SimpleReservoir1_003.pressure: 'Pressure at well #4'
            Title: Reservoir Pressure at Wellbore Location

Here, three plots are being requested. The firsts two plots will illustrate
the |CO2| leakage to the shallow aquifer and the thief zone aquifer;
the third plot will illustrate the pressures in the reservoir for the
four wellbore locations specified earlier in the control file.
*CO2_Leakage1*, *CO2_Leakage2* and *Pressure_plot* are the user defined names
of the three plots to be created: these will also be used as the filenames of the figures
saved in the output directory. ``TimeSeries`` is a keyword that instructs
the program to plot the observation data as a time series plot. The values to be
plotted (**CO2_aquifer1**, **CO2_aquifer2** and **pressure** above) have to be defined
in the control file as outputs from one of the specified component models.
Each plot will have a title corresponding to the values plotted. A user
defined title can be specified with the ``Title`` keyword (as illustrated for
the *Pressure_plot*) in the given plot section. For each aquifer
the |CO2| leakage rates for all wells will be plotted
on the same figure but on different subplot. If each observation is to be plotted
on a separate subplot, the ``subplot`` keyword with ``use`` set to ``True``
must be specified, as illustrated in the example setup. Additionally,
the ``ncols`` keyword (under ``subplot`` section) can be used to set the number
of subplot columns to use. The number of rows is controlled by the number
of different values (observations) to plot over the number of columns.
Each subplot will be given a (default) title of the variable plotted unless specified by user.
The default title names can be replaced with the user defined ones by
using the full observation name as a key and the desired title as the value
under ``subplot`` section as shown in the setup of *Pressure_plot*.

The example file described here can be found in the *examples/Control_Files* directory with
the filename *ControlFile_ex1a.yaml*. To run this example, open a command
prompt in the *examples/Control_Files* directory and run the command::

    python ../../source/openiam/openiam_cf.py --file ControlFile_ex1a.yaml

Note: use \\ on Windows and / on Mac and Linux.

Other example control files can be found in the same directory. They can be
run by replacing the file name in the above command with the user specified one.
