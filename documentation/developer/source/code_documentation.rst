******************
Code Documentation
******************

.. toctree::

Docstrings
==========

Since the integration of the component based on simulation model into the
NRAP-Open-IAM framework will not necessarily be performed by the developers
who initially developed the model, it is necessary for the model developers
to provide some level of details covering the possible use of the model.
The description of the code and logistics behind the model is provided
in comments and dosctrings written within the code. Docstrings are a
"special form" of comments. They usually occur as
the first statement in a class, method or function description.
We strongly recommend to include dosctrings in all modules written
for the NRAP-Open-IAM tool. Information provided in docstrings is used to create
the components description section in the NRAP-Open-IAM user's manual.
This means that docstrings describing the components should follow the specific format.
For examples of the desired format, the developer should consult the available
components code. Here, we consider the docstring for the CementedWellbore component
class which is the first statement one would see after opening a Python file
with the component code. This type of docstrings is referred to as a module string.
The pieces of code considered above also contain docstrings.

.. testcode::

        """
        The Cemented Wellbore component model is based on a multiphase well leakage
        model implemented in the NRAP-IAM-CS, :cite:`HARP2016150`. The model is
        built off detailed full-physics Finite Element Heat and Mass (FEHM)
        simulations, :cite:`Zyvoloski2007`. The FEHM transfer simulations are
        three-dimensional (3-D), multiphase solutions of heat and mass transfer of
        water and supercritical, liquid, and gas |CO2|. After the simulations are
        completed, the surrogate model is built based on the key input parameters
        and corresponding output parameters. The approximate (surrogate) model is
        represented by polynomials in terms of input parameters that then can be
        sampled to estimate leakage rate for wells.  Early development work can be found
        in :cite:`RN1606`.

        When using the control file interface with more than 3 shale layers,
        the *ThiefZone* keyword can be used to specify the thief zone aquifer and the
        *LeakTo* keyword can be specified to name the upper aquifer.  These values will
        default to 'aquifer1' and 'aquifer2' respectively.

        Component model input definitions:

        * **logWellPerm** [|log10| |m^2|] (-13.95 to -10.1) - logarithm of wellbore
          permeability (default -13).

        * **logThiefPerm** [|log10| |m^2|] (-13.995 to -12) - logarithm of thief zone
          permeabilty (default -12).

        * **wellRadius** [m] (0.025 to 0.25) - radius of the wellbore (default 0.05 m).

        * **initPressure** [Pa] (1.0e+5 to 5.0e+7) - initial pressure at the base of
          the wellbore (default 2.0e+7 Pa, or 20 MPa). *From linked component.*

        * **wellDepth** [m] (960 to 3200) - depth in meters from ground surface to
          top of reservoir (default 1500 m). *Linked to Stratigraphy.*

        * **depthRatio** [-] (0.3 to 0.7) - fraction of well depth to the center of
          the thief zone from the top of the reservoir (default 0.5).
          *Linked to Stratigraphy.*

        The possible outputs from the Cemented Wellbore component are leakage rates
        of |CO2| and brine to aquifer, thief zone and atmosphere. The names of the
        observations are of the form:

        * **CO2_aquifer1**, **CO2_aquifer2**, **CO2_atm** [kg/s] - |CO2| leakage rates.

        * **brine_aquifer1**, **brine_aquifer2**, **brine_atm** [kg/s] -
          for brine leakage rates.

        * **mass_CO2_aquifer1**, **mass_CO2_aquifer2** [kg] - mass of |CO2|
          leaked into aquifers.

        """

In this example of a module docstring, the first part of docstring is devoted
to the general description of the component, including references for users
interested in additional details about the model, and details covering
the use of model. The second and the most important part provides the description
of all parameters of the model, including the names, short descriptions and ranges.
The final part of dosctring contains the names of all possible observations
of the  ``simulation_model`` method. When this docstring is compiled as part of the user's guide
it looks different from above and is formatted according to the specifications
used inside the dosctring. Below is the example of how the CementedWellbore
module docstring will be compiled for the document.

----

**Cemented Wellbore Component Docstring Output**

.. automodule:: openiam.CementedWellbore

----

Since we utilize Sphinx for the compilation of the user's guide, all module
dosctrings provided with the component code have to satisfy the
reStructuredText (reST) format that Sphinx relies on. The general docstrings
conventions are described in `PEP 257 <https://www.python.org/dev/peps/pep-0257/>`_,
the Python docstrings guide.
The main difference between comments and dosctrings is that the former
explains what a given section of code is doing, while the latter describes
how a particular method can be used.
