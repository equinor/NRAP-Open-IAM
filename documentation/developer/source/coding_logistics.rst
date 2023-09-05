****************
Coding Logistics
****************

.. toctree::

Component model
===============

The NRAP-Open-IAM is built on the open-source concept to promote transparency and
to give advanced users the capability to contribute to the code. The main
framework for the NRAP-Open-IAM is written in Python 3. Python provides
the cross-platform capabilities, an extensive number of libraries for
data handling, analysis, and visualization, the flexibility to interact
with other languages that component models may use (Fortran, C++, etc.).
No coding style is enforced but conforming to
`PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_ is recommended.
The exceptions are that the line lengths should not exceed 120 characters,
and spaces (no tabs) must be used for indentations.

Underpinning the NRAP-Open-IAM is the MATK Python package (Model Analysis ToolKit)
:cite:`MATK` which provides a framework for the system model components.
The NRAP-Open-IAM is designed to be flexible with regard to how component models
are added so only the desired components of the system model need
to be included for a specific use case.

We start with describing the two main terms developer has to know
for the successful integration of a component model into the NRAP-Open-IAM framework.
The first term is *simulation model*
(sometimes referred to as a *component model*) that stands behind every building block (component)
of the system. For example, a reservoir and a collection of wellbores can constitute
one system within the NRAP-Open-IAM. The behavior of the reservoir and wellbores
within the system is described by simulation models developed
specifically for a given component (e.g., reservoir or wellbore).
Thus, simulation model is a function which when provided with the component
input parameters calculates the corresponding output (describes the component behavior).
The component code must be licensable as open source. This guide
is intended to instruct simulation code developers on the integration
of their simulation models into the NRAP-Open-IAM framework regardless
of the language used to develop the simulation model.

The second term is *component model* (``ComponentModel``) *class*.
The NRAP-Open-IAM code relies on using several classes providing the framework
for all use cases. Among these classes are ``SystemModel`` and ``ComponentModel`` classes.
Every component (class) currently available within the NRAP-Open-IAM distribution is written
as a class which inherits its main features, attributes and methods from ``ComponentModel`` class.

The ``ComponentModel`` class is written to help integrate different kinds
of simulation models into the NRAP-Open-IAM framework. It serves as a wrapper in the case
when the simulation model is developed in a language different from Python.
It is used as an interface between the NRAP-Open-IAM and simulation models already
written in Python. Two instance methods of the ``ComponentModel`` class that
we consider first are ``__init__`` and ``simulation_model``. These two methods are
examples of methods that are usually redefined within a given component. ::

    class ComponentModel(object):
        """
        NRAP-Open-IAM ComponentModel class.

        Class represents basic framework to construct component models, handle different
        types of parameters and observations.

        """
        def __init__(self, name, parent, model='',
                     model_args=[], model_kwargs={}, workdir=None):
            """
            Constructor method of ComponentModel class.

            :param name: name of component model
            :type name: str

            :param parent: the SystemModel object that the component model belongs to
            :type parent: SystemModel object

            :param model: Python function whose first argument is a dictionary
                of parameters. The function returns dictionary of model outputs.
            :type model: function or method that is called when the component is run

            :param model_args: additional optional parameters of the component
                model; by default, model_args is empty list []
            :type model_args: [float]

            :param model_kwargs: additional optional keyword arguments of
                the component model; by default, model_kwargs is empty dictionary {}.
            :type model_kwargs: dict

            :param workdir: name of directory to use for model runs (serial run case)
            :type workdir: str

            :returns: object -- ComponentModel object
            """
            self._parent = parent
            self.name = name
            self.model = model
            self.model_args = model_args
            self.model_kwargs = model_kwargs
            self.workdir = workdir

            # Parameters
            self.default_pars = OrderedDict()
            self.deterministic_pars = OrderedDict()
            self.pars = OrderedDict()
            self.composite_pars = OrderedDict()
            self.gridded_pars = OrderedDict()
            self.parlinked_pars = OrderedDict()
            self.obslinked_pars = OrderedDict()

            # Keyword arguments
            self.obs_linked_kwargs = OrderedDict()
            self.grid_obs_linked_kwargs = OrderedDict()
            self.dynamic_kwargs = OrderedDict()
            self.collection_linked_kwargs = OrderedDict()

            # Observations
            self.linkobs = OrderedDict()
            self.accumulators = OrderedDict()
            self.obs = OrderedDict()
            self.grid_obs = OrderedDict()
            self.local_obs = OrderedDict()

            # Set the working directory index for parallel runs
            self.workdir_index = 0

            # Setup how often the model method should be run
            self.run_frequency = 2

        # Other methods of the class follow the __init__ method

The ``__init__`` method of ``ComponentModel`` class defines possible
attributes of the ``ComponentModel`` class object. We will discuss them later
in the guide. The ``model`` attribute of the class object is defined by the *model*
argument provided to the ``__init__`` method. A user's component model class
derived from the ``ComponentModel`` class
usually utilizes ``__init__`` method to define the setup of the component.
It can be used to set the attributes of the component needed for the proper
functioning of its simulation model: e.g., input parameters of the component and
their bounds, properties of component observations, directory containing
any extra files. We will use snapshots of the available NRAP-Open-IAM code
to show the possible implementation of the ``__init__`` method.


Method ``__init__``
===================

In this section we consider an example of a typical ``__init__`` method in a component. ::

    class Component1(ComponentModel):

        def __init__(self, name, parent, attr1, attr2, **kwargs):
            # Set up keyword arguments of the 'model' method provided by the system model
            model_kwargs = {'time_point': 365.25, 'time_step': 365.25}
            super().__init__(
                name, parent, model=self.simulation_model, model_kwargs=model_kwargs)

            # Set default parameters of the component model
            self.add_default_par('par1', value=3.0)
            self.add_default_par('par2', value=2.0)
            self.add_default_par('par3', value=-7.0)
            self.add_default_par('par4', value=0.0)

            # Define dictionary of parameters boundaries
            self.pars_bounds = dict()
            self.pars_bounds['par1'] = [1.0, 5.5]
            self.pars_bounds['par2'] = [0.0, 10.0]
            self.pars_bounds['par3'] = [-100.0, 100.0]
            self.pars_bounds['par4'] = [-100.0, 100.0]

            # Define dictionary of temporal data limits
            # Boundaries of temporal inputs are defined by the simulation model
            self.temp_data_bounds = dict()
            self.temp_data_bounds['temp_input1'] = ['Temporal input 1', 1., 5.]
            self.temp_data_bounds['temp_input2'] = ['Temporal input 2', 1.5, 4.5]

            # Define accumulators and their initial values
            self.add_accumulator('accumulator1', sim=0.0)
            self.add_accumulator('accumulator2', sim=1.0)

            # Define additional component model attributes
            self.additional_attr1 = attr1
            self.additional_attr2 = attr2

            # Indicate how often the component should be run
            self.run_frequency = 2     # 2 is default

            # Define the conditional attribute if it is provided in kwargs
            if 'cond_attr' in kwargs:
                self.conditional_attr = kwargs['cond_attr']

Next we consider the code line by line. The first line initiates
the ``Component1`` class and specifies that it is derived from the ``ComponentModel``
class (i.e., inherits all methods and attributes from the ``ComponentModel`` class).
The second (non-blank) line shows that the constructor method ``__init__``
has several arguments: *name*, *parent*, *attr1*, *attr2*, and possibly
some extra keyword arguments whose names are not provided. ::

            model_kwargs = {'time_point': 365.25, 'time_step': 365.25}
            super().__init__(
                name, parent, model=self.simulation_model, model_kwargs=model_kwargs)

The next line indicates that the simulation model of the component is
a time-dependent function since component requires time point and time step
as its keyword arguments. Method ``super`` calls the constructor method
``__init__`` of the base class ``ComponentModel``. It provides a shortcut
to include a base class's methods without having to know the base class type
or name. The ``super`` method uses the *model_kwargs* variable as arguments
of the ``ComponentModel`` class's ``__init__`` method.
Analysis of the arguments of the ``__init__`` method shows that for this component
the ``model`` attribute of the class object (``self.model``) is defined as a method
since ``model`` argument of the ``__init__`` method has to provide a
*simulation model* of the component, i.e. model (method/function) that
transforms components input parameters into its output. ::

            self.add_default_par('par1', value=3.0)
            self.add_default_par('par2', value=2.0)
            self.add_default_par('par3', value=-7.0)
            self.add_default_par('par4', value=0.0)

The next four lines show that the ``Component1`` class object has four parameters
with names *par1*, *par2*, *par3* and *par4*. The parameters are assigned
default values of 3.0, 2.0, -7.0, and 0.0, respectively. The number
of possible model parameters is not limited by the NRAP-Open-IAM framework.
The purpose of this section of code in the ``__init__`` method is to ensure
that all parameters of the simulation model corresponding to the component
are defined even when not all parameters values are defined by users at runtime.
The values provided here are arbitrary. ::

            self.pars_bounds = dict()
            self.pars_bounds['par1'] = [1.0, 5.5]
            self.pars_bounds['par2'] = [0.0, 10.0]
            self.pars_bounds['par3'] = [-100.0, 100.0]
            self.pars_bounds['par4'] = [-100.0, 100.0]

The four lines above define the dictionary attribute ``pars_bounds`` containing
the upper and lower boundaries for each parameter of the simulation model.
This attribute is used in the method ``check_input_parameters``
of the ``ComponentModel`` class which can be redefined within the derived component
class to overwrite the default parameter bounds checks. The method ``check_input_parameters``
is called before the start of each simulation to check whether the provided input parameters
satisfy the defined boundaries. An implementation of the default
``check_input_parameters`` method is discussed later.

In many cases, the simulation model accepts time-varying inputs which have
to (possibly) satisfy some model limitations. The check for these limitations
should be implemented in the ``simulation_model`` method since these inputs change in time:
this way the limits will be rechecked at each time step as opposed to only at
instantiation of the component model object.
Similar to input parameters, the number of possible temporal inputs
of the simulation model is not limited by the NRAP-Open-IAM framework. ::

            self.temp_data_bounds = dict()
            self.temp_data_bounds['temp_input1'] = ['Temporal input 1', 1.0, 5.0]
            self.temp_data_bounds['temp_input2'] = ['Temporal input 2', 1.5, 4.5]

Each entry of the dictionary ``temp_data_bounds`` is a list of three elements
(right sides of the last two expressions above): the first is an explanatory name
of the temporal input, the last two are the lower and upper
boundaries of the input. An example of the possible implementation
of the temporal inputs check method is discussed later in the guide. ::

            self.add_accumulator('accumulator1', sim=0.0)
            self.add_accumulator('accumulator2', sim=1.0)

The two lines above allow the creation of observations wihin the model (sometimes
strictly internal to the component model) whose current value depends
on observation values the component model produced at the previous time step(s).
The accumulators are usually used in the ``simulation_model`` method for the storage
of internally calculated values needed for the ``simulation_model`` method at the next
time step (e.g., total mass calculated from flow rates at each time step). ::

            self.additional_attr1 = attr1
            self.additional_attr2 = attr2

In some cases the setup of the component requires the definition of additional
attributes specific for the given component (e.g., attributes that address
the setup of the component) that do not change during simulation and
can be utilized in the component model method or by other components.
The names ``additional_attr1`` and ``additional_attr2``
are arbitrary and can be defined to represent the intended purpose.

The following line serves as a possible example of the use of additional attributes. ::

            self.run_frequency = 2

Attribute ``run_frequency`` is a flag variable that informs the system model
of the frequency at which the component simulation model should be called. The possible values
are 2, 1, and 0. A value of 2 means that the model is called for each time point
supplied by the system, thus, the value of 2 is assigned, by default, to every
``ComponentClass`` instance and to all derived class instances.
Consequently, the line demonstrated above (``self.run_frequency = 2``)
is not required if this is an intended behavior. If the behavior (frequency of
calling the simulation model) should be different, the line should be changed to ::

            self.run_frequency = 1

or ::

            self.run_frequency = 0

A value of 1 means that the simulation model should be run only for the first
time point. A value of 0 means that the simulation model of a particular component
should not be run at all. This is useful in the situations when the component
does not have a ``simulation_model`` method, and serves as a container for parameters rather
than a source of outputs, or when the ``simulation_model`` method should not be run after
a particular time point. In the latter case some other system component
should control when the model method of a given component should be turned off. ::

            if 'cond_attr' in kwargs:
                self.conditional_attr = kwargs['cond_attr']

Some of the ``Component1`` object attributes may not be defined for all instances.
In this case the attribute can be assigned a value, for example, only if
a particular argument of the ``__init__`` method was provided.
The conditional arguments can control and define some features of
the simulation model defined by the developer. The names ``conditional_attr`` and *cond_attr*
are arbitrary and for illustration purposes only. They can be defined
to represent the intended purpose of the attribute.


Method ``simulation_model``
===========================

The ``simulation_model`` method is an instance method of the component class that either
calls the simulation model (typical if the simulation model is an external
code not written in Python), or is the simulation model itself (typical
if the simulation model is written in Python and can be implemented directly
in the component class). To avoid confusion the method is named differently
(self.simulation_model) from the instance attribute (self.model) containing reference
to the method. The ``simulation_model`` method must accept a dictionary
of input parameter values keyed by parameter names as the first argument
and return a dictionary of model results keyed by distinct observation names.
It is assumed that the arguments of the ``simulation_model`` method can be split
into time-constant and time-varying arguments.
All constant parameters defined by the user of the component model
must be passed to the ``simulation_model`` method in the dictionary described above.
The time-varying and other types of arguments (not defined by user)
can be passed as keyword arguments. The header of the ``simulation_model`` method
can include the default values of these arguments for situations
when they are not provided (see first line below; e.g., ``temp_input1=2.0``).
Sample code of a ``simulation_model`` method is provided below.
Note that this is just one of the possible implementations of the ``simulation_model``
method. The comments along with code sections provide suggestions of
``simulation_model`` method implementation. For additional examples please refer
to the source code of existing NRAP-Open-IAM components, e.g.,
*source/openiam/multisegmented_wellbore_component.py*.
The source code of the example component can be found in the folder
*examples/scripts/iam_example_component.py*. ::

    def simulation_model(self, p, temp_input1=2.0, temp_input2=3.0,
                       time_point=365.25, time_step=365.25):
        """
        :param p: input parameters of Component1 model
        :type p: dict

        :param temp_in1: the first of the two varying in time inputs of simulation_model
            method with default value of 2.0
        :type temp_in1: float

        :param temp_in2: the second of the two varying in time inputs of simulation_model
            method with default value of 2.0
        :type temp_in2: float

        :param time_point: time point in days at which the model output is
            to be calculated; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param time_step: difference between the current and previous
            time points in days; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        """
        # Obtain the default values of the parameters from dictionary
        # of default parameters
        actual_p = dict([(k,v.value) for k,v in self.default_pars.items()])

        # Update default values of parameters with the provided ones
        actual_p.update(p)

        # For the initial time point 0.0 the model should be able to return
        # the initial values of the component observations
        if time_point == 0.0:
            # Define initial values of the model observations
            # The values chosen are arbitrary and the number of possible
            # observations of model is not NRAP-Open-IAM framework limited
            out = {}
            out['obs1'] = 0.0
            out['obs2'] = 0.0
            out['obs3'] = 0.0

            # Exit method
            return out

        # Check whether the temporal inputs satisfy the model requirements
        # and/or assumptions if needed
        # The instance method check_temporal_input can use the attribute
        # temp_data_bounds defined in the __init__ method
        assumptions_satisfied = self.check_temporal_input(time, temp_input1, temp_input2)

        # The next steps depend on a particular implementation of the simulation_model method
        if assumptions_satisfied:
            # Calculate output of the component using parameters and
            # temporal keyword arguments. The signature of the model_function
            # is not defined by the NRAP-Open-IAM framework
            output = model_function(p, temp_input1, temp_input2, time_point, time_step)

            # Assign values to the component accumulators
            # acc_fun1 and acc_fun2 are replacement names for some actions performed
            # on the variable output in order to obtain accumulators values
            self.accumulators['accumulator1'].sim = acc_fun1(output)
            self.accumulators['accumulator2'].sim = acc_fun2(output)

            # Assign model observations
            # f1, f2, f3 are function names replacing some actions performed
            # on the variable output in order to obtain observations values
            out['obs1'] = f1(output)
            out['obs2'] = f2(output)
            out['obs3'] = f3(output)

            # The next type of statement is required: the model method should
            # return a dictionary with keys corresponding to the names of
            # all possible Component1 observations
            return out


Additional coding is required for implementing a ``simulation_model`` method serving
as a wrapper for the models developed in different programming languages.
Simple examples of wrappers for the models implemented in Fortran are presented
in the file *iam_simple_models.py* located in the *examples/scripts* folder of the NRAP-Open-IAM
distribution and described in the next section.


Method ``simulation_model`` as a wrapper
========================================

In some situations the developer of the component might have a code written
in language different from Python (e.g., Fortran), and rewriting the code in Python
might not be worth the efforts. In this section, we present an approach which is used for some
components in NRAP-Open-IAM and describes integration of the Fortran code. We start
with a simple example of a Fortran routine involving different types
of input parameters and results (model outputs). ::

    subroutine quad_eq_fun(a, b, c, x, N, x1, x2, y, flag) bind(C, name='quad_eq_fun')
    implicit none
    real*8, intent(in) :: a
    real*8, intent(in) :: b
    real*8, intent(in) :: c
    real*8 :: D
    real*8, intent(out) :: x1
    real*8, intent(out) :: x2
    integer, intent(in) :: N
    integer, intent(out) :: flag
    integer :: i
    real*8, parameter :: epsil = 1d-20
    real*8, dimension(1:N), intent(in) :: x
    real*8, dimension(1:N), intent(out) :: y

    D = b*b - 4.0*a*c

    if (abs(D) > epsil) then
        if (D .GT. 0.0) then
            x1 = (-b + sqrt(D))/(2.0*a)
            x2 = (-b - sqrt(D))/(2.0*a)
            flag = 1                 ! indicates two real distinct roots
        else
            x1 = -b/(2.0*a)          ! returns real part of the roots
            x2 = sqrt(-D)/(2.0*a)    ! returns complex part of the roots
            flag = 3
        endif
    else
        x1 = -b/(2.0*a)
        x2 = x1
        flag = 2
    endif

    do i = 1, N
        ! calculates parabola y-values for the entered parameters a, b, c and x-values
        y(i) = a*x(i)*x(i) + b*x(i) + c
    end do

    return
    end subroutine quad_eq_fun

First, for the Fortran function to be used in the Python code it should be
compiled into a library with the Fortran compilers. Compilation of the
Fortran code should be tested for all supported platforms.
We assume that the users working on Mac or Linux platforms should be able
to easily and properly compile the code into corresponding libraries.
For Windows users the compiled libraries will be provided. For illustration
purposes we assume that the example Fortran routine is located in file
*iam_quad_eq_fortran_model.f90*. In order to compile the above code in Windows
one can use the gfortran compiler installed and run
in a Cygwin environment as ::

    $ gfortran -c iam_quad_eq_fortran_model.f90
    $ gfortran -shared -o quad_eq_fun.dll iam_quad_eq_fortran_model.o

Option ``-c`` of the first step directs gfortran to compile the Fortran file to an object file,
rather than producing a standalone executable. This flag should be used if
the program source code consists of several files. The object files produced
by this command can later be linked together into a complete program.

The compilation of the code on Mac or Linux differs only in the second line.
In Mac the second step is ::

    $ gfortran -dynamiclib -o quad_eq_fun.dylib iam_quad_eq_fortran_model.o

For Linux the steps look like ::

    $ gfortran -fpic -c iam_quad_eq_fortran_model.f90
    $ gfortran -shared -fpic -o quad_eq_fun.so iam_quad_eq_fortran_model.o

The resulting library files for Windows, Mac and Linux will be *quad_eq_fun.dll*,
*quad_eq_fun.dylib*, and *quad_eq_fun.so*, respectively.

The following code demonstrates calling the compiled code within a ``simulation_model`` method.
The method returns the *y*-coordinates of the points on a parabola
|y = ax^2 + bx + c| also calculated inside library for defined *x*-values
and coefficients *a*, *b*, *c*. ::

    def simulation_model(p, input_array_x=None):
        '''
        Return three outputs based on the provided input parameters.

        :param p: dictionary of parameters passed to the function
        :type p: dict

        :param input_array_x: input array
        :type input_array_x: numpy.array

        :returns dictionary containing output variable and array
        '''
        # Obtain the default values of the parameters from dictionary
        # of default parameters
        actual_p = dict([(k,v.value) for k,v in self.default_pars.items()])

        # Update default values of parameters with the provided ones
        actual_p.update(p)

        a = actual_p['a']
        b = actual_p['b']
        c = actual_p['c']

        if input_array_x is None:
            input_array_x = np.linspace(-10.0, 10.0, 1000)

        # Determine the size of input array
        N = np.size(input_array_x)

        # Setup library and needed function names
        if platform == "linux" or platform == "linux2":
            # linux
             library = "quad_eq_fun.so"
        elif platform == "darwin":
            # OS X
            library = "quad_eq_fun.dylib"
        elif platform == "win32":
            # Windows...
            library = "quad_eq_fun.dll"
        functionName = "quad_eq_fun"

        # Load DLL
        external_lib = ctypes.cdll.LoadLibrary(os.path.join(os.getcwd(),library))

        # Get needed function as attribute of the library
        function = getattr(external_lib, functionName)

        # Define c classes to be used for inputs and outputs of the Fortran function
        INT = ctypes.c_int
        DOUBLE = ctypes.c_double
        NPointsArrayType = DOUBLE*N

        # Set argument types for values and pointers
        # The order should coincide with the order of arguments in the Fortran function
        function.argtypes = [ctypes.POINTER(DOUBLE),  # type for argument a
                             ctypes.POINTER(DOUBLE),  # type for argument b
                             ctypes.POINTER(DOUBLE),  # type for argument c
                             ctypes.POINTER(DOUBLE),  # type for array argument x
                             ctypes.POINTER(INT),     # type for integer argument N
                             ctypes.POINTER(DOUBLE),  # type for output x1
                             ctypes.POINTER(DOUBLE),  # type for output x2
                             ctypes.POINTER(DOUBLE),  # type for array output y
                             ctypes.POINTER(INT)]     # type for integer output flag
        function.restype = None

        # Define values of the input parameters that will be passed
        # to the Fortran function
        fun_arg_a = DOUBLE(a)
        fun_arg_b = DOUBLE(b)
        fun_arg_c = DOUBLE(c)
        fun_arg_N = INT(N)
        fun_array_x = NPointsArrayType(*input_array_x)

        out_x1 = DOUBLE()                   # initialize output variable x1
        out_x2 = DOUBLE()                   # initialize output variable x2
        out_flag = INT()                    # initialize output variable flag
        out_array_y = NPointsArrayType()    # initialize output array

        function(fun_arg_a, fun_arg_b, fun_arg_c, fun_array_x, fun_arg_N,
                 out_x1, out_x2, out_array_y, out_flag)

        # Create output dictionary
        out = dict()
        # The present simulation model is supposed to return sum and absolute
        # value of difference of the roots and y-values of the function.
        # The following calculations depend on the value of flag
        # Check flag values
        if out_flag.value==1:          # two real distinct roots
            # Extract value from the float type of output
            out['root_sum'] = abs(out_x1.value + out_x2.value)
            out['root_diff'] = abs(out_x1.value - out_x2.value)
        elif out_flag.value==2:        # single real root
            out['root_sum'] = abs(2*out_x1.value)
            out['root_diff'] = 0.0
        else:                          # two complex roots
            # Sum of roots is a doubled real part
            out['root_sum'] = abs(2*out_x1.value)
            out['root_diff'] = abs(2*out_x2.value)

        out['y'] = out_array_y[0:N]   # extract values from array type of output

        # Return output dictionary
        return out


Method ``check_input_parameters``
=================================

The purpose of the ``check_input_parameters`` method is to ensure that the parameters supplied
to the ``simulation_model`` method satisfy the model assumptions and/or limitations. Below
we provide the default ``check_input_parameters`` method
of the ``ComponentModel`` class. ::

    def check_input_parameters(self, p):
        """
        Check whether input parameters satisfy the specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        # Import of logging package can be done at the module level
        # above the definition of ComponentModel class
        import logging
        logging.debug(
            'Input parameters of component {name} are {p}.'.format(
                name=self.name, p=p))

        if hasattr(self, 'pars_bounds'):
            for key, val in p.items():
                if key in self.pars_bounds:
                    if ((val < self.pars_bounds[key])or(val > self.pars_bounds[key])):
                        logging.warning(
                            'Parameter {key} is out of boundaries.'.format(key=key))
                else:
                    logging.warning(('Parameter {key} is not recognized as '+
                                     'component {name} input parameter.').format(
                                         key=key, name=self.name))
        else:
            logging.debug(('Component {name} does not define boundaries of '+
                           'its model parameters.').format(name=self.name))

Note that there is no ``print`` function used in this example. NRAP-Open-IAM utilizes
the python package ``logging`` for the functionality needed to provide the user
with different kinds of messages related to the performance and setup of the system model.
The single argument of the ``check_input_parameters`` method is a dictionary of
simulation model parameters *p*, the same dictionary provided to the ``simulation_model`` method.
Since parameters do not vary in time the method is called only once per simulation.
We recommend to use a separate method for time-varying model inputs,
with possible (among others) name ``check_temporal_inputs``. The method
``check_temporal_inputs`` should be setup to be called for each time step
within the ``simulation_model`` method of the component. The default
``check_input_parameters`` method assumes that the parameters boundaries
are defined in the ``pars_bounds`` attribute of the ``ComponentModel`` class object.
If the attribute ``pars_bounds`` is not defined by the developer, a message
will be printed to the log to indicate this.

The method can be redefined within the class derived from the ``ComponentModel``
class. The following code snapshot provides an example of the method
modified according to the needs of ``LookupTableReservoir`` component class. ::

    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        debug_msgs = ['Checking input parameters...',
                      'Input parameters {}'.format(p)]
        for ind in [0, 1]:
            logging.debug(debug_msgs[ind])

        if not self.linked_to_intr_family:
            err_msg = ''.join(['Application cannot proceed further. ',
                               'Lookup Table Reservoir component is ',
                               'not linked to any interpolator family.'])
            logging.error(err_msg)
            raise LinkError(err_msg)

        # Save interpolators names
        if self.intr_names is None:
            self.intr_names = {}
        for intr_nm, intr_obj in self._parent.interpolators[self.intr_family].items():
            self.intr_names[intr_obj.index] = intr_nm

        # For lookup table reservoir component we need to make sure that the
        # signature created with input parameters coincide with the signature
        # of one of the interpolators used by component

        # Create signature based on default parameter values
        param_signature = {k: v.value for k, v in self.default_pars.items()}

        # Check whether there are any parameters not belonging to the signature
        for key in p:
            if key not in param_signature:
                msg = ''.join([
                    'Parameter {key} not recognized as ',
                    'a LookupTableReservoir input parameter.']).format(key=key)
                logging.warning(msg)

        # Update default signature with values of input parameters p
        param_signature.update(p)

        # Extract index from updated dictionary
        index = int(param_signature.pop('index'))
        if index != -2:
            if index not in self.intr_names:
                err_msg = ''.join([
                    'Value {} of index parameter does not correspond ',
                    'to any of the linked interpolators.']).format(index)
                logging.error(err_msg)
                raise ValueError(err_msg)
        else:
            # Check for the same signature among all connected interpolators
            signature_found = False
            for interpr in self._parent.interpolators[self.intr_family].values():
                # Compare signature of interpolator with the provided input parameters
                if interpr.signature == param_signature:
                    signature_found = True
                    break

            if not signature_found:
                err_msg = ''.join([
                    'Signature of input parameters do not coincide with ',
                    'signatures of connected interpolators {}.']).format(param_signature)
                logging.error(err_msg)
                raise ValueError(err_msg)

First, the method checks whether the component is linked
to the interpolators, needed for the proper execution of the ``simulation_model`` method,
and only after that checks whether the input parameters satisfy the model requirements.
The next example of the ``check_input_parameters`` method taken from
the ``MultisegmentedWellbore`` component class is close to the default one
but represents a modified version of the original method and accounts
for the varying number of possible (and similar) model parameters. ::

    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        debug_msg = 'Input parameters of component {} are {}'.format(self.name, p)
        logging.debug(debug_msg)

        for key, val in p.items():
            warn_msg = ''.join([
                'Parameter {} of MultisegmentedWellbore component {} ',
                'is out of boundaries.']).format(key, self.name)
            if (key[0:5] == 'shale' and key[-9:] == 'Thickness'):
                if (val < self.pars_bounds['shaleThickness'][0]) or (
                        val > self.pars_bounds['shaleThickness'][1]):
                    logging.warning(warn_msg)
                continue
            if key[0:7] == 'aquifer' and key[-9:] == 'Thickness':
                if (val < self.pars_bounds['aquiferThickness'][0]) or (
                        val > self.pars_bounds['aquiferThickness'][1]):
                    logging.warning(warn_msg)
                continue
            if key[0:7] == 'logWell' and key[-4:] == 'Perm':
                if (val < self.pars_bounds['logWellPerm'][0]) or (
                        val > self.pars_bounds['logWellPerm'][1]):
                    logging.warning(warn_msg)
                continue
            if key[0:6] == 'logAqu' and key[-4:] == 'Perm':
                if ((val < self.pars_bounds['logAquPerm'][0]) or
                        (val > self.pars_bounds['logAquPerm'][1])):
                    logging.warning(warn_msg)
                continue
            if key[0:3] == 'aqu' and key[-18:] == 'BrineResSaturation':
                if ((val < self.pars_bounds['aquBrineResSaturation'][0]) or
                        (val > self.pars_bounds['aquBrineResSaturation'][1])):
                    logging.warning(warn_msg)
                continue
            if key in self.pars_bounds:
                if ((val < self.pars_bounds[key][0]) or
                        (val > self.pars_bounds[key][1])):
                    logging.warning(warn_msg)


Method ``reset``
================

The method ``reset`` is called once before each simulation. Its purpose is
to reset all the needed attributes, parameters, accumulators and observations
of the particular component to some initial state. ::

    def reset(self):
        """
        Reset parameters, observations and accumulators.

        Parameters, observations and accumulators are reset to their
        initial/default values at the beginning of each new simulation.
        """
        pass

In the default method of the ``ComponentModel`` class there is only one statement:
``pass``. The current default method serves as a placeholder for the method that
can be redefined within derived component classes. Due to the inherent multitude of different
features possibly needed by component models, it is not practicable to write a method common
for all components; so the "proper" implementation of the method is left to the developer.
In the current version of the NRAP-Open-IAM, the ``SystemModel`` instance method
``single_step_model`` sets the values of all accumulators of the system
components to zero before each simulation: the accumulators are assumed to accumulate
sum-like quantities (e.g., cumulative masses and/or volumes).
For some accumulators the initial value of zero is assumed mainly out of
convenience. The redefinition of the ``reset`` method is needed
if, for example, the accumulator is supposed to keep track of product-like
quantities and the initial value should be 1, or for some other reason an initial
value of zero does not work. Additionally, the ``reset`` method can be used
for setting other component variables that need to be reinitialized before each simulation.
Setting the accumulator or other attribute values within the ``reset`` method
is straightforward. Referring to the same accumulators we used for the example
of the method ``_init_`` above, we replace the ``pass`` statement with appropriate
statements for this case. Note that we assume that the initial value
of ``accumulator1`` is 0, and the initial value of ``accumulator2`` is 1, thus,
only the ``accumulator2`` has to be reinitialized within the method. ::

    def reset(self):
        """
        Reset parameters, observations and accumulators.

        Parameters, observations and accumulators are reset to their
        initial/default values at the beginning of each new simulation.
        """
        self.accumulators['accumulator2'].sim = 1.0
        # self.accumulators['accumulator1'].sim = 0.0


Connecting components
=====================

There are many defaults methods of the ``ComponentModel`` class that are not meant
to be reimplemented within the derived component class. The main purpose of these
methods is to provide means to connect the components within the system model
and setup the parameters and observations of each component model.
By providing the examples illustrating the functionality of the available methods,
we want to show the capabilities of NRAP-Open-IAM which might be limited or
limit the development of new component models. Knowing what connections
can be created between models helps to make sure that the new model
is consistent with the available framework and can be
merged seemlessly. On the other hand, many of these methods were added or modified
based on the feedback of active developers, that is, new development is possible
after a justified proposal and review.

Attribute ``component_models`` of a ``SystemModel`` class object
is an ordered dictionary containing references to all component models
involved in the simulation. The order in which the components are arranged
in the system model is important and represents the order in which
the corresponding ``simulation_model`` methods are called during the simulation. Note that
the component ``simulation_model`` method should be developed in a way that would allow
it to be called at each time step provided by system model. After all component models
are run for a given time point, the control is returned back to the system model.
The instance method ``reorder_component_models`` of the ``SystemModel`` class changes
the order of the components in the system model after they have been added
to the container, and can be used to fix the order components models that have
been added out of execution order. ::

    def reorder_component_models(self, order):
        """
        Reorder execution order of component models.

        :param order: list of component model names in desired execution order
        :type order: lst(str)
        """
        self.component_models = OrderedDict((k, self.component_models[k]) for k in order)

There are many methods of the component model class
used to define the connections between the components that

* determine parameters of the model that the user can control (modify values)
  and specify observations that the user can analyze;

* determine which of the observations of a given component can provide
  input parameters of the next component;

* define which of the parameters of a given component model can be defined
  (calculated) in terms of the parameters of another component model;

* are used either at the script writing stage or are more likely to be used
  within the class derived from ``ComponentModel``.


Parameters
----------

We start with the description of the methods that allow the addition of parameters
and observations to the component models.

----

.. automethod:: openiam.iam_model_classes.ComponentModel.add_par

----

The ``add_par`` method can be used to add deterministic or stochastic parameter
by utilizing different values of the input argument *vary*. To add deterministic
parameter one would use ::

    cm.add_par(name=par_name, value=par_value, vary=False)

Here, *cm* is a name of variable containing a reference to the ``ComponentModel``
class instance.

By default the added parameter is assumed to be stochastic, thus, to add
stochastic parameter one would use ::

    cm.add_par(name=name, min=min_value, max=max_value, value=value)

Deterministic parameters of the component models can be accessed
by referring to ``deterministic_pars`` attribute of the component model class
instance, while the stochastic parameters of the component are kept in the ``pars``
attribute. The stochastic parameters are also added to and tracked by the
system model to which the component belongs. Similar to the component model,
the parameters can also be accessed through the ``pars`` attribute
of the system model. The only difference is that the system model ``pars``
attribute contain all stochastic parameters of the system,
i.e. stochastic parameters of all components.

----

.. automethod:: openiam.iam_model_classes.ComponentModel.add_default_par

----

The next method ``add_default_par`` is often used in the ``__init__`` method
of the component model class. The primary purpose of this method is to
ensure that all parameters of the model are defined during the simulation
even when not explicitly defined by the user. This is a convention
that all existing component models of NRAP-Open-IAM rely on. If the parameter
of a given component model is not defined by the user then the model itself
should take care of the default value of the parameter, for example,
assigning the value of the parameter not determined by the user within
the ``simulation_model`` method of the derived class.

As mentioned above the common approach to the default parameters is
to add them in ``__init__`` method as follows (the values of the parameters
are chosen arbitrarily for illustration) ::

    self.add_default_par('par_a', value=2.71)
    self.add_default_par('par_b', value=3.14)

Usually the two lines above would be followed by the definition
of the component attribute ``pars_bounds`` which would describe the lower and upper boundaries
of each model parameter (see notes above and example below). ::

    # Define dictionary of boundaries
    self.pars_bounds = dict()
    self.pars_bounds['par_a'] = [1.0, 15.0]
    self.pars_bounds['par_b'] = [-100.0, 5.0]

Note that the default value of the parameters can be either on the boundary
or inside the interval specified by ``pars_bounds``.

Definition of parameters default values allows utilizing them by other
components. Other component models can use the default parameter value
by referring to the ``default_pars`` attribute of the component model class instance.
For example, ::

    cm.add_default_par(par_name)
    print(cm.default_pars[par_name])

Observations
------------

The following methods allows adding observations that will be tracked
by the system model and will be available for different kinds of analysis.
Component models responses not explicitly added using these methods are ignored
and are not available for the processing once the simulation is complete.

----

.. automethod:: openiam.iam_model_classes.ComponentModel.add_obs

----

Note that the method ``add_obs`` allows adding only scalar observations.
The references to the observations of the component models are kept
in the identically named attributes ``obs`` of the component model and system model.
If the component model is supposed to return a structured (or "gridded") observation
(array, matrix) then the observation should be added with the method ``add_grid_obs``.

----

.. automethod:: openiam.iam_model_classes.ComponentModel.add_grid_obs

----

Due to the possibly large size, the structured observations are not tracked
by the system model (i.e. not available by reference to ``obs[obs_name]``
attribute of the system model) but the simulated structured observations are saved
in the compressed *.npz* format files with a filename defined by a pattern: component name,
gridded observation name, simulation number, and time index all joined
by undescore symbol **_**. For example, if the component with name *well* has
a structured observation with name *rate* then the file which contains the mentioned
observation evaluated at time point with index 10 will be named:

* *well_rate_sim_0_time_10.npz* for a single forward run of a system model and

* *well_rate_sim_1_time_10.npz*, *well_rate_sim_2_time_10.npz*, ...
  for multiple stochastic simulations of a system model.

The reference to the properties of the structured observations are kept in
the ``grid_obs`` attribute of the component model class instance. One can see the names
of the gridded observations added to the output of the given component
by accessing the mentioned attribute: ::

    cm.add_grid_obs('grid_obs1', constr_type='array')
    print(cm.grid_obs.keys())
    print(cm.grid_obs['grid_obs1'])

If only a single element of the structured (gridded) observation (array) is needed, then
method *add_local_obs* should be used. Since in this case the observation is a
scalar (not an array) the system model keeps track of the observation.

----

.. automethod:: openiam.iam_model_classes.ComponentModel.add_local_obs

----

Since the local observation is derived from the
structured observation, the component model additionally keeps the index
(or tuple of indices) of the element in an array (or matrix). In order
to get the index of the local observation in the structured observation
one has to know the name of the structured observation
from which the local observation is derived (name of the structured observation
is provided by the given component model) and the name of the local
observation assigned by the user. For example, ::

    # Add local observation
    cm.add_local_obs(loc_obs_name, grid_obs_name)

    # Print the index of the local observation
    print(cm.local_obs[grid_obs_name][loc_obs_name])


----

.. automethod:: openiam.iam_model_classes.ComponentModel.add_obs_to_be_linked

----

If the observation of a given component model is to be used as input
for other model it should be specifically added as such with the method
``add_obs_to_be_linked``. For example, ::

    cm.add_obs_to_be_linked(obs_name)

If the observation in addition should be tracked by the system model
then it also should be added with the method ``add_obs``, i.e. ::

    cm.add_obs(obs_name)               # to be tracked by system model
    cm.add_obs_to_be_linked(obs_name)  # to be used as input for other component

One important thing to mention before we continue with the next methods is that
the data (parameters, keyword arguments, observations) passed between models
need to have consistent units. There are some assumptions about the units of
the data that comes into and out of the component models. So here is a list of
units used by the NRAP-Open-IAM.

- Pressure is assumed to be in units of Pascals (Pa).
- Time is assumed to be in days (d).
- Length type parameters (distance, width, etc.) are assumed to be in units of meters (m).
- Fluxes of fluids are assumed to be in units of kilograms per second (kg/s).
- Masses are assumed to be in units of kilograms (kg).
- Viscosities are assumed to be in units of Pascal seconds (Pa s).

Linked parameters and observations
----------------------------------

The following methods allow adding component model parameters that depend
on other parameters and/or observations of the same or different component model.
Method ``add_par_linked_to_par`` adds a parameter that has the same properties and
value as an already defined parameter.

----

.. automethod:: openiam.iam_model_classes.ComponentModel.add_par_linked_to_par

----

In general, any of the model parameters can be of any allowed type. For example,
in one scenario a particular parameter can be setup as stochastic,
in another it can be setup as deterministic. Consider an example of the
code utilizing the ``add_par_linked_to_par`` method. As in the previous
examples, the values of the parameters are defined arbitrarily. ::

    # Assume we have component model cm1 defined somewhere above in the code
    # as instance of Component1 class with name 'cmpnt1'
    cm1.add_par('par_1', value=4.5, min=2.0, max=5.0)  # added parameter is stochastic
    cm1.add_par('par_2', value=1.8, vary=False)  # added parameter is deterministic
    cm1.add_default_par('par_3', value=1.5)      # added parameter is default

    # Now assume there is another component model cm2 also defined as
    # an instance of Component1 class with name 'cmpnt2'. We also assume that cm2
    # have parameters defined as dependent on the parameters of component model cm1
    # The corresponding dictionary should be used for each type of parameters of
    # component model cm1
    cm2.add_par_linked_to_par('par_1', cm1.pars['par_1'])
    cm2.add_par_linked_to_par('par_2', cm1.deterministic_pars['par_2'])

    # The names of the parameters of component cm2 and linked parameters of component cm1
    # should not necessary be the same if this is what is intended
    cm2.add_par_linked_to_par('par_3', cm1.default_pars['par_4'])
    cm2.add_par('par_4', value=2.5)                # added parameter is default

Method ``add_par_linked_to_obs`` adds a parameter that obtains its value from the
output of some other component. The observation that the added parameter
depends on should be added with the ``add_obs_to_be_linked`` method.

----

.. automethod:: openiam.iam_model_classes.ComponentModel.add_par_linked_to_obs

----

Consider the following example on the use of the method. ::

    # Add observation of component cm1 to be used as parameter of component cm2
    cm1.add_obs_to_be_linked('obs_1')
    # Add parameter of component cm2: parameter name of component cm2 does not
    # necessarily coincide with the name of observation returned by component model cm1
    cm2.add_par_linked_to_obs('par_1', cm1.linkobs['obs_1'])

Since the work of NRAP-Open-IAM is based on the assumption that parameters
of the component model are constant in time, the use of method
``add_par_linked_to_obs`` is appropriate only in situations when
the observation that is linked to the parameter does not vary in time.

Method ``add_composite_par`` adds a parameter whose value is determined by an
expression which may contain references to parameters of the same and/or other components.

----

.. automethod:: openiam.iam_model_classes.ComponentModel.add_composite_par

----

The following piece of code contain several examples utilizing
``add_composite_par`` method. ::

    # Assume we have component model cm1 defined somewhere above in the code
    # as instance of Component1 class with name 'cmpnt1'
    cm1.add_par('par_1', value=4.5, min=2.0, max=5.0) # added parameter is stochastic
    cm1.add_par('par_2', value=1.8, vary=False)  # added parameter is deterministic
    cm1.add_default_par('par_3', value=1.5)      # added parameter is default

    # Now assume there is another component model cm2 also defined as
    # an instance of Component1 class with name 'cmpnt2'.
    cm2.add_par_linked_to_par('par_1', value=2.5, min=2.1, max=4.9, vary=True)
    cm2.add_par_linked_to_par('par_2', min=2.0, max=5.0, vary=True)

    # We define par_3 of component cm2 as a sum
    # of the first three parameters of component cm1
    cm3.add_composite_par('par_3', expr='cmpnt1.par_1+cmpnt1.par_2+cmpnt1.par_3')

    # We define par_4 of component cm2 as a difference
    # of the first parameters of component cm1 and component cm2
    cm3.add_composite_par('par_4', expr='cmpnt1.par_1-cmpnt2.par_1')

Note that the expression for each added parameter directly utilizes the name of the parameter
which consists of the parental component name and parameter name
as defined during the component setup separated by a dot **.**. It does not require
knowing the parameter type (stochastic, deterministic, default). There is
a different way to write an expression for composite parameters which utilizes
the variables containing references to the corresponding components and names
of the parameters used in the expression. For example, line ::

    cm3.add_composite_par('par_3', expr='cmpnt1.par_1+cmpnt1.par_2+cmpnt1.par_3')

can be replaced with ::

    cm3.add_composite_par('par_3', expr='+'.join([cm1.pars['par_1'].name,
                                                  cm1.deterministic_pars['par_2'].name,
                                                  cm1.default_pars['par_3'].name]))

Then line ::

    cm3.add_composite_par('par_4', expr='cmpnt1.par_1-cmpnt2.par_1')

can be replaced with ::

    cm3.add_composite_par('par_4', expr='-'.join([cm1.pars['par_1'].name,
                                                  cm2.pars['par_1'].name]))

This method does not require knowing the name of the parental component
but rather the name of the variable that keeps the reference
to the corresponding component and the type of parameters involved in the
expression for the composite parameter.

Keyword arguments
-----------------

We discussed previously that inputs to the component model can be
of two main types: constant in time, scalar numerical parameters and
(possibly) varying in time model arguments. If the component's ``simulation_model``
method requires one or more of the later types, the model arguments have
to be added to the component model using one of the methods discussed below.
If the argument of the ``simulation_model`` method does not change with time:
for example, cannot be defined as a parameter of the model but might change
from one setup of the component to another, the simplest way to define it
is to make it an argument of the constructor method ``__init__``.
For example, see the definition above of *time_step* as both an argument
of ``__init__`` and ``simulation_model`` methods. Adding of ``'time_step'`` key
(and/or ``'time_point'``, ``'time_index'`` keys) to the dictionary attribute
``model_kwargs`` of the component tells the system model that these arguments
of the component ``simulation_model`` method are the same arguments provided by the system model.
If the attribute ``model_kwargs`` of the component model does not contain keys
``'time_step'``, ``'time_point'`` and ``'time_index'`` the component model will not be
aware of the values provided by the system model and will have to utilize
the default values provided with the definition of the ``__init__`` method. ::

    # Inside code of the __init__ method
    self.model_kwargs = {'time_step': 365.25, 'time_point': 365.25}

...
::

    # Inside script setting component model
    cm = Component1(name='cm1', parent=sm)  # sm is a system model cm belongs to
    cm.model_kwargs['time_step'] = 365.25

Argument *time_point* is an argument of the component model that changes in
time but is defined by the system model. To define the arguments that change
in time in a predetermined way, one can use method ``add_dynamic_kwarg``.

----

.. automethod:: openiam.iam_model_classes.ComponentModel.add_dynamic_kwarg

----

The main purpose of this method is to allow writing the scripts and/or tests
for a single component in the system model whose arguments otherwise
would have to depend on the output of other components. In the code example below,
we use the definition of the ``Component1`` class defined above. ::

    # Create a component
    cm = Component1(name='cm1', parent=sm)

    # Create an array of ten random real numbers between 1 and 5 that would
    # serve as an input for the model method
    t_array = 1.0 + 4.0*np.random.rand(10)
    # Add dynamic argument of the model method
    cm.add_dynamic_kwarg('temp_input1', t_array)

To satisfy possible needs of keyword arguments types of model methods,
``add_kwarg_linked_to_obs`` and ``add_kwarg_linked_to_collection`` were added
to the inventory of ``ComponentModel`` instance methods. Method
``add_kwarg_linked_to_obs`` can be used to connect keyword argument
of one component to the observation of another component.

----

.. automethod:: openiam.iam_model_classes.ComponentModel.add_kwarg_linked_to_obs

----

Observation to which keyword arguments are to be linked should be added
to the corresponding components (observation provider) as such with
``add_obs_to_be_linked`` method. Keyword arguments of the ``simulation_model`` method
can be linked to both types of observations in NRAP-Open-IAM: scalar and gridded/structured.
Depending on the type of the linked observation a corresponding set of arguments should be used.
We consider several examples illustrating the possible use of this method. ::

    # Suppose that component cm1 simulation_model method returns observation with name obs1
    # that can serve as a keyword argument for the component cm2.
    # Add observation of component cm1 to be linked to the argument
    # of the second component. Note that use of the obs_type argument is not
    # necessary when observation is scalar. Here, it is used to emphasize this fact
    cm1.add_obs_to_be_linked('obs1', obs_type='scalar')

    # Add keyword argument of cm2 linked to the observation of component cm1
    cm2.add_kwarg_linked_to_obs('input1', cm1.linkobs('obs1'))

As illustrated in the example the names of the keyword argument and of the
observation provided as input should not necessary be the same. The next
example illustrates the situation when the keyword argument is linked
to the gridded observation. Keyword argument can be linked to any part of
the gridded observation. Recall that the gridded observation should be
returned either as an array or matrix: ::

    # Suppose that component cm1 simulation_model method returns gridded observation
    # with name grid_obs1 that can serve as a keyword argument for the components
    # cm2, cm3 and cm4 linked to the different parts of the gridded observation
    # Add observation of component cm1 to be linked to the argument
    # of the second component. Note that use of the obs_type argument
    # is not necessary when observation is scalar.
    # Here, it is used to emphasize the type of observation
    cm1.add_obs_to_be_linked('grid_obs1', obs_type='grid')

    # It is a responsibility of the user to make sure that the returned
    # observation types are compatible with the format of keyword arguments
    # accepted by the subsequent components.
    # Add keyword argument of cm2 linked to the observation of component cm1
    cm2.add_kwarg_linked_to_obs('input1', cm1.linkobs('grid_obs1'))

    # Add keyword argument of cm3 linked to the several elements
    # of observation of component cm1
    cm3.add_kwarg_linked_to_obs('input1', cm1.linkobs('grid_obs1'),
                                constr_type='array', loc_ind=list(range(10)))

    # Add keyword argument of cm4 linked to a single element
    # of observation of component cm1
    cm4.add_kwarg_linked_to_obs('input2', cm1.linkobs('grid_obs1'),
                                constr_type='array', loc_ind=[0])

There are differences in three uses of the method ``add_kwarg_linked_to_obs``
that we want to discuss next. The keyword argument *input1* of the component
*cm2* ``model`` method is linked to the gridded observation with name
*grid_obs1* provided by component *cm1*. We do not use any
extra arguments of the method which shows that all observation data provided by *cm1*
is copied to the keyword argument. The keyword argument *input1* of the component *cm3*
is linked only to the part of the gridded observation with name *grid_obs1*
provided by component *cm1*. We indicate this by specifying the list of indices
of observation array at which it should be provided to the keyword argument *input1*.
Note that both components *cm2* and *cm3* should be able to accept the array-like
keyword argument *input1*. With component *cm4* the situation is slightly different.
Its keyword argument *input2* is also linked to the gridded observation *grid_obs1*
but it requests a single value of the gridded observation indicated
by a single index (in the example it is 0) of the element in observation array.
Note that although there is only a sinle index it still has to be provided in the list format.
The value of the observation array will be passed as a scalar variable rather
than array-like type as it is done for components *cm2* and *cm3*.

The following method is used to add a keyword argument created from
a collection of scalar observations provided by the same or different components.
Essentially, the observations are combined into one structure (collection)
and passed in this form to the corresponding component.

----

.. automethod:: openiam.iam_model_classes.ComponentModel.add_kwarg_linked_to_collection

----

The example below illustrates the use of the method for linking to the collection of similar observations. ::

    # Suppose that 5 components references to which are stored in a list variable sup_cm
    # can return an observation with name obs1
    # Add observation of components to be used to create a collection of observations
    for ind in range(5):
        # Option obs_type='scalar' can be omitted
        sup_cm[ind].add_obs_to_be_linked('obs1', obs_type='scalar')

    # Create a list of references to just created linked observations
    obs_collection = [sup_cm[ind].linkobs('grid_obs1') for ind in range(5)]

    # Add keyword argument of cm2 linked to the collection
    cm2.add_kwarg_linked_to_collection('input1', obs_collection)

The method is used to link to the collection of scalar observations.
We emphasize this by specifying option obs_type, non-mandatory for scalar observations.
The following simple example illustrates the situation when the collection
is created from not necessarily similar observations. ::

    # Suppose that component cm1 simulation_model method returns observation with name obs1
    # Add observation obs1 of component cm1 to be used to create
    # a collection of observations. Note that we omit option obs_type='scalar'
    # since it is not needed
    cm1.add_obs_to_be_linked('obs1')

    # Suppose that component cm2 simulation_model method returns observation with name obs2
    # Add observation obs2 of component cm2 to be used to create
    # a collection of observations. Note that we omit option obs_type='scalar'
    # since it is not needed
    cm2.add_obs_to_be_linked('obs2')

    # Create a list of references to just created linked observations
    obs_collection = [cm1.linkobs('obs1'), cm2.linkobs('obs2')]
    # Add keyword argument input1 of component cm3
    cm3.add_kwarg_linked_to_collection('input1', obs_collection)

In this example we assume that the ``simulation_model`` method of component *cm3* would accept
an argument *input1* which in this case is a 2-element array.
