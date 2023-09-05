******************
Tests and examples
******************

Tests
=====
Tests and examples can originate from the same process in development but
they fulfill very different roles. Tests are the integral part of the
code development and QA process. In general, tests are not for the
users application. There are many types of tests designed for specific
purposes at the different stages of the code development. At this point,
the test suite of NRAP-Open-IAM is made up of several types of tests:
unit tests (tests of single components), feature tests (particular functionality testing),
integration tests (whole system model tests) and installation (setup) tests.
Installation test checks user's Python environment and tool functionalility
work by starting all unit, feature and integration tests.
Each component model distributed as part of NRAP-Open-IAM has at least one test in the test suite
meant to check the observation values against expected results. We rely on ``unittest``,
a Python standard library, to facilitate test writing and execution.
It is considered a good practice to run a test suite after each update
of the code to make sure that changes do not cause unexpected changes in the output.
An approach consistent with the code testing practices of NRAP-Open-IAM
development team is described in :cite:`Reitz2016`.
Here, we state the most important rules that should be followed with all the new development.
Each test should be written for testing a particular feature or new component.
The test should be fast since with every major update the new functionality
should be tested which increases the total number of tests to be run.
The NRAP-Open-IAM tests are located in the file *iam_test.py* in the test *folder*.
Below we provide an example of the test written for a component utilizing a
Fortran library calculating the roots of the quadratic equation.
This component model method was described in detail in Chapter Coding Logistics.
Recall that the model method calculates the absolute values of the sum and
difference of two roots and values of the quadratic function for entered
coefficients and *x*-values.
We show the test which would check the solutions of the two quadratic
equations - one with both real roots and another one with complex roots -
against the known solutions. We create a test for the component as a method
of a ``ExampleTests`` class derived from ``unittest.TestCase`` class. Note that
we do not write a separate component class but rather utilize the available base class
``ComponentModel`` and create an instance of the class with model method specified by
``quad_eq_model`` (defined above in Chapter Coding Logistics). The example
provided below can be found in the subfolder *scripts* of folder *examples*
in the file *iam_simple_models.py*. ::

    class ExampleTests(unittest.TestCase):

        def test_quad_model(self):
            """
            Test work of quadratic model function and corresponding fortran library.
            """
            # Create system model
            sm = SystemModel()

            # Add component model with model function utilizing the fortran library
            # calculating the roots of quadratic equation
            qmc = sm.add_component_model('quad', model=quad_eq_model)

            # Add parameters of qmc component
            qmc.add_par('a', value=2.0)
            qmc.add_par('b', value=2.0)
            qmc.add_par('c', value=-12.0)
            # Add observations of qmc component
            qmc.add_obs('root_sum')
            qmc.add_obs('root_diff')

            # Run forward simulation
            sm.forward()
            # True roots for the defined a, b, c coefficients are -3 and 2.
            # Thus, absolute values of roots sum and difference are 1 and 5.
            true_vals = [1.0, 5.0]
            # Get simulated values of the observations
            sim_vals = [sm.obs['quad.root_sum'].sim, sm.obs['quad.root_diff'].sim]

            # Compare true and simulated values
            for tv, sv in zip(true_vals, sim_vals):
                self.assertTrue(abs((tv-sv)) < 0.001,
                                'The result is {} but should be {}'.format(sv,tv))

            # Check model output for complex roots
            # Change values of the component parameters
            qmc.pars['a'].value = 1.0
            qmc.pars['b'].value = -4.0
            qmc.pars['c'].value = 13.0

            # Run forward simulation
            sm.forward()
            # True roots for the defined a, b, c coefficients are 2+3i and 2-3i.
            # Thus, absolute values of roots sum and difference are 4 and 6.
            true_vals = [4.0, 6.0]
            sim_vals = [sm.obs['quad.root_sum'].sim, sm.obs['quad.root_diff'].sim]

            # Check results
            for tv, sv in zip(true_vals, sim_vals):
                self.assertTrue(abs((tv-sv)) < 0.001,
                                'The result is {} but should be {}'.format(sv,tv))

            return

In order to run the test one can use the following script: ::

    # Setup test runner
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stderr)
    # Create a test suite to which tests can be added
    test_suite = unittest.TestSuite()
    # Add corresponding test(s)
    test_suite.addTest(ExampleTests('test_quad_model'))
    # Execute added tests
    runner.run(test_suite)


A test runner is responsible for the execution of tests and returns the outcome to the user.
Test method name must start with word ``test``: in our example the name
of the test is ``test_quad_model``. It signals the test runner which of the methods
should be run. Additionally, each of the tests should contain the assertion statement involving
one of the following:

- ``assertEqual()`` to compare the obtained results versus the expected ones; or
- ``assertTrue()`` or assertFalse() to verify a particular condition; or
- ``assertRaises()`` to check whether a specific exception gets raised.

There are other available options (e.g., see ``unittest`` documentation
https://docs.python.org/3.6/library/unittest.html) that can be used
for specific checks. For example, the code ::

    for tv, sv in zip(true_vals, sim_vals):
        self.assertTrue(abs((tv-sv)) < 0.001,
                        'The result is {} but should be {}'.format(sv,tv))

can be replaced with ::

    for tv, sv in zip(true_vals, sim_vals):
        self.assertAlmostEqual(tv, sv, 2,
                               'The result is {} but should be {}'.format(sv,tv))

In the last code snapshot, the "2" argument in the assertAlmostEqual indicates
the number of places to which the difference of two values is rounded before comparison with zero.

At the final stages of the development of a new component the component's test
should be written and made available as part of the NRAP-Open-IAM test suite.


Examples
========
Integration of any component to the NRAP-Open-IAM framework involves not only
the development of the tests but also the creation of examples illustrating
the utility and capabilities of the new component. Examples are written to show users how
the component code works and how to interact with it. Ideally, examples should show off
all of the major features of the new development. Writing an example starts
with a basic description of the scenario. It helps to start with understanding which
additional components (if any) should be utilized. The NRAP-Open-IAM framework
allows running the system model with a single component included
by utilizing ``add_dynamic_kwarg`` method of the ``ComponentModel``. If the work
of the component can be shown without utilizing other components, it makes
sense to think whether the additional scenarios, where the component can be linked
to other component either as a source of input or output, can be created. The information
provided in the section Connecting components in Chapter Coding Logistics
can be very useful here. The examples illustrating work of the components
currently available as part of the NRAP-Open-IAM can be found in the *examples/scripts*
folder of the tool distribution. Examples often utilize the methods for
adding and creating component parameters and observations and connecting
the components within a given system model.
In general, most examples consist of the following steps:

- setup system model: simulation time, time step size or number of time points;
- add component models needed in the scenario;
- add system and components parameters and observations;
- create connections between components through models input/output
  using appropriate linking methods;
- run a chosen type of analysis: forward simulation, multiple stochastic
  simulations (Latin Hypercube sampling or Monte Carlo), parameter studies;
- collect observations from the system model;
- post-process the analysis results (if/as needed);
- plot or print the results.

The example files developed for the all components available within the NRAP-Open-IAM tool
can be found in the *examples/scripts* folder of the tool distribution. It is
strongly encouraged that all example files developed for the tool follow
the same naming convention: the file names should start with *iam_*
followed either by the list of components used for a specific scenario
and/or functionality of the tool featured in the example.
The following code illustrates an example of the system model containing
two linked component models available in NRAP-Open-IAM: reservoir model and
cemented wellbore model. Some of the components parameters are deterministic, some are
stochastic, some are composite or linked to other parameters. Keyword arguments
of the wellbore component (pressure and |CO2| saturation at the bottom of the well)
are obtained from the observations of the reservoir component. Leakage rates
of two fluids (|CO2| and brine) are calculated and shown on the produced figures.
The name of the file containing the example script is *iam_sys_reservoir_cmwell.py*
in the folder *examples/scripts*. ::

    import sys,os
    import matplotlib.pyplot as plt
    import numpy as np

    sys.path.insert(0,os.sep.join(['..','..','source']))
    from openiam import SystemModel, SimpleReservoir, CementedWellbore

    if __name__=='__main__':
        # Create system model
        num_years = 50.
        time_array = 365.25*np.arange(0.0, num_years+1)
        sm_model_kwargs = {'time_array': time_array}   # time is given in days
        sm = SystemModel(model_kwargs=sm_model_kwargs)

        # Add reservoir component
        sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))

        # Add parameters of reservoir component model
        sres.add_par('numberOfShaleLayers', value=3, vary=False)
        sres.add_par('shale1Thickness', min=500.0, max=550., value=525.0)
        sres.add_par('shale2Thickness', min=450.0, max=500., value=475.0)
        sres.add_par('shale3Thickness', value=11.2, vary=False)
        sres.add_par('aquifer1Thickness', value=22.4, vary=False)
        sres.add_par('aquifer2Thickness', value=19.2, vary=False)
        sres.add_par('reservoirThickness', value=51.2, vary=False)

        # Add observations of reservoir component model
        sres.add_obs_to_be_linked('pressure')
        sres.add_obs_to_be_linked('CO2saturation')

        # Add cemented wellbore component
        cw = sm.add_component_model_object(CementedWellbore(name='cw', parent=sm))

        # Add parameters of cemented wellbore component
        cw.add_par('logWellPerm', min=-13.9, max=-12.0, value=-12.0)

        # Add keyword arguments of the cemented wellbore component model
        cw.add_kwarg_linked_to_obs('pressure', sres.linkobs['pressure'])
        cw.add_kwarg_linked_to_obs('CO2saturation', sres.linkobs['CO2saturation'])

        # Add composite parameters of cemented wellbore component
        # Here, we illustrate two ways to define expressions for composite parameters
        # One way
        cw.add_composite_par('wellDepth', expr=sres.pars['shale1Thickness'].name+
            '+'+sres.pars['shale2Thickness'].name+
            '+'+sres.deterministic_pars['shale3Thickness'].name+
            '+'+sres.deterministic_pars['aquifer1Thickness'].name+
            '+'+sres.deterministic_pars['aquifer2Thickness'].name)
        # Second shorter (and more explicit) way
        cw.add_composite_par('depthRatio',
            expr='(sres.shale2Thickness+sres.shale3Thickness' +
            '+ sres.aquifer2Thickness + sres.aquifer1Thickness/2)/cw.wellDepth')
        cw.add_composite_par('initPressure',
            expr='sres.datumPressure + cw.wellDepth*cw.g*sres.brineDensity')

        # Add observations of the cemented wellbore component
        cw.add_obs('CO2_aquifer1')
        cw.add_obs('CO2_aquifer2')
        cw.add_obs('brine_aquifer1')
        cw.add_obs('brine_aquifer2')

        # Run forward simulation
        sm.forward()

        # Collect observations
        CO2_leakrates_aq1 = sm.collect_observations_as_time_series(cw, 'CO2_aquifer1')
        CO2_leakrates_aq2 = sm.collect_observations_as_time_series(cw, 'CO2_aquifer2')
        brine_leakrates_aq1 = sm.collect_observations_as_time_series(cw, 'brine_aquifer1')
        brine_leakrates_aq2 = sm.collect_observations_as_time_series(cw, 'brine_aquifer2')

        # Print results: CO2/brine leakage rates and pressure/saturation at the well
        print('CO2 leakage rates to aquifer 1:', CO2_leakrates_aq1, sep='\n')
        print('CO2 leakage rates to aquifer 2:', CO2_leakrates_aq2, sep='\n')
        print('Brine leakage rates to aquifer 1:', brine_leakrates_aq1, sep='\n')
        print('Brine leakage rates to aquifer 2:', brine_leakrates_aq2, sep='\n')

        # Plot CO2 and brine leakage rates along the wellbore
        plt.figure(1)
        plt.plot(sm.time_array/365.25, CO2_leakrates_aq1, color='blue',
                 linewidth=2, label='aquifer 1')
        plt.plot(sm.time_array/365.25, CO2_leakrates_aq2, color='red',
                 linewidth=2, label='aquifer 2')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.legend()
        plt.xlabel('Time, t [years]')
        plt.ylabel('Leakage rates, q [kg/s]')
        plt.title(r'Leakage of CO$_2$ to aquifer 1 and aquifer 2')

        plt.figure(2)
        plt.plot(sm.time_array/365.25, brine_leakrates_aq1, color='blue',
                 linewidth=2, label='aquifer 1')
        plt.plot(sm.time_array/365.25, brine_leakrates_aq2, color='red',
                 linewidth=2, label='aquifer 2')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.legend()
        plt.xlabel('Time, t [years]')
        plt.ylabel('Leakage rates, q [kg/s]')
        plt.title('Leakage of brine to aquifer 1 and aquifer 2')
        plt.show()
