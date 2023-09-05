import ctypes
import unittest
import os
import sys
from sys import platform
import numpy as np
import matplotlib.pyplot as plt

try:
    from openiam import SystemModel
except:
    try:
        sys.path.append(os.sep.join(['..', '..', 'source']))
        from openiam import SystemModel
    except ImportError as err:
        print('Unable to load IAM class module: '+str(err))
currentWorkDir = os.getcwd()


class ExampleTests(unittest.TestCase):

    def test_quad_model(self):
        """ Test work of quadratic model function and corresponding fortran library.

        """
        # Create system model
        sm = SystemModel()

        # Add component model with model function utilizing the fortran library
        # calculating the roots of quadratic equation
        qmc = sm.add_component_model('quad', model=quad_eq_model)
        qmc.grid_obs_keys = ['y'] # to let system know to handle grid obs differently

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
            self.assertAlmostEqual(
                tv, sv, 2, msg='The result is {} but should be {}'.format(sv, tv))

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
                            'The result is {} but should be {}'.format(sv, tv))

def simple_model(p, input_array=None):
    """
    Return output based on the provided input parameters.

    :param p: dictionary of parameters passed to the function
    :type p: dict

    :param input_array: input array
    :type input_array: numpy.array

    :returns dictionary containing 3 output variables and array
    """
    # If parameters are floats and can vary stochastically, they should be
    # passed in p, dictionary of input parameters. If model parameters are
    # arrays they should be passed as separate arguments, similar to input_array.

    # Check which parameters are provided in the input dictionary
    if 'var1' in p:
        var1 = p['var1']
    else:
        var1 = -1.0  # default values if the desired parameter value was not passed

    if 'var2' in p:
        var2 = p['var2']
    else:
        var2 = 2.0

    if 'var3' in p:
        var3 = p['var3']
    else:
        var3 = -3.0

    # Check whether input_array was passed and define output array
    if input_array is None:
        out_array = np.ones(10)
    else:
        out_array = 4*input_array - 7*var1 + np.log10(var2) - var3**2

    # Define output variables of the component model
    out_var1 = var1**2 + var2*var3
    out_var2 = var2**3 + var1*var3
    out_var3 = var3**2 + var1*var2 + var1**3

    # Define output of the function (output is also a dictionary)
    output = dict()
    output['output1'] = out_var1  # float
    output['output2'] = out_var2  # float
    output['output3'] = out_var3  # float
    output['output4'] = out_array  # array

    # Component model should return a dictionary of outputs
    return output

def quad_eq_model(p, input_array_x=None):
    '''
    Return three outputs based on the provided input parameters.

    :param p: dictionary of parameters passed to the function
    :type p: dict

    :param input_array_x: input array
    :type input_array_x: numpy.array

    :returns dictionary containing output variable and array
    '''
    if 'a' in p:
        a = p['a']
    else:
        a = 3.0

    if 'b' in p:
        b = p['b']
    else:
        b = 5.0

    if 'c' in p:
        c = p['c']
    else:
        c = 2.0

    if input_array_x is None:
        input_array_x = np.linspace(-10.0, 10.0, 1000)

    # Determine the size of input array
    N = np.size(input_array_x)

    # Setup library and needed function names
    if platform in ("linux", "linux2"):
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

    # Define c classes to be used for inputs and outputs of the fortran function
    INT = ctypes.c_int
    DOUBLE = ctypes.c_double
    NPointsArrayType = DOUBLE*N

    # Set argument types for values and pointers
    # The order should coincide with the order of arguments in the fortran function
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

    # Define values of the input parameters that will be passed to fortran function
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
    # Suppose that model returns absolute values of sum and difference of the roots and
    # y-values of the function. Set of calculations depends on the value of flag
    # Check flag values
    if out_flag.value==1:          # two real distinct roots
        # extract value from the float type of output
        out['root_sum'] = abs(out_x1.value + out_x2.value)
        out['root_diff'] = abs(out_x1.value - out_x2.value)
    elif out_flag.value==2:        # single real root
        out['root_sum'] = abs(2*out_x1.value)
        out['root_diff'] = 0.0
    else:                    # two complex roots
        out['root_sum'] = abs(2*out_x1.value)       # sum of roots is a doubled real part
        out['root_diff'] = abs(2*out_x2.value)      # difference of roots is a doubled complex part

    out['y'] = out_array_y[0:N]        # extract values from array type of output

    # Return output dictionary
    return out


def fortran_based_model(p, input_array=None):
    '''
    Return two outputs based on the provided input parameters.

    :param p: dictionary of parameters passed to the function
    :type p: dict

    :param input_array: input array
    :type input_array: numpy.array

    :returns dictionary containing output variable and array
    '''
    # Check what parameters are provided in the input dictionary
    # float input variable
    if 'var1' in p:
        var1 = p['var1']
    else:
        var1 = 1.0

    if 'var2' in p:
        var2 = p['var2']
    else:
        var2 = 3.0

    if input_array is None:
        input_array = np.array([3.4, 7.8, 6.5, 8.4, 7.6, 2.1, 6.4, 6.2, 6.5])

    # Determine the size of input array
    N = np.size(input_array)

    # Setup library and needed function names
    if platform in ("linux", "linux2"):
        # linux
         library = "example_rom.so"
    elif platform == "darwin":
        # OS X
        library = "example_rom.dylib"
    elif platform == "win32":
        # Windows...
        library = "example_rom.dll"
    functionName = "example_rom"

    # Define c classes
    INT = ctypes.c_int
    DOUBLE = ctypes.c_double
    NPointsArrayType = DOUBLE*N

    # Load DLL
    external_lib = ctypes.cdll.LoadLibrary(os.path.join(os.getcwd(),library))

    # Get needed function as attribute of the library
    function = getattr(external_lib, functionName)

    # Set argument types for value and pointers
    function.argtypes = [ctypes.POINTER(DOUBLE), ctypes.POINTER(DOUBLE),
                         ctypes.POINTER(INT),
                         ctypes.POINTER(DOUBLE),
                         ctypes.POINTER(DOUBLE)]
    function.restype = None

    # Define values of the input parameters that will be passed to fortran function
    fun_N = INT(N)
    fun_var1 = DOUBLE(var1+var2)
    fun_array1 = NPointsArrayType(*input_array)

    out_var1 = DOUBLE()                # initialize output variable
    out_array1 = NPointsArrayType()    # initialize output array

    function(fun_var1, fun_array1, fun_N, out_var1, out_array1)

    # Create output dictionary
    out = dict()
    out['out_var1'] = out_var1.value   # extract value from the float type of output
    out['out_array1'] = out_array1[0:N]   # extract values from array type of output

    # Return output dictionary
    return out

def test_ROM():
    """
    Test the correct work of fortran code compiled into library.
    """
    # Setup library and needed function names
    if platform in ("linux", "linux2"):
        # linux
         library = "example_rom.so"
    elif platform == "darwin":
        # OS X
        library = "example_rom.dylib"
    elif platform == "win32":
        # Windows...
        library = "example_rom.dll"
    functionName = "example_rom"

    # Define c classes
    INT = ctypes.c_int
    DOUBLE = ctypes.c_double

    N = 5
    FivePointsArrayType = DOUBLE*N

    # Load DLL
    external_lib = ctypes.cdll.LoadLibrary(os.path.join(os.getcwd(),library))

    # Get needed function as attribute of the library
    function = getattr(external_lib, functionName)

    # Set argument types for value and pointers
    function.argtypes = [ctypes.POINTER(DOUBLE), ctypes.POINTER(DOUBLE),
                         ctypes.POINTER(INT),
                         ctypes.POINTER(DOUBLE),
                         ctypes.POINTER(DOUBLE)]
    function.restype = None

    # Define values of the input parameters
    N = INT(5)
    var1 = DOUBLE(13.5)
    array1 = FivePointsArrayType(1.1, 4.5, 8.1, 2.7, 9.5)

    out_var1 = DOUBLE(3.6)
    out_array1 = FivePointsArrayType()
    function(var1, array1, N, out_var1, out_array1)

    # Return output array
    return out_array1[0:5]


if __name__=='__main__':
    # Test set up example
    res = test_ROM()
    print('Output array of example ROM:', res)

    # Test fortran_based_model function
    parameters = dict()
    parameters['var1'] = 3.5
    parameters['var2'] = 7.5
    input_array = np.array([3.4, 7.8, 6.5, 8.4])

    out = fortran_based_model(parameters, input_array)
    print('============ Fortran dll based simple model ====================')
    print(out['out_array1'])
    print(out['out_var1'])

    # Test Fortran based quadratic equation solver function
    parameters = dict()
    a = [1.5, 2.0, 2.0]
    b = [6.0, 4.0, 5.5]
    c = [2.0, 2.0, 7.0]
    colors = ['red', 'blue', 'orange']
    line_labels = ['2 real roots', '1 real root', '2 complex roots']
    input_array = np.linspace(-5.,2.,1000)
    plt.figure(1)

    print('======= Fortran dll based quadratic equation solver model =========')
    for k in range(3):
        parameters['a'] = a[k]
        parameters['b'] = b[k]
        parameters['c'] = c[k]
        print('============== {} case ============ '.format(k))
        out = quad_eq_model(parameters, input_array)
        print(out['root_sum'])
        print(out['root_diff'])
        print(out['y'])
        plt.plot(input_array, out['y'], '-', color=colors[k], label=line_labels[k])
    plt.plot([-5., 2.],[0.0, 0.0], '-k')
    plt.legend()

    # Test simple_model function
    parameters = dict()
    parameters['var1'] = 3.5
    parameters['var2'] = 0.5
    parameters['var3'] = 1.5
    input_array = np.array([3.4, 7.8, 6.5, 8.4, 0.4, 3.5, 8.7])
    out = simple_model(parameters, input_array)
    print('========================== Simple model ========================')
    print(out['output1'])
    print(out['output2'])
    print(out['output3'])
    print(out['output4'])

    print('=================== Test models with test suite ================')

    # Setup test runner
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stderr)
    # Create a test suite to which tests can be added
    test_suite = unittest.TestSuite()
    # Add corresponding test(s)
    test_suite.addTest(ExampleTests('test_quad_model'))
    # Execute added tests
    runner.run(test_suite)
