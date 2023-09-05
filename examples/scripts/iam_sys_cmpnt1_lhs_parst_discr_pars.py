'''
Example illustrates system model containing only simple_model from iam_simple_models.py
'''

import sys
import os
import matplotlib.pyplot as plt

sys.path.insert(0, os.sep.join(['..', '..', 'source']))
from openiam import SystemModel


if __name__ == '__main__':
    # For multiprocessing in Spyder
    __spec__ = None

    from iam_simple_models import simple_model
    # Create system model
    sm = SystemModel()

    # Add a component model
    tm = sm.add_component_model(name='tm', model=simple_model)

    # Add parameters of the component model
    num_values = 11
    par_values = list(range(20, 20+num_values))
    weights = num_values*[1./num_values]
    tm.add_par('var1', value=25,
               discrete_vals=(par_values, weights))
    tm.add_par('var2', value=4, min=2, max=12)
    tm.add_par('var3', value=10, vary=False)

    # Add observations
    tm.add_obs('output1')
    tm.add_obs('output2')
    tm.add_obs('output3')
    n2 = 2

    s1 = sm.lhs(siz=n2*num_values, seed=345)   # create LHS sample set

    s2 = sm.parstudy(nvals=[num_values, n2])   # create parstudy sample set

    # print('Sample set 1 values\n', s1.samples.values)

    # print('Sample set 2 values\n', s2.samples.values)

    plt.figure()
    plt.plot(s1.samples.values[:, 0], s1.samples.values[:, 1], 'ks', label='LHS')
    plt.plot(s2.samples.values[:, 0], s2.samples.values[:, 1], 'ro', label='parstudy')
    plt.title('Comparison of LHS and parstudy samples')
    plt.xlabel('var1')
    plt.ylabel('var2')
    plt.legend()

    # Run model using values in samples for parameter values
    results1 = s1.run(cpus=1, verbose=False)
    results2 = s2.run(cpus=1, verbose=False)

    plt.figure()
    colors = ['blue', 'red', 'green']
    for ind in range(3):
        plt.semilogy(range(1, n2*num_values+1), results1[:, ind], 'o',
                     color=colors[ind], label='output {}'.format(ind+1))
    plt.title('LHS')
    plt.xlabel('Realization number')
    plt.ylabel('Observations')
    plt.ylim(100, 1.0e+5)
    plt.legend()

    plt.figure()
    for ind in range(3):
        plt.semilogy(range(1, n2*num_values+1), results2[:, ind], 'o',
                     color=colors[ind], label='output {}'.format(ind+1))
    plt.title('Parameter study')
    plt.xlabel('Realization number')
    plt.ylabel('Observations')
    plt.ylim(100, 1.0e+5)
    plt.legend()
