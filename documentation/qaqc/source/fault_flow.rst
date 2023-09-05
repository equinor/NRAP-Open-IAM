.. fault_flow:
.. include:: ../../common/replace_math.rst

Fault Flow Model
================

Class documentation
--------------------
.. autoclass:: openiam.FaultFlow

Unittests
----------
.. autoclass:: iam_test.Tests
   :members: test_fault_flow
   :noindex:

Additional QA documentation
---------------------------
In addition to the mentioned test several scripts and control file examples were
developed to compare the standalone code Fault Flo with its integrated version
of Fault Flow component in NRAP-Open-IAM. Although not all results
can be recreated in NRAP-Open-IAM: specifically, the ones involving the stochastic
simulations, the results of deterministic runs were identical to the ones produced
by standalone code.

.. bibliography:: ../../bibliography/project.bib
    :filter: docname in docnames
