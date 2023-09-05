.. seal_horizon:
.. include:: ../../common/replace_math.rst

Seal Horizon Model
==================

Class documentation
--------------------
.. autoclass:: openiam.SealHorizon

Unittests
----------
.. autoclass:: iam_test.Tests
   :members: test_seal_horizon
   :noindex:

Additional QA documentation
---------------------------
Additional information about the standalone Seal Flux model behind the Seal Horizon component
can be found in the model user guide :cite:`Lindner2022` available for download
:download:`here <docs/seal_flux_user_manual.pdf>`.
In addition to the test mentioned above several scripts and control file examples were
developed to compare the standalone code Seal Flux with its integrated version
of Seal Horizon component in NRAP-Open-IAM. Although not all results
can be recreated in NRAP-Open-IAM: specifically, the ones involving the stochastic
simulations, the results of deterministic runs were identical to the ones produced
by standalone code. The developed tests did not provide full comparison of the
capabilities of both tools as some of Seal Flux capabilities
are not fully integrated into NRAP-Open-IAM.

.. bibliography:: ../../bibliography/project.bib
    :filter: docname in docnames
