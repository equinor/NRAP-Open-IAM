.. open_wellbore:
.. include:: ../../common/replace_math.rst

Open Wellbore Model
===================

Class documentation
--------------------
.. autoclass:: openiam.OpenWellbore

Unittests
----------
.. autoclass:: iam_test.Tests
   :members: test_coupled_reservoir_open_carbonate_forward, test_coupled_reservoir_open_wellbore_forward
   :noindex:

Additional QA documentation
---------------------------
The verification of the Open Wellbore component is documented in the technical
report :cite:`BaconEtAl2021` available for download
:download:`here <docs/open_wellbore_report.pdf>`. The component was
updated with additional lookup lables to improve the accuracy and extend
the possible range of input parameters.

.. bibliography:: ../../bibliography/project.bib
    :filter: docname in docnames
