.. cemented_wellbore:
.. include:: ../../common/replace_math.rst

Cemented Wellbore Leakage Model
===============================

Class documentation
--------------------
.. autoclass:: openiam.CementedWellbore

Unittests
----------
.. autoclass:: iam_test.Tests
   :members: test_coupled_reservoir_cemented_wellbore_forward, test_lhs, test_openiam_cf, test_parstudy
   :noindex:

Additional QA documentation
---------------------------
The development and verification of the Cemented Wellbore ROM is documented
in :cite:`HARP2016150`. In April 2020, the Cemented Wellbore was updated
with significant improvements to its accuracy using a multiple ROM approach.
Click :download:`here <docs/wellbore_ROM_report_draft_20200916.pdf>` to view
a report detailing this update and its QA.

.. bibliography:: ../../bibliography/project.bib
    :filter: docname in docnames
