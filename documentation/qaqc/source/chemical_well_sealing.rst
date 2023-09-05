.. chemical_well_sealing:
.. include:: ../../common/replace_math.rst

Chemical Well Sealing
=====================

Class documentation
--------------------
.. autoclass:: openiam.ChemicalWellSealing

Unittests
----------
.. autoclass:: iam_test.Tests
   :members: test_chemical_sealing_coupled_simple_reservoir, test_chemical_sealing_not_seal_forward, test_chemical_sealing_seal_forward
   :noindex:

Additional QA documentation
---------------------------
The Chemical Well Sealing model is described in :cite:`WalshEtAl2013`, :cite:`WalshEtAl2014a`,
and :cite:`IyerEtAl2017`, and was calibrated using experimental data presented
in :cite:`WalshEtAl2013`, :cite:`WalshEtAl2014a`, and :cite:`WalshEtAl2014b`.

.. bibliography:: ../../bibliography/project.bib
    :filter: docname in docnames
