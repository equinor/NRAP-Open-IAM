.. toctree::


.. _testsuite:

**********
Test Suite
**********

.. toctree::

Execution
=========

The NRAP-Open-IAM package comes with an automated test suite located
in the *test* directory.
The test suite can be run from the command land as::

    python iam_test.py

It is also run automatically after installation when running::

    python openiam_setup_tests.py

in the *setup* directory.
The test suite is also run automatically when developer's push changesets
(upload code) to the GitLab development code repository
(refer to :ref:`qaqc_dev` for details).

Test Documentation
==================

The following provides a description of each test in the test suite.

.. autoclass:: iam_test.Tests
   :members:
   :exclude-members: setUp, shortDescription
