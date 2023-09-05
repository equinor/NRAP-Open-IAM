.. toctree::

.. _qaqc_dev:

************************
QA development practices
************************

NRAP-Open-IAM utilizes the following modern code development practices
to ensure QAQC during development:

* Code versioning using `Git <https://git-scm.com/>`_
* Online development code repository `<https://gitlab.com/NRAP/nrap-open-iam-dev>`_
  (private repository for developers)
* Issue tracking using GitLab (`<https://docs.gitlab.com/ee/user/project/issues/>`_)
* :ref:`testsuite` included with software package
* Continuous integration using GitLab: Installation and test suite are tested
  on GitLab cloud resources each time developers push code to the development
  repository. If either one fails, the developer is notified.
* Git branching to allow multiple branches of development prior to merging
  into the *master* branch.
