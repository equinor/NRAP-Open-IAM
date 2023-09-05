************
Introduction
************

.. toctree::

This document describes the software design and functional requirements
of the Phase II Integrated Assessment Model (IAM), referred to as the
NRAP-Open-IAM in the current document. This design process adapted
a use-case-driven approach :cite:`RosenbergScott1999` to design a general purpose
open source version of the IAM, built upon the development effort of
Phase I NRAP-IAM-CS :cite:`PawarEtAl2016`, :cite:`StaufferEtAl2015`. Based on this
general-purpose design, a core-functionality prototype version of the NRAP-Open-IAM
was developed and tested.

In NRAP Phase II researchers are developing an integrated assessment
model that will incorporate workflows for containment assurance,
monitoring design, post-injection site care, assessment of model concordance to measured field data,
evaluation of the performance of mitigation alternatives, and updating of
probabilistic assessments as new data become available. The design goals
are to create a flexible framework for the integrate assessment model
that will be easily maintainable, flexible, and extendable.

The document is intended to outline the NRAP-Open-IAM design structure and functionality
and give guidance to the developers of reduced order models (ROMs),
monitoring modules, utility applications and algorithm developers
for interfacing and integrating their models into the NRAP-Open-IAM. The scope of the functionalities
of the NRAP-Open-IAM include risk quantification, risk management and support
of iterative risk assessment processes. The system components include
the geophysical and geochemical model components for risk quantification,
analysis components for risk management and parameter passing workflows,
procedures and mechanisms allowing the communications between
the different system components to support iterative risk assessment processes.
