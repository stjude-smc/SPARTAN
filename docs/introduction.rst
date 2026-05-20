.. _introduction:

Getting Started
===============

Documentation for **Single-molecule Platform for Automated, Real-Time Analysis (SPARTAN)**.

Authorship and version
----------------------

**Daniel Terry**

St Jude Children’s Research Hospital

262 Danny Thomas Pl, Memphis, TN 38105

**November 19, 2025** (version **3.11.0**)

License
-------

.. note::

   This software is licensed under **CC BY-NC-SA 4.0**. For users affiliated with
   commercial or for-profit organizations, use requires a specific licensing
   agreement. Contact the Office of Technology Licensing at St. Jude
   (`technology.licensing@stjude.org <mailto:technology.licensing@stjude.org>`_)
   for terms and licensing.

Project links and contact
-------------------------

* **Repository and issues:** `SPARTAN on GitHub <https://github.com/stjude-smc/SPARTAN>`_
* **Collaboration / inquiries:** `Scott.Blanchard@stjude.org <mailto:Scott.Blanchard@stjude.org>`_

Overview
--------

This documentation describes software for the analysis of single-molecule
fluorescence and **FRET** (smFRET) data. Written in **MATLAB**, it includes tools
for extracting traces from movies, selecting traces, applying corrections,
hidden Markov modeling, simulation, and visualization. The apparent FRET
efficiency is often denoted :math:`E` (dimensionless, typically in
:math:`[0,1]`). Further background appears in **Juette and Terry** *et al.*, **(2016)**


:Title:
    Single-molecule imaging of non-equilibrium molecular ensembles on the millisecond timescale
:Journal:
    *Nat. Methods* **16**, 669–676
:DOI:
    `doi: 10.1038/nmeth.3769`_


Recommended system requirements
-------------------------------

* **MATLAB** (source build): **2019b** or later, with Image Processing,
  Statistics, Curve Fitting, Optimization, and Distributed Computing toolboxes.
* **Operating system:** Microsoft Windows 10 (macOS and Linux may work but are
  not the primary supported platforms).
* **CPU:** 64-bit Intel-compatible x86 with at least **four** physical cores.
* **Memory:** **32 GB** RAM.
* **Storage:** **~100 GB** free on a modern solid-state drive.

Version numbering
-----------------

Versions use **MAJOR.MINOR.REVISION**:

* **Major:** may break compatibility (file formats, core syntax, major features).
* **Minor:** new features and bug fixes; largely backward compatible.
* **Revision:** small bug-fix–only updates.
