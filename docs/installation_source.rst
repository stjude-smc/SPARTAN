.. _installation_source:

Installation from source code
=============================

The source build requires **MATLAB** but exposes more functionality and allows
custom code that calls SPARTAN functions.

Prerequisites
-------------

Install 64-bit MATLAB with the required toolboxes (see :doc:`introduction`).

Download and path
-----------------

#. Download the latest SPARTAN source:
   `github.com/stjude-smc/SPARTAN <https://github.com/stjude-smc/SPARTAN>`_
#. Unzip to a stable location (e.g. ``C:\SPARTAN``).
#. In MATLAB, open **HOME → Set Path**:

   * Click **Default** to remove non-default entries.
   * **Add with Subfolders** → select the SPARTAN source folder.
   * **Save**, then **Close**.

Optional: parallel pool
-----------------------

Open the parallel pool indicator (lower left) → **Parallel Preferences**. Set
the shutdown period to **3000** minutes so the pool stays warm. On systems with
**≤ 32 GB** RAM, limit workers to **at most 4**.

Optional: ebFRET
----------------

Download `ebFRET 1.1.1 <https://github.com/ebfret/ebfret-gui/releases/>`_ and
add its ``src`` folder to the MATLAB path so **batchKinetics** can use it as an
additional modeling method.

Launch
------

From the MATLAB prompt:

.. code-block:: matlab

 spartan

Known issues (source)
---------------------

.. warning::

   Some multi-processor machines stall after ``parfor`` loops. If so, set
   ``constants.enable_parfor`` to ``false`` in ``cascadeConstants.m``.

Other toolboxes or scripts on the path can conflict with SPARTAN. Reset the path
so only one SPARTAN tree is present, and keep ``Documents/MATLAB`` free of stray
``.m`` files when possible.

If problems persist, restart MATLAB; then contact
`Scott.Blanchard@stjude.org <mailto:Scott.Blanchard@stjude.org>`_ or file an
issue on GitHub.
