.. _quickstart:

SPARTAN Quick Start Guide
=========================

This walkthrough reproduces plot panels similar to **Supplementary Fig. 1a–d** from *Single-molecule imaging of non-equilibrium molecular ensembles on the millisecond timescale*.

Overview
--------

Follow the steps below using the bundled example data and **Cornell SPARTAN** (compiled build). Menu paths refer to the **spartan** main window unless noted.

Step 1 — Install the program
------------------------------

#. Extract **SPARTAN example data.zip** to a fast local drive.
#. Extract **Supplementary Software 1.zip** and run ``Cornell_SPARTAN.exe``. Follow the on-screen installer.

   .. note::

      The compiled quick-start build targets **64-bit Windows 7** (or compatible).

Step 2 — Start the program
----------------------------

#. Launch **Cornell SPARTAN** from the Start menu.
#. Click **Browse…** and point to the **example data** directory.

Step 3 — Extract traces from the example movie
------------------------------------------------

#. In the main menu, choose **Extract fluorescence traces from movies** (``gettraces``).
#. In **gettraces**, use **Open Movie** and select ``sCMOS beads for alignment.tif`` from the example folder.
#. Under **Software Alignment**, choose **Iterative Closest Points**. After peaks are detected, use **Save…** to write ``align.mat``; then **Load** the same file so later movies reuse it.
#. **Open Movie** again and select ``sCMOS pre-steady-state cognate tRNA selection 10 ms.tif``.
#. Click **Save Traces**. This creates ``sCMOS pre-steady-state cognate tRNA selection 10 ms.rawtraces`` beside the movie.
#. Close **gettraces**.

Step 4 — Pre-steady-state analysis pipeline
-------------------------------------------

#. From the main menu, open **Pre-steady-state analysis pipeline** (``rtdgui``).
#. Keep defaults and click **Run**. When prompted, select ``sCMOS pre-steady-state cognate tRNA selection 10 ms.rawtraces``, then **Cancel** to finish file selection.

``rtdgui`` filters traces, idealizes, post-synchronizes, and splits **productive** vs. **non-productive** events. When processing finishes, a multi-panel window opens.

.. note::

   **Visualization layout** (columns often show **All**, **Productive**, **Non-prod.**):

   * Row 1 — post-synchronized ensemble (FRET–time) plots
   * Row 2 — per-state FRET histograms
   * Row 3 — transition density plots
   * Row 4 — state occupancy vs. time

Step 5 — Further exploration
------------------------------

**Per-trace review**

* From the main menu, open **Manually view, select, and correct traces** (``sorttraces``).
* **Open Traces File** → e.g. ``sCMOS pre-steady-state cognate tRNA selection 10 ms_auto_postSync_sel.traces``.
* Use **Previous** / **Next**; tick **Best FRET** for exemplars; **Save Selected Traces** writes a reduced file.

**Batch statistics**

* Open **Select traces according to defined criteria** (``autotrace``) → **Load Traces Files** → choose a ``.rawtraces`` file for histograms and custom selection rules.

.. seealso::

   Full feature documentation is in the distributed manual (**Cornell_SPARTAN_Documentation.pdf**) and the rest of this Sphinx guide.
