.. _rtd_pipeline:

Analysis pipeline for pre-steady-state experiments (rtdgui, rtdTool, and rtdPlots)
=====================================================================================

Real-time delivery (**RTD**) experiments need post-synchronization of asynchronous
events. **rtdgui** automates a pipeline starting from ``.rawtraces`` (typical
**gettraces** output):

#. **autotrace** filtering — excludes traces already high-FRET at :math:`t=0`.
#. **SKM** “on/off” idealization to locate first productive FRET.
#. Full idealization with a user kinetic model (default: aa-tRNA selection).
#. Optional **productive / non-productive** split by minimum dwell in a target state.
#. **makeplots**-style ensemble figures (extended).

.. note::

   **Figure 9.** **rtdgui** dialog — add ``.. figure::`` when available.

Running **rtdgui**
------------------

#. Pick a **kinetic model** (default E. coli aa-tRNA / tRNA–tRNA FRET); **Browse**
   for a custom ``.model``.
#. **Event separation** — enable to classify productive vs. non-productive using
   a minimum dwell (ms) in the chosen state.
#. **Run** → select one or more ``.rawtraces`` / ``.traces``; repeat until
   **Cancel** starts processing. Progress prints in the MATLAB **Command Window**
   when applicable.

The plot window columns correspond to outputs: with event separation, expect
**all / productive / non-productive** per file. Rows mirror **makeplots** plus
**state occupancy vs. time**.

.. note::

   **Figure 10.** Post-synchronized ensemble plots — add ``.. figure::`` when available.

Replotting saved results
------------------------

**Remake plot for already processed files** selects existing ``.traces`` outputs.
From MATLAB you can also call:

.. code-block:: matlab

   rtdPlots

Advanced / scripting
--------------------

Enable **advanced settings** to load a MATLAB script (``*.m``) that builds a
``customOpt`` structure—see ``advanced_settings_example.m`` in the **rtd2**
source folder (source build only).

The core implementation is **rtdTool**:

.. code-block:: matlab

   help rtdTool
