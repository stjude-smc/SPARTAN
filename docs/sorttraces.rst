.. _sorttraces:

Viewing traces, manual corrections, and trace selection (sorttraces)
=========================================================================

.. note::

   **Figure 5.** **sorttraces** screenshot — add ``.. figure::`` when available.

Loading ``.traces`` files
-------------------------

**File → Open** (toolbar): load ``.rawtraces`` or ``.traces``. Layout:

* **Top** — donor/acceptor intensity traces.
* **Middle** — total intensity (black) and threshold (blue); FRET forced to **0** below threshold.
* **Bottom** — FRET (blue) and idealization (red) when present.

**Navigation** (lower right): **Previous / Next**, or jump by trace index. Filename
and index appear at the top.

Info panel
----------

* **Lifetime** — donor then acceptor lifetimes before bleaching.
* **Correlation** — Pearson :math:`r` (full trace and zoomed region).
* **SNR** — total intensity vs. background σ (left) and pre-bleach signal (right).
* **Molecule location** — red circle in the field of view.

Manual bins
-----------

Assign traces to **No FRET / All FRET / Best FRET** (names editable via
right-click → **Rename Bin**). Counts show above each checkbox. Use the bin menu
to clear or advance within a bin.

Saving
------

**File → Save Selected Traces** writes one file per non-empty bin
(``*_no_fret.traces``, etc.) and removes empty-bin files. Corrections and picks
are stored in ``*_savedState.mat`` for later resume.

Idealization
------------

**sorttraces** auto-loads ``.dwt`` / ``.qub.dwt`` beside the ``.traces`` file.
**Edit → Load Idealization** / **Clear Idealization**; idealizations export with
saved bins when a template was loaded.

Corrections
-----------

* **Background:** zoom to a flat region → **Corrections → Subtract All**;
  **Reset Background** undoes.
* **FRET threshold:** drag the blue line (disable zoom); or **Edit → Set Zero
  Method → SKM** (also configurable in **gettraces**).
* **Crosstalk:** **D/A crosstalk**, **D/A2**, **A1/A2** for multi-color.
* **Scaling:** **A1 Scaling**, **A2 Scaling** so total intensity is flat across :math:`E` changes.

**Save current file as…** writes the whole corrected file. **Corrections →
Reset All** reloads originals.

VIII. Visualizing trace ensembles (**makeplots**)
-------------------------------------------------

**makeplots** summarizes ensembles: **row 1** FRET–time contours + molecule
count; **row 2** either pooled histograms (``frethistComparison`` overlays) or
per-state histograms when a ``.dwt`` exists; **row 3** transition densities,
:math:`N_t`, and mean rate. **Edit → Settings** persists options; **Edit →
Reset settings** restores defaults. **File → Export .txt** saves histogram
tables for Origin and similar tools.

.. note::

   **Figure 6.** Ensemble plots from **makeplots** — add figure when available.
