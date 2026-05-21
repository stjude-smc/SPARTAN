.. _autotrace:

Select traces according to defined criteria (autotrace)
=======================================================

**autotrace** summarizes large smFRET datasets (:math:`\gtrsim 10^4` traces per
movie): signal-to-noise, bleaching, contamination, and other **traceStat**
metrics. Define thresholds to keep molecules with desired behavior (high SNR,
long lifetime, few artifacts). This is faster and more reproducible than purely
manual selection.

.. note::

   **Figure 3.** **autotrace** screenshot — add ``.. figure::`` when available.

Loading and histograms
----------------------

**File → Open File** (toolbar) or **File → Open Directory** (all ``.rawtraces``
in a folder). Histograms show each statistic; change the plotted statistic from
the dropdown. **File → Save Statistics** exports raw numbers. Algorithm details:
``traceStat.m`` (add new statistics there). Visual guide:

.. note::

   **Figure 4.** Commonly used selection statistics (donor/acceptor lifetimes,
   MTI/SNR, Pearson :math:`r`, bleaching-step detection, pre-steady-state
   criteria). Insert figure from documentation bundle.

Selection criteria (summary)
-----------------------------

See ``autotrace/traceStat.m`` for implementation details.

* **Mean total intensity** — average donor+acceptor when donor is on.
* **Highest FRET / FRET at first frame / Average FRET** — efficiency summaries.
* **Number of FRET events / FRET lifetime** — threshold crossings (lifetime uses
  runs of **≥ 5** frames).
* **Donor lifetime / End of trace** — bleaching-related frames.
* **Correlation of Fluor.** — Pearson :math:`r` (donor vs acceptor); **Derivative** variant.
* **SNR-bg / SNR-signal / SNR-signal (FRET)** — noise definitions as in GUI.
* **SNR-bg / SNR-signal** — ratio (rarely used).
* **Background noise** — σ after donor bleaching.
* **# Cy3 Blinks** — donor dark events excluding final bleach.
* **Multi-step photobleaching** — 0 = single step, 1 = multi-step.
* **Average FRET2 / Fret2 lifetime / Highest Fret2** — second acceptor channel, if present.

Applying and saving selections
-------------------------------

Use **Standard** and **Specialized** criteria on the right. Each line is
*statistic* **operator** *value* (e.g. donor lifetime **>** 50). Press **Enter**
to apply. Histograms update for the **selected** subset only.

**File → Save Traces** writes ``*_auto.traces`` plus a ``.log`` (criteria,
counts). **View → Traces** opens **sorttraces**; **View → Contour Plot** shows
ensemble contours.

Batch mode
----------

**Batch → Batch Analysis** walks directories/subdirectories. **Batch →
Auto-detect new files** watches a folder for real-time processing.

Artifacts and corrections
-------------------------

Accurate :math:`E` requires standard corrections (see benchmark literature
`DOI:10.1038/s41592-018-0085-0 <https://doi.org/10.1038/s41592-018-0085-0>`_).
Apply in the order below; changing order can give inconsistent results.

“Flagging” FRET histograms
++++++++++++++++++++++++++

Some samples with :math:`E > 0.8` show a jump in the first **~10** frames, then a
slow rise—often from picking defaults and background bleaching. Try increasing
**Frames to average for picking** (e.g. **50**) and enable **Subtract background
trace** in **gettraces** when appropriate.

Spectral crosstalk (``crosstalkcorrect``)
++++++++++++++++++++++++++++++++++++++++++

Bleed-through elevates apparent acceptor signal (especially :math:`E < 0.5`).
``gettraces`` applies the profile value; refine with ``crosstalkcorrect``.
``calc_crosstalk`` reports the ensemble-average crosstalk without writing a new
file.

Gamma correction (``gammacorrect``)
++++++++++++++++++++++++++++++++++++

Unequal detection/QY skews brightness vs. :math:`E` (e.g. Roy *et al.*,
`Nat. Methods 2008 <https://doi.org/10.1038/nmeth.1208>`_). **gammacorrect**
estimates γ from acceptor-photobleaching steps. **Note:** SPARTAN’s γ definition
is the **inverse** of some papers’ convention.

Flat-field / local gamma (``calcLocalGamma`` / ``applyLocalGamma``)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Field-dependent efficiency broadens histograms (large sensors). Calibrate with
**~0.5 FRET** samples via ``calcLocalGamma``; apply with ``applyLocalGamma``
(produces ``*_ffcorr.traces``). If histograms worsen, skip this step.

Acceptor direct excitation (``adecorrect``)
+++++++++++++++++++++++++++++++++++++++++++

Direct acceptor excitation biases low-:math:`E` data. Use ``adecorrect``; factor
from extinction ratios at the laser line (e.g. **532 nm**). Writes ``*_ade.traces``.

Other effects
+++++++++++++

Triplet/saturation can bias :math:`E` at high power (`Pati *et al.*, 2023`).
Prefer powers below onset of intensity-dependent :math:`E`. **Drift** is not
corrected in SPARTAN; see imaging literature (e.g. `Han *et al.*
<https://doi.org/10.1186/s13628-014-0015-1>`_).

.. note::

   **Figure** — placeholder for crosstalk / gamma illustration from the Word/PDF
   bundle (``FIGURE HERE`` in source doc).

Limitations
+++++++++++

Many steps threshold :math:`E` to define acceptor “alive”. Several tools ignore
ALEX direct-excitation timing and use :math:`E` thresholds instead. Prefer
**either** the batch correction functions **or** per-trace tweaks in
**sorttraces**, not both.
