.. _gettraces:

Extracting fluorescence traces from movies (gettraces)
======================================================

Under TIRF wide-field imaging, movies contain many PSFs, each from a single
labeled particle. **gettraces** detects peaks, maps channels, sums intensities
into time traces, and computes corrected **FRET** traces (efficiency
:math:`E`). Open it from **spartan → Extract fluorescence traces from movies**.

.. note::

   **Figure 2.** Screenshot of **gettraces** — add ``.. figure::`` when available.

Select instrument profile
-------------------------

Choose an **Imaging Profile** from **Settings** at the top. Profiles encode the
setup (channels, defaults). As of SPARTAN **3.10**, new profiles are edited in
``cascadeConstants.m`` (see comments there).

Use **Customize** at the bottom of **Settings**, or the toolbar, for the active
profile. Right-click a field image for **Change field settings…** or **Change
field arrangement…**. Crosstalk and scaling: **Settings → Set crosstalk values…**
and **Settings → Set channel scaling…**.

Changes apply only to the current window; restarting MATLAB or a new
``gettraces`` session reloads ``cascadeConstants.m``.

Loading and viewing movies
--------------------------

**File → Open** or the toolbar: pick a ``.tif`` stack. By default, field images
average the **first 10 frames**; the sum image aids peak finding. Background is
interpolated and subtracted. Use the **left slider** for display scaling; use
the **slider + Play** under fields to scroll frames.

Detecting particles
-------------------

**View → Pick Peaks** (or toolbar): local maxima above threshold in the total
intensity image. Default threshold is a few **σ** above background; override via
**Settings → Customize**. Circles mark centers; **gray** circles mark pairs
rejected for being too close (overlap contamination). Tune **minimum separation**
in **Settings → Customize**. Aim for **~40%** rejected for dense loading. Counts
show **selected / (all)** molecules.

Alignment of spectral channels
-------------------------------

By default **ICP** maps donor/acceptor PSFs. After convergence, peaks are
redetected on the **aligned** sum image. Toggle ICP or load a saved alignment
from **Alignment**.

ICP needs usable signal on all channels (roughly :math:`E \in [0.2,0.8]`). If
ICP fails, record **bead** movies, pick peaks, then **Alignment → Memorize**.
The table shows transform parameters: **Deviation** (sub-pixel residual; aim
**< 0.4** px) and **Contrast** (Weber score; **> 1.4** is strong).

Extract traces, corrections, and FRET
--------------------------------------

**File → Save Traces** (or toolbar): sum pixels in a neighborhood per frame;
sizes are set under **Settings → Customize** (tradeoff signal vs. noise; **~80%**
of PSF flux is typical—see **% intensity collected**). Post-bleach background is
subtracted; **crosstalk** and **acceptor scaling** corrections are applied; FRET
is computed (zero when donor is dark). Donor blinking: threshold or **SKM** in
**Settings → Customize**. Output: ``.rawtraces`` beside the movie name.

Batch processing
----------------

**Batch → Select Directory** processes all movies with current settings;
optional subdirectory search and **skip if** ``.rawtraces`` exists.
**Batch → Auto-detect new movies** watches a folder (useful during acquisition).
Monitor **Deviation**, PSF size, and rejection fraction during automation.

Analysis settings (**Settings → Customize**)
++++++++++++++++++++++++++++++++++++++++++++

**Name** — instrument label.

**Picking Threshold Value** — minimum peak intensity.

**Use Automatic Threshold** — default uses **8 σ** of background (adjust via
**Auto picking sensitivity**).

**Integration window size (px)** — pixels summed per PSF (target **~80%**
collection).

**Integration neighborhood (px)** — region is :math:`2n+1` per side
(:math:`n=1` → :math:`3\times3`).

**Minimum separation (px)** — Euclidean cutoff for overlap rejection.

**Frames to average for picking** — default **10**; lower for fast bleaching.

**Background trace field** — saves ``bgTrace`` in metadata when enabled.

**Subtract background image** — interpolated background map (recommended).

**Subtract baseline trace** — removes frame-wise baseline drift (injections,
contaminant photobleaching).

Defaults assume **60×** objective and **13 µm** binned pixels—retune for other
geometries.

Field and arrangement dialogs
+++++++++++++++++++++++++++++++

**Field settings:** **Role** (donor/acceptor), **Name**, **Wavelength**,
**photonsPerCount** (vendor ADU→photon).

**Field arrangement:** **rows × columns** split; channel dropdown per tile
(blank = ignore).

Troubleshooting
---------------

Metrics turn **red** when out of range. Large PSFs often mean focus or dirty
optics. Edge misses suggest uneven illumination or vignetting. High **Deviation**
or asymmetric rejections point to alignment or optics. Dim channels (Cy2/Cy7)
may show low contrast even when alignment is valid. Zoom to verify alignment
when scores are marginal.
