.. _spartan_workflow:

Main menu and core analysis workflow (spartan)
==============================================

The **spartan** GUI is the main entry point for the suite. For compiled builds it
opens at launch. The current working directory is shown at the top; change it
with **Browse…**.

Typical workflow
----------------

#. Acquire wide-field movies as **TIFF** stacks.
#. Run ``gettraces`` to extract traces and compute FRET.
#. Run ``autotrace`` for statistics and subset selection.
#. Run ``makeplots`` for ensemble visualization.
#. Run ``sorttraces`` for per-trace review and selection.
#. Run ``batchKinetics`` for hidden Markov modeling (**HMM**).
#. Use **batchKinetics** simulation to test models against data.

.. note::

   **Figure 1.** SPARTAN main menu — add ``.. figure:: /_static/your_menu.png`` when
   a screenshot is available.

Acquiring TIRF smFRET data
--------------------------

Acquisition software reads camera frames and saves **TIFF** stacks. Examples:

* `µManager <https://www.micro-manager.org/>`_
* `MetaMorph <http://www.moleculardevices.com/systems/metamorph-research-imaging>`_
* `CPLC software <https://cplc.illinois.edu/software/>`_

Custom **LabVIEW** or vendor packages are fine if they output compatible stacks.
Contact the authors if a format is not read correctly.

Instrumentation and formats
+++++++++++++++++++++++++++++

**Fluorescence channel layout.** Channels are usually **side-by-side** in each
frame (e.g. Cy3 left, Cy5 right). **Sequential** stacks (Cy3 frames then Cy5)
need a new profile in ``cascadeConstants.m``.

**Movie format.** Use monochrome **TIFF** stacks. For large sCMOS movies use
`BigTIFF <http://bigtiff.org/>`_. Suggested filename pattern:
``YYMMDD sample-name experimental-condition.tif``. Prefer **≥ 50** frames,
often **1,000–3,000**. Store exposure time (seconds) in EXIF tag
**ExposureTime** (33434); missing tags break batch workflows.

**Immobilization.** SPARTAN targets **surface-tethered** particles; diffusing
molecules and **drift correction** are not supported.

**Alignment.** Align channels with beads (e.g. TetraSpeck 0.1 µm); minor
residuals are handled automatically.

**Microscope.** Stable optics, even illumination, minimal drift, and a correctly
set objective correction collar.

**Multiple cameras.** Match cameras (e.g. sequential serials), calibrate,
**synchronize** (external trigger), and align axes to limit rolling-shutter
issues on **sCMOS**.

Data-quality tips
+++++++++++++++++++

**Hardware binning.** Binning reduces size and read noise trade-offs; it lowers
maximum immobilization density.

**Immobilization density.** Estimate background particles; aim for **≥ 10×**
signal over background; keep overlap **≤ ~40%**.

**Minimum intensity.** Mean total intensity **≥ ~300** photons (**~20:1**
:math:`\mathrm{SNR}_{\mathrm{bg}}`). Very low :math:`\mathrm{SNR}_{\mathrm{bg}}`
can cause missed picks.

**Fast photobleaching.** ``gettraces`` sums the first **10** frames by default;
very fast bleaching biases detection—reduce power, improve dyes/scavenging, or
lower the summed frame count (e.g. **4**).

**First frame.** Remove artifact-prone first frame(s) in acquisition software.

**Illumination uniformity.** Global thresholding in ``gettraces`` needs **≤
~30%** drop from center to edge.
