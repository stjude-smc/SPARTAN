.. _appendix_functions:

Appendix A. Full list of MATLAB functions
=========================================

Use ``doc functionName`` in MATLAB for full help. Entries below summarize roles;
gray-listed internals are omitted here—prefer the in-app documentation for rare
hooks.

.. note::

   **Figure 11.** Schematic overview of key functions — add ``.. figure::`` when
   available.

Gettraces
---------

* ``gettraces_gui`` — GUI for movies → traces.
* ``gettraces`` — shortcut to ``gettraces_gui``.
* ``pickPeaks`` — locate PSF maxima.
* ``icpalign`` — iterative closest-point alignment.
* ``getCentroids`` — fast PSF centroids.
* ``subfield`` — extract named subfields (e.g. quadrants).
* ``translatePeaks`` — map donor peaks to expected acceptor positions.
* ``weberQuality`` — alignment contrast score.
* ``Wavelength_to_RGB`` — pseudo-color by wavelength.

Autotrace
---------

* ``autotrace`` — histogram GUI + selection export.
* ``traceStat`` — per-trace statistics.
* ``calcLifetime`` — donor bleach detection (used by ``traceStat``).
* ``medianfilter`` — median smoothing for bleaching features.
* ``rleFilter`` — zero short runs in binary masks.
* ``RLEncode`` — run-length encoding helper.
* ``pickTraces`` — apply selection rules.
* ``loadPickSaveTraces`` — batch load/pick/save.

Trace corrections
-----------------

* ``bgsub`` — zero post-bleach baseline.
* ``crosstalkcorrect`` — donor→acceptor bleed-through.
* ``gammacorrect`` — brightness / detection imbalance.
* ``scaleacceptor`` — manual acceptor scaling.
* ``correctTraces`` — ordered correction pipeline (internal).

Visualization
-------------

* ``sorttraces`` — inspect / correct / bin traces.
* ``forOrigin`` — text export for Origin.
* ``makeplots`` — contour, occupancy, TDP overview.
* ``makecplot``, ``cplot`` — FRET–time contours.
* ``statehist`` — per-state histograms (``makeplots`` middle row).
* ``tdplot``, ``tplot`` — transition density plots.
* ``frethistComparison`` — overlay FRET histograms across files.
* ``avgFretTime`` — mean :math:`E(t)` trends.
* ``occtime`` — state occupancy vs. time.
* ``percentTime`` — fractional dwell from ``.dwt``.
* ``dwellhist`` — lifetime / Siegworth plots.
* ``lifetime_exp`` — exponential fits to ``dwellhist`` data.
* ``transitionsPerSecond`` — global transition rate metric.

Hidden Markov modeling
----------------------

* ``batchKinetics`` — HMM GUI (fit, idealize, simulate).
* ``skm``, ``idealize``, ``forward_viterbi`` — SKM / Viterbi idealization.
* ``forwardBackward`` — forward–backward probabilities.
* ``milOptimize`` — maximum interval likelihood.
* ``mplOptimize`` — maximum point likelihood.
* ``BWoptimize`` — Baum–Welch EM.
* ``runEbFRET`` — external ebFRET driver.
* ``simulate`` — synthetic traces.
* ``simulateMovie`` — synthetic wide-field movies.
* ``loadDWT``, ``saveDWT`` — dwell-time I/O.
* ``dwtToIdl``, ``idlToDwt`` — dwell list ↔ idealization.

Third-party interchange
-----------------------

* ``forHammy`` / ``hammyToDWT`` — HaMMy formats.
* ``forvbFRET`` / ``vbFRET_dwt`` — vbFRET formats.
* ``forQuB`` — QuB text import.
* ``tracesToTxt``, ``tracesToMat`` — generic exports.
* ``cy5ForQuB`` — Cy5 channel helper for QuB.
* ``qub_loadTree``, ``qub_saveTree`` — QuB tree binaries.

Utilities
---------

* ``getFile``, ``getFiles``, ``getFileGroups`` — UI file pickers.
* ``haranfilter`` — anti-correlation-preserving smoother.
* ``combineDatasets`` — merge trace files.
* ``loadTraces``, ``saveTraces`` — binary trace I/O.
* ``sizeTraces`` — probe trace counts/lengths quickly.
* ``resizeTraces`` — truncate to common length.
* ``rebintraces`` — downsample traces (legacy helper).
* ``avgfret`` — file-wide average FRET.

Classes
-------

* ``Movie``, ``Movie_STK``, ``Movie_TIFF`` — movie readers.
* ``MovieParser`` — detection + trace export pipeline.
* ``Traces``, ``TracesFret``, ``TracesFret4`` — trace containers.
* ``TraceListViewer`` — scrolling trace UI (``batchKinetics``).
* ``QubModel``, ``QubModelViewer`` — kinetic model object + editor.

Pre-steady-state / RTD
----------------------

* ``rtdgui`` — RTD pipeline UI.
* ``rtdPlots`` — replot saved RTD outputs.
* ``rtdTool`` — scripted backend.
* ``postSyncTraces`` — model-based post-synchronization.
* ``simplePostSync`` — threshold post-sync.
* ``stateMin`` — productive vs. non-productive classifier.
* ``stateOccupancy`` — state occupancy vs. time.
