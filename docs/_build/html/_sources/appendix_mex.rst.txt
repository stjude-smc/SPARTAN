.. _appendix_mex:

Appendix C. Compiled functions (.mex files)
===========================================

Some routines ship as **MEX** binaries (C/C++) for speed. They are
**platform-specific**; sources and libs live under ``binary/`` in the
distribution. Except ``qub_loadTree`` / ``qub_saveTree``, each compiled entry has
a pure-MATLAB fallback named without the trailing ``x``.

* ``qub_loadTree``, ``qub_saveTree`` — QuB **qubtree** I/O (``.qmf``). Sources:
  ``treestruct.cpp`` plus QuB drop-ins in ``qubtree/``;
  `QuB sources <https://qub.mandelics.com/sources.html>`_.

* ``forward_viterbix`` — Viterbi idealization (`Rabiner, 1989
  <https://doi.org/10.1109/5.18626>`_); used by ``idealize`` / ``skm``. Code under
  ``HMM/``.

* ``gillespiex`` — Gillespie SSA for ``simulate``.

* ``forwardBackwardx`` — Forward–backward probabilities for Baum–Welch / MPL.
