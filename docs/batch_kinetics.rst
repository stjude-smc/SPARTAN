.. _batch_kinetics:

Hidden Markov modeling (batchKinetics)
======================================

Single-molecule interpretation usually requires **(i)** choosing a kinetic model
and **(ii)** idealizing each trace (hidden state vs. time). **batchKinetics**
implements HMM-related tools. Related software:

* `QuB <http://www.qub.buffalo.edu/>`_
* `HaMMy <http://bio.physics.illinois.edu/HaMMy.asp>`_
* `ebFRET <http://ebfret.github.io/>`_
* `SMART <https://simtk.org/home/smart>`_
* `FRET community software list <https://fret.community/software/>`_

.. note::

   **Figure 7.** **batchKinetics** screenshot — add ``.. figure::`` when available.

Inputs: traces and a model
----------------------------

Load ``.traces`` with **File → Load traces** (toolbar). The **Files to analyze**
list selects the active movie; traces appear in blue, idealizations in red if a
``.dwt`` exists. Scrollbars navigate traces and time zoom.

**File → Load Model** or **New Model** opens the graphical **QubModel** editor:
states (boxes), classes (colors), rates (labels, in :math:`\mathrm{s}^{-1}`).
Click a rate to fix/optimize it; right-click states for **properties** and
connectivity. Save as ``*.model`` (MATLAB ``load``-able ``QubModel`` object).

Analysis methods
----------------

Choose the **Method** in **Analysis Settings** (and **Edit → Analysis Settings**
for tolerances).

* **SKM** (`Qin, Biophys. J. 2004 <https://doi.org/10.1016/S0006-3495(04)74217-4>`_) —
  fast segmental *k*-means; pair with **MIL** for rates.
* **Baum–Welch** (`Rabiner, Proc. IEEE 1989 <https://doi.org/10.1109/5.18626>`_;
  `Qin *et al.*, Biophys. J. 2000 <https://doi.org/10.1016/S0006-3495(00)76441-1>`_) —
  EM on transition probabilities; combine with **MIL** for physical rates.
* **MPL** — optimizes rates with ``fmincon`` (slower, flexible constraints).
* **MIL** (`Qin *et al.*, Biophys. J. 1996 <https://doi.org/10.1016/S0006-3495(96)79568-1>`_) —
  uses dwell statistics for rate fits (often after SKM/Baum–Welch).
* **ebFRET** (`van de Meent *et al.* <https://doi.org/10.1016/j.bpj.2013.12.055>`_;
  `Bronson *et al.* <https://doi.org/10.1016/j.bpj.2009.09.031>`_) — Bayesian
  model discovery (slow; needs only state count).
* **HMJP** (`Kilic *et al.* <https://doi.org/10.1016/j.bpj.2020.12.022>`_) —
  accounts for frame integration when rates approach the camera rate.

Run the **toolbar** action on the current file. Except **MIL**, outputs save a
``.dwt`` beside the input. Enable **Update Model Parameters** to write optimized
rates back into the model object.

Results and plots
-----------------

After a run, **Plot** / toolbar actions include **dwellhist**, **makeplots**,
**percentTime**, **occtime**, **transitionsPerSecond**, etc.

Simulation
----------

**Action → Simulate Traces** draws continuous-time trajectories (**Gillespie**),
bins them to the camera interval, and adds Gaussian / optional Poisson noise.
Wide-field mode can overlay PSFs on a background movie. Outputs share a base
name: ``*.sim.log``, ``*.sim.dwt``, ``*.traces``, optional ``*.tif``.

Key simulation knobs include **integration time**, **mean photons/frame**,
**SNR**, **shot noise**, **Apparent Gamma**, bleaching options, and (optional)
movie **grid / density / PSF size / alignment** parameters—tune using **autotrace**
statistics from real data when possible.
