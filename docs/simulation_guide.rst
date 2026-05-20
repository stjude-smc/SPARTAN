.. _simulation_guide:

SPARTAN Simulation Guide
========================

**Step-by-step guide to simulating data with SPARTAN**

.. rubric:: Author

**Daniel Terry** — St. Jude Children’s Research Hospital, Memphis, TN

**March 31, 2020**

Summary
-------

This guide explains how to simulate single-molecule fluorescence traces in SPARTAN. Simulations test whether a kinetic model is **physically plausible** and **consistent with data**. Prefer analytic theory when your parameters allow it; use stochastic simulation when you need explicit trajectories and noise.

Software and background: `Scott C. Blanchard lab software <https://www.scottcblanchardlab.com/software>`_.

.. note::

   This guide assumes SPARTAN is installed and you are familiar with **batchKinetics** basics.

Determine model observation parameters (FRET states)
----------------------------------------------------

Each conformational state maps to an **apparent FRET** value :math:`E`.

* **From data:** Gaussian fits to experimental FRET histograms give mean :math:`E` per state.
* **From structure:** estimate dye–dye distance :math:`r` (Å) and convert to :math:`E` using Förster theory with :math:`R_0` for your dye pair.

.. math::

   E = \frac{1}{1 + (r/R_0)^6}

Here :math:`E` is apparent FRET efficiency, :math:`r` is inter-dye distance, and :math:`R_0` is the distance where :math:`E = 0.5`. For **LD555–LD655**, :math:`R_0 \approx 58\ \mathrm{Å}` is a reasonable first guess; refine using spectra and quantum yields on your construct (*Principles of Fluorescence Spectroscopy*, Lakowicz).

Structural :math:`E` values are **approximate** (linkers, dye orientation, etc.). Implicit-solvent models (e.g. `SMOG <http://smog-server.org/>`_) can crudely probe linker effects but remain approximate.

Determine model kinetic parameters
----------------------------------

Specify **rate constants** for allowed transitions between states.

.. important::

   **Degenerate states:** several *physical* states may share one **apparent FRET** (“class”). They differ in kinetics but look identical in FRET—broad or multi-modal dwell histograms are a clue. **Classes** must be ordered by increasing FRET in **batchKinetics**; the first state (black background) is the **donor dark** state (:math:`E=0`).

Set **initial state probabilities** :math:`p_0` (sum to **1**). Equilibrium models can auto-compute :math:`p_0`; non-equilibrium experiments need priors or measured occupancies.

Building a model in SPARTAN
----------------------------

Open ``batchKinetics`` in MATLAB.

* Click **New Model** (toolbar). States appear as **boxes** (index + **class** color = apparent FRET). **Double-click** a state to edit means / classes (simulation ignores **stdev** and **fix** flags for simulation-only use).
* **Right-click** empty space → **New State** to add states.
* **Arrows** show allowed transitions; numeric labels are **rate constants** (:math:`\mathrm{s}^{-1}`). **Right-click** a state → **Connect to…** to wire the model. Click a rate to edit; set **both** rates to **0** to remove a link.
* Ensure :math:`p_0` sums to **1**. For equilibrium, **right-click** background → **calculate and apply equilibrium p0**. **Save Model** when done.

Simulate smFRET traces
----------------------

Click **Simulate Data** (toolbar). Key fields include:

* **Traces** — number of independent trajectories.
* **Frames** — trace length (often **1000–3000**).
* **Integration time (ms)** — camera exposure (**~10–100** ms typical).
* **Signal-to-background SNR** — match **autotrace** when reproducing experiments.
* **Intensity** — mean total photons (donor + acceptor).
* **Intensity stdev** — optional illumination non-uniformity.
* **Donor / acceptor lifetime** — underlying photobleaching timescales (may differ from apparent lifetimes in **autotrace**).

Save a screenshot of settings, click **Ok**, and choose an output filename.

Evaluating simulated traces
---------------------------

**batchKinetics** does not archive intermediate dwell lists. Idealize simulated ``.traces`` (e.g. **SKM**) using the same or a fresh model.

* **percentTime** (toolbar) — steady-state occupancies; **Edit → Copy Values** for numbers.
* **Dwell-time** tools — compare model-implied kinetics to experiment.

Because many noise sources are omitted, simulated FRET histograms are often **narrower** than experiment (optics, cameras, blinking, background, etc.).

Tips and caveats (blinking)
---------------------------

Blinking is a **separate** Markov process from conformational dynamics; a **single** HMM rarely captures both faithfully.

A **single** dark state wired to all bright states implies **full randomization** via dark visits—often wrong for **intramolecular** FRET. Prefer **one dark state per bright FRET state** if you must mimic fast acceptor dark events; avoid simulating blinking unless it is the scientific target.

.. note::

   **Photon-emission series (advanced)** — placeholder for a future section.

Algorithm description
---------------------

**Continuous-time trajectories** use the `Gillespie algorithm <https://en.wikipedia.org/wiki/Gillespie_algorithm>`_:

#. Sample the starting state from :math:`p_0`.
#. For each exit, draw exponential dwells :math:`\tau_{ij} = 1/k_{ij}` and take the **shortest** event; advance time and state; repeat until the trace length is reached.

**Time-averaging:** ideal FRET traces are sampled at the camera interval; fast kinetics can **blur** multiple dwells into one frame.

**Noise model:** ideal intensities use :math:`I_A = E \cdot \mathrm{MTI}`, :math:`I_D = (1-E)\cdot \mathrm{MTI}` (with optional **gamma** scaling on donor). **Poisson** shot noise is drawn per frame, then **Gaussian** background noise is added to reach the target SNR.

.. warning::

   The simulator does **not** include orientation effects, complex photophysics, slow baseline drift, or per-trace optical variation—use synthetic data as a controlled approximation, not a full experimental duplicate.
