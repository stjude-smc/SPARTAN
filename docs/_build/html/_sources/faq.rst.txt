.. _faq:

Troubleshooting FAQ
==========================

.. _faq-emccd:

How can EMCCD data be analyzed in SPARTAN?
------------------------------------------

SPARTAN does not require a separate EMCCD pipeline. Analysis is the same as for
sCMOS or other cameras: wide-field **TIFF** movies → **gettraces** → **autotrace**
→ **sorttraces** / **makeplots** → **batchKinetics**, as in :doc:`spartan_workflow`.

**Export movies.** Save monochrome TIFF stacks from your acquisition software

**Channel layout.** Adjust ``photonsPerCount`` considering what the EM gain and gain settings are,
we right click on the gettraces on each channel to adjust what the photons per count is.

**Workflow entry.** Open **spartan** → **Extract fluorescence traces from movies**,
choose your profile, open the ``.tif`` stack, align channels (beads if needed),
**Save Traces**, then continue with :doc:`quickstart` / :doc:`autotrace`.

If your EMCCD data use a non-standard layout (e.g. separate TIFF files per channel),
note that the repo already supports **dual TIFF merge** via profile
``Dual TIFF merge (L-R, 2-color)`` in ``cascadeConstants.m``—mention that in the FAQ
if it matches your lab’s export pattern.