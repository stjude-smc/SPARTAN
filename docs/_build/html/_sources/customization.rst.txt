.. _customization:

Customization
=============

Beyond the GUIs, nearly every routine is callable from scripts. Start from
``examples.m`` in the SPARTAN source: each block demonstrates command-line
patterns with inline comments.

Central configuration: ``cascadeConstants.m``
---------------------------------------------

Most tunables live in **one** file for traceability. Sections include:

**Global**
   Shared defaults across tools.

**Gettraces**
   Imaging profiles and GUI-hidden defaults; add profiles carefully so existing
   entries keep working.

**Autotrace**
   Default selection rules and **traceStat** parameters.

**Makeplots**
   Contour length, styling, and options also used by helpers like
   ``frethistComparison``.

**Other**
   Miscellaneous flags—including disabling **parallel** pools.

Some niche behaviors still require edits inside individual functions.

Extending **autotrace**
-----------------------

Add fields to the **traceStat** output struct (and its description struct); the
GUI picks them up automatically.
