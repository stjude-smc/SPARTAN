.. _appendix_file_formats:

Appendix B. File formats
========================

Image stacks (movies)
---------------------

Save wide-field movies as **TIFF** stacks; for large **sCMOS** datasets prefer
`BigTIFF <https://www.loc.gov/preservation/digital/formats/fdd/fdd000328.shtml>`_.

Binary traces (``.traces`` / ``.rawtraces``)
--------------------------------------------

Binary layout for fast I/O; load in MATLAB via ``loadTraces``. Pseudocode:

.. code-block:: text

   struct TracesFile {
     uint32 zero = 0
     char[4] signature = "TRCS"
     uint16 version = 5
     uint8 dataType = enum( 9 = single )
     uint8 nChannels = C
     uint32 nTraces = N
     uint32 nFrames = M
     uint16 chNameLen = ...
     char[*] chNames = channel names, delimited by char(31)
     dataType[M] time = acquisition times (ms)
     dataType[N×M] channel1 = donor
     dataType[N×M] channel2 = acceptor
     dataType[N×M] channel3 = FRET
     MetadataField root = fileMetadata + traceMetadata tree
   }

``dataType`` encodes MATLAB classes (``char``, integers, ``single``, ``double``,
``logical``, ``cell``, ``struct``, …). Metadata pages follow:

.. code-block:: text

   struct MetadataField {
     uint8 dataType
     uint32 fieldSize
     uint8 nameLen
     char[*] fieldName
     uint8 ndim
     uint32 dataSize[...]
     dataType[...] contents
   }

Packed structs repeat **fieldName / isPacked / MetadataField** children. Char
arrays concatenate with ``char(31)`` delimiters. **fileMetadata** holds
experiment-wide fields (e.g. wavelength list ``[532 640]``); **traceMetadata**
stores per-trace entries (IDs, ``donor_x``, ``donor_y``, …).

Dwell-time text (``.dwt``)
--------------------------

QuB-compatible idealizations. Header example::

 Segment: 1 Dwells: 5 Sampling(ms): 25 Start(ms): 0 ClassCount: N µ1 σ1 … µN σN

Each dwell line: **class index** (file uses **0-based** classes; MATLAB loaders
may shift by **1**) and **duration (ms)**. Offsets in **Start(ms)** align idealizations
with FRET arrays. Specification:
`QuB DWT format <https://qub.mandelics.com/m/DWT.html>`_.

Kinetic models (``.model``)
----------------------------

Version ≥3 saves ``QubModel`` objects as ``.mat`` files with a ``.model``
extension. Format reference:
`MAT-file format (MathWorks) <https://www.mathworks.com/help/pdf_doc/matlab/matfile_format.pdf>`_.

QuB import text (``.qub.txt``)
------------------------------

Concatenated FRET samples for QuB import; segment boundaries are assigned inside
QuB.

QuB_Tree (``.qmf``)
--------------------

Legacy QuB model trees are **deprecated** as of SPARTAN **3.10**. Documentation:
`QuB manual <https://qub.mandelics.com/m/qubdoc/>`_.
