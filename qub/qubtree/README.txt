/* Copyright 1998-2011 Research Foundation State University of New York */

/* This file is part of QuB.                                            */

/* QuB is free software; you can redistribute it and/or modify          */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */

/* QuB is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */

/* You should have received a copy of the GNU General Public License,   */
/* named LICENSE.txt, in the QuB program directory.  If not, see        */
/* <http://www.gnu.org/licenses/>.                                      */

http://www.qub.buffalo.edu/

QuB is a graphical environment for hidden-Markov analysis of time series, especially ion channel recordings.
This package contains the source and package files needed to compile QuB.  You will also need:

 * Delphi 5 -- it may or may not work in later versions of Delphi, with minor revision
 * Python 2.5, pythonwin (with COM extensions), Python4Delphi
 * Microsoft Visual Studio 2003 -- another compiler may work for libraries, but please read the notes below
 * boost++ 1_45_0, specifically boost/graph and boost/ublas
 * NSIS -- if you want to build installers
 * msvcr60.dll, msvcrp71.dll and msvcrt71.dll (the MSVC runtime for qub*.dll),
 * cc3260.dll and cc3260mt.dll (the Borland C runtime, for CVODE.dll).

QuB's interface is written in Delphi pascal.  Some algorithms and file handling are done in C/C++ libraries.
Data acquisition uses plugin dlls, some of which are Delphi and some C/C++.  A few features are written
partially in Python, and limited Python scripting is available.


== The libraries and how they fit together ==

There are several mentions of qubtree and its python bindings (_qubtree)
because in addition to filing, it helps pass structured data between the three languages.

QuB.exe
|- cvode.dll
|- python25.dll
|  |- _qubfits.dll
|  |- _qubtree.dll
|  |- qubopt.dll
|- qubacqdx.dll
|- qubacqni.dll
|- qubacqnm.dll
|- qubacqsc.dll
|- qubfits.dll
|  |- _qubfits.dll
|  |  |- _qubtree.dll
|  |- qubtree.dll
|  |- qubfits plugin.dll
|- qublib.dll
|  |- qubtree.dll
|- qubopt.dll
|  |- qubtree.dll
|  |- _qubtree.dll
|  |- qublib.dll
|  |- maxill.dll
|- qubtree.dll


== Contents of the folders ==

  CVODE              (b) third-party source for cvode.dll (Dynamic-clamp differential systems)
  include            (i) C/C++ header files common to some qub* libraries
  Lib                (a) Utilities for QUB.exe (Delphi source)
  maxill             (c) MIL algorithm extended for varying stimulus
  QUB                (d) QUB.exe Delphi source and project files
  QUBQCQNM           (c) QUB.exe acquisition plugin (dll) for National Instruments MX series
  qubcommon          (i) C/C++ source common to some qub* libraries
  QUBDAQ             (d) QUB.exe acquisition plugins for Sound Card, National Instruments, Data Translation
  qubdoc             (h) QUB/python and qubtree documentation
  qubfits            (c) Curve fitting dll
  _qubfits           (c) qubfits python bindings
  qubfits plugin     (c) Sample plugin for qubfits, with additional curves
  qublib             (c) QUB.exe analysis DLL (MPL, AMP, filtering)
  qubopt             (c) QUB.exe analysis DLL (SKM, MIL, utilities)
  QuBSuite           (p) QUB.exe, dll's, installers, python scripts, presets, samples
  qubtree            (c) qubtree file manipulation
  _qubtree           (c) python bindings for qubtree
  QuB Tree Editor    (w) GUI for editing qubtree files (py2exe/wxPython)
  tpmath1            (a) third-party source for Turbo Pascal Math (numerical recipes)

(a) -- files that are included in QUB.exe
(b) -- library compiled with Borland C
(c) -- library compiled with MS Visual C++
(d) -- Delphi project and source
(h) -- html documentation
(i) -- files that are included in (c) libraries
(p) -- Program folder, with Python source
(w) -- wxPython using qubtree and _qubtree, packaged with py2exe, dead-end.


== Notes on the compiled libraries ==

QuB and its libraries output object files to c:\out and subdirectories.  You'll have to create them
or edit the object paths.

qublib, qubopt and maxill use code from boost++.  You may have to edit
include and library paths to locate your installation.

Some dlls expect python 2.5 to be installed at C:\python25.  If elsewhere, you'll need to edit
the library and include paths for _qubfits, _qubtree, qubopt, and customcurve.

Most of the provided project files output .dll and .lib to folder QuBSuite.  They are configured
to look for each other's .lib files there.

The following libraries all use the "Multithreaded DLL" runtime library (/MD).  If you use a different
compiler, they should all share a runtime library:
  qubtree, _qubtree, qubfits, _qubfits, qubopt, qublib


== qubtree ==

qubtree is pervasive.  It was originally conceived as an extensible file structure, a tree with
typed data at each named node.  The original Delphi implementation "QFS" handles models and data,
with extensive support for editing and undo.

The C/C++ implementation "qubtree" introduced automatic reference-counting and thread-safety.  It presents
a bare-bones C API, with object-oriented wrappers for C++, Delphi, and Python.  With the reference counting
you can pass structures and return detailed results across languages without lifecycle issues.  Consequently,
preferences, models, data, and nearly everything else touched by multiple languages has a qubtree representation.
Many of these are documented in qubdoc\trees.

The Delphi wrapper is QUB\QUB_Node.pas.  It uses Delphi's COM facilities to automate reference counts and type magic.
The Python wrapper is its own dll, _qubtree.dll.  with Python helper file qubtree.py.

