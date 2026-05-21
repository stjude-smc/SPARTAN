SPARTAN documentation (Sphinx)
==============================

Local setup (once)
------------------
  python3 -m venv .venv
  source .venv/bin/activate          # Windows: .venv\Scripts\activate
  pip install -r requirements.txt

Build HTML
----------
  make html
  # or: sphinx-build -b html . _build/html

Preview locally (recommended)
-----------------------------
  python3 -m http.server --directory _build/html
  Open http://127.0.0.1:8000 in a browser.

Regenerate narrative text from Word
------------------------------------
  python3 import_docx.py
  (Reads ../arxiv/SPARTAN Documentation 20240416/SPARTAN Documentation.docx)

  python3 import_docx_quickstart.py
  (Reads .../SPARTAN Quick Start Guide.docx → quickstart.rst)

  python3 import_docx_kinetic.py
  (Reads .../SPARTAN simulation guide.docx → simulation_guide.rst)

  Each run overwrites the target .rst (commit or branch before regenerating if
  you have hand edits).

Notes
-----
  Figures embedded only in Word are not extracted automatically; place assets in
  _static/ and reference them with .. image:: in the .rst files.

  Equations use MathJax (may load from the network in the browser unless you
  configure a local mathjax_path in conf.py).

RST styling (BIASD-style patterns)
----------------------------------
  - Page label: .. _my_page:
  - Sections: underline with ===, ---, +++ (see arxiv/docs/*.rst)
  - Links: `Title <https://example.com>`_
  - Code names: ``functionName``, ``file.m``
  - Code blocks: .. code-block:: matlab (or bash, text, python)
  - Math: :math:`E` inline or .. math:: block
  - Papers: :Title: / :Journal: / :DOI: field lists (see references.rst)
  - Callouts: .. note::, .. warning::, .. seealso::

Folder structure (what each part does)
--------------------------------------
  index.rst
    Master document: home page text plus the main .. toctree:: that defines
    which pages exist and in what order (sidebar navigation).

  conf.py
    Sphinx configuration: project name, theme, extensions (e.g. mathjax),
    html_static_path, etc.

  *.rst
    One file per page (except index.rst, which is special). The filename (without
    .rst) is what you list in the toctree.

  _static/
    Images, PDFs, etc. referenced from .rst (e.g. .. figure:: /_static/foo.png).

  _build/
    Generated HTML; safe to delete; recreated by "make html". Listed in .gitignore.

  requirements.txt
    Python packages for Sphinx (local venv or Read the Docs).

  Makefile / make.bat
    Shortcuts; "make html" runs sphinx-build.

  import_docx.py
    Optional: regenerates a fixed set of section .rst files from the Word doc.
    It does not discover files automatically; only index.rst plus your .rst files
    define the full site map.

  README.txt
    This file (not part of the built HTML unless you add it to a toctree).

Order of pages in the sidebar (important)
-------------------------------------------
  Sphinx does NOT sort .rst files alphabetically and does NOT put new pages
  "last" by default. The sidebar order is exactly the order of entries under
  .. toctree:: in index.rst (top to bottom). To put a new page first, add its
  name as the first line under the toctree; to put it between two existing
  pages, insert that line between those two names.

  You can have more than one .. toctree:: on index.rst (each can have its own
  :caption:), e.g. "User guide" vs "Developer notes". Each toctree lists its
  pages in the order you want for that group.

Adding a new sidebar entry (a new documentation page)
-----------------------------------------------------
  The Read the Docs theme lists each top-level toctree entry in the left
  sidebar (sometimes called a "tab"). To add a page for a new feature:

  1. Create a new file under this folder, e.g. my_feature.rst

  2. Give it a title and optional reference label at the top, for example:

       .. _my_feature:

       My feature
       ==========

       Intro text...

  3. Open index.rst and add the document name (no .rst suffix) to the
     .. toctree:: block, in the position you want in the sidebar (see "Order
     of pages" above):

       .. toctree::
          :maxdepth: 2
          :caption: User guide

          introduction
          ...
          my_feature

  4. Rebuild: make html

  Link from another page with :doc:`my_feature` or :ref:`my_feature` (if you set .. _my_feature:). Subsections inside the .rst file appear nested under
  that entry when :maxdepth: is 2 or higher.

  Do not edit import_docx.py for hand-written pages; import_docx only overwrites
  the section files it knows about. Keep standalone guides as separate .rst
  files and list them only in index.rst.


Publishing on GitHub Pages 
---------------------------------------
Live preview (after workflow + Pages are enabled on the branch):
  https://kiliczeliha.github.io/SPARTAN/

Sources live in this docs/ folder. HTML under docs/_build/ is generated locally
and in CI; do not commit docs/_build/ (see repo .gitignore).

Updates: push to branch my-feature on your fork; Actions rebuilds the site.
Manual rebuild: GitHub → Actions → "Deploy Documentation to GitHub Pages" → Run workflow.

Unpublish or freeze:
  - Take site offline: fork Settings → Pages → Unpublish.
  - Stop automatic updates: edit .github/workflows/deploy-documentation.yml
    (e.g. workflow_dispatch only) or disable/delete that workflow.