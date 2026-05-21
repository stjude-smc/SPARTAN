# -*- coding: utf-8 -*-
# SPARTAN documentation configuration for Sphinx.

import os

extensions = [
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'

project = 'SPARTAN'
copyright = 'SPARTAN documentation authors'
author = 'SPARTAN'

version = '1.0'
release = '1.0'

language = 'en'

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'README.txt']

pygments_style = 'sphinx'
todo_include_todos = False

html_theme = 'sphinx_rtd_theme'
html_title = 'SPARTAN documentation'
html_short_title = 'SPARTAN'
html_static_path = ['_static']

htmlhelp_basename = 'SPARTANdoc'
html_baseurl = "https://kiliczeliha.github.io/SPARTAN/"


latex_documents = [
    (master_doc, 'SPARTAN.tex', 'SPARTAN Documentation',
     author, 'manual'),
]

man_pages = [
    (master_doc, 'spartan', 'SPARTAN Documentation',
     [author], 1)
]

texinfo_documents = [
    (master_doc, 'SPARTAN', 'SPARTAN Documentation',
     author, 'SPARTAN', 'Single-molecule analysis software.',
     'Miscellaneous'),
]
