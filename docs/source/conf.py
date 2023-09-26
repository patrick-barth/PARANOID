# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sphinx_rtd_theme

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PARANOiD'
copyright = '2022, Patrick Barth'
author = 'Patrick Barth'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx_rtd_theme'
]

templates_path = ['_templates']
exclude_patterns = []
master_doc = 'index'
smartquotes = False

# -- Include images ----------------------------------------------------------

rst_prolog = """
.. |IGV overview wig| image:: images/IGV-wig-overview.png
    :width: 800
    :alt: Showing an overview of cross-link events via IGV
.. |IGV zoom wig| image:: images/IGV-wig-zoomed.png
    :width: 800
    :alt: Zoom into cross-link events via IGV
"""


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
