import os
import sys
sys.path.insert(0, os.path.abspath('../src'))

# Project information
project = 'spindynam'
copyright = '2026, Antigravity and the spindynam developers'
author = 'Antigravity'
release = '0.1.0'

# General configuration
extensions = [
    'sphinx_immaterial',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'myst_parser',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# HTML output options
html_static_path = ['_static']

# --- Overridden by user preference ---
html_theme = 'sphinx_immaterial'

html_theme_options = {
    "icon": {
        "repo": "fontawesome/brands/github",
    },
    "site_url": "https://github.com/ddmms/spin_dybamics",
    "repo_url": "https://github.com/ddmms/spin_dybamics",
    "repo_name": "spindynam",
    "font": {
        "text": "Roboto",
        "code": "Roboto Mono",
    },
    "palette": [
        {
            "media": "(prefers-color-scheme: light)",
            "scheme": "default",
            "primary": "indigo",
            "accent": "light-blue",
            "toggle": {
                "icon": "material/lightbulb-outline",
                "name": "Switch to dark mode",
            },
        },
        {
            "media": "(prefers-color-scheme: dark)",
            "scheme": "slate",
            "primary": "indigo",
            "accent": "light-blue",
            "toggle": {
                "icon": "material/lightbulb",
                "name": "Switch to light mode",
            },
        },
    ],
    "features": [
        "navigation.expand",
        "navigation.sections",
        "navigation.top",
        "search.share",
        "toc.sticky",
    ],
}

# Required for Material theme indexing stability
html_last_updated_fmt = ""

# --- Overridden by user preference ---

# Napoleon settings for NumPy/Google style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_use_mathjax = True
