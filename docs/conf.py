from importlib.metadata import PackageNotFoundError, version as resolve_version
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

project = "Asymmetric Uncertainty"
author = "Caden Gobat"
copyright = "2022, Caden Gobat"
version = resolve_version("asymmetric_uncertainty")

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.intersphinx",
]

autosectionlabel_prefix_document = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

root_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "navbar_start": ["navbar-logo"],
    "navbar_center": [],
    "navbar_end": [
        "theme-switcher",
        "navbar-icon-links",
    ],
    "navbar_persistent": [
        "search-button",
    ],
    "navbar_align": "right",
    "github_url": "https://github.com/cgobat/asymmetric_uncertainty",
    "show_prev_next": False,
    "show_nav_level": 2,
}
html_sidebars = {
    "index": [],
    "**": [
        "sidebar-collapse",
        "sidebar-nav-bs",
    ],
}
