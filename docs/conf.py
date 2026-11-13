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
trim_footnote_reference_space = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

root_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_favicon = "assets/favicon.png"
html_theme = "alabaster"
html_theme_options = {
    "github_user": "cgobat",
    "github_repo": "asymmetric_uncertainty",
    "github_button": True,
    "github_type": "star",
    "page_width": "80%",
}
html_sidebars = {
    "**": [
        "about.html",
        "navigation.html",
        "relations.html",
        "searchfield.html",
    ],
}
