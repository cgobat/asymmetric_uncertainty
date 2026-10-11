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
]

autosectionlabel_prefix_document = True
root_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "navigation_depth": 4,
}
