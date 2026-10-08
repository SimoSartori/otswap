"""Sphinx configuration for the otswap documentation.

The C++ reference comes from Doxygen XML through Breathe, the Python reference
from the type stubs through sphinx-autoapi. The stubs carry the full
docstrings and are what an IDE shows, and reading them needs no compiled
module, so the documentation builds without building otswap.

    pip install -r docs/requirements.txt
    sphinx-build -W --keep-going docs docs/_build/html
"""

import pathlib
import subprocess
import tomllib

HERE = pathlib.Path(__file__).parent
ROOT = HERE.parent

project = "otswap"
author = "Simone Sartori"
copyright = "2026, Simone Sartori"
release = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
version = release

# Doxygen runs as part of the Sphinx build, so there is one command to remember
# and no way to publish a C++ reference built from a stale header.
subprocess.run(["doxygen", "Doxyfile"], cwd=HERE, check=True)

extensions = [
    "autoapi.extension",
    "breathe",
    "myst_parser",
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx_design",
]

breathe_projects = {"otswap": str(HERE / "doxygen" / "xml")}
breathe_default_project = "otswap"
breathe_domain_by_extension = {"h": "cpp"}

# The stubs only: each .pyi is read in place of the .py of the same module.
# No pages are generated; python.rst places every entry by hand.
autoapi_type = "python"
autoapi_dirs = [str(ROOT / "python" / "otswap")]
autoapi_file_patterns = ["*.pyi", "*.py"]
autoapi_generate_api_docs = False
autoapi_add_toctree_entry = False
autoapi_options = ["members", "undoc-members", "show-inheritance"]
autoapi_member_order = "bysource"

napoleon_numpy_docstring = True
napoleon_google_docstring = False

# Files and folders of the repository, linked on GitHub: the example data are
# not copied into the site.
REPOSITORY = "https://github.com/SimoSartori/otswap"
extlinks = {
    "repo-file": (f"{REPOSITORY}/blob/main/%s", "%s"),
    "repo-tree": (f"{REPOSITORY}/tree/main/%s", "%s"),
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
}

myst_enable_extensions = ["amsmath", "attrs_block", "colon_fence", "dollarmath"]
myst_heading_anchors = 3

templates_path = []
exclude_patterns = ["_build", "doxygen", "requirements.txt"]

html_theme = "furo"
html_title = f"otswap {release}"
html_static_path = []
