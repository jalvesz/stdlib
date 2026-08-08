from __future__ import annotations

import subprocess
import sys
from pathlib import Path

project = 'Fortran-lang stdlib'
copyright = '2026, fortran-lang contributors'
author = 'fortran-lang contributors'
release = 'experimental'

extensions = [
    'myst_parser',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'sphinx_fortran_domain',
]

source_suffix = {
    '.md': 'markdown',
    '.rst': 'restructuredtext',
}

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'pydata_sphinx_theme'
html_title = 'Fortran-lang stdlib'
html_static_path = ['media']
html_favicon = 'media/favicon.ico'
suppress_warnings = [
    'misc.highlighting_failure',
    'myst.domains',
    'myst.header',
    'myst.xref_missing',
]
todo_include_todos = True

myst_heading_anchors = 3
myst_enable_extensions = ['colon_fence', 'deflist', 'substitution']

DOCS_DIR = Path(__file__).resolve().parent
REPO_ROOT = DOCS_DIR.parent
STDLIB_FPM_SRC_DIR = REPO_ROOT / 'stdlib-fpm' / 'src'
FORTRAN_PILOT_SOURCES = [
    STDLIB_FPM_SRC_DIR / 'stdlib_kinds.f90',
    STDLIB_FPM_SRC_DIR / 'stdlib_string_type.f90',
    STDLIB_FPM_SRC_DIR / 'stdlib_strings.f90',
]


def preprocess_fortran_sources() -> None:
    subprocess.run(
        [sys.executable, 'config/fypp_deployment.py', '--deploy_stdlib_fpm'],
        check=True,
        cwd=REPO_ROOT,
    )


preprocess_fortran_sources()

fortran_lexer = 'ford'
fortran_doc_chars = ['>', '!']
fortran_sources = [
    str(path)
    for path in FORTRAN_PILOT_SOURCES
]
fortran_sources_exclude = []
