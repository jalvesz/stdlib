from __future__ import annotations

import subprocess
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

html_theme = 'furo'
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
GENERATED_FORTRAN_DIR = DOCS_DIR / '_generated' / 'fortran'
FORTRAN_PILOT_SOURCES = {
    REPO_ROOT / 'src' / 'core' / 'stdlib_kinds.fypp': GENERATED_FORTRAN_DIR / 'core' / 'stdlib_kinds.f90',
    REPO_ROOT / 'src' / 'strings' / 'stdlib_string_type.fypp': GENERATED_FORTRAN_DIR / 'strings' / 'stdlib_string_type.f90',
    REPO_ROOT / 'src' / 'strings' / 'stdlib_strings.fypp': GENERATED_FORTRAN_DIR / 'strings' / 'stdlib_strings.f90',
}


def preprocess_fortran_sources() -> None:
    version = (REPO_ROOT / 'VERSION').read_text().strip().split('.')
    GENERATED_FORTRAN_DIR.mkdir(parents=True, exist_ok=True)
    for source, destination in FORTRAN_PILOT_SOURCES.items():
        destination.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                'fypp',
                '-I',
                str(REPO_ROOT / 'include'),
                f'-DPROJECT_VERSION_MAJOR={version[0]}',
                f'-DPROJECT_VERSION_MINOR={version[1]}',
                f'-DPROJECT_VERSION_PATCH={version[2]}',
                str(source),
                str(destination),
            ],
            check=True,
            cwd=REPO_ROOT,
        )


preprocess_fortran_sources()

fortran_lexer = 'regex'
fortran_doc_chars = ['>', '!']
fortran_sources = [
    str(path.relative_to(DOCS_DIR))
    for path in FORTRAN_PILOT_SOURCES.values()
]
fortran_sources_exclude = []
