"""Generate docs pages from README.md and the notebooks.

Writes docs/index.md and docs/notebooks/*.md (with extracted output images).
Notebooks are not executed: the stored outputs are rendered as-is.
Run before `zensical build`.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import nbformat
from nbconvert import MarkdownExporter

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
NB_SRC = ROOT / "notebooks"
NB_OUT = DOCS / "notebooks"
REPO_BLOB = "https://github.com/ToledoEM/PubMatrixPython/blob/main/"


def build_index() -> None:
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    # docs/foo.md -> foo.md (docs_dir is the site root)
    text = re.sub(r"\]\(docs/([^)]+)\)", r"](\1)", text)
    # Files outside docs_dir -> link to GitHub
    text = re.sub(r"\]\((LICENSE\.md)\)", rf"]({REPO_BLOB}\1)", text)
    (DOCS / "index.md").write_text(text, encoding="utf-8")


def build_notebooks() -> None:
    if NB_OUT.exists():
        shutil.rmtree(NB_OUT)
    NB_OUT.mkdir(parents=True)
    exporter = MarkdownExporter()
    for nb_path in sorted(NB_SRC.glob("*.ipynb")):
        nb = nbformat.read(nb_path, as_version=4)
        files_dir = f"{nb_path.stem}_files"
        body, resources = exporter.from_notebook_node(
            nb,
            resources={"output_files_dir": files_dir, "unique_key": nb_path.stem},
        )
        for rel, data in resources.get("outputs", {}).items():
            out = NB_OUT / rel
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_bytes(data)
        (NB_OUT / f"{nb_path.stem}.md").write_text(body, encoding="utf-8")


def main() -> None:
    build_index()
    build_notebooks()
    print(f"Generated {DOCS / 'index.md'} and notebook pages in {NB_OUT}")


if __name__ == "__main__":
    main()
