#!/usr/bin/env python3
"""
Import SPARTAN simulation guide.docx into a single simulation_guide.rst.

Run from the docs/ directory:

    python3 import_docx_kinetic.py

Reuses parsers from import_docx.py (same folder).
"""
from __future__ import annotations

import os

import import_docx as idoc

DOCX_NAME = "SPARTAN simulation guide.docx"
OUTPUT_RST = "simulation_guide.rst"
PAGE_TITLE = "SPARTAN Simulation Guide"


def _flatten_segments(
    segments: list[tuple[str, list]],
    page_title: str,
) -> list[tuple[str | None, str, int | None]]:
    """One RST page: each Word H1 becomes a section; first H1 uses canonical title."""
    flat: list[tuple[str | None, str, int | None]] = []
    for i, (title, body) in enumerate(segments):
        h1 = page_title if i == 0 else title
        flat.append((None, h1, 1))
        flat.extend(body)
    return flat


def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    docx = os.path.normpath(
        os.path.join(
            here,
            "..",
            "arxiv",
            "SPARTAN Documentation 20240416",
            DOCX_NAME,
        )
    )
    if not os.path.isfile(docx):
        raise SystemExit(f"Missing docx: {docx}")

    blocks = idoc.parse_docx(docx)
    segments = idoc.merge_duplicate_h1_segments(idoc.split_segments(blocks))
    if not segments:
        raise SystemExit("No content segments found in docx.")

    flat = _flatten_segments(segments, PAGE_TITLE)
    rst = idoc.blocks_to_rst(flat)
    out_path = os.path.join(here, OUTPUT_RST)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(rst)
    print("Wrote", out_path)


if __name__ == "__main__":
    main()
