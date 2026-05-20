#!/usr/bin/env python3
"""
Import narrative text from SPARTAN Documentation.docx into section .rst files.

Requires Python 3.6+ and no extra packages. Run from repo root or docs/:

    python3 import_docx.py

Source path is resolved relative to this file.
"""
from __future__ import annotations

import os
import re
import zipfile
from xml.etree import ElementTree as ET

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

# Canonical page titles (sidebar / H1) — match user TOC
RST_TITLES = {
    "introduction.rst": "Getting Started",
    "installation_compiled.rst": "Installation of the compiled version",
    "installation_source.rst": "Installation from source code",
    "spartan_workflow.rst": "Main menu and core analysis workflow (spartan)",
    "gettraces.rst": "Extracting fluorescence traces from movies (gettraces)",
    "autotrace.rst": "Select traces according to defined criteria (autotrace)",
    "sorttraces.rst": (
        "VII. Viewing traces, manual corrections, and trace selection (sorttraces)"
    ),
    "batch_kinetics.rst": "Hidden Markov modeling (batchKinetics)",
    "rtd_pipeline.rst": (
        "X. Analysis pipeline for pre-steady-state experiments "
        "(rtdgui, rtdTool, and rtdPlots)"
    ),
    "customization.rst": "Customization",
    "appendix_functions.rst": "Appendix A. Full list of MATLAB functions",
    "appendix_file_formats.rst": "Appendix B. File formats",
    "appendix_mex.rst": "Appendix C. Compiled functions (.mex files)",
    "references.rst": "References",
}

RST_FILES_ORDER = list(RST_TITLES.keys())


def _para_text(p: ET.Element) -> str:
    parts: list[str] = []
    for t in p.findall(f".//{W}t"):
        if t.text:
            parts.append(t.text)
        if t.tail:
            parts.append(t.tail)
    return "".join(parts)


def _para_style(p: ET.Element) -> str | None:
    pPr = p.find(f"{W}pPr")
    if pPr is None:
        return None
    ps = pPr.find(f"{W}pStyle")
    if ps is None:
        return None
    return ps.get(f"{W}val")


def _heading_level(style: str | None) -> int | None:
    if not style:
        return None
    m = re.match(r"heading\s*(\d+)\s*$", style, re.I)
    if m:
        return int(m.group(1))
    if style.lower() in ("title", "subtitle"):
        return 1
    return None


def _escape_rst_line(s: str) -> str:
    """Light escaping for lines that could break inline markup."""
    if re.match(r"^[#=*\-+.\d\\].*", s) and not s.startswith(" "):
        # avoid accidental directive / list starts
        if s.startswith((". ", "- ", "+ ")):
            s = "\\" + s[0] + s[1:]
    return s


def _underline(title: str, char: str) -> str:
    return char * max(len(title), 3)


def _style_to_rst_heading(text: str, level: int) -> str:
    text = text.strip()
    if not text:
        return ""
    text = _escape_rst_line(text)
    if level <= 1:
        u = "="
    elif level == 2:
        u = "-"
    elif level == 3:
        u = "^"
    else:
        u = '"'
    return f"{text}\n{_underline(text, u)}\n\n"


def _is_toc_style(style: str | None) -> bool:
    if not style:
        return False
    return style.upper().startswith("TOC")


def parse_docx(path: str) -> list[tuple[str | None, str, int | None]]:
    """
    Return list of (style_name, text, heading_level_or_None).
    Skips empty paragraphs and TOC lines.
    """
    out: list[tuple[str | None, str, int | None]] = []
    with zipfile.ZipFile(path) as zf:
        root = ET.fromstring(zf.read("word/document.xml"))
    body = root.find(f"{W}body")
    if body is None:
        return out
    for p in body.findall(f"{W}p"):
        style = _para_style(p)
        if _is_toc_style(style):
            continue
        text = _para_text(p).strip()
        if not text:
            continue
        hl = _heading_level(style)
        out.append((style, text, hl))
    return out


def split_segments(blocks: list[tuple[str | None, str, int | None]]):
    """Split on Heading1 into segments: list of (h1_text, blocks_without_h1_line)."""
    segments: list[tuple[str, list[tuple[str | None, str, int | None]]]] = []
    current_title: str | None = None
    current: list[tuple[str | None, str, int | None]] = []

    for style, text, hl in blocks:
        if hl == 1:
            if current_title is not None:
                segments.append((current_title, current))
            current_title = text.strip()
            current = []
        else:
            if current_title is None:
                current_title = "I. Introduction"
            current.append((style, text, hl))
    if current_title is not None:
        segments.append((current_title, current))
    return segments


def merge_duplicate_h1_segments(
    segments: list[tuple[str, list[tuple[str | None, str, int | None]]]],
):
    """Word sometimes repeats the same Heading1; merge bodies in order."""
    out: list[tuple[str, list[tuple[str | None, str, int | None]]]] = []
    for title, body in segments:
        if out and title == out[-1][0]:
            out[-1] = (title, out[-1][1] + body)
        else:
            out.append((title, list(body)))
    return out


def split_ii_and_iii(blocks: list[tuple[str | None, str, int | None]]):
    """After II., Word may run III. in the same segment without a new H1."""
    for i, (st, text, hl) in enumerate(blocks):
        if not re.match(r"^III\.\s*Installation", text):
            continue
        part_a = list(blocks[:i])
        # Title may be glued to body: "III. ... versionThe source"
        m = re.match(
            r"^(III\.\s*Installation[^\n]*?)([A-Z].*)$", text, re.DOTALL
        )
        if m:
            part_b = [(None, m.group(2), None)] + list(blocks[i + 1 :])
        else:
            part_b = list(blocks[i:])
        return part_a, part_b
    return blocks, []


def blocks_to_rst(blocks: list[tuple[str | None, str, int | None]]) -> str:
    parts: list[str] = []
    for style, text, hl in blocks:
        if hl and hl >= 1:
            parts.append(_style_to_rst_heading(text, hl))
        else:
            t = text.replace("\r", "").strip()
            if not t:
                continue
            parts.append(_escape_rst_line(t) + "\n\n")
    return "".join(parts).rstrip() + "\n"


def merge_reference_tail(segments: list[tuple[str, list]]):
    """Move a bogus trailing H1 (numeric citation) into References body."""
    if len(segments) < 2:
        return segments
    last_title, last_body = segments[-1]
    if last_title.lower().startswith("references"):
        return segments
    # e.g. last segment title starts with digit
    if re.match(r"^\d+\.", last_title):
        ref_title, ref_body = segments[-2]
        if ref_title.lower().startswith("references"):
            merged = list(ref_body) + [(None, last_title, None)] + list(last_body)
            return segments[:-2] + [(ref_title, merged)]
    return segments


def build_rst_page(filename: str, body_blocks: list) -> str:
    title = RST_TITLES[filename]
    header = f"{title}\n{_underline(title, '=')}\n\n"
    body = blocks_to_rst(body_blocks)
    return header + body


def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    docx = os.path.join(
        here,
        "..",
        "arxiv",
        "SPARTAN Documentation 20240416",
        "SPARTAN Documentation.docx",
    )
    docx = os.path.normpath(docx)
    if not os.path.isfile(docx):
        raise SystemExit(f"Missing docx: {docx}")

    blocks = parse_docx(docx)
    segments = split_segments(blocks)
    segments = merge_duplicate_h1_segments(segments)
    segments = merge_reference_tail(segments)

    # Map segments to output files (order matches Word H1 sequence, with II split)
    out_map: dict[str, list] = {}

    # 0: I
    t, b = segments[0]
    out_map["introduction.rst"] = b

    # 1: II (+ III inline split)
    t, b = segments[1]
    if not t.startswith("II."):
        raise SystemExit(f"Expected II., got {t!r}")
    a, iii = split_ii_and_iii(b)
    out_map["installation_compiled.rst"] = a
    out_map["installation_source.rst"] = iii
    # 2: IV ...
    rest_keys = RST_FILES_ORDER[3:]  # after intro, compiled, source
    for j, key in enumerate(rest_keys):
        seg_j = 2 + j
        if seg_j >= len(segments):
            raise SystemExit(f"Missing segment for {key}")
        t, b = segments[seg_j]
        out_map[key] = b

    for fname in RST_FILES_ORDER:
        rst = build_rst_page(fname, out_map[fname])
        path = os.path.join(here, fname)
        with open(path, "w", encoding="utf-8") as f:
            f.write(rst)
        print("Wrote", path)


if __name__ == "__main__":
    main()
