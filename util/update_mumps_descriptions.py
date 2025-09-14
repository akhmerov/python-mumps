#!/usr/bin/env python3
r"""
Single-shot updater that injects parameter descriptions for MUMPS arrays
(ICNTL, CNTL, INFO, RINFO, RINFOG) into ``src/mumps/enums.py`` using the
official MUMPS Users' Guide as the primary source of truth.

Main workflow
- Resolve a local copy of the Users' Guide as plain text (preferred), falling
    back to the Markdown copy if present, or download the PDF and convert it to
    text (requires pdfminer.six) when needed.
- Parse the Control and Information parameter sections and extract a snippet per
    parameter, preserving original lines and adding provenance (page number when
    available, source filename, and 1-based line-range) to the generated comments.
- Inject or update these snippets in ``src/mumps/enums.py`` under stable markers
    and update the embedded "MUMPS LICENSE" section if available locally.

Extraction heuristics (robust to TXT/MD formatting)
- Sections are detected by matching body headers (e.g., "6 Control parameters",
    "6.1 Integer control parameters", "7 Information parameters", "7.1", "7.2").
    Scan ranges per array are clamped to those sections to avoid incidental
    mentions elsewhere in the guide.
- Parameter headings are detected line-by-line with a BOL regex allowing leading
    spaces: ``^(CNTL|ICNTL|INFO|RINFO|RINFOG)\(\d+\)``. Lines that mention multiple
    indices (e.g., "RINFOG(7), RINFOG(8) and RINFOG(9)") are expanded so each index
    receives the same snippet.
- For each candidate heading, the snippet spans from that heading line up to the
    next heading line, but never past the end of the array's section (clamped to
    the section boundary).
- Sequential matching heuristic with quality preference:
        - Prefer full description blocks (non-bullet headings) over one-line bullet
            summaries when both exist for the same (ARRAY, index).
        - Among candidates, prefer the one with the longest span up to the next
            heading, as a proxy for descriptive richness.
        - Tie-breaker: earliest position in the document. This still prevents a
            later incidental re-mention (like a cross-reference) from overriding the
            true definition.

Page numbers and cleaning
- The Users' Guide TXT may contain page numbers or form feeds; the preprocessor
    maps each line to its current page number (simple integer-only line) and
    records the page on which a description begins in the snippet header, e.g.:
    ``# === Begin MUMPS snippet: ICNTL(4) page 72 from userguide_5.8.1.txt:4035-4042 ===``
- Within the snippet body, page-number lines and immediately surrounding blank
    lines are removed; leading/trailing blank lines are also trimmed. Snippet text
    is then left-dedented by the common leading spaces.

Marker format (begin/end comments)
- Begin: ``# === Begin MUMPS snippet: <ARRAY>(<index>) [page N] [from <file>:<a>-<b>] ===``
- End:   ``# === End MUMPS snippet ===``
    Both are treated as opaque markers by this script; it will replace existing
    blocks, collapse duplicates, and tolerate optional "page N" before the
    optional "from ..." suffix.

Usage:
        pixi run python util/update_mumps_descriptions.py

Notes:
- License injection is best-effort and skipped if ``mumps_sources/LICENSE`` is
    not available. Version annotation is derived from the environment variable
    ``MUMPS_VERSION``, from ``pixi list`` output, or parsed from the Users' Guide.
- If you rely on PDF-to-text conversion, install pdfminer.six in your
    environment (conda-forge: "pdfminer.six").
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]
ENUMS_PATH = ROOT / "src" / "mumps" / "enums.py"
SRC_DIR = ROOT / "mumps_sources"
LICENSE_PATH = SRC_DIR / "LICENSE"
# Users' Guide local copies (preferred resolution order)
USERGUIDE_TXT_GLOB = "userguide_*.txt"
USERGUIDE_MD_GLOB = "userguide_*.md"
USERGUIDE_PDF_NAME_TMPL = "userguide_{ver}.pdf"

# Note: legacy regex utilities removed as they were unused.

# Canonical lengths used to bound array index expansions from ranges like ICNTL(52-55)
LEN_CNTL = 15
LEN_ICNTL = 60
LEN_INFO = 80
LEN_INFOG = 80
LEN_RINFO = 40
LEN_RINFOG = 40

# Canonical array names and derived helpers
ARRAY_NAMES: tuple[str, ...] = ("ICNTL", "CNTL", "INFO", "INFOG", "RINFO", "RINFOG")
LEN_BY_ARRAY: dict[str, int] = {
    "ICNTL": LEN_ICNTL,
    "CNTL": LEN_CNTL,
    "INFO": LEN_INFO,
    "INFOG": LEN_INFOG,
    "RINFO": LEN_RINFO,
    "RINFOG": LEN_RINFOG,
}
ARRAYS_RX = "|".join(ARRAY_NAMES)


# -------------
# Fetch helpers
# -------------


def _get_version_from_pixi() -> Optional[str]:
    try:
        p = subprocess.run(
            ["pixi", "list", "mumps-seq", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        data = json.loads(p.stdout)
        if isinstance(data, dict):
            for k in ("version", "versions", "pkg_version"):
                if k in data:
                    return data[k]
            if "mumps-seq" in data and isinstance(data["mumps-seq"], dict):
                return data["mumps-seq"].get("version")
        if isinstance(data, list) and data:
            first = data[0]
            if isinstance(first, dict) and "version" in first:
                return first["version"]
    except Exception:
        return None
    return None


def _fetch_mumps_sources_if_needed() -> None:
    """Deprecated: kept for license retrieval only.

    Fetch a MUMPS source archive to populate mumps_sources/ with LICENSE if
    it's missing. Parameter extraction no longer depends on Fortran sources.
    """
    need_license = not LICENSE_PATH.exists()
    if not need_license:
        return

    version = os.environ.get("MUMPS_VERSION") or _get_version_from_pixi()
    if not version:
        print(
            "[update_enums_inline] Could not determine MUMPS version (set MUMPS_VERSION or install pixi). Skipping source fetch."
        )
        return

    url = f"https://mumps-solver.org/MUMPS_{version}.tar.gz"
    SRC_DIR.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as td:
        tmp_tar = Path(td) / f"MUMPS_{version}.tar.gz"
        tmp_extract = Path(td) / f"MUMPS_{version}"
        try:
            print(f"[update_enums_inline] Downloading {url} ...")
            urllib.request.urlretrieve(url, tmp_tar)
        except Exception as e:
            print(f"[update_enums_inline] Download failed: {e}")
            return
        try:
            with tarfile.open(tmp_tar, "r:gz") as tar:
                tar.extractall(path=tmp_extract)
        except Exception as e:
            print(f"[update_enums_inline] Extraction failed: {e}")
            return
        # license-like files
        found_any = False
        for lp in ("LICENSE", "LICENSE.txt", "COPYING", "COPYRIGHT"):
            for p in tmp_extract.glob(f"**/{lp}*"):
                if p.is_file():
                    target = SRC_DIR / (
                        "LICENSE" if lp.startswith("LICENSE") else p.name
                    )
                    try:
                        shutil.copy2(p, target)
                        print(
                            f"[update_enums_inline] Copied license {p.name} -> {target}"
                        )
                        found_any = True
                    except Exception as e:
                        print(f"[update_enums_inline] Copy failed for {p}: {e}")
        if not found_any:
            print(
                "[update_enums_inline] Warning: no license-like files found in archive."
            )


def _prune_sources_dir() -> None:
    """Keep only files required by the project in mumps_sources/."""
    if not SRC_DIR.exists():
        return
    keep_names = {"LICENSE", "LICENSE.txt", "COPYING", "COPYRIGHT"}
    for p in SRC_DIR.iterdir():
        try:
            if p.is_file() and p.name not in keep_names:
                p.unlink()
                print(f"[update_enums_inline] Removed unneeded file: {p}")
            elif p.is_dir():
                shutil.rmtree(p)
                print(f"[update_enums_inline] Removed directory: {p}")
        except Exception as e:
            print(f"[update_enums_inline] Failed to remove {p}: {e}")


# ------------------------
# Users' Guide resolution
# ------------------------


def _parse_version_from_userguide_text(text: str) -> Optional[str]:
    """Attempt to parse a version like '5.8.1' from the Users' Guide header.

    Looks for patterns like 'MUMPS 5.8.1' in the first ~200 lines.
    """
    head = "\n".join(text.splitlines()[:200])
    m = re.search(r"MUMPS\s+([0-9]+\.[0-9]+\.[0-9]+)", head)
    return m.group(1) if m else None


def _find_local_userguide_text() -> Optional[Path]:
    """Find a local Users' Guide text or markdown file in the repository root.

    Preference: *.txt, then *.md; returns the first match in lexicographic order
    (stable when only one version is present).
    """
    # Prefer explicit text files
    txt_candidates = sorted(ROOT.glob(USERGUIDE_TXT_GLOB))
    if txt_candidates:
        return txt_candidates[0]
    md_candidates = sorted(ROOT.glob(USERGUIDE_MD_GLOB))
    if md_candidates:
        return md_candidates[0]
    return None


def _download_userguide_pdf(version: str) -> Optional[Path]:
    """Attempt to download the Users' Guide PDF for a given version.

    Tries a list of candidate URLs; returns local path if successful.
    """
    candidates = [
        f"https://mumps-solver.org/doc/userguide_{version}.pdf",
        f"http://graal.ens-lyon.fr/MUMPS/userguide_{version}.pdf",
        f"https://graal.ens-lyon.fr/MUMPS/userguide_{version}.pdf",
    ]
    SRC_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = SRC_DIR / USERGUIDE_PDF_NAME_TMPL.format(ver=version)
    for url in candidates:
        try:
            print(f"[update_enums_inline] Trying Users' Guide URL: {url}")
            urllib.request.urlretrieve(url, pdf_path)
            if pdf_path.exists() and pdf_path.stat().st_size > 0:
                print(f"[update_enums_inline] Downloaded Users' Guide -> {pdf_path}")
                return pdf_path
        except Exception as e:
            print(f"[update_enums_inline] Failed to download from {url}: {e}")
    return None


def _pdf_to_text(pdf_path: Path) -> Optional[str]:
    """Convert a PDF to text, preferring the external `pdftotext -layout` tool.

    Fallback to pdfminer.six when pdftotext is not available.
    Returns the extracted text or None if conversion fails.
    """
    # Try pdftotext (Poppler)
    try:
        # Use stdout output with '-' and enforce UTF-8 encoding
        res = subprocess.run(
            [
                "pdftotext",
                "-layout",
                "-enc",
                "UTF-8",
                str(pdf_path),
                "-",
            ],
            check=True,
            capture_output=True,
        )
        data = res.stdout
        if isinstance(data, bytes):
            try:
                return data.decode("utf-8", errors="ignore")
            except Exception:
                return data.decode("latin1", errors="ignore")
        return str(data)
    except FileNotFoundError:
        print(
            "[update_enums_inline] pdftotext not found; attempting pdfminer.six as a fallback."
        )
    except subprocess.CalledProcessError as e:
        print(
            f"[update_enums_inline] pdftotext failed (exit {e.returncode}); attempting pdfminer.six."
        )
    except Exception as e:
        print(f"[update_enums_inline] pdftotext error: {e}; attempting pdfminer.six.")

    # Fallback: pdfminer.six
    try:
        from pdfminer.high_level import extract_text  # type: ignore
    except Exception:
        print(
            "[update_enums_inline] pdfminer.six not available; cannot convert PDF to text."
        )
        return None
    try:
        text = extract_text(str(pdf_path))
        return text
    except Exception as e:
        print(f"[update_enums_inline] PDF->text conversion failed: {e}")
        return None


def _resolve_userguide_text() -> tuple[Optional[Path], Optional[str]]:
    """Resolve a usable Users' Guide text source.

    Returns (path, text) such that at least one is non-None.
    Priority:
    1) Local userguide_*.txt
    2) Local userguide_*.md
    3) Download PDF for resolved version and convert to text
    """
    # 1/2: Local files in repo
    p = _find_local_userguide_text()
    if p and p.suffix.lower() in {".txt", ".md"}:
        try:
            return p, read_text(p)
        except Exception as e:
            print(f"[update_enums_inline] Failed to read local guide {p}: {e}")

    # 3: Try to download PDF using known or discovered version
    version = os.environ.get("MUMPS_VERSION") or _get_version_from_pixi()
    # Try to discover version from a nearby markdown if present
    if not version:
        md = next(iter(sorted(ROOT.glob(USERGUIDE_MD_GLOB))), None)
        if md:
            try:
                md_text = read_text(md)
                version = _parse_version_from_userguide_text(md_text)
            except Exception:
                pass
    if not version:
        print(
            "[update_enums_inline] Could not resolve MUMPS version; cannot download Users' Guide PDF."
        )
        return None, None

    pdf = _download_userguide_pdf(version)
    if not pdf:
        return None, None
    txt = _pdf_to_text(pdf)
    if not txt:
        return None, None
    # Write out a cached text copy under project root for future runs
    cached = ROOT / f"userguide_{version}.txt"
    try:
        write_text(cached, txt)
        print(f"[update_enums_inline] Cached Users' Guide text -> {cached}")
    except Exception as e:
        print(f"[update_enums_inline] Failed to cache Users' Guide text: {e}")
    return cached, txt


@dataclass
class ParamChunk:
    array: str
    index: int
    source_file: Optional[str]
    start_line: Optional[int]
    end_line: Optional[int]
    start_page: Optional[int]
    snippet: str

    @property
    def begin_marker(self) -> str:
        meta_parts = []
        if self.start_page:
            meta_parts.append(f"page {self.start_page}")
        if self.source_file and self.start_line and self.end_line:
            meta_parts.append(
                f"from {self.source_file}:{self.start_line}-{self.end_line}"
            )
        meta = (" " + " ".join(meta_parts)) if meta_parts else ""
        return f"# === Begin MUMPS snippet: {self.array}({self.index}){meta} ==="

    @property
    def end_marker(self) -> str:
        return "# === End MUMPS snippet ==="

    def make_comment_block(self) -> List[str]:
        lines = [self.begin_marker]
        for ln in self.snippet.splitlines():
            lines.append(f"# {ln}".rstrip())
        lines.append(self.end_marker)
        return lines


# ------------
# tiny helpers
# ------------


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore")


def write_text(path: Path, data: str) -> None:
    path.write_text(data, encoding="utf-8")


# --------------------
# shared text utilities
# --------------------


def build_page_map(lines: List[str]) -> List[Optional[int]]:
    """Given document lines, return page number per line by tracking simple page-number lines.

    A page line is like '106' on its own (optionally surrounded by spaces). Form feeds are ignored here;
    callers can treat them as page breaks if needed.
    """
    page_at_line: List[Optional[int]] = [None] * len(lines)
    current_page: Optional[int] = None
    page_re = re.compile(r"^\s*(\d{1,4})\s*$")
    for i, ln in enumerate(lines):
        m = page_re.match(ln)
        if m:
            try:
                current_page = int(m.group(1))
            except Exception:
                pass
        page_at_line[i] = current_page
    return page_at_line


def clean_and_dedent_block(slice_lines: List[str]) -> List[str]:
    """Remove standalone page-number/formfeed lines and surrounding blanks, trim edges, left-dedent.

    This mirrors the logic used for parameter snippets and for section 5.9/8 contextual blocks.
    """
    page_re2 = re.compile(r"^\s*(\d{1,4})\s*$")
    ff_re2 = re.compile(r"^\f\s*$")

    def is_blank(s: str) -> bool:
        return s.strip() == ""

    def is_page_line(s: str) -> bool:
        return bool(page_re2.match(s)) or bool(ff_re2.match(s))

    n = len(slice_lines)
    to_drop: set[int] = set()
    i = 0
    while i < n:
        if is_page_line(slice_lines[i]):
            a = i - 1
            while a >= 0 and is_blank(slice_lines[a]):
                to_drop.add(a)
                a -= 1
            to_drop.add(i)
            b = i + 1
            while b < n and is_blank(slice_lines[b]):
                to_drop.add(b)
                b += 1
            i = b
            continue
        i += 1
    cleaned = [slice_lines[k] for k in range(n) if k not in to_drop]
    # Trim leading/trailing blanks
    start = 0
    while start < len(cleaned) and is_blank(cleaned[start]):
        start += 1
    end = len(cleaned)
    while end > start and is_blank(cleaned[end - 1]):
        end -= 1
    block = cleaned[start:end]

    def leading_spaces(s: str) -> int:
        n = 0
        while n < len(s) and s[n] == " ":
            n += 1
        return n

    non_empty = [ln for ln in block if ln.strip()]
    common = min((leading_spaces(ln) for ln in non_empty), default=0)
    if common > 0:
        block = [ln[common:] if len(ln) >= common else ln for ln in block]
    return block


def resolve_userguide_path_text_pagemap() -> (
    tuple[Optional[Path], Optional[str], List[Optional[int]]]
):
    """Resolve userguide text via existing helper and return (path, text, page_map)."""
    ug_path, ug_text = _resolve_userguide_text()
    if not ug_text:
        return ug_path, None, []
    ug_lines = ug_text.splitlines()
    return ug_path, ug_text, build_page_map(ug_lines)


def find_body_section_span(
    lines_all: List[str], start_rx: str, end_rx_list: List[str]
) -> tuple[int, int]:
    """Find body section span starting at header matching start_rx, skipping TOC entries.

    Heuristic identical to the nested _find_section_span used previously.
    """
    start_pat = re.compile(start_rx)
    end_pats = [re.compile(rx) for rx in end_rx_list]
    candidates: list[int] = []
    for i, ln in enumerate(lines_all):
        if start_pat.match(ln):
            candidates.append(i)
    if not candidates:
        return -1, -1

    def looks_like_toc(s: str) -> bool:
        if re.search(r"\.[\s\.]*\d\s*$", s):
            return True
        if re.search(r"\s\d{1,4}\s*$", s) and (". . ." in s or ".." in s):
            return True
        return False

    good = [i for i in candidates if not looks_like_toc(lines_all[i])]
    start = good[0] if good else candidates[-1]
    end = len(lines_all)
    for j in range(start + 1, len(lines_all)):
        ln = lines_all[j]
        for ep in end_pats:
            if ep.match(ln):
                end = j
                break
        if end != len(lines_all):
            break
    return start, end


def make_snippet_comment_block(
    title: str,
    src_name: Optional[str],
    page: Optional[int],
    a: int,
    b: int,
    slice_lines: List[str],
) -> List[str]:
    meta_parts = []
    if page:
        meta_parts.append(f"page {page}")
    if src_name:
        meta_parts.append(f"from {src_name}:{a+1}-{b}")
    meta = (" " + " ".join(meta_parts)) if meta_parts else ""
    begin = f"# === Begin MUMPS snippet: {title}{meta} ==="
    end = "# === End MUMPS snippet ==="
    out = [begin]
    for ln in slice_lines:
        out.append(f"# {ln}".rstrip())
    out.append(end)
    return out


# (Removed several unused comment helpers during cleanup.)


# --------------------
# extract param chunks
# --------------------


def extract_param_chunks() -> tuple[List[ParamChunk], set[tuple[str, int]]]:
    """Extract parameter description snippets from the MUMPS Users' Guide.

    Parsing strategy (robust against formatting in TXT/MD):
    - Scan the whole document for lines that BEGIN with ARRAY(index) where ARRAY in
      {ICNTL, CNTL, INFO, RINFO, RINFOG} (ignoring leading spaces).
    - Handle lines that declare multiple indices for the same array at once, e.g.:
        "RINFOG(7), `RINFOG(8) and `RINFOG(9) ..." by attributing the same snippet
        to indices 7, 8, and 9.
    - Build snippets from a heading line up to the next heading line.
    - Deduplicate by first occurrence of (ARRAY, index).
    - Preserve source line numbers and dedent common leading spaces.
    """
    path, text = _resolve_userguide_text()
    if not text:
        raise FileNotFoundError(
            "Could not resolve a Users' Guide text; place userguide_<ver>.txt/.md in repo root or install pdfminer.six for PDF conversion."
        )
    lines = text.splitlines()

    # Preprocess page numbers: map each line to the current page number if known
    page_at_line: List[Optional[int]] = build_page_map(lines)
    page_re = re.compile(r"^\s*(\d{1,4})\s*$")
    ff_re = re.compile(r"^\f\s*$")  # form feed on its own line

    # Section boundaries (robust to minor spacing variations)
    def _find_first(pattern: str) -> Optional[int]:
        rx = re.compile(pattern)
        for i, ln in enumerate(lines):
            if rx.match(ln):
                return i
        return None

    def _find_from(offset: int, pattern: str) -> Optional[int]:
        rx = re.compile(pattern)
        for i in range(offset, len(lines)):
            if rx.match(lines[i]):
                return i
        return None

    # Prefer body headers (fully anchored to end of line) over TOC entries
    i_ctrl = _find_first(r"^\s*6\s+Control\s+parameters\s*$") or -1
    i_info_body = _find_first(r"^\s*7\s+Information\s+parameters\s*$") or -1

    i_icntl = _find_from(
        max(i_ctrl, 0), r"^\s*6\.1\s+Integer\s+control\s+parameters\s*$"
    )
    i_cntl = _find_from(
        max(i_ctrl, 0), r"^\s*6\.2\s+Real/complex\s+control\s+parameters\s*$"
    )
    i_info = (
        i_info_body
        if i_info_body != -1
        else _find_first(r"^\s*7\s+Information\s+parameters\s*$")
    )
    i_info1 = _find_from(
        max(i_info or 0, 0),
        r"^\s*7\.1\s+Information\s+local\s+to\s+each\s+processor\s*$",
    )
    i_info2 = _find_from(
        max(i_info or 0, 0),
        r"^\s*7\.2\s+Information\s+available\s+on\s+all\s+processors\s*$",
    )
    i_err = _find_first(r"^\s*8\s+Error\s+and\s+warning\s+diagnostics\s*$")
    nlines = len(lines)
    # Note: We anchor to body headers (not TOC) to define robust section bounds for scanning.

    # Define scan ranges per array to avoid incidental mentions elsewhere in the guide
    ranges_by_array: dict[str, List[tuple[int, int]]] = {
        name: [] for name in ARRAY_NAMES
    }
    # ICNTL in 6.1
    if i_icntl is not None:
        end = (
            i_cntl if i_cntl is not None else (i_info if i_info is not None else nlines)
        )
        ranges_by_array["ICNTL"].append((i_icntl, end))
    # CNTL in 6.2
    if i_cntl is not None:
        end = i_info if i_info is not None else nlines
        ranges_by_array["CNTL"].append((i_cntl, end))
    # INFO and RINFO in 7.1 (local info)
    if i_info1 is not None:
        end = i_info2 if i_info2 is not None else nlines
        ranges_by_array["INFO"].append((i_info1, end))
        ranges_by_array["RINFO"].append((i_info1, end))
    # RINFOG in 7.2 (global info)
    if i_info2 is not None:
        end = i_err if i_err is not None else nlines
        ranges_by_array["RINFOG"].append((i_info2, end))
        ranges_by_array["INFOG"].append((i_info2, end))

    # Primary heading detection at BOL, allowing leading spaces
    # Allow ranges like ICNTL(52-55) and single indices, both at BOL for headings and anywhere on the line for expansion
    head_bol = re.compile(
        rf"^(?P<pre>\s*(?:[•\-*]\s*)?)(?P<arr>{ARRAYS_RX})\s*\(\s*(?P<idxtok>\d+(?:\s*-\s*\d+)?)\s*\)"
    )
    token_anywhere = re.compile(
        rf"\b(?P<arr>{ARRAYS_RX})\s*\(\s*(?P<idxtok>\d+(?:\s*-\s*\d+)?)\s*\)"
    )
    reserved_re = re.compile(r"\b(reserved|not\s+used)\b", re.IGNORECASE)

    # Gather candidate headings with position and all indices on that line for the same array
    # Store whether the line is a bullet-summary (prefix contains '•') to de-prioritize it
    raw_heads: List[tuple[int, str, List[int], bool]] = []
    reserved_keys: set[tuple[str, int]] = set()
    # Track accepted heading positions per array (used to compute per-array boundaries)
    positions_by_array: dict[str, List[int]] = {name: [] for name in ARRAY_NAMES}

    def _effective_indent(s: str) -> int:
        k = 0
        while k < len(s) and s[k] == " ":
            k += 1
        if re.match(r"^\s*[\u2022\-*]\s+", s):  # bullet '•' or '-' or '*'
            return k + 1
        return k

    # Helper to expand an idxtok like "7" or "52-55" into concrete indices and clamp within array logical bounds
    def _expand_indices(arr: str, idxtok: str) -> List[int]:
        def clamp(i: int, lo: int, hi: int) -> int:
            return max(lo, min(hi, i))

        lo, hi = 1, LEN_BY_ARRAY.get(arr, 1000)
        idxtok = idxtok.strip()
        if "-" in idxtok:
            a, b = idxtok.split("-", 1)
            try:
                ia = clamp(int(a.strip()), lo, hi)
                ib = clamp(int(b.strip()), lo, hi)
            except Exception:
                return []
            if ia > ib:
                ia, ib = ib, ia
            # Avoid absurd expansions
            if ib - ia > 1000:
                return []
            return list(range(ia, ib + 1))
        try:
            iv = clamp(int(idxtok), lo, hi)
            return [iv]
        except Exception:
            return []

    def _in_any_range(i: int, ranges: List[tuple[int, int]]) -> bool:
        for a, b in ranges:
            if a <= i < b:
                return True
        return False

    # Scan per array to respect section bounds
    for arr_name, arr_ranges in ranges_by_array.items():
        if not arr_ranges:
            # If section bounds were not detected, fall back to whole document
            arr_ranges = [(0, nlines)]
        # Remember last accepted heading for this array to allow consecutive same-indent one-liners
        last_acc_pos: Optional[int] = None
        last_acc_indent: Optional[int] = None
        for i in range(nlines):
            if not _in_any_range(i, arr_ranges):
                continue
            ln = lines[i]
            m = head_bol.match(ln)
            if not m:
                continue
            # Exclude false headings such as "ICNTL(26)=1" that appear within other sections;
            # if an '=' follows immediately after the token, it's a reference, not a heading.
            suffix = ln[m.end() :]
            if re.match(r"^\s*=", suffix):
                continue
            arr = m.group("arr").upper()
            if arr != arr_name:
                continue
            # Collect all indices (including ranges) for same ARR appearing on this line
            idxs: List[int] = []
            for mm in token_anywhere.finditer(ln):
                if mm.group("arr").upper() == arr:
                    idxs.extend(_expand_indices(arr, mm.group("idxtok")))
            if not idxs:
                continue

            # Apply indentation heuristic unless this heading includes index 1 for the array.
            # First entries (e.g., INFO(1), RINFO(1)) often start blocks where indentation may not
            # strictly decrease; in that case we bypass the rule to avoid filtering out true starts.
            cur_eff = _effective_indent(ln)
            accepted_as: Optional[str] = None  # 'true' (indent-decrease) or 'sequence'
            if 1 not in idxs:
                lookback_ok = True
                for j in range(i - 1, max(i - 6, -1), -1):
                    prev = lines[j]
                    if prev.strip() == "":
                        continue
                    if page_re.match(prev) or ff_re.match(prev):
                        continue
                    prev_eff = _effective_indent(prev)
                    if not (cur_eff < prev_eff):
                        lookback_ok = False
                    break  # only compare to the nearest significant previous line
                if not lookback_ok:
                    # Sequence rule (INFO-like lists only): accept if this heading immediately follows the last
                    # accepted heading for the same array, with the same indentation.
                    if arr_name in {"INFO", "RINFO", "RINFOG"} and (
                        last_acc_pos is not None
                        and i == last_acc_pos + 1
                        and last_acc_indent is not None
                        and cur_eff == last_acc_indent
                    ):
                        accepted_as = "sequence"
                    else:
                        continue
                else:
                    accepted_as = "true"
            else:
                accepted_as = "true"  # first entries bypass indentation rule
            # Keep unique, sorted indices for determinism
            idxs = sorted(dict.fromkeys(idxs))
            is_bullet = "•" in (m.group("pre") or "")
            raw_heads.append((i, arr, idxs, is_bullet))
            # Count boundaries: always count true headings; for sequence acceptance, only for INFO-like arrays
            if accepted_as == "true" or (
                accepted_as == "sequence" and arr_name in {"INFO", "RINFO", "RINFOG"}
            ):
                positions_by_array[arr_name].append(i)
            last_acc_pos = i
            last_acc_indent = cur_eff
            # If the heading line declares these as reserved/not used, remember to drop their snippets
            if reserved_re.search(ln):
                for k in idxs:
                    reserved_keys.add((arr, k))
    # We keep positions of all heading-like lines and defer selection until we compute spans.

    if not raw_heads:
        return [], reserved_keys

    # Expand to candidate positions per key (arr, idx)
    candidates_by_key: dict[tuple[str, int], List[tuple[int, bool]]] = {}
    for pos, arr, idxs, is_bullet in raw_heads:
        for idx in idxs:
            candidates_by_key.setdefault((arr, idx), []).append((pos, is_bullet))

    # Prepare accepted headings by array with their indices to compute robust spans
    # We skip same-index re-mentions when advancing to the next boundary.
    accepted_heads_by_arr: dict[str, List[tuple[int, List[int]]]] = {
        k: [] for k in positions_by_array.keys()
    }
    for pos, arr, idxs, _is_bullet in raw_heads:
        accepted_heads_by_arr[arr].append((pos, idxs))
    for arr in accepted_heads_by_arr:
        accepted_heads_by_arr[arr].sort(key=lambda t: t[0])

    def _next_end_pos(arr: str, idx: int, pos: int) -> int:
        # find the next accepted heading for this array that introduces a different index than current
        lst = accepted_heads_by_arr.get(arr, [])
        for p2, idxs2 in lst:
            if p2 <= pos:
                continue
            # If the next heading mentions indices other than the current idx, we stop there
            # (e.g., INFO(13) -> INFO(14)). If it only re-mentions the same idx (e.g., ICNTL(6)=...), skip it.
            if any(j != idx for j in idxs2):
                return p2
        return nlines

    # Select best position per key: prefer non-bullet headings (full blocks) and
    # longer spans over one-liner bullet summaries; tie-breaker: earliest position.
    selected: List[tuple[int, str, int]] = []
    for (arr, idx), cand_list in candidates_by_key.items():

        def span_of(p: int, _arr: str = arr, _idx: int = idx) -> int:
            return _next_end_pos(_arr, _idx, p) - p

        # If any non-bullet candidates exist for this key, ignore bullet ones.
        any_non_bullet = any(not is_b for (_p, is_b) in cand_list)
        filtered = (
            [(p, is_b) for (p, is_b) in cand_list if not is_b]
            if any_non_bullet
            else list(cand_list)
        )

        # Among remaining candidates, prefer longer spans; tie-breaker: earliest position.
        filtered.sort(key=lambda t: (-span_of(t[0]), t[0]))

        # Additional guard: if the top candidate has a tiny span (likely a summary) and a longer one
        # exists, pick the longer one.
        if filtered:
            best_pos = filtered[0][0]
            selected.append((best_pos, arr, idx))

    # Sort selected by position to compute snippet boundaries
    selected.sort(key=lambda t: t[0])

    # Helpers to clean snippet content: remove page-number lines and surrounding blank lines
    def _clean_snippet(slice_lines: List[str]) -> List[str]:
        return clean_and_dedent_block(slice_lines)

    # Build chunks
    chunks: List[ParamChunk] = []
    for pos, arr, idx in selected:
        end_pos = _next_end_pos(arr, idx, pos)
        # Clamp to end of the array's section range containing this position
        arr_ranges = ranges_by_array.get(arr, []) or [(0, nlines)]
        for a, b in arr_ranges:
            if a <= pos < b:
                if end_pos > b:
                    end_pos = b
                break
        slice_lines = lines[pos:end_pos]
        # Remove page numbers and surrounding blanks from snippet body
        slice_lines = _clean_snippet(slice_lines)

        # Compute provenance (1-based)
        start_line_no = pos + 1
        end_line_idx = end_pos - 1
        while end_line_idx >= pos and lines[end_line_idx].strip() == "":
            end_line_idx -= 1
        end_line_no = (end_line_idx + 1) if end_line_idx >= pos else (pos + 1)

        # Dedent common leading spaces
        def leading_spaces(s: str) -> int:
            n = 0
            while n < len(s) and s[n] == " ":
                n += 1
            return n

        non_empty = [ln for ln in slice_lines if ln.strip()]
        common = min((leading_spaces(ln) for ln in non_empty), default=0)
        snippet_lines = [
            ln[common:] if common > 0 and len(ln) >= common else ln
            for ln in slice_lines
        ]
        snippet = "\n".join(snippet_lines).rstrip() + "\n"

        # Skip reserved/not-used entries entirely
        if (arr, idx) in reserved_keys:
            continue

        chunks.append(
            ParamChunk(
                array=arr,
                index=idx,
                source_file=path.name if path else "userguide.txt",
                start_line=start_line_no,
                end_line=end_line_no,
                start_page=page_at_line[pos] if 0 <= pos < len(page_at_line) else None,
                snippet=snippet,
            )
        )
    return chunks, reserved_keys


# ------------------
# enums.py mutation
# ------------------


def inject_license(enums_src: str, license_text: str) -> str:
    begin = "# MUMPS LICENSE"
    end = "# END MUMPS LICENSE"
    if begin not in enums_src or end not in enums_src:
        return enums_src  # be conservative
    pre, rest = enums_src.split(begin, 1)
    _, post = rest.split(end, 1)
    commented = "\n".join(
        [begin, *[f"# {ln}".rstrip() for ln in license_text.splitlines()], "", end]
    )
    return pre + commented + post


def inject_mumps_version(enums_src: str, version: Optional[str]) -> str:
    """Replace the placeholder line containing 'MUMPS VERSION' with the
    resolved version, preserving indentation. If version is None/empty,
    return the source unchanged.

    Expected header appears near the top:
        # The script was last run using
        # MUMPS VERSION
    We'll replace the second line with, e.g.:
        # MUMPS VERSION: 5.6.2
    """
    if not version:
        return enums_src
    lines = enums_src.splitlines(True)
    ver_re = re.compile(r"^(?P<indent>\s*)#?\s*MUMPS\s+VERSION\b.*$")
    replaced = False
    for i, line in enumerate(lines):
        m = ver_re.match(line)
        if m:
            indent = m.group("indent")
            lines[i] = f"{indent}# MUMPS VERSION: {version}\n"
            replaced = True
            break
    if not replaced:
        # Try to locate the header line and insert after it
        header_re = re.compile(r"^(?P<indent>\s*)#\s*The script was last run using\s*$")
        for i, line in enumerate(lines):
            m = header_re.match(line)
            if m:
                indent = m.group("indent")
                lines.insert(i + 1, f"{indent}# MUMPS VERSION: {version}\n")
                replaced = True
                break
    return "".join(lines)


def find_section_bounds(
    lines: List[str], section_title: str
) -> Optional[tuple[int, int]]:
    # locate '# <TITLE> members' and use the 3-line banner convention
    hdr_idx = None
    for i, line in enumerate(lines):
        if line.strip() == f"# {section_title} members":
            hdr_idx = i
            break
    if hdr_idx is None:
        return None
    start = hdr_idx + 3
    end = len(lines)
    for j in range(start, len(lines)):
        t = lines[j].strip()
        if (
            t
            in {
                "# CNTL members",
                "# ICNTL members",
                "# INFO members",
                "# RINFO members",
                "# RINFOG members",
            }
            and j != hdr_idx + 1
        ):
            end = j - 1
            break
        if t == "# Array wrappers with basic validation":
            end = j - 1
            break
    return (start, max(start, end))


def find_array_class_bounds(lines: List[str], array: str) -> Optional[tuple[int, int]]:
    """Find the class <ARRAY>(ParamArray|RawArray) region: returns (start, end_exclusive)."""
    class_re = re.compile(
        rf"^\s*class\s+{re.escape(array)}\((?:ParamArray|RawArray)\):\s*$"
    )
    start = None
    for i, ln in enumerate(lines):
        if class_re.match(ln):
            start = i
            break
    if start is None:
        return None
    # End at next top-level class or end of file or next section banner
    for j in range(start + 1, len(lines)):
        if lines[j].startswith("class ") and not lines[j].startswith("    "):
            return (start, j)
        t = lines[j].strip()
        if t in {
            "# CNTL members",
            "# ICNTL members",
            "# INFO members",
            "# RINFO members",
            "# RINFOG members",
            "# INFOG members",
        }:
            return (start, j)
    return (start, len(lines))


def find_param_decorator_line(
    lines: List[str], array: str, index: int
) -> Optional[int]:
    bounds = find_array_class_bounds(lines, array)
    if not bounds:
        return None
    start, end = bounds
    # Match decorator possibly with other keyword args, we will later update/add page kw
    deco_re = re.compile(rf"^\s*@param\(.*index\s*=\s*{index}[^)]*\)\s*$")
    for i in range(start + 1, end):
        if deco_re.match(lines[i]):
            return i
    return None


def remove_existing_snippet_block_above(
    lines: List[str], line_idx: int, begin_marker: str, end_marker: str
) -> int:
    i = line_idx - 1
    begin_idx = None
    end_idx = None
    while i >= 0 and line_idx - i <= 300:
        if lines[i].strip() == begin_marker:
            begin_idx = i
            break
        if not (lines[i].strip().startswith("#") or lines[i].strip() == ""):
            break
        i -= 1
    if begin_idx is not None:
        for j in range(begin_idx, line_idx):
            if lines[j].strip() == end_marker:
                end_idx = j
                break
        if end_idx is None:
            return line_idx
        del lines[begin_idx : end_idx + 1]
        line_idx -= end_idx + 1 - begin_idx
        if begin_idx < len(lines) and lines[begin_idx].strip() == "":
            del lines[begin_idx]
            line_idx -= 1
        return begin_idx
    return line_idx


def find_all_existing_snippet_blocks(
    lines: List[str], array: str, index: int, end_marker: str
) -> List[Tuple[int, int]]:
    """Return all existing snippet blocks for ARRAY(index) as (begin,end) pairs."""
    pat = re.compile(
        rf"^#\s===\sBegin\sMUMPS\ssnippet:\s{re.escape(array)}\({index}\)(?:\s+page\s+\d+)?(?:\s+from.*)?\s===\s*$"
    )
    blocks: List[Tuple[int, int]] = []
    i = 0
    while i < len(lines):
        if pat.match(lines[i].strip()):
            b = i
            j = i + 1
            while j < len(lines) and lines[j].strip() != end_marker:
                j += 1
            if j < len(lines):
                blocks.append((b, j))
                i = j + 1
                continue
        i += 1
    return blocks


def insert_block(lines: List[str], insert_at: int, comment_block: List[str]) -> None:
    if insert_at > 0 and lines[insert_at - 1].strip() != "":
        lines.insert(insert_at, "\n")
        insert_at += 1
    for ln in comment_block:
        lines.insert(insert_at, ln + "\n")
        insert_at += 1
    if insert_at < len(lines) and lines[insert_at].strip() != "":
        lines.insert(insert_at, "\n")


def _update_or_add_page_kw_to_decorator_line(src: str, page: Optional[int]) -> str:
    """Given a line like '@param(index=4)' (possibly with spaces or other kwargs),
    return a version that includes/updates 'page=<page>' if page is provided.
    Conservatively only modifies the content inside the parentheses.
    """
    if page is None:
        return src
    m = re.match(r"^(\s*@param\()(.*?)(\)\s*)$", src)
    if not m:
        return src
    pre, args, post = m.group(1), m.group(2), m.group(3)
    # Split by commas, keep simple key=val pairs
    parts = [p.strip() for p in args.split(",") if p.strip()]
    has_page = False
    for i, p in enumerate(parts):
        if re.match(r"^page\s*=", p):
            parts[i] = f"page={int(page)}"
            has_page = True
            break
    if not has_page:
        parts.append(f"page={int(page)}")
    new_args = ", ".join(parts)
    return (
        f"{pre}{new_args}{post}\n"
        if not src.endswith("\n")
        else f"{pre}{new_args}{post}"
    )


def update_enums_inline() -> None:
    # Best-effort: fetch a LICENSE if missing (parameter extraction no longer depends on sources)
    _fetch_mumps_sources_if_needed()
    _prune_sources_dir()

    if not ENUMS_PATH.exists():
        raise FileNotFoundError(f"Missing enums file: {ENUMS_PATH}")
    if not LICENSE_PATH.exists():
        print(
            f"[update_enums_inline] MUMPS LICENSE not found at {LICENSE_PATH}; skipping license injection."
        )

    chunks, reserved_keys = extract_param_chunks()
    enums_src = read_text(ENUMS_PATH)
    license_text = read_text(LICENSE_PATH) if LICENSE_PATH.exists() else None

    # Resolve version (env takes precedence over pixi discovery)
    # Try to resolve version from env/pixi or parse from guide text
    resolved_version = os.environ.get("MUMPS_VERSION") or _get_version_from_pixi()
    if not resolved_version:
        # Parse from any local Users' Guide
        ug_path, ug_text = _resolve_userguide_text()
        if ug_text:
            resolved_version = _parse_version_from_userguide_text(ug_text)

    # license injection
    updated_src = inject_license(enums_src, license_text) if license_text else enums_src
    # version injection (best-effort; no-op if version unknown)
    updated_src = inject_mumps_version(updated_src, resolved_version)
    lines = updated_src.splitlines(True)

    # ----------------------------------------------------------
    # Insert full-text contextual snippets for Subsection 5.9 and Section 8
    # ----------------------------------------------------------
    # We extract the entire body of subsection 5.9 (Error analysis) and the
    # entire Section 8 (Error and warning diagnostics) from the Users' Guide
    # and inject them as commented blocks into enums.py for convenient in-source
    # reference. We place 5.9 just before the RINFOG members banner, and 8 just
    # before the INFO members banner. If such blocks already exist (detected by
    # their unique begin markers), we replace their content.

    # Reuse shared helpers defined above

    def _insert_or_replace_block_after_banner(
        file_lines: list[str],
        header_title: str,
        begin_marker_prefix: str,
        block_lines: list[str],
    ) -> None:
        """Insert or replace a block right AFTER the full 'XXX members' banner.

        Banner is considered to be the header line '# XXX members' and its underline
        line '# -----...'. We insert after both lines (and any immediate blank line).
        If an existing block with the same begin marker exists anywhere, we replace it
        in place instead of inserting anew.
        """
        # If a block with a begin line starting with begin_marker_prefix exists, replace it.
        begin_idx = None
        end_idx = None
        for i, ln in enumerate(file_lines):
            if ln.strip().startswith(begin_marker_prefix):
                begin_idx = i
                # find matching end
                for j in range(i + 1, len(file_lines)):
                    if file_lines[j].strip() == "# === End MUMPS snippet ===":
                        end_idx = j
                        break
                break
        if begin_idx is not None and end_idx is not None:
            # Remove existing block so we can reinsert at canonical position
            del file_lines[begin_idx : end_idx + 1]
            # Remove a single trailing blank line if present at that position
            if begin_idx < len(file_lines) and file_lines[begin_idx].strip() == "":
                del file_lines[begin_idx]
            # Continue to compute insert_at for reinsertion
        # Otherwise, find banner line '# <TITLE> members' and insert AFTER banner
        banner = f"# {header_title} members"
        hdr_line = None
        for i, ln in enumerate(file_lines):
            if ln.strip() == banner:
                hdr_line = i
                break
        if hdr_line is not None:
            # Move to the underline line if present, then beyond it
            insert_at = hdr_line + 1
            if insert_at < len(file_lines) and re.match(
                r"^\s*#\s*-{3,}\s*$", file_lines[insert_at]
            ):
                insert_at += 1
            # Skip a single blank line if present to keep formatting tidy
            if insert_at < len(file_lines) and file_lines[insert_at].strip() == "":
                insert_at += 1
        else:
            # Fallback 1: insert before the class <HEADER_TITLE>(
            class_re = re.compile(rf"^\s*class\s+{re.escape(header_title)}\(")
            insert_at = None  # type: ignore
            for i, ln in enumerate(file_lines):
                if class_re.match(ln):
                    insert_at = i
                    break
            if insert_at is None:
                # Fallback 2: if there are already parameter snippets for this header (e.g. '# === Begin MUMPS snippet: RINFOG('),
                # insert just before the first such snippet so our contextual block sits in front of the section.
                snippet_hdr_re = re.compile(
                    rf"^\s*#\s===\sBegin\sMUMPS\ssnippet:\s{re.escape(header_title)}\(\d+\).*$"
                )
                for i, ln in enumerate(file_lines):
                    if snippet_hdr_re.match(ln):
                        insert_at = i
                        break
            if insert_at is None:
                # Fallback 3: insert at end of file
                insert_at = len(file_lines)
        for bl in block_lines:
            file_lines.insert(insert_at, bl + "\n")
            insert_at += 1
        if insert_at < len(file_lines) and file_lines[insert_at].strip() != "":
            file_lines.insert(insert_at, "\n")

    # Perform extraction and insertion
    ug_path, ug_text, page_map = resolve_userguide_path_text_pagemap()
    if ug_text:
        ug_lines = ug_text.splitlines()
        # Subsection 5.9 span: start at '5.9 ' and end at next '5.<digit+>' or start of section 6
        s59, e59 = find_body_section_span(
            ug_lines,
            r"^\s*5\.9\s+.*$",
            [r"^\s*5\.[1-9]\d*\s+.*$", r"^\s*6\s+.*$"],
        )
        if s59 != -1 and e59 != -1 and e59 > s59:
            s59_lines = clean_and_dedent_block(ug_lines[s59:e59])
            p59 = page_map[s59] if 0 <= s59 < len(page_map) else None
            block59 = make_snippet_comment_block(
                "SUBSECTION(5.9)",
                ug_path.name if ug_path else None,
                p59,
                s59,
                e59,
                s59_lines,
            )
            _insert_or_replace_block_after_banner(
                lines, "RINFOG", "# === Begin MUMPS snippet: SUBSECTION(5.9)", block59
            )

        # Section 8 span: start at '8 Error and warning diagnostics' to the next top-level section
        # Find Section 8 header (skip TOC via earlier heuristic) then scan for next top-level section header
        s8, _tmp = find_body_section_span(
            ug_lines,
            r"^\s*8\s+Error\s+and\s+warning\s+diagnostics\s*$",
            [r"$^"],  # no-op end; we'll compute end ourselves
        )
        # Compute end as next line matching a top-level section header: optional FF, small leading indent (flush-left header),
        # digits, at least one space, then title text starting with a letter; number > 8. This avoids stopping on indented
        # numeric table rows inside Section 8 (e.g., "18  RHS loc").
        e8 = -1
        if s8 != -1:
            # Prefer an exact match for the next section header to avoid accidental matches inside tables.
            exact_next = re.compile(r"^\s*\f?\s*9\s+Calling\s+MUMPS\s+from\s+C\s*$")
            for k in range(s8 + 1, len(ug_lines)):
                if exact_next.match(ug_lines[k]):
                    e8 = k
                    break
            # Fallback: a generic top-level section header
            if e8 == -1:
                top_hdr_fallback = re.compile(r"^\s{0,3}\f?\s*(\d{1,2})\s+[A-Za-z].+$")
                for k in range(s8 + 1, len(ug_lines)):
                    m2 = top_hdr_fallback.match(ug_lines[k])
                    if m2:
                        try:
                            num = int(m2.group(1))
                            if num > 8:
                                e8 = k
                                break
                        except Exception:
                            pass
            if e8 == -1:
                e8 = len(ug_lines)
        if s8 != -1 and e8 != -1 and e8 > s8:
            s8_lines = clean_and_dedent_block(ug_lines[s8:e8])
            p8 = page_map[s8] if 0 <= s8 < len(page_map) else None
            block8 = make_snippet_comment_block(
                "SECTION(8)",
                ug_path.name if ug_path else None,
                p8,
                s8,
                e8,
                s8_lines,
            )
            _insert_or_replace_block_after_banner(
                lines, "INFO", "# === Begin MUMPS snippet: SECTION(8)", block8
            )

    # Pre-pass: sanitize any malformed or duplicate snippet markers left by prior runs
    def _sanitize_snippet_markers(lines: List[str]) -> List[str]:
        """Remove malformed snippet markers and collapse duplicates by key.

        Streaming algorithm:
        - Outside of any block: copy lines as-is except stray End markers.
        - On Begin: start accumulating a block; if another Begin appears before End, drop the previous open block and start anew.
        - On End: if there is an open block, close it and emit it only if this (array, index) hasn't been seen yet; otherwise discard.
        - If EOF with an open block: discard it.
        """
        # Accept any standard Begin marker and capture its full title (e.g., 'ICNTL(4)', 'SECTION(8)', 'SUBSECTION(5.9)').
        # We deduplicate by this title to keep a single block per unique snippet.
        begin_pat = re.compile(
            r"^#\s===\sBegin\sMUMPS\ssnippet:\s(?P<title>.+?)\s===\s*$"
        )
        end_pat = re.compile(r"^#\s===\sEnd\sMUMPS\ssnippet\s===\s*$")

        new_lines: List[str] = []
        seen: set[str] = set()
        acc_block: Optional[List[str]] = None
        acc_key: Optional[str] = None

        for line in lines:
            m = begin_pat.match(line.strip())
            if m:
                # start a new block; drop any previously open one
                acc_block = [line]
                acc_key = m.group("title").strip()
                continue

            if end_pat.match(line.strip()):
                if acc_block is None or acc_key is None:
                    # stray end marker: skip
                    continue
                # close current block
                acc_block.append(line)
                if acc_key not in seen:
                    new_lines.extend(acc_block)
                    seen.add(acc_key)
                # reset
                acc_block = None
                acc_key = None
                continue

            # regular line
            if acc_block is not None:
                acc_block.append(line)
            else:
                new_lines.append(line)

        # drop any unterminated block at EOF
        return new_lines

    lines = _sanitize_snippet_markers(lines)

    # Remove any existing blocks for reserved/not-used indices
    END_MARK = "# === End MUMPS snippet ==="
    for arr, idx in sorted(reserved_keys):
        blocks = find_all_existing_snippet_blocks(lines, arr, idx, END_MARK)
        # remove from the bottom up to preserve indices
        for b, e in reversed(blocks):
            del lines[b : e + 1]
            # also remove a single extra blank line if present at the position
            if b < len(lines) and lines[b].strip() == "":
                del lines[b]

    for ch in chunks:
        array = ch.array
        index = ch.index
        comment_block = ch.make_comment_block()
        # 1) Find and replace existing blocks; if multiple exist, collapse to one
        blocks = find_all_existing_snippet_blocks(lines, array, index, ch.end_marker)
        if blocks:
            # Replace first block, remove the rest
            b0, e0 = blocks[0]
            begin_line = lines[b0]
            indent = begin_line[: len(begin_line) - len(begin_line.lstrip())]
            replacement = [indent + ln + "\n" for ln in comment_block]
            lines[b0 : e0 + 1] = replacement
            # Offsets have changed; remove remaining blocks by searching again
            # until only one remains
            while True:
                extra = find_all_existing_snippet_blocks(
                    lines, array, index, ch.end_marker
                )
                if len(extra) <= 1:
                    break
                # remove the last one to minimize shifting
                rb, re_ = extra[-1]
                del lines[rb : re_ + 1]
            # After snippet replacement, try to update/add page kw to existing decorator
            if array in ARRAY_NAMES:
                deco_line = find_param_decorator_line(lines, array, index)
                if deco_line is not None:
                    # Ensure a blank line separates snippet end and decorator
                    end_after = b0 + len(replacement)
                    if end_after < len(lines) and lines[end_after].strip() != "":
                        lines.insert(end_after, "\n")
                    lines[deco_line] = _update_or_add_page_kw_to_decorator_line(
                        lines[deco_line], ch.start_page
                    )
                else:
                    # Insert a minimal placeholder after the snippet within the class bounds
                    bounds = find_array_class_bounds(lines, array)
                    if bounds is not None:
                        # Find insertion point: just after the snippet end (e0 was replaced; find location of this snippet again)
                        # We look for the last line of the replacement we just inserted
                        insert_at = b0 + len(replacement)
                        # Ensure insertion is inside class; if not, append near class end
                        cstart, cend = bounds
                        insert_at = max(insert_at, cstart + 1)
                        insert_at = min(insert_at, cend)
                        # Determine indent from class body (use 4 spaces)
                        indent = "    "
                        # Build placeholder
                        page_kw = (
                            f", page={int(ch.start_page)}" if ch.start_page else ""
                        )
                        placeholder = [
                            "\n"
                            if (
                                insert_at < len(lines)
                                and lines[insert_at].strip() != ""
                            )
                            else "",
                            f"{indent}@param(index={index}{page_kw})\n",
                            f"{indent}class param_{index}:\n",
                            f'{indent}    "Placeholder parameter, not processed yet."\n',
                            "\n",
                        ]
                        for ln in reversed([x for x in placeholder if x != ""]):
                            lines.insert(insert_at, ln)
            continue

        # 2) Prefer inserting above @param(index=K) in ICNTL/CNTL class
        if array in ARRAY_NAMES:
            deco_line = find_param_decorator_line(lines, array, index)
            if deco_line is not None:
                insert_at = remove_existing_snippet_block_above(
                    lines, deco_line, ch.begin_marker, ch.end_marker
                )
                insert_block(lines, insert_at, comment_block)
                # Ensure a blank line between inserted snippet and the decorator
                if deco_line < len(lines) and deco_line > 0:
                    # After insertion, the decorator index may have shifted; find it again
                    deco_line = (
                        find_param_decorator_line(lines, array, index) or deco_line
                    )
                    # Insert a blank line if the line before decorator is not blank and not an end marker
                    if deco_line > 0 and lines[deco_line - 1].strip() != "":
                        lines.insert(deco_line, "\n")
                # update page kw on decorator line
                lines[deco_line] = _update_or_add_page_kw_to_decorator_line(
                    lines[deco_line], ch.start_page
                )
                continue
            else:
                # Missing decorator: insert snippet at a reasonable spot in class and add placeholder
                bounds = find_array_class_bounds(lines, array)
                if bounds is not None:
                    cstart, cend = bounds
                    # Insert snippet near end of class but before next section
                    insert_at = cend
                    # Create the snippet block first
                    insert_block(lines, insert_at, comment_block)
                    # After insertion, compute new class end offset
                    # Simplify by re-finding bounds
                    new_bounds = find_array_class_bounds(lines, array) or bounds
                    cstart, cend = new_bounds
                    # Insert minimal placeholder directly after the snippet we just inserted
                    # Find position of end marker we just added
                    # Scan backwards from cend to find the end marker line
                    end_marker_line = None
                    for i in range(cend - 1, cstart, -1):
                        if lines[i].strip() == ch.end_marker:
                            end_marker_line = i
                            break
                    insert_at2 = (end_marker_line + 1) if end_marker_line else cend
                    indent = "    "
                    page_kw = f", page={int(ch.start_page)}" if ch.start_page else ""
                    placeholder = [
                        "\n"
                        if (insert_at2 < len(lines) and lines[insert_at2].strip() != "")
                        else "",
                        f"{indent}@param(index={index}{page_kw})\n",
                        f"{indent}class param_{index}:\n",
                        f'{indent}    "Minimal placeholder; update tokens and docs as needed."\n',
                        "\n",
                    ]
                    for ln in reversed([x for x in placeholder if x != ""]):
                        lines.insert(insert_at2, ln)
                    continue

        # 3) Fallback: append at end of the corresponding '<ARRAY> members' section
        bounds = find_section_bounds(lines, array)
        if bounds is not None:
            _, end = bounds
            insert_block(lines, end, comment_block)
            continue
        # 4) Last resort: append at file end
        if lines and lines[-1] and not lines[-1].endswith("\n"):
            lines[-1] = lines[-1] + "\n"
        insert_block(lines, len(lines), comment_block)

    write_text(ENUMS_PATH, "".join(lines))


if __name__ == "__main__":
    update_enums_inline()
