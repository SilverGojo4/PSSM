# src/preprocess/smp_parser.py
"""
CDD .smp ASN.1 parser utilities.

This module provides functions to parse CDD .smp files and extract
scaled PSSM matrices.

Design Principles
-----------------
- Keep parsing logic self-contained and reusable
- No pipeline I/O responsibilities
- Return structured data objects for downstream processing
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

MISSING_SENTINEL = -32768

AA_ORDER_20 = list("ARNDCQEGHILKMFPSTWYV")

# From large-scale voting inference on CDD .smp blocks
CDD_ROW_TO_AA: Dict[int, str] = {
    1: "A",
    3: "C",
    4: "D",
    5: "E",
    6: "F",
    7: "G",
    8: "H",
    9: "I",
    10: "K",
    11: "L",
    12: "M",
    13: "N",
    14: "P",
    15: "Q",
    16: "R",
    17: "S",
    18: "T",
    19: "V",
    20: "W",
    22: "Y",
}

# Inverse map for quick lookup
CDD_AA_TO_ROW = {aa: row for row, aa in CDD_ROW_TO_AA.items()}


# ------------------------------ Helpers ------------------------------


def _extract_first(pattern: str, text: str) -> Optional[str]:
    m = re.search(pattern, text, re.DOTALL)
    return m.group(1).strip() if m else None


def _extract_int(pattern: str, text: str) -> Optional[int]:
    s = _extract_first(pattern, text)
    return int(s) if s else None


def _extract_bool(pattern: str, text: str) -> Optional[bool]:
    s = _extract_first(pattern, text)
    if s is None:
        return None
    return s == "TRUE"


def _extract_scores(block_text: str) -> List[int]:
    scores_block = re.search(
        r"finalData\s*\{.*?scores\s*\{\s*(.*?)\s*\}",
        block_text,
        re.DOTALL,
    )
    if not scores_block:
        raise ValueError("Cannot find finalData.scores block in ASN.1 text")

    raw = scores_block.group(1)
    scores = list(map(int, re.findall(r"-?\d+", raw)))

    if not scores:
        raise ValueError("Scores block found but no integers parsed")

    return scores


def _extract_query_sequence(block_text: str) -> Optional[str]:
    query_seq = _extract_first(r'seq-data\s+ncbieaa\s+"([^"]+)"', block_text)
    if not query_seq:
        return None
    return query_seq.replace("\n", "").replace(" ", "")


def split_pssm_blocks(asn1_text: str) -> List[str]:
    parts = re.split(r"PssmWithParameters\s*::=\s*\{", asn1_text)
    return ["PssmWithParameters ::= {" + p for p in parts[1:]]


# ------------------------------ Data Model ------------------------------


@dataclass
class ParsedCDDProfile:
    cd_id: Optional[str]
    title: Optional[str]

    num_rows: int
    num_cols: int
    by_row: bool

    query_seq: Optional[str]
    scaling_factor: Optional[int]

    # 20 x L matrix in AA_ORDER_20 order
    pssm_aa20: np.ndarray


# ------------------------------ Core Parsing ------------------------------


def _reshape_matrix(
    scores: List[int], num_rows: int, num_cols: int, by_row: bool
) -> np.ndarray:
    expected = num_rows * num_cols

    if len(scores) < expected:
        raise ValueError(f"Scores length too short: {len(scores)} < {expected}")

    if len(scores) > expected:
        scores = scores[:expected]

    arr = np.array(scores, dtype=np.int32)

    if by_row:
        return arr.reshape((num_rows, num_cols))

    # byRow == FALSE means column-major layout
    return arr.reshape((num_cols, num_rows)).T


def _convert_to_scaled_float(
    matrix_raw: np.ndarray, scaling_factor: Optional[int]
) -> np.ndarray:
    mat = matrix_raw.astype(float)
    mat[mat == MISSING_SENTINEL] = np.nan

    if scaling_factor and scaling_factor != 0:
        mat = mat / scaling_factor

    return mat


def _extract_cd_id(title: Optional[str]) -> Optional[str]:
    if not title:
        return None
    m = re.match(r"(cd\d+)", title)
    return m.group(1) if m else None


def _extract_pssm_aa20(matrix_scaled: np.ndarray) -> np.ndarray:
    """
    Extract 20 AA rows from full CDD PSSM matrix based on fixed mapping.
    Output shape: (20, num_cols)
    """
    missing = [aa for aa in AA_ORDER_20 if aa not in CDD_AA_TO_ROW]
    if missing:
        raise ValueError(f"CDD mapping incomplete, missing AAs: {missing}")

    rows = []
    for aa in AA_ORDER_20:
        row_idx = CDD_AA_TO_ROW[aa]
        if row_idx >= matrix_scaled.shape[0]:
            raise ValueError(
                f"Row index out of bounds: aa={aa} row={row_idx} "
                f"num_rows={matrix_scaled.shape[0]}"
            )
        rows.append(row_idx)

    return matrix_scaled[rows, :]


def parse_single_pssm_block(block_text: str) -> ParsedCDDProfile:
    title = _extract_first(r'title\s+"([^"]+)"', block_text)

    num_rows = _extract_int(r"numRows\s+(\d+)", block_text)
    num_cols = _extract_int(r"numColumns\s+(\d+)", block_text)
    by_row = _extract_bool(r"byRow\s+(TRUE|FALSE)", block_text)

    if num_rows is None or num_cols is None:
        raise ValueError("Missing numRows/numColumns in block")

    if by_row is None:
        by_row = True

    scaling_factor = _extract_int(r"scalingFactor\s+(\d+)", block_text)
    query_seq = _extract_query_sequence(block_text)

    scores = _extract_scores(block_text)
    matrix_raw = _reshape_matrix(
        scores, num_rows=num_rows, num_cols=num_cols, by_row=by_row
    )
    matrix_scaled = _convert_to_scaled_float(matrix_raw, scaling_factor)

    cd_id = _extract_cd_id(title)
    pssm_aa20 = _extract_pssm_aa20(matrix_scaled)

    return ParsedCDDProfile(
        cd_id=cd_id,
        title=title,
        num_rows=num_rows,
        num_cols=num_cols,
        by_row=by_row,
        query_seq=query_seq,
        scaling_factor=scaling_factor,
        pssm_aa20=pssm_aa20,
    )


def parse_smp_file(path: str | Path) -> List[ParsedCDDProfile]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    text = path.read_text(encoding="utf-8", errors="ignore")
    blocks = split_pssm_blocks(text)

    if not blocks:
        raise ValueError(f"No PssmWithParameters blocks found in {path}")

    return [parse_single_pssm_block(b) for b in blocks]
