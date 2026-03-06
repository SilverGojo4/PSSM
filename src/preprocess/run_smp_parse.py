# src/preprocess/run_smp_parse.py
"""
CDD .smp PSSM Extraction Stage (Branch B)

This stage parses CDD .smp ASN.1 profile files directly and extracts
domain-level PSSM matrices in unified AA column order.

Pipeline Role
-------------
- Input: CD-Search top hit table (query_id -> PSSM_ID + title)
- Input: Local CDD .smp files (from CDD database package)
- Output: domain-level PSSM matrix TSV files (20xL + biochemical features)

Output Directory Layout
-----------------------
results/pssm/matrices/smp/
    ├── {query_id}_{pssm_id}_hseq_with_gap.tsv
    ├── smp_metadata.tsv
    └── smp_parse_error.log
"""

# ============================== Standard Library Imports ==============================
import os
import time
from pathlib import Path
from typing import List, Optional

# ============================== Third-Party Imports ==============================
import numpy as np
import pandas as pd

# ============================== Local Imports ==============================
from src.preprocess.smp_parser import ParsedCDDProfile, parse_smp_file

# ============================== Constants ==============================
AA_ORDER = list("GAILVMFWPCSTYNQHKRDE")

HYDROPHOBIC = {"A", "I", "L", "V", "M", "F", "W", "P", "C"}
CHARGED = {"H", "D", "E", "K", "R"}
POLAR = {"S", "T", "Y", "N", "Q"}


# ============================== Helper Functions ==============================
def _resolve_smp_filename_from_title(title: str) -> str:
    """
    Resolve .smp filename from CD-Search title field.

    Example
    -------
    title:
        "NF041153, organophos_OPH, organophosphate hydrolase OPH..."
    resolved:
        "NF041153.smp"
    """
    if not title or not isinstance(title, str):
        raise ValueError("Empty or invalid title field (cannot resolve .smp filename)")

    base = title.split(",")[0].strip()
    if not base:
        raise ValueError(f"Cannot resolve .smp filename from title: {title}")

    return f"{base}.smp"


def _resolve_smp_path(smp_root_dir: str, title: str) -> Path:
    """
    Resolve .smp file path from title-derived filename.

    Common CDD extraction layouts include:
        <root>/<name>.smp
        <root>/smps/<name>.smp
    """
    smp_fname = _resolve_smp_filename_from_title(title)

    p1 = Path(smp_root_dir) / smp_fname
    p2 = Path(smp_root_dir) / "smps" / smp_fname

    if p1.exists():
        return p1
    if p2.exists():
        return p2

    raise FileNotFoundError(f".smp file not found: {smp_fname} under: {smp_root_dir}")


def _select_profile_block(
    profiles: List[ParsedCDDProfile], pssm_id: int
) -> ParsedCDDProfile:
    """
    Select the correct PSSM block from a .smp file.

    Strategy
    --------
    - If only one block exists: use it
    - Otherwise:
        prefer block whose title contains "cd{pssm_id}"
        else fallback to the first block
    """
    if len(profiles) == 1:
        return profiles[0]

    key = f"cd{pssm_id}"

    for p in profiles:
        if p.title and key in p.title:
            return p

    # fallback (still valid for many .smp files)
    return profiles[0]


def _convert_pssm_to_dataframe(
    pssm_aa20: np.ndarray, query_seq: Optional[str] = None
) -> pd.DataFrame:
    """
    Convert (20 x L) matrix into long-form TSV with AA columns + biochemical features.

    pssm_aa20 is assumed to be in AA_ORDER_20 from parser:
        "ARNDCQEGHILKMFPSTWYV"

    We must reorder into global AA_ORDER:
        "GAILVMFWPCSTYNQHKRDE"
    """
    parser_order = list("ARNDCQEGHILKMFPSTWYV")
    parser_index = {aa: i for i, aa in enumerate(parser_order)}

    missing = [aa for aa in AA_ORDER if aa not in parser_index]
    if missing:
        raise ValueError(f"Missing AA in parser output: {missing}")

    L = pssm_aa20.shape[1]

    # build reordered matrix (L x 20)
    reordered = np.zeros((L, 20), dtype=float)
    for j, aa in enumerate(AA_ORDER):
        reordered[:, j] = pssm_aa20[parser_index[aa], :]

    df = pd.DataFrame(reordered, columns=AA_ORDER)

    # Position and Residue columns
    df.insert(0, "Position", list(range(1, L + 1)))

    if query_seq and len(query_seq) == L:
        df.insert(1, "Residue", list(query_seq))
    else:
        df.insert(1, "Residue", ["-"] * L)

    # biochemical features (positive-only rule consistent with Branch A)
    def _sum_positive(row, aa_set):
        return sum(float(row[aa]) for aa in AA_ORDER if aa in aa_set and row[aa] > 0)

    po_vals = []
    hy_vals = []
    ch_vals = []
    absch_vals = []
    abshypo_vals = []

    for _, r in df.iterrows():
        po = _sum_positive(r, POLAR)
        hy = _sum_positive(r, HYDROPHOBIC)
        ch = _sum_positive(r, CHARGED)

        po_vals.append(po)
        hy_vals.append(hy)
        ch_vals.append(ch)
        absch_vals.append(abs(hy - ch))
        abshypo_vals.append(abs(hy - po))

    df["Po"] = po_vals
    df["Hy"] = hy_vals
    df["Ch"] = ch_vals
    df["|Hy-Ch|"] = absch_vals
    df["|Hy-Po|"] = abshypo_vals

    return df


# ============================== Main Stage Function ==============================
def run_smp_parse(**kwargs) -> None:
    """
    Pipeline Stage: Extract PSSM matrices directly from CDD .smp files (Branch B).
    """
    start_time = time.time()

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║              [ CDD .smp PSSM Extraction Stage Started ]              ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    cdsearch_table = kwargs.get("smp_cdsearch_table")
    smp_root_dir = kwargs.get("smp_root_dir")
    output_dir = kwargs.get("smp_matrix_output_dir")

    if not cdsearch_table:
        raise ValueError("smp_cdsearch_table must be provided")
    if not smp_root_dir:
        raise ValueError("smp_root_dir must be provided")
    if not output_dir:
        raise ValueError("smp_matrix_output_dir must be provided")

    if not os.path.exists(cdsearch_table):
        raise FileNotFoundError(f"CD-search table not found: {cdsearch_table}")

    os.makedirs(output_dir, exist_ok=True)

    meta_path = os.path.join(output_dir, "smp_metadata.tsv")
    log_path = os.path.join(output_dir, "smp_parse_error.log")

    print(f"📄 CD-Search Table      : {cdsearch_table}")
    print(f"📂 SMP Root Dir         : {smp_root_dir}")
    print(f"📁 Output Matrix Dir    : {output_dir}")
    print(f"🧾 Metadata Summary     : {meta_path}")
    print(f"🐛 Error Log            : {log_path}\n")

    df_cd = pd.read_csv(cdsearch_table, sep="\t")

    required_cols = {"query_id", "PSSM_ID", "title"}
    if not required_cols.issubset(df_cd.columns):
        raise ValueError(f"CD-search table missing required columns: {required_cols}")

    total = len(df_cd)
    print(f"🔬 Total domain hits to process: {total}\n")

    summary_rows = []
    success = skipped = failed = 0

    # ===================== Parse Loop =====================
    for idx, row in df_cd.iterrows():
        query_id = row["query_id"]
        pssm_id = int(row["PSSM_ID"])
        title = str(row["title"])

        out_name = f"{query_id}_{pssm_id}_hseq_with_gap.tsv"
        out_path = os.path.join(output_dir, out_name)

        progress = (idx + 1) / total * 100  # type: ignore

        print("─" * 70)
        print(
            f"▶️  [{idx+1:3d}/{total:<3d} | {progress:5.1f}% ] {query_id} (PSSM {pssm_id})"  # type: ignore
        )
        print("─" * 70)

        if os.path.exists(out_path):
            skipped += 1
            summary_rows.append(
                {
                    "query_id": query_id,
                    "PSSM_ID": pssm_id,
                    "status": "skipped",
                    "note": "already exists",
                }
            )
            print(f"   ⏩ Skip (already extracted): {out_name}\n")
            continue

        try:
            smp_path = _resolve_smp_path(smp_root_dir, title)
            smp_fname = smp_path.name

            profiles = parse_smp_file(smp_path)
            selected = _select_profile_block(profiles, pssm_id)

            matrix_df = _convert_pssm_to_dataframe(
                selected.pssm_aa20,
                query_seq=selected.query_seq,
            )

            matrix_df.to_csv(out_path, sep="\t", index=False)

            success += 1
            summary_rows.append(
                {
                    "query_id": query_id,
                    "PSSM_ID": pssm_id,
                    "smp_file": smp_fname,
                    "status": "success",
                    "note": str(smp_path),
                }
            )

            print(f"   ✅ SMP parsed successfully -> {out_name}")
            print(f"   📄 SMP file used: {smp_fname}\n")

        except Exception as e:
            failed += 1
            summary_rows.append(
                {
                    "query_id": query_id,
                    "PSSM_ID": pssm_id,
                    "status": "error",
                    "note": str(e),
                }
            )

            with open(log_path, "a") as log:
                log.write(f"[{query_id} | PSSM {pssm_id}] {e}\n")

            print(f"   ❌ SMP parsing failed: {e}\n")

    # ===================== Save Metadata Summary =====================
    pd.DataFrame(summary_rows).to_csv(meta_path, sep="\t", index=False)

    # ===================== Summary =====================
    elapsed = time.time() - start_time

    print("\n" + "╔" + "═" * 70 + "╗")
    print("║" + " " * 22 + "CDD SMP Extraction Summary" + " " * 21 + "║")
    print("╚" + "═" * 70 + "╝\n")

    print("📊 Summary Statistics")
    print("─" * 70)
    print(f"Total hits     : {total}")
    print(f"Successful     : {success}")
    print(f"Skipped        : {skipped}")
    print(f"Failed         : {failed}\n")

    print("📄 Output Files")
    print("─" * 70)
    print(f"📁 Matrix Dir         : {output_dir}")
    print(f"🧾 Metadata Summary   : {meta_path}")
    print(f"🐛 Error Log          : {log_path}")
    print(f"⏱️  Elapsed time      : {elapsed:.2f} seconds\n")

    print("✅  CDD .smp extraction stage completed successfully!\n")
