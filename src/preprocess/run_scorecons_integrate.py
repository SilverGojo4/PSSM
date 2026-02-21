# src/preprocess/run_scorecons_integrate.py
"""
Scorecons Conservation Integration Stage

Project alignment-based Scorecons conservation scores
back to full-length protein residue positions.

Updated Design (Branch-aware Integrated Output)
-----------------------------------------------
This stage no longer writes to a separate conservation/reconstruct folder.

Instead, it integrates Scorecons scores directly into reconstructed
full-length PSSM tables and writes the updated integrated tables into:

    results/pssm/integrated/<branch>/

Inputs
------
- CD-search alignment table: cdsearch_top_hits_detailed.tsv
- Reconstructed full-length PSSM tables: results/pssm/reconstruct/<branch>/*.tsv
- Scorecons server outputs: results/scorecons/*.txt

Outputs
-------
- Integrated PSSM + Scorecons tables: results/pssm/integrated/<branch>/*.tsv
- Error log: <integrated_output_dir>/scorecons_integrate_error.log

Important
---------
This stage must be executed BEFORE ConSurf integration stage.
"""

# ============================== Standard Library Imports ==============================
import os
import time

# ============================== Third-Party Imports ==============================
import numpy as np
import pandas as pd
from Bio import SeqIO

# ============================== Constants ==============================
AA_ORDER = list("GAILVMFWPCSTYNQHKRDE")
PSSM_NUMERIC_COLS = AA_ORDER + ["Po", "Hy", "Ch", "Hy+Ch-Po", "|Hy-Ch|"]


# ============================== Helper Functions ==============================
def _parse_scorecons_txt(path: str) -> pd.DataFrame:
    """
    Parse Scorecons server output file.

    Expected format
    ---------------
    - Header lines (first 9 lines) are ignored
    - Data lines contain:
        <score> # <two residues>

    Returns
    -------
    DataFrame with columns:
        aln_index, conservation, res_1, res_2
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Scorecons file not found: {path}")

    records = []

    with open(path) as f:
        lines = f.readlines()

    data_lines = lines[9:]  # skip header

    aln_index = 1
    for line in data_lines:
        line = line.strip()
        if not line:
            continue

        try:
            left, right = line.split("#", 1)
            score = float(left.strip())
            residues = right.strip()

            records.append(
                {
                    "aln_index": aln_index,
                    "conservation": score,
                    "res_1": residues[0] if len(residues) > 0 else None,
                    "res_2": residues[1] if len(residues) > 1 else None,
                }
            )
            aln_index += 1

        except Exception:
            continue

    if not records:
        raise ValueError(f"No valid Scorecons scores found: {path}")

    return pd.DataFrame(records)


def _load_alignment_table(path: str) -> pd.DataFrame:
    """
    Load cdsearch_top_hits_detailed.tsv.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"CD-search table not found: {path}")
    return pd.read_csv(path, sep="\t")


def _load_original_sequences(fasta_path: str) -> dict:
    """
    Load original sequences from FASTA.

    Returns
    -------
    dict
        {seq_id: sequence_string}
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Input FASTA not found: {fasta_path}")
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}


def _restore_numeric_types(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """
    Restore numeric dtypes:
    - If a column contains only integer-like values (or NA), cast to Int64
    - Otherwise keep float

    Notes
    -----
    Pandas will often upcast integer columns to float when NA exists.
    This helper restores stable dtypes to avoid integer -> float drift.
    """
    for c in cols:
        if c not in df.columns:
            continue

        series = pd.to_numeric(df[c], errors="coerce")

        if series.isna().all():
            continue

        non_na = series.dropna()

        # safer integer-like check (avoid float precision issues)
        if np.isclose(non_na, non_na.round()).all():
            df[c] = series.round().astype("Int64")
        else:
            df[c] = series.astype(float)

    return df


def _load_pssm_reconstruct_table(
    pssm_reconstruct_dir: str, query_id: str
) -> pd.DataFrame:
    """
    Load reconstructed full-length PSSM table.

    Expected path:
        <pssm_reconstruct_dir>/<query_id>.tsv
    """
    fpath = os.path.join(pssm_reconstruct_dir, f"{query_id}.tsv")
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"PSSM reconstruct file not found: {fpath}")

    df = pd.read_csv(fpath, sep="\t")
    df = _restore_numeric_types(df, PSSM_NUMERIC_COLS)

    return df


def _resolve_scorecons_path(scorecons_dir: str, query_id: str) -> str:
    """
    Resolve Scorecons txt file path.

    Expected naming:
        <scorecons_dir>/<query_id>.txt
    """
    fpath = os.path.join(scorecons_dir, f"{query_id}.txt")
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"Scorecons file not found: {fpath}")
    return fpath


def _integrate_scorecons(
    aln_row: pd.Series,
    scorecons_df: pd.DataFrame,
    base_table: pd.DataFrame,
) -> pd.DataFrame:
    """
    Integrate alignment-based Scorecons scores into full-length table.

    Rules
    -----
    - Query gap positions are not projectable
    - Only query residues receive scores
    - Unaligned positions remain NaN

    Notes
    -----
    We only use qseq + qstart mapping (query coordinate system),
    and project Scorecons aln_index back to full-length residue positions.
    """
    scorecons_df = scorecons_df.set_index("aln_index")

    qseq = str(aln_row["qseq"])
    orig_pos = int(aln_row["qstart"])  # 1-based
    aln_index = 1

    out = base_table.copy()

    title = str(aln_row["title"]).split(",")[0].strip()
    col_name = f"Conservation ({title})"

    if col_name in out.columns:
        out = out.drop(columns=[col_name])

    out.insert(3, col_name, pd.NA)  # type: ignore

    for q_char in qseq:
        if q_char != "-":
            if aln_index in scorecons_df.index:
                score = scorecons_df.loc[aln_index, "conservation"]
                out.loc[out["Position"] == orig_pos, col_name] = score

            orig_pos += 1

        aln_index += 1

    return out


# ============================== Main Stage ==============================
def run_scorecons_reconstruct(**kwargs) -> None:
    """
    Stage: integrate Scorecons conservation scores into reconstructed PSSM tables.

    Required Args
    -------------
    pssm_fasta_path
    pssm_cdsearch_table
    pssm_reconstruct_dir
    scorecons_dir
    pssm_integrated_output_dir
    """
    start = time.time()

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║           [ Scorecons Integration Stage Started ]                   ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    fasta_path = kwargs.get("pssm_fasta_path")
    cdsearch_path = kwargs.get("pssm_cdsearch_table")

    scorecons_dir = kwargs.get("scorecons_dir")
    integrated_output_dir = kwargs.get("pssm_integrated_output_dir")

    pssm_reconstruct_dir = kwargs.get("pssm_reconstruct_dir")

    if not fasta_path:
        raise ValueError("pssm_fasta_path must be provided")
    if not cdsearch_path:
        raise ValueError("pssm_cdsearch_table must be provided")
    if not scorecons_dir:
        raise ValueError("scorecons_dir must be provided")
    if not integrated_output_dir:
        raise ValueError("pssm_integrated_output_dir must be provided")

    # ------------------------------------------------------
    # Resolve reconstruct dir if not explicitly provided
    # ------------------------------------------------------
    if not pssm_reconstruct_dir:
        branch = kwargs.get("branch", "psiblast")

        if branch not in ("psiblast", "smp"):
            raise ValueError(
                "pssm_reconstruct_dir is not provided, and branch is invalid. "
                "Branch must be 'psiblast' or 'smp'."
            )

        pssm_reconstruct_dir = os.path.join(
            os.path.dirname(os.path.dirname(integrated_output_dir)),
            "reconstruct",
            branch,
        )

    if not os.path.exists(pssm_reconstruct_dir):
        raise FileNotFoundError(
            f"PSSM reconstruct directory not found: {pssm_reconstruct_dir}"
        )

    os.makedirs(integrated_output_dir, exist_ok=True)

    log_path = os.path.join(integrated_output_dir, "scorecons_integrate_error.log")

    print(f"📂 FASTA Path             : {fasta_path}")
    print(f"📄 CD-Search Table        : {cdsearch_path}")
    print(f"📁 PSSM Reconstruct Dir   : {pssm_reconstruct_dir}")
    print(f"📁 Scorecons Dir          : {scorecons_dir}")
    print(f"📁 Integrated Output Dir  : {integrated_output_dir}")
    print(f"🐛 Error Log              : {log_path}\n")

    aln_table = _load_alignment_table(cdsearch_path)
    original_seqs = _load_original_sequences(fasta_path)

    total = len(aln_table)
    success = 0
    skipped = 0
    failed = 0

    for idx, row in aln_table.iterrows():
        query_id = row["query_id"]

        progress = (idx + 1) / total * 100  # type: ignore

        print("─" * 70)
        print(f"▶️  [{idx+1:3d}/{total:<3d} | {progress:5.1f}% ] {query_id}")  # type: ignore
        print("─" * 70)

        out_path = os.path.join(integrated_output_dir, f"{query_id}.tsv")

        if os.path.exists(out_path):
            skipped += 1
            print(f"   ⏩ Skip (already integrated): {out_path}\n")
            continue

        try:
            if query_id not in original_seqs:
                raise KeyError(f"Query sequence not found in FASTA: {query_id}")

            base_table = _load_pssm_reconstruct_table(pssm_reconstruct_dir, query_id)
            scorecons_path = _resolve_scorecons_path(scorecons_dir, query_id)
            scorecons_df = _parse_scorecons_txt(scorecons_path)

            merged = _integrate_scorecons(
                row,
                scorecons_df,
                base_table,
            )

            merged = _restore_numeric_types(merged, PSSM_NUMERIC_COLS)
            merged.to_csv(out_path, sep="\t", index=False)

            print(f"   ✅ Scorecons integrated -> {out_path}\n")
            success += 1

        except Exception as e:
            failed += 1
            print(f"   ❌ Failed: {e}\n")

            with open(log_path, "a") as log:
                log.write(f"[{query_id}] {e}\n")

    elapsed = time.time() - start

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║                 Scorecons Integration Summary                       ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    print(f"Total proteins : {total}")
    print(f"Successful     : {success}")
    print(f"Skipped        : {skipped}")
    print(f"Failed         : {failed}")
    print(f"🐛 Error log   : {log_path}")
    print(f"⏱️ Elapsed time: {elapsed:.2f} seconds\n")

    print("✅  Scorecons integration completed!\n")
