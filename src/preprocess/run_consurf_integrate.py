# src/preprocess/run_consurf_integrate.py
"""
ConSurf Conservation Integration Stage

Integrate residue-level ConSurf conservation scores into
integrated full-length PSSM tables.

Updated Design (Branch-aware Integrated Output)
-----------------------------------------------
This stage no longer reads/writes separate conservation/reconstruct folders.

Instead:
- Input: integrated tables from results/pssm/integrated/<branch>/
- Input: raw ConSurf grades files from results/consurf/
- Output: updated integrated tables (same directory)

Design principles
-----------------
- No alignment
- No gaps
- Variant-aware
- Window-based local sequence mapping
- Reference sequence (UniProt) is absolute coordinate system
"""

# ============================== Standard Library Imports ==============================
import os
import time
from glob import glob

# ============================== Third-Party Imports ==============================
import numpy as np
import pandas as pd
from Bio import SeqIO

# ============================== Constants ==============================
AA_ORDER = list("GAILVMFWPCSTYNQHKRDE")
PSSM_NUMERIC_COLS = AA_ORDER + ["Po", "Hy", "Ch", "|Hy-Ch|", "|Hy-Po|"]


# ============================== Helper Functions ==============================
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


def _load_base_integrated_table(path: str) -> pd.DataFrame:
    """
    Load integrated TSV table (PSSM + Scorecons).

    Important
    ---------
    We must restore numeric dtypes after loading.
    Otherwise integer PSSM matrices will be auto-cast to float when NA exists.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Integrated table not found: {path}")

    df = pd.read_csv(path, sep="\t")
    df = _restore_numeric_types(df, PSSM_NUMERIC_COLS)

    return df


def _parse_consurf_grades(path: str) -> pd.DataFrame:
    """
    Parse ConSurf grades file.

    Returns
    -------
    DataFrame with columns:
        pos, aa, score
    """
    records = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            try:
                pos = int(parts[0])
                aa = parts[1]
            except Exception:
                continue

            score = None
            for token in parts[2:]:
                try:
                    score = float(token)
                    break
                except Exception:
                    continue

            if score is None:
                continue

            records.append(
                {
                    "pos": pos,
                    "aa": aa,
                    "score": score,
                }
            )

    if not records:
        raise ValueError(f"No valid ConSurf scores found: {path}")

    return pd.DataFrame(records)


def _build_consurf_sequence(consurf_df: pd.DataFrame) -> str:
    """
    Build AA sequence string from ConSurf parsed table.
    """
    return "".join(consurf_df.sort_values("pos")["aa"].tolist())


def _parse_mutation(query_id: str):
    """
    Parse mutation info from query_id if present.

    Expected format:
        UniProt_A123B
    """
    try:
        mut = query_id.split("_")[1]
        from_aa = mut[0]
        to_aa = mut[-1]
        pos = int(mut[1:-1])
        return pos, from_aa, to_aa
    except Exception:
        return None, None, None


def _restore_reference_sequence(seq: str, mut_pos, mut_from):
    """
    Restore original residue in mutated sequence if mutation exists.

    ConSurf grades are often computed on the reference sequence.
    """
    if mut_pos is None:
        return seq

    idx = mut_pos - 1
    if 0 <= idx < len(seq):
        return seq[:idx] + mut_from + seq[idx + 1 :]

    return seq


# ============================== Core Mapping Logic ==============================
def _find_best_local_mapping(
    ref_seq: str,
    consurf_seq: str,
    window: int = 25,
    max_mismatch: int = 5,
):
    """
    Find best local matching region between two sequences.

    Returns
    -------
    ref_start : int
        0-based start index in ref_seq
    consurf_start : int
        0-based start index in consurf_seq
    mismatch : int
        mismatch count in best window
    """
    best = None

    def scan(seq_a, seq_b, label):
        nonlocal best

        for i in range(len(seq_a) - window + 1):
            anchor = seq_a[i : i + window]

            for j in range(len(seq_b) - window + 1):
                compare = seq_b[j : j + window]

                mismatch = sum(a != b for a, b in zip(anchor, compare))

                if mismatch <= max_mismatch:
                    candidate = (i, j, mismatch, label)

                    if best is None or mismatch < best[2]:
                        best = candidate

    scan(ref_seq, consurf_seq, "ref_in_consurf")
    scan(consurf_seq, ref_seq, "consurf_in_ref")

    if best is None:
        raise ValueError("No reliable local mapping found")

    ref_i, consurf_i, mismatch, label = best  # type: ignore

    if label == "ref_in_consurf":
        ref_start = ref_i
        consurf_start = consurf_i
    else:
        ref_start = consurf_i
        consurf_start = ref_i

    return ref_start, consurf_start, mismatch


def _integrate_consurf(
    base_table: pd.DataFrame,
    consurf_df: pd.DataFrame,
    ref_start: int,
    consurf_start: int,
) -> pd.DataFrame:
    """
    Integrate ConSurf scores into integrated table.

    Adds / overwrites column:
        Evolutionary conservation
    """
    out = base_table.copy()

    col_name = "Evolutionary conservation"

    # Always overwrite if rerun
    out[col_name] = pd.NA

    for _, row in consurf_df.iterrows():
        consurf_pos = int(row["pos"]) - 1
        score = row["score"]

        ref_pos = ref_start + (consurf_pos - consurf_start) + 1

        if 1 <= ref_pos <= len(out):
            out.loc[out["Position"] == ref_pos, col_name] = score

    # move to 4th column
    cols = out.columns.tolist()
    cols.remove(col_name)
    cols.insert(3, col_name)
    out = out[cols]

    return out


def _find_consurf_grades_file(consurf_dir: str, query_id: str) -> str:
    """
    Find ConSurf grades file for query_id.

    Search priority:
    - query_id-based
    - UniProt-based

    Example patterns:
        {query_id}_*_consurf_grades.txt
        {uniprot}_*_consurf_grades.txt
    """
    uniprot = query_id.split("_")[0]

    patterns = [
        f"{query_id}_*_consurf_grades.txt",
        f"{uniprot}_*_consurf_grades.txt",
    ]

    for p in patterns:
        matches = glob(os.path.join(consurf_dir, p))
        if matches:
            return sorted(matches)[0]

    raise FileNotFoundError(f"No ConSurf grades file found for {query_id}")


# ============================== Main Stage ==============================
def run_consurf_integrate(**kwargs) -> None:
    """
    Stage: integrate ConSurf evolutionary conservation scores into integrated PSSM tables.

    Required Args
    -------------
    pssm_fasta_path
    pssm_integrated_dir
    consurf_dir
    """
    start = time.time()

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║          [ ConSurf Conservation Integration Stage Started ]          ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    fasta_path = kwargs.get("pssm_fasta_path")
    integrated_dir = kwargs.get("pssm_integrated_dir")
    consurf_dir = kwargs.get("consurf_dir")

    if not fasta_path:
        raise ValueError("pssm_fasta_path must be provided")
    if not integrated_dir:
        raise ValueError("pssm_integrated_dir must be provided")
    if not consurf_dir:
        raise ValueError("consurf_dir must be provided")

    if not os.path.exists(integrated_dir):
        raise FileNotFoundError(f"Integrated directory not found: {integrated_dir}")
    if not os.path.exists(consurf_dir):
        raise FileNotFoundError(f"ConSurf directory not found: {consurf_dir}")

    log_path = os.path.join(integrated_dir, "consurf_integrate_error.log")

    print(f"📂 FASTA Path           : {fasta_path}")
    print(f"📁 Integrated Dir       : {integrated_dir}")
    print(f"📁 ConSurf Dir          : {consurf_dir}")
    print(f"🐛 Error Log            : {log_path}\n")

    original_seqs = _load_original_sequences(fasta_path)

    integrated_files = sorted(glob(os.path.join(integrated_dir, "*.tsv")))
    integrated_files = [p for p in integrated_files if not p.endswith("_metadata.tsv")]

    total = len(integrated_files)

    if total == 0:
        raise ValueError(f"No integrated TSV files found in: {integrated_dir}")

    success = 0
    failed = 0

    for idx, path in enumerate(integrated_files):
        query_id = os.path.basename(path).replace(".tsv", "")
        progress = (idx + 1) / total * 100

        print("─" * 70)
        print(f"▶️  [{idx+1:3d}/{total:<3d} | {progress:5.1f}% ]  {query_id}")
        print("─" * 70)

        try:
            if query_id not in original_seqs:
                raise KeyError(f"Query sequence not found in FASTA: {query_id}")

            base_table = _load_base_integrated_table(path)
            consurf_path = _find_consurf_grades_file(consurf_dir, query_id)
            consurf_df = _parse_consurf_grades(consurf_path)

            ref_seq = original_seqs[query_id]
            mut_pos, mut_from, _ = _parse_mutation(query_id)
            restored_seq = _restore_reference_sequence(ref_seq, mut_pos, mut_from)

            consurf_seq = _build_consurf_sequence(consurf_df)

            ref_start, consurf_start, mismatch = _find_best_local_mapping(
                restored_seq,
                consurf_seq,
            )

            merged = _integrate_consurf(
                base_table,
                consurf_df,
                ref_start,
                consurf_start,
            )

            merged = _restore_numeric_types(merged, PSSM_NUMERIC_COLS)
            merged.to_csv(path, sep="\t", index=False)

            print(
                f"   ✅ mapped (ref_start={ref_start+1}, "
                f"consurf_start={consurf_start+1}, mismatch={mismatch})"
            )
            print(f"   🧬 integrated -> {path}\n")

            success += 1

        except Exception as e:
            failed += 1
            print(f"   ❌ Failed: {e}\n")

            with open(log_path, "a") as log:
                log.write(f"[{query_id}] {e}\n")

    elapsed = time.time() - start

    print("\n══════════════════════════════════════════════════════════════════════")
    print("ConSurf Integration Summary")
    print("──────────────────────────────────────────────────────────────────────")
    print(f"Total       : {total}")
    print(f"Successful  : {success}")
    print(f"Failed      : {failed}")
    print(f"🐛 Error log: {log_path}")
    print(f"Elapsed     : {elapsed:.2f} sec")
    print("══════════════════════════════════════════════════════════════════════\n")
