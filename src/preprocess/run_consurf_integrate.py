"""
ConSurf Conservation Integration Stage

Integrate residue-level ConSurf conservation scores into
successfully reconstructed full-length conservation tables.

Design principles
-----------------
- No alignment
- No gaps
- Variant-aware
- Offset-free
- Window-based sequence mapping
- Reference sequence is the absolute coordinate system
"""

# ============================== Standard Library Imports ==============================
import os
import time
from glob import glob

# ============================== Third-Party Imports ==============================
import pandas as pd
from Bio import SeqIO


# ============================== Helper Functions ==============================
def _load_original_sequences(fasta_path: str) -> dict:
    """Return {seq_id: sequence} dict."""
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Input FASTA not found: {fasta_path}")
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}


def _load_base_reconstruct_table(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Reconstruct table not found: {path}")
    return pd.read_csv(path, sep="\t")


def _parse_consurf_grades(path: str) -> pd.DataFrame:
    """
    Parse ConSurf grades file.

    Return columns:
        pos (int)   - position in ConSurf sequence (1-based)
        aa (str)    - amino acid
        score (float)
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
    """Build linear ConSurf sequence string."""
    return "".join(consurf_df.sort_values("pos")["aa"].tolist())


def _parse_mutation(query_id: str):
    """
    Example:
        P01857_L125C  -> (125, 'L', 'C')
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
    """Replace mutated residue back to reference amino acid."""
    if mut_pos is None:
        return seq

    idx = mut_pos - 1
    if 0 <= idx < len(seq):
        return seq[:idx] + mut_from + seq[idx + 1 :]
    return seq


def _find_consurf_start(
    full_seq: str,
    consurf_seq: str,
    max_mismatch: int = 10,
):
    """
    Find where ConSurf sequence maps onto full-length sequence.

    Returns:
        start_index (0-based)
        mismatch_count
    """

    L = len(consurf_seq)

    best_start = None
    best_mismatch = float("inf")

    for i in range(len(full_seq) - L + 1):
        window = full_seq[i : i + L]

        mismatch = sum(1 for a, b in zip(window, consurf_seq) if a != b)

        if mismatch < best_mismatch:
            best_mismatch = mismatch
            best_start = i

    if best_start is None:
        raise ValueError("Unable to locate ConSurf region")

    return best_start, best_mismatch


def _integrate_consurf_with_start(
    base_table: pd.DataFrame,
    consurf_df: pd.DataFrame,
    start_index: int,
) -> pd.DataFrame:
    """
    Write evolutionary conservation scores into full-length table
    and move the column to the 4th position.
    """

    out = base_table.copy()

    col_name = "Evolutionary conservation"

    # create column
    out[col_name] = pd.NA

    # fill values
    for _, row in consurf_df.iterrows():
        consurf_pos = int(row["pos"])  # 1-based
        score = row["score"]

        full_pos = start_index + consurf_pos

        out.loc[out["Position"] == full_pos, col_name] = score

    # move column to 4th position
    cols = out.columns.tolist()

    cols.remove(col_name)
    cols.insert(3, col_name)  # index 3 = 4th column

    out = out[cols]

    return out


def _find_consurf_grades_file(consurf_dir: str, query_id: str) -> str:
    parts = query_id.split("_")
    uniprot_id = parts[0]

    variant_pattern = os.path.join(
        consurf_dir,
        f"{query_id}_*_consurf_grades.txt",
    )
    matches = glob(variant_pattern)
    if matches:
        return sorted(matches)[0]

    wt_pattern = os.path.join(
        consurf_dir,
        f"{uniprot_id}_*_consurf_grades.txt",
    )
    matches = glob(wt_pattern)
    if matches:
        return sorted(matches)[0]

    raise FileNotFoundError(f"No ConSurf grades file for {query_id}")


# ============================== Main Stage ==============================
def run_consurf_integrate(**kwargs):
    start = time.time()

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║          [ ConSurf Conservation Integration Stage Started ]          ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    fasta_path = kwargs.get("conservation_fasta_path")
    reconstruct_dir = kwargs.get("conservation_reconstruct_dir")
    consurf_root = kwargs.get("conservation_dir")

    grades_dir = os.path.join(consurf_root, "consurf")  # type: ignore
    output_dir = os.path.join(consurf_root, "reconstruct")  # type: ignore

    os.makedirs(output_dir, exist_ok=True)

    original_seqs = _load_original_sequences(fasta_path)  # type: ignore
    reconstruct_files = sorted(glob(os.path.join(reconstruct_dir, "*.tsv")))  # type: ignore

    total = len(reconstruct_files)
    success = 0
    failed = 0

    for idx, path in enumerate(reconstruct_files):
        query_id = os.path.basename(path).replace(".tsv", "")

        print("─" * 70)
        print(f"▶️  [{idx+1}/{total}]  {query_id}")
        print("─" * 70)

        try:
            base_table = _load_base_reconstruct_table(path)
            consurf_path = _find_consurf_grades_file(grades_dir, query_id)
            consurf_df = _parse_consurf_grades(consurf_path)

            ref_seq = original_seqs[query_id]

            mut_pos, mut_from, _ = _parse_mutation(query_id)
            restored_seq = _restore_reference_sequence(ref_seq, mut_pos, mut_from)

            consurf_seq = _build_consurf_sequence(consurf_df)

            start_index, mismatch = _find_consurf_start(restored_seq, consurf_seq)

            merged = _integrate_consurf_with_start(
                base_table,
                consurf_df,
                start_index,
            )

            out_path = os.path.join(output_dir, f"{query_id}.tsv")
            merged.to_csv(out_path, sep="\t", index=False)

            print(
                f"   ✅ ConSurf start = {start_index + 1} " f"(mismatch={mismatch})\n"
            )
            success += 1

        except Exception as e:
            print(f"   ❌ Failed: {e}\n")
            failed += 1

    elapsed = time.time() - start

    print("\n══════════════════════════════════════════════════════════════════════")
    print("ConSurf Integration Summary")
    print("──────────────────────────────────────────────────────────────────────")
    print(f"Total       : {total}")
    print(f"Successful  : {success}")
    print(f"Failed      : {failed}")
    print(f"Elapsed     : {elapsed:.2f} sec")
    print("══════════════════════════════════════════════════════════════════════\n")
