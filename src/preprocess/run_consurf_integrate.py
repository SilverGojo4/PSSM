"""
ConSurf Conservation Integration Stage

Integrate residue-level ConSurf conservation scores into
successfully reconstructed full-length conservation tables.

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
import pandas as pd
from Bio import SeqIO


# ============================== Helper Functions ==============================
def _load_original_sequences(fasta_path: str) -> dict:
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Input FASTA not found: {fasta_path}")
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}


def _load_base_reconstruct_table(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Reconstruct table not found: {path}")
    return pd.read_csv(path, sep="\t")


def _parse_consurf_grades(path: str) -> pd.DataFrame:
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
    return "".join(consurf_df.sort_values("pos")["aa"].tolist())


def _parse_mutation(query_id: str):
    try:
        mut = query_id.split("_")[1]
        from_aa = mut[0]
        to_aa = mut[-1]
        pos = int(mut[1:-1])
        return pos, from_aa, to_aa
    except Exception:
        return None, None, None


def _restore_reference_sequence(seq: str, mut_pos, mut_from):
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

    Returns:
        ref_start (0-based)
        consurf_start (0-based)
        mismatch
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

    # A in B
    scan(ref_seq, consurf_seq, "ref_in_consurf")

    # B in A
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
    out = base_table.copy()

    col_name = "Evolutionary conservation"
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

            merged.to_csv(
                os.path.join(output_dir, f"{query_id}.tsv"),
                sep="\t",
                index=False,
            )

            print(
                f"   ✅ mapped (ref_start={ref_start+1}, "
                f"consurf_start={consurf_start+1}, mismatch={mismatch})\n"
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
