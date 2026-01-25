# src/preprocess/run_conservation_reconstruct.py
"""
Full-length Conservation Score Reconstruction Stage
(Project alignment-based Scorecons results back to protein positions)
"""

# ============================== Standard Library Imports ==============================
import os
import time

# ============================== Third-Party Imports ==============================
import pandas as pd
from Bio import SeqIO


# ============================== Helper Functions ==============================
def _parse_scorecons_txt(path: str) -> pd.DataFrame:
    """
    Parse Scorecons Server output file into a structured DataFrame.
    """

    if not os.path.exists(path):
        raise FileNotFoundError(f"Scorecons file not found: {path}")

    records = []

    with open(path, "r") as f:
        lines = f.readlines()

    # Scorecons output always uses the first 9 lines as header
    data_lines = lines[9:]

    aln_index = 1

    for line in data_lines:
        line = line.strip()
        if not line:
            continue

        try:
            left, right = line.split("#", 1)

            conservation = float(left.strip())
            residues = right.strip()

            # --- split residues ---
            res1 = residues[0] if len(residues) > 0 else None
            res2 = residues[1] if len(residues) > 1 else None

            records.append(
                {
                    "aln_index": aln_index,
                    "conservation": conservation,
                    "res_1": res1,
                    "res_2": res2,
                }
            )

            aln_index += 1

        except Exception:
            continue

    if not records:
        raise ValueError(f"No valid conservation scores found: {path}")

    return pd.DataFrame(records)


def _load_alignment_table(path: str) -> pd.DataFrame:
    """Load cdsearch_top_hits_detailed.tsv."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"CD-search table not found: {path}")
    return pd.read_csv(path, sep="\t")


def _load_original_sequences(fasta_path: str) -> dict:
    """Return {seq_id: sequence} dict."""
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Input FASTA not found: {fasta_path}")
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}


def _load_pssm_reconstruct_table(pssm_dir: str, query_id: str) -> pd.DataFrame:
    """Load reconstructed full-length PSSM table."""
    fpath = os.path.join(pssm_dir, f"{query_id}.tsv")
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"PSSM reconstruct file not found: {fpath}")
    return pd.read_csv(fpath, sep="\t")


def _reconstruct_full_conservation(
    original_seq: str,
    aln_row: pd.Series,
    conservation_df: pd.DataFrame,
    base_table: pd.DataFrame,
) -> pd.DataFrame:
    """
    Project alignment-based conservation scores back to full-length protein.

    Rules
    -----
    - Alignment column where query == '-' cannot be projected
    - Conservation score is written only when query residue exists
    - Unaligned positions remain NaN
    """

    # Ensure alignment-index lookup
    conservation_df = conservation_df.set_index("aln_index")

    qseq = aln_row["qseq"]
    hseq = aln_row["hseq"]

    orig_pos = aln_row["qstart"]  # 1-based
    aln_index = 1

    out = base_table.copy()

    # extract short title (before first comma)
    title_raw = aln_row["title"]
    short_title = title_raw.split(",")[0].strip()
    col_name = f"Conservation ({short_title})"

    # Insert as 4th column (after Alignment)
    out.insert(3, col_name, pd.NA)  # type: ignore

    for q_char, h_char in zip(qseq, hseq):

        # Case 1: query residue exists → projectable
        if q_char != "-":
            if aln_index in conservation_df.index:
                score = conservation_df.loc[aln_index, "conservation"]
                out.loc[out["Position"] == orig_pos, col_name] = score

            orig_pos += 1

        # Case 2: query gap → cannot be projected
        # (alignment insertion)
        else:
            pass

        aln_index += 1

    return out


# ============================== Main Stage Function ==============================
def run_conservation_reconstruct(**kwargs):
    """
    Stage: reconstruct full-length conservation score table using alignment mapping.
    """
    start = time.time()
    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║        [ Full-length Conservation Reconstruction Stage Started ]     ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    fasta_path = kwargs.get("conservation_fasta_path")
    cdsearch_path = kwargs.get("conservation_cdsearch_table")
    pssm_reconstruct_dir = kwargs.get("pssm_reconstruct_dir")
    conservation_root = kwargs.get("conservation_dir")
    scorecons_dir = os.path.join(conservation_root, "scorecons")  # type: ignore
    output_dir = os.path.join(conservation_root, "reconstruct")  # type: ignore

    os.makedirs(output_dir, exist_ok=True)  # type: ignore

    print(f"📂 FASTA Path            : {fasta_path}")
    print(f"📄 CD-Search Table       : {cdsearch_path}")
    print(f"📁 PSSM Reconstruct Dir  : {pssm_reconstruct_dir}")
    print(f"📁 Conservation ROOT     : {conservation_root}")
    print(f"📁 Output Directory      : {output_dir}\n")

    # Load data
    aln_table = _load_alignment_table(cdsearch_path)  # type: ignore
    original_seqs = _load_original_sequences(fasta_path)  # type: ignore

    total = len(aln_table)
    success = 0
    failed = 0

    # ===================== Reconstruction Loop =====================
    for idx, row in aln_table.iterrows():
        query_id = row["query_id"]

        print("─" * 70)
        progress = (idx + 1) / total * 100  # type: ignore
        print(f"▶️  [{idx+1:3d} / {total:<3d} | {progress:5.1f}% ]  {query_id}")  # type: ignore
        print("─" * 70)

        try:
            # Load base full-length table (PSSM reconstructed)
            base_table = _load_pssm_reconstruct_table(pssm_reconstruct_dir, query_id)  # type: ignore

            # Load conservation scorecons result
            scorecons_path = os.path.join(scorecons_dir, f"{query_id}.txt")
            conservation_df = _parse_scorecons_txt(scorecons_path)

            # Original sequence
            original_seq = original_seqs[query_id]

            # Reconstruct
            full = _reconstruct_full_conservation(
                original_seq,
                row,
                conservation_df,
                base_table,
            )

            # Save
            out_path = os.path.join(output_dir, f"{query_id}.tsv")  # type: ignore
            full.to_csv(out_path, sep="\t", index=False)

            print(f"   ✅ Conservation reconstructed for {query_id}\n")
            success += 1

        except Exception as e:
            print(f"   ❌ Reconstruction failed for {query_id}: {e}\n")
            failed += 1

    # ===================== Summary =====================
    elapsed = time.time() - start
    print("\n" + "╔" + "═" * 70 + "╗")
    print("║" + " " * 20 + "Conservation Reconstruction Summary" + " " * 19 + "║")
    print("╚" + "═" * 70 + "╝\n")

    print("📊 Summary Statistics")
    print("─" * 70)
    print(f"Total proteins : {total}")
    print(f"Successful     : {success}")
    print(f"Failed         : {failed}")
    print(f"⏱️ Elapsed time: {elapsed:.2f} seconds\n")

    print("✅  Conservation reconstruction stage completed!\n")
