"""
Scorecons Conservation Reconstruction Stage

Project alignment-based Scorecons conservation scores
back to full-length protein residue positions.
"""

# ============================== Standard Library Imports ==============================
import os
import time

# ============================== Third-Party Imports ==============================
import pandas as pd
from Bio import SeqIO


# ============================== Helper Functions ==============================
def _parse_scorecons_txt(path: str) -> pd.DataFrame:
    """Parse Scorecons Server output file."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Scorecons file not found: {path}")

    records = []

    with open(path) as f:
        lines = f.readlines()

    data_lines = lines[9:]  # header

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
    if not os.path.exists(path):
        raise FileNotFoundError(f"CD-search table not found: {path}")
    return pd.read_csv(path, sep="\t")


def _load_original_sequences(fasta_path: str) -> dict:
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Input FASTA not found: {fasta_path}")
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}


def _load_pssm_reconstruct_table(pssm_dir: str, query_id: str) -> pd.DataFrame:
    fpath = os.path.join(pssm_dir, f"{query_id}.tsv")
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"PSSM reconstruct file not found: {fpath}")
    return pd.read_csv(fpath, sep="\t")


def _reconstruct_full_scorecons(
    original_seq: str,
    aln_row: pd.Series,
    scorecons_df: pd.DataFrame,
    base_table: pd.DataFrame,
) -> pd.DataFrame:
    """
    Project alignment-based Scorecons scores back to full-length protein.

    Rules
    -----
    - Query gap positions are not projectable
    - Only query residues receive scores
    - Unaligned positions remain NaN
    """

    scorecons_df = scorecons_df.set_index("aln_index")

    qseq = aln_row["qseq"]
    hseq = aln_row["hseq"]

    orig_pos = aln_row["qstart"]  # 1-based
    aln_index = 1

    out = base_table.copy()

    title = aln_row["title"].split(",")[0].strip()
    col_name = f"Conservation ({title})"

    out.insert(3, col_name, pd.NA)  # type: ignore

    for q_char, h_char in zip(qseq, hseq):

        if q_char != "-":
            if aln_index in scorecons_df.index:
                score = scorecons_df.loc[aln_index, "conservation"]
                out.loc[out["Position"] == orig_pos, col_name] = score
            orig_pos += 1

        aln_index += 1

    return out


# ============================== Main Stage ==============================
def run_scorecons_reconstruct(**kwargs):
    """
    Stage: reconstruct full-length Scorecons conservation table.
    """
    start = time.time()

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║        [ Scorecons Reconstruction Stage Started ]                   ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    fasta_path = kwargs.get("conservation_fasta_path")
    cdsearch_path = kwargs.get("conservation_cdsearch_table")
    pssm_reconstruct_dir = kwargs.get("pssm_reconstruct_dir")
    conservation_root = kwargs.get("conservation_dir")

    scorecons_dir = os.path.join(conservation_root, "scorecons")  # type: ignore
    output_dir = os.path.join(conservation_root, "reconstruct")  # type: ignore

    os.makedirs(output_dir, exist_ok=True)

    print(f"📂 FASTA Path           : {fasta_path}")
    print(f"📄 CD-Search Table      : {cdsearch_path}")
    print(f"📁 PSSM Reconstruct Dir : {pssm_reconstruct_dir}")
    print(f"📁 Scorecons Dir        : {scorecons_dir}")
    print(f"📁 Output Dir           : {output_dir}\n")

    aln_table = _load_alignment_table(cdsearch_path)  # type: ignore
    original_seqs = _load_original_sequences(fasta_path)  # type: ignore

    total = len(aln_table)
    success = 0
    failed = 0

    for idx, row in aln_table.iterrows():
        query_id = row["query_id"]

        print("─" * 70)
        print(f"▶️  [{idx+1:3d}/{total:<3d}] {query_id}")  # type: ignore
        print("─" * 70)

        try:
            base_table = _load_pssm_reconstruct_table(pssm_reconstruct_dir, query_id)  # type: ignore
            scorecons_path = os.path.join(scorecons_dir, f"{query_id}.txt")
            scorecons_df = _parse_scorecons_txt(scorecons_path)

            original_seq = original_seqs[query_id]

            full = _reconstruct_full_scorecons(
                original_seq,
                row,
                scorecons_df,
                base_table,
            )

            out_path = os.path.join(output_dir, f"{query_id}.tsv")
            full.to_csv(out_path, sep="\t", index=False)

            print(f"   ✅ Scorecons reconstructed\n")
            success += 1

        except Exception as e:
            print(f"   ❌ Failed: {e}\n")
            failed += 1

    elapsed = time.time() - start

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║                 Scorecons Reconstruction Summary                    ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    print(f"Total proteins : {total}")
    print(f"Successful     : {success}")
    print(f"Failed         : {failed}")
    print(f"⏱️ Elapsed time: {elapsed:.2f} seconds\n")

    print("✅  Scorecons reconstruction completed!\n")
