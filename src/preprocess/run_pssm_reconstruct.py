# src/preprocess/run_pssm_reconstruct.py
"""
Full-length PSSM Reconstruction Stage
"""

# ============================== Standard Library Imports ==============================
import os
import time

# ============================== Third-Party Imports ==============================
import pandas as pd
from Bio import SeqIO

# ============================== Constants ==============================
AA_ORDER = list("GAILVMFWPCSTYNQHKRDE")  # same order as extraction


# ============================== Helper Functions ==============================
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


def _load_domain_pssm_matrix(
    matrix_dir: str, query_id: str, pssm_id: int
) -> pd.DataFrame:
    """Load domain PSSM matrix extracted earlier."""
    fname = f"{query_id}_{pssm_id}_hseq_with_gap.tsv"
    fpath = os.path.join(matrix_dir, fname)

    if not os.path.exists(fpath):
        raise FileNotFoundError(f"PSSM matrix not found: {fpath}")

    return pd.read_csv(fpath, sep="\t")


def _reconstruct_full_pssm(
    original_seq: str, aln_row: pd.Series, domain_pssm: pd.DataFrame
) -> pd.DataFrame:
    """
    Create full-length L x (features) PSSM matrix using alignment mapping.
    All unmatched positions filled with NaN.

    Conservation:
        "aligned"   = aligned to PSSM row
        "unaligned" = no alignment match
    """
    L = len(original_seq)

    # Prepare empty output table
    columns = (
        ["Position", "Residue", "Conservation"]
        + AA_ORDER
        + ["Po", "Hy", "Ch", "Hy+Ch-Po", "|Hy-Ch|"]
    )
    out = pd.DataFrame(index=range(1, L + 1), columns=columns)

    out["Position"] = range(1, L + 1)
    out["Residue"] = list(original_seq)
    out["Conservation"] = "0 (unaligned)"  # default

    # Alignment strings
    qseq = aln_row["qseq"]
    hseq = aln_row["hseq"]

    orig_pos = aln_row["qstart"]  # original sequence index pointer (1-based)
    pssm_row_index = 0  # domain_pssm row pointer (0-based)

    for q_char, h_char in zip(qseq, hseq):
        q_is_gap = q_char == "-"
        h_is_gap = h_char == "-"

        # Case 1: both aligned → copy PSSM row & mark aligned
        if not q_is_gap and not h_is_gap:
            if pssm_row_index < len(domain_pssm):
                pssm_row = domain_pssm.iloc[pssm_row_index]

                # --------- ★ 這裡做強制整數轉換 ★ ---------
                numeric_values = []
                for val in pssm_row[
                    AA_ORDER + ["Po", "Hy", "Ch", "Hy+Ch-Po", "|Hy-Ch|"]
                ]:
                    try:
                        numeric_values.append(int(float(val)))  # e.g. "3.0" → 3
                    except:
                        numeric_values.append(val)

                out.loc[orig_pos, columns[3:]] = numeric_values
                out.loc[orig_pos, "Conservation"] = "1 (aligned)"

                # out.loc[orig_pos, columns[3:]] = pssm_row[
                #     AA_ORDER + ["Po", "Hy", "Ch", "Hy+Ch-Po", "|Hy-Ch|"]
                # ]
                # out.loc[orig_pos, "Conservation"] = "1 (aligned)"

            pssm_row_index += 1
            orig_pos += 1

        # Case 2: residue present in original seq but gap in hseq
        elif not q_is_gap and h_is_gap:
            orig_pos += 1

        # Case 3: gap in qseq but residue in hseq
        elif q_is_gap and not h_is_gap:
            pssm_row_index += 1

        # Case 4: both gaps — do nothing

    return out


# ============================== Main Stage Function ==============================
def run_pssm_reconstruct(**kwargs):
    """
    Stage: reconstruct full-length PSSM matrix for each domain-hit protein.
    """
    start = time.time()
    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║              [ Full-length PSSM Reconstruction Stage Started ]       ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    fasta_path = kwargs.get("pssm_fasta_path")  # original FASTA
    cdsearch_path = kwargs.get("pssm_cdsearch_table")  # cdsearch_top_hits_detailed.tsv
    matrix_dir = kwargs.get("pssm_matrix_dir")  # domain pssm_matrices
    output_dir = kwargs.get("pssm_reconstruct_output")  # target folder

    if not os.path.exists(output_dir):  # type: ignore
        os.makedirs(output_dir, exist_ok=True)  # type: ignore

    print(f"📂 FASTA Path        : {fasta_path}")
    print(f"📄 CD-Search Results : {cdsearch_path}")
    print(f"📁 Matrix Directory  : {matrix_dir}")
    print(f"📁 Output Directory  : {output_dir}\n")

    # Load data
    aln_table = _load_alignment_table(cdsearch_path)  # type: ignore
    original_seqs = _load_original_sequences(fasta_path)  # type: ignore

    total = len(aln_table)
    success = 0
    failed = 0

    # ===================== Reconstruction Loop =====================
    for idx, row in aln_table.iterrows():
        query_id = row["query_id"]
        pssm_id = int(row["PSSM_ID"])

        print("─" * 70)
        progress = (idx + 1) / total * 100  # type: ignore
        print(
            f"▶️  [{idx+1:3d} / {total:<3d} | {progress:5.1f}% ]  {query_id} (PSSM {pssm_id})"  # type: ignore
        )
        print("─" * 70)

        try:
            # Load domain PSSM for this query
            domain_pssm = _load_domain_pssm_matrix(matrix_dir, query_id, pssm_id)  # type: ignore

            # Original sequence
            original_seq = original_seqs[query_id]

            # Reconstruct
            full = _reconstruct_full_pssm(original_seq, row, domain_pssm)

            # Save
            out_path = os.path.join(output_dir, f"{query_id}_full_pssm.tsv")  # type: ignore
            full.to_csv(out_path, sep="\t", index=False)

            print(f"   ✅ Full-length PSSM reconstructed for {query_id}\n")
            success += 1

        except Exception as e:
            print(f"   ❌ Reconstruction failed for {query_id}: {e}\n")
            failed += 1

    # ===================== Summary =====================
    elapsed = time.time() - start
    print("\n" + "╔" + "═" * 70 + "╗")
    print("║" + " " * 22 + "Full-length PSSM Reconstruction Summary" + " " * 18 + "║")
    print("╚" + "═" * 70 + "╝\n")

    print("📊 Summary Statistics")
    print("─" * 70)
    print(f"Total proteins : {total}")
    print(f"Successful     : {success}")
    print(f"Failed         : {failed}")
    print(f"⏱️ Elapsed time: {elapsed:.2f} seconds\n")

    print("✅  Full-length PSSM reconstruction stage completed!\n")
