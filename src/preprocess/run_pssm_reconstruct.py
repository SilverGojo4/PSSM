# src/preprocess/run_pssm_reconstruct.py
"""
Full-length PSSM Reconstruction Stage

This stage projects domain-level PSSM matrices back to full-length protein
coordinates using CD-Search alignment information.

Directory layout (under pssm_root_dir):

domain_psiblast/
├── pssm_matrices/        # input (domain-level)
├── pssm_reconstruct/     # output (full-length)
└── pssm_reconstruct_error.log
"""

# ============================== Standard Library Imports ==============================
import os
import time

# ============================== Third-Party Imports ==============================
import pandas as pd
from Bio import SeqIO

# ============================== Constants ==============================
AA_ORDER = list("GAILVMFWPCSTYNQHKRDE")


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
    """Load domain-level PSSM matrix."""
    fname = f"{query_id}_{pssm_id}_hseq_with_gap.tsv"
    fpath = os.path.join(matrix_dir, fname)

    if not os.path.exists(fpath):
        raise FileNotFoundError(f"PSSM matrix not found: {fpath}")

    return pd.read_csv(fpath, sep="\t")


def _reconstruct_full_pssm(
    original_seq: str, aln_row: pd.Series, domain_pssm: pd.DataFrame
) -> pd.DataFrame:
    """
    Reconstruct full-length PSSM matrix using alignment mapping.

    Rules
    -----
    - Only aligned residues receive PSSM values
    - Unaligned positions remain NaN
    """

    L = len(original_seq)

    columns = (
        ["Position", "Residue", "Alignment"]
        + AA_ORDER
        + ["Po", "Hy", "Ch", "Hy+Ch-Po", "|Hy-Ch|"]
    )

    out = pd.DataFrame(index=range(1, L + 1), columns=columns)

    out["Position"] = range(1, L + 1)
    out["Residue"] = list(original_seq)
    out["Alignment"] = "0 (unaligned)"

    qseq = aln_row["qseq"]
    hseq = aln_row["hseq"]

    orig_pos = aln_row["qstart"]
    pssm_row_index = 0

    for q_char, h_char in zip(qseq, hseq):
        q_gap = q_char == "-"
        h_gap = h_char == "-"

        if not q_gap and not h_gap:
            if pssm_row_index < len(domain_pssm):
                row = domain_pssm.iloc[pssm_row_index]

                numeric_values = []
                for val in row[AA_ORDER + ["Po", "Hy", "Ch", "Hy+Ch-Po", "|Hy-Ch|"]]:
                    try:
                        numeric_values.append(int(float(val)))
                    except Exception:
                        numeric_values.append(val)

                out.loc[orig_pos, columns[3:]] = numeric_values
                out.loc[orig_pos, "Alignment"] = "1 (aligned)"

            pssm_row_index += 1
            orig_pos += 1

        elif not q_gap and h_gap:
            orig_pos += 1

        elif q_gap and not h_gap:
            pssm_row_index += 1

    return out


# ============================== Main Stage Function ==============================
def run_pssm_reconstruct(**kwargs):
    """
    Stage: reconstruct full-length PSSM matrices using alignment projection.
    """
    start = time.time()

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║              [ Full-length PSSM Reconstruction Stage Started ]       ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    pssm_root = kwargs.get("pssm_root_dir")
    fasta_path = kwargs.get("pssm_fasta_path")
    cdsearch_path = kwargs.get("pssm_cdsearch_table")

    if not pssm_root:
        raise ValueError("pssm_root_dir must be provided")

    matrix_dir = os.path.join(pssm_root, "pssm_matrices")
    output_dir = os.path.join(pssm_root, "pssm_reconstruct")
    log_path = os.path.join(pssm_root, "pssm_reconstruct_error.log")

    os.makedirs(output_dir, exist_ok=True)

    print(f"📁 PSI-BLAST Root Dir : {pssm_root}")
    print(f"📂 PSSM Matrices     : {matrix_dir}")
    print(f"📁 Reconstruct Out   : {output_dir}")
    print(f"🐛 Error Log         : {log_path}")
    print(f"📂 FASTA Path        : {fasta_path}")
    print(f"📄 CD-Search Table   : {cdsearch_path}\n")

    aln_table = _load_alignment_table(cdsearch_path)  # type: ignore
    original_seqs = _load_original_sequences(fasta_path)  # type: ignore

    total = len(aln_table)
    success, failed = 0, 0

    # ===================== Reconstruction Loop =====================
    for idx, row in aln_table.iterrows():
        query_id = row["query_id"]
        pssm_id = int(row["PSSM_ID"])

        progress = (idx + 1) / total * 100  # type: ignore

        print("─" * 70)
        print(
            f"▶️  [{idx+1:3d}/{total:<3d} | {progress:5.1f}% ] "  # type: ignore
            f"{query_id} (PSSM {pssm_id})"
        )
        print("─" * 70)

        try:
            domain_pssm = _load_domain_pssm_matrix(matrix_dir, query_id, pssm_id)

            original_seq = original_seqs[query_id]

            full = _reconstruct_full_pssm(original_seq, row, domain_pssm)

            out_path = os.path.join(output_dir, f"{query_id}.tsv")
            full.to_csv(out_path, sep="\t", index=False)

            success += 1
            print(f"   ✅ Reconstructed successfully for {query_id}\n")

        except Exception as e:
            failed += 1
            print(f"   ❌ Reconstruction failed for {query_id}: {e}\n")

            with open(log_path, "a") as log:
                log.write(f"[{query_id}] {e}\n")

    # ===================== Summary =====================
    elapsed = time.time() - start

    print("\n" + "╔" + "═" * 70 + "╗")
    print("║" + " " * 22 + "PSSM Reconstruction Summary" + " " * 20 + "║")
    print("╚" + "═" * 70 + "╝\n")

    print("📊 Summary Statistics")
    print("─" * 70)
    print(f"Total proteins : {total}")
    print(f"Successful     : {success}")
    print(f"Failed         : {failed}")
    print(f"🐛 Error log   : {log_path}")
    print(f"⏱️  Elapsed    : {elapsed:.2f} seconds\n")

    print("✅  Full-length PSSM reconstruction stage completed!\n")
