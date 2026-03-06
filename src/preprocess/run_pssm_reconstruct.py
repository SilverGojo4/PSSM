# src/preprocess/run_pssm_reconstruct.py
"""
Full-length PSSM Reconstruction Stage

This stage projects domain-level PSSM matrices back to full-length protein
coordinates using CD-Search alignment information.

Updated Design (Branch-aware)
-----------------------------
This stage now supports decoupled input/output paths:

Input:
    --pssm_matrix_dir (directory containing domain-level *.tsv matrices)

Output:
    --pssm_reconstruct_output_dir (directory containing reconstructed *.tsv)

Error Log:
    <pssm_reconstruct_output_dir>/pssm_reconstruct_error.log

Legacy Support:
    --pssm_root_dir is still supported for backward compatibility.
"""

# ============================== Standard Library Imports ==============================
import os
import time
from typing import Optional, Tuple

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


def _resolve_paths(
    *,
    pssm_root_dir: Optional[str],
    pssm_matrix_dir: Optional[str],
    pssm_reconstruct_output_dir: Optional[str],
) -> Tuple[str, str]:
    """
    Resolve matrix input dir and reconstruction output dir.

    Priority:
    - New interface: pssm_matrix_dir + pssm_reconstruct_output_dir
    - Legacy interface: pssm_root_dir/pssm_matrices + pssm_root_dir/pssm_reconstruct
    """
    if pssm_matrix_dir and pssm_reconstruct_output_dir:
        return pssm_matrix_dir, pssm_reconstruct_output_dir

    if pssm_root_dir:
        matrix_dir = os.path.join(pssm_root_dir, "pssm_matrices")
        out_dir = os.path.join(pssm_root_dir, "pssm_reconstruct")
        return matrix_dir, out_dir

    raise ValueError(
        "Must provide either:\n"
        "  (1) pssm_matrix_dir + pssm_reconstruct_output_dir\n"
        "or\n"
        "  (2) pssm_root_dir (legacy mode)"
    )


def _load_domain_pssm_matrix(
    matrix_dir: str, query_id: str, pssm_id: int
) -> pd.DataFrame:
    """
    Load domain-level PSSM matrix.

    Expected file naming convention:
        {query_id}_{pssm_id}_hseq_with_gap.tsv

    Example:
        P00001_231391_hseq_with_gap.tsv
    """
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
    - Output uses UniProt absolute coordinate system (1..L)
    """
    L = len(original_seq)

    columns = (
        ["Position", "Residue", "Alignment"]
        + AA_ORDER
        + ["Po", "Hy", "Ch", "|Hy-Ch|", "|Hy-Po|"]
    )

    out = pd.DataFrame(index=range(1, L + 1), columns=columns)

    out["Position"] = range(1, L + 1)
    out["Residue"] = list(original_seq)
    out["Alignment"] = "0 (unaligned)"

    qseq = aln_row["qseq"]
    hseq = aln_row["hseq"]

    orig_pos = int(aln_row["qstart"])
    pssm_row_index = 0

    for q_char, h_char in zip(qseq, hseq):
        q_gap = q_char == "-"
        h_gap = h_char == "-"

        if not q_gap and not h_gap:
            if pssm_row_index < len(domain_pssm):
                row = domain_pssm.iloc[pssm_row_index]

                values = row[AA_ORDER + ["Po", "Hy", "Ch", "|Hy-Ch|", "|Hy-Po|"]].values
                out.loc[orig_pos, columns[3:]] = values
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

    New Directory Layout (recommended):
        input:  results/pssm/matrices/psiblast/*.tsv
        output: results/pssm/reconstruct/psiblast/*.tsv

    Legacy Layout (still supported):
        input:  <pssm_root_dir>/pssm_matrices/*.tsv
        output: <pssm_root_dir>/pssm_reconstruct/*.tsv
    """
    start = time.time()

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║              [ Full-length PSSM Reconstruction Stage Started ]       ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    pssm_root_dir = kwargs.get("pssm_root_dir")

    fasta_path = kwargs.get("pssm_fasta_path")
    cdsearch_path = kwargs.get("pssm_cdsearch_table")

    pssm_matrix_dir = kwargs.get("pssm_matrix_dir")
    reconstruct_output_dir = kwargs.get("pssm_reconstruct_output_dir")

    if not fasta_path:
        raise ValueError("pssm_fasta_path must be provided")
    if not cdsearch_path:
        raise ValueError("pssm_cdsearch_table must be provided")

    matrix_dir, output_dir = _resolve_paths(
        pssm_root_dir=pssm_root_dir,
        pssm_matrix_dir=pssm_matrix_dir,
        pssm_reconstruct_output_dir=reconstruct_output_dir,
    )

    os.makedirs(output_dir, exist_ok=True)

    log_path = os.path.join(output_dir, "pssm_reconstruct_error.log")

    print(f"📂 PSSM Matrix Dir     : {matrix_dir}")
    print(f"📁 Reconstruct Out Dir : {output_dir}")
    print(f"🐛 Error Log           : {log_path}")
    print(f"📂 FASTA Path          : {fasta_path}")
    print(f"📄 CD-Search Table     : {cdsearch_path}\n")

    aln_table = _load_alignment_table(cdsearch_path)
    original_seqs = _load_original_sequences(fasta_path)

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
            if query_id not in original_seqs:
                raise KeyError(f"Query sequence not found in FASTA: {query_id}")

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
