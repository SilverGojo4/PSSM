# src/preprocess/run_pssm_features.py
"""
PSSM Score Matrix Extraction Stage

This stage extracts domain-level PSSM score matrices from ASCII PSI-BLAST output
and computes biochemical group features (Po / Hy / Ch).
All outputs are organized under the domain_psiblast root directory.
"""

# ============================== Standard Library Imports ==============================
import os
import time
from glob import glob

# ============================== Third-Party Library Imports ==============================
import pandas as pd

# ============================== Constants ==============================
AA_ORDER = list("GAILVMFWPCSTYNQHKRDE")


# ============================== Helper Function ==============================
def _extract_score_matrix(pssm_file: str) -> pd.DataFrame:
    """
    Extract 20×L score matrix and biochemical group features from ASCII PSSM.
    """
    with open(pssm_file) as f:
        lines = f.readlines()

    start_idx = None
    for i, line in enumerate(lines):
        if "Last position-specific scoring matrix computed" in line:
            start_idx = i + 1
            break

    if start_idx is None:
        raise ValueError(f"Matrix start not found in {pssm_file}")

    HYDROPHOBIC = {"A", "I", "L", "V", "M", "F", "W", "P", "C"}
    CHARGED = {"H", "D", "E", "K", "R"}
    POLAR = {"S", "T", "Y", "N", "Q"}

    rows = []

    for line in lines[start_idx:]:
        s = line.strip()
        if not s or not s[0].isdigit():
            continue

        parts = s.split()
        if len(parts) < 22:
            continue

        pos = int(parts[0])
        aa = parts[1]
        scores = list(map(int, parts[2:22]))

        row = {"Position": pos, "Residue": aa}

        for sym, sc in zip(AA_ORDER, scores):
            row[sym] = sc

        po_sum = sum(sc for sym, sc in zip(AA_ORDER, scores) if sym in POLAR and sc > 0)
        hy_sum = sum(
            sc for sym, sc in zip(AA_ORDER, scores) if sym in HYDROPHOBIC and sc > 0
        )
        ch_sum = sum(
            sc for sym, sc in zip(AA_ORDER, scores) if sym in CHARGED and sc > 0
        )

        row.update(
            {
                "Po": po_sum,
                "Hy": hy_sum,
                "Ch": ch_sum,
                "Hy+Ch-Po": (hy_sum + ch_sum) - po_sum,
                "|Hy-Ch|": abs(hy_sum - ch_sum),
            }
        )

        rows.append(row)

    return pd.DataFrame(
        rows,
        columns=["Position", "Residue"]
        + AA_ORDER
        + ["Po", "Hy", "Ch", "Hy+Ch-Po", "|Hy-Ch|"],
    )


# ============================== Main Stage Function ==============================
def run_pssm_extract_matrix(**kwargs):
    """
    Stage: extract domain-level PSSM matrices from PSI-BLAST output.

    Directory layout:
        domain_psiblast/
        ├── pssm_profiles/     (input)
        ├── pssm_matrices/     (output)
        └── pssm_extract_error.log
    """
    start_time = time.time()

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║                [ PSSM Matrix Extraction Stage Started ]              ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    pssm_root = kwargs.get("pssm_root_dir")
    if not pssm_root:
        raise ValueError("pssm_root_dir must be provided")

    pssm_profiles_dir = os.path.join(pssm_root, "pssm_profiles")
    matrix_dir = os.path.join(pssm_root, "pssm_matrices")
    log_path = os.path.join(pssm_root, "pssm_extract_error.log")

    os.makedirs(matrix_dir, exist_ok=True)

    pssm_files = sorted(glob(os.path.join(pssm_profiles_dir, "*.pssm")))
    total = len(pssm_files)

    if total == 0:
        raise ValueError(f"No .pssm files found in {pssm_profiles_dir}")

    print(f"📁 PSI-BLAST Root Dir : {pssm_root}")
    print(f"📂 PSSM Profiles     : {pssm_profiles_dir}")
    print(f"📁 PSSM Matrices     : {matrix_dir}")
    print(f"🐛 Error Log         : {log_path}\n")

    success, failed = 0, 0

    # ===================== Extraction Loop =====================
    for idx, pssm_file in enumerate(pssm_files, start=1):
        domain_id = os.path.basename(pssm_file).replace(".pssm", "")
        progress = (idx / total) * 100

        print("─" * 70)
        print(f"▶️  [{idx:3d} / {total:<3d} | {progress:5.1f}% ]  {domain_id}")
        print("─" * 70)

        try:
            df = _extract_score_matrix(pssm_file)
            out_path = os.path.join(matrix_dir, f"{domain_id}.tsv")
            df.to_csv(out_path, sep="\t", index=False)

            success += 1
            print(f"   ✅ Matrix extracted successfully for {domain_id}\n")

        except Exception as e:
            failed += 1
            print(f"   ❌ Extraction failed for {domain_id}: {e}\n")

            with open(log_path, "a") as log:
                log.write(f"[{domain_id}] {e}\n")

    # ===================== Summary =====================
    elapsed = time.time() - start_time

    print("\n" + "╔" + "═" * 70 + "╗")
    print("║" + " " * 25 + "PSSM Extraction Summary" + " " * 22 + "║")
    print("╚" + "═" * 70 + "╝\n")

    print("📊 Summary Statistics")
    print("─" * 70)
    print(f"Total PSSM files : {total}")
    print(f"Successful       : {success}")
    print(f"Failed           : {failed}")
    print(f"🐛 Error log      : {log_path}")
    print(f"⏱️  Elapsed time  : {elapsed:.2f} seconds\n")

    print("✅  PSSM matrix extraction stage completed successfully!\n")
