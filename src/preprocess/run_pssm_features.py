# src/preprocess/run_pssm_features.py
"""
PSSM Score Matrix Extraction Stage
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
    """Extract only the 20xL score matrix (A..V) from ASCII PSSM."""
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
        hych_minus_po = (hy_sum + ch_sum) - po_sum
        hy_minus_ch_abs = abs(hy_sum - ch_sum)

        row.update(
            {
                "Po": po_sum,
                "Hy": hy_sum,
                "Ch": ch_sum,
                "Hy+Ch-Po": hych_minus_po,
                "|Hy-Ch|": hy_minus_ch_abs,
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
    Stage: Extract score matrices from all .pssm files.
    """
    start_time = time.time()
    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║                [ PSSM Matrix Extraction Stage Started ]              ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    pssm_dir = kwargs.get(
        "pssm_profiles_dir"
    )  # e.g. results/domain_psiblast/pssm_profiles
    output_dir = kwargs.get("pssm_matrix_output_dir")  # e.g. results/domain_psiblast
    os.makedirs(output_dir, exist_ok=True)  # type: ignore

    matrix_dir = os.path.join(output_dir, "pssm_matrices")  # type: ignore
    os.makedirs(matrix_dir, exist_ok=True)

    pssm_files = sorted(glob(os.path.join(pssm_dir, "*.pssm")))  # type: ignore
    total = len(pssm_files)

    print(f"📂 PSSM Input Directory : {pssm_dir}")
    print(f"📁 Output Directory     : {output_dir}")
    print(f"📄 Total PSSM files to extract : {total}\n")

    success, failed = 0, 0

    # ===================== Extraction Loop =====================
    for idx, p in enumerate(pssm_files, start=1):
        domain_id = os.path.basename(p).replace(".pssm", "")
        progress = (idx / total) * 100
        print("─" * 70)
        print(f"▶️  [{idx:3d} / {total:<3d} | {progress:5.1f}% ]  {domain_id}")
        print("─" * 70)

        try:
            df = _extract_score_matrix(p)
            out_path = os.path.join(matrix_dir, f"{domain_id}.tsv")
            df.to_csv(out_path, sep="\t", index=False)
            success += 1
            print(f"   ✅ Matrix extracted successfully for {domain_id}\n")

        except Exception as e:
            failed += 1
            print(f"   ❌ Extraction failed for {domain_id}\n")
            with open(os.path.join(output_dir, "pssm_extract_error.log"), "a") as log:  # type: ignore
                log.write(f"[{domain_id}] {e}\n")

    # ===================== Summary Output =====================
    elapsed = time.time() - start_time
    print("\n" + "╔" + "═" * 70 + "╗")
    print("║" + " " * 25 + "PSSM Extraction Summary" + " " * 22 + "║")
    print("╚" + "═" * 70 + "╝\n")

    print("📄 Output Files")
    print("─" * 70)
    print(f"📁 PSSM Matrices     : {matrix_dir}")
    print(
        f"🐛 Error Log (if any): {os.path.join(output_dir, 'pssm_extract_error.log')}\n"  # type: ignore
    )

    print("📊 Summary Statistics")
    print("─" * 70)
    print(f"Total PSSM files : {total}")
    print(f"Successful       : {success}")
    print(f"Failed           : {failed}\n")

    print(f"⏱️  Elapsed time : {elapsed:.2f} seconds")
    print("✅  PSSM Matrix extraction stage completed successfully!\n")
