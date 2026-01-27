"""
Final conservation-based feature masking stage

Input:
    results/evolutionary_conservation/reconstruct/*.tsv

Output:
    results/conservation_filtered/*.tsv

All residue positions are retained.
PSSM profile features are masked (set to NA)
if conservation criteria are not satisfied.
"""

# ============================== Standard Library ==============================
import os
import time
from glob import glob

# ============================== Third Party ==============================
import pandas as pd


# ============================== Helper ==============================
def _detect_conservation_column(df: pd.DataFrame) -> str:
    for col in df.columns:
        if col.startswith("Conservation (") and col.endswith(")"):
            return col
    raise ValueError("No Conservation (xxxx) column found")


def _detect_pssm_columns(df: pd.DataFrame, cons_col: str) -> list[str]:
    """
    All columns after the conservation column
    are treated as PSSM-derived features.
    """
    cols = df.columns.tolist()
    cons_idx = cols.index(cons_col)
    return cols[cons_idx + 1 :]


def _mask_pssm_features(
    df: pd.DataFrame,
    ec_min: float,
    ec_max: float,
    cons_min: float,
    cons_max: float,
) -> pd.DataFrame:

    cons_col = _detect_conservation_column(df)
    pssm_cols = _detect_pssm_columns(df, cons_col)

    evo = df["Evolutionary conservation"]
    cons = df[cons_col]

    valid_mask = (
        evo.notna()
        & cons.notna()
        & (evo >= ec_min)
        & (evo <= ec_max)
        & (cons >= cons_min)
        & (cons <= cons_max)
    )

    # mask PSSM features where condition is NOT satisfied
    df.loc[~valid_mask, pssm_cols] = pd.NA

    return df


# ============================== Main ==============================
def run_conservation_filter(**kwargs):
    start = time.time()

    reconstruct_dir = kwargs.get("conservation_reconstruct_dir")
    base_path = kwargs.get("base_path")

    # final output directory
    results_root = os.path.join(base_path, "results")  # type: ignore
    output_dir = os.path.join(results_root, "conservation_filtered")
    os.makedirs(output_dir, exist_ok=True)

    ec_min = float(kwargs.get("ec_min"))  # type: ignore
    ec_max = float(kwargs.get("ec_max"))  # type: ignore
    cons_min = float(kwargs.get("cons_min"))  # type: ignore
    cons_max = float(kwargs.get("cons_max"))  # type: ignore

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║     [ Conservation-Based PSSM Feature Masking Started ]              ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    files = sorted(glob(os.path.join(reconstruct_dir, "*.tsv")))  # type: ignore

    total, success, failed = len(files), 0, 0

    for i, path in enumerate(files):
        query_id = os.path.basename(path).replace(".tsv", "")
        print("─" * 70)
        print(f"▶️  [{i+1}/{total}] {query_id}")
        print("─" * 70)

        try:
            df = pd.read_csv(path, sep="\t")

            masked = _mask_pssm_features(
                df,
                ec_min=ec_min,
                ec_max=ec_max,
                cons_min=cons_min,
                cons_max=cons_max,
            )

            out_path = os.path.join(output_dir, f"{query_id}.tsv")
            masked.to_csv(out_path, sep="\t", index=False)

            kept = masked.iloc[:, :].notna().sum().sum()
            print(f"   ✅ PSSM masked where necessary\n")
            success += 1

        except Exception as e:
            print(f"   ❌ Failed: {e}\n")
            failed += 1

    elapsed = time.time() - start

    print("\n══════════════════════════════════════════════════════════════════════")
    print("Conservation Feature Masking Summary")
    print("──────────────────────────────────────────────────────────────────────")
    print(f"Total       : {total}")
    print(f"Successful  : {success}")
    print(f"Failed      : {failed}")
    print(f"Elapsed     : {elapsed:.2f} sec")
    print("══════════════════════════════════════════════════════════════════════\n")
