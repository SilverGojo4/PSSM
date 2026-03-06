"""
Final conservation-based feature masking stage

Updated Design (Branch-specific Filtered Output)
------------------------------------------------
This stage operates on integrated PSSM tables:

Input:
    results/pssm/integrated/<branch>/*.tsv

Output:
    results/pssm/filtered/<branch>/*.tsv

All residue positions are retained.
PSSM profile features are masked (set to NA)
if conservation criteria are not satisfied.

Strict Mode
-----------
If required conservation columns are missing, the file is skipped
and the error is recorded into a log file.
"""

# ============================== Standard Library Imports ==============================
import os
import time
from glob import glob

# ============================== Third Party Imports ==============================
import numpy as np
import pandas as pd

# ============================== Constants ==============================
AA_ORDER = list("GAILVMFWPCSTYNQHKRDE")
FEATURE_COLS = ["Po", "Hy", "Ch", "|Hy-Ch|", "|Hy-Po|"]
PSSM_NUMERIC_COLS = AA_ORDER + FEATURE_COLS


# ============================== Helper Functions ==============================
def _detect_scorecons_column(df: pd.DataFrame) -> str:
    """
    Detect Scorecons conservation column.

    Expected format:
        Conservation (XXXX)

    Returns
    -------
    str
        Detected Scorecons column name.
    """
    matches = [
        col
        for col in df.columns
        if col.startswith("Conservation (") and col.endswith(")")
    ]

    if len(matches) == 0:
        raise ValueError("No Scorecons column found (expected: Conservation (xxxx))")

    if len(matches) > 1:
        raise ValueError(
            f"Multiple Scorecons columns found: {matches}. "
            "Ambiguous integration result."
        )

    return matches[0]


def _validate_required_columns(df: pd.DataFrame) -> tuple[str, list[str]]:
    """
    Validate required columns exist.

    Required
    --------
    - Position
    - Evolutionary conservation
    - Conservation (xxxx)
    - PSSM columns (AA_ORDER + FEATURE_COLS)

    Returns
    -------
    cons_col : str
        Scorecons column name
    pssm_cols : list[str]
        PSSM feature columns that will be masked
    """
    if "Position" not in df.columns:
        raise ValueError("Missing required column: Position")

    if "Evolutionary conservation" not in df.columns:
        raise ValueError("Missing required column: Evolutionary conservation")

    cons_col = _detect_scorecons_column(df)

    missing = [c for c in PSSM_NUMERIC_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required PSSM feature columns: {missing}")

    return cons_col, PSSM_NUMERIC_COLS


def _restore_numeric_types(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """
    Restore numeric dtypes:
    - If integer-like -> Int64
    - Otherwise float

    Notes
    -----
    Needed because pandas auto-upcasts int columns into float when NA exists.
    """
    for c in cols:
        if c not in df.columns:
            continue

        series = pd.to_numeric(df[c], errors="coerce")

        if series.isna().all():
            continue

        non_na = series.dropna()

        if np.isclose(non_na, non_na.round()).all():
            df[c] = series.round().astype("Int64")
        else:
            df[c] = series.astype(float)

    return df


def _mask_pssm_features(
    df: pd.DataFrame,
    ec_min: float,
    ec_max: float,
    cons_min: float,
    cons_max: float,
) -> pd.DataFrame:
    """
    Mask PSSM features where conservation thresholds are not satisfied.
    """
    cons_col, pssm_cols = _validate_required_columns(df)

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

    out = df.copy()
    out.loc[~valid_mask, pssm_cols] = pd.NA

    out = _restore_numeric_types(out, pssm_cols)

    return out


# ============================== Main Stage ==============================
def run_conservation_filter(**kwargs) -> None:
    """
    Stage: mask PSSM profile features by conservation thresholds.

    Required Args
    -------------
    pssm_integrated_dir
    pssm_filtered_output_dir
    ec_min
    ec_max
    cons_min
    cons_max
    """
    start = time.time()

    integrated_dir = kwargs.get("pssm_integrated_dir")
    filtered_dir = kwargs.get("pssm_filtered_output_dir")

    if not integrated_dir:
        raise ValueError("pssm_integrated_dir must be provided")
    if not filtered_dir:
        raise ValueError("pssm_filtered_output_dir must be provided")

    ec_min = float(kwargs.get("ec_min"))  # type: ignore
    ec_max = float(kwargs.get("ec_max"))  # type: ignore
    cons_min = float(kwargs.get("cons_min"))  # type: ignore
    cons_max = float(kwargs.get("cons_max"))  # type: ignore

    if not os.path.exists(integrated_dir):
        raise FileNotFoundError(f"Integrated directory not found: {integrated_dir}")

    os.makedirs(filtered_dir, exist_ok=True)

    log_path = os.path.join(filtered_dir, "conservation_filter_error.log")

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║     [ Conservation-Based PSSM Feature Masking Started ]              ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    print(f"📥 Integrated Dir : {integrated_dir}")
    print(f"📤 Filtered Dir   : {filtered_dir}")
    print(f"🐛 Error Log      : {log_path}")
    print(f"📌 EC range       : {ec_min} ~ {ec_max}")
    print(f"📌 CONS range     : {cons_min} ~ {cons_max}\n")

    files = sorted(glob(os.path.join(integrated_dir, "*.tsv")))
    files = [p for p in files if not p.endswith("_metadata.tsv")]

    total = len(files)

    if total == 0:
        print(f"⚠️ No TSV files found in: {integrated_dir}\n")
        return

    success = 0
    failed = 0

    for i, path in enumerate(files):
        query_id = os.path.basename(path).replace(".tsv", "")
        progress = (i + 1) / total * 100

        print("─" * 70)
        print(f"▶️  [{i+1:3d}/{total:<3d} | {progress:5.1f}% ] {query_id}")
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

            out_path = os.path.join(filtered_dir, f"{query_id}.tsv")
            masked.to_csv(out_path, sep="\t", index=False)

            print(f"   ✅ Filtered -> {out_path}\n")
            success += 1

        except Exception as e:
            failed += 1
            print(f"   ❌ Failed: {e}\n")

            with open(log_path, "a") as log:
                log.write(f"[{query_id}] {e}\n")

    elapsed = time.time() - start

    print("\n══════════════════════════════════════════════════════════════════════")
    print("Conservation Filtering Summary")
    print("──────────────────────────────────────────────────────────────────────")
    print(f"Total       : {total}")
    print(f"Successful  : {success}")
    print(f"Failed      : {failed}")
    print(f"🐛 Error log: {log_path}")
    print(f"Elapsed     : {elapsed:.2f} sec")
    print("══════════════════════════════════════════════════════════════════════\n")
