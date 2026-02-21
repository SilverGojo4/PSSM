# src/postprocess/run_mutation_site_screen.py

"""
Mutation Site Recall Benchmarking (Grid Search Version)

Extended:
    - Standard Recall (all mutation sites)
    - Polar-only Recall

Polar residues:
    {"S", "T", "Y", "N", "Q"}

Rules:
    - If no GT mutation sites → Recall = NA
    - If no Polar GT mutation sites → Polar Recall = NA
"""

# ============================== Standard Library ==============================
import os
import time
from glob import glob

# ============================== Third Party ==============================
import pandas as pd

# ============================== Constants ==============================
POLAR = {"S", "T", "Y", "N", "Q"}


# ============================== Helpers ==============================


def _extract_selected_sites(df, hychpo_min, abs_hych_min):
    required_cols = ["Position", "Residue", "Hy+Ch-Po", "|Hy-Ch|"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: '{col}'")

    valid_mask = df["Hy+Ch-Po"].notna() & df["|Hy-Ch|"].notna()

    pass_mask = (
        valid_mask & (df["Hy+Ch-Po"] >= hychpo_min) & (df["|Hy-Ch|"] >= abs_hych_min)
    )

    selected_df = df.loc[pass_mask, ["Position", "Residue"]]

    return [
        f"{row['Residue']}{int(row['Position'])}" for _, row in selected_df.iterrows()
    ]


def _load_known_mutation_sites(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Known mutation sites TSV not found: {path}")

    df = pd.read_csv(path, sep="\t")

    if "ID" not in df.columns or "Known Mutation Sites" not in df.columns:
        raise ValueError("Known mutation sites TSV missing required columns.")

    df["ID"] = df["ID"].astype(str)
    return df


# ============================== Main ==============================


def run_mutation_site_screen(**kwargs):

    start = time.time()

    base_path = kwargs.get("base_path")
    branch = kwargs.get("branch")
    hychpo_max = int(kwargs.get("hychpo_max"))  # type: ignore
    abs_hych_max = int(kwargs.get("abs_hych_max"))  # type: ignore
    known_path = kwargs.get("known_mutation_sites_tsv")

    if not base_path:
        raise ValueError("base_path must be provided.")
    if not branch:
        raise ValueError("branch must be provided.")

    input_dir = os.path.join(base_path, "results", "pssm", "filtered", branch)

    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"Filtered directory not found: {input_dir}")

    known_df = _load_known_mutation_sites(known_path)

    known_dict = {
        row["ID"]: set(
            [
                x.strip()
                for x in str(row["Known Mutation Sites"]).split(",")
                if x.strip()
            ]
        )
        for _, row in known_df.iterrows()
    }

    # Output base
    output_base = os.path.join(
        base_path, "results", "analysis", "mutation_site_screen", branch
    )

    threshold_grid_dir = os.path.join(output_base, "threshold_grid")
    os.makedirs(threshold_grid_dir, exist_ok=True)

    files = sorted(glob(os.path.join(input_dir, "*.tsv")))
    recall_rows = []

    print("\n╔════════════════════════════════════════════════════════════╗")
    print("║   Mutation Site Recall Benchmarking (Grid Search)        ║")
    print("╚════════════════════════════════════════════════════════════╝\n")

    for hychpo_min in range(1, hychpo_max + 1):
        for abs_hych_min in range(1, abs_hych_max + 1):

            print(f"Running: Hy+Ch-Po >= {hychpo_min}, |Hy-Ch| >= {abs_hych_min}")

            per_protein_rows = []

            total_tp_all = 0
            total_known_all = 0

            total_tp_polar = 0
            total_known_polar = 0

            for path in files:
                query_id = os.path.basename(path).replace(".tsv", "")
                df = pd.read_csv(path, sep="\t")

                selected_sites = set(
                    _extract_selected_sites(df, hychpo_min, abs_hych_min)
                )

                known_all = known_dict.get(query_id, set())
                known_polar = {s for s in known_all if s[0] in POLAR}

                # ===================== ALL RECALL =====================
                known_n = len(known_all)

                if known_n > 0:
                    tp_all = len(selected_sites & known_all)
                    recall_all = tp_all / known_n
                    total_tp_all += tp_all
                    total_known_all += known_n
                else:
                    recall_all = pd.NA

                # ===================== POLAR RECALL =====================
                polar_known_n = len(known_polar)

                if known_n == 0:
                    polar_recall = pd.NA
                    polar_status = "NO_GT"
                elif polar_known_n == 0:
                    polar_recall = pd.NA
                    polar_status = "NO_POLAR_GT"
                else:
                    tp_polar = len(selected_sites & known_polar)
                    polar_recall = tp_polar / polar_known_n
                    polar_status = "OK"

                    total_tp_polar += tp_polar
                    total_known_polar += polar_known_n

                per_protein_rows.append(
                    {
                        "ID": query_id,
                        "N Selected Sites": len(selected_sites),
                        "Selected Sites": ",".join(sorted(selected_sites)),
                        "Known Mutation Sites": ",".join(sorted(known_all)),
                        "Recall": recall_all,
                        "Known Polar Mutation Sites": ",".join(sorted(known_polar)),
                        "N Known Polar Sites": polar_known_n,
                        "Polar Recall": polar_recall,
                        "Polar Recall Status": polar_status,
                    }
                )

            # ===================== MICRO RECALL =====================
            micro_recall_all = (
                total_tp_all / total_known_all if total_known_all > 0 else pd.NA
            )

            micro_recall_polar = (
                total_tp_polar / total_known_polar if total_known_polar > 0 else pd.NA
            )

            threshold_dir = os.path.join(
                threshold_grid_dir,
                f"hychpo_{hychpo_min}__abshych_{abs_hych_min}",
            )
            os.makedirs(threshold_dir, exist_ok=True)

            per_protein_df = pd.DataFrame(per_protein_rows)
            per_protein_df.to_csv(
                os.path.join(threshold_dir, "per_protein_predictions.tsv"),
                sep="\t",
                index=False,
            )

            recall_rows.append(
                {
                    "Hy+Ch-Po": hychpo_min,
                    "|Hy-Ch|": abs_hych_min,
                    "Micro Recall": micro_recall_all,
                    "Micro Recall (Polar)": micro_recall_polar,
                }
            )

    recall_summary_df = pd.DataFrame(recall_rows)
    recall_summary_df = recall_summary_df.sort_values(
        by="Micro Recall", ascending=False
    )

    recall_summary_df.to_csv(
        os.path.join(output_base, "recall_summary.tsv"),
        sep="\t",
        index=False,
    )

    elapsed = time.time() - start

    print("\n════════════════════════════════════════════════════════════")
    print(f"Output directory: {output_base}")
    print(f"Elapsed: {elapsed:.2f} sec")
    print("════════════════════════════════════════════════════════════\n")
