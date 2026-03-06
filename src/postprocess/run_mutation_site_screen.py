"""
Mutation Site Recall Benchmarking (Ultra Fast Parallel Grid Search)
"""

# ============================== Standard Library ==============================
import itertools
import os
import time
from glob import glob
from multiprocessing import Pool, cpu_count
from statistics import mean, median

# ============================== Third Party ==============================
import pandas as pd
from tqdm import tqdm

# ============================== Constants ==============================
POLAR = {"S", "T", "Y", "N", "Q"}

# ============================== Global Shared Objects ==============================
PROTEIN_TABLES = None
KNOWN_DICT = None


# ============================== Initializer ==============================


def _init_worker(protein_tables, known_dict):
    global PROTEIN_TABLES
    global KNOWN_DICT

    PROTEIN_TABLES = protein_tables
    KNOWN_DICT = known_dict


# ============================== Core Computation ==============================


def _extract_selected_sites(table, x, y, z, score_thr, feature, feature_thr):

    score = x * table["Hy"] + y * table["Ch"] + z * table["Po"]

    score_mask = score > score_thr

    if feature == "HyCh":
        feature_mask = table["HyCh"] > feature_thr
    else:
        feature_mask = table["HyPo"] > feature_thr

    mask = score_mask & feature_mask

    residues = table["Residue"][mask]
    positions = table["Position"][mask]

    return set(residues + positions)


def _evaluate_model(params):

    x, y, z, score_thr, feature, feature_thr = params

    total_tp_all = 0
    total_known_all = 0

    total_tp_polar = 0
    total_known_polar = 0

    selected_counts = []

    for query_id, table in PROTEIN_TABLES.items():  # type: ignore

        selected_sites = _extract_selected_sites(
            table,
            x,
            y,
            z,
            score_thr,
            feature,
            feature_thr,
        )

        selected_counts.append(len(selected_sites))

        known_all = KNOWN_DICT.get(query_id, set())  # type: ignore
        known_polar = {s for s in known_all if s[0] in POLAR}

        if known_all:
            tp_all = len(selected_sites & known_all)
            total_tp_all += tp_all
            total_known_all += len(known_all)

        if known_polar:
            tp_polar = len(selected_sites & known_polar)
            total_tp_polar += tp_polar
            total_known_polar += len(known_polar)

    micro_recall_all = total_tp_all / total_known_all if total_known_all > 0 else pd.NA

    micro_recall_polar = (
        total_tp_polar / total_known_polar if total_known_polar > 0 else pd.NA
    )

    return {
        "x": x,
        "y": y,
        "z": z,
        "score_thr": score_thr,
        "feature": feature,
        "feature_thr": feature_thr,
        "Micro Recall": micro_recall_all,
        "Micro Recall Polar": micro_recall_polar,
        "Mean N Selected Sites": mean(selected_counts),
        "Median N Selected Sites": median(selected_counts),
        "Formula Complexity": abs(x) + abs(y) + abs(z),
    }


# ============================== Helpers ==============================


def _load_known_mutation_sites(path):

    df = pd.read_csv(path, sep="\t")

    if "ID" not in df.columns or "Known Mutation Sites" not in df.columns:
        raise ValueError("Known mutation sites TSV missing required columns.")

    df["ID"] = df["ID"].astype(str)

    return df


def _build_known_dict(df):

    return {
        row["ID"]: {
            x.strip() for x in str(row["Known Mutation Sites"]).split(",") if x.strip()
        }
        for _, row in df.iterrows()
    }


# ============================== Main ==============================


def run_mutation_site_screen(**kwargs):

    start = time.time()

    base_path = kwargs.get("base_path")
    branch = kwargs.get("branch")
    known_path = kwargs.get("known_mutation_sites_tsv")

    x_min = int(kwargs.get("x_min"))  # type: ignore
    x_max = int(kwargs.get("x_max"))  # type: ignore
    y_min = int(kwargs.get("y_min"))  # type: ignore
    y_max = int(kwargs.get("y_max"))  # type: ignore
    z_min = int(kwargs.get("z_min"))  # type: ignore
    z_max = int(kwargs.get("z_max"))  # type: ignore

    score_thr_max = int(kwargs.get("score_thr_max"))  # type: ignore
    feature_thr_max = int(kwargs.get("feature_thr_max"))  # type: ignore

    n_jobs = kwargs.get("n_jobs")

    if n_jobs is None:
        n_jobs = max(cpu_count() - 2, 1)

    print(f"\nUsing {n_jobs} CPU cores\n")

    input_dir = os.path.join(base_path, "results", "pssm", "filtered", branch)  # type: ignore

    known_df = _load_known_mutation_sites(known_path)
    known_dict = _build_known_dict(known_df)

    output_dir = os.path.join(
        base_path,  # type: ignore
        "results",
        "analysis",
        "mutation_site_screen",
        branch,  # type: ignore
    )

    os.makedirs(output_dir, exist_ok=True)

    files = sorted(glob(os.path.join(input_dir, "*.tsv")))

    print("\n╔════════════════════════════════════════════════════════════╗")
    print("║ Mutation Site Recall Benchmarking (Ultra Fast Grid Search)║")
    print("╚════════════════════════════════════════════════════════════╝\n")

    # ===================== LOAD DATA (NUMPY OPTIMIZED) =====================

    protein_tables = {}

    for path in files:

        query_id = os.path.basename(path).replace(".tsv", "")

        df = pd.read_csv(path, sep="\t")

        protein_tables[query_id] = {
            "Hy": df["Hy"].to_numpy(),
            "Ch": df["Ch"].to_numpy(),
            "Po": df["Po"].to_numpy(),
            "HyCh": df["|Hy-Ch|"].to_numpy(),
            "HyPo": df["|Hy-Po|"].to_numpy(),
            "Residue": df["Residue"].astype(str).to_numpy(),
            "Position": df["Position"].astype(str).to_numpy(),
        }

    # ===================== GRID =====================

    valid_x = [x for x in range(x_min, x_max + 1) if x != 0]
    valid_y = [y for y in range(y_min, y_max + 1) if y != 0]
    valid_z = [z for z in range(z_min, z_max + 1) if z != 0]

    grid = itertools.product(
        valid_x,
        valid_y,
        valid_z,
        range(0, score_thr_max + 1),
        ["HyCh", "HyPo"],
        range(0, feature_thr_max + 1),
    )

    total_models = (
        len(valid_x)
        * len(valid_y)
        * len(valid_z)
        * (score_thr_max + 1)
        * 2
        * (feature_thr_max + 1)
    )

    print(f"Total models: {total_models:,}\n")

    summary_rows = []

    # ===================== PARALLEL EXECUTION =====================

    with Pool(
        processes=n_jobs,
        initializer=_init_worker,
        initargs=(protein_tables, known_dict),
    ) as pool:

        for result in tqdm(
            pool.imap_unordered(_evaluate_model, grid, chunksize=100),
            total=total_models,
            desc="Grid Search Progress",
            dynamic_ncols=True,
        ):
            summary_rows.append(result)

    # ===================== DATAFRAME =====================

    df = pd.DataFrame(summary_rows)

    df = df.sort_values(
        by=[
            "Micro Recall",
            "Micro Recall Polar",
            "Mean N Selected Sites",
            "Formula Complexity",
        ],
        ascending=[False, False, True, True],
    )

    df["Rank Overall"] = range(1, len(df) + 1)

    df = df.sort_values(
        by=[
            "Micro Recall Polar",
            "Micro Recall",
            "Mean N Selected Sites",
        ],
        ascending=[False, False, True],
    )

    df["Rank Polar"] = range(1, len(df) + 1)

    # ===================== SAVE FULL GRID =====================

    output_path = os.path.join(output_dir, "grid_summary.tsv")
    df.to_csv(output_path, sep="\t", index=False)

    # ===================== SPLIT RESULTS BY feature_thr =====================
    for thr, df_sub in df.groupby("feature_thr"):

        df_sub_sorted = df_sub.sort_values(
            by=[
                "Micro Recall",
                "Micro Recall Polar",
                "Mean N Selected Sites",
                "Formula Complexity",
            ],
            ascending=[False, False, True, True],
        )

        sub_output = os.path.join(output_dir, f"grid_summary_feature_thr_{thr}.tsv")

        df_sub_sorted.to_csv(sub_output, sep="\t", index=False)

    elapsed = time.time() - start

    print("\n════════════════════════════════════════════════════════════")
    print(f"Output file : {output_path}")
    print(f"Total models: {len(df)}")
    print(f"Elapsed     : {elapsed:.2f} sec")
    print("════════════════════════════════════════════════════════════\n")
