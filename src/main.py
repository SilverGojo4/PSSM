# src/main.py
"""
Project Main Entry Point
"""

# ============================== Standard Library Imports ==============================
import argparse
import importlib
import os
import sys

# ============================== Project Root Path Setup ==============================
BASE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))

# Fix Python Path
if BASE_PATH not in sys.path:
    sys.path.append(BASE_PATH)

# ============================== Stage Configuration ==============================
SUPPORTED_STAGES = {
    "cdsearch": {
        "title": "Run CD-Search Alignment and PSSM-ID Extraction",
        "import_path": "src.preprocess.run_cdsearch.run_cdsearch_pipeline",
    },
    "domain_psiblast": {
        "title": "Run Domain-level PSI-BLAST (using extracted domain FASTAs)",
        "import_path": "src.preprocess.run_domain_psiblast.run_domain_psiblast",
    },
    "pssm_features": {
        "title": "Compute Hy/Ch/Po features from PSSM profiles",
        "import_path": "src.preprocess.run_pssm_features.run_pssm_extract_matrix",
    },
    "pssm_reconstruct": {
        "title": "Reconstruct full-length Lx20 PSSM matrix using alignment",
        "import_path": "src.preprocess.run_pssm_reconstruct.run_pssm_reconstruct",
    },
    "scorecons_reconstruct": {
        "title": "Reconstruct full-length conservation scores using Scorecons alignment",
        "import_path": "src.preprocess.run_scorecons_integrate.run_scorecons_reconstruct",
    },
    "consurf_integrate": {
        "title": "Integrate ConSurf residue-level conservation scores",
        "import_path": "src.preprocess.run_consurf_integrate.run_consurf_integrate",
    },
    "conservation_filter": {
        "title": "Mask PSSM features by conservation thresholds",
        "import_path": "src.postprocess.run_conservation_filter.run_conservation_filter",
    },
}


# ============================== Pipeline Dispatcher ==============================
def dispatch_stage(args: argparse.Namespace) -> None:
    """
    Dispatch execution to the appropriate pipeline stage using lazy import.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing stage-specific options.
    """
    stage = args.stage.lower()
    if stage not in SUPPORTED_STAGES:
        available = ", ".join(SUPPORTED_STAGES.keys())
        raise ValueError(f"Unknown stage '{stage}'. Available stages: {available}.")

    # Dynamic import for the selected stage
    stage_info = SUPPORTED_STAGES[stage]
    module_path, func_name = stage_info["import_path"].rsplit(".", 1)
    module = importlib.import_module(module_path)
    stage_func = getattr(module, func_name)

    # Execute stage
    extra_args = vars(args)
    extra_args.pop("stage", None)
    stage_func(base_path=BASE_PATH, **extra_args)


# ============================== Main Entry ==============================
def main():
    """
    Main CLI entry point for the pipeline.
    Parses CLI arguments and routes execution to the selected pipeline stage.
    """
    available_stage_lines = [
        f"  - {stage:<15} {info['title']}" for stage, info in SUPPORTED_STAGES.items()
    ]
    available_stages_text = "\n".join(available_stage_lines)
    example_stage = list(SUPPORTED_STAGES.keys())[0]
    example_command = f"  python main.py --stage {example_stage}"

    parser = argparse.ArgumentParser(
        description=(
            "Protein Analysis Pipeline - Unified Execution Interface\n\n"
            "This script provides a centralized entry point to run different stages of the workflow.\n\n"
            "Available stages:\n"
            f"{available_stages_text}\n\n"
            "Example:\n"
            f"{example_command}\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # -------------------- General Options --------------------
    parser.add_argument(
        "--stage",
        type=str,
        choices=list(SUPPORTED_STAGES.keys()),
        required=True,
        help="Pipeline stage to run. Choose one of: "
        + ", ".join(SUPPORTED_STAGES.keys()),
    )

    # -------------------- CD-Search Args --------------------
    parser.add_argument(
        "--cdsearch_input_fasta",
        type=str,
        help="Input FASTA file used for CD-Search and PSI-BLAST (shared across stages).",
    )
    parser.add_argument(
        "--cdsearch_output_dir",
        type=str,
        help="Directory for CD-Search results.",
    )
    parser.add_argument(
        "--cdsearch_cdd_db",
        type=str,
        help="Path to local CDD BLAST database (used for both CD-Search and PSI-BLAST).",
    )
    parser.add_argument(
        "--cdsearch_cddid_tbl",
        type=str,
        help="Path to cddid_all.tbl mapping file for CDD PSSM IDs.",
    )

    # -------------------- Domain PSI-BLAST Args --------------------
    parser.add_argument(
        "--domain_fasta_dir",
        type=str,
        default=os.path.join(BASE_PATH, "results", "cdsearch_results", "domains_fasta"),
        help="Directory containing extracted domain FASTA files from CD-Search.",
    )
    parser.add_argument(
        "--psiblast_output_dir",
        type=str,
        default=os.path.join(BASE_PATH, "results", "domain_psiblast"),
        help="Output directory for domain-level PSI-BLAST results.",
    )
    parser.add_argument(
        "--psiblast_blast_db",
        type=str,
        default=os.path.join(BASE_PATH, "blastdb", "cdd"),
        help="Protein BLAST database path for PSI-BLAST.",
    )
    parser.add_argument(
        "--psiblast_threads",
        type=int,
        default=8,
        help="Number of threads for PSI-BLAST.",
    )

    # -------------------- PSSM Features Args --------------------
    parser.add_argument(
        "--pssm_root_dir",
        type=str,
        help=(
            "Root directory for domain-level PSI-BLAST results "
            "(e.g., results/domain_psiblast). "
            "This directory must contain pssm_profiles/."
        ),
    )

    # -------------------- PSSM Reconstruction Args --------------------
    parser.add_argument(
        "--pssm_fasta_path",
        type=str,
        help="Original input FASTA (full-length sequences).",
    )

    parser.add_argument(
        "--pssm_cdsearch_table",
        type=str,
        help="CD-Search result table: cdsearch_top_hits_detailed.tsv",
    )

    # -------------------- Conservation Reconstruction Args --------------------
    parser.add_argument(
        "--conservation_fasta_path",
        type=str,
        help="Original input FASTA (full-length sequences).",
    )

    parser.add_argument(
        "--conservation_cdsearch_table",
        type=str,
        help="CD-Search result table: cdsearch_top_hits_detailed.tsv",
    )

    parser.add_argument(
        "--pssm_reconstruct_dir",
        type=str,
        help="Directory containing reconstructed full-length PSSM tables.",
    )

    parser.add_argument(
        "--conservation_dir",
        type=str,
        help="Directory containing Scorecons server results (*.txt).",
    )

    parser.add_argument(
        "--conservation_reconstruct_dir",
        type=str,
        help="Directory containing reconstructed conservation tables (Scorecons output).",
    )

    # -------------------- Conservation Filtering Args --------------------
    parser.add_argument(
        "--ec_min",
        type=float,
        help="Lower bound of evolutionary conservation score",
    )

    parser.add_argument(
        "--ec_max",
        type=float,
        help="Upper bound of evolutionary conservation score",
    )

    parser.add_argument(
        "--cons_min",
        type=float,
        help="Lower bound of domain conservation score",
    )

    parser.add_argument(
        "--cons_max",
        type=float,
        help="Upper bound of domain conservation score",
    )

    # -------------------- Parse & Dispatch --------------------
    args = parser.parse_args()

    try:
        dispatch_stage(args)
        sys.exit(0)

    except Exception as e:
        print(f"\n[❌ Pipeline Error] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
