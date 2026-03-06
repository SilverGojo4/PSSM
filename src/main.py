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
    "smp_parse": {
        "title": "Extract domain-level PSSM matrices directly from CDD .smp profiles",
        "import_path": "src.preprocess.run_smp_parse.run_smp_parse",
    },
    "scorecons_reconstruct": {
        "title": "Integrate Scorecons conservation into integrated PSSM tables",
        "import_path": "src.preprocess.run_scorecons_integrate.run_scorecons_reconstruct",
    },
    "consurf_integrate": {
        "title": "Integrate ConSurf evolutionary conservation into integrated PSSM tables",
        "import_path": "src.preprocess.run_consurf_integrate.run_consurf_integrate",
    },
    "conservation_filter": {
        "title": "Mask PSSM features by conservation thresholds",
        "import_path": "src.postprocess.run_conservation_filter.run_conservation_filter",
    },
    "mutation_site_screen": {
        "title": "Screen mutation candidate sites using PSSM thresholds",
        "import_path": "src.postprocess.run_mutation_site_screen.run_mutation_site_screen",
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

    stage_info = SUPPORTED_STAGES[stage]
    module_path, func_name = stage_info["import_path"].rsplit(".", 1)
    module = importlib.import_module(module_path)
    stage_func = getattr(module, func_name)

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
        f"  - {stage:<20} {info['title']}" for stage, info in SUPPORTED_STAGES.items()
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

    parser.add_argument(
        "--branch",
        type=str,
        choices=["psiblast", "smp", "all"],
        default="all",
        help=(
            "Which PSSM branch to run (used in integrated stages).\n"
            "Options:\n"
            "  psiblast : Branch A (PSI-BLAST)\n"
            "  smp      : Branch B (.smp parsing)\n"
            "  all      : Run both branches (default)"
        ),
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

    # -------------------- Domain PSI-BLAST Args --------------------
    parser.add_argument(
        "--domain_fasta_dir",
        type=str,
        default=os.path.join(
            BASE_PATH,
            "results",
            "cdsearch",
            "domains_fasta",
            "hseq_with_gap",
        ),
        help="Directory containing extracted domain FASTA files from CD-Search.",
    )
    parser.add_argument(
        "--psiblast_output_dir",
        type=str,
        default=os.path.join(
            BASE_PATH,
            "results",
            "pssm",
            "profiles",
            "psiblast",
        ),
        help="Output directory for domain-level PSI-BLAST results.",
    )
    parser.add_argument(
        "--psiblast_blast_db",
        type=str,
        default=os.path.join(BASE_PATH, "blastdb", "cdd", "Cdd"),
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
            "Legacy root directory for PSI-BLAST results "
            "(contains pssm_profiles/ and pssm_matrices/)."
        ),
    )

    parser.add_argument(
        "--pssm_profiles_dir",
        type=str,
        help="Directory containing .pssm profiles (input for PSSM matrix extraction).",
    )

    parser.add_argument(
        "--pssm_matrix_output_dir",
        type=str,
        help="Output directory for extracted PSSM matrices (*.tsv).",
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

    parser.add_argument(
        "--pssm_matrix_dir",
        type=str,
        help="Directory containing extracted domain-level PSSM matrices (*.tsv).",
    )

    parser.add_argument(
        "--pssm_reconstruct_output_dir",
        type=str,
        help="Output directory for reconstructed full-length PSSM tables (*.tsv).",
    )

    # -------------------- SMP Parsing Args (Branch B) --------------------
    parser.add_argument(
        "--smp_cdsearch_table",
        type=str,
        help="CD-Search result table used to resolve PSSM_ID for SMP parsing.",
    )

    parser.add_argument(
        "--smp_root_dir",
        type=str,
        help="Root directory containing CDD .smp files (e.g., blastdb/cdd or blastdb/cdd/smps).",
    )

    parser.add_argument(
        "--smp_matrix_output_dir",
        type=str,
        help="Output directory for SMP-derived PSSM matrices (*.tsv).",
    )

    # -------------------- Integrated Conservation Stages Args --------------------
    parser.add_argument(
        "--scorecons_dir",
        type=str,
        help="Directory containing Scorecons output files (*.txt). (e.g., results/scorecons)",
    )

    parser.add_argument(
        "--consurf_dir",
        type=str,
        help="Directory containing ConSurf grades files (*_consurf_grades.txt). (e.g., results/consurf)",
    )

    parser.add_argument(
        "--pssm_integrated_output_dir",
        type=str,
        help="Output directory for integrated tables (Scorecons stage writes here).",
    )

    parser.add_argument(
        "--pssm_integrated_dir",
        type=str,
        help="Directory containing integrated tables (ConSurf stage reads/writes here).",
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

    parser.add_argument(
        "--pssm_filtered_output_dir",
        type=str,
        help="Output root directory for conservation-filtered TSV files (default: results/pssm/filtered).",
    )

    # -------------------- Backward Compatibility (Legacy Args) --------------------
    parser.add_argument(
        "--conservation_reconstruct_dir",
        type=str,
        help="(Legacy) reconstruct directory containing integrated TSV files.",
    )

    # -------------------- Mutation Site Screening Args (Extended Grid Version) --------------------
    parser.add_argument(
        "--x_min",
        type=int,
        help="Minimum weight for Hy in score = x*Hy + y*Ch + z*Po",
    )

    parser.add_argument(
        "--x_max",
        type=int,
        help="Maximum weight for Hy in score = x*Hy + y*Ch + z*Po",
    )

    parser.add_argument(
        "--y_min",
        type=int,
        help="Minimum weight for Ch in score = x*Hy + y*Ch + z*Po",
    )

    parser.add_argument(
        "--y_max",
        type=int,
        help="Maximum weight for Ch in score = x*Hy + y*Ch + z*Po",
    )

    parser.add_argument(
        "--z_min",
        type=int,
        help="Minimum weight for Po in score = x*Hy + y*Ch + z*Po",
    )

    parser.add_argument(
        "--z_max",
        type=int,
        help="Maximum weight for Po in score = x*Hy + y*Ch + z*Po",
    )

    parser.add_argument(
        "--score_thr_max",
        type=int,
        help="Maximum threshold for score filtering (score > threshold).",
    )

    parser.add_argument(
        "--feature_thr_max",
        type=int,
        help="Maximum threshold for secondary feature filtering (|Hy-Ch| or |Hy-Po|).",
    )

    parser.add_argument(
        "--known_mutation_sites_tsv",
        type=str,
        help="TSV file containing ground truth mutation sites (columns: ID, Known Mutation Sites).",
    )

    parser.add_argument(
        "--n_jobs",
        type=int,
        default=None,
        help=(
            "Number of CPU cores for parallel grid search.\n" "Default: cpu_count()-2"
        ),
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
