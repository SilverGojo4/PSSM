# Protein Analysis Pipeline (PSSM)

![Pipeline Flowchart](docs/pssm_workflow.jpg)

This repository implements a **branch-aware** protein residue feature pipeline based on:

- **CDD CD-Search / RPS-BLAST** (domain detection)
- **PSSM matrix reconstruction** (domain → full-length)
- **Scorecons conservation** integration (external)
- **ConSurf evolutionary conservation** integration (external)
- **Conservation-based feature masking**

A key design update of this project is that **two independent PSSM branches are supported**:

- **Branch A (`psiblast`)**: PSI-BLAST generated `.pssm` profiles
- **Branch B (`smp`)**: direct parsing of CDD `.smp` profile matrices

Both branches are processed under the same unified folder layout and scripts.

## Project Structure

```
PSSM/
├── src/
│   ├── main.py                        # Unified pipeline entry point
│   ├── preprocess/
│   │   ├── run_cdsearch.py
│   │   ├── run_domain_psiblast.py
│   │   ├── run_pssm_features.py
│   │   ├── run_pssm_reconstruct.py
│   │   ├── run_smp_parse.py
│   │   ├── run_scorecons_integrate.py
│   │   └── run_consurf_integrate.py
│   │
│   └── postprocess/
│       ├── run_conservation_filter.py
│       └── run_mutation_site_screen.py
│
├── scripts/
│   ├── setup_cdd.sh                   # One-time setup for CDD RPS-BLAST DB
│   ├── cdsearch.sh                    # Stage 1: CD-Search alignment
│   ├── pssm.sh                        # Stage 2+3: PSI-BLAST PSSM build (Branch A)
│   ├── conservation.sh                # Stage 4: Scorecons + ConSurf integration (branch-aware)
│   ├── conservation_filter.sh         # Stage 5: Conservation-based masking (branch-aware)
│   └── mutation_site_screen.sh        # Stage 6: Mutation Site Screening & Recall Benchmarking
│
├── results/
│   ├── cdsearch/
│   ├── pssm/
│   │   ├── profiles/
│   │   ├── matrices/
│   │   ├── reconstruct/
│   │   ├── integrated/
│   │   └── filtered/
│   ├── scorecons/                     # External Scorecons outputs (*.txt)
│   └── consurf/                       # External ConSurf grades (*_consurf_grades.txt)
│
├── blastdb/
│   └── cdd/                           # Local CDD database + extracted .smp profiles
│
├── env/
│   └── pssm.yml                       # Conda environment definition
│
└── docs/
    └── pssm_workflow.jpg              # Pipeline diagram
```

## Getting Started

### Environment Setup

```bash
conda env create -f env/pssm.yml
conda activate pssm
```

Environment includes `biopython`, `pandas`, and NCBI `blast+`.

## One-Time Setup: Install CDD Database (RPS-BLAST + .smp Profiles)

```bash
bash scripts/setup_cdd.sh <PROJECT_DIR>
```

This step prepares:

- `blastdb/cdd/Cdd.*` (RPS-BLAST database)
- extracted `.smp` profiles for Branch B parsing

## Input FASTA Requirements

The pipeline requires a **protein FASTA file** as input.

Each FASTA record must follow a strict naming rule so that
wild-type and mutant sequences can be tracked consistently.

### FASTA Header Naming Rules

- **Wild-type sequence**: `{UniProt_ID}`
- **Mutant sequence**: `{UniProt_ID}_{mutation}`

Where:

- `{UniProt_ID}` is the UniProt accession ID
- `{mutation}` follows standard notation (e.g., `A53L`)

### Examples

```fasta
>D4Z2G1
MSLGAKPFGEKKFIEIKGRRMAYIDEGTGDPILFQHGNPTSSYLWRNIMPHCAGLGRLIACDLIGMGDSDKLDPSGPERYAYAEHRDYLDALWEALDLGDRVVLVVHDWGSALGFDWARRHRERVQGIAYMEAIAMPIEWADFPEQDRDLFQAFRSQAGEELVLQDNVFVEQVLPGLILRPLSEAEMAAYREPFLAAGEARRPTLSWPRQIPIAGTPADVVAIARDYAGWLSESPIPKLFINAEPGALTTGRMRDFCRTWPNQTEITVAGAHFIQEDSPDEIGAAIAAFVRRLRPA

>D4Z2G1_A53L
MSLGAKPFGEKKFIEIKGRRMAYIDEGTGDPILFQHGNPTSSYLWRNIMPHCLGLGRLIACDLIGMGDSDKLDPSGPERYAYAEHRDYLDALWEALDLGDRVVLVVHDWGSALGFDWARRHRERVQGIAYMEAIAMPIEWADFPEQDRDLFQAFRSQAGEELVLQDNVFVEQVLPGLILRPLSEAEMAAYREPFLAAGEARRPTLSWPRQIPIAGTPADVVAIARDYAGWLSESPIPKLFINAEPGALTTGRMRDFCRTWPNQTEITVAGAHFIQEDSPDEIGAAIAAFVRRLRPA
```

# Pipeline Overview (Branch-Aware)

This project is designed as a modular stage-based pipeline.

Two PSSM branches are supported:

| Branch | Name       | Source                               |
| -----: | ---------- | ------------------------------------ |
|      A | `psiblast` | PSI-BLAST generated `.pssm` profiles |
|      B | `smp`      | parsed from CDD `.smp` matrices      |

Both branches eventually converge into the same downstream structure:

```
results/pssm/reconstruct/<branch>/
results/pssm/integrated/<branch>/
results/pssm/filtered/<branch>/
```

# Stage 1 – CD-Search Alignment (RPS-BLAST)

This stage detects conserved domains for each query protein sequence.

### Run

```bash
bash scripts/cdsearch.sh <PROJECT_DIR> <INPUT_FASTA>
```

### Output

| File / Folder                                     | Description                      |
| ------------------------------------------------- | -------------------------------- |
| `results/cdsearch/cdsearch_all_hits_detailed.tsv` | all domain hits                  |
| `results/cdsearch/cdsearch_top_hits_detailed.tsv` | top domain hit per query         |
| `results/cdsearch/domains_fasta/`                 | extracted domain FASTA fragments |
| `results/cdsearch/cdsearch_metadata.tsv`          | execution summary                |

# Stage 2 + 3 – Branch A: Domain PSI-BLAST + PSSM Matrix Construction (`psiblast`)

This branch runs PSI-BLAST against the CDD database for each extracted domain fragment,
then converts `.pssm` profiles into residue-wise feature matrices.

### Run

```bash
bash scripts/pssm.sh <PROJECT_DIR> <INPUT_FASTA>
```

### Output (Branch A)

```
results/pssm/profiles/psiblast/
results/pssm/matrices/psiblast/
results/pssm/reconstruct/psiblast/
```

# Stage 2B – Branch B: Parse `.smp` Profiles (`smp`)

This branch extracts domain-level PSSM matrices directly from CDD `.smp` files.

### Run (via main.py)

```bash
python src/main.py --stage smp_parse \
  --smp_cdsearch_table results/cdsearch/cdsearch_top_hits_detailed.tsv \
  --smp_root_dir blastdb/cdd \
  --smp_matrix_output_dir results/pssm/matrices/smp
```

### Output (Branch B)

```
results/pssm/matrices/smp/
```

# Stage 3 – Full-Length PSSM Reconstruction (Both Branches)

This stage reconstructs a full-length residue-wise matrix (L × 20)
from domain-level PSSM matrices using CD-search alignment coordinates.

Output is written into branch-specific folders.

### Output

```
results/pssm/reconstruct/psiblast/
results/pssm/reconstruct/smp/
```

# Stage 4 – Conservation Integration (Scorecons + ConSurf)

This stage integrates two external conservation sources:

- **Scorecons** → domain conservation (alignment-based)
- **ConSurf** → evolutionary conservation (residue-level)

Both are integrated into reconstructed PSSM tables and stored in:

```
results/pssm/integrated/<branch>/*.tsv
```

### Required External Inputs

#### Scorecons outputs

Place Scorecons `.txt` outputs into:

```
results/scorecons/*.txt
```

Each filename must match the query FASTA header:

```
<query_id>.txt
```

#### ConSurf outputs

Place ConSurf grades outputs into:

```
results/consurf/*_consurf_grades.txt
```

The pipeline will automatically match ConSurf grades files using:

- exact `query_id`
- fallback `UniProt_ID`

## Run (Branch-aware integration pipeline)

This script automatically loops over branches.

```bash
bash scripts/conservation.sh <PROJECT_DIR> <INPUT_FASTA> [BRANCH]
```

Where:

- `BRANCH` = `psiblast | smp | all` (default: `all`)

### Output

```
results/pssm/integrated/psiblast/
results/pssm/integrated/smp/
```

# Stage 5 – Conservation-Based PSSM Feature Masking

This stage masks PSSM feature columns (sets to `NA`)
if conservation thresholds are not satisfied.

Filtering rule:

```text
EC_MIN   ≤ Evolutionary conservation ≤ EC_MAX
CONS_MIN ≤ Conservation (Scorecons)  ≤ CONS_MAX
```

This stage is strict:

- If Scorecons column is missing → file is skipped
- If Evolutionary conservation column is missing → file is skipped

## Run (Branch-aware)

```bash
bash scripts/conservation_filter.sh <PROJECT_DIR> <EC_MIN> <EC_MAX> <CONS_MIN> <CONS_MAX> [BRANCH]
```

Where:

- `BRANCH` = `psiblast | smp | all` (default: `all`)

### Output

```
results/pssm/filtered/psiblast/
results/pssm/filtered/smp/
```

# Stage 6 – Mutation Site Screening & Recall Benchmarking

This stage evaluates mutation prediction rules via **large-scale grid
search**.

Score function:

```text
score = x * Hy + y * Ch + z * Po
```

Selection rule:

```text
score > score_thr
AND
feature > feature_thr
```

Where feature is:

```text
|Hy-Ch|  (HyCh)
|Hy-Po|  (HyPo)
```

## Run

```bash
bash scripts/mutation_site_screen.sh   <PROJECT_DIR>   <X_MIN>   <X_MAX>   <Y_MIN>   <Y_MAX>   <Z_MIN>   <Z_MAX>   <SCORE_THR_MAX>   <FEATURE_THR_MAX>   <KNOWN_MUTATION_TSV>   <N_JOBS>   [BRANCH]
```

Example:

```bash
bash scripts/mutation_site_screen.sh   ~/PSSM   -10   10   -10   10   -10   10   10   10   known_mutation_sites.tsv   120   all
```

## Recall Rules

### Standard Recall

    Recall = TP / N_known_sites

If no ground-truth mutation sites exist:

    Recall = NA

### Polar Recall

Restricted to mutation sites where residue ∈ {S, T, Y, N, Q}.

Special cases:

| Condition            | Polar Recall | Status      |
| -------------------- | ------------ | ----------- |
| No GT mutation sites | NA           | NO_GT       |
| No polar GT sites    | NA           | NO_POLAR_GT |
| Polar GT exists      | Computed     | OK          |

Proteins with NA Polar Recall are excluded from Polar micro recall
calculation.

## Output Structure

```text
results/analysis/mutation_site_screen/<branch>/
├── grid_summary.tsv
├── grid_summary_feature_thr_0.tsv
├── grid_summary_feature_thr_1.tsv
├── grid_summary_feature_thr_2.tsv
└── grid_summary_feature_thr_3.tsv
```

`grid_summary.tsv` contains the complete grid search results.

Each row includes:

- scoring parameters
- recall metrics
- selected site statistics

Ranking priority:

1.  Micro Recall
2.  Polar Recall
3.  Mean Selected Sites
4.  Formula Complexity

# Full Pipeline (Recommended)

```bash
bash scripts/setup_cdd.sh <PROJECT_DIR>

bash scripts/cdsearch.sh <PROJECT_DIR> <INPUT_FASTA>

bash scripts/pssm.sh <PROJECT_DIR> <INPUT_FASTA>

# Conservation integration (Scorecons + ConSurf)
bash scripts/conservation.sh <PROJECT_DIR> <INPUT_FASTA> all

# Conservation filtering
bash scripts/conservation_filter.sh <PROJECT_DIR> 0.15 1 -1.5 10 all

# Mutation site screening (grid search recall benchmarking)
bash scripts/mutation_site_screen.sh <PROJECT_DIR> -10 10 -10 10 -10 10 10 10 data/known_mutation_sites.tsv 32 all
```

# Final Outputs

## Integrated Feature Tables

After Stage 4:

```
results/pssm/integrated/<branch>/*.tsv
```

Each file contains:

- Position
- Amino acid
- PSSM matrix columns (20 AA)
- biochemical group features (Po/Hy/Ch)
- Scorecons conservation column
- ConSurf evolutionary conservation column

## Model-Ready Filtered Tables

After Stage 5:

```
results/pssm/filtered/<branch>/*.tsv
```

These are the final recommended tables for:

- downstream statistical analysis
- ML model training
- mutation site screening

## Mutation Screening Benchmark Results

After Stage 6:

```
results/analysis/mutation_site_screen/<branch>/
├── grid_summary.tsv
├── grid_summary_feature_thr_0.tsv
├── grid_summary_feature_thr_1.tsv
├── grid_summary_feature_thr_2.tsv
└── grid_summary_feature_thr_3.tsv
```

This stage provides:

- Large-scale grid search benchmarking
- Polar-specific mutation performance analysis
- Transparent threshold sensitivity evaluation

# Notes

## Numeric Type Stability (int vs float)

Some PSSM matrices are integer-valued, while others may be float-valued.

This pipeline preserves numeric stability by restoring column dtypes
after loading/saving tables.

- integer-like columns → `Int64`
- float-like columns → `float`

This prevents unwanted `int → float` drift during integration stages.

# License

Internal research pipeline. Modify freely for lab use.
