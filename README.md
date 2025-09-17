# 🧬 Protein Analysis Pipeline

A modular, CLI-based bioinformatics pipeline for:

- Fetching **UniProt sequences** (wild-type and optionally mutated),
- Running **PSI-BLAST** against NCBI’s **Conserved Domain Database (CDD)** to generate **PSSM** and **domain hit** annotations.

## 📁 Project Structure

```
AMPscope/
├── src/
│   ├── main.py                  # Unified pipeline entry point
│   └── preprocess/
│       ├── fetch_uniprot.py     # UniProt sequence & mutation handler
│       └── run_pssm.py          # PSI-BLAST + PSSM generator
├── scripts/
│   ├── fetch_uniprot.sh         # Shell wrapper for UniProt stage
│   ├── run_pssm.sh              # Shell wrapper for PSI-BLAST stage
│   └── setup_cdd.sh             # One-time setup for CDD RPS-BLAST DB
├── data/
│   ├── raw/                     # Input protein lists
│   ├── processed/               # FASTA files and UniProt metadata
├── results/                     # PSSM & CDD domain hits
└── blastdb/
    └── cdd/                     # RPS-BLAST CDD database
```

## 🚀 Getting Started

### 🧱 Environment Setup

This project uses a conda environment named `pssm`. You can recreate it using the provided YAML file:

```bash
conda env create -f env/pssm.yml
```

Activate it:

```bash
conda activate pssm
```

The environment includes:

- `biopython`
- `pandas`
- `blast+` (for `psiblast`, `makeprofiledb`)

Make sure to update the environment file if dependencies change.

## 🧬 Stage 1: Fetch UniProt Sequences + Apply Mutations

### 🔧 Input Format

Supported formats: `.xlsx`, `.csv`, `.tsv`

| UniProt_ID | Mutation |
| ---------- | -------- |
| P69905     | A17V     |
| P68871     | K11E     |

- `Mutation` is optional.
- Mutation format follows standard notation (e.g., `A17V` = Ala at pos 17 mutated to Val).

### ▶️ Run

```bash
bash scripts/fetch_uniprot.sh <PROJECT_DIR>
```

### 📝 Output

- `proteins_wt.fasta`: Wild-type FASTA sequences
- `proteins_mut.fasta`: Mutated sequences (if any)
- `uniprot_metadata.csv`: Summary table (WT, mutation, mutation status)

---

## 🧪 Stage 2: Run PSI-BLAST to Generate PSSM

This step runs **PSI-BLAST** for each input sequence (wild-type or mutated), and collects:

- `.pssm` profiles
- Domain hit annotations from **CDD**

### ▶️ Run

```bash
bash scripts/run_pssm.sh <PROJECT_DIR>
```

### 📝 Output

- `results/pssm_raw/`: All `.pssm` and `.tsv` per sequence
- `results/cdd_hits.tsv`: Aggregated domain hits

## 🧰 One-Time: Setup CDD Database for RPS-BLAST

### 💾 Prerequisites

- Download the `.tar` archive of the **CDD profiles** (e.g. from NCBI FTP)
- Place it as: `blastdb/cdd.tar`

### ▶️ Run

```bash
bash scripts/setup_cdd.sh <PROJECT_DIR>
```

This will:

- Extract `.smp` domain profiles
- Create `Cdd.pn` and build RPS-BLAST DB via `makeprofiledb`

## 🧠 Advanced Usage

Run each stage individually via CLI:

```bash
# Run uniprot fetch stage
python src/main.py \
  --stage uniprot_fetch \
  --uniprot_input_path data/raw/test_protein_list.xlsx \
  --uniprot_output_wt_fasta data/processed/proteins_wt.fasta \
  --uniprot_output_mut_fasta data/processed/proteins_mut.fasta \
  --uniprot_output_csv data/processed/uniprot_metadata.csv

# Run PSI-BLAST
python src/main.py \
  --stage pssm_search \
  --pssm_input_fasta data/processed/proteins_mut.fasta \
  --pssm_cdd_db blastdb/cdd/Cdd \
  --pssm_output_dir results/pssm_raw \
  --pssm_hits results/cdd_hits.tsv
```

## 📊 Output Example

Example line from `cdd_hits.tsv`:

```
qseqid      sseqid      evalue   bitscore   pident   length  qstart  qend  sstart  send
P69905|A17V cd08953     1e-05     57.1      42.0      90      5       88     3      86
```

## ❗ Troubleshooting

- ❌ **Missing UniProt ID**: Make sure `UniProt_ID` column exists (case-insensitive).
- ❌ **Mutation failed**: Check format (e.g., `K23R`, not `23K→R`)
- ❌ **CDD not found**: Ensure `.smp` files exist and `makeprofiledb` is installed.
