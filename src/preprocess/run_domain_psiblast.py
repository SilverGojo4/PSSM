# src/preprocess/run_domain_psiblast.py
"""
Domain-level PSI-BLAST Stage
"""

# ============================== Standard Library Imports ==============================
import os
import subprocess
import time
from glob import glob

# ============================== Third-Party Library Imports ==============================
import pandas as pd


# ============================== Helper Function ==============================
def run_psiblast_single(fasta_path, db, pssm_dir, threads=8):
    """
    Run PSI-BLAST on one domain FASTA.
    Generates a .pssm file and temporarily a .tsv hit file (auto-deleted on success).
    """
    base = os.path.basename(fasta_path).replace(".fasta", "")
    pssm_out = os.path.join(pssm_dir, f"{base}.pssm")
    hits_out = os.path.join(pssm_dir, f"{base}_hits.tsv")

    # Skip if already completed
    if os.path.exists(pssm_out):
        print(f"⏩ Skip {base} (already done)")
        return "skipped", base

    cmd = [
        "psiblast",
        "-query",
        fasta_path,
        "-db",
        db,
        "-num_iterations",
        "3",
        "-evalue",
        "0.001",
        "-num_threads",
        str(threads),
        "-out_ascii_pssm",
        pssm_out,
        "-out",
        hits_out,
        "-outfmt",
        "6 qseqid sseqid evalue bitscore pident length qstart qend sstart send",
    ]

    try:
        subprocess.run(cmd, check=True, stderr=subprocess.PIPE)
        if os.path.exists(hits_out):
            os.remove(hits_out)
        return "success", base

    except subprocess.CalledProcessError as e:
        log_file = os.path.join(os.path.dirname(pssm_dir), "psi_error.log")
        with open(log_file, "a") as log:
            log.write(f"[{base}] PSI-BLAST failed\n{e.stderr.decode('utf-8')}\n")
        return "error", base


# ============================== Main Stage Function ==============================
def run_domain_psiblast(**kwargs):
    """
    Pipeline Stage: Run PSI-BLAST for all domain FASTA files.
    """
    start_time = time.time()
    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║                [ Domain-level PSI-BLAST Stage Started ]              ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    domain_fasta_dir = kwargs.get("domain_fasta_dir")
    blast_db = kwargs.get("psiblast_blast_db")
    output_dir = kwargs.get("psiblast_output_dir")
    threads = kwargs.get("psiblast_threads", 8)

    if not os.path.exists(domain_fasta_dir):  # type: ignore
        raise FileNotFoundError(f"❌ Domain FASTA folder not found: {domain_fasta_dir}")
    os.makedirs(output_dir, exist_ok=True)  # type: ignore

    # Subdirectory for .pssm files
    pssm_dir = os.path.join(output_dir, "pssm_profiles")  # type: ignore
    os.makedirs(pssm_dir, exist_ok=True)

    # Find all domain FASTA files
    fasta_files = sorted(glob(os.path.join(domain_fasta_dir, "*.fasta")))  # type: ignore
    total = len(fasta_files)
    print(f"📂 Domain FASTA directory : {domain_fasta_dir}")
    print(f"📁 Output Directory       : {output_dir}")
    print(f"🧬 BLAST Database         : {blast_db}")
    print(f"🧵 Threads used           : {threads}\n")
    print(f"🔬 Total domain queries to process: {total}\n")

    summary = []
    success = skipped = failed = 0

    # ===================== Run PSI-BLAST for each domain =====================
    for idx, fasta in enumerate(fasta_files, start=1):
        base = os.path.basename(fasta).replace(".fasta", "")
        progress = (idx / total) * 100

        print("─" * 70)
        print(f"▶️  [{idx:3d} / {total:<3d} | {progress:5.1f}% ]  {base}")
        print("─" * 70)

        status, base = run_psiblast_single(fasta, blast_db, pssm_dir, threads)
        summary.append({"domain_id": base, "status": status})

        if status == "success":
            success += 1
            print(f"   ✅ PSI-BLAST completed successfully for {base}\n")
        elif status == "skipped":
            skipped += 1
            print(f"   ⏩ Skipped (already processed)\n")
        else:
            failed += 1
            print(f"   ❌ PSI-BLAST failed for {base}\n")

    # ===================== Save Metadata Summary =====================
    meta_path = os.path.join(output_dir, "psi_metadata.tsv")  # type: ignore
    pd.DataFrame(summary).to_csv(meta_path, sep="\t", index=False)

    # ===================== Summary Output =====================
    elapsed = time.time() - start_time
    print("\n" + "╔" + "═" * 70 + "╗")
    print("║" + " " * 26 + "Domain PSI-BLAST Summary" + " " * 20 + "║")
    print("╚" + "═" * 70 + "╝\n")

    print("📄 Output Files")
    print("─" * 70)
    print(f"📁 PSSM Profiles     : {pssm_dir}")
    print(f"🧾 Metadata Summary  : {meta_path}")
    print(f"🐛 Error Log (if any): {os.path.join(output_dir, 'psi_error.log')}\n")  # type: ignore

    print("📊 Summary Statistics")
    print("─" * 70)
    print(f"Total domains  : {total}")
    print(f"Successful     : {success}")
    print(f"Skipped        : {skipped}")
    print(f"Failed         : {failed}\n")

    print(f"⏱️  Elapsed time : {elapsed:.2f} seconds")
    print("✅  Domain PSI-BLAST stage completed successfully!\n")
