# src/preprocess/run_cdsearch.py
"""
CD-Search Alignment (RPS-BLAST) Stage
"""

# ============================== Standard Library Imports ==============================
import os
import subprocess
import sys
import time

# ============================== Third-Party Library Imports ==============================
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ============================== Project Root Path Setup ==============================
BASE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
if BASE_PATH not in sys.path:
    sys.path.append(BASE_PATH)


# ============================== Helper Function ==============================
def execute_cdsearch_single(seq_record, cdd_db: str, output_dir: str):
    """Run RPS-BLAST for a single protein sequence."""
    seq_id = seq_record.id.replace("|", "_")
    intermediate_dir = os.path.join(output_dir, "intermediate")
    os.makedirs(intermediate_dir, exist_ok=True)

    fasta_tmp = os.path.join(intermediate_dir, f"{seq_id}.fasta")
    hits_out = os.path.join(intermediate_dir, f"{seq_id}.tsv")

    if os.path.exists(hits_out):
        print(f"⏩ Skip {seq_id} (already completed)")
        return hits_out, "skipped", ""

    SeqIO.write(seq_record, fasta_tmp, "fasta")

    cmd = [
        "rpsblast",
        "-query",
        fasta_tmp,
        "-db",
        os.path.join(cdd_db, "Cdd"),
        "-outfmt",
        "6 qseqid sseqid pident evalue bitscore qstart qend sstart send",
        "-evalue",
        "0.001",
        "-num_threads",
        "8",
    ]

    try:
        with open(hits_out, "w") as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
        if os.path.getsize(hits_out) == 0:
            os.remove(hits_out)
            return None, "no_hit", "no domain alignment"
        return hits_out, "success", ""
    except subprocess.CalledProcessError as e:
        with open(os.path.join(output_dir, "cdsearch_error.log"), "a") as log:
            log.write(f"{seq_id}\t{e}\n")
        return None, "error", str(e)
    finally:
        if os.path.exists(fasta_tmp):
            try:
                os.remove(fasta_tmp)
            except Exception:
                pass


# ============================== Main Pipeline Function ==============================
def run_cdsearch_pipeline(**kwargs) -> None:
    """Run CD-Search (RPS-BLAST) against CDD for all sequences."""
    start_time = time.time()
    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║                 [ CD-Search Alignment Stage Started ]                ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    try:
        input_fasta = kwargs.get("cdsearch_input_fasta")
        output_dir = kwargs.get("cdsearch_output_dir")
        cdd_db = kwargs.get("cdsearch_cdd_db")

        if not os.path.exists(input_fasta):  # type: ignore
            raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")
        os.makedirs(output_dir, exist_ok=True)  # type: ignore

        print(f"📂 Input FASTA      : {input_fasta}")
        print(f"📁 Output Directory : {output_dir}")
        print(f"🧬 CDD Database     : {cdd_db}\n")

        seq_records = {rec.id: rec for rec in SeqIO.parse(input_fasta, "fasta")}
        seq_list = list(seq_records.items())
        total = len(seq_list)
        print(f"🔬 Total sequences to process: {total}\n")

        results, metadata = [], []
        success = skipped = failed = nohit = 0

        # ===================== Run RPS-BLAST for each sequence =====================
        for idx, (seq_id, seq_record) in enumerate(seq_list, start=1):
            progress = (idx / total) * 100
            print("─" * 70)
            print(f"▶️  [{idx:3d} / {total:<3d} | {progress:5.1f}% ]  {seq_id}")
            print("─" * 70)

            hits_file, status, note = execute_cdsearch_single(
                seq_record, cdd_db, output_dir  # type: ignore
            )

            if status == "success" and hits_file:
                df = pd.read_csv(
                    hits_file,
                    sep="\t",
                    header=None,
                    names=[
                        "qseqid",
                        "PSSM_ID",
                        "pident",
                        "evalue",
                        "bitscore",
                        "qstart",
                        "qend",
                        "sstart",
                        "send",
                    ],
                )
                df["query_id"] = seq_id
                results.append(df)
                success += 1
                print(f"   ✅ Alignment success for {seq_id}\n")
            elif status == "skipped":
                skipped += 1
                print(f"   ⏩ Skipped (already processed)\n")
            elif status == "no_hit":
                nohit += 1
                print(f"   ⚠️  No conserved domain detected\n")
            else:
                failed += 1
                print(f"   ❌ RPS-BLAST error encountered\n")

            metadata.append(
                {
                    "seq_id": seq_id,
                    "alignment": "yes" if status == "success" else "no",
                    "note": note,
                }
            )

        # ===================== Combine and Process Results =====================
        hits_out_path = top_hits_path = domain_fasta_dir = meta_out_path = None
        if results:
            all_hits = pd.concat(results, ignore_index=True).sort_values(
                by=["bitscore", "evalue"], ascending=[False, True]
            )
            hits_out_path = os.path.join(output_dir, "cdsearch_all_hits.tsv")  # type: ignore
            all_hits.to_csv(hits_out_path, sep="\t", index=False)

            top_hits = all_hits.groupby("query_id", as_index=False).first()
            domain_fasta_dir = os.path.join(output_dir, "domains_fasta")  # type: ignore
            os.makedirs(domain_fasta_dir, exist_ok=True)
            top_hits["domain_seq"] = None

            for _, row in top_hits.iterrows():
                qid, start, end = row["query_id"], int(row["qstart"]), int(row["qend"])
                pssm_id = row["PSSM_ID"].replace("gnl|CDD|", "")
                domain_seq = seq_records[qid].seq[start - 1 : end]
                top_hits.at[_, "domain_seq"] = str(domain_seq)  # type: ignore
                SeqIO.write(
                    SeqRecord(
                        domain_seq,
                        id=f"{qid}_{pssm_id}",
                        description=f"domain_fragment {start}-{end}",
                    ),
                    os.path.join(domain_fasta_dir, f"{qid}_{pssm_id}.fasta"),
                    "fasta",
                )

            top_hits_path = os.path.join(output_dir, "cdsearch_top_hits.tsv")  # type: ignore
            top_hits.to_csv(top_hits_path, sep="\t", index=False)

        meta_out_path = os.path.join(output_dir, "cdsearch_metadata.tsv")  # type: ignore
        pd.DataFrame(metadata).to_csv(meta_out_path, sep="\t", index=False)

        # ===================== Pretty Console Output =====================
        elapsed = time.time() - start_time
        print("╔" + "═" * 70 + "╗")
        print("║" + " " * 24 + "CD-Search Alignment Summary" + " " * 19 + "║")
        print("╚" + "═" * 70 + "╝\n")

        print("📄 Output Files")
        print("─" * 70)
        if hits_out_path:
            print(f"✅ Alignment results : {hits_out_path}")
        if top_hits_path:
            print(f"⭐ Top domain hits   : {top_hits_path}")
        if domain_fasta_dir:
            print(f"📁 Domain FASTAs     : {domain_fasta_dir}")
        print(f"🧾 Metadata summary  : {meta_out_path}\n")

        print("📊 Summary Statistics")
        print("─" * 70)
        print(f"Total sequences : {total}")
        print(f"Successful hits : {success}")
        print(f"No hits         : {nohit}")
        print(f"Skipped         : {skipped}")
        print(f"Failed          : {failed}\n")

        print(f"⏱️  Elapsed time : {elapsed:.2f} seconds")
        print("✅  CD-Search stage completed successfully!\n")

    except Exception as e:
        print(f"\n❌ Critical error in run_cdsearch_pipeline(): {e}")
        raise
