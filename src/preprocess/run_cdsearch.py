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
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
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
    xml_out = os.path.join(intermediate_dir, f"{seq_id}.xml")

    if os.path.exists(xml_out):
        print(f"⏩ Skip {seq_id} (already completed)")
        return xml_out, "skipped", ""

    SeqIO.write(seq_record, fasta_tmp, "fasta")

    cmd = [
        "rpsblast",
        "-query",
        fasta_tmp,
        "-db",
        os.path.join(cdd_db, "Cdd"),
        "-outfmt",
        "5",  # XML format (contains alignment)
        "-evalue",
        "0.001",
        "-num_threads",
        "8",
        "-out",
        xml_out,
    ]

    try:
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        if os.path.getsize(xml_out) == 0:
            os.remove(xml_out)
            return None, "no_hit", "no domain alignment"
        return xml_out, "success", ""
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


def parse_rpsblast_xml(xml_path: str, seq_id: str):
    """Parse XML to extract detailed alignment info."""
    results = []
    try:
        with open(xml_path) as f:
            blast_records = NCBIXML.parse(f)
            for record in blast_records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        qseq = hsp.query
                        hseq = hsp.sbjct
                        midline = hsp.match

                        alignment_map = "".join(
                            ["M" if c == "|" else "-" for c in midline]
                        )

                        results.append(
                            {
                                "query_id": seq_id,
                                "PSSM_ID": alignment.hit_id.replace("gnl|CDD|", ""),
                                "title": alignment.hit_def,
                                "bitscore": hsp.bits,
                                "evalue": hsp.expect,
                                "pident": (
                                    sum(c1 == c2 for c1, c2 in zip(qseq, hseq))
                                    / len(qseq)
                                )
                                * 100,
                                "qstart": hsp.query_start,
                                "qend": hsp.query_end,
                                "sstart": hsp.sbjct_start,
                                "send": hsp.sbjct_end,
                                "qseq": qseq,
                                "hseq": hseq,
                                "midline": midline,
                                "alignment_map": alignment_map,
                            }
                        )
        return results
    except Exception as e:
        print(f"⚠️ Error parsing XML for {seq_id}: {e}")
        return []


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

        all_results = []
        metadata = []
        success = skipped = failed = nohit = 0

        # ===================== Run RPS-BLAST for each sequence =====================
        for idx, (seq_id, seq_record) in enumerate(seq_list, start=1):
            progress = (idx / total) * 100
            print("─" * 70)
            print(f"▶️  [{idx:3d} / {total:<3d} | {progress:5.1f}% ]  {seq_id}")
            print("─" * 70)

            xml_file, status, note = execute_cdsearch_single(
                seq_record, cdd_db, output_dir  # type: ignore
            )

            if status == "success" and xml_file:
                parsed_hits = parse_rpsblast_xml(xml_file, seq_id)
                if parsed_hits:
                    all_results.extend(parsed_hits)
                    success += 1
                    print(f"   ✅ Alignment success for {seq_id}\n")
                else:
                    nohit += 1
                    print(f"   ⚠️  No domains parsed from XML\n")
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
        detailed_out = os.path.join(output_dir, "cdsearch_all_hits_detailed.tsv")  # type: ignore
        top_out = os.path.join(output_dir, "cdsearch_top_hits_detailed.tsv")  # type: ignore
        meta_out = os.path.join(output_dir, "cdsearch_metadata.tsv")  # type: ignore
        align_dir = os.path.join(output_dir, "alignment_blocks")  # type: ignore
        fasta_dir = os.path.join(output_dir, "domains_fasta")  # type: ignore
        query_dir = os.path.join(fasta_dir, "query")
        query_with_gap_dir = os.path.join(fasta_dir, "query_with_gap")
        hseq_with_gap_dir = os.path.join(fasta_dir, "hseq_with_gap")
        os.makedirs(align_dir, exist_ok=True)
        os.makedirs(fasta_dir, exist_ok=True)
        os.makedirs(query_dir, exist_ok=True)
        os.makedirs(query_with_gap_dir, exist_ok=True)
        os.makedirs(hseq_with_gap_dir, exist_ok=True)

        if all_results:
            df = pd.DataFrame(all_results)
            df.to_csv(detailed_out, sep="\t", index=False)

            # Get top hit for each query
            top_hits = (
                df.sort_values(by=["bitscore", "evalue"], ascending=[False, True])
                .groupby("query_id", as_index=False)
                .first()
            )

            for _, row in top_hits.iterrows():
                qid = row["query_id"]
                pssm_id = row["PSSM_ID"]
                qseq = row["qseq"]
                hseq = row["hseq"]
                midline = row["midline"]
                qstart, qend = int(row["qstart"]), int(row["qend"])

                # --- Write alignment text file ---
                align_txt = os.path.join(align_dir, f"{qid}_{pssm_id}.txt")
                with open(align_txt, "w") as f:
                    f.write(f">{qid} vs PSSM{pssm_id}\n")
                    f.write(f"Query  {qstart:4d}  {qseq}  {qend}\n")
                    f.write(f"        {midline}\n")
                    f.write(
                        f"Sbjct  {int(row['sstart']):4d}  {hseq}  {int(row['send'])}\n"
                    )

                # --- Write domain FASTA (remove gaps) ---
                domain_seq = qseq.replace("-", "")
                query_fragment = seq_records[qid].seq[qstart - 1 : qend]
                SeqIO.write(
                    SeqRecord(
                        query_fragment,
                        id=f"{qid}_{pssm_id}_query_fragment",
                        description=f"domain_fragment {qstart}-{qend}",
                    ),
                    os.path.join(query_dir, f"{qid}_{pssm_id}_query_fragment.fasta"),
                    "fasta",
                )

                # --- Write query WITH gaps (aligned qseq) ---
                SeqIO.write(
                    SeqRecord(
                        seq=Seq(qseq),
                        id=f"{qid}_{pssm_id}_query_with_gap",
                        description=f"aligned_query_with_gap {qstart}-{qend}",
                    ),
                    os.path.join(
                        query_with_gap_dir, f"{qid}_{pssm_id}_query_with_gap.fasta"
                    ),
                    "fasta",
                )

                # --- Write query WITH gaps (aligned hseq) ---
                SeqIO.write(
                    SeqRecord(
                        seq=Seq(hseq),
                        id=f"{qid}_{pssm_id}_hseq_with_gap",
                        description=f"aligned_hseq_with_gap sstart={row['sstart']} send={row['send']}",
                    ),
                    os.path.join(
                        hseq_with_gap_dir, f"{qid}_{pssm_id}_hseq_with_gap.fasta"
                    ),
                    "fasta",
                )

            top_hits.to_csv(top_out, sep="\t", index=False)

        pd.DataFrame(metadata).to_csv(meta_out, sep="\t", index=False)

        # ===================== Pretty Console Output =====================
        elapsed = time.time() - start_time
        print("╔" + "═" * 70 + "╗")
        print("║" + " " * 24 + "CD-Search Alignment Summary" + " " * 19 + "║")
        print("╚" + "═" * 70 + "╝\n")

        print("\n📄 Output Files")
        print("─" * 70)
        print(f"✅ Detailed hits    : {detailed_out}")
        print(f"⭐ Top domain hits  : {top_out}")
        print(f"📁 Alignment blocks : {align_dir}")
        print(f"📁 Domain FASTAs    : {fasta_dir}")
        print(f"🧾 Metadata summary : {meta_out}\n")

        print("📊 Summary Statistics")
        print("─" * 70)
        print(f"Total sequences : {total}")
        print(f"Successful hits : {success}")
        print(f"No hits         : {nohit}")
        print(f"Skipped         : {skipped}")
        print(f"Failed          : {failed}\n")

        print(f"⏱️  Elapsed time : {elapsed:.2f} seconds")
        print("✅  Accurate CD-Search stage completed successfully!\n")

    except Exception as e:
        print(f"\n❌ Critical error in run_cdsearch_pipeline(): {e}")
        raise
