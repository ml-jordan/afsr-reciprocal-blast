import os
import re
import pandas as pd
from Bio import SeqIO

# Input and Output locations - change as required for your application
forward_dir = "" #path to results directory from forward_blast.py
reverse_dir = "" #path to results directory from reverse_blast.py
source_fasta_dir = "source_fastas" #path to source fasta directory previously defined in forward_blast.py
query_id = "EOY49372.1" #ID of the original query protein - worth checking in the fasta file used to generate the blast database used in reverse_blast.py in case of multiple names for same protein.
summary_output = "output.tsv" #Path to output TSV reporting all of the reciprocal best hits
rbh_fasta_output = "output.faa" # Path to output .faa with the sequences of all of the reciprocal best hits for further analysis.

cols = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

summary_rows = []
rbh_records = []

for filename in os.listdir(forward_dir):
    if not filename.endswith("_blast.tsv"):
        continue

    accession_full = filename.replace("_blast.tsv", "")
    match = re.match(r"(GCA_\d+\.\d+)", accession_full)
    if not match:
        print(f"Could not extract accession from: {accession_full}")
        continue
    accession = match.group(1)

    forward_path = os.path.join(forward_dir, filename)

    # Dynamically find matching reverse BLAST file
    reverse_file_match = next(
        (f for f in os.listdir(reverse_dir)
         if f.startswith(accession) and f.endswith("_reverse.tsv")),
        None
    )
    if not reverse_file_match:
        print(f"Reverse BLAST result not found for accession {accession}, skipping.")
        continue

    reverse_path = os.path.join(reverse_dir, reverse_file_match)
    source_fasta_path = os.path.join(source_fasta_dir, f"{accession_full}_fixed.faa")

    if not os.path.exists(forward_path) or not os.path.exists(reverse_path):
        print(f"Missing forward or reverse for {accession_full}, skipping.")
        continue

    try:
        fwd_df = pd.read_csv(forward_path, sep='\t', header=None, names=cols)
        rev_df = pd.read_csv(reverse_path, sep='\t', header=None, names=cols)

        if fwd_df.empty or rev_df.empty:
            continue

        # Top forward hit
        top_fwd = fwd_df.loc[fwd_df['bitscore'].idxmax()]
        top_hit_id = top_fwd['sseqid']

        # Reverse: check top reverse hit for that hit ID
        reverse_hits = rev_df[rev_df['qseqid'] == top_hit_id]
        if reverse_hits.empty:
            continue

        top_rev = reverse_hits.loc[reverse_hits['bitscore'].idxmax()]
        rev_top_id = top_rev['sseqid']

        if rev_top_id == query_id:
            combined = {
                "genome": accession_full,
                **{f"fwd_{col}": top_fwd[col] for col in cols},
                **{f"rev_{col}": top_rev[col] for col in cols}
            }
            summary_rows.append(combined)

            # Add sequence to output FASTA
            fasta_records = SeqIO.to_dict(SeqIO.parse(source_fasta_path, "fasta"))
            if top_hit_id in fasta_records:
                rbh_records.append(fasta_records[top_hit_id])
            else:
                print(f"ID {top_hit_id} not found in {source_fasta_path}")

    except Exception as e:
        print(f"Error processing {accession_full}: {e}")

# Write output files
pd.DataFrame(summary_rows).to_csv(summary_output, sep='\t', index=False)
print(f"\n RBH summary saved: {summary_output}")

SeqIO.write(rbh_records, rbh_fasta_output, "fasta")
print(f"✅ RBH sequences saved: {rbh_fasta_output}")

