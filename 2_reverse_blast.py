import os
import subprocess
import pandas as pd
import urllib.request
import gzip
from Bio import SeqIO
import re
 
# Input/output directory/file paths - change as necessary
blast_results_dir = 'forward_blast_results_directory' #Replace with path to the directory containing the results from forward_blast.py
source_fasta_dir = 'source_fastas' #Replace with path to the directory containing the protein fastas for each proteome queried.
reference_db = os.path.join() #Path to the BLAST database constructed from the proteome of the source of the original query.
reverse_results_dir = 'output_directory' #replace with path to the desired output directory for the results from this script.
info_table_path = 'assembly_summary.tsv' #replace with the filtered assembly summary tsv containing all of the genomes included in the analysis.

# Load filtered assembly summary file with ftp_path info
info_df = pd.read_csv(info_table_path, sep='\t')

# Generate output directories if not present.
os.makedirs(reverse_results_dir, exist_ok=True)
os.makedirs(source_fasta_dir, exist_ok=True)

for file in os.listdir(blast_results_dir):
    if not file.endswith("_blast.tsv"):
        continue

    db_prefix = file.replace("_blast.tsv", "")
    blast_result_path = os.path.join(blast_results_dir, file)

    # Extract accession from filename (e.g., GCA_008064995.1)
    accession_match = re.match(r'(GCA_\d+\.\d+)', db_prefix)
    if not accession_match:
        print(f"Could not extract accession from filename: {file}")
        continue
    accession = accession_match.group(1)

    source_fasta_path = os.path.join(source_fasta_dir, f"{db_prefix}_fixed.faa")

    print(f"\n Reverse BLAST for: {blast_result_path}")

    # Download FASTA if not already available
    if not os.path.isfile(source_fasta_path):
        print(f"FASTA file not found for {db_prefix}, trying to fetch using accession {accession}...")
        row = info_df[info_df['accession'] == accession]  # Adjust if your column is named differently
        if not row.empty:
            ftp_url = row.iloc[0]['ftp_path'].rstrip('/')
            filename = ftp_url.split('/')[-1]  # e.g. GCA_008064995.1_ASM806499v1
            download_url = f"{ftp_url}/{filename}_protein.faa.gz"
            downloaded_fasta_gz = os.path.join(source_fasta_dir, f"{filename}_protein.faa.gz")

            try:
                urllib.request.urlretrieve(download_url, downloaded_fasta_gz)
                print(f"Downloaded: {download_url}")

                with gzip.open(downloaded_fasta_gz, 'rt') as f_in, open(source_fasta_path, 'w') as f_out:
                    f_out.write(f_in.read())
                os.remove(downloaded_fasta_gz)
            except Exception as e:
                print(f"Failed to download or extract FASTA for {accession}: {e}")
                continue
        else:
            print(f"No FTP path found for accession {accession} in info file.")
            continue

    # Parse BLAST results
    try:
        df = pd.read_csv(blast_result_path, sep='\t', header=None)
        subject_ids = set(df[1])  # Column 1 = sseqid
        if not subject_ids:
            print("No hits found. Skipping.")
            continue
    except Exception as e:
        print(f"Failed to read BLAST results: {e}")
        continue

    # Extract matching sequences
    reverse_query_fasta = os.path.join(reverse_results_dir, f"{db_prefix}_reverse_query.faa")

    try:
        records = list(SeqIO.parse(source_fasta_path, "fasta"))
        matched = [r for r in records if r.id in subject_ids]
        if not matched:
            print("No matching sequences found in source FASTA.")
            continue
        SeqIO.write(matched, reverse_query_fasta, "fasta")
    except Exception as e:
        print(f"Error preparing reverse query FASTA: {e}")
        continue

    # Run BLAST
    reverse_output_path = os.path.join(reverse_results_dir, f"{db_prefix}_reverse.tsv")

    try:
        subprocess.run([
            'blastp',
            '-query', reverse_query_fasta,
            '-db', reference_db,
            '-out', reverse_output_path,
            '-outfmt', '6',
            '-evalue', '1e-5',
            '-num_threads', '2'
        ], check=True)
        print(f"Reverse BLAST written: {reverse_output_path}")
    except subprocess.CalledProcessError as e:
        print(f"BLAST failed: {e}")
