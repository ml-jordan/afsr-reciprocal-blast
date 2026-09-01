import os
import pandas as pd
import subprocess

# Defining input/output directories - change as necessary.
input_file = 'genomes.tsv'  # Input file with 'ftp_path' column which links  - usually a modifed version of the assembly summary file from NCBI FTP, the script assumes it is saved as a TSV.
query_fasta = 'query.fasta'                    # FASTA file containing protein sequence of query.
blast_output_dir = 'output_directory'                # Location of the forward BLAST output
tmp_blastdb_dir = 'tmp_blastdbs'                  # temporary location holding the blast fatabases for each of the genomes being searched.
#Generating output directory
os.makedirs(blast_output_dir, exist_ok=True)

# Load input file defining the set of genomes to search
df = pd.read_csv(input_file, sep='\t')
if 'ftp_path' not in df.columns:
    raise ValueError("Missing 'ftp_path' column in the input file")

# Get unique assembly ids from the FTP path
def extract_filename_from_ftp(ftp_path):
    return str(ftp_path).rstrip('/').split('/')[-1]

df['db_name'] = df['ftp_path'].apply(extract_filename_from_ftp)
unique_ids = df['db_name'].dropna().unique()

print(f"Found {len(unique_ids)} unique database names to process.")

# Main BLAST loop
for i, db_name in enumerate(unique_ids, 1):
    db_path = os.path.join(tmp_blastdb_dir, db_name)
    blast_output_path = os.path.join(blast_output_dir, f"{db_name}_blast.tsv")

    print(f"\n [{i}/{len(unique_ids)}] Processing: {db_name}")

    # Check if BLAST DB exists
    if not os.path.exists(f"{db_path}.pin"):  # Check for one key file
        print(f"BLAST DB not found for {db_name}, expected at {db_path}.pin. Skipping.")
        continue

    # Skip if output already exists
    if os.path.exists(blast_output_path):
        print(f"Output already exists for {db_name}, skipping.")
        continue

    # Run BLAST
    try:
        subprocess.run([
            'blastp',
            '-query', query_fasta,
            '-db', db_path,
            '-out', blast_output_path,
            '-outfmt', '6',
            '-evalue', '1e-5',
            '-num_threads', '2'
        ], check=True)

        print(f"BLAST completed: {blast_output_path}")

    except subprocess.CalledProcessError as e:
        print(f"BLAST failed for {db_name} with error: {e}")
    except Exception as e:
        print(f"Unexpected error for {db_name}: {e}")
