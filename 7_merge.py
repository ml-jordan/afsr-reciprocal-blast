import pandas as pd

# Load the reciprocal_hits file
hits = pd.read_csv("reciprocal_hits.tsv", sep="\t") #change for your reciprocal hit output from reciprocal_best_hit_filter.py

# Trim the genome column to match the accession format
hits['accession'] = hits['genome'].str.extract(r'^(GCA_\d+\.\d+)')

# Rename the rev_qseqid column to query_hits
hits = hits.rename(columns={'rev_qseqid': 'query_hits'})

# Load the filtered assembly summary tsv file
file_path = "" #path to the file.
file = pd.read_csv(file_path, sep="\t")

# Merge on the accession column
merged = file.merge(hits[['accession', 'query_hits']], on='accession', how='left')

# Save the merged file (optional)
merged.to_csv("output.tsv", sep="\t", index=False) #update to path/filename as required.

