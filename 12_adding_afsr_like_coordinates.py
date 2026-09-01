#Import python packages
import pandas as pd
import os
import glob

# Define required files/paths
csv_file = "260827_locus_tag_data_as.csv" #csv file listing afsr-like proteins with genbank locus tag from 11_mapping_afsr_like_locus_tag.py
gff_folder = '/Users/jordanm/Documents/260423_antismashdb_afsr_like/whole_antismash/final version/gff_files' #path to folder containing refseq version GFF files.
output_file = "260827_output_with_coordinates.csv" #output file. 

# Read in csv file.
df = pd.read_csv(csv_file)

df.columns = df.columns.str.strip()
df["locus_tag"] = df["locus_tag"].astype(str).str.strip()
df["gbrs_paired_asm"] = df["gbrs_paired_asm"].astype(str).str.strip()

# define output columns
df["start"] = None #gene start
df["end"] = None #gene end
df["sequence_id"] = None #NC_, NZ_ format molecule ID - enables matching with antismash data using NCBI accession field in that data.

# Read in and index the gff files.
gff_files = glob.glob(os.path.join(gff_folder, "*.gff")) + \
            glob.glob(os.path.join(gff_folder, "*.GFF"))

gff_map = {
    os.path.splitext(os.path.basename(f))[0]: f
    for f in gff_files
}

# Parse the function
def find_coordinates(gff_path, target_locus):
    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue

            seq_id = parts[0]
            start = parts[3]
            end = parts[4]
            attributes = parts[8]

            # parse attributes
            attrs = {}
            for item in attributes.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attrs[k] = v

            # match locus_tag or old_locus_tag
            if (
                attrs.get("locus_tag") == target_locus or
                attrs.get("old_locus_tag") == target_locus
            ):
                return start, end, seq_id

    return None, None, None

# Main loop for running the function.
for idx, row in df.iterrows():

    asm_id = row["gbrs_paired_asm"]
    locus = row["locus_tag"]

    gff_path = gff_map.get(asm_id)

    if not gff_path:
        print(f"⚠ Missing GFF: {asm_id}")
        continue

    start, end, seq_id = find_coordinates(gff_path, locus)

    if start is None:
        print(f"⚠ Not found: {locus} in {asm_id}")

    df.at[idx, "start"] = start
    df.at[idx, "end"] = end
    df.at[idx, "sequence_id"] = seq_id

# Save the output to a csv file
df.to_csv(output_file, index=False)

print(f"Saved: {output_file}")
