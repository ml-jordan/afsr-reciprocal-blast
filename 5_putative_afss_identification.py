from pathlib import Path
import pandas as pd
from Bio import SeqIO

# Defer BCBio import (Had a previous isue with the environment)
try:
    from BCBio import GFF
    bcbio_available = True
except ImportError:
    bcbio_available = False

if not bcbio_available:
    raise ImportError("BCBio.GFF is not installed. Install using 'pip install bcbio-gff'.")

# Define input/output
#Input is the .tsv dataframe including a column listing the AfsR hit ('AfsR_hit' for example) and a column listing rhe genome assembly ('accession')
input_file = "PATH_TO_INPUT_DATAFRAME"
#Need a directory containing GFF file for each genome assembly.
gff_dir = Path("PATH_TO_GFF_FILES_OF_ASSEMBLY")
#The script saves the output as the input plus a column listing the AfsS hit IDs
output_file = "output.tsv"

# Read input
df = pd.read_csv(input_file, sep="\t")
df["fwd_sseqid"] = df["AfsR_hit"].fillna("").astype(str).str.strip()

# Track hits
afss_hits = []

# Helper: recursively extract CDS features
def extract_cds_features(feature):
    cds_list = []
    if feature.type == "CDS":
        cds_list.append(feature)
    if hasattr(feature, "sub_features"):
        for sub in feature.sub_features:
            cds_list.extend(extract_cds_features(sub))
    return cds_list

# Process each genome
total = len(df)
for idx, row in df.iterrows():
    accession = row["accession"]
    orthologue_id = row["fwd_sseqid"]

    print(f"[{idx+1}/{total}] Processing {accession}...", end=" ")

    if not orthologue_id or orthologue_id.lower() == "nan":
        afss_hits.append("")
        print("No orthologue.")
        continue

    matching_files = list(gff_dir.glob(f"{accession}*.gff"))
    if not matching_files:
        afss_hits.append("")
        print("No GFF found.")
        continue

    gff_path = matching_files[0]
    hit_found = ""

    with open(gff_path) as in_handle:
        for rec in GFF.parse(in_handle):
            # Collect all CDS features
            cds_features = []
            for parent in rec.features:
                cds_features.extend(extract_cds_features(parent))

            # Find orthologue feature
            orthologue_feature = None
            for feature in cds_features:
                name = feature.qualifiers.get("Name", [""])[0].strip()
                if name == orthologue_id:
                    orthologue_feature = feature
                    break

            if orthologue_feature is None:
                continue  # not in this record

            ortho_strand = orthologue_feature.strand

            # Sort CDS features based on strand direction
            if ortho_strand == 1:
                cds_sorted = sorted(cds_features, key=lambda f: f.location.start)
            else:
                cds_sorted = sorted(cds_features, key=lambda f: f.location.end, reverse=True)

            # Find orthologue index
            ortho_index = None
            for i, f in enumerate(cds_sorted):
                name = f.qualifiers.get("Name", [""])[0].strip()
                if name == orthologue_id:
                    ortho_index = i
                    break

            # Check immediate downstream CDS
            if ortho_index is not None and ortho_index + 1 < len(cds_sorted):
                neighbor = cds_sorted[ortho_index + 1]

                if neighbor.strand == ortho_strand:
                    neighbor_start = int(neighbor.location.start)
                    neighbor_end = int(neighbor.location.end)
                    length = abs(neighbor_end - neighbor_start)

                    if length < 600:
                        neighbor_name = neighbor.qualifiers.get("Name", [""])[0].strip()
                        desc = neighbor.qualifiers.get("product", ["hypothetical protein"])[0]
                        hit_found = f"{neighbor_name} ({desc}) [downstream]"

            break  # stop after first record processed

    afss_hits.append(hit_found)
    print("Hit found." if hit_found else "No hit.")

# Save results
df["afss_hit_new"] = afss_hits
df.to_csv(output_file, sep="\t", index=False)
print(f"\n Output saved to {output_file}")
