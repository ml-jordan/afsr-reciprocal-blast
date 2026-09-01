#Import python packages
import os
import pandas as pd

#Define required paths
csv_file = "260826_afsr_like_proteins.csv" #Output from 8_afsr_like_identification.R. .csv listing afsr_like proteins along with accession.
gff_folder = "/Users/jordanm/Documents/250516_rec_BLAST/gff_files" #Path to a file containing genbank GFF files.
output_file = "260826_afsr_like_locus_tag.csv" #output file.

#Read in csv_file
df = pd.read_csv(csv_file)

#Prepare results
results = []
missing_gff = []
unmatched_proteins = []

#Cache parsed GFFs to avoid re-reading
gff_cache = {}

def parse_gff(gff_path):
    """Parse GFF and return mapping of protein_id -> locus_tag"""
    mapping = {}
    
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            
            attributes = parts[8]
            
            attr_dict = {}
            for item in attributes.split(";"):
                if "=" in item:
                    key, val = item.split("=", 1)
                    attr_dict[key] = val
            
            protein_id = attr_dict.get("protein_id")
            locus_tag = attr_dict.get("locus_tag")
            
            if protein_id and locus_tag:
                mapping[protein_id] = locus_tag
    
    return mapping


#Main function loop for matching the locus tags
for _, row in df.iterrows():
    protein = row["X0"]
    accession = row["accession"]
    
    # Find matching GFF
    gff_file = None
    for file in os.listdir(gff_folder):
        if accession in file and file.endswith(".gff"):
            gff_file = os.path.join(gff_folder, file)
            break
    
    if not gff_file:
        missing_gff.append(accession)
        results.append({"protein": protein, "accession": accession, "locus_tag": None})
        continue
    
    # Parse GFF once
    if gff_file not in gff_cache:
        gff_cache[gff_file] = parse_gff(gff_file)
    
    mapping = gff_cache[gff_file]
    
    # Match protein
    locus_tag = mapping.get(protein)
    
    if not locus_tag:
        unmatched_proteins.append(protein)
    
    results.append({
        "protein": protein,
        "accession": accession,
        "locus_tag": locus_tag
    })


# Save output file
out_df = pd.DataFrame(results)
out_df.to_csv(output_file, index=False)

# Write out logs for missing gff files or unmatched proteins - had a few unmatched to manually match so worth including.
with open("missing_gff.txt", "w") as f:
    for acc in sorted(set(missing_gff)):
        f.write(acc + "\n")

with open("unmatched_proteins.txt", "w") as f:
    for prot in sorted(set(unmatched_proteins)):
        f.write(prot + "\n")

print("Done.")
print(f"Missing GFFs: {len(set(missing_gff))}")
print(f"Unmatched proteins: {len(set(unmatched_proteins))}")
