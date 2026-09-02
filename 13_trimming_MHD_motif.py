#Script used (with required alterations) to trim alignments to the MHD motif in order to visualise these.
from Bio import SeqIO

# Input and output filenames
input_file = "input.fasta" #aligned fasta. 
output_file = "output.fasta" #just the MHD motif region

# Trim region: in this case I want to keep characters 1944 - 1953, requires a manual check (I use one of the sequences to base this on).
start = 1943  # Remember that Python uses 0-based indexing, so character 1944 is in Python terms character 1943.
end = 1953    # Need to remember that Python is exclusive!

# Read, trim, and write
with open(output_file, "w") as out_handle:
    for record in SeqIO.parse(input_file, "fasta"):
        trimmed_seq = record.seq[start:end]
        record.seq = trimmed_seq
        record.description += f" [trimmed {start+1}-{end}]"
        SeqIO.write(record, out_handle, "fasta")

print(f"Trimmed sequences saved to {output_file}")

