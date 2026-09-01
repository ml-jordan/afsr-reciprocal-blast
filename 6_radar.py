import requests
import time
import os

base_url = "https://www.ebi.ac.uk/Tools/services/rest/radar" #basic url for the RADAR tool at the time of the analysis.
fasta_file = "input.faa" #protein fasta containing putative AfsS sequences - i.e. the small proteins encoded downstream of AfsR.

def run_radar(seq):
    params = {"email": "EMAIL_ADDRESS", "sequence": seq} #need to replace EMAIL_ADDRESS with your email address!
    r = requests.post(f"{base_url}/run/", data=params)
    r.raise_for_status()
    return r.text.strip()

def get_result(job_id):
    # Wait until job completes
    while True:
        status = requests.get(f"{base_url}/status/{job_id}").text.strip()
        if status == "FINISHED":
            break
        elif status == "ERROR":
            print(f"Job {job_id} failed.")
            return None
        time.sleep(5)
    result = requests.get(f"{base_url}/result/{job_id}/out").text
    return result

# Split multi-FASTA into sequences
with open(fasta_file) as f:
    seqs = f.read().strip().split(">")[1:]  # remove first empty

for seq in seqs:
    header, *seq_lines = seq.strip().splitlines()
    sequence = "".join(seq_lines)

    # Define output filename
    out_file = f"{header}_radar.txt" #saves output file for each sequence to a .txt file.

    # Skip if already processed - key for large numbers of sequences.
    if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        print(f"Skipping {header}: result file already exists.")
        continue

    # Submit job
    try:
        job_id = run_radar(f">{header}\n{sequence}")
        print(f"Submitted {header}, job ID: {job_id}")
        result = get_result(job_id)

        if result:
            with open(out_file, "w") as out:
                out.write(result)
        else:
            print(f"No result for {header} (job {job_id})")

    except Exception as e:
        print(f"Error processing {header}: {e}")
        continue
