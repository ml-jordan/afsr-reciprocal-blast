#import python packages
from Bio import Entrez
import csv
import time
import os
from urllib.error import HTTPError, URLError

#enter email for access
Entrez.email = "email@email.com" #Update to your email.

#Defining the function for accessing the nucleotide sequences
def get_assembly_from_nuccore(nuccore_id, max_retries=3):
    for attempt in range(max_retries):
        try:
            handle = Entrez.elink(
                dbfrom="nuccore",
                db="assembly",
                id=nuccore_id
            )
            record = Entrez.read(handle)
            handle.close()

            assembly_ids = []
            for linkset in record:
                for linkdb in linkset.get("LinkSetDb", []):
                    if linkdb["LinkName"] == "nuccore_assembly":
                        assembly_ids.extend([link["Id"] for link in linkdb["Link"]])

            if not assembly_ids:
                return []

            handle = Entrez.esummary(
                db="assembly",
                id=",".join(assembly_ids)
            )
            summaries = Entrez.read(handle)
            handle.close()

            gcas = []
            for docsum in summaries["DocumentSummarySet"]["DocumentSummary"]:
                acc = docsum.get("AssemblyAccession")
                if acc:
                    gcas.append(acc)

            return gcas

        except (HTTPError, URLError, RuntimeError) as e:
            wait = 2 ** attempt
            print(f"Retry {attempt+1} for {nuccore_id} in {wait}s: {e}")
            time.sleep(wait)

        except Exception as e:
            print(f"Error with {nuccore_id}: {e}")
            return []

    print(f"Failed after retries: {nuccore_id}")
    return []


# Ensure it is restartable - this took a very long time to perfrom so essential.
def process_file(input_file, output_csv, checkpoint_every=100):
    # Load all IDs
    with open(input_file) as f:
        all_ids = [line.strip() for line in f if line.strip()]

    # Load already processed IDs (if restarting)
    processed = set()
    if os.path.exists(output_csv):
        with open(output_csv) as f:
            reader = csv.reader(f)
            next(reader, None)  # skip header
            for row in reader:
                if row:
                    processed.add(row[0])

    remaining_ids = [i for i in all_ids if i not in processed]

    print(f"Total IDs: {len(all_ids)}")
    print(f"Already processed: {len(processed)}")
    print(f"Remaining: {len(remaining_ids)}")

    # Open in append mode
    write_header = not os.path.exists(output_csv)

    with open(output_csv, "a", newline="") as out:
        writer = csv.writer(out)

        if write_header:
            writer.writerow(["RefSeq_NZ", "Assembly_GCA"])

        for idx, nz_id in enumerate(remaining_ids, 1):
            gcas = get_assembly_from_nuccore(nz_id)

            if gcas:
                for gca in gcas:
                    writer.writerow([nz_id, gca])
            else:
                writer.writerow([nz_id, ""])

            # Save checkpoint every N records
            if idx % checkpoint_every == 0:
                out.flush()
                os.fsync(out.fileno())
                print(f"Checkpoint saved at {idx} new records")

            # Progress update
            if idx % 50 == 0 or idx == len(remaining_ids):
                percent = (idx / len(remaining_ids)) * 100
                print(f"Processed {idx}/{len(remaining_ids)} ({percent:.1f}%)")

            time.sleep(0.4)  # safe rate

    print("Done.")


# Running the script
# input_file.txt = list of the nz nucleotide IDs for example in a single column, with no columns: NZ_JBEZJC010000004.1, NZ_JBEZJC010000006.1, NZ_JBEZJC010000007.1 etc.
# output.csv = output file. 

if __name__ == "__main__":
    process_file("input_file.txt", "output.csv")
