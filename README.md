# AfsR reciprocal BLAST introduction
This repository contains the code relating to reciprocal BLAST analysis of AfsR, AfsS and WblH in the phylum Actinomycetota.
This analysis was published as part of: Jordan et al., 2026 'Autoactivation of an NLR-type transcriptional regulator unlocks global control of antibiotic production in _Streptomyces_' (preprint available here: ____). 

The included scripts are intended for use for taking a protein sequence (in this case the sequence of the AfsR orthologue in _Streptomyces lividans_ - SLI_4664 or the WblH orthologue in _S. lividans_ - SLI_7059) and performing a reciprocal BLAST analysis to identify potential reciprocal best hit (RBH) orthologues across Actinomycete genomes. The time taken per protein differs based on sequence length and how well conserved it is, but it is possible to do this at quite high throughput. I have used a derivative of these scripts to perform reciprocal BLAST analysis on a wide array of _Streptomyces_ proteins - it is probably not peak efficiency but it has been a very reliable approach for our projects!

Downstream of the core reciprocal blast hit analysis additional analysis was sometimes applied. For WblH, the reciprocal best hit was taken as the 'true' orthologue. For AfsR, InterProScan was used to confirm no extraneous domains and confirm the presence of the expected annotations. This analysis was further expanded to identify additional AfsR-like proteins - reciprocal hits that are not the reciprocal best hit, located in genomes with an identified AfsR orthologue, containing all required domain annotations. Furthermore, putative AfsS orthologues were identified by screening for small proteins encoded immediately downstream of AfsR which contain repeats within their sequence according to RADAR (https://www.ebi.ac.uk/jdispatcher/pfa/radar).  


# Software versions:
The analysis was performed using the following versions:
- Python: 3.13.2
- R: 4.4.1
- Rstudio: 2024.09.0+375
- InterProScan: 5.52-86.0


# List of scripts:
- 1_forward_blast.py: Takes protein of interest and performs a forward BLAST search against a directory of fasta files containing the proteomes of the desired genomes.
- 2_reverse_blast.py: Takes the output from forward_blast.py, performs a BLAST search against a BLAST database corresponding to the proteome from which the protein of interest is derived.
- 3_reciprocal_best_hit_filter.py: Applies a reciprocal best hit filter based on bitscore, identifies proteins which were the top hit in a particular proteome for the query, and also have the query as the top hit in the original proteome.
- 4_afsr_domain_check.R: R Script which was used to filter for AfsR orthologues based on InterProScan data.
- 5_putative_afss_identification.py: Script to use GFF files to identify a small protein encoded immediately downstream of a protein of interest - used in our case to identify candidate AfsS proteins based on their size and position encoded next to AfsR.
- 6_radar.py: Script to run the RADAR webtool on a protein fasta of interest - used in our case to identify repeats within AfsS - basic automation of submitting the sequences.
- 7_merge.py: Accessory script to add an additional column to the original assembly summary dataframe listing identified reciprocal best hits.
- 8_afsr_like_identification.R: R script used to identify AfsR-like proteins in proteomes containing an AfsR orthologue.
- 9_matching_nz_to_GCF.py: Matches the NZ_ format nucleotide accessions and matches them to a RefSeq assembly. 
- 10_AfsR_like_BGC_localisation.R: R script used to process the AfsR-like proteins to look for whether they are encoded within, or close to biosynthetic gene clusters (BGCs) predicted in the antiSMASH database. Requires some of the subsequent scripts at various points.
- 11_mapping_afsr_like_locus_tag.py: Takes the AfsR-like proteins and adds the locus_tag onto the proteins from the GenBank assemblies.
- 12_adding_afsr_like_coordinates.py: Takes the AfsR-like protein, adds the coordinates of the genes to enable matching with the BGCs.
- 13_trimming_MHD_motif.py: Python script to trim an aligned fasta for the MHD motif. 


# Required files:
- Protein FASTA of sequence of interest.
- Assembly summary file from NCBI FTP filtered for genomes of interest in the .TSV format (in our case the file was downloaded from ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt and then filtered for complete Actinomycete genomes assembled to fewer than 10 contigs with _S. coelicolor_ A3(2) and _S. lividans_ 1326 added in as an internal control. 
- source_fasta folder containing protein FASTAs for the proteome of each organism of interest - can be downloaded from the NCBI FTP.
- BLAST database constructed from the proteome which the initial protein sequence comes from (in this case _S. lividans_ 1326).
- GFF files for each genome of interest (genbank assembly) if performing the AfsS identification extension.
- If using the afsr_domain_check.R or afsr_like_identification.R scripts, InterProScan output for your sequences of interest.
- For BGC analysis, you also need BGC data downloaded from antiSMASH database, a folder of RefSeq format GFF files. 


# Pipeline instructions:
Step 1. Forward BLAST using forward_blast.py
1. Using the asssembly summary file (ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt), filter for genomes of interest - in this case complete Actinomycete genomes assembled to less than 10 contigs. Download these genomes to source_fastas file.
2. Save the protein sequences from each query to a .fasta file.
3. Generate a BLAST database for each of the query fastas.
4. Correct the locations of relevent files/directories in the forward_blast.py script.
5. Run forward_blast.py.

Step 2. Reverse BLAST using reverse_blast.py
1. Use make-blast-db with the -parse_seq_ids flag to generate a blast database from the protein FASTA corresponding to the proteome of the source of your query.
2. Update the file/directory locations in the reverse_blast.py script.
3. Run reverse_blast.py. I would recommend using the caffeinate option (caffeinate python3 reverse_blast.py) if working on Mac to stop the computer from going to sleep whilst the script is running.

Step 3. Reciprocal best hit filter using reciprocal_best_hit_filter.py
1. Check the name of your original query in the source proteome - I have found in organisms where there are multiple annotations (Streptomyces venezuelae for example) this is an easy place to trip up!
2. Check that reciprocal_best_hit_filter.py is up to date with the correct locations of files/directories generated in Step 1 and Step 2.
3. Run reciprocal_best_hit_filter.py.
4. Can merge the hits into the dataframe using merge.py or perform additional downstream analyses as appropriate.

After this additional analysis can be performed - the R script afsr_domain_check.R was used to confirm the presence of specific domains within AfsR orthologues, afsr_like_identification.R to identify the AfsR-like proteins, putative_afss_identification.py to identify small proteins encoded immediately downstream of AfsR (i.e. potential AfsS orthologues) and radar.py to run RADAR on these sequences.
