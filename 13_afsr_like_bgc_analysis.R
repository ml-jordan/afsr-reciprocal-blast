##### AfsR-like data input #####
# afsr_like_of_interest = AfsR-like data with the coordiantes previously generated in 12_adding_afsr_like_coordinates.py
afsr_like_of_interest <- read.csv('/Users/jordanm/Documents/250516_rec_BLAST/260316_afsr_like/260428_mapping_afsr_like_locus/output_with_coordinates.csv', head = T)
head(afsr_like_of_interest) #check data is read in correctly.
afsr_like_genome_list <- unique(afsr_like_of_interest$sequence_id)


##### Read in antiSMASH BGC data #####
bgc_data <- read.csv('/Users/jordanm/Documents/260423_antismashdb_afsr_like/whole_antismash/actinomycete_bgcs_edited.csv', head = T)


##### Adding 20 kb flanking region to antiSMASH annotations #####
#We noticed that often there are AfsR-like proteins encoded immediately adjacent to known BGCs. 
#This suggests to us that the boundaries are often missing these. 
#An example of this is in GCA_000024025.1 where ACU73113.1 is encoded immediately adjacent to a BGC (#Region 15 - lanthipeptide-class-iii).
#To account for this, have added a 20 kb flanking region to the antiSMASH predictions for the cluster boundary. 
#20 kb was chosen as a recent paper Augustijn et al 2025 identified non-BGC associated SARP regulators if they were more than 20 kb away. 

head(bgc_data_relevant)
bgc_data_relevant <- subset(bgc_data, bgc_data$NCBI.accession %in% afsr_like_genome_list)
bgc_data_relevant$from_20_kb <- bgc_data_relevant$From - 20000
bgc_data_relevant$to_20_kb <- bgc_data_relevant$To + 20000

afsr_like_of_interest_with_BGC_data <- afsr_like_of_interest

afsr_like_of_interest_with_BGC_data$in_region <- FALSE
match_count <- 0

for (i in seq_len(nrow(bgc_data_relevant))) {
  
  new_matches <- (
    afsr_like_of_interest_two$sequence_id == bgc_data_relevant$NCBI.accession[i] &
      afsr_like_of_interest_two$start >= bgc_data_relevant$from_20_kb[i] &
      afsr_like_of_interest_two$end   <= bgc_data_relevant$to_20_kb[i]
  )
  
  match_count <- match_count + sum(new_matches)
  
  afsr_like_of_interest_with_BGC_data$in_region <- afsr_like_of_interest_with_BGC_data$in_region | new_matches
}

write.csv(afsr_like_of_interest_with_BGC_data, 'afsr_like_of_interest_with_BGC_data.csv', row.names = F)


##### Count the number in each class #####
table(afsr_like_of_interest_with_BGC_data$in_region) #Summarise the number in each category


##### Filtering dataset to generate sets which are either inside the BGC+Flanking_region or outside #####
inside_region <- subset(afsr_like_of_interest_with_BGC_data, afsr_like_of_interest_with_BGC_data$in_region == 'TRUE')
inside_region_list <- unique(inside_region$protein)

outside_region <- subset(afsr_like_of_interest_with_BGC_data, afsr_like_of_interest_with_BGC_data$in_region == 'FALSE')
outside_region_list <- unique(outside_region$protein)

write.csv(afsr_like_of_interest_two, '260514_afsr_like_BGC_colocalisation.csv', row.names = F)


##### Write to FASTA file #####
library(Biostrings)

# read in fasta containing all AfsR-like proteins
total_afsr_like <- readAAStringSet("/Users/jordanm/Documents/250516_rec_BLAST/260316_afsr_like/augustijn_mhd/260324_extracted_sequences_augustijn.fasta")
head(total_afsr_like)

# Need to generate vectors to denote which sequences to keep based on the headers. 

total_afsr_like_headers <- names(total_afsr_like) # fasta headers

inside_region_afsr_like <- total_afsr_like_headers[total_afsr_like_headers %in% inside_region_list]
outside_region_afsr_like <- total_afsr_like_headers[total_afsr_like_headers %in% outside_region_list]

#Write FASTA files. 
writeXStringSet(inside_region_afsr_like, 
                "inside_region_sequences.fasta")
writeXStringSet(outside_region_afsr_like, 
                "outside_region_sequences.fasta")
