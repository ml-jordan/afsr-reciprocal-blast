#Define AfsR-like proteins
#1. Get the total set from the reciprocal BLAST.
#2. Remove non-included genomes.
#3. Filter for AfsR as the reverse hit
#4. Remove AfsR RBH orthologues.
#5. Aligned length greater than 500 residues
#6. Use InterProScan data to check for the presence of the domains. 
#7. Remove those in genomes which lack AfsR RBH orthologue. 
#8. Remove any duplicate proteins. 
#9. Add into the original AfsR dataframe the by-genome count.



#Reading in the start data - a combined dataframe of all of the reverse blast data
data <- read.delim('/Users/jordanm/Documents/250516_rec_BLAST/260316_afsr_like/260206_total_reverse_blast_data.tsv', head = T)
head(data) #Check the data is read in correctly.

#From the source file column, which defines the .tsv file from whcih the reverse blast data was derived, generate a new column with just the GCA_ accession.
data <- data %>%
  mutate(accession = sub("^((?:[^_]+_){1}[^_]+).*", "\\1", source_file))
head(data) #check the data is read in correctly.

#Read in the final AfsR data
afsr_data <- read.csv('/Users/jordanm/Documents/2505 AfsR manuscript/260826_afsr_like/afsr_data.csv', head = T)
head(afsr_data) #check it is read in correctly.
afsr_total_genomes <- unique(afsr_data$accession) #get a list of the genomes included in this analysis.
afsr_rbh <- unique(afsr_data$AfsR.orthologue) #get a list of the identified AfsR orthologues.

####subset data for included genomes - the data is derived from some reciprocal blast data including worse quality genomes
data <- subset(data, data$accession %in% afsr_total_genomes) 

#Filter for S. lividans AfsR as the reverse hit
data <- subset(data, data$X1 == 'EOY49372.1')
summary(data)


#Remove RBH best hit orthologues 
data <- subset(data, !(data$X0 %in% afsr_rbh))

#Filter for greater than 500 aligned length
head(data)
#X3 column is aligned length from the BLAST search
data <- subset(data, data$X3 > 500)
summary(data)


#Beginning the check for the correct domain annotations. 
#Read in the InterPro data, this is fragmented as reusing data from other analyses.
combined_ips <- read.delim('/Users/jordanm/Documents/250516_rec_BLAST/260303_afsr_like/combined_interpro_results.tsv', head = F)
head(combined_ips) #check it is read in correctly. 
additional_ips <- read.delim('/Users/jordanm/Documents/250516_rec_BLAST/260219_afsr_like/extracted_sequences.faa.tsv', head = F)
head(additional_ips) #check it is read in correctly. 
new_ips <- read.delim('/Users/jordanm/Documents/2505 AfsR manuscript/260826_afsr_like/260826_additional_sequences.faa.tsv', head = F)
head(new_ips) #check it is read in correctly. 
ips_data <- rbind(combined_ips, additional_ips, new_ips) #concatenate interpro data to one dataframe.

afsr_like_interpro_data <- subset(ips_data, ips_data$V1 %in% afsr_like_list) #Subset interpro data for that covering the potential AfsR-like proteins.
count_afsr_like_interpro <- unique(afsr_like_interpro_data$V1)
#18483 - correct number.

#Check whether any are missing.
list_afsr_like_in_ips <- subset(afsr_like_list, !(afsr_like_list %in% count_afsr_like_interpro))
summary(list_afsr_like_in_ips)
#0 are missing


######Filter for presence of domains
#C-terminal transcriptional regulatory protein domain (PF00486)
PF00486_set <- subset(afsr_like_interpro_data, afsr_like_interpro_data$V5 == 'PF00486')
PF00486_list <- unique(PF00486_set$V1)
summary(PF00486_list)
#9532

#bacterial transcriptional activator domain (BTAD; PF03704)
PF03704_set <- subset(afsr_like_interpro_data, afsr_like_interpro_data$V5 == 'PF03704')
PF03704_list <- unique(PF03704_set$V1)
summary(PF03704_list)
#12338

#P-loop NTPase (SSF52540)
SSF52540_set <- subset(afsr_like_interpro_data, afsr_like_interpro_data$V5 == 'SSF52540')
SSF52540_list <- unique(SSF52540_set$V1)
summary(SSF52540_list)
#18461

#NB-ARC domain (PF00931)
PF00931_set <- subset(afsr_like_interpro_data, afsr_like_interpro_data$V5 == 'PF00931')
PF00931_list <- unique(PF00931_set$V1)
summary(PF00931_list)
#11201

#tetratricopeptide repeats (TPRs; G3DSA:1.25.40.10)
TPR_set <- subset(afsr_like_interpro_data, afsr_like_interpro_data$V5 == 'G3DSA:1.25.40.10')
TPR_list <- unique(TPR_set$V1)
summary(TPR_list)
#18469

##### Filter for these domains.
data <- subset(data, data$X0 %in% PF00486_list)
data <- subset(data, data$X0 %in% PF03704_list)
data <- subset(data, data$X0 %in% SSF52540_list)
data <- subset(data, data$X0 %in% PF00931_list)
data <- subset(data, data$X0 %in% TPR_list)
afsr_like_list <- unique(data$X0)
#5692 proteins.

#Need to check that there is a TPR domain downstream of the ATPase domain - this is because the BTAD annotation co-occurs with a TPR annotation.
#Crucial here you need dplyr loaded but NOT plyr.
library(dplyr)

afsr_like_interpro_data <- subset(afsr_like_interpro_data, afsr_like_interpro_data$V1 %in% afsr_like_list)
summary(unique(afsr_like_interpro_data$V1)) #Checked 5692 are presnt


downstream_check <- afsr_like_interpro_data %>%
  group_by(V1) %>%
  summarise(
    A = max(V7[V5 == 'G3DSA:1.25.40.10'], na.rm = TRUE),
    B = max(V7[V5 == 'SSF52540'], na.rm = TRUE),
    .groups = "drop"
  )

head(downstream_check)
#A is the maximum start for TPR. 
#B is the maximum start for the P-loop 
#Need A is greater than B. 

downstream_check_subset <- subset(downstream_check, downstream_check$A > downstream_check$B)
#4269

downstream_check$check <- downstream_check$A - downstream_check$B
head(downstream_check)
downstream_check_approach_two <- subset(downstream_check, downstream_check$check > 0)
#Also 4269 

downstream_check_subset_list <- unique(downstream_check_subset$V1)
downstream_check_approach_two_check <- subset(downstream_check_approach_two, !(downstream_check_approach_two$V1 %in% downstream_check_subset_list))
#Same. 

##### Now need to subset data
data <- subset(data, data$X0 %in% downstream_check_subset_list)
#Had 4275 rows, why?
check_for_duplicates <- as.data.frame(table(data$X0))
View(check_for_duplicates)
#There are enough duplicates to cause this. 
#4269 unique proteins. 



#7. Remove those in genomes which lack AfsR RBH orthologue. 
genomes_with_afsr <- subset(afsr_data, afsr_data$AfsR.orthologue != '')
genome_accessions_with_afsr <- unique(genomes_with_afsr$accession)

data <- subset(data, data$accession %in% genome_accessions_with_afsr)
head(data)
count_afsr_like <- unique(data$X0) 
#3853 proteins


######### Removing any duplicates
#Preserving the data for the highest bitscore.
data_unique <- data %>%
  group_by(X0) %>%
  slice_max(order_by = X11, n = 1, with_ties = FALSE) %>%
  ungroup()
summary(data_unique)

#Write the AfsR-like proteins to a dataframe.
write.csv(data_unique, '260826_afsr_like_proteins.csv', row.names = F)

####Generating summary data. 
afsr_like_count_data <- as.data.frame(table(data_unique$accession)) #Count per genome
head(afsr_like_count_data)
afsr_like_genome_list <- unique(afsr_like_count_data$Var1) #genomes with an AfsR-like protein
head(afsr_like_genome_list)

genomes_without_afsr_like <- subset(afsr_data, afsr_data$AfsR.orthologue != '') #Getting genomes with AfsR RBH
genomes_with_afsr_like <- subset(genomes_without_afsr_like, genomes_without_afsr_like$accession %in% afsr_like_genome_list) #Genomes with AfsR RBH and AfsR-like
genomes_without_afsr_like <- subset(genomes_without_afsr_like, !(genomes_without_afsr_like$accession %in% afsr_like_genome_list)) Genomes with AfsR RBH, no AfsR-like so should be 0.

summary(genomes_with_afsr_like)
summary(genomes_without_afsr_like)

#232 + 1264 = 1496 #check the number is correct. 

#Generating dataframe with 0 as count for those with no AfsR-like proteins. 
Var1 <- unique(genomes_without_afsr_like$accession)
head(Var1)
Freq <- c(0)
genomes_without_afsr_like_df <- data.frame(Var1, Freq)
head(genomes_without_afsr_like_df)

#Assembling count dataframe for genomes with and without AfsR-like proteins.
total_afsr_like_df <- rbind(afsr_like_count_data, genomes_without_afsr_like_df)
View(total_afsr_like_df)
sum(total_afsr_like_df$Freq)

total_afsr_like_df
head(total_afsr_like_df)

#Merge into AfsR data for AfsR-like column
colnames(total_afsr_like_df)[c(1, 2)] <- c("accession", "afsr_like_count")
head(total_afsr_like_df)

new_afsr_data <- afsr_data %>%
  left_join(total_afsr_like_df %>% select(accession, new_afsr_like_count),
            by = "accession")

View(new_afsr_data)

#Filter specifically for Streptomyces. 
streptomyces <- subset(new_afsr_data, new_afsr_data$genus == 'Streptomyces')
streptomyces_genome_list <- unique(streptomyces$accession)
sum(streptomyces$new_afsr_like_count)

#2670 in total

library(ggplot2)
ggplot(streptomyces, aes(x = afsr_like_count)) +
  geom_histogram(binwidth = 1) + 
  theme_bw() + 
  ylab('Count') + 
  xlab('Number of AfsR-like proteins per genome')


#Plot for total AfsR-like
ggplot(new_afsr_data, aes(x = afsr_like_count)) +
  geom_histogram(binwidth = 1) + 
  theme_bw() + 
  ylab('Count') + 
  xlab('Number of AfsR-like proteins per genome')







