####### Load in the BGC data
bgc_data <- read.csv('/Users/jordanm/Documents/260423_antismashdb_afsr_like/whole_antismash/actinomycete_bgcs_edited.csv', head = T)
head(bgc_data) #check it is loaded in correctly



###### Get the coordinates for AfsR-like proteins

#a. Merge the GCF in from AfsR data

locus_tag_data <- read.csv('/Users/jordanm/Documents/2505 AfsR manuscript/260826_afsr_like/adding locus tag to afsr/260826_afsr_like_locus_tag_edited.csv', head = T)
head(locus_tag_data) #check the data
afsr_data <- read.csv('/Users/jordanm/Documents/2505 AfsR manuscript/260826_afsr_like/afsr_data.csv', head = T)
head(afsr_data) #check the data

locus_tag_data <- locus_tag_data %>%
  left_join(afsr_data %>% select(accession, gbrs_paired_asm),
            by = "accession")
head(locus_tag_data)

#Saving file
write.csv(locus_tag_data, '260827_locus_tag_data_with_gcf.csv', row.names = F)


#####Filter for AfsR-like proteins in genomes with antiSMASH coverage.

#Import data which defines the GCF assembly for each NZ molecule.
nz_to_gcf <- read.csv('/Users/jordanm/Documents/2505 AfsR manuscript/260826_afsr_like/adding afsr like coordinates/nz_to_gcf_two_complete.csv', head = T)
head(nz_to_gcf)
count <- unique(nz_to_gcf$Assembly_GCA)

#### Check all are present of the genomes are present
#Need to subset the NZ to GCF for the genomes present in the locus tag and the genomes which are not present. 

afsr_like_genomes_present <- subset(nz_to_gcf, nz_to_gcf$Assembly_GCA %in% afsr_like_gcf_list)
afsr_like_genomes_present_list <- unique(afsr_like_genomes_present$Assembly_GCA)
summary(afsr_like_genomes_present_list)
#1035 genomes are present. 


##### Subset AfsR-like data for genomes included in antismash data
locus_tag_data_as <- subset(locus_tag_data, locus_tag_data$gbrs_paired_asm %in% afsr_like_genomes_present_list)
head(locus_tag_data_as)
summary(unique(locus_tag_data_as$protein))
write.csv(locus_tag_data_as, '260827_locus_tag_data_as.csv', row.names = F)
#Performed matching using the python script - 12_adding_afsr_like_coordinates.py


#There are a few without a matching coordinates - fixed these manually using NCBI website/Snapgene to visualise.
#CQR62715.1	GCA_001013905.1	sle_32540	GCF_001013905.1
#Checked identical protein group - it is WP_099053645.1
#Checked GCF_001013905.1 online for this WP. The chromosome is NZ_LN831790.1
#Locus tag of this protein is BN2145_RS37850
#Coordinates are 3866928 - 3870587

#CRK57723.1	GCA_900070365.1	nan	GCF_900070365.1
#Checked identical protein group - WP_054048549.1
#Checked GCF_900070365.1 on NCBI website. The chromosome is NZ_LN850107.1
#Locus tag: BN1701_RS12665
#Coordinates: 2860514 - 2863513

#CRK60910.1	GCA_900070365.1	nan	GCF_900070365.1
#Checked identical protein group - WP_054053577.1
#Sequence: NZ_LN850107.1
#Locus tag are BN1701_RS27515
#Coordinates are 6047024 - 6050146

#CRK61986.1	GCA_900070365.1	nan	GCF_900070365.1
#Checked identical protein group - WP_054055208.1
#Sequence: NZ_LN850107.1
#Locus tag are BN1701_RS32585
#Coordinates are 7114939 - 7117665

#CRK62208.1	GCA_900070365.1	nan	GCF_900070365.1
#Checked identical protein group - not present. Check manually.
#On the GenBank assembly it is 7357623 - 7361093
#On the RefSeq assembly is starts one codon later - the coordinates are 7357626 - 7361093 (use this)
#NZ_LN850107.1
#Locus tag is BN1701_RS33635

#WRZ77839.1	GCA_035917275.1	OG251_40300	GCF_035917275.1
#Checked identical protein group - not present. Check manually.
#On the GenBank assembly it is 881551 - 886902 on a massive 'plasmid' - CP108509.1
# So There is an issue here where the gene is split. 
# From 881551..>883836 there is a pseudogene of labelled as SAV_2336 N-terminal domain-related protein
# From 883945..886902 there is a 'SARP-family global regulator AfsR' 
# Use the original coordinates 881551 - 886902
# NZ_CP108509.1

#fWSF26707.1	GCA_036237215.1	OG261_25495	GCF_036237215.1
#Checked identical protein group - Not present. 
# On the genabnk assembly (CP108349.1) the coordinates are 5747450 - 5750341
#Refeq is NZ_CP108349.1
#Locus tag is OG261_RS41300 - suggests pseudogene. Same start codon (5747450..>5748925)
#Use original coordinates.

#WSQ45894.1	GCA_041435455.1	OG345_24355	GCF_041435455.1
#Checked identical protein group - Not present. 
#Molecule:CP108531.1
#Coordinates: 5428246 - 5431137
#RefSeq:NZ_CP108531.1
#Seems to be a pseudogenisation suggested. 
#There is a pseudogene with the same start but a different stop. 
#5428246..>5429721 - no stop codon. I think the genbank is correct - an incomplete pseudogene would finish at the same spot as another gene (i.e. actually one gene).

#### Check the manually added molecules in the antismash database
#GCF_001013905.1 = NZ_LN831790.1 = Streptomyces leeuwenhoekii C34 = DSM 42122 = NRRL B-24963 = Correct! (33 BGCs detected.)
#GCF_900070365.1 = NZ_LN850107.1 = Alloactinosynnema sp. L-07 = Correct! (36 BGCs detected)
#GCF_035917275.1 = NZ_CP108509.1 = Streptomyces sp. NBC_01237 = Correct! (10 BGCs detected)
#GCF_036237215.1 = NZ_CP108349.1 = Streptomyces sp. NBC_01358 = Correct1 (23 BGCs detected)
#GCF_041435455.1 = NZ_CP108531.1 = Streptomyces sp. NBC_01220 = Correct! (23 BGCs detected)


###Check the locus tags are accurate
afsr_like_with_coordinates <- read.csv('/Users/jordanm/Documents/2505 AfsR manuscript/260826_afsr_like/adding afsr like coordinates/260827_output_with_coordinates_edited.csv', head = T)
#Checked manually that a subset of these were correct using NCBI website. 


#####6. Add the 20 kb flanking regions onto the BGCs. 
afsr_like_NC_set <- unique(afsr_like_with_coordinates$sequence_id)
bgc_data_afsr_like_set <- subset(bgc_data, bgc_data$NCBI.accession %in% afsr_like_NC_set)

bgc_data_afsr_like_set$Twenty_kb_flank_start <- bgc_data_afsr_like_set$From - 20000
bgc_data_afsr_like_set$Twenty_kb_flank_end <- bgc_data_afsr_like_set$To + 20000
head(bgc_data_afsr_like_set)


######Check for overlap in these. 
#AfsR-like protein complete within the flanking region. 

completely_inside <- afsr_like_with_coordinates %>%
  left_join(
    bgc_data_afsr_like_set,
    by = c("sequence_id" = "NCBI.accession"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    inside_region = start >= Twenty_kb_flank_start &
      end <= Twenty_kb_flank_end
  ) %>%
  group_by(across(all_of(names(afsr_like_with_coordinates)))) %>%
  summarise(
    inside_region = any(inside_region, na.rm = TRUE),
    .groups = "drop"
  )

head(completely_inside)

table(completely_inside$inside_region)
#FALSE  TRUE 
#1867  1333

1333/(1333+1867)
#0.4165625
#41.66%
