#Script to take the identified AfsR orthologues and confirm the presence of the desired domain annotations:
# Transcriptional regulatory protein, C terminal (PF00486)
# Bacterial transcriptional activator domain (PF03704)
# P-loop containing nucleoside triphosphate hydrolases (SSF52540)
# NB-ARC (PF00931)
# Tetratricopeptide repeat domain (G3DSA:1.25.40.10)
#Additionally as the BTAD annotation in the SARP domain also contains a Tetratricopeptide repeat annotation - confirm there is a C-terminal TPR annotation 
#by looking for one starting downstream of the P-loop NTPase domain. 

##### AfsR data input ######
#Starting from a version of the assembly summary tsv containing a column 'afsr_hit' which lists the ID of identified putative AfsR proteins. 
total_afsr_data <- read.delim('', head = T) #update to path
head(total_afsr_data) #check data input
total_afsr_list <- unique(total_afsr_data$afsr_hit) #list of AfsR orthologues to screen
total_afsr_list <- subset(total_afsr_list, total_afsr_list != '0') #remove the placeholder 0 which was input in the dataframe in previous analysis.


##### Read in the InterProScan data ######
interpro_data <- read.delim('interpro.tsv', head = F) #update to interpro location.
head(interpro_data) #check it has been read in correctly.
interpro_afsr <- subset(interpro_data, interpro_data$V1 %in% total_afsr_list) #subset for AfsR hits

putative_afsr_hits <- unique(interpro_afsr$V1)
summary(putative_afsr_hits)
#Check all putative AfsR orthologues are present

###### Check presence of domains ######

#HTH 
HTH <- subset(interpro_afsr, interpro_afsr$V5 == 'PF00486')
HTH_list <- unique(HTH$V1)
# 1501 have this domain 

#BTAD
BTAD <- subset(interpro_afsr, interpro_afsr$V5 == 'PF03704')
BTAD_list <- unique(BTAD$V1)
#1551 have this domain

#NB-ARC
NB_ARC <- subset(interpro_afsr, interpro_afsr$V5 == 'PF00931')
NB_ARC_list <- unique(NB_ARC$V1)
#1542 have this domain

#P-loop
P_loop <- subset(interpro_afsr, interpro_afsr$V5 == 'SSF52540')
P_loop_list <- unique(P_loop$V1)
#1551 have this domain

#TPR 
TPR <- subset(interpro_afsr, interpro_afsr$V5 == 'G3DSA:1.25.40.10')
TPR_list <- unique(TPR$V1)
#1551 have this domain


###### Additionally want to check for TPR downstream of P-loop NTPase #######

#Crucial here you need dplyr loaded but NOT plyr.
library(dplyr)

head(interpro_afsr)

TPR_downstream <- interpro_afsr %>%
  group_by(V1) %>%
  summarise(
    A = max(V7[V5 == 'G3DSA:1.25.40.10'], na.rm = TRUE),
    B = max(V7[V5 == 'SSF52540'], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(A > B)
TPR_downstream_list <- unique(TPR_downstream$V1)
#1551 have this domain organisation.

###### Filtering for AfsR sequences with these domains - generating the 'stringent_set'
#As all have the TPR, P-loop NTPase, BTAD and the TPR downstream of P-loop, just filtering for HTH and NB-ARC.
stringent_set <- subset(interpro_afsr, interpro_afsr$V1 %in% HTH_list)
stringent_set <- subset(stringent_set, stringent_set$V1 %in% NB_ARC_list)
stringent_list <- unique(stringent_set$V1)
#1496 AfsR orthologues in this stringent set.


######Updating data frame ######
#Using the total_afsr_data dataframe.
head(total_afsr_data) #checking it's unchanged

total_afsr_data$AfsR_hit_stringent <- total_afsr_data$afsr_hit
total_afsr_data$AfsR_hit_stringent[!total_afsr_data$AfsR_hit_stringent %in% stringent_list] <- ""
head(total_afsr_data) #checking the new column is ehre correct.

stringent_AfsR_orthologues <- unique(total_afsr_data$AfsR_hit_stringent) #just to double check.
original_AfsR_hits <- unique(total_afsr_data$afsr_hit)


#####Write output to file #####
write.csv(afsr_data, 'output.csv', row.names = F) #Write the new version of the dataframe to a .csv file.
