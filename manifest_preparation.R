library(tidyverse)

sdrf <- read_delim(file = "broad.mit.edu_PRAD.Genome_Wide_SNP_6.sdrf.txt", delim = "\t", col_names = T)

sdrf$sample <- as.numeric(str_sub(sdrf$`Comment [TCGA Barcode]`, 14, 15))

#the 10 below indicates a "Blood Derived Normal" sample

sdrf <- sdrf %>% filter(sample == 10)

sdrf <- sdrf %>% select(`Comment [TCGA Barcode]`, "Derived Array Data Matrix File_1")

sdrf <- sdrf %>% rename(patient = `Comment [TCGA Barcode]`, filename = "Derived Array Data Matrix File_1")

sdrf$patient = substr(sdrf$patient,1,nchar(sdrf$patient)-16)

ancestry <- read_delim(file = "ancestry_from_paper_Table_S1.txt", delim = "\t", col_names = T)

#change tumor_type and consensus_ancestry depending on cancer and population being studied

ancestry <- ancestry %>% filter(tumor_type == "PRAD" & consensus_ancestry == "eur")

ancestry <- ancestry %>% select(patient)

sdrf <- inner_join(ancestry, sdrf)

sdrf <- sdrf %>% select(filename)

manifest <- read_delim(file = "gdc_manifest.2020-12-20.txt", delim = "\t", col_names = T)

manifest <- inner_join(manifest, sdrf)

write.table(manifest, file = "gdc_manifest_prad.txt", sep = "\t", row.names = F, quote = F)