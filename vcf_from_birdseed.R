library(tidyverse)

library(vroom)

files <- fs::dir_ls(path = "./tcga_birdseeds", glob = "*.txt")

confs <- files %>%
  map(vroom, col_names = FALSE, skip = 2, col_select = 3) %>%
  reduce(cbind)

confs <- Map(function(x) {x[x >= 0.1] <- NA;x}, confs)

confs <- as.data.frame(confs)

confs <- as_tibble(confs)

colnas <- map(confs, ~sum(is.na(.)))

colnas <- as.data.frame(colnas)

colnas <- as_tibble(colnas)

rownas <- apply(confs, MARGIN = 1, function(x) sum(is.na(x)))

rownas <- as.data.frame(rownas)

rownas <- as_tibble(rownas)

calls <- files %>%
  map(vroom, col_names = FALSE, skip = 2, col_select = 2) %>%
  reduce(cbind)

calls <- as.data.frame(calls)

#change 370 below to reflect the number of birdseed.data.txt files in the directory

calls <- calls %>% setNames(1:370)

calls <- as_tibble(calls)

#change 370 below to reflect the number of birdseed.data.txt files in the directory

colnas <- colnas %>% setNames(1:370)

calls <- calls %>% add_row(colnas)

#the 45330 below should stay the same regardless of TCGA cancer data set but check and change if needed, calculated as 0.05 * 906601 (which is the total number of SNPs on the array and should remain constant)

calls <- calls %>% select(which(calls[906601,]<=45330))

calls <- calls %>% slice(-906601)

calls <- calls %>% add_column(rownas)

#the 18 below was calculated as 0.05*370, change the 370 based on the number of birdseed.data.txt files in the directory and re-calculate the 18 accordingly

calls <- calls %>% filter(rownas<=18)

calls <- calls %>% select(-rownas)

#pick any birdseed.data.txt file from the directory you are in and substitute as file in the vroom function below

snps <- vroom(file = "./tcga_birdseeds/MESNE_p_TCGAb_401_02_03_04_05_N_GenomeWideSNP_6_H11_1486824.birdseed.data.txt", delim = "\t", col_names = FALSE, skip = 2, col_select = 1)

snps <- snps %>% setNames("snp")

snps <- snps %>% add_column(rownas)

#the 18 below was calculated as 0.05*370, change the 370 based on the number of birdseed.data.txt files in the directory and re-calculate the 18 accordingly

snps <- snps %>% filter(rownas<=18)

snps <- snps %>% select(-rownas)

# starts here - code to generate ID linker file to link this TCGA cancer data set to other data sets/types for the same TCGA cancer

files <- as_tibble(files)

#change 370 below to reflect the number of birdseed.data.txt files in the directory

order <- as_tibble(1:370)

nmissing <- as_tibble(as.numeric(t(colnas)))

header <- cbind(order, files, nmissing)

colnames(header) <- c("order","filename","nmissing")

#the 45330 below should stay the same regardless of TCGA cancer data set but check and change if needed, calculated as 0.05 * 906601 (which is the total number of SNPs on the array and should remain constant)

header <- header %>% filter(nmissing<=45330)

header <- header %>% select(-nmissing)

write.table(header, file = "tcga_prad_id_linker.txt", sep = "\t", row.names = F, quote = F)

# ends here - code to generate ID linker file to link this TCGA cancer data set to other data sets/types for the same TCGA cancer

#to get the file GenomeWideSNP_6.na35.annot.csv go to http://www.affymetrix.com/support/technical/byproduct.affx?product=genomewidesnp_6, click GenomeWideSNP_6 Annotations, CSV format, Release 35 (313 MB, 4/30/15), download, and unzip

affy <- vroom(file = "./affy/GenomeWideSNP_6.na35.annot.csv", delim = ",", skip = 18)

anno <- affy %>% select(`Probe Set ID`, `dbSNP RS ID`, Chromosome, `Physical Position`, Strand, `Allele A`, `Allele B`)

colnames(anno) <- c("snp","ID","#CHROM","POS","str","a1","a2")

anno <- inner_join(snps, anno, by = "snp")

calls <- Map(function(x) {x[x == 0] <- "0/0";x}, calls)

calls <- as.data.frame(calls)

calls <- as_tibble(calls)

calls <- Map(function(x) {x[x == 1] <- "0/1";x}, calls)

calls <- as.data.frame(calls)

calls <- as_tibble(calls)

calls <- Map(function(x) {x[x == 2] <- "1/1";x}, calls)

calls <- as.data.frame(calls)

calls <- as_tibble(calls)

annocalls <- bind_cols(anno, calls)

annocalls <- annocalls %>% filter(`#CHROM` != "---" & `#CHROM` != "MT" & `#CHROM` != "X" & `#CHROM` != "Y" & str != "---")

annocalls <- annocalls %>% mutate(a1 = case_when(str == "-" & a1 == "A" ~ "T",
                                            str == "-" & a1 == "T" ~ "A",
                                            str == "-" & a1 == "C" ~ "G",
                                            str == "-" & a1 == "G" ~ "C",
                                            str == "+" & a1 == "A" ~ "A",
                                            str == "+" & a1 == "T" ~ "T",
                                            str == "+" & a1 == "C" ~ "C",
                                            str == "+" & a1 == "G" ~ "G"
                                            ))

annocalls <- annocalls %>% mutate(a2 = case_when(str == "-" & a2 == "A" ~ "T",
                                            str == "-" & a2 == "T" ~ "A",
                                            str == "-" & a2 == "C" ~ "G",
                                            str == "-" & a2 == "G" ~ "C",
                                            str == "+" & a2 == "A" ~ "A",
                                            str == "+" & a2 == "T" ~ "T",
                                            str == "+" & a2 == "C" ~ "C",
                                            str == "+" & a2 == "G" ~ "G"
                                            ))

annocalls <- annocalls %>% select(-c(snp,str))

annocalls <- annocalls %>% rename(REF = a1)

annocalls <- annocalls %>% rename(ALT = a2)

annocalls <- annocalls %>% relocate(ID, .before = REF)

annocalls <- annocalls %>% add_column("QUAL" = ".", .after = "ALT")

annocalls <- annocalls %>% add_column("FILTER" = ".", .after = "QUAL")

annocalls <- annocalls %>% add_column("INFO" = ".", .after = "FILTER")

annocalls <- annocalls %>% add_column("FORMAT" = "GT", .after = "INFO")

annocalls <- annocalls %>% arrange(POS)

annocalls <- annocalls %>% filter(!(REF == "A" & ALT == "T"))

annocalls <- annocalls %>% filter(!(REF == "T" & ALT == "A"))

annocalls <- annocalls %>% filter(!(REF == "C" & ALT == "G"))

annocalls <- annocalls %>% filter(!(REF == "G" & ALT == "C"))

annocalls_test <- annocalls %>% distinct()

annocalls <- annocalls %>% unite("dup", `#CHROM`:POS, remove = FALSE)

annocalls <- annocalls[!(duplicated(annocalls["dup"]) | duplicated(annocalls["dup"], fromLast = TRUE)), ]

annocalls <- annocalls %>% select(-dup)

write.table(annocalls, file = "tcga_prad_eur_vcf.txt", sep = "\t", row.names = F, quote = F)

