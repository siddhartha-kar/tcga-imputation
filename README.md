
### **notes**

---

This cookbook applies to TCGA single cancer data sets though it could easily be extended to the pan-cancer data set.  The steps are presented in order and maintaining this order is usually important.

### **manifest file**

---

Visit the NIH/NCI GDC Legacy Archive <https://portal.gdc.cancer.gov/legacy-archive/search/f> and log in.  To do this you will need (1) an eRA Commons ID and password and (2) approved access to the controlled access level of TCGA via dbGAP.  The dbGAP study accession number for TCGA is currently phs000178.v11.p8.

Under the "Files" tab (on the left hand side) select Simple nucleotide variation (Data Category), Genotypes (Data Type), Genotyping array (Experimental Strategy), TXT (Data Format), Affymetrix SNP Array 6.0 (Platform), and controlled (Access Level). Under the "Cases" tab select TCGA (Cancer Program) and TCGA-PRAD (Project).  PRAD stands for prostate adenocarcinoma, which is the cancer that this example will focus on.  For other TCGA cancer abbreviations please visit: <https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations>.  Now download the manifest using the button on the right hand side.

### **sdrf file**

---

Right hand side of the same page > File Name > click any file name > Download button at the bottom of new page that opens up > download and extract the tar.gz file and go two levels into the directory and extract the mage-tab tar.gz file.  In it you will find the sdrf file, which is cancer specific.

### **ancestry file**

---

Download Table S1 from this paper: <https://pubmed.ncbi.nlm.nih.gov/32396860/> and prepare a tab-delimited text file from the first tab/sheet of the spreadsheet and save it as ``ancestry_from_paper_Table_S1.txt``. Delete the first/caption row in this tab/sheet and fix the header for columns I to M so that they are in the format Admixture_AFR, and so on.  Keep the remaining column headers as is.

### **run the script ``manifest_preparation.R``**

---

Please see additional comments in the script.  Samples used to obtain germline genotypes are "blood derived normals" and ancestry assignment is based on a consensus (majority) of five methods applied to this data set.  For details please see the paper linked above.

Input: requires the following three files to be in the working directory -

1. sdrf file: ``broad.mit.edu_PRAD.Genome_Wide_SNP_6.sdrf.txt``
2. ancestry file: ``ancestry_from_paper_Table_S1.txt``
3. manifest file: ``gdc_manifest.2020-12-20.txt``

Output: ``gdc_manifest_prad.txt`` - this is the final manifest file to be used in subsequent steps.

```{r, eval = FALSE}

#this is manifest_preparation.R

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

````
### **client and token files**

---

1. Log in to the GDC Legacy Archive <https://portal.gdc.cancer.gov/legacy-archive/search/f> and download the token file (link at the upper right hand corner under your name).  Rename the file to a simpler ``gdc-user-token.txt``

2. Download and unzip the GDC Data Transfer Tool from here: <https://gdc.cancer.gov/access-data/gdc-data-transfer-tool>.  The tool file is named ``gdc-client``

### **command line step - download data set**

---

cd into the directory where you will download the TCGA data set.  This directory must contain the following three files -

1. token file: ``gdc-user-token.txt``
2. client file: ``gdc-client``
3. manifest file: ``gdc_manifest_prad.txt``

Secure the token using chmod 600.

> chmod 600 gdc-user-token.txt

IMPORTANT: The files that you are about to download are from the controlled access level of TCGA so please use a secure research data storage facility.

> ./gdc-client download -m gdc_manifest_prad.txt -t gdc-user-token.txt

### **run the script ``vcf_from_birdseed.R``**

---

The working directory must contain a sub-directory named "tcga_birdseeds" containing all the downloaded birdseed.data.txt files.

Download ``GenomeWideSNP_6.na35.annot.csv`` from: <http://www.affymetrix.com/support/technical/byproduct.affx?product=genomewidesnp_6>.  Move this into a sub-directory named "affy" in the working directory.

Please see additional comments in the script -- there are several comments and these are important.

This script:

* Sets genotypes with birdseed confidence score >= 0.1 as missing (threshold based on Affymetrix recommendation).
* Keeps SNPs with > 95% call rate and samples with > 95% call rate.
* Aligns REF (Allele A or a1) and ALT (Allele B or a2) alleles to be on the + (forward) strand as per Affymetrix SNP Array 6.0 annotation.
* Removes ambiguous SNPs (SNPs with A/T or C/G alleles).
* Removes all copies of duplicate SNPs based on chromosome and position.
* Generates a preliminary VCF file (with the extension .txt) containing SNPs from all chromosomes, sorted by position.

Outputs from the script -

1. TCGA ID linker file (to be used in downstream analyses): ``tcga_prad_id_linker.txt``
2. preliminary VCF file with .txt extension: ``tcga_prad_eur_vcf.txt``

```{r, eval = FALSE}

#this is vcf_from_birdseed.R

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

```

### **command line step - add metadata to VCF file**

---

cd into directory containing ``tcga_prad_eur_vcf.txt`` and add the standard VCF header (metadata).

> { printf \'##fileformat=VCFv4.0\\n\'; cat tcga_prad_eur_vcf.txt; } > tcga_prad_eur.vcf

### **command line step - vcftools**

---

``tcga_prad_eur.vcf`` must be in the vcftools directory.

The VCFtools tar.gz was downloaded from here: <https://vcftools.github.io/man_latest.html> and built using instructions available here: <https://github.com/vcftools/vcftools>

1. None of the samples had inbreeding coefficient, F > +/- 0.2 or abs(F) > 0.2 so all were retained based on this filter.  Guide to interpretation: <https://sites.google.com/a/broadinstitute.org/ricopili/preimputation-qc#TOC-Technical-Details>.  This was done using visual inspection of the ``tcga_prad_eur_vcftools.het`` file output from the first command below but an R script may have to be introduced if there are samples to be filtered out.

> vcftools \-\-vcf tcga_prad_eur.vcf \-\-het \-\-out tcga_prad_eur_vcftools

2. None of the samples had relatedness_phi (kinship coefficient) > 0.0884 so all were retained based on this filter.  Guide to interpretation: <http://people.virginia.edu/~wc9c/KING/manual.html>.  This was done using visual inspection of the ``tcga_prad_eur_vcftools.relatedness2`` file output from the second command below but an R script may have to be introduced if there are samples to be filtered out.

> vcftools \-\-vcf tcga_prad_eur.vcf \-\-relatedness2 \-\-out tcga_prad_eur_vcftools

3. The third command removes SNPs with minor allele frequency < 0.005 (0.5%) and Hardy-Weinberg equilibrium exact test P < 1e-6.

> vcftools \-\-vcf tcga_prad_eur.vcf \-\-maf 0.005 \-\-hwe 1e-6 \-\-out tcga_prad_eur_vcftools \-\-recode \-\-recode-INFO-all

Output from vcftools: ``tcga_prad_eur_vcftools.recode.vcf``

### **command line step - bcftools and checkVCF**

---

* Download the checkVCF bundle from <http://qbrc.swmed.edu/zhanxw/software/checkVCF/checkVCF-20140116.tar.gz> and extract.
* Download the bcftools tar.gz from <http://www.htslib.org/download/> and extract.  The page contains installation instructions.

For the bcftools command to work, ``tcga_prad_eur_vcftools.recode.vcf`` and ``hs37d5.fa`` must be in the bcftools directory.  ``hs37d5.fa`` is part of the checkVCF bundle and must be copied over from there.

The bcftools command coupled with ``hs37d5.fa`` helps align the REF and ALT alleles with corresponding data from 1000 Genomes.

> ./bcftools norm \-\-check-ref ws -f hs37d5.fa tcga_prad_eur_vcftools.recode.vcf -o tcga_prad_eur_bcftools.checkref.vcf

checkVCF performs a series of QC/sanity checks on the VCF file output by bcftools - ``tcga_prad_eur_bcftools.checkref.vcf``.  It should yield 0's for all checks at this stage but enumerate all SNPs with ALT allele frequency > 0.5

> python checkVCF.py -r hs37d5.fa -o test tcga_prad_eur_bcftools.checkref.vcf

### **command line step - vcftools and htslib/bgzip**

---

Use vcftools to split ``tcga_prad_eur_bcftools.checkref.vcf`` by chromsome (code from: <https://gist.github.com/obenshaindw/c1afbedb0e317c1483e0>).  ``tcga_prad_eur_bcftools.checkref.vcf`` must be in the vcftools directory for this.

> seq 1 22 | xargs -n1 -P4 -I {} vcftools \-\-vcf tcga_prad_eur_bcftools.checkref.vcf \-\-chr {} \-\-recode \-\-recode-INFO-all \-\-out tcga_prad_eur.chr{}

Use bgzip to compress all the chromosome-level VCF files to vcf.gz files.  The VCF files must be in the htslib directory for bgzip to work.

* bgzip (<http://www.htslib.org/doc/bgzip.html>) is a utility in HTSlib.
* HTSlib was downloaded from: <http://www.htslib.org/download/>
* Additional (and useful!) installation instructions were found at: <https://github.com/samtools/htslib/blob/develop/INSTALL>

> ls *.vcf | xargs -n1 ./bgzip

### **imputation server**

---

Use the Michigan Imputation Server <https://imputationserver.sph.umich.edu/index.html#!> with the following options:

* Datatype: unphased
* Build: hg19
* Reference Panel: apps@1000g-phase-3-v5 (hg19)
* Population: eur
* Phasing: eagle
* Mode: imputation
* Rsq filter: 0.3

SNPs with alleles that differ from or do not match the alleles in the reference panel are filtered out.  These are listed in the ``snps-excluded.txt`` file on the server.

SNPs on the array but not in the reference panel are not used for phasing and imputation but they are not filtered out.  These SNPs are listed in the ``typed-only.txt`` file on the server.
