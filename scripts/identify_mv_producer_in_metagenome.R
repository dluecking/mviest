# author: dlueckin
# Tue Jan 11 09:54:50 2022 

  
# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
suppressMessages(library(dplyr))

# # working directory -------------------------------------------------------
# setwd("~/Documents/promotion/misc/test_coverage_identification")

# read arguments ----------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 2){
    print("Usage:\nRscript identify_mv_producer_in_metagenome.R <sample> <output_path>")
    quit()
}
sample = args[1]
output_path = args[2]

r1_cov_path = paste0("results/", sample,"/mappings/meta_1/r1_vs_meta.sorted.bam.coverage")
r2_cov_path = paste0("results/", sample,"/mappings/meta_2/r2_vs_meta.sorted.bam.coverage")

# load data ---------------------------------------------------------------

if(file.info(r1_cov_path)$size == 0L | file.info(r2_cov_path)$size == 0L){
    print("r1 or r2 empty: aborting...")
    write("No mv producers identified.", output_path)
    quit()
}

r1_cov <- fread(r1_cov_path)
r2_cov <- fread(r2_cov_path)

# merge data
r1_cov$pos <- paste0(r1_cov$V1, "_", r1_cov$V2)
r2_cov$pos <- paste0(r2_cov$V1, "_", r2_cov$V2)

cov <- data.table(pos = unique(c(r1_cov$pos, r2_cov$pos)),
                  cov = 0)

cov$cov <- r1_cov$V3[match(cov$pos, r1_cov$pos)] + r2_cov$V3[match(cov$pos, r2_cov$pos)]


# identify contigs with high coverage -------------------------------------

sd <- sd(cov$cov, na.rm = TRUE)
avg <- mean(cov$cov, na.rm = TRUE)

mv_producing_contigs <- cov %>% 
    filter(!is.na(cov)) %>% 
    filter(cov >= (avg + 2 * sd)) %>% 
    select(pos) %>% 
    mutate(contig_id = str_remove(pos, "\\_\\d*$")) %>% 
    select(contig_id) %>% 
    unique()


# save identified contigs -------------------------------------------------

info <- paste0(nrow(mv_producing_contigs), " contigs showed coverage higher than two times standart deviation above the mean.")
write(info, file = output_path)
fwrite(mv_producing_contigs, file = output_path, sep = "\t", append = TRUE)










