# author: dlueckin
# Thu Apr 28 10:07:06 2022 

  
# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(ggplot2)
suppressMessages(library(dplyr))


# arguments ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

if(length(args) != 2){
    print("Usage:\nstats_summary.R <sample_name> <output_file>")
    quit()
}

SAMPLE = args[1]
OUTPUT_FILE  = args[2]
RESULTS_DIR = paste0("results/", SAMPLE, "/")

# local
# SAMPLE = "sargossamock"
# OUTPUT_FILE = paste0(RESULTS_DIR, "stats_summary_", SAMPLE, ".txt")


# general -----------------------------------------------------------------

write("# general", OUTPUT_FILE, append = TRUE)
write(paste0("Sample: ", SAMPLE), OUTPUT_FILE, append = TRUE)
write(paste0("Date: ", date()), OUTPUT_FILE, append = TRUE)



# parameter ---------------------------------------------------------------
write("# parameter", OUTPUT_FILE, append = TRUE)
system(paste0("head -n 14 config.yaml >> ", OUTPUT_FILE))
write("", OUTPUT_FILE, append = TRUE)


# input stats -------------------------------------------------------------
write("# input stats", OUTPUT_FILE, append = TRUE)
CONTIGS_DIR <- paste0(RESULTS_DIR, "contigs/")

vc <- system(paste0("grep '>' ", CONTIGS_DIR, "virome.contigs.fasta | wc -l"), intern = TRUE)
vcf <- system(paste0("grep '>' ", CONTIGS_DIR, "virome.contigs.filtered.fasta | wc -l"), intern = TRUE)
write(paste0("Number of virome contigs: ", vc), OUTPUT_FILE, append = TRUE)
write(paste0("Number of filtered virome contigs: ", vcf), OUTPUT_FILE, append = TRUE)


if(file.exists(paste0(CONTIGS_DIR, "metagenome.contigs.fasta"))){
    mc <- system(paste0("grep '>' ", CONTIGS_DIR, "metagenome.contigs.fasta | wc -l"), intern = TRUE)
    mcf <- system(paste0("grep '>' ", CONTIGS_DIR, "metagenome.contigs.filtered.fasta | wc -l"), intern = TRUE)
    write(paste0("Number of metagenome contigs: ", mc), OUTPUT_FILE, append = TRUE)
    write(paste0("Number of filtered metagenome contigs: ", mcf), OUTPUT_FILE, append = TRUE)
}

r1 <- system(paste0(" echo $(zcat ", RESULTS_DIR, "reads/virome.reads1.fastq.gz | wc -l)/4|bc"), intern = TRUE)
r2 <- system(paste0(" echo $(zcat ", RESULTS_DIR, "reads/virome.reads2.fastq.gz | wc -l)/4|bc"), intern = TRUE)
write(paste0("Number of virome reads: ", r1, " +  ", r2), OUTPUT_FILE, append = TRUE)
write("", OUTPUT_FILE, append = TRUE)


# results -----------------------------------------------------------------

potential_mv_producers <- system(paste0("cat ", RESULTS_DIR, "potential_mv_producer_within_metagenome.tsv | wc -l"), intern = TRUE)
if(potential_mv_producers == "1"){
    potential_mv_producers <- "None"
}
write(paste0("Number of potential mv producers: ", potential_mv_producers), OUTPUT_FILE, append = TRUE)
write("", OUTPUT_FILE, append = TRUE)


# viromeQC enrichment score -----------------------------------------------

mvome_reads_enrichment_score_df <- fread(paste0(RESULTS_DIR, "viromeQC_mvome_", SAMPLE, ".tsv"))
virome_reads_enrichment_score_df <- fread(paste0(RESULTS_DIR, "viromeQC_virome_", SAMPLE, ".tsv"))

mvome_reads_enrichment_score <- mvome_reads_enrichment_score_df$'total enrichmnet score'[1]
virome_reads_enrichment_score <- virome_reads_enrichment_score_df$'total enrichmnet score'[1]

write(paste0("Mvome enrichment score (reads): ", mvome_reads_enrichment_score), OUTPUT_FILE, append = TRUE)
write(paste0("Virome enrichment score (reads): ", virome_reads_enrichment_score), OUTPUT_FILE, append = TRUE)









