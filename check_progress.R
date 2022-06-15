# author: dlueckin
# Wed Jun  1 14:51:02 2022 

  
# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
suppressMessages(library(dplyr))

# working directory -------------------------------------------------------
setwd("~/projects/mvome_pipeline/mv_signals/")

args = commandArgs(trailingOnly=TRUE)

if(length(args) > 0){
    if(args[1] == "long"){
        PRINT_ALL <- TRUE

    }
}else{
    PRINT_ALL <- FALSE
}

# list samples ------------------------------------------------------------

samples <- list.files("results/")

df <- data.table("sample" = samples,
                 "completed" = "",
                 "progress" = "",
                 "v_contigs" = FALSE,
                 "v_contigs_filtered" = FALSE,
                 "m_contigs" = FALSE,
                 "m_contigs_filtered" = FALSE,
                 "reads" = FALSE,
                 "reads_mvome" = FALSE,
                 "dvf" = FALSE,
                 "vs2" = FALSE,
                 "checkv" = FALSE,
                 "virus_summary" = FALSE,
                 "platon" = FALSE,
                 "m_mv_producer" = FALSE,
                 "m_metagenome" = FALSE,
                 "summary_png" = FALSE,
                 "summary_tsv" = FALSE,
                 "viromeQC_mvome" = FALSE,
                 "viromeQC_virome" = FALSE,
                 "stats" = FALSE)

for(i in 1:nrow(df)){
    sample = df$sample[i]
    df$v_contigs[i] <- file.exists(paste0("results/", sample, "/contigs/virome.contigs.fasta"))
    df$v_contigs_filtered[i] <- file.exists(paste0("results/", sample, "/contigs/virome.contigs.filtered.fasta"))
    df$m_contigs[i] <- file.exists(paste0("results/", sample, "/contigs/metagenome.contigs.fasta"))
    df$m_contigs_filtered[i] <- file.exists(paste0("results/", sample, "/contigs/metagenome.contigs.filtered.fasta"))
    
    df$reads[i] <- file.exists(paste0("results/", sample, "/reads/virome.reads1.fastq.gz"))
    df$reads_mvome[i] <- file.exists(paste0("results/", sample, "/reads/mvome.reads1.fq.gz"))
    
    df$dvf[i] <- file.exists(paste0("results/", sample, "/dvf/virome.contigs.filtered.fasta_gt2000bp_dvfpred.txt"))
    df$vs2[i] <- file.exists(paste0("results/", sample, "/vs2-pass1/final-viral-score.tsv"))
    df$checkv[i] <- file.exists(paste0("results/", sample, "/checkv/quality_summary.tsv"))
    df$virus_summary[i] <- file.exists(paste0("results/", sample, "/vs2_final_report.tsv"))
    
    df$platon[i] <- file.exists(paste0("results/", sample, "/platon_out/platon_plasmids.tsv"))
    
    df$m_mv_producer[i] <- length(list.files(paste0("results/", sample, "/mappings/mv_producer/")))/2 == 
        length(list.files("data/mv_producer/genomes/", pattern = "fna$"))
    df$m_metagenome[i] <- file.exists(paste0("results/", sample, "/mappings/meta_1/r1_vs_meta.sorted.bam.coverage"))
    
    df$summary_png[i] <- file.exists(paste0("results/", sample, "/contig_summary_", sample,".png"))
    df$summary_tsv[i] <- file.exists(paste0("results/", sample, "/contig_summary_", sample,".tsv"))
    
    df$viromeQC_mvome[i] <- file.exists(paste0("results/", sample, "/viromeQC_mvome_", sample, ".tsv"))
    df$viromeQC_virome[i] <- file.exists(paste0("results/", sample, "/viromeQC_virome_", sample, ".tsv"))
    
    df$stats[i] <- file.exists(paste0("results/", sample, "/stats_summary_", sample, ".txt"))
    
    df$progress[i] <- paste0(rowSums(df[i, -c("sample", "progress", "completed")]), "/", length(df)-3)
    if(df$progress[i] == paste0(length(df)-3, "/", length(df)-3))
        df$completed[i] <- "Yes"
    else
        df$completed[i] <- "No"
}

to_print <- df %>%
    select(sample, completed, progress)

if(PRINT_ALL){
    print(df)
}else{
    print(to_print)
}

