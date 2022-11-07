# author: dlueckin
# Thu Jun 16 13:45:23 2022 

  
# libraries ---------------------------------------------------------------
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpmisc))
suppressMessages(library(patchwork))
suppressMessages(library(RColorBrewer))

# handling arguments ------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)

if(length(args) != 7){
    print("Usage: mviest_plot.R {input.kaiju_mvome_table} {input.kaiju_virome_table} {input.viromeQC_mvome} {input.viromeQC_virome} {input.contig_summary} {input.mvome_reads_vs_metagenome_rpkm} {output}")
    quit()
}

kaiju_table_mvome <- args[1]
kaiju_table_virome <- args[2]
viromeQC_mvome <- args[3]
viromeQC_virome <- args[4]
contig_summary <- args[5]
mvome_reads_vs_metagenome_scafstats <- args[6]
outfile <- args[7]


# kaiju_table_mvome <- '4976_AA/kaiju/mvome_summary.tsv'
# kaiju_table_virome <- '4976_AA/kaiju/virome_summary.tsv'
# viromeQC_mvome <- '4976_AA/viromeQC/mvome_QC.tsv'
# viromeQC_virome <- '4976_AA/viromeQC/virome_QC.tsv'
# contig_summary <- '4976_AA/contig_summary_4976_AA.tsv'
# mvome_reads_vs_metagenome_scafstats <- "4976_AA/mappings/mvome_positive_vs_metagenome/scafstats.txt"
# outfile <- '4976_AA/4976_AA_mviest_plot_local.png'

SAMPLE <- unlist(strsplit(outfile, "\\/"))[2]

# plot 1 - viral vs non-viral ratio ---------------------------------------

NO_OF_TOTAL_CONTIGS <- nrow(fread(contig_summary))

if(file.size(viromeQC_mvome) > 0){
    NO_OF_READS_MVOME <- fread(viromeQC_mvome) %>% 
        select(Reads) %>% 
        top_n(1) %>% 
        unlist()
}else
    NO_OF_READS_MVOME <- 0

if(file.size(viromeQC_virome) > 0){
    NO_OF_READS_VIROME <- fread(viromeQC_virome) %>% 
        select(Reads) %>% 
        top_n(1) %>% 
        unlist()
}else
    NO_OF_READS_VIROME <- 0

if(file.size(viromeQC_mvome)){
    VIROMEQC_SCORE_MVOME <- fread(viromeQC_mvome) %>% 
        select(`total enrichmnet score`) %>% 
        top_n(1) %>% 
        unlist()
}else
    VIROMEQC_SCORE_MVOME <- 0

if(file.size(viromeQC_virome)){
    VIROMEQC_SCORE_VIROME <- fread(viromeQC_virome) %>% 
        select(`total enrichmnet score`) %>% 
        top_n(1) %>% 
        unlist()
}else
    VIROMEQC_SCORE_VIROME <- 0


mini_df <- data.table(sample = c(SAMPLE, SAMPLE),
                      label = c('viral', 'non-viral'), 
                      no_of_reads = unlist(c(NO_OF_READS_VIROME, NO_OF_READS_MVOME)),
                      viromeQC_scores = unlist(c(VIROMEQC_SCORE_VIROME, VIROMEQC_SCORE_MVOME)))

mini_df$viromeQC_scores <- round(mini_df$viromeQC_scores, 2)

mini_df <- mini_df %>% 
    mutate(viromeQC_label_position = (cumsum(mini_df$no_of_reads) - 0.5 * mini_df$no_of_reads) / sum(mini_df$no_of_reads) )


plot_1 <- ggplot(mini_df, aes(fill=label, y=no_of_reads, x = sample)) +
    geom_bar(position="fill", stat="identity", color = "black", alpha = 0.8) +
    scale_fill_manual(values = c("non-viral" = "#009FB7",
                                 "viral" = "#FED766")) +
    geom_text(aes(y = viromeQC_label_position, label = paste(format(viromeQC_scores, nsmall = 2), "\nviromeQC score")), size = 4) +
    theme_minimal() +
    ylab("fraction of reads [%]") +
    xlab(SAMPLE) +
    ggtitle(SAMPLE, subtitle = paste0("no of contigs = ", NO_OF_TOTAL_CONTIGS, ",\nno of mapped reads = ", sum(NO_OF_READS_MVOME, NO_OF_READS_VIROME))) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())


# kaiju plots -------------------------------------------------------------

if(file.size(kaiju_table_mvome) > 0){
    kiaju_mvome_df <- fread(kaiju_table_mvome) %>% 
        arrange(desc(percent)) %>% 
        filter(percent >= 1) %>% 
        mutate(sample = SAMPLE) %>% 
        mutate(origin = 'non-viral')
}else
    kiaju_mvome_df <- data.table()

if(file.size(kaiju_table_virome) > 0){
    kiaju_virome_df<- fread(kaiju_table_virome) %>%
        arrange(desc(percent)) %>% 
        filter(percent >= 1) %>%  
        mutate(sample = SAMPLE) %>% 
        mutate(origin = 'viral')
}else
    kiaju_virome_df <- data.table()

kaiju_plot_df <- rbind(kiaju_virome_df, kiaju_mvome_df)

getPalette <- colorRampPalette(brewer.pal(11, 'PiYG'))
kaiju_plot_df$color <- getPalette(nrow(kaiju_plot_df))
kaiju_plot_df$color[kaiju_plot_df$taxon_name == 'Viruses'] <- "#FED766"
kaiju_plot_df$color[kaiju_plot_df$taxon_name == 'unclassified'] <- "grey"
colors <- kaiju_plot_df$color
names(colors) <- kaiju_plot_df$taxon_name

kaiju_plot <- ggplot(kaiju_plot_df, aes(x = origin, y = percent, fill = taxon_name)) +
    geom_bar(position="fill", stat="identity", color = "black", alpha = 0.8) +
    theme_minimal() +
    scale_fill_manual(values = colors) +
    ylab("fraction of reads [%]") +
    xlab("") +
    ggtitle("", subtitle = 'reads classified by kaiju')

# mv producer table  -------------------------------------------------------

if(file.size(mvome_reads_vs_metagenome_scafstats) > 0){
    mv_prod_df <- fread(mvome_reads_vs_metagenome_scafstats) %>% 
        arrange(desc(assignedReads)) %>% 
        select('#name', 'assignedReads') %>% 
        head(10)
}else
    mv_prod_df <- data.table("V1" == "no scafstats calculated.")


fwrite(mv_prod_df, str_replace(outfile, pattern = '\\_mviest\\_plot\\.png', '_potential_mv_producer.tsv'))


# patchwork plots ---------------------------------------------------------

final_plot <- plot_1 +  kaiju_plot 


# save plots --------------------------------------------------------------

ggsave(paste0(outfile), plot = final_plot, width = 25, height = 15, units = "cm")


# save data ---------------------------------------------------------------

summary_df <- data.table(
    'sample' = SAMPLE,
    'no_of_contigs' = NO_OF_TOTAL_CONTIGS,
    'no_of_reads_viral' = NO_OF_READS_VIROME,
    'no_of_reads_non-viral' = NO_OF_READS_MVOME,
    'viromeQC_score_viral' = VIROMEQC_SCORE_VIROME,
    'viromeQC_score_non-viral' = VIROMEQC_SCORE_MVOME,
    'mviest_mvome_ratio' = 0
)
summary_df$mviest_mvome_ratio = NO_OF_READS_MVOME / sum(NO_OF_READS_VIROME, NO_OF_READS_MVOME)

fwrite(summary_df, str_replace(outfile, pattern = '\\_mviest\\_plot\\.png', '_mviest_summary.tsv'))
