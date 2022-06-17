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

# kaiju_table_mvome <- 'results/sargasso/kaiju/mvome_summary.tsv'
# kaiju_table_virome <- 'results/sargasso/kaiju/virome_summary.tsv'
# viromeQC_mvome <- 'results/sargasso/viromeQC/mvome_QC.tsv'
# viromeQC_virome <- 'results/sargasso/viromeQC/virome_QC.tsv'
# contig_summary <- 'results/sargasso/contig_summary_sargasso.tsv'
# mvome_reads_vs_metagenome_scafstats <- "results/sargasso/mappings/mvome_positive_vs_metagenome/scafstats.txt"
# outfile <- 'results/sargasso/sargasso_mviest_plot.png'

SAMPLE <- unlist(strsplit(outfile, "\\/"))[2]

# plot 1 - viral vs non-viral ratio ---------------------------------------

NO_OF_TOTAL_CONTIGS <- nrow(fread(contig_summary))
NO_OF_READS_MVOME <- fread(viromeQC_mvome) %>% 
    select(Reads) %>% 
    top_n(1)
NO_OF_READS_VIROME <- fread(viromeQC_virome) %>% 
    select(Reads) %>% 
    top_n(1)
VIROMEQC_SCORE_MVOME <- fread(viromeQC_mvome) %>% 
    select(`total enrichmnet score`) %>% 
    top_n(1)
VIROMEQC_SCORE_VIROME <- fread(viromeQC_virome) %>% 
    select(`total enrichmnet score`) %>% 
    top_n(1)

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

kiaju_mvome_df <- fread(kaiju_table_mvome) %>% 
    arrange(desc(percent)) %>% 
    top_n(5, wt = percent) %>% 
    mutate(sample = SAMPLE) %>% 
    mutate(origin = 'non-viral')


kiaju_virome_df<- fread(kaiju_table_virome) %>%
    arrange(desc(percent)) %>% 
    top_n(5, wt = percent) %>% 
    mutate(sample = SAMPLE) %>% 
    mutate(origin = 'viral')

kaiju_plot_df <- rbind(kiaju_virome_df, kiaju_mvome_df)

kaiju_plot_df$color <- brewer.pal(n = nrow(kaiju_plot_df), name = 'PiYG')
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
    ggtitle("", subtitle = 'i)                                                ii)')

# mv producer plots -------------------------------------------------------

mv_prod_df <- fread(mvome_reads_vs_metagenome_scafstats) %>% 
    arrange(desc(assignedReads)) %>% 
    select('#name', 'assignedReads')

mv_prod_df <-  mv_prod_df %>% 
    head(5)

kaiju_plot <- kaiju_plot +
    annotate(geom = 'table',
             x=3,
             y=0.75,
             label=list(mv_prod_df), hjust = 0, vjust = 0)


# patchwork plots ---------------------------------------------------------

final_plot <- plot_1 + kaiju_plot


# save plots --------------------------------------------------------------

ggsave(paste0(outfile), plot = final_plot, width = 25, height = 15, units = "cm")

