# title:
# author: domi
# date: Mon Aug  9 14:50:11 2021

# libraries ---------------------------------------------------------------
suppressMessages(library(dplyr))
library(data.table, quietly = TRUE)
library(stringr, quietly = TRUE)
library(seqinr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
# library(patchwork, quietly = TRUE)


# handling arguments ------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)

if(length(args) != 1){
    print("Usage:\nsummarize_contig_info.R <sample_name>")
    quit()
}

sample = args[1]
# local
output_dir = paste0("results/", sample, "/")


# load contig_df ----------------------------------------------------------

contig_df <- fread(paste0(output_dir, "/contig_summary_", sample, ".tsv"))


# plot 1: categorization --------------------------------------------------
labels <- c("long_contig", "mv", "viral", "plasmid")

plot_df <- data.frame("sample" = rep(sample, 4),
                      "label" = labels,
                      "count" = 0,
                      "fraction" = 0)


for(label in labels){
    plot_df$count[plot_df$label == label] <- contig_df %>%
        filter(final_label == label) %>%
        nrow()
}
rm(label)

plot_df$fraction <- plot_df$count/sum(plot_df$count)

plot_1 <- ggplot(plot_df, aes(fill=label, y=fraction, x=sample)) +
    geom_bar(position="fill", stat="identity", color = "black", alpha = 0.9) +
    scale_fill_manual(values = c("long_contig" = "#E6E6EA",
                                 "mv" = "#009FB7",
                                 "viral" = "#FED766",
                                 "plasmid" = "#FE4A49")) +
    theme_minimal() +
    ylab("fraction [%]") +
    xlab("") +
    ggtitle(paste0("MV fraction for ", sample), subtitle = paste0("n = ", nrow(contig_df))) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

rm(plot_df)

# THIS IS A shortcut, since I don't want to include the COG stuff righ now. 
ggsave(paste0(output_dir, "/contig_summary_", sample, ".png"), plot = plot_1, width = 10, height = 10, units = "cm")


# # gather COG information --------------------------------------------------
# cog_df_path <- paste0(output_dir, "/eggnog/proteins.emapper.annotations")
# cog_df <- fread(cog_df_path)
# cog_df$contig_id <- str_replace(cog_df$`#query`, "\\_\\d*$", "")
# 
# cog_df$final_label <- contig_df$final_label[match(cog_df$contig_id, contig_df$contig_id)]
# 
# rm(cog_df_path)
# 
# 
# # open COG translation
# cog_translation <- fread("data/COG_translation_list.csv", header = FALSE)
# names(cog_translation) <- c("Category", "Translation")
# extra_cat <- data.table("Category" = "",
#                         "Translation" = "not assigned")
# cog_translation <- rbind(cog_translation, extra_cat)
# rm(extra_cat)
# 
# 
# # COG plots ---------------------------------------------------------------
# plotCOG <- function(label, color, cog_df, cog_translation){
#     # gather data
#     local_df <- cog_df %>%
#         filter(final_label == label)
# 
#     local_plot_df <- data.frame("category" = cog_translation$Category,
#                                 "translation" = cog_translation$Translation,
#                                 "count" = 0,
#                                 "fraction" = 0)
# 
#     for(cat in cog_translation$Category){
#         count = local_df %>%
#             filter(COG_category == cat) %>%
#             nrow()
# 
#         local_plot_df$count[local_plot_df$category == cat] <- count
#     }
#     rm(cat)
# 
#     # plot
#     plot <- ggplot(local_plot_df, aes(x = category, y = count)) +
#         geom_bar(stat = "identity", fill = color, color = "black", alpha = 0.9) +
#         theme_minimal() +
#         ggtitle("", subtitle = paste0(label, " (n = ", nrow(contig_df[final_label == label]), "), COG hits = ", sum(local_plot_df$count))) +
#         theme(legend.position = "none")
#     return(plot)
# }
# 
# aux_df  <- data.frame("label" = c("long_contig", "mv", "viral", "plasmid"),
#                       "color" = c("#E6E6EA", "#009FB7", "#FED766", "#FE4A49"))
# 
# 
# plot_2_lc <- plotCOG(aux_df$label[1], aux_df$color[1],
#                      cog_df = cog_df,
#                      cog_translation = cog_translation)
# plot_2_mv <- plotCOG(aux_df$label[2], aux_df$color[2],
#                      cog_df = cog_df,
#                      cog_translation = cog_translation)
# 
# plot_2_vi <- plotCOG(aux_df$label[3], aux_df$color[3],
#                      cog_df = cog_df,
#                      cog_translation = cog_translation)
# 
# plot_2_pl <- plotCOG(aux_df$label[4], aux_df$color[4],
#                      cog_df = cog_df,
#                      cog_translation = cog_translation)
# 
# plot_2 <- plot_2_lc / plot_2_mv / plot_2_vi / plot_2_pl



# calculate total pR1SE hits ----------------------------------------------

# total_pR1SE_hits <- sum(contig_df$pR1SE_hits, na.rm = TRUE)



# combine and save plots --------------------------------------------------

# 
# final_plot <- plot_1 + plot_2 +
#     plot_layout(widths = c(1, 5)) +
#     plot_annotation(
#         title = paste0("pore summary for: ", sample),
#         subtitle = paste0("number of contigs: ", nrow(contig_df), "\n", 
#                           "pR1SE hits: ", total_pR1SE_hits)
#     )
# ggsave(paste0(output_dir, "/pore_summary.png"), plot = final_plot, width = 20, height = 20, units = "cm")

