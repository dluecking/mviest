# title:
# author: domi
# date: Tue Jul 20 10:15:23 2021

# libraries ---------------------------------------------------------------
library(data.table, quietly = TRUE)
suppressMessages(library(dplyr))
library(stringr, quietly = TRUE)
library(seqinr, quietly = TRUE)


# handling arguments ------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)

if(length(args) != 2){
    print("Usage:\nsummarize_contig_info.R <sample_name> <long_contig_cutoff>")
    quit()
}

sample = args[1]
long_contig_cutoff  = as.numeric(args[2])
# local
# sample = "mixed_positive_test"
output_dir = paste0("results/", sample, "/")



# initialize contig_df ----------------------------------------------------
sequence_file <- paste0(output_dir, "contigs/virome.contigs.filtered.fasta")
sequences <- seqinr::read.fasta(sequence_file, seqtype = "DNA")

contig_df <- data.table("contig_id" = getName(sequences),
                        "length" = getLength(sequences),
                        "circular" = NA,
                        "plasmid" = NA,
                        "vs2" = NA,
                        "dvf" = NA,
                        "mv" = NA,
                        "long_contig" = NA, 
                        "final_label" = NA)

rm(sequence_file, sequences)


# get circularity info ----------------------------------------------------
# circ_df_path <- paste0(output_dir, "circularity.tsv")
# circ_df <- fread(circ_df_path)
# names(circ_df) <- c("id", "circular")

# contig_df$circular <- circ_df$circular[match(contig_df$contig_id, circ_df$id)]

# rm(circ_df, circ_df_path)
contig_df$circular <- NA


# gather platon information -----------------------------------------------
platon_df_path <- paste0(output_dir, "/platon_out/platon_plasmids.tsv")
platon_df <- fread(platon_df_path)
platon_df$plasmid <- "TRUE"

contig_df$plasmid <- platon_df$plasmid[match(contig_df$contig_id, platon_df$ID)]

rm(platon_df, platon_df_path)


# gather vs2 pipeline information ---------------------------------------
vs_df_path <- paste0(output_dir, "/vs2_final_report.tsv")
vs_df <- fread(vs_df_path)
vs_df$short_id <- str_replace(string = vs_df$contig_id, pattern = "\\|\\|.*$", replacement = "")

contig_df$vs2 <- vs_df$label[match(contig_df$contig_id, vs_df$short_id)]
contig_df$vs2[contig_df$vs2 == "true viral"] <- "TRUE"
contig_df$vs2[contig_df$vs2 != "TRUE"] <- NA

rm(vs_df, vs_df_path)


# gather dvf information --------------------------------------------------
dvf_path <- paste0(output_dir, "/dvf/virome.contigs.filtered.fasta_gt1bp_dvfpred.txt")
dvf_df <- fread(dvf_path)

dvf_df$viral <- NA
dvf_df$viral[dvf_df$score > 0.95] <- "TRUE"

contig_df$dvf <- dvf_df$viral[match(contig_df$contig_id, dvf_df$name)]

rm(dvf_df, dvf_path)


# assign final label ------------------------------------------------------
contig_df <- contig_df %>% 
    mutate(final_label = replace(final_label, length>=long_contig_cutoff, "long_contig"))

# plasmids are plasmids
contig_df$final_label[contig_df$plasmid  == "TRUE"] <- "plasmid"
# anything that is viral become viral, including previously assigned plasmids
contig_df$final_label[contig_df$vs2 == "TRUE" | contig_df$dvf == "TRUE"] <- "viral"
# rest is mv
contig_df$final_label[is.na(contig_df$final_label)] <- "mv"


# write contig_DF ---------------------------------------------------------
fwrite(contig_df, paste0(output_dir, "/contig_summary_", sample, ".tsv"))

