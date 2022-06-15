# author: dlueckin
# Wed Nov 17 11:48:59 2021 
# this script lists all reads in data/virome/reads as well as all contigs in data/virome/contigs and
# generates an overview with the structure id - reads - contigs and saves this overview to data/viromes/viromes_table.tsv

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)

# initialize data.table
OUTFILE = "data/viromes/viromes_table.tsv"



# list contig files 
contig_files <- list.files("data/viromes/contigs/")

# list read files
R1_files <- list.files("data/viromes/reads", pattern = "1.fastq.gz")
R2_files <- list.files("data/viromes/reads", pattern = "2.fastq.gz")


# all IDs I find
ids <- str_remove(contig_files, "\\_.*")
ids <- c(ids, str_remove(R1_files, "\\_.*"))
ids <- c(ids, str_remove(R2_files, "\\_.*"))
ids <- unique(ids)


overview_dt <- data.table(id = ids,
                          contigs = as.character(),
                          read1 = as.character(),
                          read2 = as.character())
overview_dt$id <- ids

# fill in the rest of the information
for(i in 1:nrow(overview_dt)){
    current_id <- overview_dt$id[i]
    
    if(paste0(current_id, "_1.fastq.gz") %in% R1_files){
        overview_dt$read1[i] <- paste0(current_id, "_1.fastq.gz")
        
        
    }
    
    
}




# collect ID
fwrite(overview_dt, file = OUTFILE, sep = "\t")