# author: dlueckin
# Wed Jan 12 11:38:22 2022 

  
# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(ggplot2)



# handling of arguments ---------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 1 | length(args) > 2){
    print("Usage:\nRscript mv_profiles.R <init/info/update> <path/to/results>")
    quit()
}
MODE = args[1]
RESULTS_DIR = args[2]

if(MODE == "update" & is.na(RESULTS_DIR)){
    print("Missing results directory.")
    print("Usage:\nRscript mv_profiles.R <init/info/update> <path/to/results>")
    quit()
}

# paths -------------------------------------------------------------------

LIST_OF_MV_PRODUCER_PATH = "data/mv_producer/list_of_mv_producer.list"
PROFILE_SUMMARY_PATH = "data/mv_producer/profiles_summary.tsv"
PROFILE_DIR_PATH = "data/mv_producer/profiles/"
PROFILE_LOG_PATH = "data/mv_producer/profile_log.log"
GENOME_DIR_PATH = "data/mv_producer/genomes/"

LOG_MODE = "" # either append or create


# initialize --------------------------------------------------------------

if(MODE == "init" | MODE == "initialize"){
    
    # check that no profile summary exists 
    if(file.exists(PROFILE_SUMMARY_PATH)){
        print("Profile summary found!")
        print(paste0("'", PROFILE_SUMMARY_PATH, " exists"))
        
        print("Overwrite and intialize new summary? <y/n> ")
        answer <- readLines("stdin", n = 1)
        if(answer != "y"){
            print("Aborting...")
            quit()
        }
    }
    
    # check that no profiles exist
    if(length(dir(all.files=TRUE)) != 0){
        print("Profile directory not empty:")
        print(paste0("ls ", PROFILE_DIR_PATH))
        system(paste0("ls ", PROFILE_DIR_PATH))
        
        print("Overwrite and intialize new profiles? <y/n> ")
        answer <- readLines("stdin", n = 1)
        if(answer != "y"){
            print("Aborting...")
            quit()
        }
    }
    
    # check that no list_of_producers exists
    if(file.exists(LIST_OF_MV_PRODUCER_PATH)){
        print("List of mv producers exists!")
        
        print( "Overwrite and intialize list? <y/n> ")
        answer <- readLines("stdin", n = 1)
        if(answer != "y"){
            print("Aborting...")
            quit()
        }
    }
    
    # is there a log file? If yes, append or create new?
    if(file.exists(PROFILE_LOG_PATH)){
        print("Log file for profiles exists and we are in <init> mode!")
        print("Create new log or append to existing? <append/create> ")
        user_input <- readLines("stdin", n = 1)
        
        if(user_input == "create" | user_input == "append"){
            LOG_MODE = user_input
        }else{
            print("Not a valid answer. Aborting...")
            quit()
        }
    } else{
        LOG_MODE = "create"
    }
    
    # create new list of producers
    list_of_mv_producer <- list.files(GENOME_DIR_PATH, pattern = "\\.fna$")
    write(list_of_mv_producer, LIST_OF_MV_PRODUCER_PATH)
    
    # initialize datatable for each genome
    for(genome in list_of_mv_producer){
        genome_ID <- tools::file_path_sans_ext(genome)
        genome_path = paste0(GENOME_DIR_PATH, genome_ID, ".fna")
        profile_path = paste0(PROFILE_DIR_PATH, genome_ID, ".profile")
        
        genome_length <- as.numeric(system(paste0("cat ", genome_path, " | grep  -v '>' | wc -m"), intern = TRUE))
        
        profile_DT <- data.table("genome_ID" = genome_ID,
                                 "pos" = c(1:genome_length),
                                 "cov" = 0)
        fwrite(profile_DT, file = profile_path, sep = "\t")
    }
    

    # initialize profiles_summary
    profile_summary_DT <- data.table("genome_ID" = list_of_mv_producer,
                                     "updated_n_times" = 0)                     # I could add more stats here!
    fwrite(profile_summary_DT, PROFILE_SUMMARY_PATH, sep = "\t")
    
    # log message
    log_message = paste0(date(), " - initialized profiles for ", length(list_of_mv_producer), " genomes.")
    if(LOG_MODE == "create"){
        write(log_message, PROFILE_LOG_PATH)
    }
    if(LOG_MODE == "append"){
        write(log_message, PROFILE_LOG_PATH, append = TRUE)
    }
}


# update ------------------------------------------------------------------

if(MODE == "update"){
    SAMPLES <- list.files(RESULTS_DIR)
    NO_OF_SAMPLES <- length(SAMPLES)
    print(paste0("Results from ", NO_OF_SAMPLES, " samples detected:"))
    for(s in SAMPLES){
        print(s)
    }
    
    # get user input
    print("Provide a sample ID (or multiple, seperated by a comma, e.g. 'test1,med4,fish2', or 'all' to updated everything) for which you want to update the profiles.")
    USER_PROVIDED_SAMPLES <- str_split(readLines("stdin", n = 1), ",")
    
    if(USER_PROVIDED_SAMPLES[1] == "all"){
        USER_PROVIDED_SAMPLES = SAMPLES


    }
    # iterate over all samples
    for(sample in USER_PROVIDED_SAMPLES){
        # catch if user input is not an existing sample
        if(!sample %in% SAMPLES){
            print(paste0(sample, " is not among the detected samples. Aborting..."))
            quit()
        }        
        
        # list all mappings
        sample_mappings <- list.files(paste0(RESULTS_DIR, sample, "/mappings/mv_producer"), pattern = "coverage$")
        
        # check which profiles would be required for update
        required_profiles <- str_extract(sample_mappings, "vs\\_.*\\.fna") # extracts the reference name of the mapping
        required_profiles <- str_remove(required_profiles, "vs\\_") # removes the "vs_"
        required_profiles <- unique(required_profiles)
        
        profile_summary <- fread(PROFILE_SUMMARY_PATH)
        
        # check if all the profiles are existing
        for(profile in required_profiles){
            if(! profile %in% profile_summary$genome_ID){
                print(paste0("No profile is found for ", profile))
                print("Initialize the profiles again. Aborting...")
                quit()
            }
        }
        
        # update profiles for each mapping
        for(mapping in sample_mappings){
            # read the new mapping
            local_mapping <- fread(paste0(RESULTS_DIR, sample, "/mappings/mv_producer/", mapping))
            
            # which profile do I update?
            profile_to_update <- str_extract(mapping, "vs\\_.*\\.fna")
            profile_to_update <- str_remove(profile_to_update, "vs\\_") 
            profile_to_update <- str_remove(profile_to_update, "\\.fna")
            
            # read profile
            profile_path <- paste0(PROFILE_DIR_PATH,  profile_to_update, ".profile")
            profile <- fread(profile_path)
            
            # update profile
            profile$cov <- profile$cov + local_mapping$V3[match(profile$pos, local_mapping$V2)]
            
            # save the profile
            fwrite(profile, file = profile_path, sep = "\t")
            
            # update summary, WTF HOW IS THAT THE ONLY WAY TO DO this, why is there no +=?
            profile_summary$updated_n_times[profile_summary$genome_ID == paste0(profile_to_update, ".fna")] = profile_summary$updated_n_times[profile_summary$genome_ID == paste0(profile_to_update, ".fna")] + 1
            
        }
        # save the profile summary
        fwrite(profile_summary, PROFILE_SUMMARY_PATH)
        
        # logging
        log_message = paste0(date(), " - updated profiles for sample '", sample, "'.")
        write(log_message, PROFILE_LOG_PATH, append = TRUE)
    }
    print("Done.")
}


# info --------------------------------------------------------------------

if(MODE == "info"){
    print("##################")
    print("-- 1. Summary: ---")
    
    if(file.exists(PROFILE_SUMMARY_PATH)){
        profile_summary_DT <- fread(PROFILE_SUMMARY_PATH)
        print(profile_summary_DT)
    } else {
        print("No summary file found.")
    }
    
    print("-- 2. Log: -------")
    
    if(file.exists(PROFILE_LOG_PATH)){
        log <- readLines(PROFILE_LOG_PATH)
        print(log)    
    } else {
        print("No log file found.")
    }
}

if(MODE != "init" && MODE !="info" && MODE != "update"){
    print("Usage: Rscript mv_profiles.R <init/info/update>")
    quit()
}