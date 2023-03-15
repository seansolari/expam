require(dplyr)


# loading expam output files ----------------------------------------

summary_table_names <- c("classified.csv", "split.csv")
find_modes <- c("phy", "tax")
col_modes <- list("cumulative" = c("Cumulative Classified Count", "Cumulative Split Count"),
                  "raw" = c("Raw Classified Count", "Raw Split Count"))


validate_mode <- function(mode, valid_modes) {
    if (!(mode %in% valid_modes)) {
        stop("Unknown mode: %s", mode)
    }
}

get_sample_cols <- function(sample_table_path, cols) {
    readr::read_csv(sample_table_path, show_col_types = F) %>%
        rename(clade = "...1") %>%
        select(c("clade", cols))
}

find_expam_output <- function(base_path, parent_folder, is_summary = T) {
    # find parent folder of results
    subdirs <- grep(
        sprintf("%s$", parent_folder),
        list.dirs(path = base_path, full.names = T, recursive = T),
        value = TRUE
    )
    tables <- unlist(sapply(
        subdirs,
        function(p) list.files(p, include.dirs = F, full.names = T),
        USE.NAMES = F
    ))
    
    # only return valid files
    summary_mask <- function(f) { basename(f) %in% summary_table_names }
    sample_mask <- function(f) !summary_mask(f) & (basename(f) != "raw")
    masker <- if (is_summary) summary_mask else sample_mask
    
    # mask appropriate files
    tables[masker(tables)]
}

get_samples <- function(base, results_mode = "phy", col_mode = "raw") {
    validate_mode(results_mode, find_modes)
    validate_mode(col_mode, names(col_modes))
    
    # function to collect column from csv
    col_getter <- function(fpath) get_sample_cols(fpath, col_modes[[col_mode]])
    # collect all paths to tables
    table_paths <- find_expam_output(base, results_mode, is_summary = F)
    
    data <- sapply(table_paths, col_getter, USE.NAMES = T, simplify = F)
    names(data) <- gsub("\\.csv", "", basename(names(data)))
    
    data
}

get_summaries <- function(base, results_mode = "phy") {
    validate_mode(results_mode, find_modes)
        
    # group classified.csv and split.csv per sample
    summary_files <- find_expam_output(base, results_mode, is_summary = T)
    grouped_summaries <- split(
        summary_files,
        stringr::str_match(summary_files, ".+\\/(.+)\\/phy")[,2]
    )
    
    # open tables
    lapply(grouped_summaries,
           function(fvec) {
               names(fvec) <- gsub("\\.csv", "", basename(fvec))
               lapply(fvec,
                      function(f) readr::read_csv(f, show_col_types = F)%>%
                          rename(clade = "...1")
               )
           }
    )
}

merge_samples <- function(sample_list) {
    lapply(1:length(sample_list),
           function(i) {
               sample_tibl <- sample_list[[i]] %>%
                   summarise(clade = clade, count = rowSums(Filter(is.numeric, sample_list[[i]])))
               names(sample_tibl) <- c("clade", names(sample_list)[[i]])
               
               sample_tibl
           }) %>%
        purrr::reduce(dplyr::full_join, by = "clade") %>%
        mutate_if(is.numeric, ~replace(., is.na(.), 0))
}

merge_summaries <- function(summary_list) {
    lapply(summary_list,
          function(table_list) do.call(merge_summary_splits_counts, table_list)) %>%
        purrr::reduce(dplyr::full_join, by = "clade") %>%
        mutate_if(is.numeric, ~replace(., is.na(.), 0))
}

merge_summary_splits_counts <- function(classified, split) {
    # count and split matrices have the same structure, aside from `unclassified` row
    clades_equal = all(classified$clade[2:nrow(classified)] == split$clade)
    if (!clades_equal) stop("Compromised table structure!")
    
    classified[2:nrow(classified),2:ncol(classified)] <- classified[2:nrow(classified),2:ncol(classified)] + split[,2:ncol(split)]
    classified
}


# Data loading ----------------------------------------

read_expam_classified <- function(counts_path)
  read.csv(counts_path, header = T, row.names = 1)


# Stats transforms ----------------------------------------

#` 
#`  counts is a matrix of (clades, samples).
clr <- function(counts, offset = 1) {
  log_counts <- log2(counts + offset)
  sweep(log_counts, 2, colSums(log_counts) / nrow(counts))
}



