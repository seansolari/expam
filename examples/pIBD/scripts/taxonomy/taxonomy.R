# Constants -----------------------------------------

all_ranks <- c("root", "top", "superkingdom", "kingdom", "subkingdom",
               "superclade", "clade", "subclade", "superphylum", "phylum",
               "subphylum", "superclass", "class", "subclass", "superorder",
               "order", "suborder", "superfamily", "family", "subfamily",
               "supergenus", "genus", "subgenus", "species group",
               "species subgroup", "superspecies", "species", "subspecies",
               "serogroup", "serotype", "superstrain", "strain", "substrain",
               "no rank")


# Auxilliary -----------------------------------------

ravel <- function(t) t[1:length(t)]

map_with_col <- function(query, target, values) {
  if (is.matrix(query)) query <- ravel(query)
  
  mapped <- match(
    query,
    target,
    nomatch = NA)
  
  not_na <- !is.na(mapped)
  mapped[not_na] <- values[mapped[not_na]]
  
  mapped
}

map_table <- function(query, target, values) {
  m <- map_with_col(query, target, values)
  # Return shaped matrix.
  matrix(m, nrow = nrow(query), ncol = ncol(query))
}

matrix_ind_to_row_ind <- function(i, nrows) {
  row_ind <- i %% nrows
  if (row_ind == 0) row_ind <- nrows
  
  row_ind
}


# Expam taxonomic metadata functions -----------------------------------------

expam_config_paths <- function(db_path)
  list(
    accession_ids = sprintf("%s/phylogeny/accession_ids.csv", db_path),
    taxa_rank     = sprintf("%s/phylogeny/taxa_rank.csv", db_path),
    taxid_lineage = sprintf("%s/phylogeny/taxid_lineage_repaired.csv", db_path)
  )

load_taxonomy_files <- function(db_path) {
  file_locations <- expam_config_paths(db_path)
  list(
    accesion_ids  = readr::read_csv(file_locations$accession_ids,
                                    col_names = c("genome", "accession", "taxid")),
    taxa_rank     = readr::read_csv(file_locations$taxa_rank,
                                    col_names = c("taxa", "taxid", "rank")),
    taxid_lineage = readr::read_csv(file_locations$taxid_lineage,
                                    col_names = F)
  )
}


# Pre-processing -----------------------------------------

get_unique_lineages <- function(taxid_lineage) {
  taxa_cols <- ncol(taxid_lineage)
  
  # Order taxonomic lineages alphabetically and by in rank order.
  taxid_lineage <- taxid_lineage %>%
    dplyr::arrange(dplyr::across(dplyr::num_range("X", 2:taxa_cols)))
  
  
  # Take unique taxonomic lineage, i.e. those that are not subsets
  # of the previous entry in an alphabetically ordered table.
  n_taxids        <- nrow(taxid_lineage)
  is_duplicate    <- apply(
    taxid_lineage[1:n_taxids,2:taxa_cols] == taxid_lineage[lag(1:n_taxids),2:taxa_cols],  # Compare with previous row.
    1,                                                                                    # Apply along rows.
    function(r) all(r, na.rm = T))
  is_duplicate[1] <- F  # First row is always unique.
  
  
  # Return unique values.
  taxid_lineage[!is_duplicate,]
}


# Taxonomy Table creation -----------------------------------------

make_lineage_table <- function(.data, propagate = F) {
  taxa_table <- as.matrix(get_unique_lineages(.data$taxid_lineage)[,2:11])
  
  
  taxid_table <- map_table(taxa_table, .data$taxa_rank$taxa, .data$taxa_rank$taxid)
  rank_table <- map_table(taxa_table, .data$taxa_rank$taxa, .data$taxa_rank$rank)
  
  
  # Create a dataframe where columns are ranks
  # and each row is a taxonomic lineage.
  #
  
  # Get a list of active ranks.
  active_ranks = all_ranks[all_ranks %in% unique(ravel(rank_table))]
  
  taxid_data <- t(vapply(
    1:nrow(rank_table),
    function(i) map_with_col(active_ranks, rank_table[i,], taxid_table[i,]),
    numeric(length(active_ranks))))
  colnames(taxid_data) <- active_ranks
  
  if (propagate) {
    taxid_data <- propagate_assignments(taxid_data)
  }
  
  taxid_data
}

propagate_assignments <- function(taxid_data) {
  for (j in ncol(taxid_data):2) {
    na_mask <- is.na(taxid_data[,(j-1)])
    
    if (sum(na_mask) > 0) taxid_data[na_mask,(j-1)] <- taxid_data[na_mask, j]
  }
  
  taxid_data
}


# Assigning/Converting taxonomic labels -----------------------------------------

assign_at_rank <- function(query, taxid_data, tax_rank) {
  # Get range of ranks at or below supplied rank.
  rank_index <- which(all_ranks == tax_rank)
  if( length(rank_index) == 0 )
  {
    taxid_matches <- NA
    print("Warning, rank not known!")
  }
  else
  {
    rank_index <- rank_index[1]
    valid_ranks <- all_ranks[rank_index:length(all_ranks)]
    
    # Get subset of data from these ranks that are active.
    active_ranks <- valid_ranks[valid_ranks %in% colnames(taxid_data)]
    taxid_subset <- taxid_data[,active_ranks]
    
    # Search for requested taxids.
    taxid_matches <- match(query, taxid_subset, nomatch = NA)
    mask_na <- !is.na(taxid_matches)
    taxid_row_inds <- vapply(
      taxid_matches[mask_na],
      function(v) matrix_ind_to_row_ind(v, nrow(taxid_subset)),
      numeric(1),
      USE.NAMES = F
    )
    taxid_matches[mask_na] <- taxid_subset[taxid_row_inds, tax_rank]
  }
  
  taxid_matches
}

assign_species <- function(taxids, taxid_data)
  assign_at_rank(taxids, taxid_data, "species")

assign_genus <- function(taxids, taxid_data)
  assign_at_rank(taxids, taxid_data, "genus")

