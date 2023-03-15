## Cleaning and pre-processing *expam* output

> **_NOTE_:** This workbook assumes you have run *expam* on the samples and have output ready to be analysed.

``` R
source("scripts/preprocessing/preprocessing.R")
```

Functions provided by `preprocessing.R`

**get_samples** - The separate *expam* outputs from each batch of samples was placed in a single parent directory. This function returns an R list of outputs from each batch, and feeds in to `merge_samples`.

**merge_samples** - Takes a list of *expam* outputs and merges the sample counts into a single table.

**merge_summaries** - Merges the *expam* run summary files from separate batches.

> **_NOTE:_** These helper scripts are only necessary as we separated samples into batches, and are not required if you processed all metagenomic sequence files in a single run.

### Importing datasets

``` R
# input files
batch_base <- "data/expam_output/formatted"
metadata_file <- "data/metadata.csv"

# mapping metagenome file names to sample names
rename_file <- "data/expam_output/rename_metagenomes.csv"

# output files
out_base <- "data/expam_output/clean"
cumulative_out <- paste(out_base, "cumulative.csv", sep = "/")
```

Read csv files

``` R
# load cumulative counts at each clade
samples <- get_samples(
    batch_base,
    results_mode = "phy",
    col_mode = "cumulative")

# merge batches
counts <- merge_samples(samples)

# metadata
meta <- readr::read_csv(metadata_file)
rename <- readr::read_csv(rename_file, col_names = c("old", "new"))
```

Align sample and clade counts with metadata

``` R
colmap <- function(query, keys, values)
{
    idx <- match(query, keys)
    mask <- !is.na(idx)
    idx[mask] <- values[idx[mask]]
    return(idx)
}

# map names with rename file
clean_names <- colmap(colnames(counts), rename$old, rename$new)
clean_names[1] <- "clade"   # keep clade names

# remove any that don't have corresponding entries
col_mask <- !is.na(clean_names)
renamed_counts <- counts[,col_mask]
colnames(renamed_counts) <- clean_names[col_mask]

# align output with corresponding metadata
sample_names <- colmap(colnames(renamed_counts), meta$metagenome, meta$sample)
sample_names[1] <- "clade"  # keep clade names
col_mask <- !is.na(sample_names)

X_tbl <- renamed_counts[,col_mask]
colnames(X_tbl) <- sample_names[col_mask]

# there should be 231 samples remaining
dim(X_tbl)
rm(renamed_counts)

# save cumulative counts to disk
readr::write_csv(X_tbl, cumulative_out)
```

### Employ counts cutoff for spurious read assignments

Employ '5 counts-per-million reads in the sample' cutoff, for each sample

``` R
CPM <- 5

# load summary files and get million counts per sample
expam_summary <- merge_summaries(get_summaries(batch_base, results_mode = "phy"))
million <- colSums(expam_summary[,-1])
# map sample names
million_mapped_names <- colmap(names(million), rename$old, rename$new)
mask <- !is.na(million_mapped_names)
million <- million[mask]
names(million) <- colmap(million_mapped_names[mask], meta$metagenome, meta$sample)

rm(expam_summary)
hist(million)

# reorder million according to order of samples in `X_tbl`
million <- million[match(colnames(X_tbl[,-1]), names(million))]

# get counts matrix as counts per million
X_cpm <- as.matrix(X_tbl[,-1]) / (million / 1e+06)
mask <- (X_cpm < CPM) & (X_cpm > 0)

X_mat <- as.matrix(X_tbl[,-1])
rownames(X_mat) <- X_tbl %>% pull(clade)
X_mat[mask] <- 0

# now that false classifications have been zeroed out,
# remove any null rows
mask <- apply(X_mat, MARGIN = 1, FUN = function(r) any(r > 0))
X_mat <- X_mat[mask,]

# save cutoff to disk
matrix_file_name <- paste(out_base, sprintf("cumulative_cpm%d.csv", CPM), sep = "/")
write.csv(X_mat, matrix_file_name, quote = F)
```

**Next step...**

[Integrating microbiome and host data using sparse Partial Least Squares](./3_bulk_sPLS.md)


