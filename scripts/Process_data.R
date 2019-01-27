# R script for sequence clean up

# Be sure to install these packages before running this script
# They can be installed either with the intall.packages() function
# Or with the 'Packages' pane in RStudio
# Load general-use packages
library("dplyr")
library("tidyr")
library("knitr")
library("ggplot2")
# This package allows for the easy inclusion of literature citations in our Rmd
# More info here: https://github.com/crsh/citr
# And here:
# http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html
library("citr")
# These are the primary packages well use to clean and analyze the data
# This package needs to be installed from bioconductor -- it's not on CRAN
# See info here: https://benjjneb.github.io/dada2/dada-installation.html
library("dada2")
# This to export a fasta of our final denoised sequence variants
library("seqinr")
# To install this you have to install from GitHub
# See more info here: https://github.com/leffj/mctoolsr
# Run this -- install.packages("devtools")
# And then this -- devtools::install_github("leffj/mctoolsr")
library("mctoolsr")
# And this to visualize our results
# It also needs to be installed from bioconductor
library("phyloseq")
# NOTE: Much of the following follows the DADA2 tutorials available here:
# https://benjjneb.github.io/dada2/tutorial.html
# Accessed October 19, 2017
# Set the base path for our input data files
path <- "data/raw_data"

# Sort ensures samples are in order
filenames_forward_reads <- sort(list.files(path, pattern = ".fastq"))

# Extract sample names, assuming filenames have format: SAMPLENAME.fastq
sample_names <- sapply(strsplit(filenames_forward_reads, "\\."), `[`, 1)

# Specify the full path to each of the filenames_forward_reads
filenames_forward_reads <- file.path(path, filenames_forward_reads)

# Plots the quality profiles of all eighteen samples
# My data only had 18 files and thus samples
plotQualityProfile(filenames_forward_reads[1:18])

# Place filtered files in filtered/ subdirectory
# Note this will fail if the directory doesn't exist
filter_path <- file.path("output", "filtered")
filtered_reads_path <- file.path(filter_path,
                                 paste0(sample_names,
                                        "_filt.fastq.gz"))

# See ?filterAndTrim for details on the parameters
# See here for adjustments for 454 data:
# https://benjjneb.github.io/dada2/
#     faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data
filtered_output <- filterAndTrim(fwd = filenames_forward_reads,
                                 filt = filtered_reads_path,
                                 maxLen = 600, # trim the sequence to 600
                                 maxN = 0, # discard any seqs with Ns
                                 maxEE = 3, # allow w/ up to 3 expected errors
                                 truncQ = 2, # cut off if quality gets this low
                                 rm.phix = TRUE,
                                 compress = TRUE,
                                 multithread = FALSE)

# Produce nicely-formatted markdown table of read counts
# Before/after trimming
kable(filtered_output,
      col.names = c("Reads In",
                    "Reads Out"))

# This build error models from each of the samples
errors_forward_reads <- learnErrors(filtered_reads_path,
                                    multithread = FALSE)

# Quick check to see if error models match data
# (black lines match black points) and are generally decresing left to right
plotErrors(errors_forward_reads,
           nominalQ = TRUE)

# Get rid of any duplicated sequences
dereplicated_forward_reads <- derepFastq(filtered_reads_path,
                                         verbose = TRUE)

# Name the derep-class objects by the sample names
names(dereplicated_forward_reads) <- sample_names

# Parameters adjusted based on recommendations for 454 data here:
# https://benjjneb.github.io/dada2/
#     faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data
dada_forward_reads <- dada(dereplicated_forward_reads,
                           err = errors_forward_reads,
                           HOMOPOLYMER_GAP_PENALTY = -1, # reduce penalty bc 454
                           BAND_SIZE = 32) # performs local alignments bc indels
# Check dada results
dada_forward_reads

# Produce the 'site by species matrix'
sequence_table <- makeSequenceTable(dada_forward_reads)

# Quick check to look at distribution of trimmed and denoised sequences
hist(nchar(getSequences(sequence_table)),
     main = "Histogram of final sequence variant lengths",
     xlab = "Sequence length in bp")

# Check for and remove chimeras
sequence_table_nochim <- removeBimeraDenovo(sequence_table,
                                            method = "consensus",
                                            multithread = FALSE,
                                            verbose = TRUE)

# What percent of our reads are non-chimeric?
non_chimeric_reads <- round(sum(sequence_table_nochim) / sum(sequence_table),
                            digits = 4) * 100

# Build a table showing how many sequences remain at each step of the pipeline
get_n <- function(x) sum(getUniques(x)) # make a quick function
track <- cbind(filtered_output, # already has 2 columns
               sapply(dada_forward_reads, get_n),
               rowSums(sequence_table),
               rowSums(sequence_table_nochim))

# Add nice meaningful column names
colnames(track) <- c("Input",
                     "Filtered",
                     "Denoised",
                     "Sequence Table",
                     "Non-chimeric")

# Set the proper rownames
rownames(track) <- sample_names

# Produce nice markdown table of progress through the pipeline
kable(track)

# Assigns taxonomy to each sequence variant based on a supplied training set
# Made up of known sequences
taxa <- assignTaxonomy(sequence_table_nochim,
                       "data/training/rdp_train_set_16.fa.gz",
                       multithread = FALSE,
                       tryRC = TRUE) # also check with seq reverse compliments

# Show the results of the taxonomy assignment
unname(taxa)

# We want to export the cleaned, trimmed, filtered, denoised sequence variants
# So that we can build a phylogeny - we'll build the phylogeny outside of R
# But we need the fasta file to do so. We keep the names of each sequence as the
# Sequence itself (which is rather confusing), because that's how DADA2 labels
# It's columns (e.g. 'species')
# Function taken from https://github.com/benjjneb/dada2/issues/88
export_taxa_table_and_seqs <- function(sequence_table_nochim,
                                       file_seqtab,
                                       file_seqs) {
  seqtab_t <- as.data.frame(t(sequence_table_nochim)) # transpose to data frame
  seqs <- row.names(seqtab_t) # extract rownames
  row.names(seqtab_t) <- seqs # set rownames to sequences
  outlist <- list(data_loaded = seqtab_t)
  mctoolsr::export_taxa_table(outlist, file_seqtab) # write out an OTU table
  seqs <- as.list(seqs)
  seqinr::write.fasta(seqs, row.names(seqtab_t), file_seqs) # write out fasta
}

# Actually run the function, with the names of the files we want it to create
# And where to put them
export_taxa_table_and_seqs(sequence_table_nochim,
                           "output/sequence_variants_table.txt",
                           "output/sequence_variants_seqs.fa")

# Next we want to read in the metadata file so we can add that in too
# This is not a csv file, so we have to use a slightly different syntax
# Here the `sep = "\t"` tells the function that the data are tab-delimited
# And the `stringsAsFactors = FALSE` tells it not to assume that things are
# Categorical variables
metadata_in <- read.table(paste0("data/metadata/",
                                 "SraRunTable.txt"),
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = FALSE,
                          row.names = 7) # sets sample IDs to row names
# The run or 7th column had the ERR numbers for the sequences

# Construct phyloseq object (straightforward from dada2 outputs)
phyloseq_obj <- phyloseq(otu_table(sequence_table_nochim,
                                   taxa_are_rows = FALSE), # sample-spp matrix
                         sample_data(metadata_in), # metadata for each sample
                         tax_table(taxa)) # taxonomy for each sequence variant

# Save phyloseq object to Rdata file
save(phyloseq_obj, file = "output/phyloseq_obj.Rdata")

# Melt phyloseq object
melted_phyloseq <- psmelt(phyloseq_obj)

# save melted phyloseq object
save(melted_phyloseq, file = "output/melted_phyloseq.Rdata")
