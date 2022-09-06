#!/usr/bin/env Rscript

require(fastbaps)
require(ape)

# Load in arguments
args <- commandArgs(trailingOnly = TRUE)
aln_file <- args[1]
output <- args[2]
levels <- as.numeric(args[3])
tree_file <- args[4]
threads <- as.numeric(args[5])

# load FASTA
tryCatch(
    expr = {
        sparse_data <- fastbaps::import_fasta_sparse_nt(aln_file)
        sparse_data <- fastbaps::optimise_prior(sparse_data, type = "baps")
        tree <- ape::read.tree(tree_file)
        multi <- fastbaps::multi_level_best_baps_partition(sparse_data,
            h = tree, levels = levels, n.cores = threads)
    },
    error = function(e) {
        print(e)
        message("Problem loading fasta file – not generating subclusters")
        snp_data <- fastbaps:::import_fasta_to_vector_each_nt()
        names <-  gsub("^>", "", snp_data$seq.names)
        null_clusters <- array(1, c(len(names), levels))
        colnames(null_clusters) <- paste0("Level ", seq(1, levels))
        multi <- cbind(data.frame(Isolates = snp_data$seq.names), null_clusters)
    }
)

# write results
write.table(multi, file = output, row.names = FALSE, col.names = TRUE,
            sep = ",", quote = FALSE)
