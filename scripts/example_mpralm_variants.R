# Argument parser
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = "Process BCALM variant data")
parser$add_argument("--counts", type = "character", required = TRUE, help = "Path to the counts file")
parser$add_argument("--output", type = "character", required = TRUE, help = "Path to the output file")
parser$add_argument("--output-plot", type = "character", required = FALSE, help = "Path to the output file")

args <- parser$parse_args()

suppressPackageStartupMessages(library(mpra))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

# read in the data
counts_df <- read.table(args$counts, header = TRUE)
colnames(counts_df)[1] <- c("ID")
counts_df <- counts_df %>% column_to_rownames(var = "ID")

dna_var <- counts_df[, grepl("dna", colnames(counts_df))]
colnames(dna_var) <- gsub("dna_", "", colnames(dna_var))
rna_var <- counts_df[, grepl("rna", colnames(counts_df))]
colnames(rna_var) <- gsub("rna_", "", colnames(rna_var))

# create the MPRASet object
mpraset <- MPRASet(
    DNA = dna_var,
    RNA = rna_var,
    eid = rownames(dna_var),
    eseq = NULL,
    barcode = NULL
)

# create the design matrix
design <- data.frame(
    intcpt = 1,
    alleleB <- grepl("ALT", colnames(dna_var))
)
# create the block vector
block_vector <- rep(1:(ncol(dna_var) / 2), each = 2)

# run the mpralm analysis
mpralm_allele_fit <- mpralm(
    object = mpraset,
    design = design,
    aggregate = "none",
    normalize = TRUE,
    block = block_vector,
    model_type = "corr_groups",
    plot <- FALSE
)

mpra_variants <- topTable(mpralm_allele_fit, coef = 2, number = Inf, confint = TRUE)


if (!is.null(args$output_plot)) {
    p <- ggplot(mpra_variants, aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point(alpha = 0.5) +
        geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
        geom_point(data = subset(mpra_variants, adj.P.Val < 0.01), aes(x = logFC, y = -log10(adj.P.Val)), color = "red") +
        labs(x = "log2 fold change", y = "-log10(p-value)") +
        theme_minimal()

    ggsave(filename = args$output_plot, plot = p, width = 8, height = 6)
}


names <- c("ID", colnames(mpra_variants))
mpra_variants$ID <- rownames(mpra_variants)
mpra_variants <- mpra_variants[, names]

gzfile_output <- gzfile(args$output, "w")
write.table(mpra_variants, gzfile_output, row.names = FALSE, sep = "\t", quote = FALSE)
close(gzfile_output)
