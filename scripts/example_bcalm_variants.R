# Argument parser
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = "Process BCALM variant data")
parser$add_argument("--counts", type = "character", required = TRUE, help = "Path to the counts file")
parser$add_argument("--map", type = "character", required = TRUE, help = "Path to the map file")
parser$add_argument("--output", type = "character", required = TRUE, help = "Path to the output file")
parser$add_argument("--output-plot", type = "character", required = FALSE, help = "Path to the output file")

args <- parser$parse_args()

# Load the required libraries
suppressPackageStartupMessages(library(BCalm))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))


# read in the data
COUNTS <- read.table(args$counts, header = TRUE)
# TODO change this colnames to be more general
colnames(COUNTS) <- c("Barcode", "name", "dna_count_1", "rna_count_1", "dna_count_2", "rna_count_2", "dna_count_3", "rna_count_3")

MAP <- read.table(args$map, header = TRUE)

var_df <- create_var_df(COUNTS, MAP)

dna_var <- create_dna_df(var_df)
rna_var <- create_rna_df(var_df)

# create the MPRASet object
mpraset <- MPRASet(DNA = dna_var, RNA = rna_var, eid = row.names(dna_var), barcode = NULL)

# TODO make number of replicates also general
nr_reps <- 3
bcs <- ncol(dna_var) / nr_reps
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mpraset)))
block_vector <- rep(1:nr_reps, each = bcs)
mpralm_fit_var <- mpralm(object = mpraset, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = FALSE, block = block_vector)

top_var <- topTable(mpralm_fit_var, coef = 2, number = Inf)

if (!is.null(args$output_plot)) {
    png(filename = args$output_plot, width = 800, height = 600)

    plot(top_var$logFC, -log10(top_var$adj.P.Val),
        pch = ".", cex = 3,
        xlab = "log2 fold change",
        ylab = "-log10(p-value)"
    )

    abline(2, 0, col = "red", lty = 2)

    idx <- top_var$adj.P.Val < 0.01

    points(top_var$logFC[idx], -log10(top_var$adj.P.Val[idx]), col = "red", pch = ".", cex = 3)

    dev.off()
}

names <- c("ID", colnames(top_var))
top_var$ID <- rownames(top_var)
top_var <- top_var[, names]

write.table(top_var, args$output, row.names = FALSE, sep = "\t", quote = FALSE)
