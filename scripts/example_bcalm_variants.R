# Argument parser
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = "Process BCALM variant data")
parser$add_argument("--counts", type = "character", required = TRUE, help = "Path to the counts file")
parser$add_argument("--map", type = "character", required = TRUE, help = "Path to the map file")
parser$add_argument("--output", type = "character", required = TRUE, help = "Path to the output file")
parser$add_argument("--output-plot", type = "character", required = FALSE, help = "Path to the output file")
parser$add_argument("--normalize", type = "logical", default = TRUE, help = "Whether to normalize the data (TRUE or FALSE)")

args <- parser$parse_args()

# Load the required libraries
suppressPackageStartupMessages(library(BCalm))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))


# read in the data
message("Reading input files...")

message("counts file: ", args$counts)
counts_df <- read.table(args$counts, header = TRUE, sep = "\t", fill = TRUE, c("", "NA", "N/A"))
colnames(counts_df)[1:2] <- c("Barcode", "name")

message("map file: ", args$map)
variant_map <- read.table(args$map, header = TRUE)

var_df <- create_var_df(counts_df, variant_map)

dna_var <- create_dna_df(var_df)
rna_var <- create_rna_df(var_df)

# create the MPRASet object
message("Creating MPRASet object...")
mpraset <- MPRASet(DNA = dna_var, RNA = rna_var, eid = row.names(dna_var), barcode = NULL)

message("Fit the model...")
nr_reps <- as.integer((ncol(counts_df) - 2) / 2)
bcs <- ncol(dna_var) / nr_reps
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mpraset)))
block_vector <- rep(1:nr_reps, each = bcs)
mpralm_fit_var <- mpralm(
    object = mpraset, design = design, aggregate = "none",
    normalize = args$normalize, model_type = "corr_groups", plot = FALSE, block = block_vector
)

message("Get variant statistics...")
mpra_variants <- topTable(mpralm_fit_var, coef = 2, number = Inf, confint = TRUE)

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
