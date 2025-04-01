suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = "Process BCALM variant data")
parser$add_argument("--counts", type = "character", required = TRUE, help = "Path to the counts file")
parser$add_argument("--labels", type = "character", required = TRUE, help = "Path to the labels file")
parser$add_argument("--test-label", type = "character", required = TRUE, help = "Name of the test group")
parser$add_argument("--control-label", type = "character", required = TRUE, help = "Name of the control group")
parser$add_argument("--percentile",
    type = "double", default = 0.975,
    help = "Percentile of control to test on. Default is 0.975"
)
parser$add_argument("--output", type = "character", required = TRUE, help = "Path to the output file")
parser$add_argument("--output-vulcano-plot", type = "character", required = FALSE, help = "Path to store the vulcano plot")
parser$add_argument("--output-density-plot", type = "character", required = FALSE, help = "Path to store the density plot")
parser$add_argument("--normalize", type = "logical", default = TRUE, help = "Whether to normalize the data (TRUE or FALSE)")

args <- parser$parse_args()

# Load the required libraries
suppressPackageStartupMessages(library(BCalm))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))


# read in the data
counts_df <- read.table(args$counts, header = TRUE, sep = "\t", fill = TRUE, c("", "NA", "N/A"))
colnames(counts_df)[1:2] <- c("Barcode", "name")

dna_elem <- create_dna_df(counts_df, id_column_name = "name")
rna_elem <- create_rna_df(counts_df, id_column_name = "name")

labels <- read.table(args$labels, header = FALSE, sep = "\t", col.names = c("name", "label"))

labels_vec <- as.vector(labels$label)
names(labels_vec) <- labels$name
# Use only these labels of the sequences that remained after filtering
labels_vec <- labels_vec[rownames(dna_elem)]


# create the MPRASet object
cat("Creating MPRASet object...\n")
mpraset <- MPRASet(DNA = dna_elem, RNA = rna_elem, eid = row.names(dna_elem), barcode = NULL, label = labels_vec)

nr_reps <- as.integer((ncol(counts_df) - 2) / 2)
bcs <- ncol(dna_elem) / nr_reps
block_vector <- rep(1:nr_reps, each = bcs)

cat("Fit elements...\n")
fit_elem <- fit_elements(object = mpraset, normalize = args$normalize, block = block_vector, plot = FALSE)

toptab_element <- topTable(fit_elem, coef = 1, number = Inf)
percentile <- args$percentile

if (!is.null(args$output_density_plot)) {
    cat("Plot density elements...\n")
    toptab_element_label <- toptab_element %>%
        rownames_to_column(var = "name") %>%
        left_join(labels, by = "name") %>%
        column_to_rownames(var = "name")

    percentile_up <- quantile(toptab_element_label$logFC[toptab_element_label$label == args$control_label], percentile)
    up_label <- paste(percentile, "th percentile of negative controls", sep = "")

    percentile_down <- quantile(toptab_element_label$logFC[toptab_element_label$label == args$control_label], 1 - percentile)
    down_label <- paste(1 - percentile, "th percentile of negative controls", sep = "")


    density_plot <- ggplot(toptab_element_label, aes(x = logFC, fill = label, y = after_stat(density))) +
        geom_histogram(alpha = 0.5, position = "identity", binwidth = 0.1) +
        geom_density(alpha = 0.2, adjust = 1) +
        labs(x = "log2 fold change", y = "Density") +
        xlim(c(min(toptab_element_label$logFC), max(toptab_element_label$logFC))) +
        geom_vline(aes(xintercept = percentile_up, color = up_label), linetype = "dashed", linewidth = 1) +
        geom_vline(aes(xintercept = percentile_down, color = down_label), linetype = "dashed", linewidth = 1) +
        scale_color_manual(
            values = setNames(c("green", "orange"), c(up_label, down_label)),
            guide = guide_legend(override.aes = list(linetype = "dashed"))
        ) +
        theme_minimal()

    ggsave(filename = args$output_density_plot, plot = density_plot, width = 8, height = 6)
}

# # Re-evaluate
tr <- mpra_treat(fit_elem, percentile, neg_label = args$control_label)
mpra_element <- topTreat(tr, coef = 1, number = Inf)

# Make volcano plot with cutoff of FDR < 0.01
if (!is.null(args$output_vulcano_plot)) {
    cat("Plot vulcano...\n")
    p <- ggplot(mpra_element, aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point(alpha = 0.5) +
        geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
        geom_point(data = subset(mpra_element, adj.P.Val < 0.01), aes(x = logFC, y = -log10(adj.P.Val)), color = "red") +
        labs(x = "log2 fold change", y = "-log10(p-value)") +
        theme_minimal()

    ggsave(filename = args$output_vulcano_plot, plot = p, width = 8, height = 6)
}


names <- c("ID", colnames(mpra_element))
mpra_element$ID <- rownames(mpra_element)
mpra_element <- mpra_element[, names]

cat("Write output to file...\n")
gzfile_output <- gzfile(args$output, "w")
write.table(mpra_element, gzfile_output, row.names = FALSE, sep = "\t", quote = FALSE)
close(gzfile_output)
