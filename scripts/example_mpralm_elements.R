# Argument parser
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
parser$add_argument("--output-vulcano-plot", type = "character", required = FALSE, help = "Path to stroe the vulcano plot")
parser$add_argument("--output-density-plot", type = "character", required = FALSE, help = "Path to stroe the density plot")

args <- parser$parse_args()

suppressPackageStartupMessages(library(mpra))
suppressPackageStartupMessages(library(BCalm))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))


# read in the data
counts_df <- read.table(args$counts, header = TRUE)
colnames(counts_df)[1] <- c("ID")
counts_df <- counts_df %>% column_to_rownames(var = "ID")

dna_elem <- counts_df[, grepl("dna", colnames(counts_df))]
colnames(dna_elem) <- gsub("dna_", "", colnames(dna_elem))
rna_elem <- counts_df[, grepl("rna", colnames(counts_df))]
colnames(rna_elem) <- gsub("rna_", "", colnames(rna_elem))

labels <- read.table(args$labels, header = FALSE, sep = "\t", col.names = c("name", "label"))

labels_vec <- as.vector(labels$label)
names(labels_vec) <- labels$name
# Use only these labels of the sequences that remained after filtering
labels_vec <- labels_vec[rownames(dna_elem)]


# create the MPRASet object
mpraset <- MPRASet(
    DNA = dna_elem,
    RNA = rna_elem,
    eid = rownames(dna_elem),
    eseq = NULL,
    barcode = NULL,
)

# create the design matrix
design <- model.matrix(~1, data = data.frame(sample = 1:ncol(dna_elem)))


# run the mpralm analysis
fit_elem <- mpralm(
    object = mpraset,
    design = design,
    aggregate = "none",
    normalize = TRUE,
    model_type = "indep_groups",
    plot <- FALSE
)

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


# Re-evaluate
# tr <- treat(fit_elem, lfc = percentile_up)
fit_elem$label <- labels_vec
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
