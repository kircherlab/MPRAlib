library(BCalm)
library(dplyr)
library(ggplot2)
library(tidyr)

# read in the data
COUNTS <- read.table("test/test_bc_variant_counts.tsv.gz", header = T)
colnames(COUNTS) <- c("Barcode", "name", "dna_count_1", "rna_count_1", "dna_count_2", "rna_count_2", "dna_count_3", "rna_count_3")

MAP <- read.table("test/test_variant_map_fix.tsv.gz", header = T)

var_df <- create_var_df(COUNTS, MAP)

dna_var <- create_dna_df(var_df)
rna_var <- create_rna_df(var_df)

# create the MPRASet object
mpraset <- MPRASet(DNA = dna_var, RNA = rna_var, eid = row.names(dna_var), barcode = NULL)

nr_reps <- 3
bcs <- ncol(dna_var) / nr_reps
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mpraset)))
block_vector <- rep(1:nr_reps, each = bcs)
mpralm_fit_var <- mpralm(object = mpraset, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = FALSE, block = block_vector)

top_var <- topTable(mpralm_fit_var, coef = 2, number = Inf)

plot(top_var$logFC, -log10(top_var$adj.P.Val),
    pch = ".", cex = 3,
    xlab = "log2 fold change",
    ylab = "-log10(p-value)"
)

abline(2, 0, col = "red", lty = 2)

idx <- top_var$adj.P.Val < 0.01

points(top_var$logFC[idx], -log10(top_var$adj.P.Val[idx]), col = "red", pch = ".", cex = 3)


names <- c("ID", colnames(toptab))
toptab$ID <- rownames(toptab)
toptab <- toptab[, names]

write.table(toptab, "test/test_bc_variant_bcalm.tsv.gz", row.names = FALSE, sep = "\t", quote = FALSE)
