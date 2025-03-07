library(BCalm)
library(dplyr)
library(ggplot2)
library(tidyr)

# read in the data
COUNTS <- read.table("test/test_bc_element_counts.tsv.gz", header = T)

dna_elem <- create_dna_df(COUNTS, id_column_name = "oligos")
rna_elem <- create_rna_df(COUNTS, id_column_name = "oligos")

# BcLabelMPRASetExample <- MPRASet(DNA = dna_elem, RNA = rna_elem, eid = row.names(dna_elem), barcode = NULL, label=LabelExample))

# create the MPRASet object
mpraset <- MPRASet(DNA = dna_elem, RNA = rna_elem, eid = row.names(dna_elem), barcode = NULL)

nr_reps <- 3
bcs <- ncol(dna_elem) / nr_reps
block_vector <- rep(1:nr_reps, each = bcs)

mpralm_fit_elem <- fit_elements(object = mpraset, normalize = TRUE, block = block_vector, plot = FALSE)


toptab_element <- topTable(mpralm_fit_elem, coef = 1, number = Inf)
toptab_element %>% head()


pos_ctrl <- toptab_element[grep("^Positive_", rownames(toptab_element)), ]
neg_ctrl <- toptab_element[grep("^Negative_", rownames(toptab_element)), ]

test <- toptab_element[!(rownames(toptab_element) %in% c(rownames(pos_ctrl), rownames(neg_ctrl))), ]

hist(toptab_element$logFC,
    freq = F,
    ylim = c(0, 1.7), xlab = "log2FC", main = ""
)
lines(density(neg_ctrl$logFC), col = "red", lwd = 1.5)
lines(density(test$logFC), col = "blue", lwd = 1.5)
lines(density(pos_ctrl$logFC), col = "green", lwd = 1.5)


thr <- quantile(neg_ctrl$logFC, 0.9)

# Re-evaluate
tr <- treat(mpralm_fit_elem, lfc = thr)
mpra_element <- topTreat(tr, coef = 1, number = Inf)

# Make volcano plot with cutoff of FDR < 0.01
plot(mpra_element$logFC, -log10(mpra_element$adj.P.Val),
    pch = ".", cex = 3,
    xlab = "log2 fold change",
    ylab = "-log10(p-value)"
)
abline(2, 0, col = "red", lty = 2)
idx <- mpra_element$adj.P.Val < 0.01
points(mpra_element$logFC[idx], -log10(mpra_element$adj.P.Val[idx]),
    col = "red",
    pch = ".", cex = 3
)


names <- c("ID", colnames(mpra_element))
mpra_element$ID <- rownames(mpra_element)
mpra_element <- mpra_element[, names]

write.table(mpra_element, "test_mpralm_out.tsv.gz", row.names = FALSE, sep = "\t", quote = FALSE)
