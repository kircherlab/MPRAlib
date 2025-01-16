library(mpra)
library(dplyr)
library(ggplot2)
library(tidyr)


# read in the data
DNA <- read.table("test_dna.tsv.gz", header = T, row.names = "ID")
RNA <- read.table("test_rna.tsv.gz", header = T, row.names = "ID")


# create the MPRASet object
mpraset <- MPRASet(
    DNA = DNA,
    RNA = RNA,
    eid = rownames(DNA),
    eseq = NULL,
    barcode = NULL
)

# create the design matrix
design <- model.matrix(~1, data = data.frame(sample = 1:3))


# run the mpralm analysis
fit_element <- mpralm(
    object = mpraset,
    design = design,
    aggregate = "none",
    normalize = TRUE,
    model_type = "indep_groups",
    plot = TRUE
)

toptab_element <- topTable(fit_element, coef = 1, number = Inf)

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
tr <- treat(fit_element, lfc = thr)
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
mpra_element <- toptab[, names]

write.table(mpra_element, "test_mpralm_out.tsv.gz", row.names = FALSE, sep = "\t", quote = FALSE)
