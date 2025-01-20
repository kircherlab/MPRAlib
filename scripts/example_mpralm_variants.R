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
design <- data.frame(
    intcpt = 1,
    alleleB = grepl("ALT", colnames(DNA))
)
# create the block vector
block_vector <- rep(1:3, each = 2)

# run the mpralm analysis
mpralm_allele_fit <- mpralm(
    object = mpraset,
    design = design,
    aggregate = "none",
    normalize = TRUE,
    block = block_vector,
    model_type = "corr_groups",
    plot = TRUE
)

toptab <- topTable(mpralm_allele_fit, coef = 2, number = Inf, confint = TRUE)

toptab %>% head()

analyze_var <- function(var, A1, A2) {
    dna_exp <- assay(norm_counts, 1)[var, ]
    rna_exp <- assay(norm_counts, 2)[var, ]
    n_rep <- length(dna_exp) / 2
    ratio <- as.numeric(log2(rna_exp / dna_exp))

    variant <- data.frame(ratio,
        allele = rep(c(A1, A2), n_rep),
        sample = as.character(rep(1:n_rep, each = 2))
    )
    variant$allele <- factor(variant$allele, levels = c(A1, A2))

    g <- variant %>%
        ggplot(aes(allele, ratio, color = sample)) +
        geom_point() +
        geom_line(aes(group = sample, color = sample)) +
        # geom_text(aes(label=ifelse(allele=="C",sample,'')),col="black",hjust=2,vjust=0) +
        scale_color_discrete(
            name = "Sample",
            breaks = c("1", "2", "3"),
            labels = c("Rep 1", "Rep 2", "Rep 3")
        ) +
        ggtitle(var) +
        ylab("log2(RNA/DNA)") +
        theme(
            text = element_text(size = 15),
            plot.title = element_text(hjust = 0.5),
            legend.title = element_text(hjust = 0.5)
        )
    return(g)
}
var <- "NC_000003.11:187930322:TT:T"
norm_counts <- normalize_counts(mpraset, block = block_vector)
analyze_var(var, "TT", "T")

plot(toptab$logFC, -log10(toptab$adj.P.Val),
    pch = ".", cex = 3,
    xlab = "log2 fold change",
    ylab = "-log10(p-value)"
)

abline(2, 0, col = "red", lty = 2)

idx <- toptab$adj.P.Val < 0.01

points(toptab$logFC[idx], -log10(toptab$adj.P.Val[idx]), col = "red", pch = ".", cex = 3)


names <- c("ID", colnames(toptab))
toptab$ID <- rownames(toptab)
toptab <- toptab[, names]

write.table(toptab, "test_mpralm_out.tsv.gz", row.names = FALSE, sep = "\t", quote = FALSE)
