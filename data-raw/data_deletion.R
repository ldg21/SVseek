library(SVseek)

path <- system.file("extdata", package = "SVseek")
path <- file.path(path, "bams", "data_deletion_hg38.bam")

genes <- GRanges("chr3:41199505-41240443:+")

gr <- loadBAM(path, genes)

names(gr) <- as.integer(factor(names(gr), levels = unique(names(gr))))
gr$qname <- names(gr)

pattern <- "I|D|N|P|X" ## only allow M/S/H/=
exclude <- names(gr)[grep(pattern, gr$cigar)]
gr <- gr[!names(gr) %in% exclude]

breakpoint <- detectBreakpoints(gr, genes)[[1]]
vaf <- calculateVAF(gr, breakpoint, "qiaseq")
weighted_vaf <- calculateWeightedAverage(vaf)

supporting_reads <- getSupportingReads(gr, breakpoint, "qiaseq")
nonsupporting_reads <- getNonSupportingReads(gr, breakpoint, "qiaseq")

data_deletion <- list(
    gr = gr,
    breakpoint = breakpoint,
    genes = genes,
    vaf = vaf,
    supporting_reads = supporting_reads,
    nonsupporting_reads = nonsupporting_reads,
    weighted_vaf = weighted_vaf
)

usethis::use_data(data_deletion, overwrite = TRUE)
