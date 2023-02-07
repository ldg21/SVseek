test_that("Fusion bam loading", {
    path <- system.file("extdata", package = "SVseek")
    path <- file.path(path, "bams", "data_fusion_hg19.bam")
    skip_if_not(file.exists(path))
    current <- loadBAM(path, data_fusion$genes)
    names(current) <- as.integer(factor(names(current), levels = unique(names(current))))
    current$qname <- names(current)
    pattern <- "I|D|N|P|X" ## only allow M/S/H/=
    exclude <- names(current)[grep(pattern, current$cigar)]
    current <- current[!names(current) %in% exclude]
    target <- data_fusion$gr
    expect_equal(current, target)
})

test_that("Deletion bam loading", {
    path <- system.file("extdata", package = "SVseek")
    path <- file.path(path, "bams", "data_deletion_hg38.bam")
    skip_if_not(file.exists(path))
    current <- loadBAM(path, data_deletion$genes)
    names(current) <- as.integer(factor(names(current), levels = unique(names(current))))
    current$qname <- names(current)
    pattern <- "I|D|N|P|X" ## only allow M/S/H/=
    exclude <- names(current)[grep(pattern, current$cigar)]
    current <- current[!names(current) %in% exclude]
    target <- data_deletion$gr
    expect_equal(current, target)
})

test_that("Test partial bounds", {
    path <- system.file("extdata", package = "SVseek")
    path <- file.path(path, "bams", "data_fusion_hg19.bam")
    skip_if_not(file.exists(path))
    expect_length(loadBAM(path, GRanges("chr11:128665292-128665418")), 3)
})

test_that("Test empty bounds", {
    path <- system.file("extdata", package = "SVseek")
    path <- file.path(path, "bams", "data_fusion_hg19.bam")
    skip_if_not(file.exists(path))
    expect_length(loadBAM(path, GRanges("chr1:1-10")), 0)
})

test_that("Test incorrect inputs", {
    expect_error(loadBAM("notafile.bam"))
    expect_error(loadBAM(""))
})
