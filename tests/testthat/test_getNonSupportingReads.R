test_that("Deletion nonsupporting reads", {
    current <- getNonSupportingReads(data_deletion$gr, data_deletion$breakpoint, "qiaseq")
    target <- data_deletion$nonsupporting_reads
    expect_equal(current, target)
})

test_that("Fusion nonsupporting reads", {
    current <- getNonSupportingReads(data_fusion$gr, data_fusion$breakpoint, "qiaseq")
    target <- data_fusion$nonsupporting_reads
    expect_equal(current, target)
})
