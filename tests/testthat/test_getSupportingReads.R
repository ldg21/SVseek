test_that("Deletion supporting reads", {
    current <- getSupportingReads(data_deletion$gr, data_deletion$breakpoint, "qiaseq")
    target <- data_deletion$supporting_reads
    expect_equal(current, target)
})

test_that("Fusion supporting reads", {
    current <- getSupportingReads(data_fusion$gr, data_fusion$breakpoint, "qiaseq")
    target <- data_fusion$supporting_reads
    expect_equal(current, target)
})
