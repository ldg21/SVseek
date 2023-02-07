test_that("Fusion VAF calculation", {
    current <- calculateVAF(data_fusion$gr, data_fusion$breakpoint, "qiaseq")
    target <- data_fusion$vaf
    expect_equal(current, target)
})

test_that("Deletion VAF calculation", {
    current <- calculateVAF(data_deletion$gr, data_deletion$breakpoint, "qiaseq")
    target <- data_deletion$vaf
    expect_equal(current, target)
})

test_that("Predefined supporting and nonsupporting reads", {
    current <- calculateVAF(data_fusion$gr, data_fusion$breakpoint,
                            "qiaseq", supp = data_fusion$supporting_reads,
                            nonsupp = data_fusion$nonsupporting_reads)
    target <- data_fusion$vaf
    expect_equal(current, target)
})
