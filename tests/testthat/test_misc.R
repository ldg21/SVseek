test_that("Fusion supporting reads disjoint to nonsupporting reads", {
    expect_true(!all(unlist(data_fusion$supporting_reads) %in% unlist(data_fusion$nonsupporting_reads)))
})
