test_that("Fusion weighted average calculation", {
    current <- calculateWeightedAverage(data_fusion$vaf)
    target <- data_fusion$weighted_vaf
    expect_equal(current, target)
})

test_that("Deletion weighted average calculation", {
    current <- calculateWeightedAverage(data_deletion$vaf)
    target <- data_deletion$weighted_vaf
    expect_equal(current, target)
})
