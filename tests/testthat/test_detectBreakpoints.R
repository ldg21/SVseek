test_that("Fusion breakpoint detection", {
    current <- detectBreakpoints(data_fusion$gr, data_fusion$genes)[[1]]
    target <- data_fusion$breakpoint
    expect_equal(current, target)
})

test_that("Deletion breakpoint detection", {
    current <- detectBreakpoints(data_deletion$gr, data_deletion$genes)[[1]]
    target <- data_deletion$breakpoint
    expect_equal(current, target)
})
