#' Detect structural variants
#'
#' Find location of breakpoints from split reads.
#'
#' genes must be a GRanges with genomic coordinates of a single gene for
#' intragenic variants, or two genes for intergenic variants.
#'
#' top_n specifies the number of variants to return
#'
#' min_count filters variants by minimum number of supporting reads
#'
#' @param gr GRanges with reads
#' @param genes GRanges with one or two candidate genes
#' @param min_count Minimum number of supporting reads
#' @param top_n Top n variants to return
#' @param detect_reciprocal Detect reciprocal rearrangement for two candidate
#' genes (default = FALSE)
#'
#' @return GRangesList with detected variants, each defined by two breakends
#'
#' @examples detectBreakpoints(data_fusion$gr, data_fusion$genes)
#' @examples detectBreakpoints(data_fusion$gr, data_fusion$genes, 3, 1)
#'
#' @export

detectBreakpoints <- function(gr, genes, min_count = 1, top_n = Inf,
                              detect_reciprocal = FALSE) {

    validate_gr(gr)
    validate_genes(genes)

    if (length(genes) == 1) {

        detect_reciprocal <- FALSE
        genes <- c(genes[1], genes[1])

    }

    ## extract breakpoints from split reads
    bps <- extractBreakpoints(gr, genes)

    ## exclude reads with mates split at different breakpoints
    exclude <- intersect(names(bps$R1), names(bps$R2))
    tmp_R1 <- sapply(bps$R1[exclude], paste, collapse = " ")
    tmp_R2 <- sapply(bps$R2[exclude], paste, collapse = " ")
    exclude <- exclude[which(tmp_R1 != tmp_R2)]
    breakpoints <- c(
        bps$R1[!names(bps$R1) %in% exclude],
        bps$R2[!names(bps$R2) %in% exclude])

    if (detect_reciprocal) {

        breakpoints_reversed_reads <- detectBreakpoints(gr, rev(genes), top_n,
                                                    min_count, FALSE)
        breakpoints <- c(breakpoints, breakpoints_reversed_reads)

    }

    breakpoints <- getTopBreakpoints(breakpoints, top_n, min_count)

    return(breakpoints)

}
