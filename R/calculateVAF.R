#' Calculate variant allele fraction
#'
#' Calculate variant allele fraction by counting supporting and nonsupporting
#' reads for each structural variant.
#'
#' The breakpoints argument can be a GRanges with a single variant, defined by
#' two breakends, or GRangesList with multiple variants.
#'
#' @param gr GRanges with reads
#' @param breakpoints GRanges of 2 breakends representing a single variant or
#' GRangesList of multiple variants
#' @param technology Technology of reads ("qiaseq" | "other")
#' @param min_overlap minimum overlap for each side of breakend (default = 19)
#' @param min_overlap_flanking minimum overlap for flanking region
#' (set to min_overlap if not specified)
#' @param min_overlap_variant minimum overlap for variant region
#' (set to min_overlap if not specified)
#' @param supp supporting reads (optional)
#' @param nonsupp nonsupporting reads (optional)
#'
#' @return breakpoints annotated with number of supporting reads,
#' nonsupporting reads, and VAF
#'
#' @examples calculateVAF(data_fusion$gr, data_fusion$breakpoint, "qiaseq")
#'
#' @export

calculateVAF <- function(gr, breakpoints, technology, min_overlap = 19,
                         min_overlap_flanking = min_overlap,
                         min_overlap_variant = min_overlap,
                         supp = NULL, nonsupp = NULL) {

    validate_gr(gr)
    validate_breakpoints(breakpoints)

    if (is(breakpoints, "GRangesList")) {

        out <- endoapply(breakpoints, function(bp) {
            calculateVAF(gr, bp, technology, min_overlap,
                         min_overlap_flanking, min_overlap_variant,
                         supp, nonsupp)
        })
        mcols(out)$w_avg <- calculateWeightedAverage(out)
        return(out)

    }

    if (is.null(nonsupp)) {

        nonsupp <- getNonSupportingReads(gr, breakpoints, technology,
                                         min_overlap, min_overlap_flanking,
                                         min_overlap_variant)
    }
    if (is.null(supp)) {

        supp <- getSupportingReads(gr, breakpoints, technology,
                                   min_overlap, min_overlap_flanking,
                                   min_overlap_variant, FALSE)
    }

    supp_split <- getSupportingReads(gr, breakpoints, technology,
                                        min_overlap, min_overlap_flanking,
                                        min_overlap_variant, TRUE)

    breakpoints$nonsupp <- lengths(nonsupp)
    breakpoints$supp <- lengths(supp)
    breakpoints$supp_split <- lengths(supp_split)
    breakpoints$total <- breakpoints$supp + breakpoints$nonsupp
    breakpoints$vaf <- breakpoints$supp / breakpoints$total

    return(breakpoints)

}
