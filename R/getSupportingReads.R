#' Get supporting reads
#'
#' This function is a wrapper for selectReads which sets a default value for
#' overlap sizes if not specified.
#'
#' @param gr GRanges with reads
#' @param breakpoint breakpoint
#' @param technology Technology of reads ("qiaseq" | "other")
#' @param min_overlap minimum overlap for each side of breakend (default = 19)
#' @param min_overlap_flanking minimum overlap for flanking region
#' (set to min_overlap if not specified)
#' @param min_overlap_variant minimum overlap for variant region
#' (set to min_overlap if not specified)
#' @param split_only get split supporting reads only
#'
#' @return List of supporting reads by breakpoint
#'
#' @examples getSupportingReads(data_fusion$gr, data_fusion$breakpoint, "qiaseq")
#'
#' @export

getSupportingReads <- function(gr, breakpoint, technology,
                               min_overlap = 19,
                               min_overlap_flanking = min_overlap,
                               min_overlap_variant = min_overlap,
                               split_only = FALSE) {

    selectReads(gr, breakpoint, technology, TRUE,
                min_overlap_flanking, min_overlap_variant, split_only)

}
