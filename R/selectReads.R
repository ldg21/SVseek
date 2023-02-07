#' Select supporting or nonsupporting reads
#'
#' Get the names of reads that are supporting or nonsupporting
#' for a query breakpoint.
#'
#' The function returns either supporting or nonsupporting reads,
#' specified by the supporting boolean.
#'
#' TODO: Refactor function to return both supporting and nonsupporting at once.
#'
#' For technology = qiaseq, only reads in the R12 direction are returned.
#' R1 is adjacent to the gene specific primer, thus R1 and R2 map to the
#' flanking and variant region, respectively.
#'
#' A value for the minimum overlap in flanking and variant regions is required.
#'
#' @param gr GRanges with reads
#' @param breakpoint breakpoint
#' @param technology Technology of reads ("qiaseq" | "other")
#' @param supporting TRUE for supporting reads, FALSE for nonsupporting reads
#' @param min_overlap_flanking minimum overlap for flanking region
#' @param min_overlap_variant minimum overlap for variant region
#' @param split_only select split supporting reads only (ignored for
#' supporting = FALSE)
#'
#' @return List of supporting or nonsupporting reads by breakpoint
#'
#' @examples selectReads(data_fusion$gr, data_fusion$breakpoint, "qiaseq", TRUE,
#' 19, 19)
#'
#' @export

selectReads <- function(gr, breakpoint, technology, supporting,
                        min_overlap_flanking, min_overlap_variant,
                        split_only = FALSE) {

    validate_gr(gr)
    validate_breakpoint(breakpoint)

    bp_1 <- breakpoint[1]
    bp_2 <- breakpoint[2]

    max_overlap <- 1e3

    ## consider 1Mb upstream or downstream of breakpoint 1
    bp_1_upstream <- flank(bp_1, max_overlap, start = TRUE)
    bp_1_upstream <- reduce(c(bp_1_upstream, bp_1))
    bp_1_dnstream <- flank(bp_1, max_overlap, start = FALSE)
    bp_1_dnstream_rev <- invertStrand(bp_1_dnstream)

    ## consider 1Mb upstream or downstream of breakpoint 2
    bp_2_upstream <- flank(bp_2, max_overlap, start = TRUE)
    bp_2_upstream <- reduce(c(bp_2_upstream, bp_2))
    bp_2_dnstream <- flank(bp_2, max_overlap, start = FALSE)
    bp_2_dnstream_rev <- invertStrand(bp_2_dnstream)

    ## if alternative regions overlap (e.g. CTNNB1 deletion),
    ## restrict to unambiguous regions
    bp_1_dnstream_rev <- setdiff(bp_1_dnstream_rev, bp_2_upstream)
    ## w <- width(bp_1_dnstream_rev)
    ## if (w < max_overlap)
    ##     warning(paste("bp 1 downstream region restricted to", w, "nt"))
    bp_2_dnstream_rev <- setdiff(bp_2_dnstream_rev, bp_1_upstream)
    ## w <- width(bp_2_dnstream_rev)
    ## if (w < max_overlap)
    ##    warning(paste("bp 2 downstream region restricted to", w, "nt"))

    if (supporting) {

        bp_1_flanking <- bp_1_upstream
        bp_1_variant <- bp_2_upstream
        bp_2_flanking <- bp_2_upstream
        bp_2_variant <- bp_1_upstream

    } else {

        bp_1_flanking <- bp_1_upstream
        bp_1_variant <- bp_1_dnstream_rev
        bp_2_flanking <- bp_2_upstream
        bp_2_variant <- bp_2_dnstream_rev

    }

    ## non-supporting and supporting reads must satisfy overlap requirements
    ## for non-supporting reads, both mates must be non-split
    ## for supporting reads, mates can be
    ## - non-split and fully contained in overlap region, or
    ## - split and exactly match breakpoint

    tmp <- extractNonSplitReads(gr)
    R1 <- tmp$R1
    R2 <- tmp$R2

    if (supporting) {

        reads <- separateReads(gr)
        genes <- c(bp_1, invertStrand(bp_2))
        tmp <- extractBreakpoints(gr, genes)
        i <- which(sapply(tmp$R1, paste, collapse = " ") == paste(bp_1, bp_2))
        supp_split_R1 <- names(tmp$R1)[i]
        R1 <- c(R1, reads$R1[supp_split_R1])
        i <- which(sapply(tmp$R2, paste, collapse = " ") == paste(bp_1, bp_2))
        supp_split_R2 <- names(tmp$R2)[i]
        R2 <- c(R2, reads$R2[supp_split_R2])

    }

    bp_1_R12 <- overlapsRegions(R1, R2, bp_1_flanking, bp_1_variant,
                                min_overlap_flanking, min_overlap_variant,
                                supporting, max_overlap)
    bp_1_R21 <- overlapsRegions(R2, R1, bp_1_flanking, bp_1_variant,
                                min_overlap_flanking, min_overlap_variant,
                                supporting, max_overlap)
    bp_2_R12 <- overlapsRegions(R1, R2, bp_2_flanking, bp_2_variant,
                                min_overlap_flanking, min_overlap_variant,
                                supporting, max_overlap)
    bp_2_R21 <- overlapsRegions(R2, R1, bp_2_flanking, bp_2_variant,
                                min_overlap_flanking, min_overlap_variant,
                                supporting, max_overlap)

    if (technology == "qiaseq") {

        reads_bp_1 <- bp_1_R12
        reads_bp_2 <- bp_2_R12

    } else {

        reads_bp_1 <- c(bp_1_R12, bp_1_R21)
        reads_bp_2 <- c(bp_2_R12, bp_2_R21)

    }

    if (supporting && split_only) {

        supp_split <- union(supp_split_R1, supp_split_R2)
        reads_bp_1 <- intersect(reads_bp_1, supp_split)
        reads_bp_2 <- intersect(reads_bp_2, supp_split)

    }

    return(list(breakpoint_1 = reads_bp_1, breakpoint_2 = reads_bp_2))

}
