#' Extract non-split reads
#'
#' @param gr GRanges with reads
#'
#' @return List with R1, R2

extractNonSplitReads <- function(gr)
{

    validate_gr(gr)

    reads <- separateReads(gr)
    R1 <- reads$R1
    R2 <- reads$R2

    ## select non-split reads
    R1 <- R1[lengths(R1) == 1L]
    R2 <- R2[lengths(R2) == 1L]

    return(list(R1 = R1, R2 = R2))

}

#' Extract split reads
#'
#' @param gr GRanges with reads
#'
#' @return List with R1, R2

extractSplitReads <- function(gr)
{

    validate_gr(gr)

    gr_strand <- as.character(strand(gr))
    right_clipped_regex <- "^\\d+M\\d+(S|H)$"
    left_clipped_regex <- "^\\d+(S|H)\\d+M$"

    mcols(gr)$seg <- NA_integer_
    i <- which(gr_strand == "+" & grepl(right_clipped_regex, mcols(gr)$cigar))
    mcols(gr)$seg[i] <- 1L
    i <- which(gr_strand == "+" & grepl(left_clipped_regex, mcols(gr)$cigar))
    mcols(gr)$seg[i] <- 2L
    i <- which(gr_strand == "-" & grepl(right_clipped_regex, mcols(gr)$cigar))
    mcols(gr)$seg[i] <- 2L
    i <- which(gr_strand == "-" & grepl(left_clipped_regex, mcols(gr)$cigar))
    mcols(gr)$seg[i] <- 1L

    reads <- separateReads(gr)
    R1 <- reads$R1
    R2 <- reads$R2

    ## select split reads
    R1 <- R1[lengths(R1) == 2L]
    tmp <- unlist(R1)
    names(tmp) <- mcols(tmp)$qname
    valid <- intersect(
        names(tmp)[which(mcols(tmp)$seg == 1)],
        names(tmp)[which(mcols(tmp)$seg == 2)])
    R1 <- R1[names(R1) %in% valid]
    R2 <- R2[lengths(R2) == 2L]
    tmp <- unlist(R2)
    names(tmp) <- mcols(tmp)$qname
    valid <- intersect(
        names(tmp)[which(mcols(tmp)$seg == 1)],
        names(tmp)[which(mcols(tmp)$seg == 2)])
    R2 <- R2[names(R2) %in% valid]

    return(list(R1 = R1, R2 = R2))

}

#' Separate reads
#'
#' @param gr GRanges with reads
#'
#' @return List with R1, R2

separateReads <- function(gr) {

    validate_gr(gr)

    tmp <- gr[gr$R1]
    R1 <- split(tmp, names(tmp))
    tmp <- gr[gr$R2]
    R2 <- split(tmp, names(tmp))

    return(list(R1 = R1, R2 = R2))

}

#' Get top breakpoints from countBreakpoints GRangesList
#'
#' @param breakpoints GRangesList with breakpoints
#' @param top_n Top n breakpoints to return
#' @param min_count Minimum number of required reads
#'
#' @return GRangesList with top breakpoints

getTopBreakpoints <- function(breakpoints, top_n, min_count) {

    validate_breakpoints(breakpoints)

    if (length(breakpoints) == 0L) return(breakpoints)

    breakpoints <- countBreakpoints(breakpoints)
    breakpoints <- breakpoints[order(mcols(breakpoints)$count,
                                     decreasing = TRUE)]
    breakpoints <- breakpoints[mcols(breakpoints)$count >= min_count]
    if (is.finite(top_n)) breakpoints <- head(breakpoints, top_n)
    names(breakpoints) <- NULL
    mcols(breakpoints) <- NULL

    return(breakpoints)

}

#' Count unique breakpoints in breakpoints GRangesList
#'
#' Collapses breakpoints and adds metadata column with count for each
#' breakpoint
#'
#' @param breakpoints GRangesList
#'
#' @return GRangesList of breakpoints annotated with counts

countBreakpoints <- function(breakpoints) {

    validate_breakpoints(breakpoints)

    seqlevels_temp <- seqlevels(breakpoints)
    seqlengths_temp <- seqlengths(breakpoints)
    bp_n <- tapply(
        names(breakpoints),
        sapply(breakpoints, function(bp) paste(bp[1], bp[2])),
        function(x) length(unique(x)))
    bps <- strsplit(names(bp_n), " ")
    bp_1 <- as(GRanges(sapply(bps, "[", 1)), "GRangesList")
    bp_2 <- as(GRanges(sapply(bps, "[", 2)), "GRangesList")
    breakpoints <- suppressWarnings(pc(bp_1, bp_2))
    mcols(breakpoints)$count <- as.integer(bp_n)
    seqlevels(breakpoints) <- seqlevels_temp
    seqlengths(breakpoints) <- seqlengths_temp

    return(breakpoints)

}
