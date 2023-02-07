#' Extract breakpoints from aligned reads
#'
#' @param gr GRanges with reads
#' @param genes GRanges with one or two candidate genes
#'
#' @return GRangesList with breakpoints

extractBreakpoints <- function(gr, genes)
{

    gene_1 <- genes[1]
    gene_2 <- genes[2]

    gene_1_inv <- invertStrand(gene_1)
    gene_2_inv <- invertStrand(gene_2)

    split_reads <- extractSplitReads(gr)

    out <- list()

    for (read in c("R1", "R2")) {

        reads <- split_reads[[read]]

        tmp <- unlist(reads)
        names(tmp) <- mcols(tmp)$qname
        seg_1 <- tmp[which(mcols(tmp)$seg == 1)]
        seg_1 <- seg_1[match(names(reads), names(seg_1))]
        seg_2 <- tmp[which(mcols(tmp)$seg == 2)]
        seg_2 <- seg_2[match(names(reads), names(seg_2))]

        ## consider split reads in forward orientation
        ## (overlapping both genes in sense)
        reads_fwd <- reads[seg_1 %over% gene_1 & seg_2 %over% gene_2]
        breakpoints_fwd <- endoapply(reads_fwd, extractBreakpointPerRead,
                                     gene_1, gene_2, "fwd")
        breakpoints_fwd <- breakpoints_fwd[lengths(breakpoints_fwd) != 0L]

        ## consider split reads in reverse orientation
        ## (overlapping both genes in antisense)
        reads_rev <- reads[seg_1 %over% gene_2_inv & seg_2 %over% gene_1_inv]
        breakpoints_rev <- endoapply(reads_rev, extractBreakpointPerRead,
                                     gene_1, gene_2, "rev")
        breakpoints_rev <- breakpoints_rev[lengths(breakpoints_rev) != 0L]

        out[[read]] <- c(breakpoints_fwd, breakpoints_rev)

    }

    return(out)

}

#' Extract breakpoint from split read
#'
#' @param r split read
#' @param gene_1 candidate gene 1
#' @param gene_2 candidate gene 2
#' @param orientation read orientation ("fwd" | "rev")
#'
#' @return GRanges with breakpoint

extractBreakpointPerRead <- function(r, gene_1, gene_2, orientation) {

    ## extract segments expected to map to gene_1 and gene_2
    if (orientation == "fwd") {

        seg_gene_1 <- r[which(mcols(r)$seg == 1)]
        seg_gene_2 <- r[which(mcols(r)$seg == 2)]

    } else if (orientation == "rev") {

        seg_gene_1 <- invertStrand(r[which(mcols(r)$seg == 2)])
        seg_gene_2 <- invertStrand(r[which(mcols(r)$seg == 1)])

    }

    if (gene_1 == gene_2) {

        flank_gene <- flank(gene_1, 1)
        flank_seg_gene_1 <- flank(seg_gene_1, 1)
        flank_seg_gene_2 <- flank(seg_gene_2, 1)
        w_1 <- width(range(flank_seg_gene_1, flank_gene))
        w_2 <- width(range(flank_seg_gene_2, flank_gene))
        if (w_1 > w_2) return(GRanges())

    }

    ## adjust coordinates in case of ambiguous alignments at the
    ## end/start of seg_gene_1/2 due to identical reference bases
    r_qwidth <- sapply(strsplit(mcols(r)$cigar, "M|S|H"),
                       function (x) sum(as.integer(x)))
    stopifnot(length(unique(r_qwidth)) == 1L)
    qwidth <- r_qwidth[1]
    awidth <- sum(width(r))
    d <- awidth - qwidth
    if (d > 0L) {

        seg_gene_2 <- resize(seg_gene_2, width(seg_gene_2) - d, fix = "end")

    }

    ## for gene 1 segment last pos marks breakpoint
    bp_1 <- flank(seg_gene_1, -1, FALSE)
    ## for gene 2 segment first read pos marks breakpoint
    bp_2 <- flank(seg_gene_2, -1, TRUE)
    ## invert strand by convention
    bp_2 <- invertStrand(bp_2)

    names(bp_1) <- NULL
    names(bp_2) <- NULL

    bp <- GRanges(c(bp_1, bp_2))
    names(bp) <- names(r)

    return(bp)

}
