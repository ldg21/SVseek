#' Check read flags
#'
#' @param gr GRanges with reads
#'
#' @return boolean value

checkFlags <- function(gr) {

    ## check all reads paired
    all(bitwAnd(gr$flag, 0x1)) &
    ## check no unmapped reads
    all(!bitwAnd(gr$flag, 0x4)) &
    ## check no unmapped mates
    all(!bitwAnd(gr$flag, 0x8)) &
    ## check no multi-mappers
    all(!bitwAnd(gr$flag, 0x200)) &
    ## check no duplicates
    all(!bitwAnd(gr$flag, 0x400)) &
    ## check strand directions are logical
    all(gr$R1 | gr$R2)

}

#' Validate reads
#'
#' @param gr GRanges with reads
#'
#' @return void

validate_gr <- function(gr) {

    # Check class
    if (!"GRanges" %in% class(gr))
        warning("gr must be GRanges object.")

    # Check flags
    if (!checkFlags(gr))
        warning("Invalid flag.")

    # Check seqlevels
    if (!all(seqlevels(gr) == seqlevels(keepStandardChromosomes(gr))))
        warning("Seqlevels contain non standard chromosomes.")

    if(!all(names(gr) == mcols(gr)$qname))
        warning("Names and qname do not match.")

    # Check metadata
    # actual_metadata <- sapply(mcols(gr), class)
    # expected_metadata <- c(
    #     qname = ,
    #     flag = "integer",
    #     rname = "factor",
    #     pos = "integer",
    #     qwidth = "integer",
    #     mapq = "integer",
    #     cigar = "character",
    #     mrnm = "factor",
    #     mpos = "integer",
    #     isize = "integer",
    #     R1 = "logical",
    #     R2 = "logical",
    #     seg = "integer"
    # )

    # if (!all(actual_metadata == expected_metadata[names(actual_metadata)]))
    #     warning("Metadata columns are not in the right form.")

}

#' Validate breakpoint GRanges
#'
#' @param breakpoint GRanges object
#'
#' @return void

validate_breakpoint <- function(breakpoint) {

    ## Check class
    if (!inherits(breakpoint, "GRanges"))
        warning("Breakpoint must be GRanges object.")

    ## Check length
    if (!length(breakpoint) == 2)
        warning("Breakpoint must be of length 2.")

    ## Check strand
    if (!all(strand(breakpoint) %in% c("+", "-")))
        warning("Breakpoint has no strand information.")

    ## Check seqlevels
    if (!all(seqlevels(breakpoint) ==
             seqlevels(keepStandardChromosomes(breakpoint))))
        warning("Seqlevels contain non standard chromosomes.")

}

#' Validate breakpoints GRanges or GRangesList
#'
#' @param breakpoints GRanges or GRangesList object
#'
#' @return void

validate_breakpoints <- function(breakpoints) {

    ## Check class, use validate_breakpoint if GRanges
    if (inherits(breakpoints, "GRanges"))
        return(validate_breakpoint(breakpoints))

    if (!inherits(breakpoints, "GRangesList"))
        warning("Breakpoints must be GRanges or GRangesList.")

    ## Check seqlevels
    if (!all(seqlevels(breakpoints) ==
             seqlevels(keepStandardChromosomes(breakpoints))))
        warning("Seqlevels contain non standard chromosomes.")

    ## Validate each internal breakpoint
    sapply(breakpoints, validate_breakpoint)

}

#' Validate candidate genes
#'
#' @param genes GRanges with candidate genes
#'
#' @return void

validate_genes <- function(genes) {

    if (!inherits(genes, "GRanges"))
        warning("Genes must be GRanges object.")
    if (length(genes) != 1 & length(genes) != 2)
        warning("Genes must have length of 1 or 2.")

}

validate_vaf <- function(vaf) {

    if (!inherits(vaf, "GRanges"))
        warning("VAF must be GRanges object.")

    ## Check metadata
    actual_metadata <- sapply(mcols(vaf), class)
    expected_metadata <- c(
        nonsupp = "integer",
        supp = "integer",
        supp_split = "integer",
        total = "integer",
        vaf = "numeric"
    )

    if (!all(actual_metadata == expected_metadata[names(actual_metadata)]))
        warning("Metadata columns are not in the right form.")

}
