#' Parse a BAM file
#'
#' loadBAM takes an input BAM file and an optional bounding range to get
#' aligned read segments and metadata from the file.
#'
#' The function removes reads with unmapped mates, reads marked as
#' duplicates, and alignments with mapping quality zero. In the output,
#' metadata columns R1 and R2 indicate whether alignments are for
#' sequence reads 1 or 2, respectively. The output only includes alignments
#' for seqlevels in GenomeInfoDb::standardChromosomes().
#'
#' @param input_bam Input bam file
#' @param bounds GRanges with bounding ranges (optional)
#'
#' @return GRanges with reads
#'
#' @export

loadBAM <- function(input_bam, bounds = GRanges()) {

    bounds <- reduce(bounds)

    # fields from bam
    what <- c(
        "qname",
        "flag",
        "rname",
        "pos",
        "qwidth",
        "mapq",
        "cigar",
        "mrnm",
        "mpos",
        "isize")
    param <- ScanBamParam(
        flag = scanBamFlag(isUnmappedQuery = FALSE),
        what = what,
        which = bounds)
    ga <- readGAlignments(input_bam, use.names = TRUE, param = param)
    gr <- granges(ga, use.mcols = TRUE)

    ## remove reads with unmapped mates
    gr <- gr[!bitwAnd(gr$flag, 0x8)]

    ## remove reads marked as duplicates
    gr <- gr[!bitwAnd(gr$flag, 0x400)]

    ## remove reads with mapq 0
    gr <- gr[which(gr$mapq > 0)]

    ## flag alignment as first or second read
    mcols(gr)$R1 <- as.logical(bitwAnd(gr$flag, 0x40))
    mcols(gr)$R2 <- as.logical(bitwAnd(gr$flag, 0x80))

    ## remove unnecessary seqlevels
    seqlevels(gr) <- standardChromosomes(gr)

    return(gr)

}
