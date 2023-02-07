#' Get reads that overlap regions of interest
#'
#' @param R1 read 1
#' @param R2 read 2
#' @param region_1 region 1
#' @param region_2 region 2
#' @param min_overlap_1 overlap for region 1
#' @param min_overlap_2 overlap for region 2
#' @param supporting TRUE for supporting reads, FALSE for nonsupporting reads
#' @param max_overlap maximum considered overlap
#'
#' @return list with overlapping reads

overlapsRegions <- function(R1, R2, region_1, region_2, min_overlap_1,
                            min_overlap_2, supporting, max_overlap) {

    R1_ext <- pc(R1, flank(R1, max_overlap, start = FALSE))
    R2_ext <- pc(R2, flank(R2, max_overlap, start = FALSE))
    overlapping <- intersect(
        names(which(sum(width(pintersect(R1_ext, region_1))) >= min_overlap_1)),
        names(which(sum(width(pintersect(R2_ext, region_2))) >= min_overlap_2)))

    if (supporting) {

        ## supporting non-split reads must be contained in overlap regions
        R1_ns <- unlist(R1[lengths(R1) == 1L])
        names(R1_ns) <- mcols(R1_ns)$qname
        R2_ns <- unlist(R2[lengths(R2) == 1L])
        names(R2_ns) <- mcols(R2_ns)$qname
        noncontained <- union(
            names(R1_ns)[width(pintersect(R1_ns, region_1)) < width(R1_ns)],
            names(R2_ns)[width(pintersect(R2_ns, region_2)) < width(R2_ns)])
        overlapping <- setdiff(overlapping, noncontained)

    }

    return(overlapping)

}
