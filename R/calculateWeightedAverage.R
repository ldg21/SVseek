#' Calculate weighted average variant allele fraction (VAF)
#'
#' @param vaf GRanges or GRangesList with variants and VAFs
#'
#' @return weighted average VAF
#'
#' @examples calculateWeightedAverage(data_fusion$vaf)
#'
#' @export

calculateWeightedAverage <- function(vaf) {

    if (is(vaf, "GRangesList")) {

        out <- sapply(vaf, function(v) {
            calculateWeightedAverage(v)
        })
        return(out)

    }

    validate_vaf(vaf)

    ## if VAF is only calculated for one variant, return that value
    if (any(vaf$total == 0)) {

        return(vaf$vaf[vaf$total != 0])

    }

    ## calculate weighted average VAF
    w <- c(vaf[1]$total, vaf[2]$total) / (vaf[1]$total + vaf[2]$total)
    w_ave <- w[1] * vaf[1]$vaf + w[2] * vaf[2]$vaf
    names(w_ave) <- NULL
    return(w_ave)

}
