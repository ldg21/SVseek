#' QIAseq ES Example Data
#'
#' @format List of intermediate objects
#'
#' \describe{
#' \item{\code{gr}}{GRanges from parseBAM}
#' \item{\code{breakpoints}}{Top breakpoint from detectBreakpoints}
#' \item{\code{genes}}{GRanges with candidate genes}
#' \item{\code{vaf}}{Breakpoint with VAF from calculateVAF}
#' \item{\code{supporting_reads}}{String list of supporting reads}
#' \item{\code{nonsupporting_reads}}{String list of nonsupporting reads}
#' }
#'
#' @examples data_fusion$gr

"data_fusion"

#' QIAseq HB Example Data
#'
#' @format List of intermediate objects
#'
#' \describe{
#' \item{\code{gr}}{GRanges from parseBAM}
#' \item{\code{breakpoints}}{Top breakpoint from detectBreakpoints}
#' \item{\code{genes}}{GRanges with candidate genes}
#' \item{\code{vaf}}{Breakpoint with VAF from calculateVAF}
#' \item{\code{supporting_reads}}{String list of supporting reads}
#' \item{\code{nonsupporting_reads}}{String list of nonsupporting reads}
#' }
#'
#' @examples data_deletion$gr

"data_deletion"
