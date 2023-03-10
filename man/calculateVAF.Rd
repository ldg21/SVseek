% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateVAF.R
\name{calculateVAF}
\alias{calculateVAF}
\title{Calculate variant allele fraction}
\usage{
calculateVAF(
  gr,
  breakpoints,
  technology,
  min_overlap = 19,
  min_overlap_flanking = min_overlap,
  min_overlap_variant = min_overlap,
  supp = NULL,
  nonsupp = NULL
)
}
\arguments{
\item{gr}{GRanges with reads}

\item{breakpoints}{GRanges of 2 breakends representing a single variant or
GRangesList of multiple variants}

\item{technology}{Technology of reads ("qiaseq" | "other")}

\item{min_overlap}{minimum overlap for each side of breakend (default = 19)}

\item{min_overlap_flanking}{minimum overlap for flanking region
(set to min_overlap if not specified)}

\item{min_overlap_variant}{minimum overlap for variant region
(set to min_overlap if not specified)}

\item{supp}{supporting reads (optional)}

\item{nonsupp}{nonsupporting reads (optional)}
}
\value{
breakpoints annotated with number of supporting reads,
nonsupporting reads, and VAF
}
\description{
Calculate variant allele fraction by counting supporting and nonsupporting
reads for each structural variant.
}
\details{
The breakpoints argument can be a GRanges with a single variant, defined by
two breakends, or GRangesList with multiple variants.
}
\examples{
calculateVAF(data_fusion$gr, data_fusion$breakpoint, "qiaseq")

}
