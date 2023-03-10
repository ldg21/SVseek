% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractBreakpoints.R
\name{extractBreakpointPerRead}
\alias{extractBreakpointPerRead}
\title{Extract breakpoint from split read}
\usage{
extractBreakpointPerRead(r, gene_1, gene_2, orientation)
}
\arguments{
\item{r}{split read}

\item{gene_1}{candidate gene 1}

\item{gene_2}{candidate gene 2}

\item{orientation}{read orientation ("fwd" | "rev")}
}
\value{
GRanges with breakpoint
}
\description{
Extract breakpoint from split read
}
