% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_estimate_parameters.R
\name{approximate_gene_drop}
\alias{approximate_gene_drop}
\title{Approximate Gene Dropout as a gamma}
\usage{
approximate_gene_drop(data_summaries, plot = FALSE)
}
\arguments{
\item{data_summaries}{an R object that has been output by the package's
compute_data_summaries function.}

\item{plot}{a TRUE/FALSE statement for the output of a plot to observe a
histogram of dropout}
}
\value{
A vector of length two, where the first number is the shape of the
 gamma distribution and the second number is the rate of the gamma
 distribution
}
\description{
This function estimates the parameters for the gamma distribution
 of gene dropout.
}
\details{
Prior to estimating the simulation parameters, it is important to run
 the compute_data_summaries function to build an object that is in the right
 format for the following estimation functions to properly compute. The
 'ggplot2' package is required for the plotting component of this function to
 work. Please install it.
}
\examples{
\donttest{clean_expr_data <- filter_counts()
data_summaries <- compute_data_summaries(clean_expr_data)
gene_drop_params <- approximate_gene_drop(data_summaries)}

}
