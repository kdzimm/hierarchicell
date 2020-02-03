#' @title hierarchicell: A package for simulating compositional, cell-type
#'   specific, and hierarchical single-cell expression data
#'
#' @description The hierarchicell package estimates important parameters from
#'   single-cell RNAseq expression counts before simulating expression data that
#'   is compositional, cell-type specific, and hierarchical in nature. With the
#'   simulated data, hierarchicell is then able to compute power calculations
#'   under a variety of conditions. The package consists of four important
#'   categories of functions: data loading and cleaning, empirical estimation of
#'   distributions, simulating expression data, and computing power.
#'
#' @section Data loading and filtering function:
#'
#'   The data loading and cleaning function is very basic, but data input is
#'   critical to the package working correctly. If no input data is given, the
#'   package default data will be used for simulation and power calculations.
#'   For more detailed information see: \code{\link{filter_counts}}
#'
#' @section Empirical estimation of distributions:
#'
#'   The most fundamental component of this package is in the estimation of the
#'   simulation parameters. The functions to estimate parameters for the
#'   simulation estimate the empirical distributions for library size, dropout
#'   rate, and global gene means and model the hierarchical variance structure
#'   of the input data. For more detailed information see:
#'   \code{\link{empirical_estimation}}
#'
#' @section Simulating expression data:
#'
#'   With the parameters estimated, the package can simulate data under a
#'   variety of pre-determined conditions. These conditions include foldchange,
#'   number of genes, number of samples (i.e., independent experimental units),
#'   and the mean number of cells per individual. For more detailed information
#'   see: \code{\link{simulate_count_matrix}}
#'
#' @section Computing type I error:
#'
#'   With the parameters estimated, the package can compute type 1 error rates
#'   for a number of different tools under a variety of pre-determined
#'   conditions. These conditions include number of genes, number of samples
#'   (i.e., independent experimental units), and the mean number of cells per
#'   individual. For more detailed information see: \code{\link{compute_error}}
#'
#'
#' @section Computing power:
#'
#'   With the parameters estimated, the package can compute power under a
#'   variety of pre-determined conditions. These conditions include foldchange,
#'   number of genes, number of samples (i.e., independent experimental units),
#'   and the mean number of cells per individual. For more detailed information
#'   see: \code{\link{compute_power}}
#'
#'
#'
#' @docType package
#' @name hierarchicell

NULL

