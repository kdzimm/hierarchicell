#' @title Filter gene expression matrix
#'
#' @rdname filter_counts
#'
#' @description This function loads and filters the input data for the
#'   subsequent steps. The data loading and cleaning function is very basic, but
#'   data input is critical to the package working correctly. If no input data
#'   is given, the package defaults to using normal mucosal cells data for the
#'   simulation and power calculations (see \code{\link{mucosal_cells}}).
#'
#' @details Input data should be formatted as follows:
#'
#'   \tabular{cccccc}{ \strong{Cell_ID} \tab \strong{Individual_ID} \tab
#'   \strong{Gene1} \tab \strong{Gene2} \tab \strong{Gene3} \tab \strong{...}
#'   \cr Cell1_Ind1 \tab Ind1 \tab 12 \tab 24 \tab 0 \tab ...\cr Cell2_Ind1 \tab
#'   Ind1 \tab 11 \tab 2 \tab 0 \tab ... \cr Cell3_Ind1 \tab Ind1 \tab 10 \tab 0
#'   \tab 0 \tab ... \cr Cell4_Ind1 \tab Ind1 \tab 0 \tab 124 \tab 10 \tab ...
#'   \cr Cell1_Ind2 \tab Ind2 \tab 9 \tab 37 \tab 18 \tab ... \cr Cell2_Ind2
#'   \tab Ind2 \tab 0 \tab 29 \tab 0 \tab ... }
#'
#'   Where the unique cell identifier is in column one and the sample identifier
#'   is in column two with the remaining columns all being genes.
#'
#' @note Data should be \strong{only for cells of the specific cell-type} you
#'   are interested in simulating or computing power for. Data should also
#'   contain as many unique sample identifiers as possible. If you are inputing
#'   data that has less than 5 unique values for sample identifier (i.e.,
#'   independent experimental units), then the empirical estimation of the
#'   inter-individual heterogeneity is going to be very unstable. Finding such a
#'   dataset will be difficult at this time, but, over time (as experiments grow
#'   in sample size and the numbers of publically available single-cell RNAseq
#'   datasets increase), this should improve dramatically.
#'
#' @param expr a data.frame where the unique cell identifier is in column one
#'   and the sample identifier is in column two with the remaining columns all
#'   being genes.
#' @param gene_thresh the mean expression threshold for retaining genes.
#'   Defaults to 0.
#' @param cell_thresh the mean expression threshold for retaining cells.
#'   Defaults to 5.
#' @return a data.frame that has filtered out cells with mean count < 5 and
#'   genes with mean count = 0
#' @examples
#'
#' n_genes <- 10
#' n_cells <- 10
#'
#' make_data <- function(x){
#'  mu_random <- round(rgamma(n=1, shape=1, rate=0.001),0)
#'  size_random <- runif(n=1, min=0, max=3)
#'  rnbinom(n_cells, size=size_random, mu=mu_random)
#' }
#'
#' expr_dat <- as.data.frame(replicate(n_genes,make_data()))
#' expr_dat$CellID <- paste0("Cell",1:n_cells)
#' expr_dat$IND <- "IND1"
#' expr_dat <- expr_dat[,c(11,12,1:10)]
#' clean_expr_data <- filter_counts(expr_dat)
#'
#' @export


filter_counts <- function(expr=hierarchicell:::mucosal_cells, gene_thresh=0, cell_thresh=5){
  if (identical(expr, hierarchicell:::mucosal_cells)){
    message("Filtering default dataset")
  } else {
    message("Filtering user input")
  }

  hierar.env <- new.env(parent = emptyenv())

  hierar.env$ids <- expr[ ,c(1,2)]
  hierar.env$genes <- expr[ ,c(-1,-2)]
  hierar.env$genes <- hierar.env$genes[ , which(apply(hierar.env$genes, 2, mean) > gene_thresh)]
  hierar.env$ids <- hierar.env$ids[which(apply(hierar.env$genes, 1, mean) > cell_thresh), ]
  hierar.env$genes <- hierar.env$genes[which(apply(hierar.env$genes, 1, mean) > cell_thresh), ]
  message("Genes and cells have been filtered, ready for estimating parameters")
  as.data.frame(cbind(hierar.env$ids, hierar.env$genes))
}



