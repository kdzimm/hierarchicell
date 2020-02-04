#' @title Computing Power for Single Cell Expression Case-Control Differential
#'   Expression Analysis
#'
#' @name compute_power
#'
#' @description Computes power for single cell data that is cell-type specifc,
#'   hierarchical, and compositonal. This function computes power using random
#'   effects to account for the correlation structure that exists among measures
#'   from cells within an individual. The power calculations will borrow
#'   information from the input data (or the package default data) to simulate
#'   data under a variety of pre-determined conditions. These conditions include
#'   foldchange, number of genes, number of samples (i.e., independent
#'   experimental units), and the mean number of cells per individual.
#'
#' @details Prior to running the \code{\link{power_hierarchicell}} function, it
#'   is important to run the \code{\link{filter_counts}} function followed by
#'   the \code{\link{compute_data_summaries}} function to build an R object that
#'   is in the right format for the following simulation function to properly
#'   work.
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

NULL

#'@title Compute Power for Single Cell Expression Case-Control Analysis
#'
#'@rdname power_hierarchicell
#'
#'@description Computes power for single cell data that is cell-type specifc,
#'  hierarchical, and compositonal. This function computes power using random
#'  effects to account for the correlation structure that exists among measures
#'  from cells within an individual. The power calculations will borrow
#'  information from the input data (or the package default data) to simulate
#'  data under a variety of pre-determined conditions. These conditions include
#'  foldchange, number of genes, number of samples (i.e., independent
#'  experimental units), and the mean number of cells per individual.
#'
#'@details Prior to running the \code{\link{power_hierarchicell}} function, it
#'  is important to run the \code{\link{filter_counts}} function followed by the
#'  \code{\link{compute_data_summaries}} function to build an R object that is
#'  in the right format for the following simulation function to properly work.
#'
#'@note Data should be \strong{only for cells of the specific cell-type} you are
#'  interested in simulating or computing power for. Data should also contain as
#'  many unique sample identifiers as possible. If you are inputing data that
#'  has less than 5 unique values for sample identifier (i.e., independent
#'  experimental units), then the empirical estimation of the inter-individual
#'  heterogeneity is going to be very unstable. Finding such a dataset will be
#'  difficult at this time, but, over time (as experiments grow in sample size
#'  and the numbers of publically available single-cell RNAseq datasets
#'  increase), this should improve dramatically.
#'
#'@param data_summaries an R object that has been output by the package's
#'  compute_data_summaries function. No default
#'
#'@param n_genes an integer. The number of genes you would like to simulate for
#'  your dataset. Too large of a number may cause memory failure and may slow
#'  the simulation down tremendously. We recommend an integer less than 40,000.
#'  Defaults to 10,000.
#'
#'@param n_per_group an integer. The number of independent samples per
#'  case/control group for simulation. Creates a balanced design, for unbalanced
#'  designs, specify n_cases and n_controls separately. If not specifying a
#'  foldchange, the number of cases and controls does not matter. Defaults to 3.
#'
#'@param n_cases an integer. The number of independent control samples for
#'  simulation. Defaults to n_per_group.
#'
#'@param n_controls an integer. The number of independent case samples for
#'  simulation. Defaults to n_per_group.
#'
#'@param cells_per_individual an integer. The mean number of cells per
#'  individual you would like to simulate. Too large of a number may cause
#'  memory failure and may slow the simulation down tremendously. We recommend
#'  an integer less than 300, but more is possible. We note that anything
#'  greater than 100, brings marginal improvements in power. Defaults to 100.
#'
#'@param pval a number. The significance threshold (alpha) to use for
#'  significance. Defaults to 0.05.
#'
#'@param foldchange an integer between 1 and 10. The amount of fold change to
#'  simulate a difference in expression between case and control groups. The
#'  foldchange changes genes in either direction, so a foldchange of 2  would
#'  cause the mean expression in cases to be either twice the amount or half the
#'  amount for any particular gene. Defaults to 1.
#'
#'@param decrease_dropout a numeric proportion between 0 and 1. The proportion
#'  by which you would like to simulate decreasing the amount of dropout in your
#'  data. For example, if you would like to simulate a decrease in the amount of
#'  dropout in your data by twenty percent, then 0.2 would be appropriate. This
#'  component of the simulation allows the user to adjust the proportion of
#'  dropout if they believe future experiments or runs will have improved
#'  calling rates (due to improved methods or improved cell viability) and
#'  thereby lower dropout rates. Defaults to 0.
#'
#'@param alter_dropout_cases a numeric proportion between 0 and 1. The
#'  proportion by which you would like to simulate decreasing the amount of
#'  dropout between case control groups. For example, if you would like to
#'  simulate a decrease in the amount of dropout in your cases by twenty
#'  percent, then 0.2 would be appropriate. This component of the simulation
#'  allows the user to adjust the proportion of dropout if they believe the
#'  stochastic expression of a gene will differ between cases and controls. For
#'  a two-part hurdle model, like MAST implements, this will increase your
#'  ability to detect differences. Defaults to 0.
#'
#'@return The estimated power under the specified conditions when using random
#'  effects to account for the correlation structure that exists among measures
#'  from cells within an individual.
#'
#'@examples
#'clean_expr_data <- filter_counts()
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'power_hierarchicell(data_summaries)
#'
#'@export

power_hierarchicell <- function(data_summaries,
                                   n_genes = 10000,
                                   n_per_group = 3,
                                   n_cases = n_per_group,
                                   n_controls = n_per_group,
                                   cells_per_individual = 100,
                                   pval = 0.05,
                                   foldchange = 1,
                                   decrease_dropout = 0,
                                   alter_dropout_cases = 0){




  if (!requireNamespace(c("MAST","SummarizedExperiment","lme4"),quietly = TRUE)){
    stop("The packages 'MAST', 'lme4', 'fitdistrplus', and \n
           'SummarizedExperiment' are required. Please install them.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"MAST\")' to see if this is the issue.",
         call. = FALSE)
  } else {

    if (n_genes < 1000){
      message("----------------------------------------------")
      message("Number of genes is less than 1,000.\nA low number of genes affects the compositional component of the simulation.")
      message("----------------------------------------------")
    }

    if (cells_per_individual < 50){
      message("----------------------------------------------")
      message("Mean number of cells per individual is less than 50.\nThe probability of complete separation will start to increase.")
      message("----------------------------------------------")
    }

    all_genes <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                         n_genes = n_genes,
                                                         n_per_group = n_per_group,
                                                         n_cases = n_cases,
                                                         n_controls = n_controls,
                                                         cells_per_individual = cells_per_individual,
                                                         foldchange = foldchange,
                                                         decrease_dropout = decrease_dropout,
                                                         alter_dropout_cases = alter_dropout_cases))

    genecounts <- as.matrix(t(all_genes[,c(-1,-2,-3)]))
    coldata <- all_genes[,1:3]
    coldata$Status <- as.factor(coldata$Status)
    genecounts <- genecounts[which(apply(genecounts, 1, mean) > 5), ]
    genecounts <- genecounts[,rownames(coldata)]
    log2counts <- log2(genecounts + 1)

    fData <- data.frame(primerid=rownames(genecounts))
    sca <- suppressMessages(MAST::FromMatrix(exprsArray=log2counts, cData=coldata, fData=fData))

    cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
    SummarizedExperiment::colData(sca)$ngeneson <- scale(cdr2)
    SummarizedExperiment::colData(sca)$Status <-
      factor(SummarizedExperiment::colData(sca)$Status)
    SummarizedExperiment::colData(sca)$DonorID <-
      factor(SummarizedExperiment::colData(sca)$DonorID)

    zlmCond <- suppressMessages(MAST::zlm(~ ngeneson + Status + (1 | DonorID),
                                          sca, method='glmer',ebayes = F,
                                          strictConvergence = FALSE))

    summaryCond <- suppressMessages(MAST::summary(zlmCond,
                                                  doLRT='StatusControl'))
    summaryDt <- summaryCond$datatable
    print(head(summaryDt))
    fcHurdle <-  summaryDt[summaryDt$contrast=='StatusControl' & summaryDt$component=='H', c(1,4)]

    fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
    signif <- ifelse(fcHurdle[,2] < pval, 1, 0)
    rate <- mean(signif)
    message(paste0("Power is: ", rate))

  }


}
