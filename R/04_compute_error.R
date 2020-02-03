#' @title Computing Type 1 Error for Single Cell Expression Case-Control Differential
#'   Expression Analysis
#'
#' @name compute_error
#'
#' @description Computes type 1 error for single cell data that is cell-type specifc,
#'   hierarchical, and compositonal. This function computes type 1 error with the
#'   single-cell differential expression analysis tool 'MAST', using random
#'   effects to account for the correlation structure that exists among measures
#'   from cells within an individual. The type 1 error calculations will borrow
#'   information from the input data (or the package default data) to simulate
#'   data under a variety of pre-determined conditions. These conditions include
#'   foldchange, number of genes, number of samples (i.e., independent
#'   experimental units), and the mean number of cells per individual.
#'
#' @details Prior to running the \code{\link{error_hierarchicell}} function, it
#'   is important to run the \code{\link{filter_counts}} function followed by
#'   the \code{\link{compute_data_summaries}} function to build an R object that
#'   is in the right format for the following simulation function to properly
#'   work.
#'
#' @note Data should be \strong{only for cells of the specific cell-type} you
#'   are interested in simulating or computing type 1 error for. Data should also
#'   contain as many unique sample identifiers as possible. If you are inputing
#'   data that has less than 5 unique values for sample identifier (i.e.,
#'   independent experimental units), then the empirical estimation of the
#'   inter-individual heterogeneity is going to be very unstable. Finding such a
#'   dataset will be difficult at this time, but, over time (as experiments grow
#'   in sample size and the numbers of publically available single-cell RNAseq
#'   datasets increase), this should improve dramatically.
#'

NULL

#'@title Compute Type 1 Error for Single Cell Expression Case-Control Analysis
#'
#'@rdname error_hierarchicell
#'
#'@description Computes type 1 error for single cell data that is cell-type
#'  specifc, hierarchical, and compositonal. This function computes type 1 error
#'  with the single-cell differential expression analysis tool 'MAST', using
#'  random effects to account for the correlation structure that exists among
#'  measures from cells within an individual. The type 1 error calculations will
#'  borrow information from the input data (or the package default data) to
#'  simulate data under a variety of pre-determined conditions. These conditions
#'  include foldchange, number of genes, number of samples (i.e., independent
#'  experimental units), and the mean number of cells per individual.
#'
#'@details Prior to running the \code{\link{error_hierarchicell}} function, it
#'  is important to run the \code{\link{filter_counts}} function followed by the
#'  \code{\link{compute_data_summaries}} function to build an R object that is
#'  in the right format for the following simulation function to properly work.
#'
#'@note Data should be \strong{only for cells of the specific cell-type} you are
#'  interested in simulating or computing type 1 error for. Data should also
#'  contain as many unique sample identifiers as possible. If you are inputing
#'  data that has less than 5 unique values for sample identifier (i.e.,
#'  independent experimental units), then the empirical estimation of the
#'  inter-individual heterogeneity is going to be very unstable. Finding such a
#'  dataset will be difficult at this time, but, over time (as experiments grow
#'  in sample size and the numbers of publically available single-cell RNAseq
#'  datasets increase), this should improve dramatically.
#'
#'@param data_summaries an R object that has been output by the package's
#'  compute_data_summaries function. No default
#'
#'@param method a name. The method for differential expression to be used for
#'  the computation of error. Possible methods include: MAST with random effects
#'  ("MAST_RE"), MAST ("MAST"), MAST with batch effect correction
#'  ("MAST_Combat"), GLM assuming a tweedie distribution ("GLM_tweedie"), GLMM
#'  assuming a tweedie distribution ("GLMM_tweedie"), generalized estimating
#'  equations ("GEE1"), ROTS ("ROTS"), Monocle ("Monocle"), DESeq2 ("DESeq2").
#'  Defaults to "MAST_RE" which is the currently recommended analysis pipeline
#'  for single-cell data. See \code{\link{de_methods}} for more details on each
#'  of the methods.
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
#'  greater than 100, brings marginal improvements in type 1 error. Defaults to
#'  100.
#'
#'@param pval a number. The significance threshold (alpha) to use for
#'  significance. Defaults to 0.05.
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
#'@return The estimated error under the specified conditions when using 'MAST'
#'  with random effects to account for the correlation structure that exists
#'  among measures from cells within an individual.
#'
#'@examples
#'clean_expr_data <- filter_counts()
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'error_hierarchicell(data_summaries)
#'
#'@export

error_hierarchicell <- function(data_summaries,
                                   method = "MAST_RE",
                                   n_genes = 10000,
                                   n_per_group = 3,
                                   n_cases = n_per_group,
                                   n_controls = n_per_group,
                                   cells_per_individual = 100,
                                   pval = 0.05,
                                   foldchange = 1,
                                   decrease_dropout = 0){
  if (method == "MAST_RE") {


    if (!requireNamespace(c("MAST","SummarizedExperiment","lme4"),quietly = TRUE)){
      stop("The packages 'MAST', 'lme4', 'fitdistrplus', and \n
           'SummarizedExperiment' are required. Please install them.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"MAST\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      if (foldchange != 1){
        message("----------------------------------------------")
        message("Foldchange is not equal to 1, you are not simulating under the null.\nFor type 1 error rates, please keep foldchange equal to 1")
        message("----------------------------------------------")
      }

      if (n_genes < 1000){
        message("----------------------------------------------")
        message("Number of genes is less than 1,000.\nA low number of genes affects the compositional component of the simulation\nand limits the stability of your type 1 error calculations")
        message("----------------------------------------------")
      }

      if (cells_per_individual < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      all_genes <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_individual = cells_per_individual,
                                                           foldchange = foldchange,
                                                           decrease_dropout = decrease_dropout))

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
      fcHurdle <- merge(summaryDt[summaryDt$contrast=='StatusControl'
                                  & summaryDt$component=='logFC', c(1,7,5,6,8)],
                        summaryDt[summaryDt$contrast=='StatusControl'
                                  & summaryDt$component=='H', c(1,4)],
                        by = 'primerid')

      fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
      signif <- ifelse(fcHurdle[,6] < pval, 1, 0)
      rate <- mean(signif)
      message(paste0("Type 1 error rate is: ", rate))

    }
  } else if (method == "MAST") {
    if (!requireNamespace(c("MAST","SummarizedExperiment"),quietly = TRUE)){
      stop("The packages 'MAST', 'fitdistrplus', and \n
           'SummarizedExperiment' are required. Please install them.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"MAST\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      if (foldchange != 1){
        message("----------------------------------------------")
        message("Foldchange is not equal to 1, you are not simulating under the null.\nFor type 1 error rates, please keep foldchange equal to 1")
        message("----------------------------------------------")
      }

      if (n_genes < 1000){
        message("----------------------------------------------")
        message("Number of genes is less than 1,000.\nA low number of genes affects the compositional component of the simulation\nand limits the stability of your type 1 error calculations")
        message("----------------------------------------------")
      }

      if (cells_per_individual < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      all_genes <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_individual = cells_per_individual,
                                                           foldchange = foldchange,
                                                           decrease_dropout = decrease_dropout))

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

      zlmCond <- suppressMessages(MAST::zlm(~ ngeneson + Status,
                                            sca, method='glm',ebayes = F))

      summaryCond <- suppressMessages(MAST::summary(zlmCond,
                                                    doLRT='StatusControl'))
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[summaryDt$contrast=='StatusControl'
                                  & summaryDt$component=='logFC', c(1,7,5,6,8)],
                        summaryDt[summaryDt$contrast=='StatusControl'
                                  & summaryDt$component=='H', c(1,4)],
                        by = 'primerid')

      fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
      signif <- ifelse(fcHurdle[,6] < pval, 1, 0)
      rate <- mean(signif)
      message(paste0("Type 1 error rate is: ", rate))
    }


  } else if (method == "MAST_Combat") {
    if (!requireNamespace(c("MAST","SummarizedExperiment","sva"),quietly = TRUE)){
      stop("The packages 'MAST', 'fitdistrplus',  \n
           'SummarizedExperiment', and 'sva' are required. Please install them.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"MAST\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      if (foldchange != 1){
        message("----------------------------------------------")
        message("Foldchange is not equal to 1, you are not simulating under the null.\nFor type 1 error rates, please keep foldchange equal to 1")
        message("----------------------------------------------")
      }

      if (n_genes < 1000){
        message("----------------------------------------------")
        message("Number of genes is less than 1,000.\nA low number of genes affects the compositional component of the simulation\nand limits the stability of your type 1 error calculations")
        message("----------------------------------------------")
      }

      if (cells_per_individual < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      all_genes <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_individual = cells_per_individual,
                                                           foldchange = foldchange,
                                                           decrease_dropout = decrease_dropout))

      genecounts <- as.matrix(t(all_genes[,c(-1,-2,-3)]))
      coldata <- all_genes[,1:3]
      coldata$Status <- as.factor(coldata$Status)
      genecounts <- genecounts[which(apply(genecounts, 1, mean) > 5), ]
      genecounts <- genecounts[,rownames(coldata)]
      log2counts <- log2(genecounts + 1)
      log2counts <- sva::ComBat(log2counts,coldata$DonorID)

      fData <- data.frame(primerid=rownames(genecounts))
      sca <- suppressMessages(MAST::FromMatrix(exprsArray=log2counts, cData=coldata, fData=fData))

      cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
      SummarizedExperiment::colData(sca)$ngeneson <- scale(cdr2)
      SummarizedExperiment::colData(sca)$Status <-
        factor(SummarizedExperiment::colData(sca)$Status)
      SummarizedExperiment::colData(sca)$DonorID <-
        factor(SummarizedExperiment::colData(sca)$DonorID)

      zlmCond <- suppressMessages(MAST::zlm(~ ngeneson + Status,
                                            sca, method='glm',ebayes = F))

      summaryCond <- suppressMessages(MAST::summary(zlmCond,
                                                    doLRT='StatusControl'))
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[summaryDt$contrast=='StatusControl'
                                  & summaryDt$component=='logFC', c(1,7,5,6,8)],
                        summaryDt[summaryDt$contrast=='StatusControl'
                                  & summaryDt$component=='H', c(1,4)],
                        by = 'primerid')

      fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
      signif <- ifelse(fcHurdle[,6] < pval, 1, 0)
      rate <- mean(signif)
      message(paste0("Type 1 error rate is: ", rate))

      }



  } else if (method == "GLM_tweedie") {
    if (!requireNamespace(c("glmmTMB"),quietly = TRUE)){
      stop("The 'glmmTMB package is required. Please install it.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"glmmTMB\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      message("This function is slow and requires a lot of memory.")

      if (foldchange != 1){
        message("----------------------------------------------")
        message("Foldchange is not equal to 1, you are not simulating under the null.\nFor type 1 error rates, please keep foldchange equal to 1")
        message("----------------------------------------------")
      }

      if (n_genes < 1000){
        message("----------------------------------------------")
        message("Number of genes is less than 1,000.\nA low number of genes affects the compositional component of the simulation\nand limits the stability of your type 1 error calculations")
        message("----------------------------------------------")
      }

      if (cells_per_individual < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      all_genes <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_individual = cells_per_individual,
                                                           foldchange = foldchange,
                                                           decrease_dropout = decrease_dropout))

      genecounts <- as.matrix(t(all_genes[,c(-1,-2,-3)]))
      coldata <- all_genes[,1:3]
      coldata$Status <- as.factor(coldata$Status)
      coldata$DonorID <- as.factor(coldata$DonorID)
      genecounts <- genecounts[which(apply(genecounts, 1, mean) > 5), ]
      genecounts <- log(sweep(genecounts,2,apply(genecounts,2,mean),'/'))
      genecounts[which(genecounts == '-Inf')] <- 0
      genecounts <- t(genecounts[,rownames(coldata)])
      allcells <- cbind(coldata,genecounts)

      fitfixed <- lapply(4:ncol(allcells),
                         function(x){glmmTMB::glmmTMB(allcells[,x] ~ Status,
                                                      data=allcells,
                                                      family=glmmTMB::tweedie,
                                                      ziformula= ~0)})
      summaries <- lapply(fitfixed, summary)
      pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
      signif <- ifelse(pvalues < pval, 1, 0)
      rate <- mean(signif)
      message(paste0("Type 1 error rate is: ", rate))

    }


  } else if (method == "GLMM_tweedie") {
    if (!requireNamespace(c("glmmTMB"),quietly = TRUE)){
      stop("The 'glmmTMB package is required. Please install it.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"glmmTMB\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      message("This function is slow and requires a lot of memory.")

      if (foldchange != 1){
        message("----------------------------------------------")
        message("Foldchange is not equal to 1, you are not simulating under the null.\nFor type 1 error rates, please keep foldchange equal to 1")
        message("----------------------------------------------")
      }

      if (n_genes < 1000){
        message("----------------------------------------------")
        message("Number of genes is less than 1,000.\nA low number of genes affects the compositional component of the simulation\nand limits the stability of your type 1 error calculations")
        message("----------------------------------------------")
      }

      if (cells_per_individual < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      all_genes <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_individual = cells_per_individual,
                                                           foldchange = foldchange,
                                                           decrease_dropout = decrease_dropout))

      genecounts <- as.matrix(t(all_genes[,c(-1,-2,-3)]))
      coldata <- all_genes[,1:3]
      coldata$Status <- as.factor(coldata$Status)
      coldata$DonorID <- as.factor(coldata$DonorID)
      genecounts <- genecounts[which(apply(genecounts, 1, mean) > 5), ]
      genecounts <- log(sweep(genecounts,2,apply(genecounts,2,mean),'/'))
      genecounts[which(genecounts == '-Inf')] <- 0
      genecounts <- t(genecounts[,rownames(coldata)])
      allcells <- cbind(coldata,genecounts)

      fitmixed <- lapply(4:ncol(allcells),
                         function(x){glmmTMB::glmmTMB(allcells[,x] ~ Status + (1 | DonorID),
                                                      data=allcells,
                                                      family=glmmTMB::tweedie,
                                                      ziformula= ~0)})
      summaries <- lapply(fitmixed, summary)
      pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)$cond[2,4]})))
      signif <- ifelse(pvalues < pval, 1, 0)
      rate <- mean(signif)
      message(paste0("Type 1 error rate is: ", rate))

    }


  } else if (method == "GEE1") {
    if (!requireNamespace(c("geepack"),quietly = TRUE)){
      stop("The 'glmmTMB package is required. Please install it.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"geepack\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      message("This function is slow and requires a lot of memory.")

      if (foldchange != 1){
        message("----------------------------------------------")
        message("Foldchange is not equal to 1, you are not simulating under the null.\nFor type 1 error rates, please keep foldchange equal to 1")
        message("----------------------------------------------")
      }

      if (n_genes < 1000){
        message("----------------------------------------------")
        message("Number of genes is less than 1,000.\nA low number of genes affects the compositional component of the simulation\nand limits the stability of your type 1 error calculations")
        message("----------------------------------------------")
      }

      if (cells_per_individual < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      all_genes <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_individual = cells_per_individual,
                                                           foldchange = foldchange,
                                                           decrease_dropout = decrease_dropout))

      genecounts <- as.matrix(t(all_genes[,c(-1,-2,-3)]))
      coldata <- all_genes[,1:3]
      coldata$Status <- as.factor(coldata$Status)
      coldata$DonorID <- as.factor(coldata$DonorID)
      genecounts <- genecounts[which(apply(genecounts, 1, mean) > 5), ]
      genecounts <- log(sweep(genecounts,2,apply(genecounts,2,mean),'/'))
      genecounts[which(genecounts == '-Inf')] <- 0
      genecounts <- t(genecounts[,rownames(coldata)])
      allcells <- cbind(coldata,genecounts)

      fitgee <- lapply(4:ncol(allcells),
                         function(x){geepack::geeglm(allcells[,x] ~ Status,
                                                     data=allcells,
                                                     family=stats::gaussian(link="identity"),
                                                     id = DonorID,
                                                     corstr="exchangeable")})
      summaries <- lapply(fitgee, summary)
      pvalues <- as.numeric(unlist(lapply(summaries, function(x){stats::coef(x)[2,4]})))
      signif <- ifelse(pvalues < pval, 1, 0)
      rate <- mean(signif)
      message(paste0("Type 1 error rate is: ", rate))



    }
  } else if (method == "ROTS") {

    if (!requireNamespace(c("ROTS","tidyr"),quietly = TRUE)){
      stop("The packages 'ROTS', 'fitdistrplus', and \n
           'tidyr' are required. Please install them.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"ROTS\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      if (foldchange != 1){
        message("----------------------------------------------")
        message("Foldchange is not equal to 1, you are not simulating under the null.\nFor type 1 error rates, please keep foldchange equal to 1")
        message("----------------------------------------------")
      }

      if (n_genes < 1000){
        message("----------------------------------------------")
        message("Number of genes is less than 1,000.\nA low number of genes affects the compositional component of the simulation\nand limits the stability of your type 1 error calculations")
        message("----------------------------------------------")
      }

      if (cells_per_individual < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      all_genes <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_individual = cells_per_individual,
                                                           foldchange = foldchange,
                                                           decrease_dropout = decrease_dropout))

      genecounts <- as.matrix(t(all_genes[,c(-1,-2,-3)]))
      coldata <- all_genes[,1:3]
      coldata$Status <- as.factor(coldata$Status)
      genecounts <- genecounts[which(apply(genecounts, 1, mean) > 5), ]
      genecounts <- genecounts[,rownames(coldata)]
      coldata$Status <- ifelse(coldata$Status == "Control",0,1)
      coldata$DonorID <- as.factor(coldata$DonorID)
      genecounts <- genecounts[,rownames(coldata)]
      results <- suppressMessages(ROTS::ROTS(data = genecounts, groups = coldata$Status, B = 1000, K = 500, seed = 4119))
      results <- suppressMessages(ROTS::summary.ROTS(results, num.genes=nrow(genecounts)))
      results <- as.data.frame(results)
      results <- stats::na.omit(results)

      nsig <- length(results$pvalue[results$pvalue<=pval])
      n.non.singular <- length(results$pvalue)
      rate <- nsig/n.non.singular
      message(paste0("Type 1 error rate is: ", rate))

      }


  } else if (method == "Monocle") {
    if (!requireNamespace(c("monocle","tidyr"),quietly = TRUE)){
      stop("The packages 'monocle', 'fitdistrplus', and \n
           'tidyr' are required. Please install them.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"monocle\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      if (foldchange != 1){
        message("----------------------------------------------")
        message("Foldchange is not equal to 1, you are not simulating under the null.\nFor type 1 error rates, please keep foldchange equal to 1")
        message("----------------------------------------------")
      }

      if (n_genes < 1000){
        message("----------------------------------------------")
        message("Number of genes is less than 1,000.\nA low number of genes affects the compositional component of the simulation\nand limits the stability of your type 1 error calculations")
        message("----------------------------------------------")
      }

      if (cells_per_individual < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      all_genes <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_individual = cells_per_individual,
                                                           foldchange = foldchange,
                                                           decrease_dropout = decrease_dropout))

      genecounts <- as.matrix(t(all_genes[,c(-1,-2,-3)]))
      coldata <- all_genes[,1:3]
      coldata$Status <- as.factor(coldata$Status)
      genecounts <- genecounts[which(apply(genecounts, 1, mean) > 5), ]
      genecounts <- genecounts[,rownames(coldata)]
      coldata$Status <- ifelse(coldata$Status == "Control",0,1)
      coldata$DonorID <- as.factor(coldata$DonorID)
      genecounts <- genecounts[,rownames(coldata)]

      features <- data.frame(rownames(genecounts),"Function")
      colnames(features) <- c("gene_short_name","Function")
      rownames(features) <- features$gene_short_name
      features <- methods::new("AnnotatedDataFrame", data = features)
      pheno <- methods::new("AnnotatedDataFrame", data = coldata)
      celldat <- monocle::newCellDataSet(genecounts,phenoData = pheno, featureData = features, expressionFamily = VGAM::negbinomial.size())
      celldat <- monocle::detectGenes(celldat, min_expr = 0.1)
      celldat <- BiocGenerics::estimateSizeFactors(celldat)
      celldat <- BiocGenerics::estimateDispersions(celldat)
      results <- monocle::differentialGeneTest(celldat, fullModelFormulaStr = "~Status")
      results <- stats::na.omit(results)

      nsig <- length(results$pval[results$pval<=pval])
      n.non.singular <- length(results$pval)
      rate <- nsig/n.non.singular
      message(paste0("Type 1 error rate is: ", rate))

      }

    }

}
