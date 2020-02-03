#' @title Empirical Estimation of Parameters
#'
#' @name empirical_estimation
#'
#' @description The most fundamental component of this package is in the
#'   estimation of the simulation parameters. The functions to estimate
#'   parameters for the simulation estimate the empirical distributions for
#'   library size, dropout rate, and global gene means and model the
#'   hierarchical variance structure of the input data.
#'
#' @details Prior to estimating the simulation parameters, it is important to
#'   run the \code{\link{filter_counts}} function to build a data.frame that is in the right
#'   format for the following estimation functions to properly compute.
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

#'@title Compute Data Summaries
#'
#'@rdname compute_data_summaries
#'
#'@description This function computes cell-wise dropout rates and library sizes
#'  on the unnormalized data before computing the gene-wise grand means,
#'  gene-wise dropout rates, inter-individual variance, and intra-individual
#'  variance on the normalized data.
#'
#'@details Prior to estimating the data summaries, it is important to run the
#'  \code{\link{filter_counts}} function to build a data.frame that is in the
#'  right format for the following estimation functions to properly compute.
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
#'@param expr a data.frame that has been output by filter_counts where the
#'  unique cell identifier is in column one and the sample identifier is in
#'  column two with the remaining columns all being genes.
#'
#'@return A data.frame of the summary data as well as two vectors for the
#'  cell-wise dropout rates and library sizes. The data.frame includes the
#'  gene-wise grand means, inter-individual standard deviations,
#'  intra-individual standard deviations, and dropout rates.
#'
#'@examples
#'clean_expr_data <- filter_counts()
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'
#'@export

compute_data_summaries <- function(expr){

  n_individuals <- length(unique(expr[,2]))

  message("Computing cell-wise summaries ... ")

  compute_drop <- function(a){length(a[which(a == 0)])/length(a)}
  cell_dropout <- apply(expr[,c(-1,-2)],1,compute_drop) #Compute cell dropout
  cell_libraries <- apply(expr[,c(-1,-2)],1,sum) #Compute library size

  message("Normalizing data ...")

  ids <- expr[,c(1,2)]
  expr <- expr[,c(-1,-2)]
  expr <- t(as.matrix(apply(expr,1,function(x){-log((x/sum(x)))}))) #Normalize
  expr[expr == 'Inf'] <- 0
  expr <- as.data.frame(cbind(ids,expr))

  message("Computing sample means and standard deviations ... ")

  computesd <- function(a){tapply(a,expr[,2],stats::sd)}
  temp.intrasd <- sapply(expr[,c(-1,-2)],computesd)
  rownames(temp.intrasd) <- paste0(rownames(temp.intrasd),"_SD") #Compute sample SDs

  computemeans <- function(x){tapply(x,expr[,2],mean)}
  temp.intrameans <- sapply(expr[,c(-1,-2)],computemeans)
  rownames(temp.intrameans) <- paste0(rownames(temp.intrameans),"_Mean") #Compute sample means

  temp.intra <- as.data.frame(t(do.call("rbind",list(temp.intrameans,temp.intrasd))))

  message("Computing final data summaries ... ")

  temp.intra$InterStD <- apply(temp.intra[,1:n_individuals],1,stats::sd) #Compute inter SD
  temp.intra$GrandMean <- apply(temp.intra[,1:n_individuals],1,mean) #Compute grand mean
  temp.intra$GrandMeanSq <- temp.intra$GrandMean**2
  temp.intra$IntraStD <- apply(temp.intra[,(1+n_individuals):(n_individuals*2)],1,stats::median) #Compute intra SD
  temp.intra$DropOut <- apply(expr[,c(-1,-2)],2,compute_drop) #Compute gene dropout
  main_summary <- as.data.frame(temp.intra)

  list(cell_libraries,cell_dropout,n_individuals,main_summary)
}


#'@title Approximate Library Size Distribution
#'
#'@rdname approximate_library_sizes
#'
#'@description This function estimates the parameters for the lognormal
#'  distribution of the library sizes.
#'
#'@details Prior to estimating the simulation parameters, it is important to run
#'  the compute_data_summaries function to build an object that is in the right
#'  format for the following estimation functions to properly compute.The
#'  'fitdistrplus' package is required for this function to work. Please install
#'  it.
#'
#'@param data_summaries an R object that has been output by the package's
#'  compute_data_summaries function.
#'
#'@return A vector of length two, where the first number is the mean of the
#'  lognormal distribution and the second number is the standard deviation of
#'  the lognormal distribution
#'
#'@examples
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'library_params <- approximate_library_sizes(data_summaries)
#'
#'@export


approximate_library_sizes <- function(data_summaries){
  if (!requireNamespace("fitdistrplus",quietly = TRUE)) {
    stop("The 'fitdistrplus' package is required for this function to work. Please install it",
         call. = FALSE)
  } else {
  cell_library <- fitdistrplus::fitdist(data_summaries[[1]], "lnorm",method = "mle")
  cell_library_mean <- cell_library$estimate[1]
  cell_library_sd <- cell_library$estimate[2]
  as.numeric(c(cell_library_mean,cell_library_sd))
  }
}

#'@title Approximate Cell Dropout Distribution
#'
#'@rdname approximate_cell_dropout
#'
#'@description This function estimates the parameters for the normal
#'  distribution of cell dropout.
#'
#'@details Prior to estimating the simulation parameters, it is important to run
#'  the compute_data_summaries function to build an object that is in the right
#'  format for the following estimation functions to properly compute. The
#'  'fitdistrplus' package is required for this function to work. Please install
#'  it.
#'
#'@param data_summaries an R object that has been output by the package's
#'  compute_data_summaries function.
#'
#'@return A vector of length two, where the first number is the mean of the
#'  normal distribution and the second number is the standard deviation of the
#'  normal distribution
#'
#'@examples
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'cell_drop_params <- approximate_cell_dropout(data_summaries)
#'
#'@export

approximate_cell_dropout <- function(data_summaries){
  if (!requireNamespace("fitdistrplus",quietly = TRUE)){
    stop("The 'fitdistrplus' package is required for this function to work. Please install it",
         call. = FALSE)
  } else {
  cell_dropout <- fitdistrplus::fitdist(data_summaries[[2]], "norm",method = "mle")
  cell_dropout_mean <- cell_dropout$estimate[1]
  cell_dropout_sd <- cell_dropout$estimate[2]
  as.numeric(c(cell_dropout_mean,cell_dropout_sd))
  }
}

#'@title Approximate Gene Mean Distribution
#'
#'@rdname approximate_gene_mean
#'
#'@description This function estimates the parameters for the gamma distribution
#'  of gene means.
#'
#'@details Prior to estimating the simulation parameters, it is important to run
#'  the compute_data_summaries function to build an object that is in the right
#'  format for the following estimation functions to properly compute.The
#'  'fitdistrplus' package is required for this function to work. Please install
#'  it
#'
#'@param data_summaries an R object that has been output by the package's
#'  compute_data_summaries function.
#'
#'@return A vector of length two, where the first number is the shape of the
#'  gamma distribution and the second number is the rate of the gamma
#'  distribution
#'
#'@examples
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'gene_mean_params <- approximate_gene_mean(data_summaries)
#'
#'@export

approximate_gene_mean <- function(data_summaries){
  if (!requireNamespace("fitdistrplus",quietly = TRUE)){
    stop("The 'fitdistrplus' package is required for this function to work. Please install it",
         call. = FALSE)
  } else {
  data_summaries <- data_summaries[[4]]
  data_summaries <- data_summaries[which(data_summaries$GrandMean > 0), ]
  gene_mean <- fitdistrplus::fitdist(data_summaries$GrandMean, "gamma",method = "mle")
  gene_mean_shape <- gene_mean$estimate[1]
  gene_mean_rate <- gene_mean$estimate[2]
  as.numeric(c(gene_mean_shape,gene_mean_rate))
  }
}

#'@title Model Gene Dropout as a Linear Function of the Grand Mean
#'
#'@rdname model_gene_drop
#'
#'@description This function models gene dropout as a linear function of the
#'  grand mean. Visualizing this fit and looking for oddities and extreme
#'  outliers may be helpful.
#'
#'@details Prior to estimating the simulation parameters, it is important to run
#'  the compute_data_summaries function to build an object that is in the right
#'  format for the following estimation functions to properly compute. The
#'  'ggplot2' package is required for the plotting component of this function to
#'  work. Please install it.
#'
#'@param data_summaries an R object that has been output by the package's
#'  compute_data_summaries function.
#'
#'@param plot a TRUE/FALSE statement for the output of a plot to observe how
#'  well gene dropout behaves as a linear function of the grand mean
#'
#'@return A plot (if plot=TRUE) and a vector of length two, where the first
#'  number is the estimate of the intercept (beta0) and the second number is the
#'  estimate of the slope (beta1). The estimate of the intercept should be very
#'  close to one.
#'
#'@examples
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'gene_drop_betas <- model_gene_drop(data_summaries)
#'
#'@export

model_gene_drop <- function(data_summaries, plot=FALSE){
  if ((plot == TRUE)) {
    if (!requireNamespace("ggplot2",quietly = TRUE)){
    stop("The 'ggplot2' package is required for plotting. Please install it",
         call. = FALSE)
    } else {
    n_individuals <- data_summaries[[3]]
    data_summaries <- data_summaries[[4]]
    data_summaries <- data_summaries[apply(data_summaries[,1:n_individuals],1,function(x){all(x != 0)}), ]

    message("Plotting gene-wise dropout")

    ggplot2::ggplot(data_summaries,ggplot2::aes(x=GrandMean,y=DropOut),
                      environment = environment()) +
      ggplot2::geom_point() +
      ggplot2::stat_smooth(method="lm",formula = y ~ x , size = 1) +
      ggplot2::ylim(c(0,1))
      ggplot2::ggsave("Gene_Dropout.pdf")

    temp <- stats::lm(data = data_summaries, DropOut ~ GrandMean)
    temp <- summary(temp)
    dropout.beta0 <- temp$coefficients[1,1]
    dropout.beta1 <- temp$coefficients[2,1]
    as.numeric(c(dropout.beta0,dropout.beta1))
    }
  } else {
    n_individuals <- data_summaries[[3]]
    data_summaries <- data_summaries[[4]]
    data_summaries <- data_summaries[apply(data_summaries[,1:n_individuals],1,function(x){all(x != 0)}), ]
    temp <- stats::lm(data = data_summaries, DropOut ~ GrandMean)
    temp <- summary(temp)
    dropout.beta0 <- temp$coefficients[1,1]
    dropout.beta1 <- temp$coefficients[2,1]
    as.numeric(c(dropout.beta0,dropout.beta1))
  }
}

#'@title Model Inter-Individual Heterogeneity as a Quadratic Function of the
#'  Grand Mean
#'
#'@rdname model_inter
#'
#'@description This function models inter-individual heterogeneity as a
#'  quadratic function of the grand mean. Visualizing this fit and looking for
#'  oddities and extreme outliers may be helpful.
#'
#'@details Prior to estimating the simulation parameters, it is important to run
#'  the compute_data_summaries function to build an object that is in the right
#'  format for the following estimation functions to properly compute. The
#'  'ggplot2' package is required for the plotting component of this function to
#'  work. Please install it.
#'
#'@param data_summaries an R object that has been output by the package's
#'  compute_data_summaries function.
#'
#'@param plot a TRUE/FALSE statement for the output of a plot to observe how
#'  well inter-individual standard deviation behaves as a quadratic function of
#'  the grand mean
#'
#'@return A plot (if plot=TRUE) and a vector of length two, where the first
#'  number is the estimate of beta1 and the second number is the estimate of
#'  beta2 (function is forced through the origin, so intercept is zero).
#'
#'@examples
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'inter_betas <- model_inter(data_summaries)
#'
#'@export

model_inter <- function(data_summaries, plot=FALSE){
  if ((plot == TRUE)) {
    if (!requireNamespace("ggplot2",quietly = TRUE)){
      stop("The 'ggplot2' package is required for plotting. Please install it",
           call. = FALSE)
    } else {
      n_individuals <- data_summaries[[3]]
      data_summaries <- data_summaries[[4]]
      data_summaries <- data_summaries[apply(data_summaries[,1:n_individuals],1,function(x){all(x != 0)}), ]

      message("Plotting inter-individual standard deviation")

      ggplot2::ggplot(data_summaries,ggplot2::aes(x=GrandMean,y=InterStD),
                      environment = environment()) +
        ggplot2::geom_point() +
        ggplot2::stat_smooth(method="lm",
                             formula = y ~ 0 + x + I(x**2), size = 1)
      ggplot2::ggsave("Inter_Sample_Heterogeneity.pdf")

      temp <- stats::lm(data = data_summaries, InterStD ~ 0 + GrandMean + GrandMeanSq)
      temp <- summary(temp)
      inter.beta1 <- temp$coefficients[1,1]
      inter.beta2 <- temp$coefficients[2,1]
      as.numeric(c(inter.beta1,inter.beta2))
    }
  } else {
    n_individuals <- data_summaries[[3]]
    data_summaries <- data_summaries[[4]]
    data_summaries <- data_summaries[apply(data_summaries[,1:n_individuals],1,function(x){all(x != 0)}), ]
    temp <- stats::lm(data = data_summaries, InterStD ~ 0 + GrandMean + GrandMeanSq)
    temp <- summary(temp)
    inter.beta1 <- temp$coefficients[1,1]
    inter.beta2 <- temp$coefficients[2,1]
    as.numeric(c(inter.beta1,inter.beta2))
  }
}

#'@title Model Intra-Individual Variance as a Quadratic Function of the Grand
#'  Mean
#'
#'@rdname model_intra
#'
#'@description This function models intra-individual standard deviation as a
#'  quadratic function of the grand mean. Visualizing this fit and looking for
#'  oddities and extreme outliers may be helpful.
#'
#'@details Prior to estimating the simulation parameters, it is important to run
#'  the compute_data_summaries function to build an object that is in the right
#'  format for the following estimation functions to properly compute. The
#'  'ggplot2' package is required for the plotting component of this function to
#'  work. Please install it.
#'
#'@param data_summaries an R object that has been output by the package's
#'  compute_data_summaries function.
#'
#'@param plot a TRUE/FALSE statement for the output of a plot to observe how
#'  well intra-individual standard deviation behaves as a quadratic function of
#'  the grand mean
#'
#'@return A plot (if plot=TRUE) and a vector of length three, where the first
#'  number is the estimate of beta0 (the intercept) and the second and third
#'  numbers are the estimates of beta1 and beta2, respectively.
#'
#'@examples
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'intra_betas <- model_intra(data_summaries)
#'
#'@export

model_intra <- function(data_summaries, plot=FALSE){
  if ((plot == TRUE)) {
    if (!requireNamespace("ggplot2",quietly = TRUE)){
      stop("The 'ggplot2' package is required for plotting. Please install it",
           call. = FALSE)
    } else {
      n_individuals <- data_summaries[[3]]
      data_summaries <- data_summaries[[4]]
      data_summaries <- data_summaries[apply(data_summaries[,1:n_individuals],1,function(x){all(x != 0)}), ]

      message("Plotting intra-individual standard deviation")

      ggplot2::ggplot(data_summaries,ggplot2::aes(x=GrandMean,y=IntraStD),
                      environment = environment()) +
        ggplot2::geom_point() +
        ggplot2::stat_smooth(method="lm",
                             formula = y ~ x + I(x**2), size = 1)
      ggplot2::ggsave("Intra_Sample_Variance.pdf")

      temp <- stats::lm(data = data_summaries, IntraStD ~ GrandMean + GrandMeanSq)
      temp <- summary(temp)
      intra.beta0 <- temp$coefficients[1,1]
      intra.beta1 <- temp$coefficients[2,1]
      intra.beta2 <- temp$coefficients[3,1]
      as.numeric(c(intra.beta0,intra.beta1,intra.beta2))
    }
  } else {
    n_individuals <- data_summaries[[3]]
    data_summaries <- data_summaries[[4]]
    data_summaries <- data_summaries[apply(data_summaries[,1:n_individuals],1,function(x){all(x != 0)}), ]
    temp <- stats::lm(data = data_summaries, IntraStD ~ GrandMean + GrandMeanSq)
    temp <- summary(temp)
    intra.beta0 <- temp$coefficients[1,1]
    intra.beta1 <- temp$coefficients[2,1]
    intra.beta2 <- temp$coefficients[3,1]
    as.numeric(c(intra.beta0,intra.beta1,intra.beta2))
  }
}
