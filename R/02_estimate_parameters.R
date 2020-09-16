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

utils::globalVariables(c("Dispersion",
                         "DonorID",
                         "DropOut",
                         "DropOutStD",
                         "GrandMean",
                         "InterStD",
                         "IntraMean",
                         "ToSep",
                         "n_individuals"))

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
#'@param type an identifier for the type of data being submitted. If it is raw
#'  counts put "Raw", if it is TPM or some other normalized counts per million
#'  then type "PerMillion" or "Norm". The program assumes data is in one of
#'  these two formats. Other normalizations (i.e., logs) that have negative
#'  values will cause the program to malfunction.
#'
#'@return A data.frame of the summary data as well as two vectors for the
#'  cell-wise dropout rates and library sizes. The data.frame includes the
#'  gene-wise grand means, inter-individual standard deviations,
#'  intra-individual standard deviations, and dropout rates.
#'
#'@examples
#'\donttest{clean_expr_data <- filter_counts()
#'data_summaries <- compute_data_summaries(clean_expr_data, type = "Norm")}
#'
#'@export

compute_data_summaries <- function(expr,
                                   type = "Norm"){

  if (type == "Raw"){

    n_individuals <- length(unique(expr[,2]))

    message("Normalizing ... ")

    ids <- expr[,c(1,2)]
    expr <- expr[,c(-1,-2)]
    expr <- t(as.matrix(apply(expr,1,function(x){(x/sum(x))*1000000}))) #RPM
    expr <- cbind(ids,expr)

    message("Removing highly correlated genes")

    reduced <- expr[,c(-1,-2)]
    genelist <- colnames(reduced)
    uncorrelatedgenes <- ids

    for (i in 1:500){
      genename <- sample(genelist,1)
      drawngene <- reduced[,genename]
      uncorrelatedgenes <- cbind(uncorrelatedgenes,drawngene)
      correlations <- abs(stats::cor(drawngene,reduced))
      genelist <- names(correlations[,which(correlations < 0.25)])
      reduced <- reduced[,genelist]
      if (ncol(reduced) < 10) break
    }

    expr <- uncorrelatedgenes
    rm(ids,uncorrelatedgenes,genelist,correlations,drawngene,genename,reduced)

    message("Computing sample means, dropout rates, and dispersion ... ")

    computevar <- function(a){tapply(a,expr[,2],function(a){stats::var(a[a != 0])})}
    temp.intravar <- sapply(expr[,c(-1,-2)],computevar)
    rownames(temp.intravar) <- paste0(rownames(temp.intravar),"_Var") #Compute sample variances

    computemeans <- function(x){tapply(x,expr[,2],function(a){mean(a[a != 0])})}
    temp.intrameans <- sapply(expr[,c(-1,-2)],computemeans)
    rownames(temp.intrameans) <- paste0(rownames(temp.intrameans),"_Mean") #Compute sample means

    computedrop <- function(x){tapply(x,expr[,2],function(a){length(a[which(a == 0)])/length(a)})}
    temp.drop <- sapply(expr[,c(-1,-2)],computedrop)
    rownames(temp.drop) <- paste0(rownames(temp.drop),"_Drop") #Compute sample dropout

    intravar <- na.omit(as.data.frame(cbind(c(temp.intravar),c(temp.intrameans)))) #Compute sample dispersion
    colnames(intravar) <- c("IntraVar","IntraMean")
    intravar$Dispersion <- (intravar$IntraMean**2)/((intravar$IntraVar) - intravar$IntraMean)

    message("Computing final data summaries ... ")

    temp.intra <- as.data.frame(t(do.call("rbind",list(temp.intrameans,temp.intravar,temp.drop))))
    temp.intra$InterStD <- as.numeric(as.character(apply(temp.intra[,1:n_individuals],1,function(a){stats::sd(a[!is.na(a)])}))) #Compute inter SD
    temp.intra$GrandMean <- as.numeric(as.character(apply(temp.intra[,1:n_individuals],1,function(a){mean(a[!is.na(a)])}))) #Compute grand mean
    temp.intra$DropOut <- as.numeric(as.character(apply(temp.intra[,(1+(n_individuals*2)):(n_individuals*3)],1,function(a){mean(a[!is.na(a)])}))) #Compute gene dropout
    temp.intra$DropOutStD <- as.numeric(as.character(apply(temp.intra[,(1+(n_individuals*2)):(n_individuals*3)],1,function(a){stats::sd(a[!is.na(a)])})))
    temp.intra <- temp.intra[ ,(ncol(temp.intra)-3):ncol(temp.intra)]
    main_summary <- as.data.frame(temp.intra)

    list(n_individuals,main_summary,intravar)

  } else {

    n_individuals <- length(unique(expr[,2]))

    message("Computing sample means, dropout rates, and dispersion ... ")

    computevar <- function(a){tapply(a,expr[,2],function(a){stats::var(a[a != 0])})}
    temp.intravar <- sapply(expr[,c(-1,-2)],computevar)
    rownames(temp.intravar) <- paste0(rownames(temp.intravar),"_Var") #Compute sample variances

    computemeans <- function(x){tapply(x,expr[,2],function(a){mean(a[a != 0])})}
    temp.intrameans <- sapply(expr[,c(-1,-2)],computemeans)
    rownames(temp.intrameans) <- paste0(rownames(temp.intrameans),"_Mean") #Compute sample means

    computedrop <- function(x){tapply(x,expr[,2],function(a){length(a[which(a == 0)])/length(a)})}
    temp.drop <- sapply(expr[,c(-1,-2)],computedrop)
    rownames(temp.drop) <- paste0(rownames(temp.drop),"_Drop") #Compute sample dropout

    intravar <- na.omit(as.data.frame(cbind(c(temp.intravar),c(temp.intrameans)))) #Compute sample dispersion
    colnames(intravar) <- c("IntraVar","IntraMean")
    intravar$Dispersion <- (intravar$IntraMean**2)/((intravar$IntraVar) - intravar$IntraMean)

    message("Computing final data summaries ... ")

    temp.intra <- as.data.frame(t(do.call("rbind",list(temp.intrameans,temp.intravar,temp.drop))))
    temp.intra$InterStD <- as.numeric(as.character(apply(temp.intra[,1:n_individuals],1,function(a){stats::sd(a[!is.na(a)])}))) #Compute inter SD
    temp.intra$GrandMean <- as.numeric(as.character(apply(temp.intra[,1:n_individuals],1,function(a){mean(a[!is.na(a)])}))) #Compute grand mean
    temp.intra$DropOut <- as.numeric(as.character(apply(temp.intra[,(1+(n_individuals*2)):(n_individuals*3)],1,function(a){mean(a[!is.na(a)])}))) #Compute gene dropout
    temp.intra$DropOutStD <- as.numeric(as.character(apply(temp.intra[,(1+(n_individuals*2)):(n_individuals*3)],1,function(a){stats::sd(a[!is.na(a)])})))
    temp.intra <- temp.intra[ ,(ncol(temp.intra)-3):ncol(temp.intra)]
    main_summary <- as.data.frame(temp.intra)

    list(n_individuals,main_summary,intravar)

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
#'@param plot a TRUE/FALSE statement for the output of a plot to observe a
#'  histogram of grand means
#'
#'@param data_summaries an R object that has been output by the package's
#'  compute_data_summaries function.
#'
#'@return A vector of length two, where the first number is the shape of the
#'  gamma distribution and the second number is the rate of the gamma
#'  distribution
#'
#'@examples
#'\donttest{clean_expr_data <- filter_counts()
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'gene_mean_params <- approximate_gene_mean(data_summaries)}
#'
#'@export

approximate_gene_mean <- function(data_summaries, plot = FALSE){
  if (!requireNamespace("fitdistrplus",quietly = TRUE)){
    stop("The 'fitdistrplus' package is required for this function to work. Please install it",
         call. = FALSE)
  } else if ((plot == TRUE)) {

    if (!requireNamespace("ggplot2",quietly = TRUE)){
      stop("The 'ggplot2' package is required for plotting. Please install it",
           call. = FALSE)
    } else {
      data_summaries <- data_summaries[[2]]
      data_summaries <- data_summaries[which(data_summaries$GrandMean > 0), ]
      data_summaries <- data_summaries[which(data_summaries$GrandMean < 10000), ]

      message("Plotting distribution of grand means")

      print(ggplot2::ggplot(data_summaries, ggplot2::aes(GrandMean))
            + ggplot2::geom_histogram(ggplot2::aes(y=..density..),fill="cornflowerblue",color = "black")
            + ggplot2::ylab("Density"))
      gene_mean <- fitdistrplus::fitdist(data_summaries$GrandMean, "gamma",method = "mle")
      gene_mean_shape <- gene_mean$estimate[1]
      gene_mean_rate <- gene_mean$estimate[2]
      as.numeric(c(gene_mean_shape,gene_mean_rate))
    }

  } else {
    data_summaries <- data_summaries[[2]]
    data_summaries <- data_summaries[which(data_summaries$GrandMean > 0), ]
    data_summaries <- data_summaries[which(data_summaries$GrandMean < 10000), ]
    gene_mean <- fitdistrplus::fitdist(data_summaries$GrandMean, "gamma",method = "mle")
    gene_mean_shape <- gene_mean$estimate[1]
    gene_mean_rate <- gene_mean$estimate[2]
    as.numeric(c(gene_mean_shape,gene_mean_rate))
  }
}

#'@title Model Dispersion as a function of Intra Means
#'
#'@rdname model_dispersion
#'
#'@description This function estimates the parameters for the gamma distribution
#'  of dispersion parameters.
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
#'@param plot a TRUE/FALSE statement for the output of a plot to observe how
#'  well gene dropout behaves as a function of the intra means
#'
#'@return A vector of length two, where the first number is the shape of the
#'  gamma distribution and the second number is the rate of the gamma
#'  distribution
#'
#'@examples
#'\donttest{clean_expr_data <- filter_counts()
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'dispersion_params <- model_dispersion(data_summaries)}
#'
#'@export

model_dispersion <- function(data_summaries, plot=FALSE){
  if ((plot == TRUE)) {
    if (!requireNamespace("ggplot2",quietly = TRUE)){
      stop("The 'ggplot2' package is required for plotting. Please install it",
           call. = FALSE)
    } else {
      data_summaries <- data_summaries[[3]]
      data_summaries <- data_summaries[which(data_summaries$Dispersion > 0),]

      message("Plotting dispersion")

      print(suppressWarnings(ggplot2::ggplot(data_summaries,ggplot2::aes(x=IntraMean,y=sqrt(Dispersion)),
                                       environment = environment()) +
                         ggplot2::geom_point() +
                         ggplot2::stat_smooth(method="glm",
                                              formula = y ~ I(1/x), size = 1,
                                              method.args = list(family = gaussian(link = "log")))))

      suppressWarnings(temp <- stats::glm(data = data_summaries, Dispersion ~ I(1/IntraMean), family = gaussian(link = "log")))
      suppressWarnings(temp <- summary(temp))
      disp.beta0 <- temp$coefficients[1,1]
      disp.beta1 <- temp$coefficients[2,1]
      as.numeric(c(disp.beta0,disp.beta1))
    }
  } else {
    data_summaries <- data_summaries[[3]]
    data_summaries <- data_summaries[which(data_summaries$Dispersion > 0),]
    temp <- stats::glm(data = data_summaries, Dispersion ~ I(1/IntraMean), family = gaussian(link = "log"))
    suppressWarnings(temp <- summary(temp))
    disp.beta0 <- temp$coefficients[1,1]
    disp.beta1 <- temp$coefficients[2,1]
    as.numeric(c(disp.beta0,disp.beta1))

  }
}


#'@title Approximate Gene Dropout as a gamma
#'
#'@rdname approximate_gene_drop
#'
#'@description This function estimates the parameters for the gamma distribution
#'  of gene dropout.
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
#'@param plot a TRUE/FALSE statement for the output of a plot to observe a
#'  histogram of dropout
#'
#'@return A vector of length two, where the first number is the shape of the
#'  gamma distribution and the second number is the rate of the gamma
#'  distribution
#'
#'@examples
#'\donttest{clean_expr_data <- filter_counts()
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'gene_drop_params <- approximate_gene_drop(data_summaries)}
#'
#'@export

approximate_gene_drop <- function(data_summaries, plot = FALSE){
    if (!requireNamespace("fitdistrplus",quietly = TRUE)){
      stop("The 'fitdistrplus' package is required for this function to work. Please install it",
           call. = FALSE)
    } else if ((plot == TRUE)) {

      if (!requireNamespace("ggplot2",quietly = TRUE)){
        stop("The 'ggplot2' package is required for plotting. Please install it",
             call. = FALSE)
      } else {
      data_summaries <- data_summaries[[2]]
      data_summaries$DropOut <- 1 - data_summaries$DropOut
      data_summaries <- data_summaries[which(data_summaries$DropOut > 0), ]
      message("Plotting distribution of dropout")

      print(ggplot2::ggplot(data_summaries, ggplot2::aes(DropOut))
            + ggplot2::geom_histogram(ggplot2::aes(y=..density..),fill="cornflowerblue",color = "black")
            + ggplot2::ylab("Density"))
      drop_gamma <- fitdistrplus::fitdist(data_summaries$DropOut, "gamma",method = "mle")
      drop_shape <- drop_gamma$estimate[1]
      drop_rate <- drop_gamma$estimate[2]
      as.numeric(c(drop_shape,drop_rate))
      }

    } else {
      data_summaries <- data_summaries[[2]]
      data_summaries$DropOut <- 1 - data_summaries$DropOut
      data_summaries <- data_summaries[which(data_summaries$DropOut > 0), ]
      drop_gamma <- fitdistrplus::fitdist(data_summaries$DropOut, "gamma",method = "mle")
      drop_shape <- drop_gamma$estimate[1]
      drop_rate <- drop_gamma$estimate[2]
      as.numeric(c(drop_shape,drop_rate))
    }
}

#'@title Model Dropout SD as a Quadratic Function of the Mean Dropout
#'
#'@rdname model_drop_sd
#'
#'@description This function models gene dropout SD as a quadratic function of the
#'  dropout mean. Visualizing this fit and looking for oddities and extreme
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
#'  well gene dropout SD behaves as a quadratic function of mean dropout
#'
#'@return A plot (if plot=TRUE) and a vector of length three, where each
#'  number is the estimate of the intercept or slope for the model.
#'
#'@examples
#'\donttest{clean_expr_data <- filter_counts()
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'gene_drop_betas <- model_drop_sd(data_summaries)}
#'
#'@export

model_drop_sd <- function(data_summaries, plot=FALSE){
  if ((plot == TRUE)) {
    if (!requireNamespace("ggplot2",quietly = TRUE)){
      stop("The 'ggplot2' package is required for plotting. Please install it",
           call. = FALSE)
    } else {
      n_individuals <- data_summaries[[1]]
      data_summaries <- data_summaries[[2]]

      message("Plotting dropout")

      print(suppressWarnings(ggplot2::ggplot(data_summaries,ggplot2::aes(x=DropOut,y=DropOutStD),
                                       environment = environment()) +
                         ggplot2::geom_point() +
                         ggplot2::stat_smooth(method = "glm", formula = y ~ x + I(x**2), fullrange = TRUE) +
                         ggplot2::ylim(c(0,1))))

      suppressWarnings(temp <- stats::glm(data = data_summaries, DropOutStD ~ DropOut + I(DropOut**2)))
      suppressWarnings(temp <- summary(temp))
      drop.beta0 <- temp$coefficients[1,1]
      drop.beta1 <- temp$coefficients[2,1]
      drop.beta2 <- temp$coefficients[3,1]
      as.numeric(c(drop.beta0,drop.beta1,drop.beta2))
    }
  } else {
    n_individuals <- data_summaries[[1]]
    data_summaries <- data_summaries[[2]]

    suppressWarnings(temp <- stats::glm(data = data_summaries, DropOutStD ~ DropOut + I(DropOut**2)))
    suppressWarnings(temp <- summary(temp))
    drop.beta0 <- temp$coefficients[1,1]
    drop.beta1 <- temp$coefficients[2,1]
    drop.beta2 <- temp$coefficients[3,1]
    as.numeric(c(drop.beta0,drop.beta1,drop.beta2))

  }
}


#'@title Model Inter-Individual Heterogeneity as a Linear Function of the Grand
#'  Mean
#'
#'@rdname model_inter
#'
#'@description This function models inter-individual heterogeneity as a log
#'  function of the grand mean. Visualizing this fit and looking for oddities
#'  and extreme outliers may be helpful.
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
#'  well inter-individual standard deviation behaves as a linear function of the
#'  grand mean
#'
#'@return A plot (if plot=TRUE) a vector of length two, where the first
#'  number is the estimate of the intercept for the log model second number is
#'  the estimate of the slope for the log model.
#'
#'@examples
#'\donttest{clean_expr_data <- filter_counts()
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'inter_betas <- model_inter(data_summaries)}
#'
#'@export

model_inter <- function(data_summaries, plot=FALSE){
  if ((plot == TRUE)) {
    if (!requireNamespace("ggplot2",quietly = TRUE)){
      stop("The 'ggplot2' package is required for plotting. Please install it",
           call. = FALSE)
    } else {
      n_individuals <- data_summaries[[1]]
      data_summaries <- na.omit(data_summaries[[2]])

      message("Plotting inter-individual standard deviation")

      print(suppressWarnings(ggplot2::ggplot(data_summaries,ggplot2::aes(x=GrandMean,y=InterStD),
                      environment = environment()) +
      ggplot2::geom_point() +
      ggplot2::stat_smooth(method="glm",
                             formula = y ~ 0 + x, size = 1)))

      suppressWarnings(temp <- stats::glm(data = data_summaries, InterStD ~ 0 + GrandMean))
      suppressWarnings(temp <- summary(temp))
      inter.beta1 <- temp$coefficients[1,1]
      as.numeric(inter.beta1)
    }
  } else {
    n_individuals <- data_summaries[[1]]
    data_summaries <- na.omit(data_summaries[[2]])
    suppressWarnings(temp <- stats::glm(data = data_summaries, InterStD ~ 0 + GrandMean))
    suppressWarnings(temp <- summary(temp))
    inter.beta1 <- temp$coefficients[1,1]
    as.numeric(inter.beta1)
  }
}
