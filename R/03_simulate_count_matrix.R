#' @title Simulating Expression Data
#'
#' @name simulate_count_matrix
#'
#' @description The simulation of single cell data that is cell-type specifc,
#'   hierarchical, and compositonal is the primary purpose of this package. This
#'   simulation will borrow information from the input data (or the package
#'   default data) to simulate data under a variety of pre-determined
#'   conditions. These conditions include foldchange, number of genes, number of
#'   samples (i.e., independent experimental units), and the mean number of
#'   cells per individual.
#'
#' @details Prior to running the \code{\link{simulate_hierarchicell}} function,
#'   it is important to run the \code{\link{filter_counts}} function followed by
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

#'@title Simulate Expression Data
#'
#'@rdname simulate_hierarchicell
#'
#'@description This function will compute a  simulation that will borrow
#'  information from the input data (or the package default data) to simulate
#'  data under a variety of pre-determined conditions. These conditions include
#'  foldchange, number of genes, number of samples (i.e., independent
#'  experimental units), and the mean number of cells per individual. The
#'  simulation incorporates information about the cell-wise dropout rates and
#'  library sizes from the unnormalized data and the gene-wise grand means,
#'  gene-wise dropout rates, inter-individual variance, and intra-individual
#'  variance from the normalized data.
#'
#'@details Prior to running the \code{\link{simulate_hierarchicell}} function,
#'  it is important to run the \code{\link{filter_counts}} function followed by
#'  the \code{\link{compute_data_summaries}} function to build an R object  that
#'  is in the right format for the following simulation function to properly
#'  work.
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
#'  the simulation down tremendously. We recommend an integer less than 100,000.
#'  Defaults to 1,000.
#'
#'@param n_per_group an integer. The number of independent samples per
#'  case/control group for simulation. Use when "binary" is specified as the
#'  outcome. Creates a balanced design, for unbalanced designs, specify n_cases
#'  and n_controls separately. If not specifying a foldchange, the number of
#'  cases and controls does not matter. Defaults to 3.
#'
#'@param n_cases an integer. The number of independent control samples for
#'  simulation. Defaults to n_per_group.
#'
#'@param n_controls an integer. The number of independent case samples for
#'  simulation. Defaults to n_per_group.
#'
#'@param cells_per_control an integer. The mean number of cells per control you
#'  would like to simulate. Too large of a number may cause memory failure and
#'  may slow the simulation down tremendously. We recommend an integer less than
#'  1,000.Defaults to 150.
#'
#'@param cells_per_case an integer. The mean number of cells per case you would
#'  like to simulate. Too large of a number may cause memory failure and may
#'  slow the simulation down tremendously. We recommend an integer less than
#'  1,000.Defaults to 150.
#'
#'@param ncells_variation_type either "Poisson", "NB", or "Fixed". Allows the
#'  number of cells per individual to be fixed at exactly the specified number
#'  of cells per individual, vary slightly with a poisson distribution with a
#'  lambda equal to the specified number of cells per individual, or a negative
#'  binomial with a mean equal to the specified number of cells and dispersion
#'  size equal to one.Defaults to "Poisson".
#'
#'@param foldchange an integer between 1 and 10. The amount of fold change to
#'  simulate a difference in expression between case and control groups. The
#'  foldchange changes genes in either direction, so a foldchange of 2  would
#'  cause the mean expression in cases to be either twice the amount or half the
#'  amount for any particular gene. Defaults to 1.
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
#'@param decrease_dropout a numeric proportion between 0 and 1. The proportion
#'  by which you would like to simulate decreasing the amount of dropout in your
#'  data. For example, if you would like to simulate a decrease in the amount of
#'  dropout in your data by twenty percent, then 0.2 would be appropriate. This
#'  component of the simulation allows the user to adjust the proportion of
#'  dropout if they believe future experiments or runs will have improved
#'  calling rates (due to improved methods or improved cell viability) and
#'  thereby lower dropout rates. Defaults to 0.
#'
#'@param tSNE_plot a TRUE/FALSE statement for the output of a tSNE plot to
#'  observe the global behavior of your simulated data. Seurat will need to be
#'  installed for this function to properly work. Defaults to FALSE.
#'
#'
#'@return A data.frame of the simulated data or potentially a pdf of a tSNE plot
#'  (if tSNE_plot=TRUE).
#'
#'@examples
#'\donttest{clean_expr_data <- filter_counts()
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'simulated_counts <- simulate_hierarchicell(data_summaries)}
#'
#'@export

simulate_hierarchicell <- function(data_summaries,
                                   n_genes = 1000,
                                   n_per_group = 3,
                                   n_cases = n_per_group,
                                   n_controls = n_per_group,
                                   cells_per_control = 150,
                                   cells_per_case = 150,
                                   ncells_variation_type = "Poisson",
                                   foldchange = 1,
                                   decrease_dropout = 0,
                                   alter_dropout_cases = 0,
                                   tSNE_plot = FALSE){

   if ((tSNE_plot == TRUE)) {
    if (!requireNamespace("Seurat",quietly = TRUE)){
      stop("The 'Seurat' package is required for plotting. Please install it",
           call. = FALSE)
    } else {

      message("Computing simulation parameters ...")

      gene_mean_shape <- approximate_gene_mean(data_summaries)[1]
      gene_mean_rate <- approximate_gene_mean(data_summaries)[2]
      message("-------------------------------------------------------")
      message("Distribution of grand means is a gamma\nwith shape: ",
              round(gene_mean_shape,2),
              " and rate: ",
              round(gene_mean_rate,2))

      gene_dropout_shape <- approximate_gene_drop(data_summaries)[1]
      gene_dropout_rate <- approximate_gene_drop(data_summaries)[2]
      message("-------------------------------------------------------")
      message("Distribution for gene-wise dropout is a gamma \n with shape: ",
              round(gene_dropout_shape,2),
              " and rate: ",
              round(gene_dropout_rate,2))

      dropoutstd_beta0 <- model_drop_sd(data_summaries)[1]
      dropoutstd_beta1 <- model_drop_sd(data_summaries)[2]
      dropoutstd_beta2 <- model_drop_sd(data_summaries)[3]
      message("-------------------------------------------------------")
      message("Function for dropout SD is:\nDropoutStD = ",
              round(dropoutstd_beta0,2)," + ",round(dropoutstd_beta1,2),"*DropOut + ",
              round(dropoutstd_beta2,2),"*(DropOut**2)")

      inter_beta0 <- 0
      inter_beta1 <- model_inter(data_summaries)[1]

      message("-------------------------------------------------------")
      message("Function for inter-individual SD is:\nInterStDev = ",
              round(inter_beta0,2)," + ",
              round(inter_beta1,2),"*GrandMean)")

      dispersion_beta0 <- model_dispersion(data_summaries)[1]
      dispersion_beta1 <- model_dispersion(data_summaries)[2]
      message("-------------------------------------------------------")
      message("Function for dispersion is:\n exp(",
              round(dispersion_beta0,2),
              " + ",
              round(dispersion_beta1,2),
              "/IntraMean)")


      message("-------------------------------------------------------")
      message("Simulating cells ...")

      if (ncells_variation_type == "Poisson") {
        ncells_per_control <- stats::rpois(n = n_controls, lambda = cells_per_control)
        ncells_per_case <- stats::rpois(n = n_cases, lambda = cells_per_case)
      } else if (ncells_variation_type == "Fixed") {
        ncells_per_control <- rep(times = n_controls, x = cells_per_control)
        ncells_per_case <- rep(times = n_cases, x = cells_per_case)
      } else if (ncells_variation_type == "NB") {
        ncells_per_control <- stats::rnbinom(n = n_controls, mu = cells_per_control, size = 1)
        ncells_per_case <- stats::rnbinom(n = n_cases, mu = cells_per_case, size = 1)
      } else {
        stop("The variation type you selected for the number
             of cells per individual is not properly specified.
             Please correct")
      }

      message("-------------------------------------------------------")
      message("Simulating expression values ... ")

      allcells <- NULL

      simulate_gene <- function(){

        grandmean <- stats::rgamma(n=1,shape=gene_mean_shape,rate=gene_mean_rate)
        stddev_of_within_means <- inter_beta1*grandmean
        fc <- ifelse(stats::rbinom(n=1, size=1, prob = 0.5) == 1, foldchange, 1/foldchange)
        prob_zero <- stats::rgamma(n=1,shape=gene_dropout_shape,rate=gene_dropout_rate)
        prob_zero <- ifelse(prob_zero > 1, stats::rgamma(n=1,shape=gene_dropout_shape,rate=gene_dropout_rate), prob_zero)
        drop.sd <- dropoutstd_beta0 + dropoutstd_beta1*prob_zero + dropoutstd_beta2*(prob_zero**2)
        drop.sd <- ifelse(drop.sd < 0, 0, drop.sd)
        prob_zero <- rnorm(n=1,mean = prob_zero, sd = drop.sd)
        prob_zero <- ifelse(prob_zero < 0, 0, prob_zero)
        prob_zero <- ifelse(prob_zero > 1, 1, prob_zero)
        prob_zero <- 1 - prob_zero

        for (i in 1:n_controls){

          controlmean <- grandmean + stats::rnorm(n=1,mean=0,sd=stddev_of_within_means)
          controlmean <- ifelse(controlmean < 0, 0.0000001, controlmean)
          control_size <- exp(dispersion_beta0 + (dispersion_beta1/controlmean))
          controlcells <- stats::rnbinom(n=ncells_per_control[i],mu=controlmean,size=control_size)
          controlcells <- ifelse(stats::rbinom(n=length(controlcells),size=1,prob=prob_zero)==1, 0, controlcells)
          names(controlcells) <- paste0("Control_",i,"_Cell_",1:ncells_per_control[i])
          allcells <- c(allcells,controlcells)
        }

        for (i in 1:n_cases){


          casemean <- (grandmean*fc) + stats::rnorm(n=1,mean=0,sd=stddev_of_within_means)
          casemean <- ifelse(casemean < 0, 0.0000001, casemean)
          case_size <- exp(dispersion_beta0 + (dispersion_beta1/casemean))
          casecells <- stats::rnbinom(n=ncells_per_case[i],mu=casemean,size=case_size)
          casecells <- ifelse(stats::rbinom(n=length(casecells),size=1,prob=prob_zero)==1, 0, casecells)
          names(casecells) <- paste0("Case_",i,"_Cell_",1:ncells_per_case[i])
          allcells <- c(allcells,casecells)
        }

        allcells
      }
      all_genes <- as.data.frame(replicate(n_genes,simulate_gene()))
      colnames(all_genes) <- paste0("Gene",1:n_genes)
      all_genes <- data.frame(all_genes)
      all_genes$ToSep <- rownames(all_genes)
      all_genes$wellKey <- rownames(all_genes)
      all_genes <- tidyr::separate(all_genes,ToSep,
                            c("Status", "Donor_Number", "Cell", "Cell_Number"), sep="_")
      all_genes$Cell_Number <- paste0("Cell_", all_genes$Cell_Number)
      all_genes$DonorID <- paste0(all_genes$Status, "_", all_genes$Donor_Number)
      all_genes <- all_genes[ ,c((n_genes + 5), (n_genes + 6), (n_genes + 1), 1:n_genes)]
      all_genes <- all_genes[which(apply(all_genes[,c(-1,-2,-3)],1,mean) > 0),]

      message("-------------------------------------------------------")
      message("Generating tSNE plot ...")

      counts <- stats::na.omit(as.matrix(t(all_genes[ ,-1:-3])))
      pheno <- all_genes[ ,1:3]
      all <- Seurat::CreateSeuratObject(counts, project="All_Cells", min.cells=3)
      rownames(pheno) <- pheno[ ,1]
      pheno <- pheno[ ,-1]
      pheno$Status <- as.factor(pheno$Status)
      pheno$DonorID <- as.factor(pheno$DonorID)
      pheno <- pheno[colnames(all), ]
      all(rownames(pheno) %in% colnames(all))
      all(rownames(pheno) == colnames(all))
      all <- Seurat::AddMetaData(object = all, metadata = pheno)
      all <- Seurat::NormalizeData(all)
      all <- Seurat::FindVariableFeatures(all,do.plot=F)
      all <- Seurat::ScaleData(all)
      all <- Seurat::RunPCA(all)
      all <- Seurat::FindNeighbors(all)
      all <- Seurat::FindClusters(all)
      all <- Seurat::RunTSNE(all)

      print(Seurat::DimPlot(object = all,reduction = "tsne",group.by="DonorID"))

    }

  } else {

    message("Computing simulation parameters ...")

    gene_mean_shape <- approximate_gene_mean(data_summaries)[1]
    gene_mean_rate <- approximate_gene_mean(data_summaries)[2]
    message("-------------------------------------------------------")
    message("Distribution of grand means is a gamma\nwith shape: ",
            round(gene_mean_shape,2),
            " and rate: ",
            round(gene_mean_rate,2))

    gene_dropout_shape <- approximate_gene_drop(data_summaries)[1]
    gene_dropout_rate <- approximate_gene_drop(data_summaries)[2]
    message("-------------------------------------------------------")
    message("Distribution for gene-wise dropout is a gamma \n with shape: ",
            round(gene_dropout_shape,2),
            " and rate: ",
            round(gene_dropout_rate,2))

    dropoutstd_beta0 <- model_drop_sd(data_summaries)[1]
    dropoutstd_beta1 <- model_drop_sd(data_summaries)[2]
    dropoutstd_beta2 <- model_drop_sd(data_summaries)[3]
    message("-------------------------------------------------------")
    message("Function for dropout SD is:\nDropoutStD = ",
            round(dropoutstd_beta0,2)," + ",round(dropoutstd_beta1,2),"*DropOut + ",
            round(dropoutstd_beta2,2),"*(DropOut**2)")

    inter_beta0 <- 0
    inter_beta1 <- model_inter(data_summaries)[1]

    message("-------------------------------------------------------")
    message("Function for inter-individual SD is:\nInterStDev = ",
            round(inter_beta0,2)," + ",
            round(inter_beta1,2),"*GrandMean)")

    dispersion_beta0 <- model_dispersion(data_summaries)[1]
    dispersion_beta1 <- model_dispersion(data_summaries)[2]
    message("-------------------------------------------------------")
    message("Function for dispersion is:\n exp(",
            round(dispersion_beta0,2),
            " + ",
            round(dispersion_beta1,2),
            "/IntraMean)")


    message("-------------------------------------------------------")
    message("Simulating cells ...")

    if (ncells_variation_type == "Poisson") {
      ncells_per_control <- stats::rpois(n = n_controls, lambda = cells_per_control)
      ncells_per_case <- stats::rpois(n = n_cases, lambda = cells_per_case)
    } else if (ncells_variation_type == "Fixed") {
      ncells_per_control <- rep(times = n_controls, x = cells_per_control)
      ncells_per_case <- rep(times = n_cases, x = cells_per_case)
    } else if (ncells_variation_type == "NB") {
      ncells_per_control <- stats::rnbinom(n = n_controls, mu = cells_per_control, size = 1)
      ncells_per_case <- stats::rnbinom(n = n_cases, mu = cells_per_case, size = 1)
    } else {
      stop("The variation type you selected for the number
             of cells per individual is not properly specified.
             Please correct")
    }

    message("-------------------------------------------------------")
    message("Simulating expression values ... ")

    allcells <- NULL

    simulate_gene <- function(){

      grandmean <- stats::rgamma(n=1,shape=gene_mean_shape,rate=gene_mean_rate)
      stddev_of_within_means <- inter_beta1*grandmean
      fc <- ifelse(stats::rbinom(n=1, size=1, prob = 0.5) == 1, foldchange, 1/foldchange)
      prob_zero <- stats::rgamma(n=1,shape=gene_dropout_shape,rate=gene_dropout_rate)
      prob_zero <- ifelse(prob_zero > 1, stats::rgamma(n=1,shape=gene_dropout_shape,rate=gene_dropout_rate), prob_zero)
      drop.sd <- dropoutstd_beta0 + dropoutstd_beta1*prob_zero + dropoutstd_beta2*(prob_zero**2)
      drop.sd <- ifelse(drop.sd < 0, 0, drop.sd)
      prob_zero <- rnorm(n=1,mean = prob_zero, sd = drop.sd)
      prob_zero <- ifelse(prob_zero < 0, 0, prob_zero)
      prob_zero <- ifelse(prob_zero > 1, 1, prob_zero)
      prob_zero <- 1 - prob_zero

      for (i in 1:n_controls){

        controlmean <- grandmean + stats::rnorm(n=1,mean=0,sd=stddev_of_within_means)
        controlmean <- ifelse(controlmean < 0, 0.0000001, controlmean)
        control_size <- exp(dispersion_beta0 + (dispersion_beta1/controlmean))
        controlcells <- stats::rnbinom(n=ncells_per_control[i],mu=controlmean,size=control_size)
        controlcells <- ifelse(stats::rbinom(n=length(controlcells),size=1,prob=prob_zero)==1, 0, controlcells)
        names(controlcells) <- paste0("Control_",i,"_Cell_",1:ncells_per_control[i])
        allcells <- c(allcells,controlcells)
      }

      for (i in 1:n_cases){


        casemean <- (grandmean*fc) + stats::rnorm(n=1,mean=0,sd=stddev_of_within_means)
        casemean <- ifelse(casemean < 0, 0.0000001, casemean)
        case_size <- exp(dispersion_beta0 + (dispersion_beta1/casemean))
        casecells <- stats::rnbinom(n=ncells_per_case[i],mu=casemean,size=case_size)
        casecells <- ifelse(stats::rbinom(n=length(casecells),size=1,prob=prob_zero)==1, 0, casecells)
        names(casecells) <- paste0("Case_",i,"_Cell_",1:ncells_per_case[i])
        allcells <- c(allcells,casecells)
      }

      allcells
    }

    all_genes <- as.data.frame(replicate(n_genes,simulate_gene()))
    colnames(all_genes) <- paste0("Gene",1:n_genes)
    all_genes <- data.frame(all_genes)
    all_genes$ToSep <- rownames(all_genes)
    all_genes$wellKey <- rownames(all_genes)
    all_genes <- tidyr::separate(all_genes,ToSep,
                                 c("Status", "Donor_Number", "Cell", "Cell_Number"), sep="_")
    all_genes$Cell_Number <- paste0("Cell_", all_genes$Cell_Number)
    all_genes$DonorID <- paste0(all_genes$Status, "_", all_genes$Donor_Number)
    all_genes <- all_genes[ ,c((n_genes + 5), (n_genes + 6), (n_genes + 1), 1:n_genes)]
    all_genes <- all_genes[which(apply(all_genes[,c(-1,-2,-3)],1,mean) > 0),]


  }

  message("-------------------------------------------------------")
  message("All done!")
  as.data.frame(all_genes)

  }


#'@title Simulate Expression Data for a Continuous Measure
#'
#'@rdname simulate_hierarchicell_continuous
#'
#'@description This function will compute a  simulation that will borrow
#'  information from the input data (or the package default data) to simulate
#'  data under a variety of pre-determined conditions. These conditions include
#'  correlation between fold change and the continuous measure of interest,
#'  number of genes, number of samples (i.e., independent experimental units),
#'  and the mean number of cells per individual. The simulation incorporates
#'  information about the cell-wise dropout rates and library sizes from the
#'  unnormalized data and the gene-wise grand means, gene-wise dropout rates,
#'  inter-individual variance, and intra-individual variance from the normalized
#'  data.
#'
#'@details Prior to running the \code{\link{simulate_hierarchicell}} function,
#'  it is important to run the \code{\link{filter_counts}} function followed by
#'  the \code{\link{compute_data_summaries}} function to build an R object  that
#'  is in the right format for the following simulation function to properly
#'  work.
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
#'@param n_individuals an integer. The number of independent samples for
#'  simulation. If not specifying a foldchange, the number of cases and controls
#'  does not matter. Defaults to 3.
#'
#'@param cells_per_individual an integer. The mean number of cells per control
#'  you would like to simulate. Too large of a number may cause memory failure
#'  and may slow the simulation down tremendously. We recommend an integer less
#'  than 300, but more is possible. We note that anything greater than 100,
#'  brings marginal improvements in power. Defaults to 100.
#'
#'@param ncells_variation_type either "Poisson", "NB", or "Fixed". Allows the
#'  number of cells per individual to be fixed at exactly the specified number
#'  of cells per individual, vary slightly with a poisson distribution with a
#'  lambda equal to the specified number of cells per individual, or a negative
#'  binomial with a mean equal to the specified number of cells and dispersion
#'  size equal to one.Defaults to "Poisson".
#'
#'
#'@param rho a number between -1 and 1. The amount of correlation between fold
#'  change and the continuous measure of interest.Defaults to 1.
#'
#'@param continuous_mean A number. The mean for your continuous measure of
#'  interest. Assumes a normal distribution.Defaults to 0.
#'
#'@param continuous_sd A number. The standard deviation for your continuous
#'  measure of interest. Assumes a normal distribution.Defaults to 1.
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
#'@param tSNE_plot a TRUE/FALSE statement for the output of a tSNE plot to
#'  observe the global behavior of your simulated data. Seurat will need to be
#'  installed for this function to properly work. Defaults to FALSE.
#'
#'@return A data.frame of the simulated data or potentially a pdf of a tSNE plot
#'  (if tSNE_plot=TRUE).
#'
#'@examples
#'\donttest{clean_expr_data <- filter_counts()
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'simulated_counts <- simulate_hierarchicell_continuous(data_summaries)}
#'
#'@export

simulate_hierarchicell_continuous <- function(data_summaries,
                                   n_genes = 1000,
                                   n_individuals = 3,
                                   cells_per_individual = 100,
                                   ncells_variation_type = "Poisson",
                                   rho = 1,
                                   continuous_mean = 0,
                                   continuous_sd = 1,
                                   decrease_dropout = 0,
                                   tSNE_plot = FALSE){

    if ((tSNE_plot == TRUE)) {
      if (!requireNamespace("Seurat",quietly = TRUE)){
        stop("The 'Seurat' package is required for plotting. Please install it",
             call. = FALSE)
      } else {

        message("Computing simulation parameters ...")

        gene_mean_shape <- approximate_gene_mean(data_summaries)[1]
        gene_mean_rate <- approximate_gene_mean(data_summaries)[2]
        message("-------------------------------------------------------")
        message("Distribution of grand means is a gamma\nwith shape: ",
                round(gene_mean_shape,2),
                " and rate: ",
                round(gene_mean_rate,2))

        gene_dropout_shape <- approximate_gene_drop(data_summaries)[1]
        gene_dropout_rate <- approximate_gene_drop(data_summaries)[2]
        message("-------------------------------------------------------")
        message("Distribution for gene-wise dropout is a gamma \n with shape: ",
                round(gene_dropout_shape,2),
                " and rate: ",
                round(gene_dropout_rate,2))

        dropoutstd_beta0 <- model_drop_sd(data_summaries)[1]
        dropoutstd_beta1 <- model_drop_sd(data_summaries)[2]
        dropoutstd_beta2 <- model_drop_sd(data_summaries)[3]
        message("-------------------------------------------------------")
        message("Function for dropout SD is:\nDropoutStD = ",
                round(dropoutstd_beta0,2)," + ",round(dropoutstd_beta1,2),"*DropOut + ",
                round(dropoutstd_beta2,2),"*(DropOut**2)")

        inter_beta0 <- 0
        inter_beta1 <- model_inter(data_summaries)[1]

        message("-------------------------------------------------------")
        message("Function for inter-individual SD is:\nInterStDev = ",
                round(inter_beta0,2)," + ",
                round(inter_beta1,2),"*GrandMean)")

        dispersion_beta0 <- model_dispersion(data_summaries)[1]
        dispersion_beta1 <- model_dispersion(data_summaries)[2]
        message("-------------------------------------------------------")
        message("Function for dispersion is:\n exp(",
                round(dispersion_beta0,2),
                " + ",
                round(dispersion_beta1,2),
                "/IntraMean)")


        message("-------------------------------------------------------")
        message("Simulating cells ...")

        if (ncells_variation_type == "Poisson") {
          ncells_per_individual <- stats::rpois(n = n_individuals, lambda = cells_per_individual)
        } else if (ncells_variation_type == "Fixed") {
          ncells_per_individual <- rep(times = n_individuals, x = cells_per_individual)
        } else if (ncells_variation_type == "NB") {
          ncells_per_individual <- stats::rnbinom(n = n_individuals, mu = cells_per_individual, size = 1)
        } else {
          stop("The variation type you selected for the number
             of cells per individual is not properly specified.
             Please correct")
        }

        message("-------------------------------------------------------")
        message("Simulating expression values ... ")

        allcells <- NULL
        outcome <- rnorm(n_individuals, mean = continuous_mean, sd = continuous_sd)
        mergeoutcome <- as.data.frame(cbind(paste0("Individual_",1:n_individuals),outcome))
        colnames(mergeoutcome) <- c("DonorID","Outcome")
        scaledoutcome <- scale(outcome)

        complement <- function(y, rho) {
          x <- rnorm(length(y))
          y.perp <- stats::residuals(stats::lm(x ~ y))
          rho * stats::sd(y.perp) * y + y.perp * stats::sd(y) * sqrt(1 - rho^2)
        }

        simulate_gene <- function(){

          grandmean <- stats::rgamma(n=1,shape=gene_mean_shape,rate=gene_mean_rate)
          stddev_of_within_means <- inter_beta1*grandmean
          log2foldchange <- complement(scaledoutcome,rho)
          foldchange <- 2**log2foldchange
          prob_zero <- stats::rgamma(n=1,shape=gene_dropout_shape,rate=gene_dropout_rate)
          prob_zero <- ifelse(prob_zero > 1, stats::rgamma(n=1,shape=gene_dropout_shape,rate=gene_dropout_rate), prob_zero)
          drop.sd <- dropoutstd_beta0 + dropoutstd_beta1*prob_zero + dropoutstd_beta2*(prob_zero**2)
          drop.sd <- ifelse(drop.sd < 0, 0, drop.sd)
          prob_zero <- rnorm(n=1,mean = prob_zero, sd = drop.sd)
          prob_zero <- ifelse(prob_zero < 0, 0, prob_zero)
          prob_zero <- ifelse(prob_zero > 1, 1, prob_zero)
          prob_zero <- 1 - prob_zero

          for (i in 1:n_individuals){

            indmean <- grandmean + stats::rnorm(n=1,mean=0,sd=stddev_of_within_means)
            indmean <- ifelse(indmean < 0, 0.0000001, indmean)
            indmean <- indmean*foldchange[i]
            ind_size <- exp(dispersion_beta0 + (dispersion_beta1/indmean))
            indcells <- stats::rnbinom(n=ncells_per_individual[i],mu=indmean,size=ind_size)
            indcells <- ifelse(stats::rbinom(n=length(indcells),size=1,prob=prob_zero)==1, 0, indcells)
            names(indcells) <- paste0("Individual_",i,"_Cell_",1:ncells_per_individual[i])
            allcells <- c(allcells,indcells)
          }

          allcells
        }
        all_genes <- as.data.frame(replicate(n_genes,simulate_gene()))
        colnames(all_genes) <- paste0("Gene",1:n_genes)
        all_genes <- data.frame(all_genes)
        all_genes$ToSep <- rownames(all_genes)
        all_genes$wellKey <- rownames(all_genes)
        all_genes <- tidyr::separate(all_genes,ToSep,
                                     c("Ind", "Donor_Number", "Cell", "Cell_Number"), sep="_")
        all_genes$Cell_Number <- paste0("Cell_", all_genes$Cell_Number)
        all_genes$DonorID <- paste0(all_genes$Ind, "_", all_genes$Donor_Number)
        all_genes <- merge(all_genes,mergeoutcome,by="DonorID")
        all_genes <- all_genes[ ,c(1,(n_genes + 6),(n_genes + 5),(n_genes + 7),2:(n_genes + 1))]
        rownames(all_genes) <- all_genes$wellKey

        counts <- na.omit(as.matrix(t(all_genes[,-1:-4])))
        pheno <- all_genes[,1:4]
        pheno$Outcome <- as.numeric(as.character(pheno$Outcome))

        message("-------------------------------------------------------")
        message("Generating tSNE plot ...")
        #Create a Seurat Dataset

        all <- Seurat::CreateSeuratObject(counts, project="All_Cells", min.cells=3)
        rownames(pheno) <- pheno[,2]
        pheno <- pheno[,-2]
        pheno <- pheno[colnames(all),]
        all(rownames(pheno) %in% colnames(all))
        all(rownames(pheno) == colnames(all))
        all <- Seurat::AddMetaData(object = all, metadata = pheno)
        all <- Seurat::NormalizeData(all)
        all <- Seurat::FindVariableFeatures(all,do.plot=F)
        all <- Seurat::ScaleData(all)
        all <- Seurat::RunPCA(all)
        all <- Seurat::FindNeighbors(all)
        all <- Seurat::FindClusters(all)
        all <- Seurat::RunTSNE(all, check_duplicates=F)

        print(Seurat::FeaturePlot(object = all,features = "Outcome"))

      }

    } else {

      message("Computing simulation parameters ...")

      gene_mean_shape <- approximate_gene_mean(data_summaries)[1]
      gene_mean_rate <- approximate_gene_mean(data_summaries)[2]
      message("-------------------------------------------------------")
      message("Distribution of grand means is a gamma\nwith shape: ",
              round(gene_mean_shape,2),
              " and rate: ",
              round(gene_mean_rate,2))

      gene_dropout_shape <- approximate_gene_drop(data_summaries)[1]
      gene_dropout_rate <- approximate_gene_drop(data_summaries)[2]
      message("-------------------------------------------------------")
      message("Distribution for gene-wise dropout is a gamma \n with shape: ",
              round(gene_dropout_shape,2),
              " and rate: ",
              round(gene_dropout_rate,2))

      dropoutstd_beta0 <- model_drop_sd(data_summaries)[1]
      dropoutstd_beta1 <- model_drop_sd(data_summaries)[2]
      dropoutstd_beta2 <- model_drop_sd(data_summaries)[3]
      message("-------------------------------------------------------")
      message("Function for dropout SD is:\nDropoutStD = ",
              round(dropoutstd_beta0,2)," + ",round(dropoutstd_beta1,2),"*DropOut + ",
              round(dropoutstd_beta2,2),"*(DropOut**2)")

      inter_beta0 <- 0
      inter_beta1 <- model_inter(data_summaries)[1]

      message("-------------------------------------------------------")
      message("Function for inter-individual SD is:\nInterStDev = ",
              round(inter_beta0,2)," + ",
              round(inter_beta1,2),"*GrandMean)")

      dispersion_beta0 <- model_dispersion(data_summaries)[1]
      dispersion_beta1 <- model_dispersion(data_summaries)[2]
      message("-------------------------------------------------------")
      message("Function for dispersion is:\n exp(",
              round(dispersion_beta0,2),
              " + ",
              round(dispersion_beta1,2),
              "/IntraMean)")
      message("-------------------------------------------------------")
      message("Simulating cells ...")

      if (ncells_variation_type == "Poisson") {
        ncells_per_individual <- stats::rpois(n = n_individuals, lambda = cells_per_individual)
      } else if (ncells_variation_type == "Fixed") {
        ncells_per_individual <- rep(times = n_individuals, x = cells_per_individual)
      } else if (ncells_variation_type == "NB") {
        ncells_per_individual <- stats::rnbinom(n = n_individuals, mu = cells_per_individual, size = 1)
      } else {
        stop("The variation type you selected for the number
             of cells per individual is not properly specified.
             Please correct")
      }

      message("-------------------------------------------------------")
      message("Simulating expression values ... ")

      allcells <- NULL
      outcome <- rnorm(n_individuals, mean = continuous_mean, sd = continuous_sd)
      mergeoutcome <- as.data.frame(cbind(paste0("Individual_",1:n_individuals),outcome))
      colnames(mergeoutcome) <- c("DonorID","Outcome")
      scaledoutcome <- scale(outcome)

      complement <- function(y, rho) {
        x <- rnorm(length(y))
        y.perp <- stats::residuals(stats::lm(x ~ y))
        rho * stats::sd(y.perp) * y + y.perp * stats::sd(y) * sqrt(1 - rho^2)
      }

      simulate_gene <- function(){

        grandmean <- stats::rgamma(n=1,shape=gene_mean_shape,rate=gene_mean_rate)
        stddev_of_within_means <- inter_beta1*grandmean
        log2foldchange <- complement(scaledoutcome,rho)
        foldchange <- 2**log2foldchange
        prob_zero <- stats::rgamma(n=1,shape=gene_dropout_shape,rate=gene_dropout_rate)
        prob_zero <- ifelse(prob_zero > 1, stats::rgamma(n=1,shape=gene_dropout_shape,rate=gene_dropout_rate), prob_zero)
        drop.sd <- dropoutstd_beta0 + dropoutstd_beta1*prob_zero + dropoutstd_beta2*(prob_zero**2)
        drop.sd <- ifelse(drop.sd < 0, 0, drop.sd)
        prob_zero <- rnorm(n=1,mean = prob_zero, sd = drop.sd)
        prob_zero <- ifelse(prob_zero < 0, 0, prob_zero)
        prob_zero <- ifelse(prob_zero > 1, 1, prob_zero)
        prob_zero <- 1 - prob_zero

        for (i in 1:n_individuals){

          indmean <- grandmean + stats::rnorm(n=1,mean=0,sd=stddev_of_within_means)
          indmean <- ifelse(indmean < 0, 0.0000001, indmean)
          indmean <- indmean*foldchange[i]
          ind_size <- exp(dispersion_beta0 + (dispersion_beta1/indmean))
          indcells <- stats::rnbinom(n=ncells_per_individual[i],mu=indmean,size=ind_size)
          indcells <- ifelse(stats::rbinom(n=length(indcells),size=1,prob=prob_zero)==1, 0, indcells)
          names(indcells) <- paste0("Individual_",i,"_Cell_",1:ncells_per_individual[i])
          allcells <- c(allcells,indcells)
        }

        allcells
      }
      all_genes <- as.data.frame(replicate(n_genes,simulate_gene()))
      colnames(all_genes) <- paste0("Gene",1:n_genes)
      all_genes <- data.frame(all_genes)
      all_genes$ToSep <- rownames(all_genes)
      all_genes$wellKey <- rownames(all_genes)
      all_genes <- tidyr::separate(all_genes,ToSep,
                                   c("Ind", "Donor_Number", "Cell", "Cell_Number"), sep="_")
      all_genes$Cell_Number <- paste0("Cell_", all_genes$Cell_Number)
      all_genes$DonorID <- paste0(all_genes$Ind, "_", all_genes$Donor_Number)
      all_genes <- merge(all_genes,mergeoutcome,by="DonorID")
      all_genes <- all_genes[ ,c(1,(n_genes + 6),(n_genes + 5),(n_genes + 7),2:(n_genes + 1))]
      rownames(all_genes) <- all_genes$wellKey


    }

    message("-------------------------------------------------------")
    message("All done!")
    as.data.frame(all_genes)

  }
