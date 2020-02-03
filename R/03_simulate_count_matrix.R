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
#'   the \code{\link{compute_data_summaries}} function to build an R object
#'   that is in the right format for the following simulation function to
#'   properly work.
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
#'  Defaults to 15,000.
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
#'  an integer less than 1,000.Defaults to 150.
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
#'@param tSNE_plot a TRUE/FALSE statement for the output of a tSNE plot to
#'  observe the global behavior of your simulated data. Seurat will need to be
#'  installed for this function to properly work. Defaults to FALSE.
#'
#'@param plot_name a character string preferrably containing no spaces for the
#'  name of the pdf file where the tSNE plot will be stored. Defaults to
#'  "tSNE_Simulated_Data".
#'
#'@return A data.frame of the simulated data or potentially a pdf of a tSNE
#'  plot (if tSNE_plot=TRUE).
#'
#'@examples
#'clean_expr_data <- filter_counts()
#'data_summaries <- compute_data_summaries(clean_expr_data)
#'simulated_counts <- simulate_hierarchicell(data_summaries)
#'
#'@export

simulate_hierarchicell <- function(data_summaries,
                                   n_genes = 15000,
                                   n_per_group = 3,
                                   n_cases = n_per_group,
                                   n_controls = n_per_group,
                                   cells_per_individual = 150,
                                   foldchange = 1,
                                   decrease_dropout = 0,
                                   tSNE_plot = FALSE,
                                   plot_name="tSNE_Simulated_Data"){

  if ((tSNE_plot == TRUE)) {
    if (!requireNamespace("Seurat",quietly = TRUE)){
      stop("The 'Seurat' package is required for plotting. Please install it",
           call. = FALSE)
    } else {

      message("Computing simulation parameters ...")

      library_mean <- approximate_library_sizes(data_summaries)[1]
      library_sd <- approximate_library_sizes(data_summaries)[2]
      message("-------------------------------------------------------")
      message("Distribution of library sizes is a lognormal\nwith mean: ",
              round(library_mean,2),
              " and variance:",
              round(library_sd**2,2))

      cell_dropout_mean <- approximate_cell_dropout(data_summaries)[1]
      cell_dropout_sd <- approximate_cell_dropout(data_summaries)[2]
      message("-------------------------------------------------------")
      message("Distribution of cell-wise dropout is a normal\nwith mean: ",
              round(cell_dropout_mean,2),
              " and variance: ",
              round(cell_dropout_sd**2,2))

      gene_mean_shape <- approximate_gene_mean(data_summaries)[1]
      gene_mean_rate <- approximate_gene_mean(data_summaries)[2]
      message("-------------------------------------------------------")
      message("Distribution of grand means is a gamma\nwith shape: ",
              round(gene_mean_shape,2),
              " and rate: ",
              round(gene_mean_rate,2))

      gene_dropout_beta1 <- model_gene_drop(data_summaries)[2]
      message("-------------------------------------------------------")
      message("Function for gene-wise dropout is:\nDropout = 1 + (",
              round(gene_dropout_beta1,2),")(GrandMean)")

      inter_beta1 <- model_inter(data_summaries)[1]
      inter_beta2 <- model_inter(data_summaries)[2]
      message("-------------------------------------------------------")
      message("Function for inter-individual standard deviation is:\nInterStDev = 0 + (",
              round(inter_beta1,2),")(GrandMean) + (",
              round(inter_beta2,2),")(GrandMeanSq)")

      intra_beta0 <- model_intra(data_summaries)[1]
      intra_beta1 <- model_intra(data_summaries)[2]
      intra_beta2 <- model_intra(data_summaries)[3]
      message("-------------------------------------------------------")
      message("Function for intra-individual standard deviation is:\nIntraStDev = ",
              round(intra_beta0,2)," + (",
              round(intra_beta1,2),")(GrandMean) + (",
              round(intra_beta2,2),")(GrandMeanSq)")
      message("-------------------------------------------------------")
      message("Simulating cells ...")

      ncells_per_control <- stats::rpois(n = n_controls, lambda = cells_per_individual)
      ncells_per_case <- stats::rpois(n = n_cases, lambda = cells_per_individual)
      ncells <- sum(ncells_per_case) + sum(ncells_per_control)
      prob_zero_cell_control <- vector(mode = "list", length = n_controls)
      prob_zero_cell_case <- vector(mode="list", length = n_cases)
      message("-------------------------------------------------------")
      message("Simulating cell dropout probabilities ...")

      for (i in 1:n_controls){
        prob_zero_cell_control[[i]] <- stats::rnorm(n=ncells_per_control[i],
                                             mean=(cell_dropout_mean*(1-decrease_dropout)),
                                             sd=cell_dropout_sd)
        prob_zero_cell_control[[i]] <- ifelse(prob_zero_cell_control[[i]] < 0,
                                              0, prob_zero_cell_control[[i]])
        prob_zero_cell_control[[i]] <- ifelse(prob_zero_cell_control[[i]] > 1,
                                              1, prob_zero_cell_control[[i]])
      }

      for (i in 1:n_cases){
        prob_zero_cell_case[[i]] <- stats::rnorm(n=ncells_per_case[i],
                                          mean=(cell_dropout_mean*(1-decrease_dropout)),
                                          sd=cell_dropout_sd)
        prob_zero_cell_case[[i]] <- ifelse(prob_zero_cell_case[[i]] < 0,
                                           0, prob_zero_cell_case[[i]])
        prob_zero_cell_case[[i]] <- ifelse(prob_zero_cell_case[[i]] > 1,
                                           1, prob_zero_cell_case[[i]])
      }
      message("-------------------------------------------------------")
      message("Simulating expression values ... ")

      allcells <- NULL

      simulate_gene <- function(){

        grandmean <- stats::rgamma(n=1,shape=gene_mean_shape,rate=gene_mean_rate)

        stddev_of_within_means <- (inter_beta1*grandmean) + (inter_beta2*(grandmean**2))
        stddev_of_within_means <- ifelse(stddev_of_within_means < 0, 0, stddev_of_within_means)
        within_donor_stddev <- intra_beta0 + (intra_beta1*grandmean) + (intra_beta2*(grandmean**2))
        within_donor_stddev <- ifelse(within_donor_stddev < 0, 0, within_donor_stddev)

        for (i in 1:n_controls){

          prob_zero_gene <- 1 + (gene_dropout_beta1*grandmean)
          prob_zero_gene <- ifelse(prob_zero_gene < 0, 0, prob_zero_gene)

          controlmean <- grandmean + stats::rnorm(n=1,mean=0,sd=stddev_of_within_means)
          controlcells <- stats::rnorm(n=ncells_per_control[i],mean=controlmean,sd=within_donor_stddev)
          prob_zero_control <- prob_zero_cell_control[[i]]*prob_zero_gene
          controlcells <- ifelse(stats::rbinom(n=length(controlcells),size=1,prob=prob_zero_control)==1, 0, controlcells)
          names(controlcells) <- paste0("Control_",i,"_Cell_",1:ncells_per_control[i])
          allcells <- c(allcells,controlcells)
        }

        for (i in 1:n_cases){

          prob_zero_gene <- 1 + (gene_dropout_beta1*grandmean)
          prob_zero_gene <- ifelse(prob_zero_gene < 0, 0, prob_zero_gene)

          fc <- ifelse(stats::rbinom(n=1, size=1, prob = 0.5) == 1, foldchange, 1/foldchange)
          casemean <- (grandmean*fc) + stats::rnorm(n=1,mean=0,sd=stddev_of_within_means)
          casecells <- stats::rnorm(n=ncells_per_case[i],mean=casemean,sd=within_donor_stddev)
          prob_zero_case <- prob_zero_cell_case[[i]]*prob_zero_gene
          casecells <- ifelse(stats::rbinom(n=length(casecells),size=1,prob=prob_zero_case)==1, 0, casecells)
          names(casecells) <- paste0("Case_",i,"_Cell_",1:ncells_per_case[i])
          allcells <- c(allcells,casecells)
        }

        allcells[allcells < 0] <- 0
        allcells
      }
      all_genes <- as.data.frame(replicate(n_genes,simulate_gene()))
      rn <- rownames(all_genes)

      message("-------------------------------------------------------")
      message("Simulating cell library sizes ... ")

      all_genes <- t(apply(all_genes, 1,
                           function(a){
                             l <- length(a); b <- stats::rgamma(l,a); return(b/sum(b))
                             }))
      all_genes <- t(apply(all_genes, 1,
                           function(b){
                             round((stats::rlnorm(n=1,meanlog = library_mean,sdlog = library_sd))*b,0)
                             }))
      all_genes[is.na(all_genes)] <- 0
      rownames(all_genes) <- rn
      colnames(all_genes) <- paste0("Gene",1:n_genes)
      all_genes <- data.frame(all_genes)
      all_genes$ToSep <- rownames(all_genes)
      all_genes$wellKey <- rownames(all_genes)
      all_genes <- tidyr::separate(all_genes,ToSep,
                            c("Status", "Donor_Number", "Cell", "Cell_Number"), sep="_")
      all_genes$Cell_Number <- paste0("Cell_", all_genes$Cell_Number)
      all_genes$DonorID <- paste0(all_genes$Status, "_", all_genes$Donor_Number)
      all_genes <- all_genes[ ,c((n_genes + 5), (n_genes + 6), (n_genes + 1), 1:n_genes)]

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

      grDevices::pdf(paste0(plot_name,".pdf"))
      print(Seurat::DimPlot(object = all,reduction = "tsne",group.by="DonorID"))
      grDevices::dev.off()
    }

  } else {

    message("Computing simulation parameters ...")

    library_mean <- approximate_library_sizes(data_summaries)[1]
    library_sd <- approximate_library_sizes(data_summaries)[2]
    message("-------------------------------------------------------")
    message("Distribution of library sizes is a lognormal\nwith mean: ",
            round(library_mean,2),
            " and StDev:",
            round(library_sd,2))

    cell_dropout_mean <- approximate_cell_dropout(data_summaries)[1]
    cell_dropout_sd <- approximate_cell_dropout(data_summaries)[2]
    message("-------------------------------------------------------")
    message("Distribution of cell-wise dropout is a normal\nwith mean: ",
            round(cell_dropout_mean,2),
            " and StDev: ",
            round(cell_dropout_sd,2))

    gene_mean_shape <- approximate_gene_mean(data_summaries)[1]
    gene_mean_rate <- approximate_gene_mean(data_summaries)[2]
    message("-------------------------------------------------------")
    message("Distribution of grand means is a gamma\nwith shape: ",
            round(gene_mean_shape,2),
            " and rate: ",
            round(gene_mean_rate,2))

    gene_dropout_beta1 <- model_gene_drop(data_summaries)[2]
    message("-------------------------------------------------------")
    message("Function for gene-wise dropout is:\nDropout = 1 + (",
            round(gene_dropout_beta1,2),")(GrandMean)")

    inter_beta1 <- model_inter(data_summaries)[1]
    inter_beta2 <- model_inter(data_summaries)[2]
    message("-------------------------------------------------------")
    message("Function for inter-individual standard deviation is:\nInterStDev = 0 + (",
            round(inter_beta1,2),")(GrandMean) + (",
            round(inter_beta2,2),")(GrandMeanSq)")

    intra_beta0 <- model_intra(data_summaries)[1]
    intra_beta1 <- model_intra(data_summaries)[2]
    intra_beta2 <- model_intra(data_summaries)[3]
    message("-------------------------------------------------------")
    message("Function for intra-individual standard deviation is:\nIntraStDev = ",
            round(intra_beta0,2)," + (",
            round(intra_beta1,2),")(GrandMean) + (",
            round(intra_beta2,2),")(GrandMeanSq)")
    message("-------------------------------------------------------")
    message("Simulating cells ...")

    ncells_per_control <- stats::rpois(n = n_controls, lambda = cells_per_individual)
    ncells_per_case <- stats::rpois(n = n_cases, lambda = cells_per_individual)
    ncells <- sum(ncells_per_case) + sum(ncells_per_control)
    prob_zero_cell_control <- vector(mode = "list", length = n_controls)
    prob_zero_cell_case <- vector(mode="list", length = n_cases)
    message("-------------------------------------------------------")
    message("Simulating cell dropout probabilities ...")

    for (i in 1:n_controls){
      prob_zero_cell_control[[i]] <- stats::rnorm(n=ncells_per_control[i],
                                           mean=(cell_dropout_mean*(1-decrease_dropout)),
                                           sd=cell_dropout_sd)
      prob_zero_cell_control[[i]] <- ifelse(prob_zero_cell_control[[i]] < 0,
                                            0, prob_zero_cell_control[[i]])
      prob_zero_cell_control[[i]] <- ifelse(prob_zero_cell_control[[i]] > 1,
                                            1, prob_zero_cell_control[[i]])
    }

    for (i in 1:n_cases){
      prob_zero_cell_case[[i]] <- stats::rnorm(n=ncells_per_case[i],
                                        mean=(cell_dropout_mean*(1-decrease_dropout)),
                                        sd=cell_dropout_sd)
      prob_zero_cell_case[[i]] <- ifelse(prob_zero_cell_case[[i]] < 0,
                                         0, prob_zero_cell_case[[i]])
      prob_zero_cell_case[[i]] <- ifelse(prob_zero_cell_case[[i]] > 1,
                                         1, prob_zero_cell_case[[i]])
    }
    message("-------------------------------------------------------")
    message("Simulating expression values ... ")

    allcells <- NULL

    simulate_gene <- function(){

      grandmean <- stats::rgamma(n=1,shape=gene_mean_shape,rate=gene_mean_rate)

      stddev_of_within_means <- (inter_beta1*grandmean) + (inter_beta2*(grandmean**2))
      stddev_of_within_means <- ifelse(stddev_of_within_means < 0, 0, stddev_of_within_means)
      within_donor_stddev <- intra_beta0 + (intra_beta1*grandmean) + (intra_beta2*(grandmean**2))
      within_donor_stddev <- ifelse(within_donor_stddev < 0, 0, within_donor_stddev)

      for (i in 1:n_controls){

        prob_zero_gene <- 1 + (gene_dropout_beta1*grandmean)
        prob_zero_gene <- ifelse(prob_zero_gene < 0, 0, prob_zero_gene)

        controlmean <- grandmean + stats::rnorm(n=1,mean=0,sd=stddev_of_within_means)
        controlcells <- stats::rnorm(n=ncells_per_control[i],mean=controlmean,sd=within_donor_stddev)
        prob_zero_control <- prob_zero_cell_control[[i]]*prob_zero_gene
        controlcells <- ifelse(stats::rbinom(n=length(controlcells),size=1,prob=prob_zero_control)==1, 0, controlcells)
        names(controlcells) <- paste0("Control_",i,"_Cell_",1:ncells_per_control[i])
        allcells <- c(allcells,controlcells)
      }

      for (i in 1:n_cases){

        prob_zero_gene <- 1 + (gene_dropout_beta1*grandmean)
        prob_zero_gene <- ifelse(prob_zero_gene < 0, 0, prob_zero_gene)

        casemean <- (grandmean*foldchange) + stats::rnorm(n=1,mean=0,sd=stddev_of_within_means)
        casecells <- stats::rnorm(n=ncells_per_case[i],mean=casemean,sd=within_donor_stddev)
        prob_zero_case <- prob_zero_cell_case[[i]]*prob_zero_gene
        casecells <- ifelse(stats::rbinom(n=length(casecells),size=1,prob=prob_zero_case)==1, 0, casecells)
        names(casecells) <- paste0("Case_",i,"_Cell_",1:ncells_per_case[i])
        allcells <- c(allcells,casecells)
      }

      allcells[allcells < 0] <- 0
      allcells
    }
    all_genes <- as.data.frame(replicate(n_genes,simulate_gene()))
    rn <- rownames(all_genes)

    message("-------------------------------------------------------")
    message("Simulating cell library sizes ... ")

    all_genes <- t(apply(all_genes, 1,
                         function(a){
                           l <- length(a); b <- stats::rgamma(l,a); return(b/sum(b))
                         }))
    all_genes <- t(apply(all_genes, 1,
                         function(b){
                           round((stats::rlnorm(n=1,meanlog = library_mean,sdlog = library_sd))*b,0)
                         }))
    all_genes[is.na(all_genes)] <- 0
    rownames(all_genes) <- rn
    colnames(all_genes) <- paste0("Gene",1:n_genes)
    all_genes <- data.frame(all_genes)
    all_genes$ToSep <- rownames(all_genes)
    all_genes$wellKey <- rownames(all_genes)
    all_genes <- tidyr::separate(all_genes,ToSep,
                          c("Status", "Donor_Number", "Cell", "Cell_Number"), sep="_")
    all_genes$Cell_Number <- paste0("Cell_", all_genes$Cell_Number)
    all_genes$DonorID <- paste0(all_genes$Status, "_", all_genes$Donor_Number)
    all_genes <- all_genes[ ,c((n_genes + 5), (n_genes + 6), (n_genes + 1), 1:n_genes)]


  }

  message("-------------------------------------------------------")
  message("All done!")
  as.data.frame(all_genes)
}


