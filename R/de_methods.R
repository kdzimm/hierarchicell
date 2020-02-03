#' @title Differential Expression Methods for Type 1 Error Computation
#'
#' @name de_methods
#'
#' @description The differential expression tools available for type 1 error
#'   computation. For all of these tools, a modest number of genes is required
#'   for stable estimation of the type 1 error and for appropriate distribution
#'   of the library size spread across many genes. (Too small a number of genes
#'   and the library size will be allocated all to a small number of genes,
#'   which would be inappropriate). Further description of each method can be
#'   seen below.
#'
#' @section MAST:
#'
#'   Computes type 1 error rates using MAST (Model-based Analysis of Single-cell
#'   Transcriptomics) as implemented in the vignette. MAST implements a hurdle
#'   model to account for the large number of zero values (dropout) and
#'   bi-modally distributed single-cell data. Prior to computing differential
#'   expression, the data are log-transformed. Implementing MAST without
#'   accounting for inter-individual heterogeneity, when it is present, will
#'   show inflated type 1 error rates. For more detailed information see:
#'   https://www.ncbi.nlm.nih.gov/pubmed/26653891
#'
#' @section MAST_RE:
#'
#'   Computes type 1 error rates using MAST (Model-based Analysis of Single-cell
#'   Transcriptomics) with random effects for individual. MAST implements a
#'   hurdle model to account for the large number of zero values (dropout) and
#'   bi-modally distributed single-cell data. Prior to computing differential
#'   expression, the data are log-transformed. Implementing MAST with random
#'   effects to account for inter-individual heterogeneity, when it is present,
#'   will show well-controlled type 1 error rates. For more detailed information
#'   see: https://www.ncbi.nlm.nih.gov/pubmed/26653891
#'
#'
#' @section MAST_Batch:
#'
#'   Computes type 1 error rates using MAST (Model-based Analysis of Single-cell
#'   Transcriptomics) after computing a batch effect correction for individual.
#'   MAST implements a hurdle model to account for the large number of zero
#'   values (dropout) and bi-modally distributed single-cell data. Prior to
#'   computing differential expression, the data are log-transformed.
#'   Implementing MAST with batch effect correction to account for
#'   inter-individual heterogeneity, when it is present, will show highly
#'   inflated type 1 error rates. Note: combat will not work if number of cells
#'   is low. For more detailed information see:
#'   https://www.ncbi.nlm.nih.gov/pubmed/26653891
#'
#' @section GLM_tweedie:
#'
#'   Computes type 1 error rates using generalized linear model assuming a
#'   tweedie distribution. The tweedie distribution will help to capture the
#'   large number of zero values (dropout)in single-cell data. Prior to
#'   computing differential expression, the data are normalized using an ALR
#'   transform. Implementing a tweedie GLM without accounting for
#'   inter-individual heterogeneity, when it is present, will show inflated type
#'   1 error rates. Note: to implement this with a large number of genes takes a
#'   good chunk of time.
#'
#' @section GLMM_tweedie:
#'
#'   Computes type 1 error rates using generalized linear mixed model assuming a
#'   tweedie distribution. The tweedie distribution will help to capture the
#'   large number of zero values (dropout)in single-cell data. Prior to
#'   computing differential expression, the data are normalized using an ALR
#'   transform. Implementing a tweedie GLMM using random effects to account for
#'   inter-individual heterogeneity, when it is present, will show
#'   well-controlled type 1 error rates. Note: to implement this with a large
#'   number of genes takes a good chunk of time.
#'
#' @section GEE1:
#'
#'   Computes type 1 error rates using generalized estimating equations (GEE1)
#'   with an exchangeable correlation structure. Prior to computing differential
#'   expression, the data are normalized using an ALR transform. Implementing a
#'   GEE1 accounting for inter-individual heterogeneity, when it is present,
#'   will show well-controlled type 1 error rates as the number of individuals
#'   in the study grows. Note: to implement this with a large number of genes
#'   takes a good chunk of time.
#'
#' @section ROTS:
#'
#'   Computes type 1 error rates using ROTS (Reproducibility-Optimized
#'   Statistical Testings) as implemented in the vignette. ROTS computes tests
#'   with a modified t-statistic that is adjusted according to inherent
#'   properties of the data. Implementing ROTS, when inter-individual
#'   heterogenetiy is present, shows inflated type 1 error rates. For more
#'   detailed information see:
#'   https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005562
#'
#'
#'
#' @section Monocle:
#'
#'   Computes type 1 error rates using monocle as implemented in the vignette.
#'   Monocle computes differential expression a generalized linear model
#'   assuming a negative binomial distribution. Prior to computing differential
#'   expression, the data are normalized using the DESeq median normalization
#'   method. Implementing monocle, when inter-individual heterogenetiy is
#'   present, shows inflated type 1 error rates. For more detailed information
#'   see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5330805/


NULL
