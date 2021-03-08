# hierarchicell
An R package for simulating cell-type specific and hierarchical single-cell 
expression data.The hierarchicell package estimates important parameters from
single-cell RNAseq expression counts before simulating expression data that
is cell-type specific and hierarchical in nature. With the simulated data,
hierarchicell is then able to compute power calculations under a variety of
conditions. The package consists of four important categories of functions:
data loading and cleaning, empirical estimation of distributions,
simulating expression data, and computing type 1 error or power.

The data loading and cleaning function is very basic, but data input is
critical to the package working correctly. If no input data is given, the
package default data will be used for simulation and power calculations.

The most fundamental component of this package is in the estimation of the
simulation parameters. The functions to estimate parameters for the
simulation estimate the empirical distributions for library size, dropout
rate, and global gene means and model the hierarchical variance structure
of the input data. 

With the parameters estimated, the package can simulate data under a
variety of pre-determined conditions. These conditions include foldchange,
number of genes, number of samples (i.e., independent experimental units),
and the mean number of cells per individual.Either binary (case/control)
or continuous outcomes of interest can be simulated. 

With the parameters estimated, the package can compute type 1 error rates
for a number of different tools under a variety of pre-determined
conditions. These conditions include number of genes, number of samples
(i.e., independent experimental units), and the mean number of cells per
individual. 

With the parameters estimated, the package can compute power under a
variety of pre-determined conditions. These conditions include foldchange,
number of genes, number of samples (i.e., independent experimental units),
and the mean number of cells per individual. Power can be calculated for 
either binary or continuous outcomes of interest. 

# Installation

To install, type the following command:
```{r, echo = FALSE}

devtools::install_github("kdzimm/hierarchicell")

```

# Example

```{r}

### Load and filter a cell-type specific scRNA-seq datset of interest
# (Here, we are just loading the default dataset)

clean_expr_data <- filter_counts()

### Estimate parameters from the simulated data

data_summaries <- compute_data_summaries(clean_expr_data)

### Simulate hierarchical scRNA-seq data similar to input data

simulated_counts <- simulate_hierarchicell(data_summaries)

### Estimate power with a two-part hurdle mixed model

power_hierarchicell(data_summaries, n_genes = 100, n_per_group = 10)


```
