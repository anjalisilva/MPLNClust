
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MPLNClust

Finite Mixtures of Multivariate Poisson-Log Normal Model for Clustering
Count Data

<!-- badges: start -->
<!-- https://www.codefactor.io/repository/github/anjalisilva/MPLNClust/issues -->
<!-- [![CodeFactor](https://www.codefactor.io/repository/github/anjalisilva/mplnclust/badge)](https://www.codefactor.io/repository/github/anjalisilva/mplnclust)-->

[![GitHub
issues](https://img.shields.io/github/issues/anjalisilva/MPLNClust)](https://github.com/anjalisilva/MPLNClust/issues)
[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)
![GitHub language
count](https://img.shields.io/github/languages/count/anjalisilva/MPLNClust)
![GitHub commit activity
(branch)](https://img.shields.io/github/commit-activity/y/anjalisilva/MPLNClust/master)

<!-- https://shields.io/category/license -->
<!-- badges: end -->

## Description

`MPLNClust` is an R package for performing clustering using finite
mixtures of multivariate Poisson-log normal (MPLN) distribution proposed
by [Silva et al., 2019](https://pubmed.ncbi.nlm.nih.gov/31311497/). It
was developed for count data, with clustering of RNA sequencing data as
a motivation. However, the clustering method may be applied to other
types of count data. The package provides functions for functions for
parameter estimation via 1) an MCMC-EM framework by [Silva et al.,
2019](https://pubmed.ncbi.nlm.nih.gov/31311497/) and 2) a variational
Gaussian approximation with EM algorithm by [Subedi and Browne,
2020](https://doi.org/10.1002/sta4.310). Information criteria (AIC, BIC,
AIC3 and ICL) and slope heuristics (Djump and DDSE, if more than 10
models are considered) are offered for model selection. Also included
are functions for simulating data from this model and visualization.

## Installation

To install the latest version of the package:

``` r
require("devtools")
devtools::install_github("anjalisilva/MPLNClust", build_vignettes = TRUE)
library("MPLNClust")
```

To run the Shiny app:

``` r
MPLNClust::runMPLNClust()
```

## Overview

To list all functions available in the package:

``` r
ls("package:MPLNClust")
```

`MPLNClust` contains 14 functions.

1.  ***mplnVariational*** for carrying out clustering of count data
    using mixtures of MPLN via variational expectation-maximization
2.  ***mplnMCMCParallel*** for carrying out clustering of count data
    using mixtures of MPLN via a Markov chain Monte Carlo
    expectation-maximization algorithm (MCMC-EM) with parallelization
3.  ***mplnMCMCNonParallel*** for carrying out clustering of count data
    using mixtures of MPLN via a Markov chain Monte Carlo
    expectation-maximization algorithm (MCMC-EM) with no parallelization
4.  ***mplnDataGenerator*** for the purpose of generating simlulation
    data via mixtures of MPLN
5.  ***mplnVisualizeAlluvial*** for visualizing clustering results as
    Alluvial plots
6.  ***mplnVisualizeBar*** for visualizing clustering results as bar
    plots
7.  ***mplnVisualizeHeatmap*** for visualizing clustering results as
    heatmaps
8.  ***mplnVisualizeLine*** for visualizing clustering results as line
    plots
9.  ***AICFunction*** for model selection
10. ***AIC3Function*** for model selection
11. ***BICFunction*** for model selection
12. ***ICLFunction*** for model selection
13. ***runMPLNClust*** is the shiny implementation of *mplnVariational*
14. ***mplnVarClassification*** is an implementation for classification
    is currently under construction

Framework of ***mplnVariational*** makes it computationally efficient
and faster compared to ***mplnMCMCParallel*** or
***mplnMCMCNonParallel***. Therefore, ***mplnVariational*** may perform
better for large datasets. For more information, see details section
below. An overview of the package is illustrated below:

<div style="text-align:center">

<img src="inst/extdata/Overview_MPLNClust.png" width="800" height="450"/>

<div style="text-align:left">
<div style="text-align:left">

## Details

The MPLN distribution ([Aitchison and Ho,
1989](https://www.jstor.org/stable/2336624?seq=1)) is a multivariate log
normal mixture of independent Poisson distributions. The hidden layer of
the MPLN distribution is a multivariate Gaussian distribution, which
allows for the specification of a covariance structure. Further, the
MPLN distribution can account for overdispersion in count data.
Additionally, the MPLN distribution supports negative and positive
correlations.

A mixture of MPLN distributions is introduced for clustering count data
by [Silva et al., 2019](https://pubmed.ncbi.nlm.nih.gov/31311497/).
Here, applicability is illustrated using RNA sequencing data. To this
date, two frameworks have been proposed for parameter estimation: 1) an
MCMC-EM framework by [Silva et al.,
2019](https://pubmed.ncbi.nlm.nih.gov/31311497/) and 2) a variational
Gaussian approximation with EM algorithm by [Subedi and Browne,
2020](https://doi.org/10.1002/sta4.310).

### MCMC-EM Framework for Parameter Estimation

[Silva et al., 2019](https://pubmed.ncbi.nlm.nih.gov/31311497/) used an
MCMC-EM framework via Stan for parameter estimation. This method is
employed in functions ***mplnMCMCParallel*** and
***mplnMCMCNonParallel***.

Coarse grain parallelization is employed in *mplnMCMCParallel*, such
that when a range of components/clusters (g = 1,…,G) are considered,
each component/cluster size is run on a different processor. This can be
performed because each component/cluster size is independent from
another. All components/clusters in the range to be tested have been
parallelized to run on a separate core using the *parallel* R package.
The number of cores used for clustering is calculated using
*parallel::detectCores() - 1*. No internal parallelization is performed
for *mplnMCMCNonParallel*.

To check the convergence of MCMC chains, the potential scale reduction
factor and the effective number of samples are used. The Heidelberger
and Welch’s convergence diagnostic (Heidelberger and Welch, 1983) is
used to check the convergence of the MCMC-EM algorithm. Starting values
(argument: *initMethod*) and the number of iterations for each chain
(argument: *nIterations*) play an important role for the successful
operation of this algorithm.

### Variational-EM Framework for Parameter Estimation

[Subedi and Browne, 2020](https://doi.org/10.1002/sta4.310) proposed a
variational Gaussian approximation that alleviates challenges of MCMC-EM
algorithm. Here the posterior distribution is approximated by minimizing
the Kullback-Leibler (KL) divergence between the true and the
approximating densities. A variational-EM based framework is used for
parameter estimation. This algorithm is implemented in the function
***mplnVariational***. The parsimonious family of models implemented by
considering eigen-decomposition of covariance matrix in [Subedi and
Browne, 2020](https://doi.org/10.1002/sta4.310) is not yet available
with this package.

## Model Selection and Other Details

Four model selection criteria are offered, which include the Akaike
information criterion (AIC; Akaike, 1973), the Bayesian information
criterion (BIC; Schwarz, 1978), a variation of the AIC used by Bozdogan
(1994) called AIC3, and the integrated completed likelihood (ICL;
Biernacki et al., 2000). Slope heuristics (Djump and DDSE; Arlot et al.,
2016) could be used for model selection if more than 10 models are
considered.

Starting values (argument: *initMethod*) and the number of iterations
for each chain (argument: *nInitIterations*) play an important role to
the successful operation of this algorithm. There maybe issues with
singularity, in which case altering starting values or initialization
method may help.

## Shiny App

The Shiny app employing ***mplnVariational*** could be run and results
could be visualized:

``` r
MPLNClust::runMPLNClust()
```

<div style="text-align:center">

<img src="inst/extdata/ShinyAppMPLNClust1.png" alt="ShinyApp1" width="650" height="400"/>

<div style="text-align:left">
<div style="text-align:left">
&#10;
## Tutorials

For tutorials and plot interpretation, refer to the vignette:

``` r
browseVignettes("MPLNClust")
```

## Citation for Package

``` r
citation("MPLNClust")
```

Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2019). A
multivariate Poisson-log normal mixture model for clustering
transcriptome sequencing data. *BMC Bioinformatics*. 20(1):394.

``` r
A BibTeX entry for LaTeX users is

  @Article{,
    title = {A multivariate Poisson-log normal mixture model for clustering transcriptome sequencing data},
    author = {A. Silva and S. J. Rothstein and P. D. McNicholas and S. Subedi},
    journal = {BMC Bioinformatics},
    year = {2019},
    volume = {20},
    number = {1},
    pages = {394},
    url = {https://pubmed.ncbi.nlm.nih.gov/31311497/},
  }
```

## Package References

- [Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2019). A
  multivariate Poisson-log normal mixture model for clustering
  transcriptome sequencing data. *BMC
  Bioinformatics.*](https://pubmed.ncbi.nlm.nih.gov/31311497/)

- [Subedi, S., R.P. Browne (2020). A family of parsimonious mixtures of
  multivariate Poisson-lognormal distributions for clustering
  multivariate count data. *Stat.*
  9:e310.](https://doi.org/10.1002/sta4.310)

## Other References

- [Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log
  normal distribution.
  *Biometrika.*](https://www.jstor.org/stable/2336624?seq=1)

- [Akaike, H. (1973). Information theory and an extension of the maximum
  likelihood principle. In *Second International Symposium on
  Information Theory*, New York, NY, USA, pp. 267–281. Springer
  Verlag.](https://link.springer.com/chapter/10.1007/978-1-4612-1694-0_15)

- [Arlot, S., Brault, V., Baudry, J., Maugis, C., Michel, B. (2016).
  *capushe: CAlibrating Penalities Using Slope HEuristics*. R package
  version 1.1.1](https://CRAN.R-project.org/package=capushe)

- [Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture
  model for clustering with the integrated classification likelihood.
  *IEEE Transactions on Pattern Analysis and Machine Intelligence*
  22.](https://hal.inria.fr/inria-00073163/document)

- [Bozdogan, H. (1994). Mixture-model cluster analysis using model
  selection criteria and a new informational measure of complexity. In
  *Proceedings of the First US/Japan Conference on the Frontiers of
  Statistical Modeling: An Informational Approach: Volume 2 Multivariate
  Statistical Modeling*, pp.69–113. Dordrecht: Springer
  Netherlands.](https://link.springer.com/chapter/10.1007/978-94-011-0800-3_3)

- [Ghahramani, Z. and Beal, M. (1999). Variational inference for
  bayesian mixtures of factor analysers. *Advances in neural information
  processing systems*
  12.](https://cse.buffalo.edu/faculty/mbeal/papers/nips99.pdf)

- [Robinson, M.D., and Oshlack, A. (2010). A scaling normalization
  method for differential expression analysis of RNA-seq data. *Genome
  Biology* 11,
  R25.](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25)

- [Schwarz, G. (1978). Estimating the dimension of a model. *The Annals
  of Statistics* 6.](https://www.jstor.org/stable/2958889?seq=1)

- [Wainwright, M. J. and Jordan, M. I. (2008). Graphical models,
  exponential families, and variational inference. *Foundations and
  Trends® in Machine Learning*
  1.](https://onlinelibrary.wiley.com/doi/abs/10.1002/sta4.310)

## Maintainer

- Anjali Silva (<anjali@alumni.uoguelph.ca>).

## Contributions

`MPLNClust` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/anjalisilva/MPLNClust/issues).

## Acknowledgments

- Dr. Marcelo Ponce, SciNet HPC Consortium, University of Toronto, ON,
  Canada for all the computational support.

- This work was funded by Natural Sciences and Engineering Research
  Council of Canada, Queen Elizabeth II Graduate Scholarship, and Arthur
  Richmond Memorial Scholarship.
