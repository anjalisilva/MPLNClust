
mplnNotParallelInternal <- function(dataset, membership = "none",
  gmin = 1, gmax = 2,
  nChains = 3, nIterations = 1000,
  initMethod = "kmeans", nInitIterations = 0,
  normalize = "Yes") {

  ptm <- proc.time()

  # Performing checks
  if (typeof(dataset) != "double" & typeof(dataset) != "integer") {
    stop("Dataset type needs to be integer")
  }

  if (class(dataset) != "matrix") {
    stop("Dataset needs to be a matrix")
  }

  dataset <- removeZeroCounts(dataset = dataset)$dataset

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin")
  }

  if(is.na(nIterations)) {
    nIterations <- 1000
  }

  # if(is.na(nChains) || nChains<3) {
  #  nChains <- 3
  #  cat("Recommended number of chains is minimum 3. nChains' set to 3")}

  if(nIterations < 40) {
    stop("RStan nIterations argument should be greater than 40")}


  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if(class(membership) != "character" && length(membership) != nObservations) {
    stop("Length of membership character vector and
      sample size of dataset should match")
  }

  if(class(membership) != "character" &&
      all((diff(sort(unique(membership))) == 1) != TRUE) ) {
    stop("Cluster memberships in the membership vector
      are missing a cluster, e.g. 1, 3, 4, 5, 6 is missing cluster 2")
  }

  if(length(which(apply(dataset, 1, function(x) all(x == 0)) == TRUE)) != 0) {
    cat("\nDataset row(s)",
      c(which(apply(dataset, 1, function(x) all(x == 0)) == TRUE)),
      "will be removed as this/these contain(s) all zeros")

    if(class(membership) != "character") {
      membership <- membership[- c(which(apply(dataset, 1, function(x)
        all(x == 0)) == TRUE))]
    }

    dataset <- dataset[- c(which(apply(dataset, 1, function(x)
      all(x == 0)) == TRUE)), ]
    nObservations <- nrow(dataset)
  }

  if(class(membership) == "character") {
    membership <- "Not provided"
  }

  if (gmax > nObservations) {
    stop("gmax cannot be larger than nrow(dataset)")
  }

  # Calculating normalization factors
  if(normalize == "Yes") {
    normFactors <- log(as.vector(edgeR::calcNormFactors(as.matrix(dataset),
      method = "TMM")))
  } else if(normalize == "No") {
    normFactors <- rep(0, dimensionality)
  } else{
    stop("normFactors should be 'Yes' or 'No' ")
  }

  # Construct a Stan model
  # Stan model was developed with help of Sanjeena Dang
  # <sdang@math.binghamton.edu>
  stancode <- 'data{int<lower=1> d; // Dimension of theta
  int<lower=0> N; //Sample size
  int y[N,d]; //Array of Y
  vector[d] mu;
  matrix[d,d] Sigma;
  vector[d] normfactors;
  }
  parameters{
  matrix[N,d] theta;
  }
  model{ //Change values for the priors as appropriate
  for (n in 1:N){
  theta[n,]~multi_normal(mu,Sigma);
  }
  for (n in 1:N){
  for (k in 1:d){
  real z;
  z=exp(normfactors[k]+theta[n,k]);
  y[n,k]~poisson(z);
  }
  }
  }'

  mod <- rstan::stan_model(model_code = stancode, verbose = FALSE)
  allruns <- list()

  # Constructing non parallel code
  for (gmodel in seq_along(1:(gmax - gmin + 1))) {

    if(length(1:(gmax - gmin + 1)) == gmax) {
      clustersize <- gmodel
    } else if(length(1:(gmax - gmin + 1)) < gmax) {
      clustersize <- seq(gmin, gmax, 1)[gmodel]
    }

    if(nInitIterations != 0) {
      # cat("\nRunning initialization for G =", clustersize)
      initializeruns <- initializationRun(gmodel = clustersize,
        dataset = dataset,
        init_method = initMethod,
        init_iterations = nInitIterations,
        n_chain = nChains,
        numb_iterations = nIterations,
        initialization = NA,
        normalizefactors = normFactors,
        mod = mod)
      # cat("\nInitialization done for G =", clustersize)
      # cat("\nRunning clustering for G =", clustersize)
      allruns[[gmodel]] <- mplnCluster(dataset = dataset,
        z = NA,
        G = clustersize,
        nChains = nChains,
        nIterations = nIterations,
        initialization = initializeruns,
        normalizefac = normFactors,
        mod = mod)
      # cat("\nClustering done for G =", clustersize)
    } else if(nInitIterations == 0) {
      # cat("\nNo initialization done for G =", clustersize)
      # cat("\nRunning clustering for G =", clustersize)
      allruns[[gmodel]] <- mplnCluster(dataset = dataset,
        z = mclust::unmap(stats::kmeans(log(dataset + 1/3),
          clustersize)$cluster),
        G = clustersize,
        nChains = nChains,
        nIterations = nIterations,
        initialization = NA,
        normalizefac = normFactors,
        mod = mod)
      # cat("\nClustering done for G =", clustersize)
    }
  }

  RESULTS <- list(gmin = gmin,
    gmax = gmax,
    initalization_method = initMethod,
    all_results = allruns)

  class(RESULTS) <- "mplnNotParallelInternal"
  return(RESULTS)
}

mplnNotParallel <- function(){

  mplnNotParallelInternal(dataset = dataset, membership = membership,
    gmin = 1, gmax = 1, nChains = 3, nIterations = 1000,
    initMethod = "kmeans", nInitIterations = 1,
    normalize = "Yes")

  BIC <- ICL <- AIC <- AIC3 <- Djump <- DDSE <- k <- ll <- vector()

  for(g in seq_along(1:(gmax - gmin + 1))) {
    # save the final log-likelihood
    ll[g] <- unlist(utils::tail(allruns[[g]]$loglikelihood, n = 1))

    k[g] <- calcParameters(numberG = g, dimensionality = dimensionality)

    if (g == max(1:(gmax - gmin + 1))) { # starting model selection
      bic <- BICFunction(ll = ll,
        k = k,
        n = nObservations,
        run = allruns,
        gmin = gmin,
        gmax = gmax)

      icl <- ICLFunction(bIc = bic,
        gmin = gmin,
        gmax = gmax,
        run = parallelRun)

      aic <- AICFunction(ll = ll,
        k = k,
        run = parallelRun,
        gmin = gmin,
        gmax = gmax )

      aic3 <- AIC3Function(ll = ll,
        k = k,
        run = parallelRun,
        gmin = gmin,
        gmax = gmax)
    }
  }

  # for Djump and DDSE
  if((gmax - gmin + 1) > 10 ) {
    # adapted based on HTSCluster package 2.0.8 (25 Oct 2016)
    PMM <- parallelRun
    runs <- gmin:gmax
    gmax <- gmax
    logLike.final <- suppressWarnings(do.call("cbind",
      lapply(PMM, function(x) x$loglikelihood)))
    # gives log-likelihood for each cluster at each run
    logLike.val <- apply(logLike.final, 1, max)


    message("Note: diagnostic plots for results corresponding
      to model selection via slope heuristics (Djump and DDSE)
      should be examined to ensure that sufficiently complex
      models have been considered.")
    Kchoice <- gmin:gmax
    k <- k # number of parameters
    mat <- cbind(Kchoice, k/nObservations, k/nObservations,
      - logLike.val)
    ResCapushe <- capushe::capushe(mat, nObservations)
    DDSEmodel <- ResCapushe@DDSE@model
    Djumpmodel <- ResCapushe@Djump@model
    final <- proc.time() - ptm

    RESULTS <- list(dataset = dataset,
      dimensionality = dimensionality,
      normalization_factors = normFactors,
      gmin = gmin,
      gmax = gmax,
      initalization_method = initMethod,
      all_results = parallelRun,
      loglikelihood = ll,
      numb_of_parameters = k,
      true_labels = membership,
      ICL_all = icl,
      BIC_all = bic,
      AIC_all = aic,
      AIC3_all = aic3,
      slope_heuristics = ResCapushe,
      Djumpmodel_selected = ResCapushe@Djump@model,
      DDSEmodel_selected = ResCapushe@DDSE@model,
      total_time = final)

  } else {# end of Djump and DDSE
    final <- proc.time() - ptm
    RESULTS <- list(dataset = dataset,
      dimensionality = dimensionality,
      normalization_factors=normFactors,
      gmin = gmin,
      gmax = gmax,
      initalization_method = initMethod,
      all_results = parallelRun,
      loglikelihood = ll,
      numb_of_parameters = k,
      true_labels = membership,
      ICL_all = icl,
      BIC_all = bic,
      AIC_all = aic,
      AIC3_all = aic3,
      slope_heuristics = "Not used",
      total_time = final)
  }

  class(RESULTS) <- "MPLN"
  return(RESULTS)
  # [END]

}
