main_mpln<-function(i, y, membership, Gmin, Gmax, mod, n_chain, numb_iterations=NA, init_method=NA, init_iterations=NA, normalize=NA){
  ptm<-proc.time() 
  
  source("AIC_function.R")
  source("AIC3_function.R")
  source("BIC_function.R")
  source("calc_likelihood.R")
  source("ICL_function.R")
  source("calculate_parameters.R")
  source("calling_clustering.R")
  source("cluster_mpln.R")
  source("initializationrun.R")
  source("main_mpln.R")
  source("PackageCheck.R")
  source("stanrun.R")
  source("zvalue_calculation.R")
  
  # Loading packages needed for RStan model
  LoadCheckPkg(pckgs=c("mvtnorm","mclust","coda","capushe","edgeR",      
    "clusterGeneration","pheatmap","RColorBrewer","gplots","rstan","Rcpp","parallel"))

  
  # Calculate the number of cores
  no_cores = detectCores()-1
  
  # Initiate cluster
  cl = makeCluster(no_cores) 
  
  cat("\n Doing clusterExport")
  clusterExport(cl,c("mod", "zvalue_calculation", "calc_likelihood", "stanrun", "initializationrun", "BIC_function","ICL_function","AIC_function","AIC3_function", "calculate_parameters", "cluster_mpln", "calling_clustering"))
  
  
  cat("\n Doing clusterEvalQ")
  #other packages need to be downloaded using clusterEvalQ
  clusterEvalQ(cl, library(rstan))
  clusterEvalQ(cl, library(Rcpp))
  clusterEvalQ(cl, library(mclust))
  clusterEvalQ(cl, library(mvtnorm))
  clusterEvalQ(cl, library(edgeR))
  clusterEvalQ(cl, library(capushe))
  clusterEvalQ(cl, library(clusterGeneration))
  clusterEvalQ(cl, library(coda))
  clusterEvalQ(cl, library(parallel))
  clusterEvalQ(cl, library(parallel))
  
  if (typeof(y) != "double" & typeof(y) != "integer"){
    stop("Dataset type needs to be integer");}
  
  if (Gmax<Gmin){
    stop("Gmax cannot be less than Gmin");}
  
  if(is.na(numb_iterations)) numb_iterations <- 1000
  

  if(numb_iterations<40){
    stop("RStan numb_iterations argument should be greater than 40");}
  
  if((is.na(init_iterations) != TRUE && init_iterations == !0) && is.na(init_method) == TRUE){
    stop("Number of initialization iterations specified, but no initialization method selected");}
  
  d<-ncol(y) # Saving number of dimensions
  n<-nrow(y) # Saving number of observations
  
  if(all(is.na(membership)!=TRUE) && length(membership)!=n){
    stop("Length of membership character vector and sample size of dataset should match");}
  
  if(all(is.na(membership)!=TRUE) && all((diff(sort(unique(membership)))==1)!=TRUE) ){
    stop("Cluster memberships in the membership vector are missing a cluster, e.g. 1,3,4,5,6 is missing cluster 2");}
  
  if(length(which(apply(y, 1, function(x) all(x==0))==TRUE))!=0){
    cat("\nDataset row(s)", c(which(apply(y, 1, function(x) all(x==0))==TRUE)), "will be removed as this/these contain(s) all zeros")
    if(all(is.na(membership)==FALSE)){membership<-membership[-c(which(apply(y, 1, function(x) all(x==0))==TRUE))]}
    y<-y[-c(which(apply(y, 1, function(x) all(x==0))==TRUE)),]
    n<-nrow(y)
  }
  
  if(all(is.na(membership)==TRUE)){
    membership<-"Not provided"}
  
  if (Gmax > n){
    stop("Gmax cannot be larger than n");}
  
  if(is.na(normalize) == FALSE) {
    norm_factors<-log(as.vector(calcNormFactors(as.matrix(y), method = "TMM")))
  } else {norm_factors<-rep(0,d)}
  #cat("\nNormalize factors in main_mpln are: ",norm_factors)
  
  MPLN_parallel = function(g){
    ## ** Never use set.seed(), use clusterSetRNGStream() instead,
    # to set the cluster seed if you want reproducible results
    #clusterSetRNGStream(cl=cl, iseed=g)
    test = calling_clustering(y=y, Gmin=g, Gmax=g, n_chain=n_chain, numb_iterations=numb_iterations, init_method=init_method, init_iterations=init_iterations, norm_factors=norm_factors, mod=mod)
    return(test)
  }
  
  # empty list to save output
  parallel.Wei_2 = list()
  cat("\n Running parallel code now.")
  parallel.Wei_2 = clusterMap(cl=cl,fun=MPLN_parallel, g=Gmin:Gmax)
  cat("\n Done parallel code.")
  
  BIC<-ICL<-AIC<-AIC3<-Djump<-DDSE<-k<-ll<-vector()
  
  
  for(g in 1:(Gmax-Gmin+1)) {
    ll[g]<-unlist(tail(parallel.Wei_2[[g]]$allresults$loglikelihood, n=1)) # save the final log-likelihood
    
    k[g]<-calculate_parameters(g,y)
    
    if (g==max(1:(Gmax-Gmin+1))){ # starting model selection
      bic<-BIC_function(ll=ll,k=k, n=n, run=parallel.Wei_2, gmin=Gmin, gmax=Gmax)
      icl<-ICL_function(bIc=bic, gmin=Gmin, gmax=Gmax, run=parallel.Wei_2)
      aic<-AIC_function(ll=ll,k=k, run=parallel.Wei_2, gmin=Gmin, gmax=Gmax )
      aic3<-AIC3_function(ll=ll,k=k, run=parallel.Wei_2, gmin=Gmin, gmax=Gmax)
    }
  }
  cat("\nDone model selection for ", g)
  
  # for Djump and DDSE
  if((Gmax-Gmin+1) > 10 ) {
    if (!require(capushe)) install.packages('capushe') # loading needed package
    library(capushe)
    
    # adapted based on HTSCluster package 2.0.8 (25 Oct 2016)
    PMM<-allruns
    runs <- Gmin:Gmax
    Gmax <- Gmax
    logLike.final <- suppressWarnings(do.call("cbind", lapply(PMM, function(x) x$loglikelihood))) #gives log-likelihood for each cluster at each run
    logLike.val <- apply(logLike.final,1,max) 
    
    message("Note: diagnostic plots for results corresponding to model selection via slope heuristics (Djump and DDSE) should be examined to ensure that sufficiently complex models have been considered.")
    Kchoice <- Gmin:Gmax
    k <- k # number of parameters
    mat <- cbind(Kchoice, k/n, k/n, -logLike.val)
    ResCapushe <- capushe(mat, n)
    DDSEmodel<- ResCapushe@DDSE@model
    Djumpmodel<- ResCapushe@Djump@model
    final<-proc.time()-ptm
    
    RESULTS <- list(dataset= y,
      dimensionality = d,
      normalization_factors=norm_factors,
      gmin = Gmin,
      gmax = Gmax,
      initalization_method = init_method,
      allresults = parallel.Wei_2,
      loglikelihood = ll, 
      numbofparameters = k,
      truelabels = membership,
      ICL.all = icl,
      BIC.all = bic,
      AIC.all = aic,
      AIC3.all = aic3,
      SlopeHeuristics = ResCapushe,
      Djumpmodelselected = ResCapushe@Djump@model,
      DDSEmodelselected = ResCapushe@DDSE@model,
      totaltime = final)
    
  } else {# end of Djump and DDSE
    
    final<-proc.time()-ptm
    cat("\nSaving final time for ", g)
    
    RESULTS <- list(dataset= y,
      dimensionality = d,
      normalization_factors=norm_factors,
      gmin = Gmin,
      gmax = Gmax,
      initalization_method = init_method,
      allresults = parallel.Wei_2,
      loglikelihood = ll, 
      numbofparameters = k,
      truelabels = membership,
      ICL.all = icl,
      BIC.all = bic,
      AIC.all = aic,
      AIC3.all = aic3,
      SlopeHeuristics = "Not used",
      Djumpmodelselected = "Not used",
      DDSEmodelselected = "Not used",
      totaltime = final)
    cat("\nSaving results list for ", g)
  }
  
  class(RESULTS) <- "MPLN"
  return(RESULTS)
}
