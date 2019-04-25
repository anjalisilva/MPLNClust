
main_mpln <- function(dataset, membership, Gmin, Gmax, n_chains, n_iterations=NA, init_method=NA, n_init_iterations=NA, normalize=NA){
  
  ptm<-proc.time() 
  
  if (typeof(dataset) != "double" & typeof(dataset) != "integer"){
    stop("Dataset type needs to be integer");}
  
  dataset<-remove_zero_counts(dataset=dataset)$dataset
  
  if (Gmax<Gmin){
    stop("Gmax cannot be less than Gmin");}
  
  if(is.na(n_iterations)) n_iterations <- 1000
  
  #if(is.na(n_chains) || n_chains<3) {
  #  n_chains <- 3
  #  cat("Recommended number of chains is minimum 3. n_chains' set to 3")}
  
  if(n_iterations<40){
    stop("RStan n_iterations argument should be greater than 40");}
  
  if((is.na(n_init_iterations) != TRUE && n_init_iterations == !0) && is.na(init_method) == TRUE){
    stop("Number of initialization iterations specified, but no initialization method selected");}
  
  d<-ncol(dataset)
  n<-nrow(dataset)
  
  if(all(is.na(membership)!=TRUE) && length(membership)!=n){
    stop("Length of membership character vector and sample size of dataset should match");}
  
  if(all(is.na(membership)!=TRUE) && all((diff(sort(unique(membership)))==1)!=TRUE) ){
    stop("Cluster memberships in the membership vector are missing a cluster, e.g. 1,3,4,5,6 is missing cluster 2");}
  
  if(length(which(apply(dataset, 1, function(x) all(x==0))==TRUE))!=0){
    cat("\nDataset row(s)", c(which(apply(dataset, 1, function(x) all(x==0))==TRUE)), "will be removed as this/these contain(s) all zeros")
    if(all(is.na(membership)==FALSE)){membership<-membership[-c(which(apply(dataset, 1, function(x) all(x==0))==TRUE))]}
    dataset<-dataset[-c(which(apply(dataset, 1, function(x) all(x==0))==TRUE)),]
    n<-nrow(dataset)
  }
  
  if(all(is.na(membership)==TRUE)){
    membership<-"Not provided"}
  
  if (Gmax > n){
    stop("Gmax cannot be larger than n");}
  
  if(is.na(normalize) == FALSE) {
    norm_factors<-log(as.vector(calcNormFactors(as.matrix(dataset), method = "TMM")))
  } else {norm_factors<-rep(0,d)}
  #cat("\nNormalize factors in main_mpln are: ",norm_factors)
  
  MPLN_parallel = function(g){
    ## ** Never use set.seed(), use clusterSetRNGStream() instead,
    # to set the cluster seed if you want reproducible results
    #clusterSetRNGStream(cl=cl, iseed=g)
    test = calling_clustering(y=dataset, Gmin=g, Gmax=g, n_chain=n_chains, numb_iterations=n_iterations, init_method=init_method, init_iterations=n_init_iterations, norm_factors=norm_factors, mod=mod)
    return(test)
  }
  #cat("Done MPLN_parallel")
  
  # empty list to save output
  parallel.Wei_2 = list()
  cat("\nRunning parallel code now.")
  parallel.Wei_2 = clusterMap(cl=cl,fun=MPLN_parallel, g=Gmin:Gmax)
  cat("\nDone parallel code.")
  
  BIC<-ICL<-AIC<-AIC3<-Djump<-DDSE<-k<-ll<-vector()
  

  for(g in 1:(Gmax-Gmin+1)) {
    ll[g]<-unlist(tail(parallel.Wei_2[[g]]$allresults$loglikelihood, n=1)) # save the final log-likelihood
    
    k[g]<-calculate_parameters(g,y=dataset)
    
    if (g==max(1:(Gmax-Gmin+1))){ # starting model selection
      bic<-BIC_function(ll=ll,k=k, n=n, run=parallel.Wei_2, gmin=Gmin, gmax=Gmax)
      icl<-ICL_function(bIc=bic, gmin=Gmin, gmax=Gmax, run=parallel.Wei_2)
      aic<-AIC_function(ll=ll,k=k, run=parallel.Wei_2, gmin=Gmin, gmax=Gmax )
      aic3<-AIC3_function(ll=ll,k=k, run=parallel.Wei_2, gmin=Gmin, gmax=Gmax)
    }
  }
  
  # for Djump and DDSE
  if((Gmax-Gmin+1) > 10 ) {
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
    
    RESULTS <- list(dataset= dataset,
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
    RESULTS <- list(dataset= dataset,
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
  }
  
  class(RESULTS) <- "MPLN"
  return(RESULTS)
}
