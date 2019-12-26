# Initialization 
initializationRun <- function(gmodel, y, init_method, 
                              init_iterations, n_chain, 
                              numb_iterations, 
                              initialization = NA, 
                              normalizefactors, mod) {
  
  z <- init_runs <- list()
  logL_init <- vector()
  n <- nrow(y)
  d <- ncol(y)
  
  for(iterations in 1:init_iterations) {
    if (init_method == "kmeans" | is.na(init_method)) {
      z[[iterations]] <- unmap(kmeans(log(y + 1/3), gmodel)$cluster)
    }else if (init_method == "random") {
      if(gmodel == 1) { # generating z if g=1
        z[[iterations]] <- as.matrix(rep.int(1, times = n), 
                                     ncol = gmodel,
                                     nrow = n)
      } else { # generating z if g>1
        z_conv = 0
        while(! z_conv) { 
          # ensure that dimension of z is same as G (i.e., 
          # if one column contains all 0s, then generate z again)
          z[[iterations]] <- t(rmultinom(n, size = 1,
                               prob = rep(1 / gmodel, gmodel))) 
          if(length(which(colSums(z[[iterations]]) > 0)) == gmodel) {
            z_conv = 1
          }
        }
      }
    }else if (init_method == "medoids") {
      LoadCheckPkg(pckgs = c("cluster"))
      z[[iterations]] <- unmap(pam(log(y + 1/3), k = gmodel)$cluster)
    }else if (init_method == "clara") {
      LoadCheckPkg(pckgs = c("cluster"))
      z[[iterations]] <- unmap(clara(log(y + 1/3), k = gmodel)$cluster)
    }else if (init_method == "fanny") {
      LoadCheckPkg(pckgs = c("cluster"))
      z[[iterations]] <- unmap(fanny(log(y + 1/3), k = gmodel)$cluster)
    }
    
    init_runs[[iterations]] <- cluster_mpln(y = y, 
                                            z = z[[iterations]],
                                            G = gmodel, 
                                            n_chain = n_chain,
                                            numb_iterations = numb_iterations,
                                            initialization = "init", 
                                            normalizefac = normalizefactors, 
                                            mod = mod)
    logL_init[iterations] <- 
      unlist(tail((init_runs[[iterations]]$loglikelihood), n = 1)) 
  }
  
  initialization <- init_runs[[which(logL_init == max(logL_init, 
                                                      na.rm = TRUE))[1]]]
  return(initialization)
  # Developed by Anjali Silva
}  
