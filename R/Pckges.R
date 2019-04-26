# A function to check whether a package is available in the system, if not loads the package

# Code developed based on code provided by Dr. Marcelo Ponce, April 2019

pckges <- c("rstan",
            "Rcpp",
            "mvtnorm",
            "mclust",
            #"BiocManager",
            "clusterGeneration",
            "mclust",
            "cluster",
            "capushe",
            "coda")

otherPckgs <- c("edgeR")

NeededPackages <- function(pckges, otherPckgs="") {

        availablePckges <- .packages(all.available = TRUE)

        needTOinstall <- !(pckges %in% availablePckges)
        cat("Requested packages:")
        print(pckges)
        cat("installing...", pckges[needTOinstall], '\n')
        for (pck in pckges[needTOinstall]) {
                install.packages(pck)
        }

        needTOinstall <- !(otherPckgs %in% availablePckges)
        cat("Requested packages:")
        print(otherPckgs)
        cat("installing...", otherPckgs[needTOinstall], '\n')
        for (pck in otherPckgs[needTOinstall]) {
          # BiocManager::install(pck)
          source("https://bioconductor.org/biocLite.R")
          biocLite(pck)
        }
}

checkVersion <- function(pckges) {

        print(sessionInfo())

        for (pck in pckges){
                cat(pck, as.character(packageVersion(pck)), '\n')
        }

}


#NeededPackages(pckges,otherPckgs)

checkVersion(c(pckges,otherPckgs))

