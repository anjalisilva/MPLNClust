# Script and Functions for seting up required R-packages for MPLNClust 


# define dependencies, ie. needed to be installed for the code to run
# first packages that can be installaed from CRAN using install.packages()
pckges <- c("rstan",
                "Rcpp",
                "mvtnorm",
                "mclust",
                #"BiocManager",
                "clusterGeneration",
                'mclust',
                'cluster',
                'capushe',
                'coda')
# other packages that have to be installed in a different manner, eg. packages from BioConductor
otherPckgs <- c("edgeR")


##############################################################################################
# functions for checking whether dependencies are installed and install all needed packages


# function to install neeeded packages
NeededPackages <- function(pckges, otherPckgs="", def.mirror='https://cloud.r-project.org') {
	RverM <- as.numeric(R.Version()['major'])
	Rverm <- as.numeric(R.Version()['minor'])

        availablePckges <- .packages(all.available = TRUE)

	# deal with packages from CRAN
        needTOinstall <- !(pckges %in% availablePckges)
        cat("Requested packages:")
        print(pckges)
        cat("installing...", pckges[needTOinstall], '\n')
        for (pck in pckges[needTOinstall]) {
                install.packages(pck,repos=def.mirror)
        }

	# deal with packages from BioConductor
        needTOinstall <- !(otherPckgs %in% availablePckges)
        cat("Requested packages:")
        print(otherPckgs)
        cat("installing...", otherPckgs[needTOinstall], '\n')
        for (pck in otherPckgs[needTOinstall]) {
		print(RverM); print(Rverm)
		if ((RverM >= 3 ) && (Rverm > 5)) {
			# newer R version
                        install.packages("BiocManager", repos=def.mirror)
                        BiocManager::install(pck)
		} else {
			# older R versions...
			source("https://bioconductor.org/biocLite.R") 
			library(BiocInstaller)
			BiocInstaller::biocLite(pck)
		}
        }
}

# function to check the versions of a given set of packages...
checkVersion <- function(pckges) {

        print(sessionInfo())

	print(pckges)

        for (pck in pckges){
                cat(pck, as.character(packageVersion(pck)), '\n')
        }

}

##############################################################################################

# check and install required packages
NeededPackages(pckges,otherPckgs)

# display versions of the installed/required packages
checkVersion(c(pckges,otherPckgs))

