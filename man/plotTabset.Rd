% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{plotTabset}
\alias{plotTabset}
\title{Plot a tabset}
\usage{
plotTabset(spe.list, img)
}
\arguments{
\item{spe.list}{A named \code{list} of \code{\linkS4class{SpatialExperiment}}
objects.}

\item{img}{a raster object representing a bitmap image.}
}
\value{
None. Produces a tabset for rendering with \code{rmarkdown}.
}
\description{
Plot a tabset of colData annotations of one or more
SpatialExperiment objects over an image.
}
\examples{
    # create simulated data as described in the SpatialExperiment man page
    example("SpatialExperiment", package = "SpatialExperiment", echo = FALSE)
    spe <- spe_mol   

    # add simulated cell centroids
    s <- cbind(x = runif(20), y = runif(20))    
    spatialCoords(spe) <- s 

    # add simulated cell type and cell cycle annotation
    ct <- c("ct1", "ct2", "ct3")
    cc <- c("G1", "G2", "S", "M") 
    spe$type <- sample(ct, ncol(spe), replace = TRUE)
    spe$cycle <- sample(cc, ncol(spe), replace = TRUE)

    # create an example image
    hgrid <- hcl(0, 80, seq(50, 80, 10))
    img <- as.raster(matrix(hgrid, nrow = 4, ncol = 5))

    # plotTabset
    spe.list <- list(myseg = spe)
    plotTabset(spe.list, img)  
}
