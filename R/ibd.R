#' MERFISH mouse colon IBD dataset from Cadinu et al., 2024
#' @description Obtain the MERFISH mouse colon IBD dataset from Cadinu et al.,
#' 2024
#' @details Gut inflammation involves contributions from immune and non-immune cells,
#' whose interactions are shaped by the spatial organization of the healthy gut and its
#' remodeling during inflammation.
#' The crosstalk between fibroblasts and immune cells is an important axis in this process,
#' but our understanding has been challenged by incomplete cell-type definition and biogeography.
#'
#' To address this challenge, Cadinu et al., 2024 used multiplexed error-robust
#' fluorescence in situ hybridization
#' (MERFISH) to profile the expression of 943 genes in 1.35 million cells imaged
#' across the onset and recovery from a mouse colitis model.
#' They identified diverse cell populations, charted their
#' spatial organization, and revealed their polarization or recruitment in
#' inflammation.
#'
#' The barcoding scheme contained 990 possible barcodes; 943 of them were used to
#' code the RNAs of the genes assayed via combinatorial smFISH across different
#' stages of colitis in a mouse model; 47 of these barcodes
#' were left unassigned ("blank"), providing a direct measure of the false-positive
#' rate in MERFISH. Measurements for these 47 blank barcodes is stored in an
#' `altExp` named \code{"blank"}.
#'
#' The provided data contain cell type labels with three different degrees of
#' detail. The data have been collected before the insurgence of colitis
#' (\code{sample_type="Healthy"}) and after some time intervals
#' (3 days, 9 days, 21 days).
#'
#' @return An object of class \code{\linkS4class{SpatialExperiment}}.
#' @references Cadinu et al. (2024) Charting the cellular biogeography in
#' colitis reveals fibroblast trajectories and coordinated spatial remodeling.
#' Cell, 187(8).
#' @source \url{https://doi.org/10.5061/dryad.rjdfn2zh3}
#' @examples spe <- MouseColonIBD2024()
#' @export

MouseColonIBD2024 <- function() {
    eh <- ExperimentHub::ExperimentHub()
    recs <- AnnotationHub::query(eh, c("MERFISH", "Cadinu2024"))

    # setup progress bar
    counter <- 0
    pb <- NULL
    if (interactive()) {
        nr.items <- length(recs) - 3
        pb <- txtProgressBar(counter, nr.items, style = 3)
    }

    # (1) retrieve colData and identify reducedDims
    cdat <- .getResource(recs, "_coldata")
    if (interactive()) {
        counter <- counter + 1
        setTxtProgressBar(pb, counter)
    }
    pcs <- grep("^PC\\d{1,3}$", colnames(cdat))
    umap <- grep("^UMAP", colnames(cdat))
    umap1 <- grep("Neigh_", colnames(cdat))
    umap2 <- grep("Tier2", colnames(cdat))
    umap3 <- grep("Tier3", colnames(cdat))

    # (2) create spe object for the blanks features
    blank.spe <- SpatialExperiment(
        assays = list(
            counts = .getHDF5(recs, "ibd_blanks_counts", "counts"),
            logcounts = .getHDF5(recs, "ibd_blanks_logcounts", "logcounts")
        )
    )
    if (interactive()) {
        counter <- counter + 1
        setTxtProgressBar(pb, counter)
    }
    rownames(blank.spe) <- paste("Blank", 1:47, sep = "")
    if (interactive()) {
        counter <- counter + 1
        setTxtProgressBar(pb, counter)
    }

    # (3) create spe with all the components
    spe <- SpatialExperiment(
        assays = list(
            counts = .getHDF5(recs, "ibd_counts", "counts"),
            logcounts = .getHDF5(recs, "ibd_logcounts", "logcounts")
        ),
        colData = cdat[, -c(pcs, umap, umap1, umap2, umap3)],
        metadata = .getResource(recs, "metadata"),
        reducedDims = list(
            PCA = cdat[, pcs],
            UMAP = cdat[, umap],
            UMAP_Tier1 = cdat[, umap1],
            UMAP_Tier2 = cdat[, umap2],
            UMAP_Tier3 = cdat[, umap3]
        ),
        altExps = list(Blank = blank.spe),
        spatialCoordsNames = c("x", "y")
    )
    if (interactive()) {
        counter <- counter + 1
        setTxtProgressBar(pb, counter)
    }

    rownames(spe) <- .getResource(recs, "ibd_rownames")
    if (interactive()) {
        counter <- counter + 1
        setTxtProgressBar(pb, counter)
    }

    return(spe)
}

.getHDF5 <- function(recs, suffix, name) {
    file_path <- .getResource(recs, suffix)
    obj <- HDF5Array(file_path, name)
    return(obj)
}
