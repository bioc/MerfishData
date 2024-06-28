context("Colon IBD")

checkSPE <- function(spe) {
    # check overall object
    expect_is(spe, "SpatialExperiment")
    expect_true(nrow(spe) > 940)
    expect_true(ncol(spe) > 1300000)
    
    # check colData
    cdat.cols <- c("cell_id", "sample_id", "sample_type", "mouse_id", 
                   "technical_repeat_number", "slice_id", "fov", "n_genes",
                   "tier1", "tier2", "tier3", "leiden_neigh")
    expect_true(all(cdat.cols %in% colnames(colData(spe))))
    
    # check spatialCoords
    expect_equal(ncol(spatialCoords(spe)), 2)
    expect_true(all(colnames(spatialCoords(spe)) == c("x", "y")))
    expect_true(is.numeric(spatialCoords(spe)))
    cdat <- colData(spe)
    expect_true(ncol(cdat) == length(cdat.cols))
    expect_true(nrow(cdat) == ncol(spe))
    
    # check altExps
    expect_true("Blank" %in% altExpNames(spe))
    blank_spe <- altExp(spe, "Blank")
    expect_is(blank_spe, "SpatialExperiment")
    expect_equal(nrow(blank_spe), 47)
    
    # check reducedDims
    expect_true("PCA" %in% reducedDimNames(spe))
    expect_true("UMAP" %in% reducedDimNames(spe))
    expect_true("UMAP_Tier1" %in% reducedDimNames(spe))
    expect_true("UMAP_Tier2" %in% reducedDimNames(spe))
    expect_true("UMAP_Tier3" %in% reducedDimNames(spe))
    expect_true(ncol(reducedDim(spe, "PCA")) == 50)
    
    # check metadata
    metadata <- metadata(spe)
    expect_true(length(metadata) > 0)
    expect_true(is.list(metadata))
}

test_that("standard usage", {
  spe <- MouseColonIbdCadinu2024()
  checkSPE(spe)
})


