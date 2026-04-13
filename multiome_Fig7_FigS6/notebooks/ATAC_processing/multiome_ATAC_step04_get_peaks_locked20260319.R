library(parallel)
library(ArchR)

Sys.setenv(RHDF5_USE_FILE_LOCKING="FALSE")
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE")

addArchRGenome("mm10")
addArchRThreads(threads = 1)

args = commandArgs(trailingOnly=T)
proj_name <- args[1]

print(args)

proj <- loadArchRProject(proj_name)

#####################################
# Continue standard scATAC processing
#####################################

getAvailableMatrices(proj)
colnames(proj)

proj <- addTileMatrix(proj, force=TRUE)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", force=TRUE)
proj <- saveArchRProject(ArchRProj = proj)

# Gene scores with selected features
# Artificial black list to exclude all non variable features
chrs <- getChromSizes(proj)
var_features <- proj@reducedDims[["IterativeLSI"]]$LSIFeatures
var_features_gr <- GRanges(var_features$seqnames, IRanges(var_features$start, var_features$start + 500))
blacklist <- setdiff(chrs, var_features_gr)
proj <- addGeneScoreMatrix(proj, matrixName='GeneScoreMatrix', force=TRUE, blacklist=blacklist, subThreading=FALSE)
proj <- saveArchRProject(ArchRProj = proj)

groupBy <- "pheno_30"
proj <- addGroupCoverages(proj, maxFragmentLength=147, groupBy = groupBy, force=TRUE)
proj <- saveArchRProject(ArchRProj = proj)
proj <- addReproduciblePeakSet(proj, groupBy = groupBy)
proj <- saveArchRProject(ArchRProj = proj)

# Counts
proj <- addPeakMatrix(proj, maxFragmentLength=147, ceiling=10^9)
proj <- saveArchRProject(ArchRProj = proj)
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif_cisbp", force=TRUE)
# Save 
proj <- saveArchRProject(ArchRProj = proj)
