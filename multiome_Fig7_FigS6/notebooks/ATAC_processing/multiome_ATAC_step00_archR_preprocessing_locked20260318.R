library(parallel)
library(ArchR)

Sys.setenv(RHDF5_USE_FILE_LOCKING="FALSE")
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE")

addArchRGenome("mm10")
addArchRThreads(threads = 1)

# JR multiome experiments
experiment_path <- '../../rawdata/multiome'
sample_ids_current <- c('p489c_shRen_caer_48h_multiome_1', 'p489c_shRen_caer_48h_multiome_2')
sample_ids_all <- sample_ids_current
sample_paths <- file.path(experiment_path, sample_ids_current, 'cr-arc-results/atac_fragments.tsv.gz')

# DAC multiome experiments
experiment_path <- '../../rawdata/DAC_multiome'
sample_ids_current <- c('DACE657_mKate2_multiome')
sample_ids_all <- c(sample_ids_all, sample_ids_current)
sample_paths <- c(sample_paths, file.path(experiment_path, sample_ids_current, 'cr-arc-results/atac_fragments.tsv.gz'))

gex_directory <- "../../preprocessed/step0_simplePreprocessing"

# Create Arrow files
ArrowFiles <- createArrowFiles(c(sample_paths), sampleNames = c(sample_ids_all), addTileMat = FALSE, addGeneScoreMat = FALSE, minTSS=1, minFrags=100, maxFrags = 1000000, subThreading = F, force=T)

arrow_filenames = paste(sample_ids_all, ".arrow", sep = "")

proj_name <- "multiome_PP20251110"
proj <- ArchRProject(
  ArrowFiles = arrow_filenames, 
  outputDirectory = proj_name,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
proj@cellColData$replicate_id_simple <- gsub("caer_48h_multiome_", "rep", gsub("p489c_", "", as.character(proj@cellColData$Sample)))
proj <- saveArchRProject(ArchRProj = proj)

write.csv(proj@cellColData, file.path(gex_directory, "multiome_PP20251110_step0_atacQC.csv"), quote=F, row.names=T)

###########################################################
# Here I inspect ArchR QC in python and re-import the cellColData matrix with added annotations
###########################################################

atac_QC <- read.csv(file.path(gex_directory, "multiome_PP20251110_step0_atacQC_custom.csv"), header=T, stringsAsFactors=F)
valid_cells <- atac_QC[,"custom_QC"] == "True" & atac_QC[,"in_gex"] == "True"
proj <- proj[rownames(proj@cellColData) %in% atac_QC[valid_cells, "simplified_name"],]
print(identical(rownames(proj@cellColData), atac_QC[valid_cells, "simplified_name"]))
proj <- saveArchRProject(ArchRProj = proj)

###########################################################
# Add doublet scores
###########################################################

proj <- addTileMatrix(proj, force=TRUE)
doubScores <- addDoubletScores(
    input = getArrowFiles(proj),
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)
ds <- unlist(sapply(doubScores, "[[", "doubletScore"))
names(ds)<- gsub(".*?\\.", "", names(ds))
de <- unlist(sapply(doubScores, "[[", "doubletEnrich"))
names(de)<- gsub(".*?\\.", "", names(de))

proj <- addCellColData(proj, ds, "DoubletScore", names(ds), force=TRUE)
proj <- addCellColData(proj, de, "DoubletEnrichment", names(de), force=TRUE)

atac_QC = proj@cellColData
proj <- saveArchRProject(ArchRProj = proj)

# I don't save the filtered object, and only record what would've been filtered
proj2 <- filterDoublets(ArchRProj = proj)
atac_QC$archr_doublet <- ifelse(rownames(atac_QC) %in% rownames(proj2), "False", "True")
write.csv(atac_QC, file.path(gex_directory, "multiome_PP20251110_step0_atacQC.csv"), quote=F, row.names=T)
