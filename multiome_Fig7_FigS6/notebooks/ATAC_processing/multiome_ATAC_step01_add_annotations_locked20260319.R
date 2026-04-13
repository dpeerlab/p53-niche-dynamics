###########################################################
# Start here after finishing processing gene expression.
# This is important because cells can still be lost during SCRAN normalization
# Thus this step ensures that both modalities will have the same set of cells
# I also make a copy of the base project to consider only a subset of conditions
# This subset of cells is defined during processing of scRNA-seq at the UMAP/seacells steps
###########################################################

library(parallel)
library(ArchR)

Sys.setenv(RHDF5_USE_FILE_LOCKING="FALSE")
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE")

addArchRGenome("mm10")
addArchRThreads(threads = 1)

proj_name_original <- "multiome_PP20251110"
proj <- loadArchRProject(proj_name_original)

seacells_path <- "../intermediate_output/step3_seacells/Reyes_multiome_consolidated__SEACells.csv"
seacells <- read.csv(seacells_path, header = T)
rownames(seacells) <- seacells$X
seacells <- seacells[rownames(proj@cellColData),]
print(identical(rownames(proj@cellColData), rownames(seacells)))

selected_columns = c("DC_gastric_progenitor", "SEACell", "pheno_10", "pheno_30", "batch", "prefix", "condition_id", "condition", "cell_type", "DC_bin")

for(my_column in selected_columns){
	proj <- addCellColData(ArchRProj = proj, data = seacells[,my_column], name = my_column, cells = rownames(seacells), force=T)
}
proj <- saveArchRProject(proj)
