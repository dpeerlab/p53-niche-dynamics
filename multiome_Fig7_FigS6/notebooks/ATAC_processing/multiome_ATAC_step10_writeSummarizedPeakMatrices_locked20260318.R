library(parallel)
library(ArchR)
library(reshape2)
library(Matrix.utils)

library(BSgenome.Mmusculus.UCSC.mm10)
mm10 <- BSgenome.Mmusculus.UCSC.mm10

Sys.setenv(RHDF5_USE_FILE_LOCKING="FALSE")
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE")

addArchRGenome("mm10")
addArchRThreads(threads = 1)

args <- commandArgs(trailingOnly=T)
proj_name <- args[1]
proj <- loadArchRProject(proj_name)

motif_output <- file.path(proj_name, "MotifData", "motifmatchr")
peak_output <- file.path(proj_name, "PeakData")
motif_instances <- file.path(motif_output, "MotifInstances")
if(!dir.exists(file.path(proj_name, "MotifData"))){
    dir.create(file.path(proj_name, "MotifData"))
}
if(!dir.exists(motif_output)){
    dir.create(motif_output)
}
if(!dir.exists(motif_instances)){
    dir.create(motif_instances)
}
if(!dir.exists(peak_output)){
    dir.create(peak_output)
}
if(!dir.exists(file.path(proj_name, "TileData"))){
    dir.create(file.path(proj_name, "TileData"))
}

####################################
# Extract peaks
####################################

peaks = getPeakSet(proj)
peaks$region = paste(as.character(seqnames(peaks)), ":", as.character(ranges(peaks)), sep = "")
chr_order <- sort(seqlevels(peaks))
reordered_features <- list()
for(chr in chr_order)
    reordered_features[[chr]] = peaks[seqnames(peaks) == chr]
for(chr in chr_order)
    reordered_features[[chr]] = which(seqnames(peaks) == chr)
reordered_features <- Reduce("c", reordered_features)

peaks = peaks[reordered_features,]

get_bed_file <- function(my_peaks){
   starts = start(ranges(my_peaks))
   ends = end(ranges(my_peaks))
   seqname = seqnames(my_peaks)
   return(data.frame(seqname, starts, ends, paste(seqname, paste(starts, ends, sep = "-"), sep = ":")))
}

########################################################
# Write cross-compatible summarized peak matrices
########################################################

print("Writing cross-compatible summarized peak matrices")

peak_matrix <- getMatrixFromProject(proj,  useMatrix="PeakMatrix")
pm <- assay(peak_matrix, "PeakMatrix", useDimNames=T)

coldata = colnames(proj)
potential_clusters = coldata[grep("seacell\\d+", coldata, perl=T)]
print(potential_clusters)

for(column_name in potential_clusters){
	pm_agg <- aggregate.Matrix(t(pm), colData(peak_matrix)[,column_name])
	normBy <- tapply(colData(peak_matrix)$ReadsInTSS, colData(peak_matrix)[,column_name], sum)
	pm_agg <- apply(pm_agg, 2, "/", normBy) * 10^4
	colnames(pm_agg) = peaks$region
	write.csv(data.frame(t(pm_agg)), file.path(proj_name, "PeakData", sprintf("PeakMatrix_%s.csv", column_name)), quote=F)
}
