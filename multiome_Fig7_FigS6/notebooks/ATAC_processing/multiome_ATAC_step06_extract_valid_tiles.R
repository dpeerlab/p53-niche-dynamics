library(parallel)
library(ArchR)

Sys.setenv(RHDF5_USE_FILE_LOCKING="FALSE")
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE")

addArchRGenome("mm10")
addArchRThreads(threads = 1)

args <- commandArgs(trailingOnly = T)
#args <- c("multiomeK2shp53shRenIndependentHVG3000_PP20231223", "seacell_tile_100_mode_insertion")
ArchRProj <- args[1]
grouping <- args[2]
threshold <- 0.05
extension <- 20

proj <- loadArchRProject(ArchRProj)
bigwig_dir <- file.path(getOutputDirectory(proj), "GroupBigWigs", grouping)
print(bigwig_dir)

bw_files <- list.files(bigwig_dir)
bw_files <- bw_files[grep("\\.bw", bw_files, perl=TRUE)]
bw_files <- bw_files[grep("fragLen-1-147", bw_files)]

for(i in seq_along(bw_files)){
    my_file = bw_files[i]
    print(my_file)
    tryCatch(
        {
            current_bw <- rtracklayer::import.bw(file.path(bigwig_dir, my_file))
            current_bw <- current_bw[current_bw$score > threshold]
            if(my_file == bw_files[1]){
                result <- current_bw
            }else{
                result <- c(result, current_bw[!(seq(length(current_bw)) %in% queryHits(findOverlaps(current_bw, result, minoverlap=2L)))])
            }
        }, error = function(cond){
            message(paste("Problematic file:", my_file))
            message("Here's the original error message:")
            message(conditionMessage(cond))
            NA
        }
    )
}
result <- sortSeqlevels(result)
result <- sort(result)

gene_annotation <- getGeneAnnotation(proj)
result <- ArchR:::.fastAnnoPeaks(result, "mm10", gene_annotation)

dataframe_to_GRanges <- function(current_df){
    region_names <- sprintf("%s:%d-%d", current_df$seqnames, current_df$start, current_df$end)
    current_ranges <- GRanges(region_names)
    strand(current_ranges) <- current_df$strand
    columns <- colnames(current_df)
    columns <- columns[!(columns %in% c("seqnames", "start", "end", "strand"))]
    mcols(current_ranges) <- current_df[,columns]
    return(current_ranges)
}

GRanges_to_dataframe <- function(current_ranges){
    seqnames <- as.character(seqnames(current_ranges))
    start <- start(current_ranges)
    end <- end(current_ranges)
    strand <- as.character(strand(current_ranges))
    metadata <- mcols(current_ranges)
    region <- sprintf("%s:%d-%d", seqnames, start, end)
    result <- data.frame(seqnames = seqnames, start = start, end = end, region = region, strand = strand)
    for(mycol in colnames(metadata)){
        result[[mycol]] <- as.character(metadata[[mycol]])
    }
    return(result)
}

print("Writing file")
df <- GRanges_to_dataframe(result)
write.table(df, file.path(getOutputDirectory(proj), "TileData", sprintf("Tile_%s.csv", grouping)), sep = ",", row.names=F, quote=F, col.names=T)
write.table(df[,1:6], file.path(getOutputDirectory(proj), "TileData", sprintf("Tile_%s.bed", grouping)), sep = "\t", row.names=F, quote=F, col.names=T)

df$start <- df$start - extension
df$end <- df$end + extension
df$region <- sprintf("%s:%d-%d", df$seqnames, df$start, df$end)
write.table(df, file.path(getOutputDirectory(proj), "TileData", sprintf("Tile_%s_extended_%d.csv", grouping, extension)), sep = ",", row.names=F, quote=F, col.names=T)
write.table(df[,1:6], file.path(getOutputDirectory(proj), "TileData", sprintf("Tile_%s_extended_%d.bed", grouping, extension)), sep = "\t", row.names=F, quote=F, col.names=T)
