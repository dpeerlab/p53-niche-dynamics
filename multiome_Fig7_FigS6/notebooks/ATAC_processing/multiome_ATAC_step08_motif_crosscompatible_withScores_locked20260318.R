library(parallel)
library(ArchR)
library(reshape2)
library(rhdf5)

library(BSgenome.Mmusculus.UCSC.mm10)
mm10 <- BSgenome.Mmusculus.UCSC.mm10

Sys.setenv(RHDF5_USE_FILE_LOCKING="FALSE")
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE")

addArchRGenome("mm10")
addArchRThreads(threads = 1)

GRanges_to_dataframe <- function(current_ranges){
    seqnames <- as.character(seqnames(current_ranges))
    start <- start(current_ranges)
    end <- end(current_ranges)
    strand <- as.character(strand(current_ranges))
    metadata <- mcols(current_ranges)
    result <- data.frame(seqnames = seqnames, start = start, end = end, strand = strand)
    for(mycol in colnames(metadata)){
        result[[mycol]] <- as.character(metadata[[mycol]])
    }
    #result <- cbind(result, metadata)
    return(result)
}
dataframe_to_GRanges <- function(current_df){
    region_names <- sprintf("%s:%d-%d", current_df$seqnames, current_df$start, current_df$end)
    current_ranges <- GRanges(region_names)
    strand(current_ranges) <- current_df$strand
    columns <- colnames(current_df)
    columns <- columns[!(columns %in% c("seqnames", "start", "end", "strand"))]
    mcols(current_ranges) <- current_df[,columns]
    return(current_ranges)
}

proj_name <- "multiome_PP20251110"
proj <- loadArchRProject(proj_name)

motif_output <- file.path(proj_name, "MotifData")
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
write.table(peaks, gzfile(file.path(peak_output, "PeakData.tsv.gz")), quote=F, sep = "\t", row.names=F)
write.table(get_bed_file(peaks), file.path(peak_output, "PeakData.bed"), quote=F, row.names=F, col.names=F, sep = "\t")

#######################################################
# Write cross-compatible motif matrices
#######################################################

print("Writing cross-compatible motif matrices")

getMaxScorePerRegion <- function(ranges, regions){
    overlaps <- findOverlaps(ranges, regions)
    max_score <- tapply(ranges$score[queryHits(overlaps)], subjectHits(overlaps), max)
    result <- regions[as.numeric(names(max_score))]
    result$score <- max_score
    result <- GRanges_to_dataframe(result)
    result <- result[,c("region", "score")]
    return(result)
}

my_motif_set <- "Motif_cisbp"
pos <- getPositions(proj, name = my_motif_set) ## Get motif positions

########################
# Write motif instances
########################

for(my_motif in names(pos)){
    current_motifs <- pos[[my_motif]]
    df <- GRanges_to_dataframe(current_motifs)
    output_file <- file.path(motif_instances, sprintf("DATABASE_%s_MOTIF_%s.tsv.gz", my_motif_set, my_motif))
    write.table(df, gzfile(output_file), row.names = F, quote=F, sep = "\t")
}

########################
# Write p53 motif instances
########################
p53_motif_path = file.path(proj_name, "MotifData", "homer_tile_based")
p53_possible_files = list.files(p53_motif_path)
selected_p53_file = p53_possible_files[grep("_merged.bed", p53_possible_files)]
p53_motifs <- read.table(file.path(p53_motif_path, selected_p53_file), sep = "\t", header = T, stringsAsFactors=F)

p53_canonical <- p53_motifs[p53_motifs$gap == 0,]
output_file <- file.path(motif_instances, sprintf("DATABASE_%s_MOTIF_%s.tsv.gz", "homer", "Trp53_canonical"))
write.table(p53_canonical, gzfile(output_file), row.names = F, quote=F, sep = "\t")

p53_noncanonical <- p53_motifs[p53_motifs$gap != 0,]
output_file <- file.path(motif_instances, sprintf("DATABASE_%s_MOTIF_%s.tsv.gz", "homer", "Trp53_noncanonical"))
write.table(p53_noncanonical, gzfile(output_file), row.names = F, quote=F, sep = "\t")

########################
# Get maximum score per peak and write files
########################

# Import curated p53 motifs and format them to be compatible with ArchR output

p53_canonical <- p53_canonical[,c("seqnames", "start", "end", "strand", "score")]
p53_canonical <- dataframe_to_GRanges(p53_canonical)
overlaps <- findOverlaps(p53_canonical, peaks)
p53_canonical <- p53_canonical[queryHits(overlaps)]
colnames(mcols(p53_canonical)) <- "score"

p53_noncanonical <- p53_noncanonical[,c("seqnames", "start", "end", "strand", "score")]
p53_noncanonical <- dataframe_to_GRanges(p53_noncanonical)
overlaps <- findOverlaps(p53_noncanonical, peaks)
p53_noncanonical <- p53_noncanonical[queryHits(overlaps)]
colnames(mcols(p53_noncanonical)) <- "score"

# Add peak annotations to motif set
cisbp_matrix <- lapply(pos, getMaxScorePerRegion, peaks)
cisbp_matrix <- do.call("rbind", cisbp_matrix)
cisbp_matrix$motif_name <- gsub("\\.\\d+", "", rownames(cisbp_matrix))

p53_matrix_canonical <- getMaxScorePerRegion(p53_canonical, peaks)
p53_matrix_canonical$motif_name <- "Trp53_canonical"
p53_matrix_noncanonical <- getMaxScorePerRegion(p53_noncanonical, peaks)
p53_matrix_noncanonical$motif_name <- "Trp53_noncanonical"

mat1 <- cisbp_matrix
mat2 <- cisbp_matrix[-grep("Trp53", cisbp_matrix$motif_name),]
mat2 <- rbind(mat2, p53_matrix_canonical, p53_matrix_noncanonical)

write.table(data.frame(mat1), gzfile(file.path(motif_output, "PeakMotifScoreMatrix_cisbp.csv.gz")), row.names=F, quote=F, sep = ",")
write.table(data.frame(mat2), gzfile(file.path(motif_output, "PeakMotifScoreMatrix_cisbp_p53.csv.gz")), row.names=F, quote=F, sep = ",")
