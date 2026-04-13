library(GenomicRanges)
library(dplyr)

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

# Define input dataset
selected_dataset <- "multiome_PP20251110"
masking <- "unmasked"
input_path <- file.path(selected_dataset, "MotifData", "homer_tile_based")

# Define selected half site from Encode
p53_motif_file <- sprintf("TileData_homerp53_%s.bed", masking)
selected_p53_motif <- "TP53_half_HH"
motif_half <- read.table(file.path(input_path, p53_motif_file), sep = "\t", skip=1, header=F)
motif_half <- motif_half[grep(selected_p53_motif, motif_half[,4]),]
colnames(motif_half) <- c("seqnames", "start", "end", "MotifName", "score", "strand")

# These ranges are used as positive controls for known p53 motifs in canonical targets
cdkn1a_motif_ranges = GRanges(c("chr17:29090918-29090956"))
mdm2_motif_ranges = GRanges(c("chr10:117710096-117710134"))

###############################################
# In this section, I classify p53 sites in tiers
###############################################

# Correct 0-indexing from Perl to be consistent with R
motif_half$start <- motif_half$start + 1
reference_ranges <- dataframe_to_GRanges(motif_half)

# First find side-by-side pairs and remove query range from overlaps list
distal <- findOverlaps(reference_ranges, reference_ranges, maxgap = 11L)
start_distance <- abs(start(reference_ranges)[queryHits(distal)] - start(reference_ranges)[subjectHits(distal)])
adjacent_motif <- start_distance == 10
paired_adjacent <- unique(queryHits(distal)[adjacent_motif])
adjacent_pairs <- distal[(queryHits(distal) %in% paired_adjacent) & adjacent_motif]

# Next find distal pairs and remove query range from overlaps list
distal <- distal[!(queryHits(distal) %in% paired_adjacent)]
start_distance <- abs(start(reference_ranges)[queryHits(distal)] - start(reference_ranges)[subjectHits(distal)])
distal_motif <- start_distance > 10
paired_distal <- unique(queryHits(distal)[distal_motif])
distal_pairs <- distal[(queryHits(distal) %in% paired_distal) & distal_motif]

# Next find overlapping hits and remove query range from overlaps list
distal <- distal[!(queryHits(distal) %in% paired_distal)]
start_distance <- abs(start(reference_ranges)[queryHits(distal)] - start(reference_ranges)[subjectHits(distal)])
overlapping_motif <- start_distance < 10 & start_distance > 0
paired_overlapping <- unique(queryHits(distal)[overlapping_motif])
overlapping_pairs <- distal[(queryHits(distal) %in% paired_overlapping) & overlapping_motif]

# Find remaining ranges without another half-site nearby
distal <- distal[!(queryHits(distal) %in% paired_overlapping)]
start_distance <- abs(start(reference_ranges)[queryHits(distal)] - start(reference_ranges)[subjectHits(distal)])
singleton_motif <- start_distance == 0
paired_singleton <- unique(queryHits(distal)[singleton_motif])
singleton_pairs <- distal[(queryHits(distal) %in% paired_singleton) & singleton_motif]

#################################################################
# In this section, I consolidate motifs to create composite sites
#################################################################
consolidate_ranges <- function(overlaps){
    start1 <- start(reference_ranges[queryHits(overlaps)])
    start2 <- start(reference_ranges[subjectHits(overlaps)])
    end1 <- end(reference_ranges[queryHits(overlaps)])
    end2 <- end(reference_ranges[subjectHits(overlaps)])
    # For each pair of half motifs, find the earliest start and latest end to form a composite motif
    selected_start <- rep(-1, length(start1))
    selected_start[start1 <= start2] <- start1[start1 <= start2]
    selected_start[start1 > start2] <- start2[start1 > start2]
    selected_end <- rep(-1, length(end1))
    selected_end[end1 <= end2] <- end2[end1 <= end2]
    selected_end[end1 > end2] <- end1[end1 > end2]
    selected_chr <- seqnames(reference_ranges[queryHits(overlaps)])
    selected_strand <- strand(reference_ranges[queryHits(overlaps)])
    # Aggregate motif scores
    selected_score <- as.numeric(reference_ranges[queryHits(overlaps)]$score) + as.numeric(reference_ranges[subjectHits(overlaps)]$score)
    new_ranges <- GRanges(seqnames = selected_chr, ranges = IRanges(start = selected_start, end=selected_end), strand = selected_strand)
    new_ranges$score <- selected_score
    new_ranges$gap <- end(new_ranges) - start(new_ranges) + 1 - 20
    return(new_ranges)
}
print("Processing adjacent ranges")
adjacent_ranges <- consolidate_ranges(adjacent_pairs)
print(adjacent_ranges[queryHits(findOverlaps(adjacent_ranges, cdkn1a_motif_ranges))])

print("Processing distal ranges")
distal_ranges <- consolidate_ranges(distal_pairs)

print("Processing overlapping ranges")
overlapping_ranges_temp <- consolidate_ranges(overlapping_pairs)
overlapping_ranges <- reference_ranges[queryHits(overlapping_pairs)]
overlapping_ranges$gap <- overlapping_ranges_temp$gap
mcols(overlapping_ranges) <- mcols(overlapping_ranges)[,c("score", "gap")]

print("Processing singleton ranges")
singleton_ranges <- reference_ranges[queryHits(singleton_pairs)]
singleton_ranges$gap <- -10
mcols(singleton_ranges) <- mcols(singleton_ranges)[,c("score", "gap")]

#######################################################################################
# Import p53 binding data to filter non-canonical motifs
#######################################################################################

# Import ChIP-seq results from Kenzelmann-Broz et al.
dataset_path <- "/data1/lowes/reyesj3/databases/p53_ChIP_seq/KenzelmannBroz_GSE46240/GSE46240_ChIP-seq_peaks_p53_binding_sites_mm10.bed"
p53_chip <- read.table(dataset_path, header = F, stringsAsFactors = F, sep = "\t")
colnames(p53_chip) <- c("seqnames", "start", "end", "peak_name", "score", "strand")
p53_chip$start <- p53_chip$start -100
p53_chip$end <- p53_chip$end +100
p53_chip$strand <- "*"
p53_chip_ranges <- dataframe_to_GRanges(p53_chip)

# Compute overlaps with singleton and distal motif sites
overlaps <- findOverlaps(singleton_ranges, p53_chip_ranges)
print(sprintf("Original number of singleton motifs: %d", length(singleton_ranges)))
singleton_ranges <- singleton_ranges[queryHits(overlaps)]
print(sprintf("Filtered number of singleton motifs: %d", length(singleton_ranges)))

overlaps <- findOverlaps(distal_ranges, p53_chip_ranges)
print(sprintf("Original number of distal motifs: %d", length(distal_ranges)))
distal_ranges <- distal_ranges[queryHits(overlaps)]
print(sprintf("Filtered number of distal motifs: %d", length(distal_ranges)))

#######################################################################################
# Consolidate final motif list
#######################################################################################

# Excluded overlapping ranges to keep highest confidence motifs
final_ranges <- c(adjacent_ranges, distal_ranges, singleton_ranges)
print(final_ranges[queryHits(findOverlaps(final_ranges, cdkn1a_motif_ranges))])

final_ranges_df <- GRanges_to_dataframe(final_ranges)
final_ranges_df$range_name <- paste(seqnames(final_ranges), ":", start(final_ranges), "-", end(final_ranges), sep = "")
# Remove redundancies, selecting motif/strand with highest score when there are competing motifs with exactly the same site
final_ranges_df <- unique(final_ranges_df)
final_ranges_df <- slice_max(final_ranges_df, order_by= score, by=range_name)
# This last filter is for cases in which the motif has the exact same score in two strands
final_ranges_df <- final_ranges_df[!duplicated(final_ranges_df$range_name),]
#head(final_ranges)

# Write final table
final_ranges_df$MotifName <- sprintf("p53_HH_consolidated_%s_tile_based", masking)
final_ranges_df <- final_ranges_df[,c('seqnames','start','end','MotifName','score','strand','gap')]
write.table(final_ranges_df, file.path(input_path, sprintf("TileData_homerp53_%s_tile_based_HH_consolidated.bed", masking)), row.names = F, quote=F, sep="\t")

################################################################
# Finally, I write the actual full site that was found by homer
################################################################
p53_motif_file <- sprintf("TileData_homerp53_%s.bed", masking)
selected_p53_motif <- "p53_full_HH_corrected"
motif_full <- read.table(file.path(input_path, p53_motif_file), sep = "\t", skip=1, header=F)
motif_full <- motif_full[grep(selected_p53_motif, motif_full[,4]),]
motif_full$MotifName <- sprintf("p53_full_HH_corrected_%s_tile_based", masking)
colnames(motif_full) <- c("seqnames", "start", "end", "motif_id", "score", "strand", "MotifName")
motif_full$range_name <- paste(motif_full$seqnames, ":", motif_full$start, "-", motif_full$end, sep = "")
motif_full <- unique(motif_full)
motif_full <- slice_max(motif_full, order_by= score, by=range_name)
# This last filter is for cases in which the motif has the exact same score in two strands
motif_full <- motif_full[!duplicated(motif_full$range_name),]

# Correct 0-indexing from Perl to be consistent with R
motif_full$start <- motif_full$start + 1
motif_full$gap <- 0

motif_full <- motif_full[,c('seqnames','start','end','MotifName','score','strand','gap')]
write.table(motif_full, file.path(input_path, sprintf("TileData_homerp53_%s_tile_based_HH_full.bed", masking)), row.names=F, col.names=T, sep = "\t", quote=F)

######################################################
# Consolidate chimeric motifs and homer full motif scan
######################################################

# Annotate motif type
motif_full$motif_type <- "canonical"
# This annotation is not completely correct, as some non-canonical distal elements are in reality canonical (gap = 0). This doesn't affect downstream analysis, since I re-classify later on based on the gap information
final_ranges_df$motif_type <- ifelse(final_ranges_df$gap < 0, "noncanonical_singleton", "noncanonical_distal")

motif_consolidated <- rbind(motif_full, final_ranges_df)
motif_consolidated$range_name <- paste(motif_consolidated$seqnames, ":", motif_consolidated$start, "-", motif_consolidated$end, sep = "")
motif_consolidated <- slice_max(motif_consolidated, order_by= score, by=range_name)
motif_consolidated$original_name <- motif_consolidated$MotifName
motif_consolidated$MotifName <- "Trp53_custom"

# Identify motifs that are bound
consolidated_ranges <- dataframe_to_GRanges(motif_consolidated)
overlaps <- findOverlaps(consolidated_ranges, p53_chip_ranges)
consolidated_ranges$p53_ChIP_200 <- seq(length(consolidated_ranges)) %in% queryHits(overlaps)

# Write final table
motif_consolidated <- GRanges_to_dataframe(consolidated_ranges)
write.table(motif_consolidated, file.path(input_path, sprintf("TileData_homerp53_%s_tile_based_HH_custom_merged.bed", masking)), row.names=F, col.names=T, sep = "\t", quote=F)
