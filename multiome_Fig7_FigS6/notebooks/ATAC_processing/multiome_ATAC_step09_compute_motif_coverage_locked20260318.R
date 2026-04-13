################################################################
# Compute motif coverage per SEACell. Coverage is defined as the number of fragments that contain partly or fully a motif 
################################################################
library(ArchR)

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

getGroupFragmentsSingleChromosome <- function(ArrowFiles, selected_cells, current_chr, target_regions){
	result <- list()
	for(j in seq_along(ArrowFiles)){
		current_arrow <- names(ArrowFiles)[j]
		valid_cells <- selected_cells %in% cellsInArrow[[current_arrow]]
		if(sum(valid_cells) == 0){
			#print("Didn't find valid cells in this arrow file")
			next
		}
		current_cells <- selected_cells[valid_cells]
		result[[j]] <- ArchR:::.getFragsFromArrow(ArrowFiles[[current_arrow]], chr = current_chr, out = "GRanges", cellNames = selected_cells, maxFragmentLength = maxFragmentLength, minFragmentLength=minFragmentLength)
	}
	result <- result[sapply(result, length) > 0]
    all_fragments <- do.call("c", result)
    if(length(all_fragments) > 0){
        overlaps <- findOverlaps(all_fragments, target_regions)
        fragments_subset <- all_fragments[unique(queryHits(overlaps))]
    }else{
        fragments_subset <- all_fragments
    }
	return(fragments_subset)
}

getMotifCoverage <- function(ArrowFiles, selected_cells, target_regions){
	result <- list()
	for(k in seq_along(availableChr)){
		current_chr <- availableChr[k]
		#print(sprintf("Processing %s", current_chr))
		result[[k]] <- getGroupFragmentsSingleChromosome(ArrowFiles, selected_cells, current_chr, target_regions)
	}
	result <- result[sapply(result, length) > 0]
    all_fragments <- do.call("c", result)
    overlaps <- findOverlaps(all_fragments, target_regions)
    target_coverage <- table(subjectHits(overlaps))
	return(target_coverage)
}

args <- commandArgs(trailingOnly = T)
#args <-c("multiome_PP20251110", "SEACell", "multiome_PP20251110/MotifData/MotifInstances/DATABASE_homer_MOTIF_Trp53_canonical.tsv.gz", 1, 147)
proj_name <- args[1]
groupBy <- args[2]
motif_file <- args[3]
minFragmentLength <- as.numeric(args[4])
maxFragmentLength <- as.numeric(args[5])
motif_name <- gsub("\\.tsv.gz", "", basename(motif_file))

# Load ArchR project and define general metadata
proj <- loadArchRProject(proj_name)
ArrowFiles <- getArrowFiles(proj)
cellsInArrow <- split(
    rownames(getCellColData(proj)), stringr::str_split(rownames(getCellColData(proj)), pattern="\\#", simplify=TRUE)[,1]
)
availableChr <- ArchR:::.availableSeqnames(head(getArrowFiles(proj)))

# Define grouping
proj@cellColData$Sample	<- as.character(proj@cellColData$Sample)
if(length(grep(",", groupBy)) > 0){
	selected_fields <- unlist(stringr::str_split(groupBy, pattern = ","))
	selected_annotations <- proj@cellColData[,selected_fields]
	selected_annotations <- apply(selected_annotations, 2, as.character)
	groupBy <- paste(selected_fields, collapse = "_")
	chimeric_annotations = apply(selected_annotations, 1, paste, collapse = "_")
	proj <- addCellColData(proj, data = chimeric_annotations, name = groupBy, cells = rownames(proj@cellColData))
}
proj@cellColData[,groupBy] <- as.character(proj@cellColData[,groupBy])
ArchRProj <- addCellColData(ArchRProj = proj, data = proj@cellColData[,groupBy], name = groupBy, cells = rownames(proj@cellColData), force=T)

possible_groups <- unique(proj@cellColData[,groupBy])
valid_groups <- unique(proj@cellColData$cell_type)
group_subset <- c()
for(my_group in valid_groups){
	group_subset <- c(group_subset, possible_groups[grep(my_group, possible_groups)])
}

# Load motifs
motif_df <- read.table(motif_file, header = T, sep = "\t")
motif_df$motif_position <- paste(motif_df$seqnames, ":", motif_df$start, "-", motif_df$end, sep = "")
current_motifs <- dataframe_to_GRanges(motif_df)

# Get motif coverage in every group, then identify maximum coverage
result = list()
normalization_factor <- c()
for(i in seq_along(group_subset)){
    selected_group <- group_subset[i]
    if(selected_group %in% names(result)){
        next
    }
    print(sprintf("Processing %d of %d", i, length(group_subset)))
    selected_cells <- rownames(proj@cellColData)[proj@cellColData[,groupBy] == selected_group]
    motif_fragments <- getMotifCoverage(ArrowFiles, selected_cells, current_motifs)
    normalization_factor <- c(normalization_factor, sum(proj@cellColData[selected_cells,]$ReadsInTSS))
    result[[selected_group]] <- data.frame(motif_index = as.numeric(names(motif_fragments)), motif_coverage = as.numeric(motif_fragments))
}
motif_coverage <- data.frame(do.call("rbind", result))
colnames(motif_coverage) <- c("motif_index", "coverage")
motif_coverage$group_name <- gsub("\\.\\d+$", "", rownames(motif_coverage))
motif_coverage$group_index <- match(motif_coverage$group_name, group_subset)
motif_coverage <- motif_coverage[,c("group_index", "motif_index", "coverage")]

# Output path
output_path <- file.path(proj_name, "MotifCoverage", groupBy, sprintf("%s_GROUP_%s_FRAG_%d_%d", motif_name, groupBy, minFragmentLength, maxFragmentLength))
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive=T)
}
# Write group coverage tables
write.table(motif_coverage, gzfile(file.path(output_path, "coverage.tsv.gz")), row.names=F, quote=F, sep = "\t")
group_df <- data.frame(group = group_subset)
write.table(group_df, gzfile(file.path(output_path, "group.tsv.gz")), row.names=F, quote=F, sep = "\t")
write.table(motif_df[,c("seqnames", "start", "end", "score")], gzfile(file.path(output_path, "motif.tsv.gz")), row.names=F, quote=F, sep = "\t")

# Max coverage per motif
max_coverage <- tapply(motif_coverage$coverage, motif_coverage$motif_index, max)
max_coverage_df <- data.frame(motif_index = names(max_coverage), max_coverage = max_coverage)
write.table(max_coverage_df, gzfile(file.path(output_path, "max_coverage.tsv.gz")), row.names=F, quote=F, sep = "\t")

# Normalization
normalization_df <- data.frame(group = group_subset, normalization_factor = normalization_factor)
write.table(normalization_df, gzfile(file.path(output_path, "normalization_factor.tsv.gz")), row.names=F, quote=F, sep = "\t")
