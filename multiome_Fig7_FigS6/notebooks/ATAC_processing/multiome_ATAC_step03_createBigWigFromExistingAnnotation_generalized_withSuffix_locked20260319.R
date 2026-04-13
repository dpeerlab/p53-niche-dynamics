library(parallel)
library(ArchR)
#library(BRGenomics)

Sys.setenv(RHDF5_USE_FILE_LOCKING="FALSE")
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE")

addArchRGenome("mm10")
addArchRThreads(threads = 1)

args <- commandArgs(trailingOnly = T)
#args <- c("multiomeK2shp53shRenIndependentHVG3000_PP20231223", "genotype,cell_type", "500", "ReadsInTSS", "insertion", "1", "147","1", "step10")
proj_name <- args[1]
groupBy <- args[2]
tile_size <- as.numeric(args[3])
normMethod <- args[4]
agg_mode <- args[5]
minFragmentLength <- as.numeric(args[6])
maxFragmentLength <- as.numeric(args[7])
step_size <- as.numeric(args[8])
suffix <- args[9]

proj <- loadArchRProject(proj_name)

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

column_name <- sprintf("%s_tile_%d_mode_%s", groupBy, tile_size, agg_mode)
ArchRProj <- addCellColData(ArchRProj = proj, data = proj@cellColData[,groupBy], name = column_name, cells = rownames(proj@cellColData), force=T)

if(tolower(normMethod) %in% tolower(c("ReadsInTSS", "ReadsInPromoter", "nFrags"))){
    normBy <- getCellColData(ArchRProj = ArchRProj, select = normMethod)
}else{
    normBy <- NULL
}

bwDir1 <- file.path(getOutputDirectory(ArchRProj), "GroupBigWigs")
bwDir <- file.path(getOutputDirectory(ArchRProj), "GroupBigWigs", sprintf("%s_%s", column_name, suffix))

dir.create(bwDir1, showWarnings = FALSE)
dir.create(bwDir, showWarnings = FALSE)

cellsInArrow <- split(rownames(getCellColData(ArchRProj)), stringr::str_split(rownames(getCellColData(ArchRProj)), pattern="\\#", simplify=TRUE)[,1])
availableChr <- ArchR:::.availableSeqnames(head(getArrowFiles(ArchRProj)))
chromLengths <- getChromLengths(ArchRProj)
chromSizes <- getChromSizes(ArchRProj)
reference_ranges = list()
for(k in seq(length(availableChr))){
	current_chr <- availableChr[k]
	print(sprintf("Generating reference ranges for %s", current_chr))
	start_positions <- seq(1,chromLengths[current_chr],tile_size)
	end_positions <- start_positions + tile_size - 1
	reference_ranges[[k]] <- GRanges(rep(current_chr, length(start_positions)), ranges=IRanges(start=start_positions, end=end_positions))
}
names(reference_ranges) <- availableChr

getGroupFragmentsSingleChromosome <- function(ArrowFiles, selected_cells, current_chr){
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
	return(do.call("c", result))
}

consolidate_fragments <- function(frags, agg_mode, reference_ranges){
	if(agg_mode == "insertion"){
        insertion_sites	<- c(start(frags), end(frags))
        frags <- GRanges(rep(seqnames(frags), 2), ranges = IRanges(start = insertion_sites, end=insertion_sites))
	}
    overlaps <- findOverlaps(frags, reference_ranges)
    counts <- reference_ranges[unique(subjectHits(overlaps))]
    overlaps <- findOverlaps(frags, counts)
    coverage <- table(subjectHits(overlaps))
    counts$score <- coverage
	return(counts)
}

aggregate_fragments <- function(frags, reference_ranges){
	overlaps <- findOverlaps(frags, reference_ranges)
	counts <- tapply(frags$score[queryHits(overlaps)], subjectHits(overlaps), sum)
	result <- reference_ranges[as.integer(names(counts))]
	result$score <- counts
	return(result)
}

getGroupFragments <- function(ArrowFiles, selected_cells, agg_mode){
	result <- list()
	for(k in seq_along(availableChr)){
		current_chr <- availableChr[k]
		print(sprintf("Processing %s", current_chr))
		subfragments <- getGroupFragmentsSingleChromosome(ArrowFiles, selected_cells, current_chr)
		if(length(subfragments) > 0)
			result[[k]] <- consolidate_fragments(subfragments, agg_mode, reference_ranges[[k]])
	}
	result <- result[sapply(result, length) > 0]
	return(do.call("c", result))
}

print_bw_file <- function(ArrowFiles, selected_cells, minFragmentLength, maxFragmentLength, normMethod, agg_mode, covFile){
	if(file.exists(covFile)){
        	return()
	}
	print(sprintf("Processing %s", covFile))
	frags <- getGroupFragments(ArrowFiles, selected_cells, agg_mode)
	if(length(frags) == 0){
		print("Didn't find any fragments")
		return(0)
	}
	if(tolower(normMethod) %in% c("readsintss", "readsinpromoter", "nfrags")){
		frags$score <- frags$score * 10^4 / sum(normBy[selected_cells, 1])
	}else if(tolower(normMethod) %in% c("ncells")){
		frags$score <- frags$score / length(selected_cells)
	}
	seqlengths(frags) <- chromLengths[seqlevels(frags)]
	rtracklayer::export.bw(object = frags, con = covFile)
}

ArrowFiles <- getArrowFiles(ArchRProj)
Groups <- getCellColData(ArchRProj = ArchRProj, select = column_name, drop = TRUE)
Cells <- ArchRProj$cellNames
AllCellGroups <- split(Cells, Groups)

valid_groups <- unique(names(AllCellGroups))
valid_groups <- c("progenitor", "gastric")
group_subset <- c()
for(my_group in valid_groups){
	group_subset <- c(group_subset, names(AllCellGroups)[grep(my_group, names(AllCellGroups))])
}
cellGroups <- AllCellGroups[group_subset]

fragment_start <- seq(minFragmentLength,maxFragmentLength-step_size,step_size)
fragment_end <- seq(minFragmentLength + step_size - 1,maxFragmentLength,step_size)

for(current_group in seq_along(cellGroups)){
	selected_group <- names(cellGroups)[current_group]
	selected_cells <- cellGroups[[current_group]]
	for(i in seq_along(fragment_start)){
		minFragmentLength <- fragment_start[i]
		maxFragmentLength <- fragment_end[i]
		covFile <- file.path(bwDir, paste0(make.names(selected_group), "-TileSize-",tile_size,"-fragLen-",minFragmentLength,"-",maxFragmentLength,"-normMethod-",normMethod,"-AggMode-",agg_mode,"-ArchR.bw"))
		print_bw_file(ArrowFiles, selected_cells, minFragmentLength, maxFragmentLength, normMethod, agg_mode, covFile)
	}
}
