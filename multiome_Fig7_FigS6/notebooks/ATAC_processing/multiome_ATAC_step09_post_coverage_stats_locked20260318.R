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

load_coverage_information <- function(coverage_path){
    coverage <- read.table(gzfile(file.path(coverage_path, "coverage.tsv.gz")), header = T, sep = "\t")
    group_info <- read.table(gzfile(file.path(coverage_path, "group.tsv.gz")), header = T, sep = "\t")
    motif_info <- read.table(gzfile(file.path(coverage_path, "motif.tsv.gz")), header = T, sep = "\t")
    motif_info$motif_name <- paste(motif_info$seqnames, ":", motif_info$start, "-", motif_info$end, sep = "")
    coverage$motif_position <- motif_info$motif_name[coverage$motif_index]
    coverage$group_name <- group_info$group[coverage$group_index]
    return(coverage)
}

load_coverage_information_asDataFrame <- function(coverage_path){
    coverage <- read.table(gzfile(file.path(coverage_path, "coverage.tsv.gz")), header = T, sep = "\t")
    group_info <- read.table(gzfile(file.path(coverage_path, "group.tsv.gz")), header = T, sep = "\t")
    motif_info <- read.table(gzfile(file.path(coverage_path, "motif.tsv.gz")), header = T, sep = "\t")
    motif_info$motif_name <- paste(motif_info$seqnames, ":", motif_info$start, "-", motif_info$end, sep = "")
    
    sparse_matrix <- sparseMatrix(
      i = coverage$motif_index,
      j = coverage$group_index,
      x = coverage$coverage,
      dims = c(max(coverage$motif_index), max(coverage$group_index)) # Optional: specify dimensions
    )
    df <- as.data.frame(as.matrix(sparse_matrix))
    rownames(df) <- motif_info$motif_name[1:dim(df)[1]]
    colnames(df) <- group_info$group[1:dim(df)[2]]
    
    return(df)
}

args <- commandArgs(trailingOnly = T)
#args <-c("multiome_PP20251110", "Sample_cell_type", "multiome_PP20251110/MotifData/MotifInstances/DATABASE_homer_MOTIF_Trp53_canonical.tsv.gz", 1, 147)
proj_name <- args[1]
groupBy <- args[2]
motif_file <- args[3]
minFragmentLength <- as.numeric(args[4])
maxFragmentLength <- as.numeric(args[5])
motif_name <- gsub("\\.tsv.gz", "", basename(motif_file))
output_path <- file.path(proj_name, "MotifCoverage", groupBy, sprintf("%s_GROUP_%s_FRAG_%d_%d", motif_name, groupBy, minFragmentLength, maxFragmentLength))

# Read peak data
peak_data <- read.table(gzfile(file.path(proj_name, "PeakData", "PeakData.tsv.gz")), stringsAsFactors = T, sep = "\t", header = T)
peak_ranges <- GRanges(peak_data$region)

# Read motifs
motif_df <- read.table(motif_file, header = T, sep = "\t")
motif_df$motif_position <- paste(motif_df$seqnames, ":", motif_df$start, "-", motif_df$end, sep = "")
current_motifs <- dataframe_to_GRanges(motif_df)

# Filter motifs to find only those in peaks
selected_motifs <- current_motifs[unique(queryHits(findOverlaps(current_motifs, peak_ranges, type='within')))]

# Read coverage files
motif_coverage <- load_coverage_information_asDataFrame(output_path)
normalization <- read.table(file.path(output_path, "normalization_factor.tsv.gz"), header = T)
motif_coverage_normalized <- data.frame(t(apply(motif_coverage, 1, "/", normalization$normalization_factor)) * 10000)
rownames(motif_coverage_normalized) <- rownames(motif_coverage)
colnames(motif_coverage_normalized) <- colnames(motif_coverage)
max_coverage <- apply(motif_coverage_normalized, 1, max)
max_coverage <- max_coverage[selected_motifs$motif_position]

ACCESSIBILITY_THRESHOLD = 0.05153401921440394
print(sprintf("Fraction of motifs in peaks with coverage lower than %.4f = %.2f", ACCESSIBILITY_THRESHOLD, mean(max_coverage < ACCESSIBILITY_THRESHOLD)))

ACCESSIBILITY_THRESHOLD = 0.013
print(sprintf("Fraction of motifs in peaks with coverage lower than %.4f = %.2f", ACCESSIBILITY_THRESHOLD, mean(max_coverage < ACCESSIBILITY_THRESHOLD)))

ACCESSIBILITY_THRESHOLD = 0.008946699265574468
print(sprintf("Fraction of motifs in peaks with coverage lower than %.4f = %.2f", ACCESSIBILITY_THRESHOLD, mean(max_coverage < ACCESSIBILITY_THRESHOLD)))

print(sprintf("Fraction of motifs in peaks without coverage = %.2f", mean(max_coverage == 0)))
