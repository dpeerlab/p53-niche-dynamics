library(parallel)
library(ArchR)
library(reshape2)
library(rhdf5)
library(DESeq2)
library(Matrix)

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

proj_name <- "multiome_PP20251110"
groupBy <- "Sample_cell_type"

output_path <- file.path(proj_name, "DifferentialMotif")
if(!dir.exists(output_path)){
    dir.create(output_path)
}

proj <- loadArchRProject(proj_name)
motif_coverage_path <- file.path(proj_name, "MotifCoverage", groupBy)

#######################################################################
# Import and prepare motif coverage information per group
#######################################################################

# Import canonical motifs
possible_motif_files <- list.files(motif_coverage_path)
selected_file <- possible_motif_files[grep("Trp53_canonical", possible_motif_files)]
coverage_canonical_df <- load_coverage_information_asDataFrame(file.path(motif_coverage_path, selected_file))

# Import non-canonical motifs
selected_file <- possible_motif_files[grep("Trp53_noncanonical", possible_motif_files)]
coverage_noncanonical_df <- load_coverage_information_asDataFrame(file.path(motif_coverage_path, selected_file))

# Import ChIP-seq results from Kenzelmann-Broz et al.
dataset_path <- "../../support/GSE46240_ChIP-seq_peaks_p53_binding_sites_mm10.bed"
p53_chip <- read.table(dataset_path, header = F, stringsAsFactors = F, sep = "\t")
colnames(p53_chip) <- c("seqnames", "start", "end", "peak_name", "score", "strand")
p53_chip$start <- p53_chip$start -100
p53_chip$end <- p53_chip$end +100
p53_chip$strand <- "*"
p53_chip_ranges <- dataframe_to_GRanges(p53_chip)

# Compute overlaps with canonical and non-canonical motif sites
noncanonical_ranges <- GRanges(rownames(coverage_noncanonical_df))
overlaps <- findOverlaps(noncanonical_ranges, p53_chip_ranges)
coverage_noncanonical_df$p53_ChIP_200 <- seq(dim(coverage_noncanonical_df)[1]) %in% queryHits(overlaps)

canonical_ranges <- GRanges(rownames(coverage_canonical_df))
overlaps <- findOverlaps(canonical_ranges, p53_chip_ranges)
coverage_canonical_df$p53_ChIP_200 <- seq(dim(coverage_canonical_df)[1]) %in% queryHits(overlaps)

mean(coverage_canonical_df$p53_ChIP_200)
#[1] 0.04342826
mean(coverage_noncanonical_df$p53_ChIP_200)
#[1] 0.01047904

# Import PeakData and compute overlaps with peaks
peaks <- read.table(gzfile(file.path(proj_name, "PeakData", "PeakData.tsv.gz")), header = T, sep = "\t")
peak_ranges <- GRanges(peaks$region)
overlaps <- findOverlaps(canonical_ranges, peak_ranges, type = "within")
coverage_canonical_df$within_peaks <- seq(dim(coverage_canonical_df)[1]) %in% queryHits(overlaps)
overlaps <- findOverlaps(noncanonical_ranges, peak_ranges, type = "within")
coverage_noncanonical_df$within_peaks <- seq(dim(coverage_noncanonical_df)[1]) %in% queryHits(overlaps)

# Consolidate ranges: I take all canonical motifs, but only the non-canonical motifs that have two half-sites and overlap a p53 binding site
motifs_final <- rbind(coverage_canonical_df, coverage_noncanonical_df[coverage_noncanonical_df$p53_ChIP_200,])

#######################################################################
# Get group metadata
#######################################################################
library(stringr)
sample_names <- colnames(motifs_final[,1:(dim(motifs_final)[2]-1)])
matches <- str_match(sample_names, "(.*)_(.*)")
metadata <- data.frame(matches)
colnames(metadata) <- c("sample_name", "sample_id", "cell_type")
rownames(metadata) <- sample_names

########################
# Test all batches
########################

motifs_sub <- motifs_final[, sample_names]
cell_type_filter <- metadata$cell_type %in% c("progenitor", "gastric")
batch_filter <- metadata$sample_id %in% c("p489c_shRen_caer_48h_multiome_1", "p489c_shRen_caer_48h_multiome_2", "DACE657_mKate2_multiome")
valid_columns <- cell_type_filter & batch_filter
motifs_sub <- motifs_sub[, valid_columns]
metadata_sub <- metadata[valid_columns,]

dds <- DESeqDataSetFromMatrix(countData = motifs_sub, colData = metadata_sub, design= ~ cell_type)
keep <- rowSums(counts(dds) >= 1) >= 2
dds <- dds[keep,]
dds <- DESeq(dds)
res <- data.frame(results(dds, contrast = c("cell_type", "progenitor", "gastric")))
res$p53_ChIP_200 <- motifs_final$p53_ChIP_200[match(rownames(res), rownames(motifs_final))]
res$within_peaks <- motifs_final$within_peaks[match(rownames(res), rownames(motifs_final))]

output_file <- file.path(output_path, sprintf("DESEQ20251116_FEATURES_p53_motif_canonical_noncanonicalchip200_MODE_coverage_BATCH_all_QUERY_progenitor_REFERENCE_gastric.tsv.gz"))
write.table(res, gzfile(output_file), row.names=T, quote=F, sep = '\t')

########################
# Test only JR batch
########################

motifs_sub <- motifs_final[, sample_names]
cell_type_filter <- metadata$cell_type %in% c("progenitor", "gastric")
batch_filter <- metadata$sample_id %in% c("p489c_shRen_caer_48h_multiome_1", "p489c_shRen_caer_48h_multiome_2")
valid_columns <- cell_type_filter & batch_filter
motifs_sub <- motifs_sub[, valid_columns]
metadata_sub <- metadata[valid_columns,]

dds <- DESeqDataSetFromMatrix(countData = motifs_sub, colData = metadata_sub, design= ~ cell_type)
keep <- rowSums(counts(dds) >= 1) >= 2
dds <- dds[keep,]
dds <- DESeq(dds)
res <- data.frame(results(dds, contrast = c("cell_type", "progenitor", "gastric")))
res$p53_ChIP_200 <- motifs_final$p53_ChIP_200[match(rownames(res), rownames(motifs_final))]
res$within_peaks <- motifs_final$within_peaks[match(rownames(res), rownames(motifs_final))]

output_file <- file.path(output_path, sprintf("DESEQ20251116_FEATURES_p53_motif_canonical_noncanonicalchip200_MODE_coverage_BATCH_m3_QUERY_progenitor_REFERENCE_gastric.tsv.gz"))
write.table(res, gzfile(output_file), row.names=T, quote=F, sep = '\t')
