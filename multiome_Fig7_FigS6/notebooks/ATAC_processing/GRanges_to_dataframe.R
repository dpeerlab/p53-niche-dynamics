GRanges_to_dataframe <- function(current_ranges){
    seqnames <- as.character(seqnames(current_ranges))
    start <- start(current_ranges)
    end <- end(current_ranges)
    strand <- as.character(strand(current_ranges))
    metadata <- mcols(current_ranges)
    result <- data.frame(seqnames = seqnames, start = start, end = end, strand = strand)
    for(mycol in colnames(metadata)){
        result[[mycol]] <- ifelse(isfactor(metadata[[mycol]]), as.character(metadata[[mycol]]), metadata[[mycol]])
    }
    return(result)
}
