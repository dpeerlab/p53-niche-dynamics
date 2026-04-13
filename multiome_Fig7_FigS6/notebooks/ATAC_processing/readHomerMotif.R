library(TFBSTools)
homerToPWMatrixList <- function(filename, n = 1000L, type="logratio") {
    stopifnot(file.exists(filename))

    # parse HOMER motif file
    tmp <- readLines(filename)
    g <- grep(">", tmp)
    tmp[g] <- sub("^>","",tmp[g])
    fields <- strsplit(tmp[g], "\t", fixed = TRUE)
    L <- lapply(seq_along(g), function(i) {
        cons <- fields[[i]][1]
        nm <- fields[[i]][2]
        log2cut <- round(as.numeric(fields[[i]][3]) / log(2), 2)
        s <- g[i] + 1L
        e <- if (i == length(g)) length(tmp) else g[i + 1L] - 1L
        m <- round(n * do.call(cbind, 
                               lapply(strsplit(tmp[s:e], "\t", fixed = TRUE), 
                                      as.numeric)), 0)
        rownames(m) <- c("A", "C", "G", "T")
	freqMatrix = m/n
	freqMatrix = t(t(freqMatrix)/apply(freqMatrix, 2, sum)) / 0.25
        if(type=="logratio")
		profileMatrix = log(freqMatrix)
	else
		profileMatrix = freqMatrix
        TFBSTools::PWMatrix(
          ID = nm, name = nm, profileMatrix = profileMatrix,
          tags = list(log2cut = log2cut, 
                      comment = "imported from HOMER motif file")
        )
		
    })

    result = do.call(PWMatrixList, c(L, list(use.names = TRUE)))
    names(result) = sapply(result, name)
    return(result)
}
