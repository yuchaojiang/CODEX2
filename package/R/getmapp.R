getmapp <- function(chr, ref) {
    if (chr == "X" | chr == "x" | chr == "chrX" | chr == "chrx") {
        chrtemp <- 23
    } else if (chr == "Y" | chr == "y" | chr == "chrY" | chr == "chry") {
        chrtemp <- 24
    } else {
        chrtemp <- as.numeric(mapSeqlevels(as.character(chr), "NCBI")[1])
    }
    if (length(chrtemp) == 0) 
        message("Chromosome cannot be found in NCBI Homo sapiens database!")
    mapp <- rep(1, length(ref))
    mappability <- mappability[[chrtemp]]
    mapp_ref <- mapp_ref[[chrtemp]]
    for (k in 1:length(ref)) {
        index <- countOverlaps(mapp_ref, ref[k])
        if (sum(index) > 0) {
            mapp[k] <- mean(mappability[as.logical(index)])
        }
    }
    mapp
}