qc <- function(Y, sampname, chr, ref, mapp, gc, cov_thresh, length_thresh, 
               mapp_thresh, gc_thresh) {
    # sample wise QC
    sampfilter <- apply(Y, 2, sum) > 2000
    Y_qc <- Y[, sampfilter]
    sampname_qc <- sampname[sampfilter]
    # exon wise QC
    binfiltera <- (apply(Y_qc, 1, median) > cov_thresh[1]) & (apply(Y_qc, 
        1, median) < cov_thresh[2])
    message("Excluded ", sum(1 - binfiltera), " exons due to extreme coverage.")
    binfilterb <- ((end(ref) - start(ref)) > length_thresh[1]) & ((end(ref) - 
        start(ref)) < length_thresh[2])
    message("Excluded ", sum(1 - binfilterb), 
        " exons due to extreme exonic length.")
    binfilterc <- (mapp >= mapp_thresh)
    message("Excluded ", sum(1 - binfilterc), 
        " exons due to extreme mappability.")
    binfilterd <- (gc >= gc_thresh[1] & gc <= gc_thresh[2])
    message("Excluded ", sum(1 - binfilterd), 
        " exons due to extreme GC content.")
    binfilter <- binfiltera & binfilterb & binfilterc & binfilterd
    message("After taking union, excluded ", sum(1 - binfilter), " out of ", 
        length(binfilter), " exons in QC.")
    qcmat <- cbind(rep(chr, length(ref)), start(ref), end(ref), binfilter, 
        apply(Y_qc, 1, median), binfiltera, (end(ref) - start(ref) + 1)/1000, 
        binfilterb, mapp, binfilterc, round(gc, 2), binfilterd)
    colnames(qcmat) <- c("chr", "start_bp", "end_bp", "pass", "median_depth", 
        "pass_depth", "length_kb", "pass_length", "mapp", "pass_mapp", 
        "gc", "pass_gc")
    Y_qc <- Y_qc[binfilter, ]
    gc_qc <- gc[binfilter]
    ref_qc <- ref[binfilter]
    mapp_qc <- mapp[binfilter]
    list(Y_qc = Y_qc, sampname_qc = sampname_qc, gc_qc = gc_qc,
         mapp_qc = mapp_qc, ref_qc = ref_qc, qcmat = qcmat)
}