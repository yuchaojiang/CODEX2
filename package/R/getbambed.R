getbambed <- function(bamdir, bedFile, sampname, projectname, chr) {
    if (!file.exists(bedFile))
        stop("Please check the bed file directory provided. File could not be 
            found!")
    exomtarg <- read.table(bedFile, sep = '\t')
    exomtarg <- exomtarg[exomtarg[,1] == chr,]
    ref <- IRanges(start = exomtarg[,2], end = exomtarg[,3])
    if (length(ref) == 0) 
        message('No exonic targets loaded from the bed file. Check chr style.')
    list(bamdir = bamdir, sampname = sampname, ref = ref, projectname = 
             projectname, chr = chr)
}