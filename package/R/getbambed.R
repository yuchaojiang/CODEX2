getbambed = function (bamdir, bedFile, sampname, projectname) 
{
  if (!file.exists(bedFile)) 
    stop("Please check the bed file directory provided. File could not be \n            found!")
  exomtarg <- read.table(bedFile, sep = "\t")
  ref <- GRanges(seqnames=exomtarg[,1],ranges=IRanges(start=exomtarg[,2], end=exomtarg[,3]))
  #ref <- sort(ref)
  list(bamdir = bamdir, sampname = sampname, ref = ref, projectname = projectname)
}
