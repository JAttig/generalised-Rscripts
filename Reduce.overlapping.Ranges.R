### This function takes a GRange as input, and merges overlapping Ranges
# it does return a GRange object with same mcols information
# but asssigns to a merged range the annotation from the _first_ range only



rm.overlap.ranges <- function(gr){
  reduced <- reduce(gr, with.revmap = TRUE)
  mapping <- mcols(reduced)$revmap
  subset <- unlist(lapply(mapping, "[[", 1))
  return.gr <-  gr[subset]
  return(return.gr)
}


### Another option is to remove "real" duplicates i.e. IRanges with identical start and stop position
# The function returns a GRange object with same mcols information
# but asssigns to a merged range the annotation from the _first_ range only

rm.identical.dupl.ranges <- function(gr){
  Range.df <- data.frame(starts = start(gr), ends = end(gr))
  duples <- duplicated(Range.df)
  return.gr <- gr[!duples]
  return(return.gr)
}

rm.identical.dupl.ranges.with.identical.mcols <- function(gr){
  Range.df <- as.data.frame(gr)
  duples <- duplicated(Range.df)
  return.gr <- gr[!duples]
  return(return.gr)
}


### Yet another option is to disjoin the GRange to keep unique fragments - but this also deletes the mcols()
find.unique.ranges <- function(gr){
  return.gr <- disjoin(gr)
  return(return.gr)
}
