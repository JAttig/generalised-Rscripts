#' Generate metadata profile
#' 
#' Create the metaprofile data for a set of GRanges using a coverage object.
#' @param PcovList A RleList object containing the information on the coverage for two or more datasets - positive strand. 
#' @param NcovList A RleList object containing the information on the coverage for two or more datasets - negative strand.
#' @param query A GRanges object with positions of interest. Must have a 'gene_id' column!
#' @param flanks The width of the GRanges in the query (must be the same for all GRanges). 
#' @param breaks A numeric vector of two or more unique cut points.
#' @details The coverage RleList object can be generated from any kind of NGS data. The function extracts from RleList the coverage across each entry in the GRanges object. 
#'    Since the output of this operation can be quite heavy, it subsets the regions into bins and reports the average coverage per bin.
#'    Query should not contain any out-of-bound ranges; and hence should have seqinfo (which will truncate out of seqinfo ragnes)
#' @return A List of data tables (one for each entry in the covList) reporting the binned coverage for each entry in the query.
#'
#' @import GenomicFeatures
#' @import data.table
#' @export
#' @examples createProfile(RNAseq_Coverage_List, RNAseq_Coverage_List, myGRanges, upstream=100, downstream=500, flanks=2000, min_count=10)

### function are changed to adapt to my R set-up

# Function to generate the metadata for the profiles
createProfile = function(PcovList, NcovList, query, sizeFac=NULL, upstream=NULL, downstream=NULL, 
                         antisense=FALSE, bin_width=20, min_count=10, keys=NULL){
   
   source("/Users/attigj/Documents/Jan.Crick/computational.analysis/R.sourcefunctions/sanitycheck.r")
   
   # Checks
   if( !sanitycheck(names(PcovList), names(NcovList)) ){
      print("Please make sure names in the positive and negative coverage objects are matching.")
      return(0)
   }
   if( !is.null(sizeFac) & (!sanitycheck(names(PcovList), names(sizeFac)) | !sanitycheck(names(NcovList), names(sizeFac))) ){
      print(paste0("The names in the RleLists (", length(names(PcovList)),";", length(names(NcovList)),") do not match those of the size factors (", length(names(sizeFac)),")."))
      return(0)
   }
   
   if( is.null(upstream) | is.null(downstream) ){
      print("Please specify the upstream and downstream widths from the reference position.")
      return(0)
   }
   
   if( (upstream+downstream)!=width(query[1]) ){
      print(paste0("The sum of upstream and downstream (", upstream+downstream,") must be equal to the GRanges width (", width(query[1]),")."))
      return(0)
   }
   
   #delete any row names in the GRange; these will conflict if they are not unique and are useless. .
   names(query) <- NULL
   
   # Create bins and labels   
   width = upstream+downstream
   breaks = seq(0, width, bin_width)
   lbls = seq(-upstream, downstream, bin_width)
   
   # Correct the orientation when countin on the antisense of the feature
   if(antisense)
      lbls = rev(lbls)

   lbls = data.table(bin = seq_along(lbls), labels = lbls, key="bin")
   
   # Create a data table object to extract the metadata information
   query.dt = as.data.table(as.data.frame(copy(query)))
   if(!is.null(keys)){
      tab = query.dt[, list(Tot = .N), by=keys]
      setkeyv(tab, keys)
   }
   setkey(query.dt, gene_id)
   
   if(interactive()) message("generated data.table and keyed it")
   
   # Use for progress tracking
   t = length(PcovList)
   
   # Main function
   myList = lapply(seq_along(PcovList),
                   function(i){
                      # If in interactive mode, report the progress status
                      if(interactive())
                         print(paste(paste0(i,"/",t), names(PcovList)[[i]]))
                      # Positive strand coverage
                      cov = PcovList[[i]]
                      sel = query[strand(query) == "+"]
                      seqlevels(sel, force=TRUE) <- seqlevels(cov)
                      
                      # Extract coverage values
                      pl = cov[sel ]
                      names(pl) = sel$gene_id
                      # Remove low count features
                      pl = pl[sum(pl>0)>=min_count,]
                      pl = rbindlist(lapply(pl, as.data.table), idcol = "gene_id")
                      
                      if(interactive()) message("calculated positive coverage per gene")
                      
                      # Add positional information
                      pl[, `:=`(Position, as.integer(1:width))]
                      # Negative strand coverage
                      cov = NcovList[[i]]
                      sel = query[strand(query) == "-"]
                      seqlevels(sel, force=TRUE) <- seqlevels(cov)
                      
                      # Extract coverage values
                      nl = cov[sel]
                      names(nl) = sel$gene_id
                      # Remove low count features
                      nl = nl[sum(nl>0)>=min_count,]
                      nl = rbindlist(lapply(nl, as.data.table), idcol = "gene_id")
                      # Add positional information
                      nl[, `:=`(Position, as.integer(width:1))]
                      # Merge strand coverages
                      l = rbindlist(list(pl, nl))
                      # Calculate intervals (binning) and average the values
                      l[, `:=`(bin, findInterval(Position, breaks, all.inside=TRUE))]
                      l = l[, list(value = mean(value)), by = c("gene_id", "bin")]
                      setnames(l, "bin", "Position")
                      setkey(l, gene_id)
                      # Add the metadate information
                      query.dt[, gene_id := as.character(gene_id) ]
                      setkey(query.dt, gene_id)
                      l = l[query.dt, c("gene_id", "value", "Position", keys), with=FALSE, nomatch = 0]
                      
                      ### normalise data
                      setkeyv(l, c("gene_id", keys))
                      tab2 = l[, list(Num = length(unique(gene_id))), by=keys]
                      setkeyv(tab2, keys)
                      
                      # Normalise the values by the maximum peak height and area
                      if(!is.null(sizeFac)){
                         l = l[, list(Position, value, value.sf=value/sizeFac[i], value.max=value/max(value), value.sum=value/sum(value)), by=key(l)]
                         # Set -Inf, Inf and NaN values to 0
                         l[!is.finite(value.max), value.max := 0]
                         l[!is.finite(value.sum), value.sum := 0]
                         # Collapse positions
                         setkeyv(l, c("Position", keys))
                         l = l[, list(value=sum(value), value.sf=sum(value.sf), value.max=sum(value.max), value.sum=sum(value.sum)), by=key(l)]
                         # Add the total number of entries for the specified 'keys'
                         if(!is.null(keys)){
                            setkeyv(l, keys)
                            l = l[tab2, ][tab, ]
                         }else{
                            l[tab2, ][, Tot := nrow(query.dt)]
                         }
                         # Replace bin labels with nucleotide positions
                         setkey(l, Position)
                         l = l[lbls, nomatch=0][, Position := labels][, labels := NULL]
                      }else{
                         l = l[, list(Position, value, value.max=value/max(value), value.sum=value/sum(value)), by=key(l)]
                         # Set -Inf, Inf and NaN values to 0
                         l[!is.finite(value.max), value.max := 0]
                         l[!is.finite(value.sum), value.sum := 0]
                         # Collapse positions
                         setkeyv(l, c("Position", keys))
                         l = l[, list(value=sum(value), value.max=sum(value.max), value.sum=sum(value.sum)), by=key(l)]
                         # Add the total number of entries for the specified 'keys'
                         if(!is.null(keys)){
                            setkeyv(l, keys)
                            l = l[tab2, ][tab, ]
                         }else{
                            l[tab2, ][, Tot := nrow(query.dt)]
                         }
                         # Replace bin labels with nucleotide positions
                         setkey(l, Position)
                         l = l[lbls, nomatch=0][, Position := labels][, labels := NULL]
                      }
                      l
                   })
   names(myList) = names(PcovList)
   return(myList)
}
