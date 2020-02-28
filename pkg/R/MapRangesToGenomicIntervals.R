# GenometriCorrelation project evaluating two interval markups genomewide independence. 
# (c) 2010-2020 Alexander Favorov, Loris Mularoni, 
#               Yulia Medvedeva, Harris A. Jaffee, 
#               Ekaterina V. Zhuravleva, Veronica Busa,
#               Leslie M. Cope, Andrey A. Mironov, 
#               Vsevolod J. Makeev, Sarah J. Wheelan.

#' GRangesMappingToChainViaFile
#' 
#' GRangesMappingToChainViaFile creates a \code{Chain} object based on intervals from a \code{GRanges} object. The mapping by this \code{Chain} will collapse chromosomes of the annotation into pseudogenome that are combined from tiled intervals of the \code{GRanges} object. 
#' 
#' 
#' @param ranges_to_map_to A GRanges file with non-overlapping intervals that will be converted to a chain file. Required.
#' @param chrom_suffix The suffix to be appended to all the sestination chromosome names in the mapping "default is "mapped"
#' @param out_chain_name The name of a chain file to be written in the local directory. Default is "", is calls a tmp file creation. In this case, the file will be unlinked before the function returns.
#' @param verbose Output updates while the function is running. Default FALSE
#' @param chromosomes.length is sequinfo of the mapping object is not enough for the chromosome lengths, the additional info is provided here. Default is c(). The foemat is like the seqlengths() result for a GRanges.
#' @return a Chain object that maps all the chromosomes according to GRanges
# this ia code by Veronica Busa with some modification s by Alexander Favorov

#' @import plyranges
#' @import rtracklayer
GRangesMappingToChainViaFile<-function(ranges_to_map_to,
                                    chrom_suffix = "_mapped",
                                    out_chain_name= "",
                                    chromosomes.length=c(),
																		verbose=FALSE
{
  #confirm GRanges doesn't have any overlapping intervals
  ranges_to_map_to<-reduce(ranges_to_map_to)
  #chromosomes to get lengths
	chrom_length<-seqlengths(ranges_to_map_to)
	for (name in as.character(ranges_to_map_to@seqnames)) {
		if (is.na(seqlengths(ranges_to_map_to)[name])){ 
			if(!is.na(chromosomes.length[name])) {
				seqlengths(ranges_to_map_to)[name] <- chromosomes.length[name]
			} else {
				seqlengths(ranges_to_map_to)[name] <- max(end(ranges_to_map_to %>% filter(seqnames==name)))
			}
		}
	}
	
  
  row<-1
  interval<-1
  chrom_chains<-matrix(NA, ncol=3, nrow=(length(ranges_to_map_to)))
  
  # make a chain for each chromosome in the genome
  for(chr in as.character(ranges_to_map_to@seqnames)){
    if(verbose==TRUE){print(paste("Chromosome", chr, "starting"))}
    gtf_hold<-ranges_to_map_to %>% filter(seqnames==chr)
    if(length(gtf_hold@ranges)==0){next} # in case of chromosomes without data
    length_chr<-sum(gtf_hold@ranges@width)
    first_line<-c(paste0("chain 42 ", chr, " ", seqlengths(ranges_to_map_to)[name], 
                         " * ", gtf_hold@ranges@start[1], " ", 
                         gtf_hold@ranges@start[length(gtf_hold@ranges@start)]+gtf_hold@ranges@width[length(gtf_hold@ranges@width)]-1,
                         " ", chr, chrom_suffix, " ", length_chr, " * 1 ", 
                         length_chr," ", interval), "", "")
    interval<-interval + 1
    chrom<-matrix(NA, nrow=length(gtf_hold@ranges@start), ncol=3)
    chrom[,1]<-gtf_hold@ranges@width
    if(length(gtf_hold@ranges@width)>1){ #for cases where there is only 1 interval in a chromosome
      chrom[,2]<-c(gtf_hold@ranges@start[2:length(gtf_hold@ranges@start)]-
                  (gtf_hold@ranges@width[1:(length(gtf_hold@ranges@width)-1)]+
                  gtf_hold@ranges@start[1:(length(gtf_hold@ranges@width)-1)])-1, "")
      chrom[,3]<-0
      chrom[nrow(chrom),3]<-""
      } 
    chrom_chains[row,]<-first_line
    chrom_chains[(row+1):(row+nrow(chrom)),]<-chrom
    chrom_chains[row+nrow(chrom)+1,]<-c("", "", "")
    row<-row+nrow(chrom)+2
  }
  chrom_chains<-na.omit(chrom_chains)
  if(verbose == TRUE){print("Creating chain object")}
  format_chrom_chains<-apply(chrom_chains, 1, paste, collapse = "\t")
  format_chrom_chains<-gsub("\t\t", "", format_chrom_chains)
  if(out_chain_name == ""){
    out_chain_name<-tempfile(pattern = "", fileext = ".chain")
  } else{if(verbose == TRUE){print("Saving chain object")}}
  writeLines(format_chrom_chains, con=out_chain_name)
  # have to write intermediate file to inport chain object
  chain<-import.chain(out_chain_name)
  unlink(out_chain_name)
  return(chain)
}

#' @param chrom.suffix The suffix to be appended to all the sestination chromosome names in the mapping; default is "mapped"
#' @export
MapRangesToGenomicIntervals<-function(
	where.to.map, what.to.map,
	chrom.suffix,
	chromosomes.to.proceed=NA,
	unmapped.chromosome.warning=TRUE,
	nonnormalised.mapping.warning=TRUE
)
#subgenome and are suppose to be GRanges
#in the second case, we test the lenthg equivalence
{
	is.where.gr<-inherits(where.to.map,"GRanges")
	if (!is.where.gr)
		stop("where.to.map is not GRanges. It's all lost!\n")	
	chromosome.names.where<-NA
	chromosome.names.where<-as.character(unique(seqnames(where.to.map)))

	is.what.gr<-inherits(what.to.map,"GRanges")
	if (is.what.gr)
		stop("what.to.map is not GRanges. It's all lost!\n")	
	
	#we are here, they are granges
	if (!is.na(chromosomes.to.proceed)) {
		where.to.map <- where.to.map %>% filter(seqnames %in% chromosomes.to.proceed)
		what.to.map <- what.to.map %>% filter(seqnames %in% chromosomes.to.proceed)
	}
	
	if(! identical(where.to.map,reduce(where.to.map))){
    if (nonnormalised.mapping.warning) {warning("GRanges object to map to includes overlapping intervals. Normalised.")}
		where.to.map<-reduce(where.to.map)
  }

	if(unmapped.chromosome.warning) {
		unmapped_chroms<-setdiff(what.to.map@seqnames,where.to.map@seqnqmes)
		if(length(unmapped_chroms)>0) warning(paste0("Some chromosomes, e.g. ",unmapped_chroms[1]," has no mapping,"))
	}
	
	return
	(
		GRanges(
			seqnames=seqnames,
			ranges=IRanges(start=start,end=end),
			seqlengths=seqlengths)
	)
}

