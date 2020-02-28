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
#' @param input_GRanges A GRanges file with non-overlapping intervals that will be converted to a chain file. Required.
#' @param chrom_suffix The suffix to be appended to all the sestination chromosome names in the mapping "default is "mapped"
#' @param out_chain_name The name of a chain file to be written in the local directory. Default is "", is calls a tmp file creation. In this case, the file will be unlinked before the function returns.
#' @param verbose Output updates while the function is running. Default FALSE
#' @param chromosomes.lengths is sequinfo of the mapping object is not enough for the chromosome lengths, the additional info is provided here. Default is c(). The foemat is like the seqlengths() result for a GRanges.
#' @return a Chain object that maps all the chromosomes according to GRanges
# this ia code by Veronica Busa with some modification s by Alexander Favorov

#' @import plyranges
#' @import rtracklayer
GRangesMappingToChainFile<-function(input_GRanges,
                                    chrom_suffix = "",
                                    out_chain_name= "",
                                    chromosomes.length=c(),
																		verbose=FALSE
{
  #confirm GRanges doesn't have any overlapping intervals
  if(! identical(input_GRanges,reduce(input_GRanges))){
    stop("GRanges object includes overlapping intervals. This will cause errors when trying to use the chain object.")
  }
  #create seqinfo object for all chromosomes to get lengths
  seqinf<-as(seqinft[nchar(rownames(seqinft))<6,],"Seqinfo")
  # write chain file chromosome-by-chromosome
  if(sum(input_GRanges@seqnames %in% seqinf@seqnames) == 0){
    input_GRanges@seqnames<-paste0("chr", input_GRanges@seqnames) %>% Rle()
    if(sum(input_GRanges@seqnames %in% seqinf@seqnames) == 0){
      stop("GTF chromosomes do not resemble UCSC chromosome names. Suggested format: 'chr#' or just # for chromosome name, e.g. chr1 chr10 chrM")
    }
  }
  
  row<-1
  interval<-1
  chrom_chains<-matrix(NA, ncol=3, nrow=(length(input_GRanges)))
  
  # make a chain for each chromosome in the genome
  for(chr in seqinf@seqnames){
    if(verbose==TRUE){print(paste("Chromosome", chr, "starting"))}
    gtf_hold<-input_GRanges %>% filter(seqnames==chr)
    if(length(gtf_hold@ranges)==0){next} # in case of chromosomes without data
    length_chr<-sum(gtf_hold@ranges@width)
    first_line<-c(paste0("chain 42 ", chr, " ", seqinf@seqlengths[which(seqinf@seqnames==chr)], 
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
  tmp<-""
  if(out_chain_name == ""){
    out_chain_name<-tempfile(pattern = "", fileext = ".chain")
    tmp<-out_chain_name
  } else{if(verbose == TRUE){print("Saving chain object")}}
  writeLines(format_chrom_chains, con=out_chain_name)
  # have to write intermediate file to inport chain object
  chain<-import.chain(out_chain_name)
  unlink(tmp)
  return(chain)
}


#'@export
MapRangesToGenomicIntervals<-function(
	where.to.map, what.to.map,
	chromosomes.to.proceed=NA,
	unmapped.chromosome.warning=TRUE,
	unmapped.range.warning=FALSE,
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
	chromosome.names.what<-NA
	chromosome.names.what<-as.character(unique(seqnames(what.to.map)))
	#seqnames<-c()
	#seqleninfo<-c()
	#start=c()
	#end=c()
	for (chr in chromosome.names.what)
	{
		if (!is.na(chromosomes.to.proceed) && ! chr %in% chromosomes.to.proceed) next;
		if (! chr %in% chromosome.names.where )
		{
			if (unmapped.chromosome.warning) 
				warning(paste("The chromosome",chr,"has no space to be mapped to."))
			next;
		}
		whereranges<-sort(ranges(where.to.map)[[chr]])

		if (! isNormal(whereranges))
		{
			if (nonnormalised.mapping.warning) 
				warning(paste("The chromosome",chr,"is not normalised in where.to.map; I do it for you."))
			whereranges <- asNormalIRanges(whereranges)
		}
	
		whatranges<-sort(ranges(what.to.map)[[chr]])
		# 1.22 code
		#mapping.mat<-as.matrix(findOverlaps(whatranges,whereranges))
		
		#map.result<-apply(mapping.mat,1,
		#	function(hit){
		#		ind.what<-hit[1]
		#		ind.where<-hit[2]
		#		start<-max(1,start(whatranges)[ind.what]-start(whereranges)[ind.where]+1)
		#		end<-min(end(whatranges)[ind.what]-start(whereranges)[ind.where]+1,width(whereranges)[ind.where])
		#		seqnames<-paste0(chr,':',start(whereranges)[ind.where],'-',end(whereranges)[ind.where])
		#		c(seqnames=seqnames,end=end,start=start)
		#	}
		#)
	
		#seqnames<-c(seqnames,map.result['seqnames',])
		#start<-c(start,as.integer(map.result['start',]))
		#end<-c(end,as.integer(map.result['end',]))
	}
	if(unmapped.range.warning && ( length(space(what.to.map)) > length(result) ))
		warning("Some ranges remained unmapped.\n")
	if (length(seqnames)==0)
	{
		warning("Empty mapping result.\n");
		return (GRanges())
	}
	#prepare seqlength info
	seqlenames<-unique(seqnames)
	seqlengths<-sapply(as.vector(strsplit(seqlenames,':')), #list, we need
		function(splt){
			startend<-as.integer(strsplit(splt[2],'-')[[1]])
			startend[2]-startend[1]+1	
		}
	)
	names(seqlengths)<-seqlenames
	return
	(
		GRanges(
			seqnames=seqnames,
			ranges=IRanges(start=start,end=end),
			seqlengths=seqlengths)
	)
}


#'@export
WriteChainFromGRanges<-function(where.to.mao,chainfile) {
}
