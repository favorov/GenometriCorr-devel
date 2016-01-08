#this script contain function that counts mean value for coverage on set of intervals
#input: input_GRanges object from file and intervals set as GRanges object,
#chromosomes names and length is required

library('GenomicRanges')
library('Rsamtools')
library('IRanges')
library('rtracklayer')




GenometriCorrelationWilcoxonTest <-funcion(
  query, reference, 
  chr.to.proceed=c(), 
  chromosomes.length =с(),
  map.to=NA)

{ #input parameteres description

  #various format checks for input data #TODO
  
  result.lists<-.aggregateSpaces(
    query=query,
    reference=reference,
    chr.to.proceed = chr.to.proceed, 
    chromosomes.length = chromosomes.length, 
    map.to)

  result.estimated<-.estimateTest(result.lists)

return(result.estimated)  
}

.aggregateSpaces <- function(query, reference, chr.to.proceed, chromosomes.length, map.to){
  #first we made space.correspondence, which is that we mapped to  
  if (suppressWarnings(is.na(map.to))){
    #than we consider chromnosomes intervals
    space.description <-IRanges(rep(1, length(chr.to.proceed)), chromosomes.length)
    names(space.description) <- names(chr.to.proceed)
    
    reference.space.correspondence = sapply(split(reference, seqnames(reference)),ranges)
    space.correspondence = sapply(split(query, seqnames(query)),ranges)
    
  }  else {
    #than we consider setted intervals
    spaces.names = paste(as.vector(seqnames(map.to)), start(intervals.set), end(intervals.set), sep='_')
    space.description <-ranges(map.to)
    names(space.description) <- spaces.names
    
    #space.correspondence = suppressWarnings(lapply(map.to, function(x) intersect(x, query)))
    space.correspondence = suppressWarnings(lapply(map.to, function(x) ranges(findOverlaps(ranges(query), ranges(x)),  ranges(query), ranges(x))))
                     
    names(space.correspondence) = spaces.names
    #space.correspondence = sapply(space.correspondence, ranges)
    
    #reference.space.correspondence = suppressWarnings(lapply(map.to, function(x) intersect(x, reference)))
    reference.space.correspondence = suppressWarnings(lapply(map.to, function(x) ranges(findOverlaps(ranges(reference), ranges(x)),  ranges(reference), ranges(x))))
    
    names(reference.space.correspondence) = spaces.names
    #reference.space.correspondence =sapply(reference.space.correspondence, ranges)

  }
  
  #here we get our space.correspondence and reference.space.correspondence,
  #that is object like this:
  #$chr21_40000000_48119976
  #IRanges of length 8
  #start      end   width
  #[1] 40000000 40285945  285946
  #[2] 40285951 42955560 2669610
  

  #find overlaps in map.to or chr rergions
  in.intervals.overlapping <- sapply(names(space.correspondence), function(x) findOverlaps(space.correspondence[[x]], reference.space.correspondence[[x]]))
  
  #take intersections for overlaps
  #in.intervals.overlapping.extracted.intersection <- ranges(in.intervals.overlapping, ranges(a.g), ranges(test[i]))
  in.intervals.overlapping.extracted.intersection <- sapply(names(in.intervals.overlapping), function(x) ranges(in.intervals.overlapping[[x]], space.correspondence[[x]], reference.space.correspondence[[x]]))
  
  #todo   !!!range function
  
  #now we count coverage on each intersection, it will be a vector, we need only mean value
  mean.on.intersection.byintervals<- sapply(names(in.intervals.overlapping.extracted.intersection), function(x) mean(as.vector(coverage(in.intervals.overlapping.extracted.intersection[[x]]))))
  #mean.on.intersection<- c(mean.on.intersection, mean(as.vector(coverage(ranges(intersection.for.interval)))))
  
  #now we need calculate intervals that do not intersect and calculate mean coverage on it
  #also contain 1 width intervals
  reference.space.correspondence.inverse <- sapply(names(reference.space.correspondence), function(x) c(gaps(reference.space.correspondence[[x]]), IRanges(c(start(space.description[x]), end(c(range(reference.space.correspondence[[x]])))), c(start(c(range(reference.space.correspondence[[x]]))), end(space.description[x])))) )
  
  #delete start= end from inverse set !
  
  #calculate means on inverse interval set
  in.inverse.overlapping <- sapply(names(space.correspondence), function(x) findOverlaps(space.correspondence[[x]], reference.space.correspondence.inverse[[x]]))
  in.inverse.overlapping.extracted.intersection <- sapply(names(in.inverse.overlapping), function(x) ranges(in.inverse.overlapping[[x]], space.correspondence[[x]], reference.space.correspondence.inverse[[x]]))
  mean.on.intersection.inverse<- sapply(names(in.inverse.overlapping.extracted.intersection), function(x) mean(as.vector(coverage(in.inverse.overlapping.extracted.intersection[[x]]))))
  
  #TODO
  
}

#.estimateTest(result.lists){
  
#}


"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}


##################
#input sample 


##################

#parameters for human chromosomes
human.chrom.length <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 59373566,155270560)
names(human.chrom.length) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrY", "chrX")


input.file= c('bigWigExample.bw')
##BigWig data import query
data_gr<-import.bw(con = BigWigFile(input.file),format="bigWig")
countIntervals.between<-c('True')

#reference intervals
test<-GRanges(seqnames = c("chr21", "chr21" , 'chr21'), ranges = IRanges(c(9411201,48119776, 9411101), c(9412600, 48119976, 9411151)))

#for maipping
intervals.set = GRanges(seqnames = c("chr21", "chr21", "chr21" ), ranges = IRanges(c(0, 9000000,40000000), c(100, 9411600, 48119976)))



#-----parameteres
reference = test
query = data_gr
chromosomes.length = human.chrom.length 
chr.to.proceed = c('chr21')
map.to = intervals.set 

#map.to=NA

#---------

list <- structure(NA,class="result")
list[mean.on.intersection, mean.on.Notintersection]<- GenometriCorrelationWilcoxonTest(data_gr,test, chr.to.proceed=c('chr21'), chromosomes.length = human.chrom.length, map.to = intervals.set) 

#list[mean.on.intersection, mean.on.Notintersection]<- countMeansforIntervalsSet(data_gr,test)



