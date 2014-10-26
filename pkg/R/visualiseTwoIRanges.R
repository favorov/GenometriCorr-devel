# GenometriCorrelation project evaluating two markups genomwide independence. 
# (c) 2010-2011 Alexander Favorov, Leslie Cope, Yulia Medvedeva, 
#              Loris Mularoni, Vsevolod Makeev, Sarah Wheelan.
# VisualiseTwoIRanges is a function that visualise a pair of IRanges on a chromosome in two colors
# $Id: visualiseTwoIRanges.R 1717 2012-04-18 19:29:15Z favorov $

# HJ -- pdf and close.device arguments

VisualiseTwoIRanges<-function(irA, irB, start=1, end=NA, nameA='RangesA', nameB='RangesB',
	chrom_length=NA, title=NA, pdf=NULL, close.device=NULL)
{
	pixelizeIRanges<-function(iranges,start,end,max.pixels=10000)
	{
		len=end-start+1
		if (len>max.pixels)
		{
			#less than max_pixels*bin are to cover len
			bin <- len %/% max.pixels	
			if ((len %% max.pixels) > 0) bin<-bin+1
			test.ranges<-IRanges(start=seq.int(from=1,by=bin,length.out=max.pixels),width=bin)
			overrle<-findOverlaps(test.ranges,iranges)
			masker<-tapply(width(iranges[subjectHits(overrle)]),queryHits(overrle),sum)
			mask<-rep(0,max.pixels)
			mask[as.integer(names(masker))]<-as.integer(masker)
		} else
		{
			mask<-as.vector(coverage(iranges,shift=1-start,width=len))
			bin<-1
		}
		return(mask)
	}

	#require(grDevices)

	if (!inherits(irA,"IRanges"))
		stop("The first parameter for TwoIRangesIndependence is to be IRange.")
	
	if (!inherits(irB,"IRanges"))
		stop("The second parameter for TwoIRangesIndependence is to be IRange.")

	if (is.na(chrom_length))
	{
		if (is.na(end)) 
		{
			chrom_length<-chromosomes.length.eval(irA,irB)
			end<-chrom_length
			warning("Chromosome length and the right margin of picture are evaluated rather than pre-given")
		}
		else
			chrom_length<-chromosomes.length.eval(irA,irB)
	}

	if (is.na(end))
		end<-chrom_length

	if (end<=start)
		stop("End is less or equal than start.")

	len<-end-start+1
	# pdf
	if (!is.null(pdf)) {
		if (length(grep("\\.pdf$", pdf)) == 0)
			pdf <- paste(pdf, ".pdf", sep="")
		pdf(pdf)
		if (is.null(close.device)) close.device=TRUE
	} else
	{
		if (is.null(close.device)) close.device=FALSE
	}

	par(yaxt='n')
	plot(c(start,end), c(-1, 1), type = "n", xlab="", ylab="")

	#we are ready to think about what to plot.. \
	#let's get the length af the raster we are to prepare

	pixels<-as.integer(dev.size('px')[1]*(end-start+1)/xinch(dev.size('in')[1]))	

	maskA<-pixelizeIRanges(irA,start,end,max.pixels=pixels)
	
	maskB<-pixelizeIRanges(irB,start,end,max.pixels=pixels)
		
	#to debug	
	#maskA=sample(c(0,100,10000),pixels,replace=T,c(8/10,1/10,1/10))
	#maskB=sample(c(0,100,10000),pixels,replace=T,c(8/10,1/10,1/10))

	maskA <- 1-(maskA/max(maskA))
	maskB <- 1-(maskB/max(maskB))
	#maskA <- exp(-maskA) #to be better on view, now 0 is white and everything > 0 is a bit or more red 
	#maskB <- exp(-maskB) #to be better on view, now 0 is white and everything > 0 is a bit or more blue
 
	maxA<-max(maskA)
	maxB<-max(maskB)
	
	img_len<-length(maskA)
	if (img_len != length(maskB))
	{
		stop("Image lengthes are different, something went wrong.");
	}
	if (maxA==0)  #not to divide by 0 in as.raster
		image_red<-
			as.raster(array(c(rep(1,img_len),rep(1,img_len),rep(1,img_len)),c(1,img_len,3)),max=1)
	else
		image_red<-
			as.raster(array(c(rep(maxA,img_len),maskA,maskA),c(1,img_len,3)),max=maxA)
	if (maxB==0)  #not to divide by 0 in as.raster
		image_blue<-
			as.raster(array(c(rep(1,img_len),rep(1,img_len),rep(1,img_len)),c(1,img_len,3)),max=1)
	else
		image_blue<-
			as.raster(array(c(maskB,maskB,rep(maxB,img_len)),c(1,img_len,3)),max=maxB)

	rasterImage(image_red, start, .05,end,.75)
	rasterImage(image_blue, start, -.75,end,-.05)
	text(c(start+len/2,start+len/2),c(.9,-.9),c(nameA,nameB))
	if (!is.na(title))
	{
		title(main=title)
	}

	if (close.device)
		invisible(dev.off())
}
