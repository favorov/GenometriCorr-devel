### R code from vignette source 'GenometriCorr.Rnw'

###################################################
### code chunk number 1: set_options
###################################################
options(width=50)


###################################################
### code chunk number 2: load_genometricorr
###################################################
library('GenometriCorr')


###################################################
### code chunk number 3: read_files
###################################################
library('rtracklayer')
#cpgis<-import("../extdata/UCSCcpgis_hg19.bed")
refseq<-as(import(system.file("extdata", "UCSCrefseqgenes_hg19.bed", package = "GenometriCorr")),'RangedData');
cpgis<-as(import(system.file("extdata", "UCSCcpgis_hg19.bed", package = "GenometriCorr")),'RangedData');
#refseq <- import("../extdata/UCSCrefseqgenes_hg19.bed")


###################################################
### code chunk number 4: chrom_length
###################################################
human.chrom.length<-c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,59373566,155270560)


###################################################
### code chunk number 5: chrom_names
###################################################
names(human.chrom.length)<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrY','chrX')


###################################################
### code chunk number 6: dev_off_init
###################################################
(if (length(dev.list())) dev.off())


###################################################
### code chunk number 7: chr1_picture
###################################################
VisualiseTwoIRanges(cpgis['chr1']$ranges,refseq['chr1']$ranges,nameA='CpG Islands',nameB='RefSeq Genes',chrom_length=human.chrom.length[['chr1']],title="CpGIslands and RefGenes on chr1 of Hg19",close.device=FALSE)


###################################################
### code chunk number 8: test_chromosome_123
###################################################

pn.area<-100
pn.dist<-100
pn.jacc<-100

cpgi_to_genes<-GenometriCorrelation(cpgis,refseq,chromosomes.length=human.chrom.length,chromosomes.to.proceed=c('chr1','chr2','chr3'),ecdf.area.permut.number=pn.area,mean.distance.permut.number=pn.dist,jaccard.measure.permut.number=pn.jacc,keep.distributions=TRUE,showProgressBar=FALSE)


###################################################
### code chunk number 9: report_chromosome_123
###################################################
print(cpgi_to_genes)


###################################################
### code chunk number 10: chromosome_1_graphical_report
###################################################
graphical.report(cpgi_to_genes,pdffile="CpGi_to_RefSeq_chr1_picture.pdf",show.chromosomes=c('chr1'),show.all=FALSE)


###################################################
### code chunk number 11: chromosome_1_visualize
###################################################
visualize(cpgi_to_genes,pdffile="CpGi_to_RefSeq_chr1_picture_vis.pdf",show.chromosomes=c('chr1'),show.all=FALSE)


###################################################
### code chunk number 12: read_config_file
###################################################
config<-new("GenometriCorrConfig",system.file("extdata", "template-config.ini", package = "GenometriCorr"))


###################################################
### code chunk number 13: print_config_file
###################################################
print(config)


###################################################
### code chunk number 14: change_config_file
###################################################
config$tests$ecdf.area.permut.number<-10
config$tests$mean.distance.permut.number<-10
config$tests$jaccard.measure.permut.number<-10
config$chromosomes <- "chr18"


###################################################
### code chunk number 15: read_files
###################################################
library('rtracklayer')
histones<-import(system.file("extdata", "chr18.H3K4me3.bed", package = "GenometriCorr"));
expr_genes<-import(system.file("extdata", "chr18.mRNAseq.bed", package = "GenometriCorr"));


###################################################
### code chunk number 16: run_config_file
###################################################
conf_res<-run.config(config,query=histones,reference=expr_genes)


###################################################
### code chunk number 17: simple_config_file
###################################################
quickconfig<-new("GenometriCorrConfig",system.file("extdata", "quick-config.ini", package = "GenometriCorr"))
print(quickconfig)


###################################################
### code chunk number 18: mapping_example_random
###################################################
population<-1000

chromo.length<-c(3000000)

names(chromo.length)<-c('the_chromosome')

rquery<-RangedData(ranges=IRanges(start=runif(population,1000000,2000000-10),width=c(10)),space='the_chromosome')

rref<-RangedData(ranges=IRanges(start=runif(population,1000000,2000000-10),width=c(10)), space='the_chromosome')

#create two features, they are randomly scattered in 1 000 000...2 000 000

unmapped_result<-GenometriCorrelation(rquery,rref,chromosomes.length=chromo.length,ecdf.area.permut.number=pn.area,mean.distance.permut.number=pn.dist,jaccard.measure.permut.number=pn.jacc,keep.distributions=FALSE,showProgressBar=FALSE)

#correlate them on the whole cromosome: 1...3 000 000

map_space<-RangedData(ranges=IRanges(start=c(1000000),end=c(2000000)),space='the_chromosome')

mapped_rquery<-MapRangesToGenomicIntervals(what.to.map=rquery,where.to.map=map_space)

mapped_rref<-MapRangesToGenomicIntervals(what.to.map=rref,where.to.map=map_space)

#map them into 1 000 000...2 000 000

mapped_result<-GenometriCorrelation(mapped_rquery,mapped_rref,ecdf.area.permut.number=pn.area,mean.distance.permut.number=pn.dist,jaccard.measure.permut.number=pn.jacc,keep.distributions=FALSE,showProgressBar=FALSE)

#then, correlate again

cat('Unmapped result:\n')
print(unmapped_result)

cat('Mapped result:\n')
print(mapped_result)



###################################################
### code chunk number 19: simple_config_file
###################################################
mapconfig<-new("GenometriCorrConfig",system.file("extdata", "mapping-config.ini", package = "GenometriCorr"))
print(mapconfig)


