# GenometriCorrelation project evaluating two interval markups genomewide independence. 
# (c) 2010-2014 Alexander Favorov, Loris Mularoni, Yulia Medvedeva, 
#               Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, 
#               Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
# .aggregateSpaces is finction that arregate spacein GenometriCorrResult
# it has two parameters: space.correspondence and data. 
#space.correspondence is a character array of (N,2) and data 
#

#This is list of fields copied from docs as of Mar 21 2015
#the line starting with #* is what we do in this function
#the fields with ### line after description are the fields that are computed by 
#.estimateTests after the aggregation
#
#query.population	
#Query points used in the comparisons.
#*sum over spaces
#
#reference.population	
#Reference points used in the comparisons.
#*sum over spaces
#
#relative.distances.ks.p.value	
#p-value for local independence obtained by the Kolmogorov-Smirnov test for relative distances.
###
#
#relative.distances.ecdf.deviation.area.p.value	
#p-value for local independence obtained by the permutation test for relative distances.
###
#
#relative.distances.ecdf.area.correlation	
#Has the same sign with the relative distance-based local correlation.
###
#
#projection.test.p.value	
#p-value for chromosome-scale independence obtained by the projection test.
####
#
#projection.test.lower.tail	
#If TRUE, projection test shows negative correlation, real overlap is less than the expectation.
####
#
#scaled.absolute.min.distance.sum.p.value	
#p-value for chromosome-scale null hypothesis as obtained by the permutations of the query points and the mean of the distances to the two closest reference points.
#
###
#
#scaled.absolute.min.distance.sum.lower.tail	
#If TRUE, the query points are closer to the reference points than expected (the absolute distance is lower than the expectation).
###
#
#query.reference.intersection	
#Intersection of reference and query, in bases.
#*sum over spaces
#
#query.reference.union	
#Union of reference and query, in bases.
#*sum over spaces
#
#jaccard.measure	
#Jaccard measure of query and reference overlap.
###
#
#jaccard.measure.p.value	
#The permutation-based evaluation of the p-value for the obtained Jaccard measure, given the null hypothesis of independence.
###
#
#jaccard.measure.lower.tail	
#If TRUE, then Jaccard measure is lower that the expectation (overlap less than expected)
###
#
#The additional values that are returned if keep.distributions=TRUE
#
#relative.distances.data	
#The original relative distances
#*concatenate over spaces
#
#relative.distances.ecdf.deviation.area	
#The real value of the ECDF deviation area to be compared with the permutation to obtain the p-value
##
#
#relative.distances.ecdf.deviation.area.null.list	
#The null distribution
#?????
#
#projection.test	
#List of three values: reference.length is length of a chromosome; reference.coverage is length of that chromosome covered by reference intervale, and query.hits is the number of query points that fall into the reference intervals.
#*sum each value over spaces
#
#absolute.min.distance.data	
#The distribution of query-reference distances
#*concatenate over spaces
#
#absolute.inter.reference.distance.data	
#The distribution of reference-reference distances
#*concatenate over spaces
#
#scaled.absolute.min.distance.sum	
#The value of the sum (i.e. mean) of scaled absolute distances
##
#
#scaled.absolute.min.distance.sum.null.list	
#The null distribution for the scaled absolute distances
#?????
#
#jaccard.measure.null.list	
#The null distribution of Jaccard measures in permutations
#?????
	
	
	
	awhole.space.name="awhole",
	#awhole.space.name is the space name for this operation
			awhole.space.name=awhole.space.name,
	if (awhole.space.name!="awhole") result@config$options$awhole.space.name=awhole.space.name
	awhole.space.name="awhole",
	#awhole.space.name is the space name for this operation
		result[[awhole.space.name]]<-list()
		result[[awhole.space.name]][['query.population']]<-0
		result[[awhole.space.name]][['reference.population']]<-0
		result[[awhole.space.name]][['projection.test']]<-c()
		result[[awhole.space.name]][['projection.test']][['reference.length']]<-0
		result[[awhole.space.name]][['projection.test']][['reference.coverage']]<-0
		result[[awhole.space.name]][['projection.test']][['query.hits']]<-0
		result[[awhole.space.name]][['query.reference.intersection']]<-0
		result[[awhole.space.name]][['query.reference.union']]<-0
			result[[awhole.space.name]][['scaled.absolute.min.distance.sum']]<-0
			result[[awhole.space.name]][['scaled.absolute.min.distance.sum.null.list']]<-c()
		result[[awhole.space.name]][['relative.distances.data']]<-c()
			result[[awhole.space.name]][['absolute.min.distance.data']]<-c()
	 		result[[awhole.space.name]][['absolute.inter.reference.distance.data']]<-c()
		awhole.space.name<-c()
			result[[awhole.space.name]][['relative.distances.data']]<-
						c(result[[awhole.space.name]][['relative.distances.data']],result[[space]][['relative.distances.data']])
			result[[awhole.space.name]][['absolute.min.distance.data']]<-
						c(result[[awhole.space.name]][['absolute.min.distance.data']],result[[space]][['absolute.min.distance.data']])
			result[[awhole.space.name]][['absolute.inter.reference.distance.data']]<-
						c(result[[awhole.space.name]][['absolute.inter.reference.distance.data']],result[[space]][['absolute.inter.reference.distance.data']])
					result[[awhole.space.name]][['projection.test']][[stats]]<-
							result[[awhole.space.name]][['projection.test']][[stats]]+result[[space]][['projection.test']][[stats]]
			result[[awhole.space.name]][['query.reference.intersection']]<-result[[awhole.space.name]][['query.reference.intersection']]+result[[space]][['query.reference.intersection']]
			result[[awhole.space.name]][['query.reference.union']]<-result[[awhole.space.name]][['query.reference.union']]+result[[space]][['query.reference.union']]
			result[[awhole.space.name]][['query.population']]<-result[[awhole.space.name]][['query.population']]+result[[space]][['query.population']]
			result[[awhole.space.name]][['reference.population']]<-result[[awhole.space.name]][['reference.population']]+result[[space]][['reference.population']]
				result[[awhole.space.name]][['scaled.absolute.min.distance.sum']]<-
					result[[awhole.space.name]][['scaled.absolute.min.distance.sum']]+
		the.names<- names(result[[awhole.space.name]])
			space<-awhole.space.name	
			space<-awhole.space.name	
			done_info <- sprintf('chromosome: %s ; %i of %i done',awhole.space.name,done,pb_capacity)
			space<-awhole.space.name
			space<-awhole.space.name
			space<-awhole.space.name
		for (space in c(list.of.spaces,awhole.space.name))
		for (space in c(list.of.spaces,awhole.space.name))
		for (space in c(list.of.spaces,awhole.space.name))
		for (space in c(list.of.spaces,awhole.space.name))
	#now, we want the order of records in result[[awhole.space.name]] to be the same as for all the chromosomes
		result.awhole<-result[[awhole.space.name]];
		result[[awhole.space.name]]<-list();
		# result[[awhole.space.name]] as in all chromosome results
		#result[[awhole.space.name]] as a list()
			result[[awhole.space.name]]<-list();
			result[[awhole.space.name]][[name]]<-result.awhole[[name]]
