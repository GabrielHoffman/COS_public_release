#! /bin/env Rscript
# Gabriel Hoffman
# March 14, 2017
#
# Compute enrichments between region in two BED files, where the regions in the first 
# are sorted based on a score

# For example, compute enrichment for peaks that vary across individual and 
# epiEQTL's from Grubert, et al. Cell 2015
#
# METHODS: Enrichment of Jaccard similarities in peaks that vary across
# Indiviuals. Compare Jaccard similariy in peaks that are epiQTL's in Grubert and 
# Compare to all peaks in dataset as baseline
#
# INPUT: BED file of regions, each with correspondings scores.  Enrichments are computed with 
# a second BED file based on the regions in the first file than pass a score cutoff

# Uses: regioneR

library(getopt)
spec = matrix(c(
	'bed1',  				'a', 1, "character",
	'scores1',  			'b', 1, "character",
	'bed2',  				'z', 1, "character",
	'nperm',				'p', 1, "integer",
	'confidenceInterval',	'l', 0, "logical",
	'offset',				's', 0, "integer",
	'xlim',					'x', 1, "character",
	'log',					'g', 0, "logical",
	'ylim',					'y', 1, "character",
	'nthreads',				't', 1, "integer",
	'npoints',				'n', 1, "integer",
	'colors',				'r', 1, "character",
	'outfile',				'c', 1, "character",
	'min_set_size',			'e', 1, "integer"
), byrow=TRUE, ncol=4);
opt = getopt(spec);


suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(regioneR))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))

# Set parallel options
doParallel = TRUE
ntimes = opt$nperm
mc.cores = 1
registerDoParallel(opt$nthreads)

min_set_size = ifelse( is.null(opt$min_set_size), 20, opt$min_set_size )


#############
# Read data #
#############

BED1 = toGRanges(opt$bed1)
BED2 = toGRanges(opt$bed2) 
scoreDF = read.table( opt$scores1, sep="\t", header=TRUE)
colLocalIn = read.csv( opt$colors, header=FALSE, stringsAsFactors=FALSE)

colLocal = colLocalIn[,2]
names(colLocal) = colLocalIn[,1]

# sort colLocal based on scoreDF
colLocal = colLocal[match(colnames(scoreDF), names(colLocal))]


if( any(!colnames(scoreDF) %in% names(colLocal)) ){
	stop("color names must match score names")
}


if( nrow(scoreDF) != length(BED1) ){
	stop("scores1 and bed1 must have same number of rows")
}

# compute jaccard similarity
jaccardWidth = function(A, B,...){
	# length( commonRegions(A, B) ) / length( mergeRegions(A, B) )

	U = sum(width(union(A, B)))
	I = sum(width(intersect(A, B)))

	I / (U-I)
}		

jaccardCount = function(A, B,...){
	# length( commonRegions(A, B) ) / length( mergeRegions(A, B) )


	U = length(union(A, B))
	I = length(intersect(A, B))

	I / (U-I)
}		


# metric =  jaccardWidth
metric = jaccardCount


run_enrichment_analysis = cmpfun(function( scoreDF, BED1, BED2, npoints){
	
	# cutoffs = seq(min(scoreDF), max(scoreDF), length.out=npoints)

	# increase the number of cutoffs and make them close to the max value
	# cutoffs = apply(scoreDF, 2, function(x) seq(0, max(x)*.999, length.out=npoints))
	# cutoffs = sort( cutoffs )
	cutoffs = seq(0, .99, length.out=npoints )

	# for all colums in score1
	ptList = foreach( variable = colnames(scoreDF)) %do% {
		# for each cutoff
		res = foreach( cutoff = cutoffs ) %dopar% {
			cat("\rvariable:", variable, "\tcutoff:", cutoff, "           ")

			# get intervals that satisfy cutoff
			idx_A = which( scoreDF[,variable] >= cutoff )

			# if A is empty, add interval with max score
			if( length(idx_A) < min_set_size){
				idx_A = which( scoreDF[,variable] >= max(scoreDF[,variable]) )
				observedMetric = 0
			}else{

				# May 12, 2017
				# Simpler and faster permutation scheme implemented directly rather than thru permTest
				# This samples X regions from BED1 to evaluate permutations
				observedMetric = metric( BED1[idx_A,], BED2)
			}

			n_total_peaks = length(BED1)
			permutedMetric = sapply(1:ntimes, function(i){
				idx = sample.int(n_total_peaks, length(idx_A), replace=FALSE)
				metric( BED1[idx,], BED2)
			})

			# hist( observedMetric / permutedMetric )

			pt = list(metric = list(observed = observedMetric, permuted=permutedMetric))
			
			# for subset of entires in BED1, compute overlap with BED1.  Draw permutations from BED1
			# pt = permTest(A=BED1[idx_A,], B=BED2, universe=BED1, ntimes=ntimes, randomize.function=resampleRegions, evaluate.function=metric, mc.set.seed=FALSE, mc.cores=mc.cores, doParallel=doParallel)

			pt
		}
		names(res) = format(cutoffs)
		res
	}
	names(ptList) = colnames(scoreDF)	

	# get Fold enrichement
	ratioMetric = foreach( variable = names(ptList), .combine=cbind ) %do% {
		# get enrichment ratio
		res = sapply(ptList[[variable]], function(x) median( x$metric$observed / x$metric$permuted, na.rm=TRUE))

		# format, and replace invalid entries with NA
		res = as.matrix(res)
		colnames(res) = variable
		res[is.nan(res) | is.infinite(res)] = NA
		res
	}

	ratioMetric = data.frame(cutoff = as.numeric(rownames(ratioMetric)), ratioMetric)
	rownames(ratioMetric) = c()


	# get confidence intervals
	ratioMetricCI = foreach( variable = names(ptList), .combine=rbind ) %do% {
		# get confidence interval for enrichment ratio
		resTop = sapply(ptList[[variable]], function(x) quantile( x$metric$observed / x$metric$permuted, .95, na.rm=TRUE))
		resBottom = sapply(ptList[[variable]], function(x) quantile( x$metric$observed / x$metric$permuted, .05, na.rm=TRUE))

		resTop[resTop ==0 |is.nan(resTop)] = NA
		resBottom[resBottom ==0 |is.nan(resBottom)] = NA

		res = data.frame(variable = variable, cutoff = as.numeric(names(ptList[[variable]])), top = resTop, bottom = resBottom )
		rownames(res) = c()

		# assign missing values the previous value
		while( any(is.infinite(res$top)) ){
			idx = which(is.infinite(res$top))
			if(length(idx) > 0){
				res$top[idx] = res$top[idx-1]
			}
		}

		while( any(is.infinite(res$bottom)) ){
			idx = which(is.infinite(res$bottom))
			if(length(idx) > 0){
				res$bottom[idx] = res$bottom[idx-1]
			}
		}
		res
	}

	# set ratioMetricCI to NA where ratioMetric is NA
	for( key in levels(ratioMetricCI$variable) ){
		i = is.na(ratioMetric[[key]])

		idx = which(ratioMetricCI$variable == key)
		ratioMetricCI$top[idx][i] = NA
		ratioMetricCI$bottom[idx][i] = NA
	}

	list(ratioMetric = ratioMetric, ratioMetricCI = ratioMetricCI)
})

################
# run analysis #
################

npoints = ifelse( is.null(opt$npoints), 40, opt$npoints)

res = run_enrichment_analysis( scoreDF, BED1, BED2, npoints)

#################
# write results #
#################

file = paste0(tools::file_path_sans_ext(opt$outfile), '.ratio.tsv')
write.table( res$ratioMetric, file, quote=FALSE, sep="\t", row.names=FALSE)

file = paste0(tools::file_path_sans_ext(opt$outfile), '.ratioCI.tsv')
write.table( res$ratioMetricCI, file, quote=FALSE, sep="\t", row.names=FALSE)


##############
# Make plots #
##############

if( opt$confidenceInterval ){
	ylim = range(res$ratioMetricCI[,-c(1,2)], na.rm=TRUE)
}else{
	ylim = range(res$ratioMetric[,-1], na.rm=TRUE)
}

if( !is.null(opt$log) ){
	ylim = log2(ylim)
}

if( ! is.null(opt$ylim) ){
	ylim = as.numeric(unlist(strsplit(opt$ylim, " ")))
	if( length(ylim) == 1){
		ylim = c(-ylim, ylim)
	}
}

if( ! is.null(opt$xlim) ){
	xlim =  as.numeric(unlist(strsplit(opt$xlim, " ")))
}

cutoffs = res$ratioMetric[,1]

if( !is.null(opt$log) ){

	# log2 space
	pdf(opt$outfile)

	par(pty="s") # square plot

	plot(1, type='n', xlim=xlim, ylim=ylim, xlab=expression(Peaks~with~variance~percentage>=cutoff), ylab=expression(log[2]~Fold~enrichment), xaxt='n', cex.lab=1)
	abline(h=0, col="black", lty=2, lwd=3)
	axis(1, at = seq(0, 100, length=6), label=paste0("  ", seq(0, 100, length=6), "%"))
	matplot(cutoffs*100, log2(res$ratioMetric[,-1]), pch=20, type='l', lty=1, lwd=3, col=colLocal, add=TRUE)

	legend('topleft', colnames(res$ratioMetric)[-1], col=colLocal, lty=1, lwd=3)

	if( opt$confidenceInterval ){
		for( key in levels(res$ratioMetricCI$variable) ){
	
			idx = which(res$ratioMetricCI$variable == key)
			idx2 = !is.na(res$ratioMetricCI$top[idx]) & !is.na(res$ratioMetricCI$bottom[idx])

			ytop = log2(res$ratioMetricCI$top[idx][idx2])
			ybottom = log2(res$ratioMetricCI$bottom[idx][idx2])
			x = res$ratioMetricCI$cutoff[idx][idx2] * 100

			polygon(c(rev(x), x),c(rev(ybottom),ytop),col=adjustcolor(colLocal[key], alpha.f=0.15),border=NA)
		}
	}
	dev.off()
}else{

	# original space
	pdf(opt$outfile)
	
	par(pty="s") # square plot

	plot(1, type='n', xlim=xlim, ylim=ylim, xlab=expression(Peaks~with~variance~percentage>=cutoff), ylab="Fold enrichment", xaxt='n', cex.lab=1.3)
	abline(h=1, col="black", lty=2, lwd=3)
	axis(1, at = seq(0, 100, length=6), label=paste0("  ", seq(0, 100, length=6), "%"))
	matplot(cutoffs*100, res$ratioMetric[,-1],  pch=20, type='l', lty=1, lwd=3, col=colLocal, add=TRUE)

	legend('topleft', colnames(res$ratioMetric)[-1], col=colLocal, lty=1, lwd=3)

	if( opt$confidenceInterval ){
		for( key in levels(res$ratioMetricCI$variable) ){
	
			idx = which(res$ratioMetricCI$variable == key)
			idx2 = !is.na(res$ratioMetricCI$top[idx]) & !is.na(res$ratioMetricCI$bottom[idx])

			ytop = res$ratioMetricCI$top[idx][idx2]
			ybottom = res$ratioMetricCI$bottom[idx][idx2]
			x = res$ratioMetricCI$cutoff[idx][idx2] * 100

			polygon(c(rev(x), x),c(rev(ybottom),ytop),col=adjustcolor(colLocal[key], alpha.f=0.15),border=NA)
		}
	}
	dev.off()
}

