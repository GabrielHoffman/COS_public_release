# Gabriel Hoffman
# March 21, 2017
#
# Use enrichment_analysis.R to compute enrichment of genes that vary across individuals
# to genes that are eQTL's in CommonMind


#################
# Prepare files #
#################

suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library(variancePartition))

synapseLogin()

# get variancePartition
# first syn7254046  
# second syn9771709
# public release: syn9926786 (June 1, 2017)
varPart = read.table(getFileLocation(synGet( 'syn9926786' )), header=TRUE, row.names=1, sep="\t")[,-1]
varPart = sortCols( varPart )

# get eQTL results
eQTLInfo = read.table(getFileLocation(synGet( 'syn8000326' )), header=TRUE)
eQTLInfo = unique(eQTLInfo)

# get genes common to both datasets
commonGenes = intersect(rownames(varPart), eQTLInfo$gene)

varPart = varPart[rownames(varPart) %in% commonGenes,]
eQTLInfo = eQTLInfo[eQTLInfo$gene %in% commonGenes,]

dim(varPart)
dim(eQTLInfo)

cutoff = sort(eQTLInfo$FDR)[2001]
eQTL_genes = with(eQTLInfo, gene[FDR <cutoff])
 length(eQTL_genes)

# combine Sex with Donor
varPart2 = varPart
varPart2$Donor = with(varPart2, Donor + Sex + Dx)
# varPart2$Sex = c()

geneInfo = read.table(getFileLocation(synGet( 'syn9909914' )), header=TRUE, row.names=1, sep="\t")

# drop genes on the chrX, chrY, chrM
idx = rownames(varPart) %in% rownames(geneInfo)[geneInfo$Chrom %in% c('chrX', 'chrY', 'chrM')]


# BED1
write.table( data.frame(rownames(varPart)[!idx], 1, 10), "/hpc/users/hoffmg01/psychgen_ips/run/eqtl_enrichment/varPart.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# score1
write.table( varPart2[!idx,], "/hpc/users/hoffmg01/psychgen_ips/run/eqtl_enrichment/varPart.tsv", quote=FALSE, sep="\t", row.names=FALSE)

# BED2
write.table( data.frame(eQTL_genes, 1, 10), "/hpc/users/hoffmg01/psychgen_ips/run/eqtl_enrichment/CMC_EQTL.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

cols = cbind(colnames(sortCols(varPart)), c(ggColorHue(ncol(varPart) - 1), "grey85"))
write.table(cols, "/hpc/users/hoffmg01/psychgen_ips/run/eqtl_enrichment/colors.csv", 
	quote=FALSE, row.names=FALSE, col.names=FALSE, sep=',')



#######################
# Enrichment analysis #
#######################

# Compute enrichments based on sorting of genes in $DIR/varPart.tsv.  Genes in the file are listed in $DIR/varPart.bed along 
# with (fake) coordinates.  Enrichment is computed based on the set of genes in $DIR/CMC_EQTL.bed  
# $NPERM permutations are used  
# The colors of lines corresponding to $DIR/varPart.tsv are given in $DIR/colors.csv 

DIR=/hpc/users/hoffmg01/psychgen_ips/run/eqtl_enrichment/

SCRIPT_RUN=~/scripts/psychENCODE/EpiMap/enrichment_analysis.R

NPERM=10000

$SCRIPT_RUN --bed1 $DIR/varPart.bed --scores1 $DIR/varPart.tsv --bed2 $DIR/CMC_EQTL.bed --nperm $NPERM --confidenceInterval --nthreads 60 --xlim '0 100' --log --colors $DIR/colors.csv  --min_set_size 100 --npoints 40 --outfile $DIR/enrichment_QTL_CTC_${NPERM}.pdf 


