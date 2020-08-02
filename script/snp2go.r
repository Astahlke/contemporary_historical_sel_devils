#!/usr/bin/env Rscript
.doc <- "Run SNP2GO.

    snp2go.r --output <output tsv file> --candidates <candidate list>
        --non-candidates <non-candidate list> --gtf <GTF file for the genome>
        --go <GO term file> [--window <window around genes> --fdr <fdr>]

    FORMATS:
        
        candidate / non-candidate list: Should be tab-delimited, with no
            column names. First column is scaffold name (same as GTF file),
            second column is 1-based SNP position.
        GTF file: A GTF file from Ensembl; can be gzipped
        GO term file: A tab-delimited file with column names; Should have
            Ensembl Gene ID and GO Term Accession

    DEPENDS

        igraph, optparse, SNP2GO"

library(optparse)

library(GO.db, lib="/mnt/lfs2/stah3621/bin/R")
library(goProfiles, lib="/mnt/lfs2/stah3621/bin/R")
library(hash, lib="/mnt/lfs2/stah3621/bin/R")
library(SNP2GO, lib="/mnt/lfs2/stah3621/bin/R")
library(GenomicRanges, lib="/mnt/lfs2/stah3621/bin/R")

#==============================================================================#

opts   <- list(make_option('--output'),
               make_option('--candidates'),
               make_option('--non-candidates'),
               make_option('--gtf'),
               make_option('--go'),
               make_option('--window', type='double', default=100000),
               make_option('--fdr', type='double', default=0.05))
parser <- OptionParser(option_list=opts, usage=.doc)
args   <- parse_args(parser)

#-----------------------------------------------------------------------
#
# Setup
#
#-----------------------------------------------------------------------

outfile <- args[['output']]
canfile <- args[['candidates']]
ncfile  <- args[['non-candidates']]
gtffile <- args[['gtf']]
gofile  <- args[['go']]
fdr     <- args[['fdr']]
window  <- args[['window']]

if(is.null(outfile) || is.null(canfile) || is.null(ncfile) ||
   is.null(gtffile) || is.null(gofile)) {
    print_help(parser)
    stop('You left out some arguments')
}

#-----------------------------------------------------------------------
#
# Read data
#
#-----------------------------------------------------------------------

candata <- read.csv(canfile, sep='\t', as.is=TRUE, header=FALSE)
ncdata  <- read.csv(ncfile, sep='\t', as.is=TRUE, header=FALSE)
candidates <- GRanges(seqnames=candata[, 1],
                      ranges=IRanges(candata[, 2], candata[, 2]))
noncandidates <- GRanges(seqnames=ncdata[, 1],
                         ranges=IRanges(ncdata[, 2], ncdata[, 2]))

#-----------------------------------------------------------------------
#
# Run SNP2GO
#
#-----------------------------------------------------------------------

result <- snp2go(gtf=gtffile,
                 goFile=gofile,
                 candidateSNPs=candidates,
                 noncandidateSNPs=noncandidates,
                 FDR=.05,
                 runs=10000,
                 extension=100000,
                 min.regions=1)

#-----------------------------------------------------------------------
#
# Write output file
#
#-----------------------------------------------------------------------

write.table(result$enriched, file=paste0(outfile, ".tsv"), sep='\t', row.names=FALSE,
            quote=FALSE)

result$enriched$GO.def    <- as.character(result$enriched$GO.def)
result$enriched$child.GOs <- as.character(result$enriched$child.GOs)
result$enriched$GO        <- as.character(result$enriched$GO)

# Get all enriched GO terms of GFF analysis:
gff.significant.terms <- result$enriched$GO
        
# Get the candidate SNPs and regions associated w GO terms  
gffcans.df <- data.frame()
gffregions.df <- data.frame()

for(i in 1:length(gff.significant.terms)){
	gff.cans <- as.data.frame(candidates[unlist(as.list(result$go2ranges[["candidates"]][gff.significant.terms[i]]))])
	gff.cans$GOterm <- gff.significant.terms[i]
	gffcans.df <- rbind(gff.cans, gffcans.df)

	gff.regions <- as.data.frame(
		result$regions[unlist(as.list(
		result$go2ranges[["regions"]][gff.significant.terms[i]]))])
	gff.regions$GOterm <- gff.significant.terms[i]
	gffregions.df <- rbind(gff.regions, gffregions.df)
}

write.table(gffcans.df, file = paste0(outfile,".candidateSNPs.tsv"), sep='\t', row.names=FALSE,
            quote=FALSE)
write.table(gffregions.df, file = paste0(outfile, ".GOregions.tsv", sep='\t', row.names=FALSE,
            quote=FALSE)
#-----------------------------------------------------------------------
#
# Group GO terms
#
#-----------------------------------------------------------------------
library(igraph, lib="/mnt/lfs2/stah3621/bin/R")

ch <- strsplit(as.character(result$enriched$child.GOs), ',')
ch <- lapply(ch, function(x) if(length(x) == 0) {return(NA)} else {return(x)})
elist <- cbind(rep(as.character(result$enriched$GO), 
                   times=sapply(ch, length)),
               unlist(ch))
nochildren <- elist[is.na(elist[, 2]), 1]
noparent   <- nochildren[!nochildren %in% elist[, 2]]
g <- graph.edgelist(elist[!is.na(elist[, 2]), ], directed=TRUE)
g <- add.vertices(g, length(noparent), attr=list(name=noparent))
components <- clusters(g, 'weak')
members <- components$membership
vnames <- vertex.attributes(g)$name
f <- file(paste0(outfile, '.gotree'), 'wb')
ff <- file(paste0(outfile, '.goterms'), 'wb')
cat('cluster\tlevel\tGO\tP\tFDR\tp.L\tp.G\tg\tG\tnc\tmc\tGO.def\tchild.GOs\n', file=ff)
for(cl in unique(members)) {
    terms <- vnames[members == cl]
    levels <- unlist(result$termlevel[terms])
    l <- 0
    for(i in sort(unique(levels))) {
        terms_at_level <- names(levels)[levels == i]
        for(tm in terms_at_level) {
            cat(cl, '\t', i, '\t',
                paste(result$enriched[result$enriched$GO == tm, ], collapse='\t'),
                '\n', sep='', file=ff)
            cat(cl, '\t', i, '\t', rep('\t', times=l), tm, '\n', sep='', file=f)
        }
        l <- l + 1
    }
}
close(f)
close(ff)
