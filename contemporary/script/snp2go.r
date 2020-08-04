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

# install.packages("igraph", lib="/mnt/lfs2/stah3621/bin/R")
library(GO.db, lib="/mnt/lfs2/stah3621/bin/R")
library(goProfiles, lib="/mnt/lfs2/stah3621/bin/R")
library(hash, lib="/mnt/lfs2/stah3621/bin/R")
library(SNP2GO, lib="/mnt/lfs2/stah3621/bin/R")     #load the package
library(GenomicRanges, lib="/mnt/lfs2/stah3621/bin/R")
library(igraph, lib="/mnt/lfs2/stah3621/bin/R")

suppressMessages(library(optparse))
suppressMessages(library(SNP2GO))
suppressMessages(library(igraph))

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
                 FDR=fdr,
                 runs=10000,
                 extension=window,
                 min.regions=1)

#-----------------------------------------------------------------------
#
# Write output file
#
#-----------------------------------------------------------------------

write.table(result$enriched, file=outfile, sep='\t', row.names=FALSE,
            quote=FALSE)

result$enriched$GO.def    <- as.character(result$enriched$GO.def)
result$enriched$child.GOs <- as.character(result$enriched$child.GOs)
result$enriched$GO        <- as.character(result$enriched$GO)


#-----------------------------------------------------------------------
#
# Group GO terms
#
#-----------------------------------------------------------------------

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

# I decided not to use this piece of code for now, but
# I left it in case it is useful later (yes, I know this is
# messy).
# see also unfold.tree and denom.tree for useful ideas
#par(mar=c(0, 0, 0, 0))
#plot(0, 0, xlim=c(0, 1), ylim=c(0, 1), type='n')
#h <- 1 / nrow(result$enriched)
#w <- 1 / max(sapply(unique(components$membership), 
#             function(cluster) {
#                 terms <- vnames[components$membership == cluster] 
#                 levels <- unlist(result$termlevel[terms])
#                 max(levels) - min(levels) + 1
#             }))
#y <- 1
#for(cl in unique(members)) {
#    x <- 0
#    terms <- vnames[members == cl]
#    levels <- unlist(result$termlevel[terms])
#    l <- 0
#    level0 <- min(levels)
#    leveli <- structure(1:length(levels), names=levels)
#    for(i in sort(unique(levels))) {
#        terms_at_level <- names(levels)[levels == i]
#        for(tm in terms_at_level) {
#            x <- w * l
#            points(x, y, cex=0.5, col='blue', pch=19)
#            text(x, y, cex=0.5, result$enriched$GO.def[result$enriched$GO == tm],
#                 adj=0)
#            y <- y - h
#        }
#        l <- l + 1
#    }
#}
