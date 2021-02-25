#!/usr/bin/env Rscript
#
# Run the s_lattice code on a directory of files
#
# This is identical to slattice.r (as of 03 Apr 2015), but it includes
# 95% confidence intervals. The confidence intervals assume no dominance.
#
#==============================================================================#

suppressMessages(library(optparse))
scriptdir <- dirname(gsub("--file=", "", commandArgs()[grepl("--file", commandArgs())]))
cwd <- getwd()
setwd(file.path(scriptdir, 's_lattice'))
source(file.path(scriptdir, 's_lattice', 'include.R'))
setwd(cwd)

#==============================================================================#

opts <- list(make_option('--output'),
             make_option('--ne', type='double'))
parser <- OptionParser(option_list=opts)
args <- parse_args(parser, positional_arguments=TRUE)

ne <- args$options$ne

out <- file(args$options$output, 'w')
cat('locus\ts\tlower\tupper\tp\n', file=out, sep='')
for(fname in list.files(args$args[1], full.names=TRUE)) {
    locus <- gsub('(.+)\\..+', '\\1', basename(fname))
    d <- read.csv(fname, header=TRUE, sep='\t', as.is=TRUE)
    ss <- tryCatch(estimate.s(obs=d, Ne=ne)$s,
                   error=function(e) NaN)
    if(!is.na(ss)) {
        ci <- tryCatch(find.confidence.interval(d, ne, ss, alpha=0.05, h=0.5),
                       error=function(e) c(NaN, NaN))
        p <- tryCatch(p.value(d, ne, ss, h=0.5),
                      error=function(e) NaN)
    } else {
        ci <- c(NaN, NaN)
        p <- NaN
    }
    cat(locus, '\t', ss, '\t', ci[1], '\t', ci[2], '\t', p, '\n',
        file=out, sep='')
}
close(out)

cat('Done\n')
