#!/usr/bin/env Rscript
#
# This script calculates a composite test statistic using the column
# labeled "adjp" in a set of tsv files.
#
# css.r [--min-tests (1)] --output <output file> <input files...>
#
# DETAILS
#
# Ma method (DCMS - used here):
#
#   For each statistic, calculate a ratio for each SNP:
#       log((1-p)/p), where p is the p-value (all tests combined with
#   DCMS must have a p-value). In this case I use the "adjp" column
#   as the p-value.
#
# NOTE
#
#   A previous version of this script expected files without a header,
#   and without and "adjp" column. This script is different. The expected
#   columns are:
#
#       scaffold    the scaffold / contig / chromosome
#       start       (0-based)
#       adjp        the p-value-like quantity to use
#
#   Other columns are ok, but will be ignored.
#

#==============================================================================#

library(optparse)

#==============================================================================#

options('scipen'=50)

opts = list(make_option('--output'),
            make_option('--min-tests', type='double', default=1))
parser = OptionParser(option_list=opts)
args = parse_args(parser, positional_arguments=TRUE)

outfile   = args$options[['output']]
min_tests = args$options[['min-tests']]
infiles   = args$args

indata = vector('list', length=length(infiles))
for(i in seq_along(infiles)) {
    d = read.csv(infiles[i], sep='\t', header=TRUE, as.is=TRUE)
    indata[[i]] = d[, c('scaffold', 'start', 'adjp')]
}

# First merge the datasets
snps = unique(unlist(lapply(indata,
                            function(x) paste0(x[, 'scaffold'], '-', x[, 'start']))))
id_table  = data.frame('scaffold'=gsub('-.+', '', snps),
                       'start'=as.numeric(gsub('.+-', '', snps)))

data = vector('list', length=length(indata))
for(i in seq_along(data)) {
    tmp = merge(id_table, indata[[i]][, c('scaffold', 'start', 'adjp')],
                all=TRUE, stringsAsFactors=FALSE)
    data[[i]] = tmp[order(tmp[, 'scaffold'], tmp[, 'start']), ]
}

# Check to make sure the merging step worked
scaffolds = data[[1]][, 'scaffold']
positions = data[[1]][, 'start']
if(length(data) > 1) {
    passing = sapply(data[-1], function(x) {
                     all(x[, 'scaffold'] == scaffolds) &&
                     all(x[, 'start'] == positions)
                    })
    if(!all(passing))
        stop('Merging did not work')
}

# Get the number of tests with results at each SNP
nstats = numeric(nrow(data[[1]]))
for(snpidx in seq_along(nstats)) {
    n = 0
    for(statidx in 1:length(data)) {
        if(!is.na(data[[statidx]][snpidx, 'adjp'])) {
            n = n + 1
        }
    }
    nstats[snpidx] = n
}

# Calculate the weighting factor for each statistic using only
# SNPs with the minimum number of tests.
weights = numeric(length(data))
wmat = matrix(nrow=length(data), ncol=length(data))
for(i in 1:length(data)) {
    sum_cors = 0
    for(j in 1:length(data)) {
        r = cor.test(data[[i]][nstats >= min_tests, 'adjp'],
                     data[[j]][nstats >= min_tests, 'adjp'])$estimate
        sum_cors = sum_cors + abs(r)
        wmat[i, j] = r
    }
    weights[i] = sum_cors
}

# Then calculate the composite score for each snp.
composite = numeric(nrow(data[[1]]))
meancomp = numeric(nrow(data[[1]]))
for(snpidx in 1:length(composite)) {
    if(nstats[snpidx] < min_tests) {
        composite[snpidx] = NA
        meancomp[snpidx] = NA
        next
    }
    comp = 0
    p = numeric(length(data))
    for(statidx in 1:length(data)) {
        p[statidx] = data[[statidx]][snpidx, 'adjp']
    }
    for(statidx in 1:length(data)) {
        # This only uses the statistics calculated for a particular
        # SNP to get the weight - it is debatable whether this
        # is really a good idea.
        w = sum(abs(wmat[statidx, !is.na(p)]))
        x = log((1 - p[statidx]) / p[statidx]) / w
        if(!is.na(x)) {
            comp = comp + x
        }
    }
    if(is.infinite(comp))
        stop('infinite')
    composite[snpidx] = comp
    meancomp[snpidx] = comp / sum(!is.na(p))
}

# Write the output
output = data.frame('scaffold'=scaffolds, 'start'=positions,
                    'end'=data[[1]][, 'start'] + 1, 'nstats'=nstats,
                    'composite'=composite,
                    'mean_composite'=meancomp)[nstats >= min_tests, ]
write.table(output, file=outfile, sep='\t', col.names=TRUE,
            row.names=FALSE, quote=FALSE)
write.table(data.frame('weight'=weights, 'file'=infiles),
            file=paste0(outfile, '.weights'), sep='\t', col.names=TRUE,
            row.names=FALSE, quote=FALSE)
write.table(wmat, file=paste0(outfile, '.cors'), sep='\t', col.names=FALSE,
            row.names=FALSE, quote=FALSE)
