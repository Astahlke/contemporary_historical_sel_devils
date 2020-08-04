#!/usr/bin/env Rscript
#
# Correct p-values using a genomic inflation factor.
#
#    genomic_inflation_factor_correction.r --output <output file>
#       [--fdr] [--lambda <float> (-1)] [--min-p (0.0)] [--max-p (1.0)]
#       [--double-p] <input file>
#
# Because the downstream step can't handle p-values of 0, you can set a
# min. p-value to change all the zeroes to. Similarly, --max-p sets the
# p-value for p-values that are 1 (again, bad for DCMS).
#
# INPUT FILE FORMAT
#
# Tab-delimited, with no header, and columns in this order:
#
#   scaffold, start, end, snp, score, p-value, q-value, quantile
#

#==============================================================================#

library(optparse)

#==============================================================================#

opts = list(make_option('--output'),
            make_option('--lambda', type='double', default=-1),
            make_option('--min-p', type='double', default=0),
            make_option('--max-p', type='double', default=1),
            make_option('--double-p', default=FALSE, action='store_true'),
            make_option('--fdr', default=FALSE, action='store_true'))
args = parse_args(OptionParser(option_list=opts),
                  positional_arguments=TRUE)

infile  = args$args[1]
outfile = args$options[['output']]
lambda  = args$options[['lambda']]
doubp   = args$options[['double-p']]
usefdr  = args$options[['fdr']]

data = read.csv(infile, sep='\t', header=FALSE, as.is=TRUE)

# Calculate lambda from the data
if(lambda < 0) {
    pvals = data[, 6]
    if(doubp)
        pvals = 2 * pvals
    z2 = qchisq(1 - pvals, 1, lower.tail=TRUE)
    lambda = median(z2, na.rm=TRUE) / qchisq(0.5, 1)
    adj_p = pchisq(z2 / lambda, 1, lower.tail=FALSE)
    adj_p[adj_p == 0] = args$options[['min-p']]
    adj_p[adj_p == 1] = args$options[['max-p']]
}

if(usefdr) {
    library(fdrtool)
    qvals = numeric(length(adj_p)) * NaN
    qvals[!is.na(adj_p)] = fdrtool(adj_p[!is.na(adj_p)], 'pvalue')$qval
    output = data.frame(data, adj_p, qvals)
    names(output) = c('scaffold', 'start', 'end', 'snp', 'score',
                      'p', 'q', 'quantile', 'adjp', 'fdr')
    output = output[order(output$scaffold, output$start), ]
    write.table(output, file=outfile, sep='\t', col.names=TRUE,
                row.names=FALSE, quote=FALSE)
} else {
    output = data.frame(data, adj_p)
    names(output) = c('scaffold', 'start', 'end', 'snp', 'score',
                      'p', 'q', 'quantile', 'adjp')
    output = output[order(output$scaffold, output$start), ]
    write.table(output, file=outfile, sep='\t', col.names=TRUE,
                row.names=FALSE, quote=FALSE)
}
cat('file:', infile, 'lambda:', lambda, '\n')
