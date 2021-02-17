#!/usr/bin/env Rscript
#
# Convert quantiles to p-values.
#
#    quantile_to_p.r --output <output file> <input file>
#
# This script makes sure that there are no p-values of zero or one.
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

opts = list(make_option('--output'))
args = parse_args(OptionParser(option_list=opts),
                  positional_arguments=TRUE)

infile  = args$args[1]
outfile = args$options[['output']]

data = read.csv(infile, sep='\t', header=FALSE, as.is=TRUE)

# I re-rank the quantiles here to avoid getting a value of 0,
# which would mess up the next step.
adj_p = ((length(data[, 8]) + 1) - rank(data[, 8])) / length(data[, 8])
adj_p[adj_p == 1] = mean(sort(adj_p, decreasing=TRUE)[1:2])

output = data.frame(data, adj_p)
names(output) = c('scaffold', 'start', 'end', 'snp', 'score',
                  'p', 'q', 'quantile', 'adjp')
write.table(output, file=outfile, sep='\t', col.names=TRUE,
            row.names=FALSE, quote=FALSE)
