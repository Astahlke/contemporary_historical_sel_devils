#!/usr/bin/env python
"""
    Combine annotations"
"""

#==============================================================================#

import argparse
import csv
import re
import collections
import os.path as osp

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--bedtools', default=False, action='store_true')
parser.add_argument('--bed', default=False, action='store_true')
parser.add_argument('--dir', default=False, action='store_true')
parser.add_argument('input', nargs='+')
args = parser.parse_args()

genes = collections.defaultdict(int)
gene_sets = {}
descriptions = {}
datasets = []

for fname in args.input:
    if args.dir:
        dataset = osp.basename(osp.dirname(fname)) + \
                  '.' + osp.splitext(osp.basename(fname))[0]
    else:
        dataset = osp.splitext(osp.basename(fname))[0]
    datasets.append(dataset)
    gene_sets[dataset] = set()
    with open(fname, 'r') as handle:
        if args.bedtools or args.bed:
            rdr = csv.reader(handle, delimiter='\t')
        else:
            rdr = csv.DictReader(handle, delimiter='\t')
        for row in rdr:
            # I don't know why I had this in a try block:
            #try:
            #    gene = row['gene'].split(':')[0]
            #    description = ':'.join(row['gene'].split(':')[1:])
            #    descriptions[gene] = description
            #except ValueError:
            #    gene, description = None, None
            if args.bedtools:
                rg = row[7].split(':')
            elif args.bed:
                rg = [row[3]]
            else:
                rg = row['gene'].split(':')
            gene = rg[0]
            description = ':'.join(rg[1:])
            if len(gene) == 0:
                gene = None
            else:
                descriptions[gene] = description
            if gene is not None and gene not in gene_sets[dataset]:
                genes[gene] += 1
                gene_sets[dataset].add(gene)

with open(args.output, 'w') as out:
    out.write('gene\tdescription\ttotal\t' + '\t'.join(datasets) + '\n')
    for gene, description in descriptions.iteritems():
        out.write(gene + '\t' + str(genes[gene]) + '\t' + description + '\t' +
                  '\t'.join('1' if gene in gene_sets[d] else '0' for d in datasets) +
                  '\n')
