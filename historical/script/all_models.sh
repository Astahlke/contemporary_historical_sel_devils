#!/bin/bash

### This script runs multiple paml codeml models and then parses the results 
module load paml
module load python/3.3.2

date="1003117" ## dhange this to reflect new results

dir="/mnt/home/stah3621/scratch/molevo/"
out="/mnt/home/stah3621/scratch/molevo/results/"

#for f in ${dir}addkoala/*.phylip; do
#geneID=$(echo $f | cut -f 8 -d '/' | cut -f 1 -d .)
#echo $geneID
#echo $f
## 0 script
## 1 = alignment 2 = tree
## 3 = model type  *  0:one, 1:b, 2:2 or more dN/dS ratios for branches
## 4 = fix_omega * 1:yes, 0
## 5  = omega = 1
## 6 codon_freq model * 0: all assumed equal, 1: from average 2: from average at codon positions; 3: free parameters
## 7 = out

#python ${dir}scripts/my_branch-site.py $f  ${dir}addkoala/maycollado.koala.paml.nwk 2 0 0 2 ${dir}paml_out/${geneID}.koala.posi.branchsite.out
#python ${dir}scripts/my_branch-site.py $f ${dir}addkoala/maycollado.koala.paml.nwk 2 1 1 2 ${dir}paml_out/${geneID}.koala.null.branchsite.out ## null

#done;


for f in ${dir}paml_out/*koala.null.branchsite.out; do
geneID=$(echo $f | cut -f 8 -d '/'| cut -d '.' -f 1 )
echo $geneID 
echo $f
python ${dir}scripts/parse_branchsite-null.py $f ${geneID} ${out}paml_nullBS_can_koala_${date}.tsv
done;

for r in ${dir}paml_out/*koala.posi.branchsite.out; do
geneID=$(echo $r | cut -f 8 -d '/'| cut -d '.' -f 1 )
echo $geneID
echo $r
python ${dir}scripts/parse_branchsite-posi.py $r ${geneID} ${out}paml_posiBS_can_koala_${date}.tsv
done;
