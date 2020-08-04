#!/bin/bash

### This script runs multiple paml codeml models and then parses the results 
module load paml
module load python/3.3.2

date="091817" ## dhange this to reflect new results

dir="/mnt/home/stah3621/molevo/"
out="/mnt/home/stah3621/molevo/paml_out/"

for f in ${dir}wholegenome/ENSSHAG00000${1}.phylip; do
geneID=$(echo $f | cut -f 7 -d '/' | cut -f 1 -d .)
echo $geneID
#
### 0 script
### 1 = alignment 2 = tree
### 3 = model type  *  0:one, 1:b, 2:2 or more dN/dS ratios for branches
### 4 = fix_omega * 1:yes, 0
### 5  = omega = 1
### 6 codon_freq model * 0: all assumed equal, 1: from average 2: from average at codon positions; 3: free parameters
### 7 = out

python ${dir}scripts/my_branch-site.py ${dir}wholegenome/${geneID}.phylip ${dir}maycollado.paml.nwk 2 0 0 2 ${dir}paml_out/${geneID}.posi.branchsite.out

python ${dir}scripts/my_branch-site.py ${dir}wholegenome/${geneID}.phylip ${dir}maycollado.paml.nwk 2 1 1 2 ${dir}paml_out/${geneID}.null.branchsite.out ## null
done;


for f in ${dir}paml_out/ENSSHAG00000${1}.null.branchsite.out; do
#geneID=$(echo $f | cut -f 8 -d '/'| cut -d '.' -f 1 )
echo $geneID

python ${dir}scripts/parse_branchsite-null.py ${dir}paml_out/${geneID}.null.branchsite.out ${geneID} ${out}paml_nullBS_mars_wholegenome_${date}.tsv


for r in ${dir}paml_out/ENSSHAG00000${1}.posi.branchsite.out; do
python ${dir}scripts/parse_branchsite-posi.py ${dir}paml_out/${geneID}.posi.branchsite.out ${geneID} ${out}paml_posiBS_mars_wholegenome_${date}.tsv

done;

#mkdir ${dir}HyPhyout_${date}

#p="hyphy"

for f in ${dir}wholegenome/ENSSHAG00000${1}.phylip; do
#        samp=$(echo $f | cut -d '/' -f 7 | cut -d '.' -f 1 )
echo $geneID

(echo 1; echo 1; echo 2; echo ${f}; echo ${dir}maycollado.HyPhy.nwk; echo 7; echo 2; echo d;  echo ${dir}${samp}.BSREL.out) | HYPHYMP LIBPATH=${LIBPATH} $LIBPATH/TemplateBatchFiles/BranchSiteREL.bf

echo "done with adaptive branch site random effects model"
done;
