#!/usr/bin/env bash

dir="/mnt/home/stah3621/molevo/"

module load muscle
module load emboss

cd ${dir}wholegenome 

for f in ENSSHAG*aa.fasta; do
geneID=$(echo $f | cut -f1 -d.)

echo $geneID 
muscle -in ${geneID}.aa.fasta -out ${geneID}.aa.align

grep '>' ${geneID}.aa.align > ${geneID}.reference

python /mnt/home/stah3621/scratch/molevo/scripts/reorder_fasta.py ${geneID}.aa.align ${geneID}.reference > ${geneID}.aa.align.reordered

python /mnt/home/stah3621/scratch/molevo/scripts/reorder_fasta.py ${geneID}.nuc.fasta ${geneID}.reference > ${geneID}.nuc.fasta.reordered

tranalign -asequence ${geneID}.nuc.fasta.reordered -snucleotide -bsequence ${geneID}.aa.align.reordered -outseq ${geneID}.phylip -osformat3 phylipnon

done;
