#!/usr/bin/env bash

dir="/mnt/home/stah3621/scratch/devils/molevo/addkoala/candidates"
#mkdir $diri

#cut -f 4 ${dir}/candidate_orthologtable_*.tsv | cut -f 1 -d '/' > ${dir}/koala_candidate_orthologs.txt
#cut -f 4 ${dir}/wholegenome_orthologtable_110717.tsv | cut -f 1 -d '/' > ${dir}/koala_allg_orthologs.txt

#cat ${dir}/koala_allg_orthologs.txt

#while read line; do grep $line ${dir}/../transcriptsPC.fasta.out -A1 > ${dir}/Koala_Transcripts/${line}.nuc.fasta ; done < ${dir}/koala_candidate_orthologs.txt

echo "working director is" pwd 

module load muscle
module load emboss

mkdir ${dir}/getorf

for f in ${dir}/Koala_Transcripts/Locus*nuc.fasta; do
echo $f

koalageneID=$(echo $f | cut -f 11 -d '/' | cut -f1 -d.)
ID=$(echo $koalageneID | cut -f 2 -d '_' | cut -b 1-3)
echo $koalageneID

##getorf for geneID
#getorf -sequence $f -outseq ${dir}/getorf/${koalageneID}.getorf.out

## grab the longest ORF
#python /mnt/home/stah3621/scratch/devils/molevo/addkoala/parse_getorf.py $dir/getorf/${koalageneID}.getorf.out > $dir/getorf/${koalageneID}.getorf.aa.fasta
#ls $dir/getorf/${koalageneID}.getorf.aa.fasta

devilgeneID=$(grep ${koalageneID} ${dir}/koala_devil_*_orthologtable_*.tsv | cut -f 9) # wholegenome_orthologtable_110717.tsv
echo $devilgeneID 

ls ${dir}/ensembl_seq/${devilgeneID}.nuc*
#head ${dir}/../wholegenome/${devilgeneID}.nuc*
ls ${dir}/ensembl_seq/${devilgeneID}.aa*
#head ${dir}/../wholegenome/${devilgeneID}.aa*

#cat ${dir}/ensembl_seq/${devilgeneID}.aa.fasta >  ${dir}/${devilgeneID}.aa.fasta
#cat $dir/getorf/${koalageneID}.getorf.aa.fasta >>  ${dir}/${devilgeneID}.aa.fasta
#cat ${dir}/ensembl_seq/${devilgeneID}.nuc.fasta > ${dir}/${devilgeneID}.nuc.fasta
#echo '>PCI'$ID >> ${dir}/${devilgeneID}.nuc.fasta
#tail --lines=+2  $dir/Koala_Transcripts/${koalageneID}.nuc.fasta >> ${dir}/${devilgeneID}.nuc.fasta

done

mv ENSSHAG000000* ${dir}/nucaafasta/

### copied from dataprep.sh
mkdir ${dir}/aaalign
mkdir ${dir}/reference
mkdir ${dir}/reordered

for f in ${dir}/nucaafasta/ENSSHAG*aa.fasta; do
geneID=$(echo $f | cut -f 11 -d '/' | cut -f1 -d.)

echo $geneID
muscle -in $f -out ${dir}/aaalign/${geneID}.aa.align

grep '>' ${dir}/aaalign/${geneID}.aa.align > ${dir}/reference/${geneID}.reference

python /mnt/home/stah3621/scratch/devils/molevo/scripts/reorder_fasta.py ${dir}/aaalign/${geneID}.aa.align ${dir}/reference/${geneID}.reference > ${dir}/reordered/${geneID}.aa.align.reordered

python /mnt/home/stah3621/scratch/devils/molevo/scripts/reorder_fasta.py ${dir}/nucaafasta/${geneID}.nuc.fasta ${dir}/reference/${geneID}.reference > ${dir}/reordered/${geneID}.nuc.fasta.reordered

tranalign -asequence ${dir}/reordered/${geneID}.nuc.fasta.reordered -snucleotide -bsequence ${dir}/reordered/${geneID}.aa.align.reordered -outseq ${dir}/${geneID}.phylip -osformat3 phylipnon

done

## file lists generated in R
mkdir ${dir}/allfour
mkdir ${dir}/three
mkdir ${dir}/oppwall
mkdir ${dir}/oppkoala
mkdir ${dir}/wallkoala
mkdir ${dir}/two

# R source('prepareorthologs.R') for ach of these groups 
allfour=4
three=3
two=2

for f in ${dir}/ENSSHAG*.phylip; do
geneID=$(echo $f | cut -f 10 -d '/' | cut -f1 -d.)

_seqcount=` head -n 1 ${dir}/${geneID}.phylip | awk '{print $1}'`
if [[ $_seqcount -eq $allfour ]];
then
        mv ${dir}/${geneID}.phylip ${dir}/allfour/
        
else
                echo "seqcount not 4"
        fi
if [[ $_seqcount -eq $two ]];
then
        mv ${dir}/${geneID}.phylip ${dir}/two/
        
	fi

_seqcount=` head -n 1 ${dir}/${geneID}.phylip | awk '{print $1}'`
if [[ $_seqcount -eq $three ]];
then
        mv ${dir}/${geneID}.phylip ${dir}/three/

else
                echo "seqcount not 4 or 3"
        fi
        done


for f in ${dir}/three/ENSSHAG*.phylip; do
geneID=$(echo $f | cut -f 11 -d '/' | cut -f1 -d.)

if ( grep -q "MEU" $f ) && ( grep -q "MOD" $f );
  then
	echo "MEU & MOD $f"
	mv $f $dir/oppwall/
fi

if ( grep -q "MEU" $f ) && ( grep -q "PCI" $f );
  then
        echo "MEU & PCI $f"
        mv $f $dir/wallkoala/
fi

if ( grep -q "MOD" $f ) && ( grep -q "PCI" $f );
  then
        echo "MOD & PCI $f"
        mv $f $dir/oppkoala/
fi

done
