#!/bin/bash

#SBATCH -o hostname_%j.out      # File to which STDOUT will be written
#SBATCH -e hostname_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL        # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=amandastahlke@gmail.com # Email to which notifications will be sent

# Calculate composite statistics from individual statistics based on
# ANGSD runs. Combines allele frequency change, and time series.
# This version uses ANGSD association analyses,
# analyses, spatpg on binary environment value, Mathieson and McVean's
# time series method, and allele frequency change.
#
# - FDR is calculated for adjusted p-values
# - GO term enrichment is performed on the full set of candidates,
# - Candidate SNPs are top 1% of composite stat
#
# INPUT
#
# process_rapture_2015_150609 -> rmdup_clone_filter_rapture_2015_150619 ->
# merge_rapture_2015_clone_filter_150619_150625 ->
# bowtie2_rapture_2015_allref_150904 -> aln_filter_2016-01-13 ->
#   AF change and variance: angsd_2016-01-15 -> afchange_2016-01-26
#   spatpg: angsd_2016-03-13 -> spatpg_2016-03-14a
#   M&M: angsd_2016-03-13 -> mm_2016-05-27
#
# OUTPUT
#
# Composite files for each method, population, and everything.
#
# NOTES
#
# There are two steps:
#
#   1) Calculate an adjusted p-value using either a genomic inflation
#      factor (for spatpg) or just based on the ranking
#      (for allele frequency change). 
#   2) Combine using Ma et al. (2015) DCMS statistic.
#
# I decided to use prevalence only as a covariate for spatpg.
# I rank for allele frequency change. For ANGSD, spatpg, and
# M&M, the p-values from the tests are included.
#
# The p-values and statistics are calculated for each SNP.
#
# SETTINGS
MIN_STATS=11            # Minimum number of statistics at each SNP 
MIN_STATS_AF=4          # Minimum number of AF statistics for AF composite 
MIN_STATS_TS=4          # Minimum number of time series statistics for time series 
MAX_P=0.99999           # Maximum p-value (change p = 1 to this for assoc. and spatpg)
SPATPG_MINP=0.00005     # Change spatpg p-values of 0 to this to let the DCMS algorithm work
GO_FDR=0.05             # False discovery rate for SNP2GO
#
# How to go from path2 to path1
relpath="""
import sys
from os.path import relpath
print relpath(sys.argv[1], sys.argv[2])
"""

# Get top x%
topxr="cargs=commandArgs(TRUE); tp = as.numeric(cargs[1]) / 100; "
topxr+="x = read.csv(file('stdin'), sep='\t', header=FALSE); "
topxr+="write.table(x[1:floor(nrow(x) * tp), ], file=stdout(), "
topxr+="sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE); "

PROJDIR="/mnt/lfs2/stah3621/devils/contemporary_sel"
GITDIR="/mnt/lfs2/stah3621/devils/contemporary_historical_sel_devils"
RUN="2019-2-22"
RESULTSDIR="${PROJDIR}/angsd_2019-01-18/next"
BRENDIR="/mnt/lfs2/bepstein/devilsrapture"
OUTDIR="${RESULTSDIR}/composite_stat/${RUN}"
RESULTDIR="${OUTDIR}/results"
SCRIPTDIR="${RESULTDIR}/script_copies"
LOGDIR="${RESULTDIR}/log"
WORKDIR="${RESULTDIR}/working"
AFDIR="${RESULTSDIR}/afchange/2019-01-21"
SPATDIR="${RESULTSDIR}/time_series/spatpg_2016-03-14a"
MMDIR="${RESULTSDIR}/time_series/mm_2019-02-20"
ANNFILE="${BRENDIR}/data/annotation/biomart_reduced.bed"
GTF="${GITDIR}/annotation/go/Sarcophilus_harrisii.DEVIL7.0.100.gtf.gz"
GO="${GITDIR}/annotation/go/goterms.tsv"
FAIFILE="${BRENDIR}/data/reference/sarHar1.fa.fai"

POPULATIONS="Fentonbury Forestier Freycinet Narawntapu West_Pencil_Pine"
POPULATIONS_ALL="Fentonbury Forestier Freycinet Mt_William Narawntapu West_Pencil_Pine"

gif="${PROJDIR}/brendan/script/genomic_inflation_correction.r"
q2p="${PROJDIR}/brendan/script/quantile_to_p.r"
compsnps="${PROJDIR}/brendan/script/css.r"
combiner="${PROJDIR}/brendan/script/combine_annotations.py"
snp2go="${GITDIR}/script/snp2go.r"

mkdir -p $RESULTDIR
mkdir -p $SCRIPTDIR
mkdir -p $LOGDIR
mkdir -p $WORKDIR
mkdir -p "${OUTDIR}/next"

## Make links to this scripts' output directory from the directories
## of all the precursor analyses. This is a little complicated, but I
## am trying to make it easy to follow the progress of the analysis
## by looking at the directory tree.
for d in $AFDIR $SPATDIR $MMDIR; do
    mkdir -p "${d}/next/composite_stat"
    cd "${d}/next/composite_stat"
    relname=$(python -c "$relpath" $OUTDIR ".")
    rm -f $RUN
    ln -s $relname $RUN
done

if [[ $SLURM_SUBMIT_DIR ]]; then

    module load R/3.4.3
    module load bedtools ## version?

    set -u
    set -e
    set -o pipefail

    cd $RESULTDIR

    #1 - Do adjustments

    # Get the list of allele frequency change files and convert the
    # quantiles to p-values
    affiles=""
    for pop in $POPULATIONS; do
        fname="${AFDIR}/results/${pop}/output.scores.tsv"
        fname="${AFDIR}/results/${pop}/output.scores.tsv"
        fname="${AFDIR}/results/${pop}/output.scores.tsv"
        fname="${AFDIR}/results/${pop}/output.scores.tsv"
        if [ -f $fname ]; then
            $q2p --output "${pop}.afchange.adjusted.tsv" \
                $fname \
                || { echo "adjusting ${fname} failed"; exit 1; }
            affiles+="${pop}.afchange.adjusted.tsv "
        fi
    done

    allfiles=$affiles
	echo $affiles
	echo $allfiles
    # Adjust spatpg

    tsfiles=""
    $gif --output "spatpg.adjusted.tsv" --double-p --min-p $SPATPG_MINP \
        --fdr --max-p $MAX_P \
        "${SPATDIR}/results/output.scores.tsv" \
        || { echo "adjusting spatpg failed"; exit 1; }
    allfiles+="spatpg.adjusted.tsv " 
    tsfiles+="spatpg.adjusted.tsv "
	echo $tsfiles
	echo $allfiles

    for pop in $POPULATIONS_ALL; do 
        fname="${MMDIR}/results/${pop}/output.scores.tsv"
        if [ -f $fname ]; then
            $gif --output "${pop}.mm.adjusted.tsv" \
                --fdr --max-p $MAX_P \
                $fname \
                || { echo "adjusting ${fname} failed"; exit 1; }
            tsfiles+="${pop}.mm.adjusted.tsv "
            allfiles+="${pop}.mm.adjusted.tsv "
		echo $tsfiles
		echo $allfiles
        fi
    done


    #######

    # 2 - Calculate composite score

    # Make a composite score just based on allele frequency change
   $compsnps --min-tests $MIN_STATS_AF --output "composite.snps.afchange" \
        $affiles \
        || { echo "combining AF change failed"; exit 1; }
		

    # Make a composite score just based on time series
    $compsnps --min-tests $MIN_STATS_TS --output "composite.snps.timeseries" \
        $tsfiles \djust spatpg
	|| { echo "combining time series failed"; exit 1; }

    # Make a composite score based on everything
    $compsnps --min-tests $MIN_STATS --output "composite.snps.everything" \
        $allfiles \
        || { echo "combining everything failed";  exit 1; }


    #######

    # 3 - Annotate based on mean composite score

    # Top composite scores
    # The 100.0 run is meant to get all genes within the specified
    # distance of a filter-passing SNP - might be useful for later
    # steps.
    for tp in "100.0" "1.0" "5.0" "1.0"; do
        for dist in 100000 200000; do
            q=$(awk 'BEGIN { print 1.0 - ('$tp' / 100.0) };')
            cd $RESULTDIR
            mkdir -p "annotation_top${tp}"
            cd "annotation_top${tp}"

            cut -f 1,2 $FAIFILE > "genome.txt"

            # All statistics
            tail -n +2 "../composite.snps.everything" | \
                sort -k 6rg,6 | cut -f 1-4 | \
                Rscript -e "$topxr" $tp | \
                bedtools sort > \
                "composite.snps.everything.top.bed" \
                || { echo "making bed file for all tests failed"; exit 1; }
            bedtools window -a "composite.snps.everything.top.bed" \
                -b $ANNFILE -w $dist > \
                "composite.snps.everything.top.genes.${dist}bp.txt" \
                || { echo "running bedtools for everything failed"; exit 1; }

            # Regions
            bedtools slop -b $dist -i "composite.snps.everything.top.bed" -g \
                "genome.txt" | \
                bedtools sort | \
                bedtools merge > "composite.snps.everything.top.merged.${dist}bp.bed" \
                || { echo "merging regions around top SNPs failed"; exit 1; }
            tail -n +2 "../composite.snps.everything" | cut -f 1-3,6 | \
                sort -k 4rg,4 | \
                Rscript -e "$topxr" $tp | \
                bedtools sort > \
                "composite.snps.everything.top.scores.bed" \
                || { echo "getting scores for top SNPs failed"; exit 1; }
            bedtools intersect -b "composite.snps.everything.top.scores.bed" \
                -a "composite.snps.everything.top.merged.${dist}bp.bed" -wb -wa | \
                sort -k 7rg,7 > \
                "composite.snps.everything.top.merged.scores.${dist}bp.bed" \
                || { echo "getting scores for regions failed"; exit 1; }

            # Af change only
            tail -n +2 "../composite.snps.afchange" | 
                sort -k 6rg,6 | cut -f 1-4 | \
                Rscript -e "$topxr" $tp | \
                bedtools sort > \
                "composite.snps.afchange.top.bed" \
                || { echo "making bed file for afchange tests failed"; exit 1; }
            bedtools window -a "composite.snps.afchange.top.bed" \
                -b $ANNFILE -w $dist > \
                "composite.snps.afchange.top.genes.${dist}bp.txt" \
                || { echo "running bedtools for af change failed"; exit 1; }
	
            # time series
            tail -n +2 "../composite.snps.timeseries" | \
                sort -k 6rg,6 | cut -f 1-4 | \
                Rscript -e "$topxr" $tp | \
                bedtools sort > \
                "composite.snps.timeseries.top.bed" \
                || { echo "making bed file for time series tests failed"; exit 1; }
            bedtools window -a "composite.snps.timeseries.top.bed" \
                -b $ANNFILE -w $dist > \
                "composite.snps.timeseries.top.genes.${dist}bp.txt" \
                || { echo "running bedtools for time series failed"; exit 1; }

            # GO term enrichment files
            if [ $tp != "100.0" ]; then
                # Make candidates file
                cut -f 1,3 "composite.snps.everything.top.bed" | \
                    sed 's/_random//g' | sed 's/chr.\{1,2\}_//g' | \
                    awk '{print $1 ".1\t" $2};' | sort > \
                    "composite.snps.everything.top.candidates" \
                    || { echo "making candidates file failed"; exit 1; }
                # Make non-candidates file
                cut -f 1,3 "../annotation_top100.0/composite.snps.everything.top.bed" | \
                    sed 's/_random//g' | sed 's/chr.\{1,2\}_//g' | \
                    awk '{print $1 ".1\t" $2};' | sort > \
                   "tmp" \
                    || { echo "making tmp non-candidates file failed"; exit 1; }
                comm -23 "tmp" "composite.snps.everything.top.candidates" > \
                    "composite.snps.everything.top.noncandidates" \
                    || { echo "making non-candidates file failed"; exit 1; }
                rm "tmp"
            fi

            $combiner --output "combined.genes.${dist}bp.tsv" \
                --bedtools \
                "composite.snps.afchange.top.genes.${dist}bp.txt" \
                "composite.snps.timeseries.top.genes.${dist}bp.txt" \
                || { echo "combining annotations failed"; exit 1; }
        done
    done


    # Parallel allele frequency changes
    # Find genes within 100 kb of SNPs in top 5% of at least 3 populations
    cd $RESULTDIR
    mkdir -p "parallel_top5.0"
    cd "parallel_top5.0"
    rcmd=''
    rcmd+='indir="'$AFDIR'/results"; '
    rcmd+='threshold=0.95; '
    rcmd+='populations = c("Fentonbury", "Forestier", "Freycinet", "Narawntapu", "West_Pencil_Pine"); '
    rcmd+='data = vector("list", length=length(populations)); '
    rcmd+='names(data) = populations; '
    rcmd+='topx = vector("list", length=length(populations)); '
    rcmd+='top100 = vector("list", length=length(populations)); '
    rcmd+='names(topx) = populations; '
    rcmd+='for(pop in populations) { '
    rcmd+='data[[pop]] = read.csv(file.path(indir, pop, "output.scores.tsv"), '
    rcmd+='sep="\t", header=FALSE, as.is=TRUE); '
    rcmd+='topx[[pop]] = data[[pop]][data[[pop]][, 8] >= threshold, 4]; '
    rcmd+='top100[[pop]] = data[[pop]][data[[pop]][, 8] >= 0, 4]; '
    rcmd+=' }; '
    rcmd+='counts = table(unlist(topx)); '
    rcmd+='counts_all = table(unlist(top100)); '
    rcmd+='top = names(counts)[counts >= 3]; ' 
    rcmd+='scaffold = gsub("(.+)-.+", "\\1", names(counts)[counts >= 3]); '
    rcmd+='pos = as.numeric(gsub(".+-", "", names(counts)[counts >= 3])); '
    rcmd+='result = cbind(scaffold, pos, pos + 1, names(counts)[counts >= 3]); '
    rcmd+='write.table(result, file="top.snps.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE); '
    rcmd+='scaffold = gsub("(.+)-.+", "\\1", names(counts_all)[counts_all >= 3 & !names(counts_all) %in% top]); '
    rcmd+='pos = as.numeric(gsub(".+-", "", names(counts_all)[counts_all >= 3 & !names(counts_all) %in% top])); '
    rcmd+='result = cbind(scaffold, pos, pos + 1, names(counts_all)[counts_all >= 3 & !names(counts_all) %in% top]); '
    rcmd+='write.table(result, file="nontop.snps.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE); '
    Rscript -e "$rcmd" \
        || { echo "filtering parallel snps failed"; exit 1; }
    bedtools window -a "top.snps.bed" \
        -b $ANNFILE -w 100000 > \
        "top.genes.100000bp.txt" \
        || { echo "running bedtools for parallel changes failed"; exit 1; }
    cut -f 1,3 "nontop.snps.bed" | \
        sed 's/_random//g' | sed 's/chr.\{1,2\}_//g' | \
        awk '{print $1 ".1\t" $2};' | sort > \
        "noncandidates.txt" \
        || { echo "making non-candidates file failed"; exit 1; }
    cut -f 1,3 "top.snps.bed" | \
        sed 's/_random//g' | sed 's/chr.\{1,2\}_//g' | \
        awk '{print $1 ".1\t" $2};' | sort > \
        "candidates.txt" \
        || { echo "making candidates file failed"; exit 1; }

    # M&M - p < 0.05 in three or more populations
    # Find genes within 100 kb of SNPs significant in at least 3 populations
    cd $RESULTDIR
    mkdir -p "mm_significant_3_pops"
    cd "mm_significant_3_pops"
    rcmd=''
    rcmd+='threshold=0.05; '
    rcmd+='populations = c("Fentonbury", "Forestier", "Freycinet", "Narawntapu", "West_Pencil_Pine"); '
    rcmd+='data = vector("list", length=length(populations)); '
    rcmd+='names(data) = populations; '
    rcmd+='topx = vector("list", length=length(populations)); '
    rcmd+='top100 = vector("list", length=length(populations)); '
    rcmd+='names(topx) = populations; '
    rcmd+='for(pop in populations) { '
    rcmd+='data[[pop]] = read.csv(paste0("../", pop, ".mm.adjusted.tsv"), '
    rcmd+='sep="\t", header=TRUE, as.is=TRUE); '
    rcmd+='topx[[pop]] = data[[pop]][data[[pop]][, "fdr"] < threshold, "snp"]; '
    rcmd+='top100[[pop]] = data[[pop]][data[[pop]][, "fdr"] >= 0, "snp"]; '
    rcmd+=' }; '
    rcmd+='counts = table(unlist(topx)); '
    rcmd+='counts_all = table(unlist(top100)); '
    rcmd+='top = names(counts)[counts >= 3]; ' 
    rcmd+='scaffold = gsub("(.+)-.+", "\\1", names(counts)[counts >= 3]); '
    rcmd+='pos = as.numeric(gsub(".+-", "", names(counts)[counts >= 3])); '
    rcmd+='result = cbind(scaffold, pos, pos + 1, names(counts)[counts >= 3]); '
    rcmd+='write.table(result, file="top.snps.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE); '
    rcmd+='scaffold = gsub("(.+)-.+", "\\1", names(counts_all)[counts_all >= 3 & !names(counts_all) %in% top]); '
    rcmd+='pos = as.numeric(gsub(".+-", "", names(counts_all)[counts_all >= 3 & !names(counts_all) %in% top])); '
    rcmd+='result = cbind(scaffold, pos, pos + 1, names(counts_all)[counts_all >= 3 & !names(counts_all) %in% top]); '
    rcmd+='write.table(result, file="nontop.snps.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE); '
    Rscript -e "$rcmd" \
        || { echo "filtering MM snps failed"; exit 1; }
    bedtools window -a "top.snps.bed" \
        -b $ANNFILE -w 100000 > \
        "top.genes.100000bp.txt" \
        || { echo "running bedtools for MM failed"; exit 1; }
    cut -f 1,3 "nontop.snps.bed" | \
        sed 's/_random//g' | sed 's/chr.\{1,2\}_//g' | \
        awk '{print $1 ".1\t" $2};' | sort > \
        "noncandidates.txt" \
        || { echo "making non-candidates file failed"; exit 1; }
    cut -f 1,3 "top.snps.bed" | \
        sed 's/_random//g' | sed 's/chr.\{1,2\}_//g' | \
        awk '{print $1 ".1\t" $2};' | sort > \
        "candidates.txt" \
        || { echo "making candidates file failed"; exit 1; }

    # spatpg FDR < 0.05
    # Find genes within 100 kb
    cd $RESULTDIR
    mkdir -p "spatpg_significant"
    cd "spatpg_significant"
    tail -n +2 "../spatpg.adjusted.tsv" | \
        awk '$10 < 0.05 { print $1 "\t" $2 "\t" $3 "\t" $4};' | \
        bedtools sort > "top.snps.bed" \
        || { echo "getting significant snps for spatpg failed"; exit 1; }
    tail -n +2 "../spatpg.adjusted.tsv" | \
        awk '$10 >= 0.05 { print $1 "\t" $2 "\t" $3 "\t" $4};' | \
        bedtools sort > "nontop.snps.bed" \
        || { echo "getting non-significant snps for spatpg failed"; exit 1; }
    bedtools window -a "top.snps.bed" \
        -b $ANNFILE -w 100000 > \
        "top.genes.100000bp.txt" \
        || { echo "running bedtools for parallel changes failed"; exit 1; }
    cut -f 1,3 "top.snps.bed" | \
        sed 's/_random//g' | sed 's/chr.\{1,2\}_//g' | \
        awk '{print $1 ".1\t" $2};' | sort > \
        "candidates.txt" \
        || { echo "making candidates file failed"; exit 1; }
    cut -f 1,3 "nontop.snps.bed" | \
        sed 's/_random//g' | sed 's/chr.\{1,2\}_//g' | \
        awk '{print $1 ".1\t" $2};' | sort > \
        "noncandidates.txt" \
        || { echo "making non-candidates file failed"; exit 1; }

    # GO term enrichment
   cd $RESULTDIR
   mkdir -p "go_enrichment"
    cd "go_enrichment"

    cat "../annotation_top1.0/composite.snps.everything.top.candidates" | \
    sort | uniq > "selection.composite.candidates.txt" \
        || { echo "catting composite candidates failed"; exit 1; }
    cat "../annotation_top1.0/composite.snps.everything.top.noncandidates"| \
        sort | uniq | \
        comm -23 - "selection.composite.candidates.txt" | \
        awk '{print $1 "\t" $2 "\t"}; ' | \
        sort | uniq > "selection.composite.noncandidates.txt" \
        || { echo "catting composite non-candidates failed"; exit 1;}

    $snp2go --candidates "selection.composite.candidates.txt" --non-candidates "selection.composite.noncandidates.txt" --gtf $GTF --go $GO   --output "top.composite.snps.go.100000bp" \
        || { echo "SNP2GO failed on composite ${dist}, ${tp}"; exit 1; }

else

   d=$(date '+%Y-%m-%d-%H%M')
   sfile="${SCRIPTDIR}/${d}-$(basename $0)"
   rm -f $sfile
   cp $0 $sfile
   cp $compsnps "${SCRIPTDIR}/${d}-$(basename $compsnps)"
   cp $combiner "${SCRIPTDIR}/${d}-$(basename $combiner)"
   cp $gif "${SCRIPTDIR}/${d}-$(basename $gif)"
   cp $q2p "${SCRIPTDIR}/${d}-$(basename $q2p)"
   cd $WORKDIR
   echo $WORKDIR
   
  sbatch --job-name="compstat" $sfile
fi
