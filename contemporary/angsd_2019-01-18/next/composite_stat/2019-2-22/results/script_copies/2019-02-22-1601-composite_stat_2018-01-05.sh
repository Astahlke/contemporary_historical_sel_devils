#!/bin/bash

#SBATCH -o hostname_%j.out      # File to which STDOUT will be written
#SBATCH -e hostname_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL        # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=amandastahlke@gmail.com # Email to which notifications will be sent

# # Modified by Amanda Stahlke Jan 14 to run on new data and optimize best

# Modified by Amanda Stahlke Sep 18 to exclude populations (Woolnorth and Mt William) that do not have both pre and post RAD data.

# Modified by me (Amanda Stahlke) Dec 14 2017 to run on my own directory, excluding GEMMA results 

# Original (from Brendan): 
# Calculate composite statistics from individual statistics based on
# ANGSD runs. Combines GWAS, allele frequency change, and time series.
# This version uses ANGSD association analyses, GEMMA association
# analyses, spatpg on binary environment value, Mathieson and McVean's
# time series method, and allele frequency change.
#
# Differences from 2016-06-22 script:
#
#   1 - FDR is calculated for adjusted p-values
#   2 - GEMMA p-values are calculated using just gamma rather than
#       beta X gamma; I tried using these values as p-values, but they
#       weren't really distributed right for that (min. around 0.9).
#   3 - GO term enrichment is performed on the full set of candidates,
#       not on each set separately
#   4 - The candidates are identified somewhat differently:
#       Top 1% of composite and
#       Top 5% allele frequency change in >= 3 pops and
#       FDR < 0.05 in M&M and
#       Top 1% in composite GWAS
#
# Differences from 2016-09-13 script:
#
#   1 - Uses new tumor growth rate code that calculates growth rate
#       across all tumors
#
# Differences from 2016-10-08 script
#
#   1 - Tumor growth rate calculated on log scale.
#
# INPUT
#
# process_rapture_2015_150609 -> rmdup_clone_filter_rapture_2015_150619 ->
# merge_rapture_2015_clone_filter_150619_150625 ->
# bowtie2_rapture_2015_allref_150904 -> aln_filter_2016-01-13 ->
#   GWAS: angsd_2016-01-15 -> gwas_2016-04-19_* and 2016-10-14
#   AF change and variance: angsd_2016-01-15 -> afchange_2016-01-26
#   spatpg: angsd_2016-03-13 -> spatpg_2016-03-14a
#   M&M: angsd_2016-03-13 -> mm_2016-05-27
#   GEMMA: angsd_2016-01-15 -> vcf_2016-04-25 -> various gemma runs
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
#      factor (for association and spatpg) or just based on the ranking
#      (for allele frequency change). Some of the association runs do
#      not really need much adjustment, but I adjust anyway.
#   2) Combine using Ma et al. (2015) DCMS statistic.
#
# I decided to use prevalence only as a covariate for spatpg.
#
# GEMMA is included - I ranked SNPs by the estimated size of the
# large effect, and use the ranks as pseudo-pvalues, similar to the
# way I use ranks for allele frequency change. For ANGSD, spatpg, and
# M&M, the p-values from the tests are included.
#
# The p-values and statistics are calculated for each SNP.
#
# SETTINGS
MIN_STATS=11            # Minimum number of statistics at each site # Brendan's Oct-2106  ; What should this be? 
#MIN_STATS_NOGEMMA=11    # Minimum number of statistics for non-gemma composite *Ultimately this is redundant ARS 01/04/18
MIN_STATS_AF=4          # Minimum number of AF statistics for AF composite # unchanged from Brendan's Oct-2016 script
#MIN_STATS_GWAS=9        # Minimum number of gwas statistics for GWAS composite
MIN_STATS_TS=4          # Minimum number of time series statistics for time series # unchanged from Brendan's Oct-2016 script
MAX_P=0.99999           # Maximum p-value (change p = 1 to this for assoc. and spatpg)
SPATPG_MINP=0.00005     # Change spatpg p-values of 0 to this to let the DCMS algorithm work
GO_FDR=0.05             # False discovery rate for SNP2GO
#
# To choose the optimal number of PCs, I just took the run with
# the lowest inflation factor. In some cases, it was very close, and
# probably unnecessary, but this seems objective.
#AGE_ADJUST="no_adjust"   # PCA adjustment for age at first infection association
#CC_ADJUST="no_adjust"    # PCA adjustment for case-control
#SURV_ADJUST="adjust_5"  # PCA adjustment for survival after disease
#RATE_ADJUST="no_adjust" # PCA adjustment for tumor growth rate
#
# GEMMA runs
#GEMMA_AGE_F="gemma_bslmm_2016-04-28f_age"
#GEMMA_AGE_M="gemma_bslmm_2016-04-28m_age"
#GEMMA_CC_F="gemma_bslmm_2016-05-03f_unc_cc"
#GEMMA_CC_M="gemma_bslmm_2016-05-03m_unc_cc"
#GEMMA_SURV_F="gemma_bslmm_2016-05-03f_unc_survival"
#GEMMA_SURV_M="gemma_bslmm_2016-05-03m_unc_survival"
#GEMMA_RATE_F="gemma_bslmm_2016-10-14f_unc_rate"
#GEMMA_RATE_M="gemma_bslmm_2016-10-14m_unc_rate"
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
RUN="2019-2-22"
RESULTSDIR="${PROJDIR}/angsd_2019-01-18/next"
#BREN_RESULTS_DIR="${PROJDIR}/brendan/angsd_2016-01-15/next"
BRENDIR="/mnt/lfs2/bepstein/devilsrapture"
OUTDIR="${RESULTSDIR}/composite_stat/${RUN}"
RESULTDIR="${OUTDIR}/results"
SCRIPTDIR="${RESULTDIR}/script_copies"
LOGDIR="${RESULTDIR}/log"
WORKDIR="${RESULTDIR}/working"
#GWAS_AGEDIR="${OUTDIR}/../../gwas/2016-04-19_age/results" ## modified ARS 05-19
#GWAS_SURVDIR="${OUTDIR}/../../gwas/2016-04-19_survival/results"
#GWAS_RATEDIR="${OUTDIR}/../../gwas/2016-10-14_rate/results"
# WAS_CCDIR="${OUTDIR}/../../gwas/2016-04-19_cc/results"
#GEMMADIR="${OUTDIR}/../../vcf/2016-04-25/next/gwas"
AFDIR="${RESULTSDIR}/afchange/2019-01-21"
SPATDIR="${RESULTSDIR}/time_series/spatpg_2016-03-14a"
MMDIR="${RESULTSDIR}/time_series/mm_2019-02-20"
ANNFILE="${BRENDIR}/data/annotation/biomart_reduced.bed"
GTF="${BRENDIR}/data/annotation/go_150413/Sarcophilus_harrisii.DEVIL7.0.79.gtf.gz"
GO="${BRENDIR}/data/annotation/go_150413/goterms.tsv"
FAIFILE="${BRENDIR}/data/reference/sarHar1.fa.fai"

POPULATIONS="Fentonbury Forestier Freycinet Narawntapu West_Pencil_Pine"
POPULATIONS_ALL="Fentonbury Forestier Freycinet Mt_William Narawntapu West_Pencil_Pine"

gif="${PROJDIR}/brendan/script/genomic_inflation_correction.r"
q2p="${PROJDIR}/brendan/script/quantile_to_p.r"
compsnps="${PROJDIR}/brendan/script/css.r"
combiner="${PROJDIR}/brendan/script/combine_annotations.py"
snp2go="${PROJDIR}/brendan/script/snp2go.r"
#gemmasum="${PROJDIR}/script/summarize_gemma_bslmm_gamma.r"

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

    module load R ## version? 
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

   # gwasfiles=""
    # Do the same thing to the gemma files, but also change the way
    # the p-values are calculated
#    $gemmasum --output "age_f.gemma.scores.tsv" \
#        "${GEMMADIR}/${GEMMA_AGE_F}/results/all/output/output.param.txt" \
#        || { echo "adjusting gemma failed"; exit 1; }
#    $q2p --output "age_f.gemma.adjusted.tsv" \
#        "age_f.gemma.scores.tsv" \
#        || { echo "adjusting gemma failed"; exit 1; }
#    allfiles+="age_f.gemma.adjusted.tsv "
#   gwasfiles+="age_f.gemma.adjusted.tsv "
#   $gemmasum --output "age_m.gemma.scores.tsv" \
#       "${GEMMADIR}/${GEMMA_AGE_M}/results/all/output/output.param.txt" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   $q2p --output "age_m.gemma.adjusted.tsv" \
#       "age_m.gemma.scores.tsv" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   allfiles+="age_m.gemma.adjusted.tsv "
#   gwasfiles+="age_m.gemma.adjusted.tsv "
#   $gemmasum --output "cc_f.gemma.scores.tsv" \
#       "${GEMMADIR}/${GEMMA_CC_F}/results/all/output/output.param.txt" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   $q2p --output "cc_f.gemma.adjusted.tsv" \
#       "cc_f.gemma.scores.tsv" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   allfiles+="cc_f.gemma.adjusted.tsv "
#   gwasfiles+="cc_f.gemma.adjusted.tsv "
#   $gemmasum --output "cc_m.gemma.scores.tsv" \
#       "${GEMMADIR}/${GEMMA_CC_M}/results/all/output/output.param.txt" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   $q2p --output "cc_m.gemma.adjusted.tsv" \
#       "cc_m.gemma.scores.tsv" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   allfiles+="cc_m.gemma.adjusted.tsv "
#   gwasfiles+="cc_m.gemma.adjusted.tsv "
#   $gemmasum --output "survival_f.gemma.scores.tsv" \
#       "${GEMMADIR}/${GEMMA_SURV_F}/results/all/output/output.param.txt" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   $q2p --output "survival_f.gemma.adjusted.tsv" \
#       "survival_f.gemma.scores.tsv" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   allfiles+="survival_f.gemma.adjusted.tsv "
#   gwasfiles+="survival_f.gemma.adjusted.tsv "
#   $gemmasum --output "survival_m.gemma.scores.tsv" \
#       "${GEMMADIR}/${GEMMA_SURV_M}/results/all/output/output.param.txt" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   $q2p --output "survival_m.gemma.adjusted.tsv" \
#       "survival_m.gemma.scores.tsv" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   allfiles+="survival_m.gemma.adjusted.tsv "
#   gwasfiles+="survival_m.gemma.adjusted.tsv "
#   $gemmasum --output "rate_f.gemma.scores.tsv" \
#       "${GEMMADIR}/${GEMMA_RATE_F}/results/all/output/output.param.txt" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   $q2p --output "rate_f.gemma.adjusted.tsv" \
#       "rate_f.gemma.scores.tsv" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   allfiles+="rate_f.gemma.adjusted.tsv "
#   gwasfiles+="rate_f.gemma.adjusted.tsv "
#   $gemmasum --output "rate_m.gemma.scores.tsv" \
#       "${GEMMADIR}/${GEMMA_RATE_M}/results/all/output/output.param.txt" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   $q2p --output "rate_m.gemma.adjusted.tsv" \
#       "rate_m.gemma.scores.tsv" \
#       || { echo "adjusting gemma failed"; exit 1; }
#   allfiles+="rate_m.gemma.adjusted.tsv "
#   gwasfiles+="rate_m.gemma.adjusted.tsv "

    # Adjust spatpg and association
#    $gif --output "gwas.age.adjusted.tsv" \
#        --fdr --max-p $MAX_P --fdr \
#        "${GWAS_AGEDIR}/${AGE_ADJUST}/output.scores.tsv" \
#        || { echo "adjusting association failed"; exit 1; }
#    allfiles+="gwas.age.adjusted.tsv "
#    allfiles_nogemma+="gwas.age.adjusted.tsv "
#    gwasfiles+="gwas.age.adjusted.tsv "
#    $gif --output "gwas.cc.adjusted.tsv" \
#        --fdr --max-p $MAX_P \
#        "${GWAS_CCDIR}/${CC_ADJUST}/output.scores.tsv" \
#        || { echo "adjusting association failed"; exit 1; }
#    allfiles+="gwas.cc.adjusted.tsv "
#    allfiles_nogemma+="gwas.cc.adjusted.tsv "
#    gwasfiles+="gwas.cc.adjusted.tsv "
#    $gif --output "gwas.survival.adjusted.tsv" \
#        --fdr --max-p $MAX_P \
#        "${GWAS_SURVDIR}/${SURV_ADJUST}/output.scores.tsv" \
#        || { echo "adjusting association failed"; exit 1; }
#    allfiles+="gwas.survival.adjusted.tsv "
#    allfiles_nogemma+="gwas.survival.adjusted.tsv "
#    gwasfiles+="gwas.survival.adjusted.tsv "
#    $gif --output "gwas.rate.adjusted.tsv" \
#        --fdr --max-p $MAX_P \
#        "${GWAS_RATEDIR}/${RATE_ADJUST}/output.scores.tsv" \
#        || { echo "adjusting association failed"; exit 1; }
#    allfiles+="gwas.rate.adjusted.tsv "
#    allfiles_nogemma+="gwas.rate.adjusted.tsv "
#    gwasfiles+="gwas.rate.adjusted.tsv "

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
        $tsfiles \
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

    # Combined GO term enrichment
    # I had originally hesitated to combine all the SNPs like this,
    # because there are slightly different lists of eligible SNPs for
    # each analysis, but I think it isn't too harmful.
    # - Parallel change SNPs
    # - MM FDR < 0.05 SNPs
    # - Top 1% GWAS SNPs
    # - Top 1% composite SNPs
    cd $RESULTDIR
   mkdir -p "go_enrichment"
    cd "go_enrichment"
#"../annotation_top1.0/composite.snps.everything.top.noncandidates" \
#"../annotation_top1.0/composite.snps.everything.top.candidates" | \
    cat "../mm_significant_3_pops/candidates.txt" \
        "../spatpg_significant/candidates.txt" \
        "../parallel_top5.0/candidates.txt"  | \
    sort | uniq > "selection.candidates.txt" \
        || { echo "combining selection candidates failed"; exit 1; }
    cat "../mm_significant_3_pops/noncandidates.txt" \
        "../spatpg_significant/noncandidates.txt" \
        "../parallel_top5.0/noncandidates.txt" | \
        sort | uniq | \
        comm -23 - "selection.candidates.txt" | \
        awk '{print $1 "\t" $2 "\t"}; ' | \
        sort | uniq > "selection.noncandidates.txt" \
        || { echo "combining selection non-candidates failed"; exit 1; }

        cat "../annotation_top1.0/composite.snps.everything.top.candidates" | \
    sort | uniq > "selection.composite.candidates.txt" \
        || { echo "catting composite candidates failed"; exit 1; }
    cat "../annotation_top1.0/composite.snps.everything.top.noncandidates" \
        sort | uniq | \
        comm -23 - "selection.composite.candidates.txt" | \
        awk '{print $1 "\t" $2 "\t"}; ' | \
        sort | uniq > "selection.composite.noncandidates.txt" \
        || { echo "catting composite non-candidates failed"; }   # ; exit 1;

    $snp2go --candidates "selection.candidates.txt" --non-candidates "selection.noncandidates.txt" --gtf $GTF --go $GO   --output "top.snps.go.100000bp" \
        || { echo "SNP2GO failed on ${dist}, ${tp}"; exit 1; }
    $snp2go --candidates "selection.composite.candidates.txt" --non-candidates "selection.composite.noncandidates.txt" --gtf $GTF --go $GO   --output "top.composite.snps.go.100000bp" \
        || { echo "SNP2GO failed on composite ${dist}, ${tp}"; exit 1; }

else
    #sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-$(basename $0)-${1}"
    #rm -f $sfile
    #cp $0 $sfile
    #cd $WORKDIR
    #echo $WORKDIR

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
   
 #qsub -q reg $sfile 
  sbatch --job-name="compstat" $sfile
fi
