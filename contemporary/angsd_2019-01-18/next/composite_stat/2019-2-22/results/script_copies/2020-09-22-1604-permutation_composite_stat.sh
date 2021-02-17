#!/bin/bash

#SBATCH -o hostname_%j.out      # File to which STDOUT will be written
#SBATCH -e hostname_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL        # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=amandastahlke@gmail.com # Email to which notifications will be sent


## Author: Amanda Stahlke 
## Date: Sep 22, 2020 
## Purpose: Run permutation test to examine sensitivity of DCMS

# INPUT
#
# process_rapture_2015_150609 -> rmdup_clone_filter_rapture_2015_150619 ->
# merge_rapture_2015_clone_filter_150619_150625 ->
# bowtie2_rapture_2015_allref_150904 -> aln_filter_2016-01-13 ->
#   GWAS: angsd_2016-01-15 -> gwas_2016-04-19_* and 2016-10-14
#   AF change and variance: angsd_2016-01-15 -> afchange_2016-01-26
#   spatpg: angsd_2016-03-13 -> spatpg_2016-03-14a
#   M&M: angsd_2016-03-13 -> mm_2016-05-27

# The p-values and statistics are calculated for each SNP.
#
# SETTINGS
MIN_STATS=11            # Minimum number of statistics at each site 
MIN_STATS_AF=4          # Minimum number of AF statistics for AF composite # unchanged from Brendan's Oct-2016 script
MIN_STATS_TS=4          # Minimum number of time series statistics for time series # unchanged from Brendan's Oct-2016 script
MAX_P=0.99999           # Maximum p-value (change p = 1 to this for assoc. and spatpg)
SPATPG_MINP=0.00005     # Change spatpg p-values of 0 to this to let the DCMS algorithm work


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
BRENDIR="/mnt/lfs2/bepstein/devilsrapture"
OUTDIR="${RESULTSDIR}/composite_stat/${RUN}"
RESULTDIR="${OUTDIR}/results"
SCRIPTDIR="${RESULTDIR}/script_copies"
LOGDIR="${RESULTDIR}/log"
WORKDIR="${RESULTDIR}/working"
AFDIR="${RESULTSDIR}/afchange/2019-01-21"
SPATDIR="${RESULTSDIR}/time_series/spatpg_2016-03-14a"
MMDIR="${RESULTSDIR}/time_series/mm_2019-02-20"

POPULATIONS="Fentonbury Forestier Freycinet Narawntapu West_Pencil_Pine"
POPULATIONS_ALL="Fentonbury Forestier Freycinet Mt_William Narawntapu West_Pencil_Pine"

compsnps="${PROJDIR}/brendan/script/css.r"


if [[ $SLURM_SUBMIT_DIR ]]; then

    module load R/3.4.3
    module load bedtools ## version?

    set -u
    set -e
    set -o pipefail

    cd $RESULTDIR


# Get the list of allele frequency change files and convert the
    # quantiles to p-values
    affiles=""
    for pop in $POPULATIONS; do ## no Mt. William
        fname="${AFDIR}/results/${pop}/output.scores.tsv"
        fname="${AFDIR}/results/${pop}/output.scores.tsv"
        fname="${AFDIR}/results/${pop}/output.scores.tsv"
        fname="${AFDIR}/results/${pop}/output.scores.tsv"
        if [ -f $fname ]; then    
	affiles+="${pop}.afchange.adjusted.tsv "
        fi
    done

    allfiles=$affiles
        echo "aff files are  $affiles"
        echo "all files are $allfiles"

tsfiles=""
    allfiles+="spatpg.adjusted.tsv "
    tsfiles+="spatpg.adjusted.tsv "

    for pop in $POPULATIONS_ALL; do
            tsfiles+="${pop}.mm.adjusted.tsv "
            allfiles+="${pop}.mm.adjusted.tsv "
	done

 	echo "time series files are  $tsfiles"
	echo " all input for composite state are $allfiles"

#######

    # 2 - Calculate composite score

    # Make a composite score based on everything
    $compsnps --min-tests $MIN_STATS --output "composite.snps.everything" \
        $allfiles \
        || { echo "combining everything failed";  exit 1; }


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

  sbatch --job-name="compstat_permutation" $sfile
fi


