#!/bin/bash
#
# SETTINGS
#
# BE: I also do not specify any reads filters, but
# the input data are already filtered pretty stringently (see above).
#
# ORDER OF SUBPROGRAMS: prep
#
MYDIR="/mnt/lfs2/stah3621/devils/contemporary_sel"
PROJDIR="/mnt/lfs2/bepstein/devilsrapture"
REFERENCE="${PROJDIR}/data/reference/sarHar1.fa"
OUTDIR="${MYDIR}/angsd_2019-01-18"
RESULTDIR="${OUTDIR}/results"
WORKDIR="${RESULTDIR}/working"
SCRIPTDIR="${RESULTDIR}/script_copies"
DATADIR="${RESULTDIR}/arraydata"
datafile="${DATADIR}/data"
CAPTURELOCI="${MYDIR}/array_list_141220_table.tsv"
DATABASEFILE="${MYDIR}/fundamentals_SH_20181114.txt"


if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0` prep
	all this script does is prep"
  exit 0
fi


mkdir -p $OUTDIR
mkdir -p "${OUTDIR}/next"
mkdir -p $RESULTDIR
mkdir -p $WORKDIR
mkdir -p $SCRIPTDIR
mkdir -p $DATADIR

if [[ $SLURM_SUBMIT_DIR ]]; then
	module load bedtools/2.26.0

    set pipefail -o
    set -e
    set -u

    TASK=`cut -f 1 -d " " ${datafile}.${SLURM_ARRAY_TASK_ID}`

    regionsfile="${RESULTDIR}/regions.txt"
    masterlist="${RESULTDIR}/all_bams.txt"

    case $TASK in
 

        "prep")
        # This step does the filtering and input file preparation that
        # can be done on all individuals at once.

        # Make a list of all the bam files
        ls ${MYDIR}/6.0_All_samples_combined/*bam | sort | uniq > \
	    $masterlist \
            || { echo "making master list failed"; exit 1; }

        # Create a file specifying the regions to analyze. Only
        # non-X scaffolds and only within 600 bp of the target
        # cut site. ANGSD appears to use open 1-based intervals for
        # this file. It is important to merge overlapping intervals
        # before printing out the file - otherwise ANGSD will print
        # sites in the overlap regions twice.
        grep -v "chrX" $CAPTURELOCI | tail -n +2 | \
            awk '$3 == "+" {print $8 "\t" $4-1 "\t" $4+600 "\txxx"} $3 == "-" {print $8 "\t" $4-601 "\t" $4 "\txxx"};' | \
            awk 'BEGIN {OFS="\t"}; $2 < 0 {print $1 "\t0\t" $3 "\t" $4}; $2 >= 0 {print};' | \
            bedtools sort | bedtools merge | \
            awk '{print $1 ":" $2+1 "-" $3};' > $regionsfile \
            || { echo "making regions file failed"; exit 1; }

        ;;


    esac

    rm -f "${datafile}.${SLURM_ARRAY_TASK_ID}"

else
    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-$(basename $0)-${1}"
    rm -f $sfile
    cp $0 $sfile
    cd $WORKDIR
    echo $WORKDIR
    
    i=0
    case $1 in
        "prep")
        echo "prep" "_" "_" > "${datafile}.${i}"
        i=$(($i+1))
        ;;
    esac

    sbatch --array=0-$(($i-1)) $sfile	
    echo "Submitted ${i} jobs"
fi
