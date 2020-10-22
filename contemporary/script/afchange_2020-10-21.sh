#!/bin/bash -l
#
# Calculate allele frequency change between before and after populations
# for every SNP and arc-sin transfrom a la Kelley and Hughes 2019. Also calculate variance in allele frequency.
#
# INPUT
#
# process_rapture_2015_150609 -> rmdup_clone_filter_rapture_2015_150619 ->
# merge_rapture_2015_clone_filter_150619_150625 ->
# bowtie2_rapture_2015_allref_150904 -> aln_filter_2016-01-13 ->
# angsd_2016-01-15 (simple-run)
#
# SETTINGS
#
MAX_P=0.000001  # Only use sites that are segregating with p <= MAXP
                # in at least one of the time points. (Applied in each
                # population separately.)
MIN_INDS=10     # Only use sites with at least this many individuals
MIN_MAF=0.05    # Only use sites with MAF >= MIN_MAF in at least one time point
                # in both time points. (Separately for each population.)
POPULATIONS="Freycinet Forestier Fentonbury Narawntapu West_Pencil_Pine"
WINDOW=100000   # Sliding window size
STEP=100000     # Sliding window step
#
# VARIANCE
#
# To calculate variance in allele frequency, the sites that pass the
# above filters in all five populations are used. Then the ratio of
# variance among populations before to after is reported.
#

afchange="""
import gzip
import sys
import math

sitesfile = sys.argv[1]
maf0 = sys.argv[2]
maf1 = sys.argv[3]
sites = set()
with open(sitesfile, 'rb') as handle:
    for line in handle:
        sites.add(tuple(line.strip().split('\t')))
af0 = {}
with gzip.open(maf0, 'rb') as handle:
    for line in handle:
        fields = line.strip().split('\t')
        key = (fields[0], fields[1])
        if key not in sites:
            continue
        af = float(fields[6])
        af0[key] = af
af1 = {}
with gzip.open(maf1, 'rb') as handle:
    for line in handle:
        fields = line.strip().split('\t')
        key = (fields[0], fields[1])
        if key not in sites:
            continue
        af = float(fields[6])
        af1[key] = af
for key in sorted(sites, key=lambda x: (x[0], int(x[1]))):
    sys.stdout.write('\t'.join(key) + '\t' + str((2*math.asin(math.sqrt(af1[key]))) - (2*math.asin(math.sqrt(af0[key])))) + '\n')
"""

varcalc="""
import collections
import gzip
import sys
def var(x):
    m = sum(x) / float(len(x))
    v = sum((s - m) ** 2 for s in x)
    return float(v) / (len(x) - 1)
sitesfile = sys.argv[1]
sites = set()
with open(sitesfile, 'rb') as handle:
    for line in handle:
        sites.add(tuple(line.strip().split('\t')))
afs = collections.defaultdict(list)
for fname in sys.argv[2:]:
    with gzip.open(fname, 'rb') as handle:
        for line in handle:
            fields = line.strip().split('\t')
            key = (fields[0], fields[1])
            if key not in sites:
                continue
            af = float(fields[6])
            afs[key].append(af)
for key in sorted(sites, key=lambda x: (x[0], int(x[1]))):
    sys.stdout.write('\t'.join(key) + '\t' + str(var(afs[key])) + '\n')
"""
MYDIR="/mnt/lfs2/stah3621/devils/contemporary_sel"
PROJDIR="/mnt/lfs2/bepstein/devilsrapture"
#OUTPROJDIR="/mnt/lfs2/bepstein/devilsrapture"
PREVDIR="${MYDIR}/angsd_2019-01-18"
OUTDIR="${PREVDIR}/next/afchange/2019-01-21"
RESULTDIR="${OUTDIR}/results"
LOGDIR="${RESULTDIR}/log"
WORKDIR="${RESULTDIR}/working"
SCRIPTDIR="${RESULTDIR}/script_copies"
DATADIR="${RESULTDIR}/arraydata"
FAIFILE="${PROJDIR}/data/reference/sarHar1.fa.fai"
datafile="${DATADIR}/data"

afsummary="${MYDIR}/script/summarize_value_abs.r"
varsummary="${MYDIR}/script/summarize_value_raw.r"
stat="${MYDIR}/script/window_quantile_stats.py"

angsd_path="/mnt/lfs2/stah3621/bin/angsd"
module load R

mkdir -p "${OUTDIR}/next"
mkdir -p $RESULTDIR
mkdir -p $LOGDIR
mkdir -p $WORKDIR
mkdir -p $SCRIPTDIR
mkdir -p $DATADIR

if [[ $SLURM_SUBMIT_DIR ]]; then

    #module load angsd/0.910

    set pipefail -o
    set -u
    set -e

    task=`cut -f 1 -d " " ${datafile}.${SLURM_ARRAY_TASK_ID}`
    pop=`cut -f 2 -d " " ${datafile}.${SLURM_ARRAY_TASK_ID}`

    case $task in

        "run")
            # I mostly followed the instructions here:
            # http://popgen.dk/angsd/index.php/Fst

            indirb="${PREVDIR}/results/${pop}/before/simple-run"
            indira="${PREVDIR}/results/${pop}/after/simple-run"

            suboutdir="${RESULTDIR}/${pop}"
            mkdir -p $suboutdir
            cd $suboutdir

            # Get a set of usable sites
            zcat "${indira}/output.mafs.gz" "${indirb}/output.mafs.gz" | \
                tail -n +2 | \
                awk '$8 <= '$MAX_P' && $7 >= '$MIN_MAF' && (1 - $7) >= '$MIN_MAF' {print $1 "\t" $2 };' | \
                sort | uniq > "p.sites" \
                || { echo "filtering by p-value failed"; exit 1; }
            zcat "${indira}/output.mafs.gz" | tail -n +2 | \
                awk '$9 >= '$MIN_INDS' { print $1 "\t" $2 }; ' | \
                sort > "after_inds.sites" \
                || { echo "filtering after by inds failed"; exit 1; }
            zcat "${indirb}/output.mafs.gz" | tail -n +2 | \
                awk '$9 >= '$MIN_INDS' { print $1 "\t" $2 }; ' | \
                sort > "before_inds.sites" \
                || { echo "filtering before by inds failed"; exit 1; }
            comm -12 "after_inds.sites" "before_inds.sites" > \
                "both_inds.sites" \
                || { echo "filtering by inds failed"; exit 1; }
            comm -12 "both_inds.sites" "p.sites" | sort -k 1b,1 -k 2n,2 | \
                uniq > "sites.tsv" \
                || { echo "combining filtered files failed"; exit 1; }
            rm "p.sites" "after_inds.sites" "before_inds.sites"
	;;

        "calc-afchange")    
	    indirb="${PREVDIR}/results/${pop}/before/simple-run"
            indira="${PREVDIR}/results/${pop}/after/simple-run"

            suboutdir="${RESULTDIR}/${pop}"
            mkdir -p $suboutdir
            cd $suboutdir	

	# Calculate the change in allele frequency at the usable sites
            python -c "$afchange" "sites.tsv" "${indirb}/output.mafs.gz" \
                "${indira}/output.mafs.gz" | sort -k 1b,1 -k 2n,2 > \
                "output.afchange.tsv"
            ;;

        "variance")
            suboutdir="${RESULTDIR}/variance"
            mkdir -p $suboutdir
            cd $suboutdir

            # Get common sites
            npops=$(echo $POPULATIONS | sed 's/ /\n/g' | wc -l)
            popsites=""
            infilesb=""
            infilesa=""
            for pop in $POPULATIONS; do
                popsites+="../${pop}/sites.tsv "
                infilesb+="${PREVDIR}/results/${pop}/before/simple-run/output.mafs.gz "
                infilesa+="${PREVDIR}/results/${pop}/after/simple-run/output.mafs.gz "
            done
            cat $popsites | sort -k 1b,1 -k 2n,2 | uniq -c | \
                awk '$1 == '$npops' { print $2 "\t" $3 };' > \
                "sites.tsv" \
                || { echo "Getting list of sites failed"; exit 1; }

            # Calculate variance for these sites
            python -c "$varcalc" "sites.tsv" $infilesb > "var.before.tsv" \
                || { echo "calculating before var failed"; exit 1; }
            python -c "$varcalc" "sites.tsv" $infilesa > "var.after.tsv" \
                || { echo "calculating after var failed"; exit 1; }

            # Calculate before / after
            paste "var.before.tsv" "var.after.tsv" | \
                awk '{ if($6 > 0) { print $1 "\t" $2 "\t" $3 / $6 } else { print $1 "\t" $2 "\tInf" } };' > \
                "output.var.tsv" \
                || { echo "calculating ratio failed"; exit 1; }
            ;;

        "summary-afchange")
            suboutdir="${RESULTDIR}/${pop}"
            cd $suboutdir
	echo "running afsumary pm $pop"
            $afsummary --output "output.scores.tsv" "output.afchange.tsv" \
                || { echo "summarizing ${pop} failed"; exit 1; }
            echo "running ${stat}"
            $stat --output "output.windows.tsv" \
                --window $WINDOW --step $STEP \
                "output.scores.tsv" \
                $FAIFILE \
                || { echo "calculating window stats failed for ${pop}"; exit 1; }
            echo "finished ${stat}"
            echo "output should be in ${suboutdir}"
            ;;

        "summary-variance")
            suboutdir="${RESULTDIR}/variance"
            cd $suboutdir

            $varsummary --output "output.scores.tsv" "output.var.tsv" \
                || { echo "summarizing variance failed"; exit 1; }
            $stat --output "output.windows.tsv" \
                --window $WINDOW --step $STEP \
                "output.scores.tsv" \
                $FAIFILE \
                || { echo "calculating window stats failed for ${pop}"; exit 1; }
            ;;
    esac

    rm -f "${datafile}.${SLURM_SUBMIT_DIR}"

else

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-$(basename $0)-${1}"
    rm -f $sfile
    cp $0 $sfile
    cd $WORKDIR
    echo $WORKDIR
    
    i=0
    case $1 in
        "run")
        for pop in $POPULATIONS; do
            echo "run" $pop > "${datafile}.${i}"
            i=$(($i+1))
        done
        ;;

"calc-afchange")
        for pop in $POPULATIONS; do
            echo "calc-afchange" $pop > "${datafile}.${i}"
            i=$(($i+1))
        done
        ;;


    "summary-afchange")
        for pop in $POPULATIONS; do
            echo "summary-afchange" $pop > "${datafile}.${i}"
            i=$(($i+1))
        done
        ;;

    "variance")
        echo "variance" > "${datafile}.${i}"
        i=$(($i+1))
        ;;

    "summary-variance")
        echo "summary-variance" > "${datafile}.${i}"
        i=$(($i+1))
        ;;

    esac
	sbatch --job-name="afchange" --mem=120G --array=0-$(($i-1)) $sfile
    #qsub -q short -t 0-$(($i-1)) $sfile
    echo "Submitted ${i} jobs"

fi
