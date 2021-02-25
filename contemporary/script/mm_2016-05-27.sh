#!/bin/bash
#
# Run Mathieson and McVean's (2013) time series selection method
# using the output of ANGSD to estimate allele counts.
#
# INPUT
#
# process_rapture_2015_150609 -> rmdup_clone_filter_rapture_2015_150619 ->
# merge_rapture_2015_clone_filter_150619_150625 ->
# bowtie2_rapture_2015_allref_150904 -> aln_filter_2016-01-13 ->
# angsd_2016-01-15 (simple-run)
#
# OUTPUT
#
# SETTINGS
#
# input - same as spatpg_2016-03-14a
MAX_P=0.000001  # Sites must be segregating with p <= MAX_P in at least
                # MIN_POPS datasets (each pop. / gen. combo. is a dataset)
MIN_MAC=3       # Sites must have rounded minor allele count > MIN_MAC in
                # at least MIN_POPS datasets
MIN_MAF=0.05    # Same as above, but a minor allele frequency in addition to count
MIN_POPS=5      # I added this parameter to cut down the number of SNPs
                # a bit. This was 10 in the previous version of the script.
GENERATIONS="1999-2000 2001-2002 2003-2004 2005-2006 2007-2008 2009-2010 2011-2012 2013-2014"
POPULATIONS="Fentonbury Forestier Freycinet Mt_William Narawntapu West_Pencil_Pine"
RUN="simple-run"
CHUNKSIZE=100
#
#
# Summary settings
WINDOW=100000   # Sliding window size
STEP=100000     # Sliding window step
#
# The first generation with disease:
#   Fentonbury: 2005-2006
#   Forestier:  2003-2004
#   Freycinet:  2001-2002
#   Mt William: 1999-2000
#   Narawntapu: 2007-2008
#   WPP:        2007-2008
#   Woolnorth:  never
#
# The Ne estimates from the full RAD dataset were used for Freycinet,
# Narawntapu, and West Pencil Pine. For the other populations, I had
# an NeEstimator run that didn't work, so I just set them to 35.
#
#   Freycinet   34
#   Narawntapu  37 
#   WPP         26
#
# NOTES
#
# The input for spatpg is prepared by taking the .mafs.gz file from
# ANGSD, multiplying the minor allele frequency * n inds * 2, and
# rounding to the nearest whole number to get the allele count.
#

extractcmd="""
import gzip
import os
import sys
sites = set()
ordered_sites = []
with open(sys.argv[1], 'rb') as handle:
    for line in handle:
        s = tuple(line.strip().split())
        sites.add(s); ordered_sites.append(s)
data = {}
if os.access(sys.argv[2], os.F_OK):
    with gzip.open(sys.argv[2], 'rb') as handle:
        handle.readline()
        for line in handle:
            fields = line.strip().split('\t')
            key = (fields[0], fields[1])
            if key not in sites:
                continue
            ref = round((1 -float(fields[6])) * int(fields[8]) * 2)
            alt = round(float(fields[6]) * int(fields[8]) * 2)
            data[key] = str(int(ref)) + ',' + str(int(alt))
for site in ordered_sites:
    if site in data:
        sys.stdout.write(data[site] + '\n')
    else:
        sys.stdout.write('0,0\n')
"""

PROJDIR="/mnt/lfs2/bepstein/devilsrapture"
PREVDIR="${PROJDIR}/results/aln_geno_analysis/alignments/bowtie2_rapture_2015_allref_150904/next/aln_filter/2016-01-13/next/af_est/angsd_2016-03-13"
OUTDIR="${PREVDIR}/next/time_series/mm_2016-05-27"
RESULTDIR="${OUTDIR}/results"
LOGDIR="${RESULTDIR}/log"
WORKDIR="${RESULTDIR}/working"
SCRIPTDIR="${RESULTDIR}/script_copies"
DATADIR="${RESULTDIR}/arraydata"
datafile="${DATADIR}/data"
FAIFILE="${PROJDIR}/data/reference/sarHar1.fa.fai"

summary="${PROJDIR}/script/summarize_mm.r"
stat="${PROJDIR}/script/window_quantile_stats.py"
slattice="/mnt/lfs2/bepstein/bin/slattice_ci.r"

snpfile="${RESULTDIR}/snps.txt"
genofile="geno.txt"

mkdir -p "${OUTDIR}/next"
mkdir -p $RESULTDIR
mkdir -p $LOGDIR
mkdir -p $WORKDIR
mkdir -p $SCRIPTDIR
mkdir -p $DATADIR


# Need a script to make the input files - should take a list of
# sites, a directory with the input files, and a list of generations
# in order. Should produce the spatpg output ready to use.


if [[ $PBS_O_WORKDIR ]]; then

    set -e
    set -u
    set -o pipefail

    task=`cut -f 1 -d " " ${datafile}.${PBS_ARRAYID}`
    pop=`cut -f 2 -d " " ${datafile}.${PBS_ARRAYID}`
    gen=`cut -f 3 -d " " ${datafile}.${PBS_ARRAYID}`
    chunk=`cut -f 4 -d " " ${datafile}.${PBS_ARRAYID}`

    case $pop in
        "Fentonbury") 
            first=2005
            ne=35
            ;;
        "Forestier") 
            first=2004
            ne=35
            ;;
        "Freycinet") 
            first=2001
            ne=34
            ;;
        "Mt_William") 
            first=1999  # I need to set this to the first year in the dataset, rather than the first year of infection here
            ne=35
            ;;
        "Narawntapu") 
            first=2007
            ne=37
            ;;
        "West_Pencil_Pine") 
            first=2006
            ne=26
            ;;
    esac

    case $task in

        "list-sites")
            # Get a list of segregating sites

            cd $RESULTDIR

            zcat "${PREVDIR}/results/${pop}/${gen}/${RUN}/output.mafs.gz" | \
                tail -n +2 | \
                awk '$8 <= '$MAX_P' && $7 >= '$MIN_MAF' && (1 - $7) >= '$MIN_MAF' && ($7 * $9 * 2) >= '$MIN_MAC' && ((1 - $7) * $9 * 2) >= '$MIN_MAC' { print $1 "\t" $2 };' > \
                "${snpfile}.${pop}.${gen}" \
                || { echo "Getting list of SNPs failed"; exit 1; }

            ;;

        "combine-sites")
            # combine the lists of SNPs into a single file

            cd $RESULTDIR

            sites_files=""
            for pop in $POPULATIONS; do
                for gen in $GENERATIONS; do
                    f="${snpfile}.${pop}.${gen}"
                    if [ -f $f ]; then
                        sites_files+="${f} "
                    fi
                done
            done

            cat $sites_files | sort -k 1b,1 -k 2n,2 | uniq -c | \
                awk '$1 >= '$MIN_POPS' { print $2 "\t" $3 };' > $snpfile \
                || { echo "combining sites failed"; exit 1; }

            ;;

        "extract-genotypes")
            # Extract the allele counts for each dataset.

            cd $RESULTDIR

            extracted_files=""
            python -c "$extractcmd" \
                $snpfile \
                "${PREVDIR}/results/${pop}/${gen}/${RUN}/output.mafs.gz" > \
                "extracted.${pop}.${gen}" \
                || { echo "Extracting data for ${pop} ${gen} failed"; exit 1; }

            ;;

        "input")

            mkdir -p "${RESULTDIR}/${pop}"
            cd "${RESULTDIR}/${pop}"
            
            rm -f $genofile
            extracted_files=""
            use="no"
            for gen in $GENERATIONS; do
                if [[ $(echo $gen | grep $first) ]]; then
                    use="yes"
                fi
                if [ $use = "yes" ]; then
                    f="../extracted.${pop}.${gen} "
                    extracted_files+="${f} "
                fi
            done

            paste $extracted_files > $genofile \
                || { echo "combining extracted files failed"; exit 1; }

            #rm $extracted_files

            nsnps=$(wc -l $snpfile | cut -f 1 -d " ") || \
                { echo "counting number of snps for ${pop} failed"; exit 1; }
            chunk=0
            mkdir -p "chunk${chunk}"
            cd "chunk${chunk}"
            for i in $(seq 1 $nsnps); do
                if [ $(($i%$CHUNKSIZE)) == 0 ]; then
                    chunk=$(($chunk+1))
                    mkdir -p "../chunk${chunk}"
                    cd "../chunk${chunk}" \
                        || { echo "switching to ${chunk} for ${pop} failed"; exit 1; }
                fi
                snpname=$(sed -n "${i}p" $snpfile | sed 's/\t/-/g') || \
                    { echo "getting SNP name for ${pop} ${i} from ${snpfile} failed"; exit 1; }
                outfile="${snpname}.tsv"
                echo "N	N.A" > $outfile
                sed -n "${i}p" "../${genofile}" | \
                    sed 's/\t/\n/g' | sed 's/,/\t/g' | \
                    awk '{print $1+$2 "\t" $2};' >> \
                    $outfile \
                    || { echo "making input: ${snpname} ${chunk} ${pop} failed"; exit 1; }
            done

            ;;

        "run")

            # Actually run the program

            cd "${RESULTDIR}/${pop}"

            $slattice --output "output.${chunk}"  --ne $ne \
                "chunk${chunk}" \
                || { echo "running slattice on ${pop} ${chunk} failed"; exit 1; }
            ;;

        "summary")
            # Summarize the data: p-value, qvalue, score, and quantile.

            cd "${RESULTDIR}/${pop}"

            # 1 - Put estimates from different chunks together
            nchunks=$(ls | grep "chunk.\+" | wc -l)
            rm -f "output.tsv"
            head -n 1 "output.0" > "output.tsv"
            for chunk in $(seq 0 $(($nchunks-1))); do
                tail -n +2 "output.${chunk}" >> "output.tsv" \
                    || { echo "concatenating ${chunk} snps failed"; exit 1; }
            done

            # 2 - Make the score files and the window quantile files
            $summary --output "output.scores.tsv" "output.tsv" $snpfile \
                || { echo "summarizing failed"; exit 1; }
            $stat --output "output.windows.tsv" \
                --window $WINDOW --step $STEP \
                "output.scores.tsv" \
                $FAIFILE \
                || { echo "calculating window stats failed"; exit 1; }

            ;;

    esac

    rm -f "${datafile}.${PBS_ARRAYID}"
else
    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-$(basename $0)-${1}"
    rm -f $sfile
    cp $0 $sfile
    cd $WORKDIR
    echo $WORKDIR
    
    i=0
    case $1 in
        "list-sites")
            q="tiny"
            for pop in $POPULATIONS; do
                for gen in $GENERATIONS; do
                    f="${PREVDIR}/results/${pop}/${gen}/${RUN}/output.mafs.gz "
                    if [ -f $f ]; then
                        echo "list-sites" $pop $gen > "${datafile}.${i}"
                        i=$(($i+1))
                    fi
                done
            done
            ;;
        "combine-sites")
            q="tiny"
            echo "combine-sites" > "${datafile}.${i}"
            i=$(($i+1))
            ;;
        "extract-genotypes")
            q="short"
            for pop in $POPULATIONS; do
                for gen in $GENERATIONS; do
                    f="${PREVDIR}/results/${pop}/${gen}/${RUN}/output.mafs.gz "
                    echo "extract-genotypes" $pop $gen > "${datafile}.${i}"
                    i=$(($i+1))
                done
            done
            ;;
        "input")
            q="short"
#            for pop in $POPULATIONS; do
            for pop in "Mt_William"; do
                echo "input" $pop > "${datafile}.${i}"
                i=$(($i+1))
            done
            ;;
        "run")
            q="long"
#            for pop in $POPULATIONS; do
            for pop in "Mt_William"; do
                nchunks=$(ls "${RESULTDIR}/${pop}" | grep -c "chunk")
                for j in $(seq 0 $(($nchunks-1))); do
                    echo "run" $pop "_" $j > "${datafile}.${i}"
                    i=$(($i+1))
                done
            done
            ;;
        "summary")
            q="tiny"
            for pop in $POPULATIONS; do
                echo "summary" $pop > "${datafile}.${i}"
                i=$(($i+1))
            done
            ;;
    esac

    if [ $q = "none" ]; then
        echo $sfile
        echo "run on the standalone server with PBS_ARRAYID=0"
    else
        qsub -q $q -t 0-$(($i-1))%150 $sfile
        echo "Submitted ${i} jobs"
    fi
fi
