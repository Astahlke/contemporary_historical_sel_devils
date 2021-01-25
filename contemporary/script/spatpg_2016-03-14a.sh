#!/bin/bash
#
# Run Zach Gompert's time series selection method using the output
# of ANGSD to estimate allele counts. This is a modified version of
# 2016-03-13 to test how setting the Ne prior effects estimates; in
# this case, the prior is 25-40. (Updated 08 April 2016 with a tighter
# prior.)
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
# input
MAX_P=0.000001  # Sites must be segregating with p <= MAX_P in at least
                # MIN_POPS datasets (each pop. / gen. combo. is a dataset)
MIN_MAC=3       # Sites must have rounded minor allele count > MIN_MAC in
                # at least MIN_POPS datasets
MIN_MAF=0.05    # Same as above, but a minor allele frequency in addition to count
MIN_POPS=5      # I added this parameter to cut down the number of SNPs
                # a bit. This was 10 in the previous version of the script.
GENERATIONS="1999-2000 2001-2002 2003-2004 2005-2006 2007-2008 2009-2010 2011-2012 2013-2014"
POPULATIONS="Fentonbury Forestier Freycinet Mt_William Narawntapu West_Pencil_Pine Woolnorth"
RUN="simple-run"
NCHUNKS=18
#
# spatpg
### ****** NOTE ********* ####
# Spatpg will crash if (NSTEPS - BURNIN) / THIN is not an
# integer - will try to do too many steps and go past the end
# of the matrix. I chose the Ne range to span the range I used
# for the full RAD data.
NSTEPS=100000       # number of mcmc steps (-n option)
BURNIN=10000        # burnin (-b option)
LOWNE=25            # lower bound on Ne (-l option)
UPPNE=40            # upper bound on Ne (-u option)
STDP=0.1            # std dev on prior coeff (-s option) - this is the default
STDNP=0.05          # std dev of proposal distr. (-p option) - this is the default
THINNING=10         # record 1 / $THINNING iterations
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
# NOTES
#
# The input for spatpg is prepared by taking the .mafs.gz file from
# ANGSD, multiplying the minor allele frequency * n inds * 2, and
# rounding to the nearest whole number to get the allele count.
#
# The thinning procedure looks at every pair of SNPs on a scaffold.
# If they are within WINSIZE kb, then only the one with the largest
# number of individuals is kept. If they have the same sample size,
# the first one is kept.
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

chunkcmd="""
import itertools
import random
import sys
genofile, snpfile, nchunks = sys.argv[1:4]
nchunks = int(nchunks)
split_data = [[] for i in range(nchunks)]
rand = random.Random()
chunks = range(nchunks)
with open(genofile, 'rb') as ghandle:
    with open(snpfile, 'rb') as shandle:
        for g, s in itertools.izip(ghandle, shandle):
            idx = rand.choice(chunks)
            split_data[idx].append((g, s))
for chunk in chunks:
    with open(genofile + '.' + str(chunk), 'wb') as gout:
        gout.write(sys.argv[4] + ' ' + sys.argv[5] + ' ' + str(len(split_data[chunk])) + '\n')
        with open(snpfile + '.' + str(chunk), 'wb') as sout:
            for g, s in split_data[chunk]:
                gout.write(g)
                sout.write(s)
"""

PROJDIR="/mnt/lfs2/bepstein/devilsrapture"
PREVDIR="${PROJDIR}/results/aln_geno_analysis/alignments/bowtie2_rapture_2015_allref_150904/next/aln_filter/2016-01-13/next/af_est/angsd_2016-03-13"
REFERENCE="${PROJDIR}/data/reference/sarHar1.fa"
OUTDIR="${PREVDIR}/next/time_series/spatpg_2016-03-14a"
RESULTDIR="${OUTDIR}/results"
LOGDIR="${RESULTDIR}/log"
WORKDIR="${RESULTDIR}/working"
SCRIPTDIR="${RESULTDIR}/script_copies"
DATADIR="${RESULTDIR}/arraydata"
datafile="${DATADIR}/data"
FAIFILE="${PROJDIR}/data/reference/sarHar1.fa.fai"

summary="${PROJDIR}/script/summarize_spatpg.r"
stat="${PROJDIR}/script/window_quantile_stats.py"

NGENS=$(echo $GENERATIONS | sed 's/ /\n/g' | wc -l)
NPOPS=$(echo $POPULATIONS | sed 's/ /\n/g' | wc -l)

snpfile="${RESULTDIR}/snps.txt"
poporder="${RESULTDIR}/poporder.txt"
genofile="${RESULTDIR}/geno.txt"
envfile="${RESULTDIR}/env.txt"
hdf5output="${RESULTDIR}/output.hdf5"

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

    module load hdf5/1.8.16
    module load spatpg/1.0

    set -e
    set -u
    set -o pipefail

    task=`cut -f 1 -d " " ${datafile}.${PBS_ARRAYID}`
    pop=`cut -f 2 -d " " ${datafile}.${PBS_ARRAYID}`
    gen=`cut -f 3 -d " " ${datafile}.${PBS_ARRAYID}`
    chunk=`cut -f 4 -d " " ${datafile}.${PBS_ARRAYID}`

    case $task in

        "list-sites")
            # Get a list of segregating sites for each dataset

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
            # Finally, put everything together

            cd $RESULTDIR
            
            # 1 - Make the genotype file
            rm -f $poporder
            rm -f $genofile
            rm -f $envfile

            extracted_files=""
            for pop in $POPULATIONS; do
                for gen in $GENERATIONS; do
                    echo $pop $gen >> $poporder
                    f="extracted.${pop}.${gen} "
                    extracted_files+="${f} "
                done
            done

            paste $extracted_files > $genofile \
                || { echo "combining extracted files failed"; exit 1; }

            # 2 - Make the environment file
            for pop in $POPULATIONS; do
                case $pop in
                    "Fentonbury") first=2005 ;;
                    "Forestier") first=2004 ;;
                    "Freycinet") first=2001 ;;
                    "Mt_William") first=1996 ;;
                    "Narawntapu") first=2007 ;;
                    "West_Pencil_Pine") first=2006 ;;
                    "Woolnorth") first=2020 ;;
                esac
                echo $GENERATIONS | sed 's/ /\n/g' | sed 's/-/\t/g' | \
                    awk '{ if($2 > '$first') { print "1" } else { print "-1" }};' >> $envfile \
                    || { echo "making env file for ${pop} failed"; exit 1; }
            done

            rm $extracted_files


            ;;

        "chunk")
            # Split input into chunks

            cd $RESULTDIR

            python -c "$chunkcmd" $genofile $snpfile $NCHUNKS $NPOPS $NGENS \
                || { echo "splitting into chunks failed"; exit 1; }

            ;;

        "run")

            # Actually run the program
            spatpg -g "${genofile}.${chunk}" \
                -e $envfile \
                -o "${hdf5output}.${chunk}" \
                -t $THINNING -n $NSTEPS -b $BURNIN \
                -l $LOWNE -u $UPPNE -s $STDP -p $STDNP \
                || { echo "spatpg failed on ${chunk}"; exit 1; }
            ;;

        "estpost")
            # This command is a little different - it doesn't actually run
            # anything - it just sets up a script to run on one of the
            # standalone servers.
            # For some reason estpost does not work when I run it in a job.
            # When run from an interactive session on one of the standalone servers,
            # it is fine.
            cd $RESULTDIR
            statistic="ne"
            for chunk in $(seq 0 $(($NCHUNKS-1))); do
                for statistic in "ne" "beta"; do
                    estpost -p $statistic -o \
                        "${RESULTDIR}/${statistic}.${chunk}.estimates.txt" \
                        -b 0 -h 20 -s 2 -w 1 \
                        "${hdf5output}.${chunk}" \
                        || { echo "estpost failed on ${statistic} ${chunk}"; exit 1; }
                    estpost -p $statistic -o \
                        "${RESULTDIR}/${statistic}.${chunk}.mcmcdiagnostics.txt" \
                        -b 0 -h 20 -s 3 -w 1 \
                        "${hdf5output}.${chunk}" \
                        || { echo "estpost failed on ${statistic} ${chunk}"; exit 1; }
                done
            done

            ;;

        "summary")
            # Summarize the data: p-value, qvalue, score, and quantile.
            # In this case the "score" will be the betas. The strength of
            # support is the proportion of iterations with beta > 0 (or
            # < 0 if the mean is negative, I think).
            cd $RESULTDIR

            # 1 - Put estimates from different chunks together
            #   - just use beta
            betaest="beta.estimates.txt"
            rm -f $betaest
            allsnps="snps.concatenated.txt"
            rm -f $allsnps
            for chunk in $(seq 0 $(($NCHUNKS-1))); do
                cat "beta.${chunk}.estimates.txt" >> $betaest \
                    || { echo "concatenating ${chunk} failed"; exit 1; }
                cat "snps.txt.${chunk}" >> $allsnps \
                    || { echo "concatenating ${chunk} snps failed"; exit 1; }
            done

            # 2 - Make the score files and the window quantile files
            $summary --output "output.scores.tsv" $betaest $allsnps \
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
            q="tiny"
            echo "input" > "${datafile}.${i}"
            i=$(($i+1))
            ;;
        "chunk")
            q="tiny"
            echo "chunk" > "${datafile}.${i}"
            i=$(($i+1))
            ;;
        "run")
            q="long"
            for j in $(seq 0 $(($NCHUNKS-1))); do
                echo "run" "_" "_" $j > "${datafile}.${i}"
                i=$(($i+1))
            done
            ;;
        "estpost")
            q="none"
            echo "estpost" > "${datafile}.${i}"
            i=$(($i+1))
            ;;
        "summary")
            q="tiny"
            echo "summary" > "${datafile}.${i}"
            i=$(($i+1))
            ;;
    esac

    if [ $q = "none" ]; then
        echo $sfile
        echo "run on the standalone server with PBS_ARRAYID=0"
    else
        qsub -q $q -t 0-$(($i-1)) $sfile
        echo "Submitted ${i} jobs"
    fi
fi
