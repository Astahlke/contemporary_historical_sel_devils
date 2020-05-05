#!/bin/bash
#
# Use ANGSD to estimate allele frequencies and a variety of summary
# statistics for each population and time point. This version is just
# like the 2016-01-15 version (as of 13 March 2016), except that it
# divides the data by year of birth.
#
# LAST UPDATED Nov 26, 2018 to run on slurm and incorporate additional sequencing
# 	which was prepared by Soraria. 
# INPUT
# /mnt/lfs2/stah3621/devils/contemporary_sel/6.0_All_samples_combined/
# Copied from Soraia: /mnt/lfs2/soraia/TD/5_bowtie2/6.0_All_samples_combined/
# In this case target loci have not been extracted
#
# Brendan's:
# process_rapture_2015_150609 -> rmdup_clone_filter_rapture_2015_150619 ->
# merge_rapture_2015_clone_filter_150619_150625 ->
# bowtie2_rapture_2015_allref_150904 -> aln_filter_2016-01-13
#
# Reads were aligned to the entire reference genome, then the target
# loci were extracted.
#
# OUTPUT
#
# The output is organized by population, time point, and program type.
# simple-run: genotype-likelihoods, saf files, and maf files for all
#   sites.
#
# SETTINGS
#
STARTYEAR=1999  # First year for multi-year datasets; I chose this b/c
                # it is the first year with data.
ENDYEAR=2014    # The last year with data
GENTIME=2       # Years per generation to use when producing
                # multi-year results
BAQ=2           # 2 = extended BAQ (cite Heng Li 2010 Bioinfo.); 0 = none; 1 = standard
NPROCS=1        # Number of processors for each ANGSD run
EPS="0.01"      # Initial guess at error rate
GL_MININDS=1    # Minimum individuals for genotype-likelihood steps
MINMAPQ=40      # Minimum mapping quality
MINBASEQ=25     # Minimum base call quality
AF_SNPP="1e-6"  # Maximum p-value for allele frequency step and genotype
                # likelihood step. I picked this value after a test run
                # with Freycinet 1999: values smaller than this did not
                # decrease the number of SNPs very much, and I wanted to
                # be able to detect a reasonable number of SNPs in smaller
                # samples.
#
# Note that I did not do any filtering of individuals aside from removing
# quolls and pouch young. I also do not specify any reads filters, but
# the input data are already filtered pretty stringently (see above).
#
# ORDER OF SUBPROGRAMS
#
#   prep, bamlist, genotype-likelihood, inbreeding,
#   folded-sfs, thetas, allele-frequency   
#
#   There is also a multisample-error program, but that is optional.
#   The folded-sfs and thetas programs are not necessary to get the
#   allele-frequency program to run. I had a snp-calls program
#   at one point, but then I decided to only run the genotype
#   likelihood step on sites that had a reasonable chance of being
#   a SNP.
#
#   This script does not calculate Fst. It is probably best to use
#   ngsPopGen to calculate Fst instead of following the ANGSD example -
#   something about not being able to use a joint SFS as a prior if
#   you have folded data.
#
# NOTES AND UPDATES
#
# One issue with this program is that it doesn't check to see
# if an individual survived to age two, it just checks whether the
# individual was born in the right year. I think most individuals
# that are in our dataset probably do survive to age two, but it is
# still not ideal. Year of birth is defined as being born on or after
# 01 Jan of the target year and before 01 Jan of the next year.
#
MYDIR="/mnt/lfs2/stah3621/devils/contemporary_sel"
PROJDIR="/mnt/lfs2/bepstein/devilsrapture"
PREVDIR="${PROJDIR}/results/aln_geno_analysis/alignments/bowtie2_rapture_2015_allref_150904/next/aln_filter/2016-01-13"
REFERENCE="${PROJDIR}/data/reference/sarHar1.fa"
OUTDIR="${MYDIR}/angsd_2019-01-18"
RESULTDIR="${OUTDIR}/results"
LOGDIR="${RESULTDIR}/log"
WORKDIR="${RESULTDIR}/working"
SCRIPTDIR="${RESULTDIR}/script_copies"
DATADIR="${RESULTDIR}/arraydata"
datafile="${DATADIR}/timeseries_data"
#POUCHYOUNG="${PROJDIR}/data/barcodes/rapture_2015_barcodes/pouch-young.txt"
#POPFILE="${PROJDIR}/data/barcodes/rapture_2015_barcodes/DevilsRaptureLocationsBarb3.tsv"
CAPTURELOCI="${MYDIR}/array_list_141220_table.tsv"
#SQLFILE="${MYDIR}/sql/Fundamentals_noRK.fixed.tsv"
DATABASEFILE="${MYDIR}/fundamentals_SH_20181114.txt"

#database='mysql -A --user=bepstein --host=data.ibest.uidaho.edu -pxbda19khb3 bepstein'

# getids year (at 2) pop
#function getids() {
#    d0=$(echo "${1} - 1" | bc)
#    d1=$(echo "${1} + 1" | bc)
#    yob="${d0}-12-31"
#    yob1="${d1}-01-01"
#    sql="""
#    SELECT Microchip, SiteOfFirstCapture FROM fundamentals
#    WHERE YOB > DATE('$yob') AND YOB < DATE('$yob1')
#    AND Species = 'Sarcophilus harrisii';
#    """
#    # Note that the tail -n +2 must come before grep or
#    # the first result will be lost.
#    echo $sql | $database | tail -n +2 | grep -i $2 | \
#        cut -f 1 | sort | uniq
#}

#beagle_converter="${PROJDIR}/script/angsd_beagle_like_to_vcf.py"

angsd_path="/mnt/lfs2/stah3621/bin/angsd"

if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0` TASK
        TASK options: 
		bamlist
		simple-run
	Run prep script first!"
  exit 0
fi

mkdir -p "${OUTDIR}/next"
mkdir -p $RESULTDIR
mkdir -p $LOGDIR
mkdir -p $WORKDIR
mkdir -p $SCRIPTDIR
mkdir -p $DATADIR

if [[ $SLURM_SUBMIT_DIR ]]; then
	#module load samtools/1.3.1
	module load R/3.3.2
	module load bedtools/2.26.0

    set pipefail -o
    set -e
    set -u

    TASK=`cut -f 1 -d " " ${datafile}.${SLURM_ARRAY_TASK_ID}`
    pop=`cut -f 2 -d " " ${datafile}.${SLURM_ARRAY_TASK_ID}`
    year=`cut -f 3 -d " " ${datafile}.${SLURM_ARRAY_TASK_ID}`
    chr=`cut -f 4 -d " " ${datafile}.${SLURM_ARRAY_TASK_ID}`

    case $pop in
        "Freycinet") dftdstart=2001 ;;
        "Narawntapu") dftdstart=2007 ;;
        "West_Pencil_Pine") dftdstart=2006 ;;
        "Forestier") dftdstart=2004 ;;
        "Mt_William") dftdstart=1996 ;;
        "Fentonbury") dftdstart=2005 ;;
        "Woolnorth") dftdstart=2020 ;; # I just set this to a time in the future
    esac


    regionsfile="${RESULTDIR}/regions.txt"
    popout="${RESULTDIR}/${pop}"
    yearslist="${popout}/years.txt"
    genslist="${popout}/generations.txt"
    popyearout="${RESULTDIR}/${pop}/${year}"
    bamlist="${popyearout}/bamlist.txt"
    masterlist="${RESULTDIR}/all_bams.txt"

    case $TASK in
 

        "bamlist")
        # This step makes directories for the population / year
        # combinations and also lists of bam files. I borrowed some
        # code from previous scripts. Individuals are counted in the
        # year in which they are two years old.
#        mkdir -p $popout
	
#	rm -f ${popyearout}/bamlist.txt
#        rm -f $genslist
        for firstyear in $(seq $STARTYEAR $GENTIME $ENDYEAR); do
            yearset=$(Rscript -e 'cat(seq('$firstyear', '$firstyear'+'$GENTIME'-1, 1), sep="-");')
            echo $yearset >> $genslist # all pops have the same genlist
#            popyearout="${popout}/${yearset}"
#            echo $pop $yearset >&2
#            mkdir -p $popyearout
#            rm -f "${popyearout}/tmp"
#            	for y in $(seq $firstyear 1 $(($firstyear+$GENTIME-1))); do
#		echo $y
#		grep -i $pop $DATABASEFILE | awk -v y="$y" '$6 ~ y' | cut -f 1 >> "${popyearout}/tmp" \
#		|| { echo "getting tmp list of samples for ${pop} ${y} in ${yearset} failed"; exit 1; }	
#            done
#      	     for samp in $(cat ${popyearout}/tmp); do
#		echo $samp
#		awk /$samp/ $masterlist >> "${popyearout}/bamlist.txt"
#		done 
            #rm -f "${popyearout}/tmp"
	done
        ;;


    "simple-run")
        # This steps estimates per-site allele frequency based
        # on genotype-likelihoods assuming HWE, as well as producing
        # genotype likelihood files and maf files for all sites.
        #
        # -doSaf 1: Calculate SAF assuming HWE (no inbreeding file);
        #   This is probably not ideal, but it gets around the inbreeding
        #   step. I need the saf files for calculating Fst.
        # -doMajorMinor 4 and -doMaf 2: I need these for allele frequency
        #   change and for the association step
        #
        # I also run the realSFS program here to get the overall
        # site frequency spectrum. A prior SFS can also be given,
        # but I don't do that here.
        #
        # Note: I found that -HWE_pval 1 caused freezing for some
        # of my test runs.

        subout="${popyearout}/simple-run"
        mkdir -p $subout
        cd $subout
        if [ $(wc -l $bamlist | cut -f 1 -d " ") -eq 0 ]; then
            echo "Skipping because no files - ${pop} ${year}"
            exit
        fi
	#for bam in $(cat $bamlist); do
	#	samtools index $bam
	#done

        # Experiment to see if I can get certain populations to not freeze
        $angsd_path/angsd -bam $bamlist \
            -rf $regionsfile \
            -minMapQ $MINMAPQ \
            -minQ $MINBASEQ \
            -ref $REFERENCE \
            -anc $REFERENCE \
            -baq $BAQ \
            -P $NPROCS \
            -nInd $(wc -l $bamlist | cut -f 1 -d " ") \
            -doMajorMinor 4 \
            -doMaf 2 \
            -SNP_pval 1 \
            -GL 2 \
            -doSaf 1 \
            -fold 1 \
            -out "output" \
            || { echo "estimating allele frequencies failed for ${pop} ${year}"; exit 1; }
        ;;

#    "beagle-prep")
        # Same as above, but only output SNPs, and actually produce the
        # genotype likelihood file. Then beagle is run on the results.

#        subout="${popyearout}/beagle"
#        mkdir -p $subout
#        cd $subout
#        if [ $(wc -l $bamlist | cut -f 1 -d " ") == 0 ]; then
#            echo "Skipping because no files - ${pop} ${year}"
#            exit
#        fi
#        $angsd_path/angsd -bam $bamlist \
#            -rf $regionsfile \
#            -minMapQ $MINMAPQ \
#            -minQ $MINBASEQ \
#            -ref $REFERENCE \
#            -anc $REFERENCE \
#            -baq $BAQ \
#            -P $NPROCS \
#            -nInd $(wc -l $bamlist | cut -f 1 -d " ") \
#            -doMajorMinor 4 \
#            -doMaf 2 \
#            -SNP_pval $AF_SNPP \
#            -minMaf 0.01 \
#            -minInd 10 \
#            -doGlf 2 \
#            -GL 2 \
#            -out "output" \
#            || { echo "making beagle input files failed for ${pop} ${year}"; exit 1; }
#
#        rm -f "output.mafs.gz"
#        for chr in "chr1" "chr2" "chr3" "chr4" "chr5" "chr6"; do
#            zcat "output.beagle.gz" | grep "marker\|^${chr}" | gzip > \
#                "output.${chr}.beagle.gz" \
#                || { echo "subsetting failed for ${pop} ${chr}"; exit 1; }
#        done
#
#        ;;

#    "beagle")
#        subout="${popyearout}/beagle"
#        cd $subout
#        $beagle_converter --output "output.${chr}.beagle.vcf" \
#            "output.${chr}.beagle.gz" \
#            || { echo "converting to vcf failed on ${chr} ${pop} ${year}"; exit 1; }
#        # Step 1: calculate genotype probabilities
#        beagle gl="output.${chr}.beagle.vcf" \
#            gprobs="true" \
#            out="output.unphased.${chr}" \
#            ne=50 \
#            nthreads=1 \
#            || { echo "beagle GP failed on ${chr} ${pop} ${year}"; exit 1; }
#        # Step 2: phase. Unfortunately, this doesn't produce any missing
#        # data (impute=false doesn't work unless you have a reference
#        # panel).
#        beagle gt="output.unphased.${chr}.vcf.gz" \
#            out="beagle.phased.${chr}" \
#            gprobs="true" \
#            ne=50 \
#            nthreads=1 \
#            || { echo "beagle phasing failed on ${chr} ${pop} ${year}"; exit 1; }
#        echo "Finished ${chr} ${pop} ${year}"
#        ;;
#
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

        "bamlist")
        for pop in "Freycinet" "Narawntapu" "West_Pencil_Pine" "Forestier" "Mt_William" "Fentonbury" "Woolnorth"; do
            echo "bamlist" $pop "_" > "${datafile}.${i}"
            i=$(($i+1))
        done
        ;;

        "simple-run")
        for pop in "Freycinet" "Narawntapu" "West_Pencil_Pine" "Forestier" "Mt_William" "Fentonbury" "Woolnorth"; do
	#for pop in "Woolnorth"; do
	   for year in `cat "${RESULTDIR}/${pop}/generations.txt"`; do
         #  for year in "1999-2000"; do  
	   if [[ $(wc -l "${RESULTDIR}/${pop}/${year}/bamlist.txt" | cut -f 1 -d " ") != 0 ]]; then
                    echo "simple-run" $pop $year > "${datafile}.${i}"
                    i=$(($i+1))
               fi
            done
        done
        ;;

#    "beagle-prep")
#        for pop in "Freycinet" "Narawntapu" "West_Pencil_Pine" "Forestier" "Mt_William" "Fentonbury" "Woolnorth"; do
#            for year in "before" "after"; do
#                if [ $(wc -l "${RESULTDIR}/${pop}/${year}/bamlist.txt" | cut -f 1 -d " ") != 0 ]; then
#                    echo "beagle-prep" $pop $year > "${datafile}.${i}"
#                    i=$(($i+1))
#                fi
#            done
#        done
#        ;;
#
#    "beagle")
#        for pop in "Freycinet" "Narawntapu" "West_Pencil_Pine" "Forestier" "Mt_William" "Fentonbury" "Woolnorth"; do
#            for year in "before" "after"; do
#                for chr in "chr1" "chr2" "chr3" "chr4" "chr5" "chr6"; do
#                    if [ $(wc -l "${RESULTDIR}/${pop}/${year}/bamlist.txt" | cut -f 1 -d " ") != 0 ]; then
#                        echo "beagle" $pop $year $chr > "${datafile}.${i}"
#                        i=$(($i+1))
#                    fi
#                done
#            done
#        done
#        ;;

    esac

    sbatch --mail-user=amandastahlke@gmail.com --mail-type=ALL --job-name="timeseries" --mem=120G --array=0-$(($i-1)) $sfile	
    echo "Submitted ${i} jobs"
fi
