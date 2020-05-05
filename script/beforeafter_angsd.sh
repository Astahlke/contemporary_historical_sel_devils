#!/bin/bash
#
# Use ANGSD to estimate allele frequencies and a variety of summary
# statistics for each population and time point.
#
# LAST UPDATED Nov 26, 2018 to run on slurm and incorporate additional sequencing
#       which was prepared by Soraria.
#
# INPUT
#
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
# genotype-likelihood:  initial genotype likelihood and allele frequency
#   estimate (glf.gz and mafs.gz, respectively) for site that have a
#   good chance of being SNPs.
# inbreeding: output.indF file has the inbreeding coefficients
# folded-sfs: output.sfs has the site frequency spectrum (folded)
#   after taking into account inbreeding
# thetas: population genetic stats
# allele-frequency: same as genotype-likelihood, but all sites, not
#   just SNPs, and also has a file called output.hwe.gz that has the
#   HWE p-values.
# simple-run: genotype-likelihoods, saf files, and maf files for all
#   sites.
#
# SETTINGS
#
STARTYEAR=1999  # First year for multi-year datasets; I chose this b/c
                # it is the first year with data.
ENDYEAR=2016    # The last year with data
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
# I dropped the "all" category (all years) from the results b/c I
# realized I wasn't using it and it was taking up a lot of space and
# time. I also originally ran this without a regions file, but then
# I decided it would save space and time to use only sites within 600 bp
# of the anticipated cut site. Finally, I originally had the genotype
# likelihoods in plain text, but switched to binary form after seeing
# how big the files were.
#
# I made a simple-run program to get what I needed for the downstream
# steps faster than if I went through the whole inbreeding calculation.
#

MYDIR="/mnt/lfs2/stah3621/devils/contemporary_sel"
PROJDIR="/mnt/lfs2/bepstein/devilsrapture"
#OUTPROJDIR="/mnt/lfs2/bepstein/devilsrapture"
#PREVDIR="${PROJDIR}/results/aln_geno_analysis/alignments/bowtie2_rapture_2015_allref_150904/next/aln_filter/2016-01-13"
REFERENCE="${PROJDIR}/data/reference/sarHar1.fa" # could it matter which reference we align to?
OUTDIR="${MYDIR}/angsd_2019-01-18"
RESULTDIR="${OUTDIR}/results"
LOGDIR="${RESULTDIR}/log"
WORKDIR="${RESULTDIR}/working"
SCRIPTDIR="${RESULTDIR}/script_copies"
DATADIR="${RESULTDIR}/arraydata"
datafile="${DATADIR}/beforeafter_data"
#POUCHYOUNG="${PROJDIR}/data/barcodes/rapture_2015_barcodes/pouch-young.txt"
DATEFILE="${PROJDIR}/data/barcodes/rapture_2015_barcodes/sampling_info-all.tsv"
CAPTURELOCI="${PROJDIR}/../devils/annotation/array_list_141220_table.tsv"
#SQLFILE="${MYDIR}/sql/Fundamentals_noRK.fixed.tsv"
DATABASEFILE="${MYDIR}/fundamentals_SH_20181114.txt"

#beagle_converter="${OUTPROJDIR}/script/angsd_beagle_like_to_vcf.py"
angsd_path="/mnt/lfs2/stah3621/bin/angsd"

if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0` TASK
        TASK options:
                bamlist
                simplerun
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

        module load samtools/1.3.1
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


 #   quollfile="${RESULTDIR}/quolls.txt"
 #   pyfile="${MYDIR}/pouchyoung.txt" # removed "sample" from Brendan's PY file
 #   excludefile="${RESULTDIR}/exclude.txt"
    regionsfile="${RESULTDIR}/regions.txt"
    popout="${RESULTDIR}/${pop}"
    yearslist="${popout}/years.txt"
    genslist="${popout}/generations.txt"
    popyearout="${RESULTDIR}/${pop}/${year}"
    bamlist="${popyearout}/bamlist.txt"
    masterlist="${RESULTDIR}/all_bams.txt"

    case $TASK in
 

        "prep")
        # This step does the filtering and input file preparation that
        # can be done on all individuals at once.

        # First use the database to generate a list of quolls
        #echo "select microchip from fundamentals where fundamentals.Species != 'Sarcophilus harrisii';" | \
        #    mysql -A --user=bepstein --host=data.ibest.uidaho.edu -pxbda19khb3 bepstein | \
        #    tail -n +2 | awk '{print "sample_" $1};' > \
        #    $quollfile \
        #    || { echo "getting list of quolls failed"; exit 1; }

        # Then combine the quolls and pouch young into a single file
        #awk '{print $2};' $POUCHYOUNG > $pyfile
        #cat $pyfile $quollfile | sort -k 1b,1 > $excludefile

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


        "bamlist")
        # This step makes directories for the population / year
        # combinations and also lists of bam files. I borrowed some
        # code from previous scripts. Repeat individuals are counted
        # in all years they were sampled.

	mkdir -p $popout
        awk 'BEGIN { FS="\t" }; {print $2 "-" $3};' $DATEFILE | \
            grep -i "$pop" | sort | uniq | grep -v "\-$" | \
            sed 's/.\+-//g' | grep -v "0000" > $yearslist || \
            { echo "Getting list of years failed for ${pop}"; exit 1; }
        year_files=""
        for year in `cat $yearslist`; do
            popyearout="${popout}/${year}"
            echo $pop $year >&2
            mkdir -p $popyearout
                echo $year
                grep -i $pop $DATABASEFILE | awk -v y="$year" '$6 ~ y' | cut -f 1 >> "${popyearout}/tmp"
	#grep -P "${pop}\t${year}" $DATEFILE | \
            #    awk '{print "sample_" $1};' | \
            #    sort -k 1b,1 | uniq > "${popyearout}/tmp" \
                #|| { echo "getting tmp list of samples for ${pop} ${year} failed"; exit 1; }
            	#done 
		#sort -k 1b,1 "${popyearout}/tmp" | comm -23 - $excludefile > "${popyearout}/tmp2" \
                #|| { echo "removing excluded samples from ${pop} ${year} failed"; exit 1; }
        	for samp in $(cat ${popyearout}/tmp); do
                echo $samp
                awk /$samp/ $masterlist >> "${popyearout}/bamlist.txt"
		done
            #    "${popyearout}/tmp2" > \
            #    "${popyearout}/bamlist.txt" \
                #|| { echo "making list of bam files for ${pop} ${year} failed"; exit 1; }
            year_files+="${popyearout}/bamlist.txt "
            rm -f "${popyearout}/tmp" #"${popyearout}/tmp2"
        done
        # Multiple years
        #rm -f $genslist
        #for firstyear in $(seq $STARTYEAR $GENTIME $ENDYEAR); do
        #    yearset=$(Rscript -e 'cat(seq('$firstyear', '$firstyear'+'$GENTIME'-1, 1), sep="-");')
        #    echo $yearset >> $genslist
        #    popyearout="${popout}/${yearset}"
        #    echo $pop $yearset >&2
        #    mkdir -p $popyearout
        #    rm -f "${popyearout}/tmp"
        #    for y in $(seq $firstyear 1 $(($firstyear+$GENTIME-1))); do
        #        grep -P "${pop}\t${y}" $DATEFILE | \
        #            awk '{print "sample_" $1};' | \
        #            sort -k 1b,1 | uniq >> "${popyearout}/tmp" \
        #            || { echo "getting tmp list of samples for ${pop} ${y} in ${yearset} failed"; exit 1; }
        #    done
        #    sort -k 1b,1 "${popyearout}/tmp" | \
        #        comm -23 - $excludefile > "${popyearout}/tmp2" \
        #        || { echo "removing excluded samples from ${pop} ${yearset} failed"; exit 1; }
        #    awk '{print "'$PREVDIR'/results/merged/" $1 ".bam"};' \
        #        "${popyearout}/tmp2" > \
        #        "${popyearout}/bamlist.txt" \
        #        || { echo "making list of bam files for ${pop} ${yearset} failed"; exit 1; }
        #    rm -f "${popyearout}/tmp" "${popyearout}/tmp2"
        #done
        # All years
        mkdir -p "${popout}/all"
        cat $year_files | sort | uniq > "${popout}/all/bamlist.txt" \
            || { echo "making list for all years failed"; exit 1; }

        # Before
        mkdir -p "${popout}/before"
        before_year_files=$(echo $year_files | sed 's/ /\n/g' | \
            awk 'gensub(".+/(.+)/bamlist.txt", "\\1", $1) <= '$dftdstart' {print $1};')
        cat $before_year_files | sort | uniq > "${popout}/before/bamlist.txt" \
            || { echo "making list for before years failed"; exit 1; }
        # After
        mkdir -p "${popout}/after"
        after_year_files=$(echo $year_files | sed 's/ /\n/g' | \
            awk 'gensub(".+/(.+)/bamlist.txt", "\\1", $1) > '$dftdstart' {print $1};')
        cat $after_year_files | sort | uniq > "${popout}/after/bamlist.txt" \
            || { echo "making list for after years failed"; exit 1; }

        ;;


#        "error-multisample")
#        # This step estimates error rates using the MAF spectrum
#        # See http://popgen.dk/angsd/index.php/Error_estimation
#        # I don't really need this step, but I included it b/c it
#        # might be useful later.
#        subout="${popyearout}/error-multisample"
#        mkdir -p $subout
#        cd $subout
#        echo $pop $year >&2
#        if [ $(wc -l $bamlist | cut -f 1 -d " ") == 0 ]; then
#            echo "Skipping because no files - ${pop} ${year}"
#            exit
#        fi
#        angsd -bam $bamlist \
#            -rf $regionsfile \
#            -minMapQ $MINMAPQ \
#            -minQ $MINBASEQ \
#            -ref $REFERENCE \
#            -fai "${REFERENCE}.fai" \
#            -baq $BAQ \
#            -doCounts 1 \
#            -out "output" \
#            -doError 1 \
#            -eps $EPS \
#            -doMajorMinor 2 \
#            -nThreads $NPROCS \
#            -minSites 1000 \
#            -minPhat $(Rscript -e "cat(1.0/$(wc -l $bamlist | cut -f 1 -d " "))") \
#            || { echo "calculating error rate for ${pop} ${year} failed"; exit 1; }
#        ;;
#
#
#        "genotype-likelihood")
#        # This step calculates genotype likelihoods using the GATK
#        # method and the samtools method. Note that this is not the
#        # latest GATK method.
#        # See https://github.com/fgvieira/ngsF/tree/master/examples .
#        #
#        # -doMajorMinor 4 = major allele is reference allele
#        # -doMaf 2 = assume known major allele, unknown minor
#        #
#        # Produces output-gatk.glf.gz and output-gatk.mafs.gz
#        # The mafs file could actually be used for allele frequencies,
#        # but the estimates can be improved by taking into account
#        # inbreeding (maybe).
#        #
#        # The mafs.gz column unknownEM = minor allele frequency est.
#        # and pu-EM = p-value (for whether this site is a SNP)
#        subout="${popyearout}/genotype-likelihood"
#        mkdir -p $subout
#        cd $subout
#        echo $pop $year >&2
#        if [ $(wc -l $bamlist | cut -f 1 -d " ") == 0 ]; then
#            echo "Skipping because no files - ${pop} ${year}"
#            exit
#        fi
#        angsd -bam $bamlist \
#            -rf $regionsfile \
#            -minMapQ $MINMAPQ \
#            -minQ $MINBASEQ \
#            -ref $REFERENCE \
#            -baq $BAQ \
#            -doMajorMinor 4 \
#            -doMaf 2 \
#            -nInd $(wc -l $bamlist | cut -f 1 -d " ") \
#            -SNP_pval $AF_SNPP \
#            -GL 2 \
#            -doGlf 3 \
#            -out "output-gatk" \
#            -minind $GL_MININDS \
#            -P $NPROCS \
#            || { echo "calculating GL-gatk for ${pop} ${year} failed"; exit 1; }
#        # I had planned to do samtools genotype likelihoods too, but
#        # I decided not to to save space.
#
#        ;;
#
#        "inbreeding")
#        # This step estimates per-individual inbreeding coefficients
#        # that can be used for better SFS estimation. The program it
#        # uses is not actually part of ANGSD, but is recommended by
#        # the ANGSD wiki.
#        # The first run is fast and is meant to get an intial estimate.
#        # The second run uses the initial estimate as a starting point
#        # and has more stringent convergence criteria. The ngsF manual
#        # suggests 1e-9 (for epsilon) for low coverage data, but I found
#        # this was taking a very long time, so I set it to 1e-5.
#        # The two step process was in the examples that came
#        # with ngsF, although I made the second step use the --approx_EM
#        # settings so that it would run faster.
#        # The ngsF example README also says to use -islog 1 with the output
#        # from ANGSD, but the example doesn't, nor does the script that
#        # comes with angsd-wrapper, so I think -islog 1 is meant to be
#        # used when running a different ngsPopGen program.
#        #
#        # I increased the number of iterations because some run were
#        # ending without hitting the convergence criterion.
#        subout="${popyearout}/inbreeding"
#        mkdir -p $subout
#        cd $subout
#        echo $pop $year >&2
#        if [ $(wc -l $bamlist | cut -f 1 -d " ") == 0 ]; then
#            echo "Skipping because no files - ${pop} ${year}"
#            exit
#        fi
#        zcat "../genotype-likelihood/output-gatk.glf.gz" | ngsF \
#            --n_ind $(wc -l $bamlist | cut -f 1 -d " ") \
#            --n_sites $((`zcat "../genotype-likelihood/output-gatk.mafs.gz" | wc -l`-1)) \
#            --init_values "r" \
#            --min_epsilon "0.001" \
#            --approx_EM \
#            --glf - \
#            --n_threads $NPROCS \
#            --out "initial_estimate" \
#            || { echo "calculating inbreeding first estimate for ${pop} ${year} failed"; exit 1; }
#
#        zcat "../genotype-likelihood/output-gatk.glf.gz" | ngsF \
#            --n_ind $(wc -l $bamlist | cut -f 1 -d " ") \
#            --n_sites $((`zcat "../genotype-likelihood/output-gatk.mafs.gz" | wc -l`-1)) \
#            --init_values "initial_estimate.pars" \
#            --min_epsilon "1e-5" \
#            --max_iters 4000 \
#            --approx_EM \
#            --glf - \
#            --n_threads $NPROCS \
#            --out "output" \
#            || { echo "calculating inbreeding refined estimate for ${pop} ${year} failed"; exit 1; }
#
#        ;;
#
#
#
#        "folded-sfs")
#        # This step calculates the folded site frequency spectrum for
#        # each population / year combination.
#        # -doSaf 2 calculates the site allele frequency likelihood using
#        # an inbreeding model.
#        subout="${popyearout}/folded-sfs"
#        mkdir -p $subout
#        cd $subout
#        echo $pop $year >&2
#        angsd -bam $bamlist \
#            -rf $regionsfile \
#            -minMapQ $MINMAPQ \
#            -minQ $MINBASEQ \
#            -ref $REFERENCE \
#            -baq $BAQ \
#            -doSaf 2 \
#            -indF "../inbreeding/output" \
#            -anc $REFERENCE \
#            -glf "../genotype-likelihood/output-gatk.glf.gz" \
#            -fold 1 \
#            -P $NPROCS \
#            -out "output" \
#            || { echo "angsd making folded SAF likelihood failed"; exit 1; }
#        # TODO!!! I think this needs a count of the number of individuals
#        realSFS "output.saf.idx" \
#            -underFlowProtect 1 \
#            -P $NPROCS > \
#            "output.sfs" \
#            || { echo "calculating ML estimate of SFS failed"; exit 1; }
#
#        ;;
#
#
#        "thetas")
#        # This step estimates thetaW and Tajima's D.
#        subout="${popyearout}/thetas"
#        mkdir -p $subout
#        cd $subout
#        echo $pop $year >&2
#        if [ $(wc -l $bamlist | cut -f 1 -d " ") == 0 ]; then
#            echo "Skipping because no files - ${pop} ${year}"
#            exit
#        fi
#        angsd -bam $bamlist \
#            -rf $regionsfile \
#            -minMapQ $MINMAPQ \
#            -minQ $MINBASEQ \
#            -ref $REFERENCE \
#            -baq $BAQ \
#            -P $NPROCS \
#            -out "output" \
#            -doThetas 1 \
#            -doSaf 2 \
#            -indF "../inbreeding/output" \
#            -pest "../folded-sfs/output.sfs" \
#            -anc $REFERENCE \
#            -glf "../genotype-likelihood/output-gatk.glf.gz" \
#            -fold 1 \
#            || { echo "estimating thetas failed"; exit 1; }
#        thetaStat make_bed "output.thetas.gz" # This may seg fault, but the next step checks for errors
#        thetaStat validate_bed "output.thetas.gz" \
#            || { echo "binary thetas file did not pass"; exit 1; }
#        thetaStat do_stat "output.thetas.gz" \
#            -nChr $(echo "2 *" $(wc -l $bamlist | cut -f 1 -d " ") | bc) \
#            || { echo "estimating Tajima's D failed"; exit 1; }
#        ;;
#
#
#        "allele-frequency")
#        # This step estimates the per-site allele frequency based on
#        # genotype likelihoods.
#        # -doMajorMinor 4 = use reference allele as major
#        # -doMaf 2 = Calculate frequency assuming known major allele
#        # -indF = inbreeding coefficients
#        # -GL 2 = calculate genotype likel. using GATK method
#        # -doGlf 3 = print out GLs in binary form (to save space)
#        # -HWE_pval 1 = calculate HWE p-values
#        # Actually works very similarly to the genotype-likelihood step.
#        # Also, I removed the requirement that sites are likely to
#        # be SNPs; the results will take up more space, but this will
#        # facilitate comparisons among different populations with
#        # different segregating sites.
#        if [ $(wc -l $bamlist | cut -f 1 -d " ") == 0 ]; then
#            echo "Skipping because no files - ${pop} ${year}"
#            exit
#        fi
#        angsd -bam $bamlist \
#            -rf $regionsfile \
#            -minMapQ $MINMAPQ \
#            -minQ $MINBASEQ \
#            -ref $REFERENCE \
#            -baq $BAQ \
#            -P $NPROCS \
#            -doMajorMinor 4 \
#            -doMaf 2 \
#            -indF "../inbreeding/output" \
#            -GL 2 \
#            -doGlf 3 \
#            -HWE_pval 1 \
#            -out "output" \
#            || { echo "estimating allele frequencies failed for ${pop} ${year}"; exit 1; }
#        ;;
#
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
        if [[ $(wc -l $bamlist | cut -f 1 -d " ") -eq 0 ]]; then
            echo "Skipping because no files - ${pop} ${year}"
            exit
        fi
        
	#for bam in $(cat $bamlist); do
        #        samtools index $bam
        #done

	# Experiment to see if I can get certain populations to not freeze
	## rf only analyze regions near caputure loci
	## minMapQ mapping quality from bam files

	## bam files prepared by Soraia
	## rf are regions near capture loci, no X locus
	## ancestral sequence is the reference
	## -doMajorMinor 4 = Pre specified Major using a reference; "the major allele according to the reference states if you have defined those -ref. The minor allele will be inferred based on the genotype likelihood"
	## -doMaf 2 = Frequency (fixed major unknown minor) we sum over the 3 possible minor alleles weighted by their probabilities. The allele frequency estimator from genotype likelihoods are from Kim et al 2011, Estimation of allele frequency and association mapping using next-generation sequencing data. 
	## snp_pval 1  = Remove sites with a pvalue larger
	## GL 2 = Use the (old) GATK genotype likelihood model
	## doSad 1 = perform multisample GL estimation, Calculate the Site allele frequency likelihood based on individual genotype likelihoods assuming HWE
	## fold deprecated 	
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

### This doesn't seem requisite for estimate raw allele frequency changes
#        realSFS "output.saf.idx" \
#            -maxIter 100 \
#            -P $NPROCS \
#            > "output.sfs" \
#            || { echo "estimating SFS failed for ${pop} ${year}"; exit 1; }
# 
        ;;

#    "beagle-prep")
#        # Same as above, but only output SNPs, and actually produce the
#        # genotype likelihood file. Then beagle is run on the results.
#
#        subout="${popyearout}/beagle"
#        mkdir -p $subout
#        cd $subout
#        if [ $(wc -l $bamlist | cut -f 1 -d " ") == 0 ]; then
#            echo "Skipping because no files - ${pop} ${year}"
#            exit
#        fi
#        angsd -bam $bamlist \
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
#
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

        "bamlist")
        for pop in "Freycinet" "Narawntapu" "West_Pencil_Pine" "Forestier" "Mt_William" "Fentonbury" "Woolnorth"; do
            echo "bamlist" $pop "_" > "${datafile}.${i}"
            i=$(($i+1))
        done
        ;;

#        "genotype-likelihood")
#        #for pop in "Freycinet" "Narawntapu" "West_Pencil_Pine" "Forestier" "Mt_William" "Fentonbury" "Woolnorth"; do
#            #for year in `cat "${RESULTDIR}/${pop}/years.txt" "${RESULTDIR}/${pop}/generations.txt"` "before" "after"; do
#        for pop in "Narawntapu"; do
#            for year in "2006"; do
#                if [ $(wc -l "${RESULTDIR}/${pop}/${year}/bamlist.txt" | cut -f 1 -d " ") != 0 ]; then
#                    echo "genotype-likelihood" $pop $year > "${datafile}.${i}"
#                    i=$(($i+1))
#                fi
#            done
#        done
#        ;;

        "simple-run")
        for pop in "Freycinet" "Narawntapu" "West_Pencil_Pine" "Forestier" "Mt_William" "Fentonbury"; do
        #for pop in "Forestier"; do    
	for year in `cat "${RESULTDIR}/${pop}/generations.txt"` "before" "after"; do
       	if [[ $(wc -l "${RESULTDIR}/${pop}/${year}/bamlist.txt" | cut -f 1 -d " ") != 0 ]]; then
                    echo "simple-run" $pop $year > "${datafile}.${i}"
                    i=$(($i+1))
                fi
            done
        done
	;;

        ## Re-runs with a longer time limit
        #echo "simple-run" "Forestier" "after" > "${datafile}.${i}"
        #i=1
        #echo "simple-run" "Woolnorth" "before" > "${datafile}.${i}"
        #i=2
        # Experiment to see if I can get certain populations to not freeze
        #i=0
        #echo "simple-run" "Woolnorth" "before" > "${datafile}.${i}"
        #i=1
        #i=3
        #echo "simple-run" "Woolnorth" "before" > "${datafile}.${i}"
        #i=4
        #;;

#    "beagle-prep")
#        for pop in "Freycinet" "Narawntapu" "West_Pencil_Pine" "Forestier" "Mt_William" "Fentonbury" "Woolnorth"; do
#            #for year in "before" "after"; do
#            for year in "all"; do
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
#            #for year in "before" "after"; do
#            for year in "all"; do
#                for chr in "chr1" "chr2" "chr3" "chr4" "chr5" "chr6"; do
#                    if [ $(wc -l "${RESULTDIR}/${pop}/${year}/bamlist.txt" | cut -f 1 -d " ") != 0 ]; then
#                        echo "beagle" $pop $year $chr > "${datafile}.${i}"
#                        i=$(($i+1))
#                    fi
#                done
#            done
#        done
#       ;;


#        "inbreeding")
#        for pop in "Freycinet" "Narawntapu" "West_Pencil_Pine" "Forestier" "Mt_William" "Fentonbury" "Woolnorth"; do
#            for year in `cat "${RESULTDIR}/${pop}/years.txt" "${RESULTDIR}/${pop}/generations.txt"` "before" "after"; do
#                if [ $(wc -l "${RESULTDIR}/${pop}/${year}/bamlist.txt" | cut -f 1 -d " ") != 0 ]; then
#                    echo "inbreeding" $pop $year > "${datafile}.${i}"
#                    i=$(($i+1))
#                fi
#            done
#        done
#        ;;

#        "folded-sfs")
#        for pop in "Freycinet" "Narawntapu" "West_Pencil_Pine" "Forestier" "Mt_William" "Fentonbury" "Woolnorth"; do
#            for year in `cat "${RESULTDIR}/${pop}/years.txt" "${RESULTDIR}/${pop}/generations.txt"` "before" "after"; do
#                if [ $(wc -l "${RESULTDIR}/${pop}/${year}/bamlist.txt" | cut -f 1 -d " ") != 0 ]; then
#                    echo "folded-sfs" $pop $year > "${datafile}.${i}"
#                    i=$(($i+1))
#                fi
#            done
#        done
#        ;;
#
#        "thetas")
#        for pop in "Freycinet" "Narawntapu" "West_Pencil_Pine" "Forestier" "Mt_William" "Fentonbury" "Woolnorth"; do
#            for year in `cat "${RESULTDIR}/${pop}/years.txt" "${RESULTDIR}/${pop}/generations.txt"` "before" "after"; do
#                if [ $(wc -l "${RESULTDIR}/${pop}/${year}/bamlist.txt" | cut -f 1 -d " ") != 0 ]; then
#                    echo "thetas" $pop $year > "${datafile}.${i}"
#                    i=$(($i+1))
#                fi
#            done
#        done
#        ;;
#
#        "allele-frequency")
#        for pop in "Freycinet" "Narawntapu" "West_Pencil_Pine" "Forestier" "Mt_William" "Fentonbury" "Woolnorth"; do
#            for year in `cat "${RESULTDIR}/${pop}/years.txt" "${RESULTDIR}/${pop}/generations.txt"` "before" "after"; do
#                if [ $(wc -l "${RESULTDIR}/${pop}/${year}/bamlist.txt" | cut -f 1 -d " ") != 0 ]; then
#                    echo "allele-frequency" $pop $year > "${datafile}.${i}"
#                    i=$(($i+1))
#                fi
#            done
#        done
#        ;;
#
#        "error-multisample")
#        for pop in "Freycinet" "Narawntapu" "West_Pencil_Pine" "Forestier" "Mt_William" "Fentonbury" "Woolnorth"; do
#            for year in `cat "${RESULTDIR}/${pop}/years.txt" "${RESULTDIR}/${pop}/generations.txt"` "before" "after"; do
#                if [ $(wc -l "${RESULTDIR}/${pop}/${year}/bamlist.txt" | cut -f 1 -d " ") != 0 ]; then
#                    echo "error-multisample" $pop $year > "${datafile}.${i}"
#                    i=$(($i+1))
#                fi
#            done
#        done
#        ;;

    esac

	sbatch --job-name="afchange-simplerun" --mem=120G --mail-user=amandastahlke@gmail.com --mail-type=ALL --array=0-$(($i-1)) $sfile
    #qsub -l "nodes=1:ppn=${NPROCS}" -q long -t 0-$(($i-1)) $sfile
    echo "Submitted ${i} jobs"
fi
