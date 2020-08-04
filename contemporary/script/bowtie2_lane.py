#!/usr/bin/env
#bowtie2 rapture data

from os.path import join as jp
import os
import sys

VERBOSE=False

#Function definitions:
def log(txt, out):
    if VERBOSE:
        print(txt)
    out.write(txt+'\n')
    out.flush()
    
# Setup folders and paths variables:
input = "/mnt/lfs2/soraia/TD/4_combine/Rapture_2018/PHAV5_HT30/"
sample_list = os.listdir(input)
resultsDir = '/mnt/lfs2/soraia/TD/5_bowtie2/6_Rapture2018_samples/HT30_PHAV5/'
slurm_scripts = '/mnt/lfs2/soraia/TD/4_combine/scripts/'
ref="/mnt/lfs2/soraia/TD/reference/sarHar1"


## Read in samples and put them in a list:
samples = []
for i in sample_list:
	if len(i) > 1:
		samples.append(i[:-6])
		
samples=set(samples)
        
os.system('mkdir -p %s' % resultsDir)
os.system('mkdir -p %s' % slurm_scripts)

##### Run pipeline ###
for sample in samples:
    # Set up files:
	logFile = jp(resultsDir, sample + '.log')
	logCommands = open(jp(slurm_scripts, sample + '_commands.sh'), 'w')
	# bowtie2 - will output .sam files
	cmd = ' '.join(['bowtie2 -p 5 -x', ref, '--phred33 --sensitive --end-to-end -X 900 --rg-id "rapture2018" --rg "PL:ILLUMINA" --rg "SM:' + sample+ '"', '-U', input + sample + '.fq.gz -S', resultsDir + sample + '.sam'])
	log(cmd, logCommands)
	os.system(cmd)
	#convert .sam to .bam
	cmd= ' '.join(['samtools view -bS', resultsDir + sample + '.sam', '>', resultsDir + sample + '.bam'])
	log(cmd, logCommands)
	os.system(cmd)
	#sort .bam
	#cmd= ' '.join(['samtools sort', resultsDir + sample + '.bam', resultsDir + sample + '_sorted.bam'])
	cmd= ' '.join(['samtools sort -o', jp(resultsDir, sample) + "_sorted.bam", ' -@ 50', jp(resultsDir, sample + ".bam")])
	log(cmd, logCommands)
	os.system(cmd)
	#index .bam
	cmd= ' '.join(['samtools index', resultsDir + sample + '_sorted.bam', resultsDir + sample + '_sorted.bai'])
	log(cmd, logCommands)
	os.system(cmd)
	logCommands.close()