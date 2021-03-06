# Contemporary & Historical Selection in Tasmanian devils
Code to study contemporary and historical selection in devils. Manuscript in review. 

## Raw input data
### contemporary 
- Demultiplexed sequence data has been deposited at NCBI under BioProject PRJNA306495 (http://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA306495) and BioProject PRJNA634071 (http://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA634071). 
- bamlists by population and cohort are found in 
contemporary_historical_sel_devils/contemporary/angsd_2019-01-18/results/
- Genotype likelihoods and allele frequency inputs can be found on dryad:https://doi.org/10.5061/dryad.jq2bvq872 

### historical
Orthology tables for devil, koala, oppossum, and tammar wallaby are found in 
contemporary_historical_sel_devils/historical/00-getorthologID/
- KoalaDevilOrthologs_12864_2014_6686_MOESM2_ESM.csv
- wholegenome_orthologtable.tsv

## Scripts
### contemporary
- prep_angsd_2018.sh
- timeseries_angsd.sh
- snp2go.R

### historical 
- 00-biomart_getgeneID-allorthologs.R
- 01-getseq.R
- addkoala_dataprep.sh
- all_models.sh
- dataprep.sh
- my_branch-site.py
- paml_parse_inputgeneID.sh
- parse_branchsite-null.py
- parse_branchsite-posi.py
- parse_getorf.py
- process_PAML_likelihoodratiotest.R
- reorder_fasta.py
- writetree.py 

## compare
- test_for_ongoing_selection.R
- candidatelistcomparison.R

## Key Results 
- composite stat for contemporary selection: contemporary/angsd_2019-01-18/next/composite_stat/2019-2-22/results/annotation_top1.0/composite.snps.everything.top.genes.100000bp.txt 
- combined raw paml results: 
contemporary_historical_sel_devils/historical/04-parsed_paml_out/
	- paml_nullBS_wgcan_koala_111617.tsv 
	- paml_posiBS_wgcan_koala_111617.tsv
- Likelihood ratio test for all paml results: contemporary_historical_sel_devils/historical/paml_LRT_21-1-13.csv 
- Gene set overlap results: permutation histogram 
