###
### Goal is to compare lists of candidate genes to better understand functionality etc
### 

### Set up to run locally ~/amandastahlke/Code/devils/ - ideally this would be synced with git etc
getwd()
options('scipen'=50)

source("historical/script/top_script.R")

#### contemporary selection ####
contemporarydir <- "all_composite_stat_results/2020-10-22/results/"

## main compstat results
compstat_results <- read.csv(paste0(contemporarydir, "annotation_top1.0/composite.snps.everything.top.genes.100000bp.txt"), 
                             row.names = NULL, sep = "\t", header = FALSE)
head(compstat_results)
compstat_results$ensembl_gene_id <- strtrim((compstat_results$V8), 18) 
length(unique(compstat_results$ensembl_gene_id)) 

contemporary_universe <- read.delim(paste0(contemporarydir, "annotation_top100.0/composite.snps.everything.top.genes.100000bp.txt"),
                                    header = FALSE, col.names = c("beg_scaffold", 
                                                                  "beg_chr_start", "beg_chr_stop", "compstat", 
                                                                  "end_scaffold", "end_chr_start", "end_chr_stop", 
                                                                  "annotation"))
contemporary_universe <- str_split(contemporary_universe$annotation, ":", simplify = TRUE)[,1]



### allele frequency for each population ####

infiles <- list.files("all_composite_stat_results/2020-10-22/results/annotation_top1.0/by_pop_method", 
                      pattern =  glob2rx("*.top.genes.*.txt"), full.names = TRUE)

infiles
indata = vector('list', length=length(infiles))

for(i in seq_along(infiles)) {
  d = read.csv(infiles[i], sep='\t', header=F, as.is=TRUE)
  indata[[i]] = str_split(d[,'V8'], ":", simplify = T)[,1]
}

nms2 <- str_split(str_split(infiles, "/", simplify = T)[,6], "[[.]]", simplify = T)[,1:2]
names(indata) <- paste(nms2[,1], nms2[,2], sep = ".")
              
                                     
nms <- combn( names(indata) , 2 , FUN = paste0 , collapse = "-" , simplify = FALSE )
ll <- combn(indata, 2 , simplify = FALSE )
out <- lapply( ll , function(x) length( intersect( x[[1]] , x[[2]] ) ) )
out <- setNames( out , nms )
out <- unlist(out)

out_names <- as.data.frame(str_split(names(out), "-", simplify = T))
out_names$V3 <- unlist(out)


myNames <- sort(unique(as.character(unlist(out_names[1:2]))))
myMat <- matrix("-", length(myNames), length(myNames), dimnames = list(myNames, myNames))
myMat[as.matrix(out_names[c(1,2)])] <- out_names$V3
diag(myMat) <- unlist(lapply(indata, function(x) length(x)))


## Is any one list of candidates over represented in composite list?
DCMS <- unlist(lapply(indata, function(x) length(intersect (x, compstat_results$ensembl_gene_id))))
myMat <- rbind(DCMS, myMat)
myMat <- cbind(DCMS = c(length(compstat_results$ensembl_gene_id), rep_len("-", length.out = length(infiles))), myMat)
myMat
write.csv(myMat, "candidate_list_intersections.csv", quote = F)


#### historical selection ####
allpamlresults <- read.csv("historical/paml_LRT_21-1-13.csv")
sig_paml <- read.csv("historical/sig_paml_results_21-1-13.csv")

nrow(sig_paml)
nrow(sig_paml)/nrow(allpamlresults)

sig_can_paml <- semi_join(sig_paml, compstat_results, by ="ensembl_gene_id") %>% 
  mutate(ID = "recurrent")
nrow(sig_can_paml)

#### Previous eresults
#### Espstein et al 2016 Rapid evolutionary response to a transmissible cancer in Tasmanian devils ####
epstein_candidates <- read.csv("Epstein2016candidates", header = FALSE, 
                               col.names = "ensembl_gene_id", stringsAsFactors = FALSE)

intersect(epstein_candidates$ensembl_gene_id, compstat_results$ensembl_gene_id)
setdiff(epstein_candidates$ensembl_gene_id, compstat_results$ensembl_gene_id)


#### Margres et al 2018 Few loci of large effect S5 ####
Margres_cans_full <- read.csv("mec14853-sup-0004-tables5.tsv", sep = "\t",
                              header = TRUE, stringsAsFactors = FALSE)
# View(Margres_cans_full)
Margres_cans_full <- mutate(Margres_cans_full, 
                            ensembl_gene_id = str_split(Margres_cans_full$gene, ":", simplify = TRUE)[,1])

library(GeneOverlap)

## What if we filter these by phenotype? 
Margres_cans_case_control <- dplyr::filter(Margres_cans_full, phenotype == "case_control")
go.obj <- newGeneOverlap(listA = compstat_results$ensembl_gene_id,
                         listB = Margres_cans_case_control$ensembl_gene_id,
                         genome.size = length(contemporary_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

Margres_cans_age <- dplyr::filter(Margres_cans_full, phenotype == "age")
go.obj <- newGeneOverlap(listA = compstat_results$ensembl_gene_id,
                         listB = Margres_cans_age$ensembl_gene_id,
                         genome.size = length(contemporary_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

Margres_cans_survival <- dplyr::filter(Margres_cans_full, phenotype == "survival")
go.obj <- newGeneOverlap(listA = compstat_results$ensembl_gene_id,
                         listB = Margres_cans_survival$ensembl_gene_id,
                         genome.size = length(contemporary_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

Margres_cans_female <- dplyr::filter(Margres_cans_full, sex == "female")
go.obj <- newGeneOverlap(listA = compstat_results$ensembl_gene_id,
                         listB = Margres_cans_female$ensembl_gene_id,
                         genome.size = length(contemporary_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

Margres_cans_both <- dplyr::filter(Margres_cans_full, sex == "both")
go.obj <- newGeneOverlap(listA = compstat_results$ensembl_gene_id,
                         listB = Margres_cans_both$ensembl_gene_id,
                         genome.size = length(contemporary_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

Margres_cans_male <- dplyr::filter(Margres_cans_full, sex == "male")
go.obj <- newGeneOverlap(listA = compstat_results$ensembl_gene_id,
                         listB = Margres_cans_male$ensembl_gene_id,
                         genome.size = length(contemporary_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

### Margres et al 2018 The Genomic Basis of Tumor Regression in Tasmanian Devils (Sarcophilus harrisii) ###
### Table S2: "We identified 1148 SNPs in the 11 SNP-based candidate genomic regions described above (supplementary table S2, Supplementary Material online)"
### ** Need to filter these down to those "highly differentiated" fst > 0.5
Margres_regression_snps <- read.csv("evy229_supp/Table_S2.csv") %>% 
  dplyr::filter(abs(WC_FST) > 0.5) %>% 
  dplyr::select(VEP) %>% 
  unique(.) %>% 
  rename(ensembl_gene_id=VEP)

Margres_regression_indel <- read.csv("evy229_supp/Table_S3.csv") %>% 
  dplyr::filter(abs(wcFST) > 0.5) %>% 
  dplyr::select(VEP) %>% 
  unique(.) %>% 
  rename(ensembl_gene_id=VEP)

## one weird case where two genes are listed
Margres_regression_indel <- full_join(data.frame(ensembl_gene_id = t(str_split(Margres_regression_indel[1,], "_or_", simplify = T))),
      Margres_regression_indel)

## Any interesting overlaps?
merge(compstat_results, Margres_regression_snps, by = "ensembl_gene_id") ## ENSSHAG00000006010 (JAKIMP3)
merge(compstat_results, Margres_regression_indel, by = "ensembl_gene_id") ## none

merge(sig_can_paml, Margres_regression_snps, by = "ensembl_gene_id") ## none (JAKIMP3)
merge(sig_can_paml, Margres_regression_indel, by = "ensembl_gene_id") ## none

lapply(indata, function(x) intersect (x, Margres_regression_snps$ensembl_gene_id))
lapply(indata, function(x) intersect (x, Margres_regression_indel$ensembl_gene_id))


### Fraik et al Disease swamps molecular signatures of genetic‐environmental associations to abiotic factors in Tasmanian devil (Sarcophilus harrisii) populations ###
### Supplementary Table 5. The 71 candidate genes detected pre-disease using landscape genomics analyses. ENSEMBL ID, NCBI gene function, NCBI gene ID, associated GO terms are included for each candidate gene. 
### Includes the twenty-four candidate genes that were shared between the pre and post-disease candidate gene lists are colored grey.
### Supplementary Table 6. The 81 candidate genes uniquely detected post-disease using landscape genomics analyses. Each candidate gene is identifiable with an ENSEMBL ID, and an NCBI gene function and gene ID, and includes all of the associated GO terms.

fraik_s5 <- read.csv("fraik_evo_14023_supp/fraik_evo_14023_supp_S5.csv") %>% 
  rename(ensembl_gene_id=ENSEMBL.ID) %>% 
  dplyr::select(ensembl_gene_id)

fraik_s6 <- read.csv("fraik_evo_14023_supp/fraik_evo_14023_supp_S6.csv") %>% 
  rename(ensembl_gene_id=ENSEMBL.ID) %>% 
  dplyr::select(ensembl_gene_id)

semi_join(compstat_results, fraik_s5) # 7 genes
semi_join(compstat_results, fraik_s6)

semi_join(sig_can_paml, fraik_s5) # ENSSHAG00000013654
semi_join(sig_can_paml, fraik_s6)


### Hubert

hubert_s1 <- read.csv("hubert.pone.0201838.s003.csv") %>% 
  rename(ensembl_gene_id=Ensembl.ID) %>% 
  dplyr::select(ensembl_gene_id) %>% 
  unique(.)

semi_join(compstat_results, hubert_s1)
semi_join(compstat_results, hubert_s1)
semi_join(sig_can_paml, hubert_s1) # none 

### Wright
wright <- read.csv("wright_natcomm_10.1038_s41598-017-00439-7.csv") %>% 
  rename(ensembl_gene_id=ENSEMBL.GENE.ID) %>% 
  dplyr::select(ensembl_gene_id)

semi_join(compstat_results, wright) # none
semi_join(sig_can_paml, wright) # none 

lapply(indata, function(x) intersect (x, wright$ensembl_gene_id))
