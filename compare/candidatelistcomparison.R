library(dplyr)
library(stringr)

library(ggplot2)
library(cowplot)

### Compare lists of my candidate genes
## parallel af
## spatpg 
## mm
## compstat
## paml historical

### historical selection
allpamlresults <- read.csv("all_paml_results_200707.csv", stringsAsFactors = F)

sig_can <- semi_join(allpamlresults, compstat_results,  by = "ensembl_gene_id") %>%
  filter(pval <.05) %>%
  arrange(pval) %>%
  mutate(ID="sig_can")

nonsig_can <- semi_join(allpamlresults, compstat_results,  by = "ensembl_gene_id") %>%
  filter(pval >=.05) %>%
  arrange(pval) %>%
  mutate(ID="nonsig_can")

sig_gwb <- anti_join(allpamlresults, compstat_results,  by = "ensembl_gene_id") %>%
  filter(pval < .05) %>%
  arrange(pval) %>%
  mutate(ID="sig_gwb")

nonsig_gwb <- anti_join(allpamlresults, compstat_results,  by = "ensembl_gene_id") %>%
  filter(pval >=.05) %>%
  arrange(pval) %>%
  mutate(ID="nonsig_gwb")

allpamlresults <- rbind(sig_can, nonsig_can, nonsig_gwb, sig_gwb)

sig_can_wg_paml <- filter(allpamlresults, pval < .05)

sig_can_wg_paml_genes <- filter(allpamlresults, pval < .05) %>%
  select(ensembl_gene_id)
nonsig_can_wg_paml_genes <- filter(allpamlresults, pval >= .05) %>%
  select(ensembl_gene_id)
# write.csv(allpamlresults, "all_paml_results_200707.csv", row.names = F, quote = F)
# write.csv(nonsig_can_wg_paml_genes, "nonsig_can_wg_paml_genes_200707.csv", row.names = F, quote = F)

compstat_results <- read.csv("~/Code/devils/all_composite_stat_results/2019-2-22/results/annotation_top1.0/composite.snps.everything.top.genes.100000bp.txt", 
                             row.names = NULL, sep = "\t", header = FALSE)
head(compstat_results)
compstat_results$ensembl_gene_id <- strtrim((compstat_results$V8), 18) 
length(unique(compstat_results$ensembl_gene_id)) # 186


composite_af_results <- read.csv("~/Code/devils/all_composite_stat_results/2019-2-22/results/annotation_top1.0/composite.snps.afchange.top.genes.100000bp.txt", 
                                row.names = NULL, sep = "\t", header = FALSE)
length(unique(composite_af_results$V8)) # 162


parallel_af_results <- read.csv("~/Code/devils/all_composite_stat_results/2019-2-22/results/parallel_top5.0/top.genes.100000bp.txt", 
                                    row.names = NULL, sep = "\t", header = FALSE)
length(unique(parallel_af_results$V8)) # 174

spatpg <- read.csv("~/Code/devils/all_composite_stat_results/2019-2-22/results/spatpg_significant/top.genes.100000bp.txt", 
                           row.names = NULL, sep = "\t", header = FALSE)
length(unique(spatpg$V8)) # 68

mm <- read.csv("~/Code/devils/all_composite_stat_results/2019-2-22/results/mm_significant_3_pops/top.genes.100000bp.txt", 
                       row.names = NULL, sep = "\t", header = FALSE)
head(mm)
length(unique(mm$V8)) # 49


shared_af_spatpg <- 
  strtrim(intersect(unique(parallel_af_results$V8) ,unique(spatpg$V8)), 18)
#2

shared_af_mm <- 
  strtrim(intersect(unique(parallel_af_results$V8) ,unique(mm$V8)), 18)
#13

shared_spatpg_mm <- 
  strtrim(intersect(unique(spatpg$V8) ,unique(mm$V8)), 18)
#19


omega <- ggplot(sig_can_wg_paml, aes(x= ID, y = f3, color = ID)) +
  geom_violin(show.legend = F) +
  geom_jitter(shape=16, size = 3, position=position_jitter(0.1), alpha = .5, show.legend = F) +
  scale_fill_grey() +
  scale_color_manual(values=c("black", "#b2182b")) +
  xlab("") +  
  scale_x_discrete(guide = guide_axis(n.dodge=3), labels= c("Contemporary Candidates", "Genome-wide Background")) +
  ylab("dN/dS Estimates in Site Class 2") +
  theme(text = element_text(size = 18)) +
  theme_minimal()
omega

p2ab <- ggplot(sig_can_wg_paml, aes(x= ID, y = p2a_2b, color = ID)) +
  geom_violin(show.legend = F) +
  geom_jitter(shape=16, size = 3, position=position_jitter(0.1), alpha = .5, show.legend = F) +
  scale_fill_grey() +
  scale_color_manual(values=c("black", "#b2182b")) +
  xlab("") +
  ylab("Proportion of Sites in Class 2") +
  scale_x_discrete(guide = guide_axis(n.dodge=3), labels= c("Contemporary Candidates", "Genome-wide Background")) +
  # labs(title = "Proportion of sites under positive selection per gene") +
  # ggtitle(paste(candidategene_lists[i], "vs whole genome")) +
  theme(text = element_text(size = 18)) +
  theme_minimal() 
p2ab

plot_grid(omega, p2ab, labels = c('A', 'B'))
# ggsave("PAMLresults_200707.pdf", width = 7, units = "in", height = 6, dpi = 300)

##### How do these compare to previously identified candidate genes? ####
## Main candidate gene lists are compstat top 1% and historical sel. 
sig_can
sig_gwb
compstat_results

## Externally identfied candidate genes
## 1. Margres et al 2018 tumor regression
## 2. Margres et al 2018 Few loci of large effect
## Epstein et al 2016
## 

#### Margres et al 2018 tumor regression ####
Margres_cans_tumor1 <- read.csv("Margres-regression-fst-2018.txt", sep = "\t",
                               header = TRUE, stringsAsFactors = FALSE,
                               col.names = c("ensembl_gene_id", "hgnc_symbol"))

Margres_cans_tumor2 <- read.csv("Margres-regression-2018.txt", sep = "\t",
                              header = TRUE, stringsAsFactors = FALSE,
                              col.names = c("ensembl_gene_id", "hgnc_symbol"))
Margres_cans_tumor <- rbind(Margres_cans_tumor1, Margres_cans_tumor2)

semi_join( Margres_cans_tumor, compstat_results) 
intersect( Margres_cans_tumor$ensembl_gene_id, 
           stringr::str_split(composite_af_results$V8, ":",  simplify = T)[,1])
intersect( Margres_cans_tumor$ensembl_gene_id, 
           stringr::str_split(parallel_af_results$V8, ":",  simplify = T)[,1])
intersect( Margres_cans_tumor$ensembl_gene_id, 
           stringr::str_split(spatpg$V8, ":",  simplify = T)[,1])
intersect( Margres_cans_tumor$ensembl_gene_id, 
           stringr::str_split(mm$V8, ":",  simplify = T)[,1])

semi_join( Margres_cans_tumor, sig_can_wg_paml_genes) 



#### Margres et al 2018 Few loci of large effect S5 ####
Margres_cans_full <- read.csv("~/Downloads/mec14853-sup-0004-tables5.tsv", sep = "\t",
                              header = TRUE, stringsAsFactors = FALSE)
# View(Margres_cans_full)
Margres_cans <- str_split(Margres_cans_full$gene, ":", simplify = TRUE)[,1]
Margres_cans_full <- cbind(Margres_cans_full, Margres_cans)
Margres_cans_full <- mutate(Margres_cans_full, 
                            ensembl_gene_id = str_split(Margres_cans_full$gene, ":", simplify = TRUE)[,1])

semi_join( Margres_cans_full, sig_can) %>%
  arrange(snp_rank)

semi_join( Margres_cans_full, sig_can) %>%
  arrange(snp_rank)

semi_join( Margres_cans_full, sig_can) %>%
  distinct(gene)

semi_join( Margres_cans_full, sig_gwb) %>%
  arrange(snp_rank) %>%
  select(sex) %>%
  table()

semi_join( Margres_cans_full, sig_gwb) %>%
  filter(snp_rank<=10)

semi_join(Margres_cans_full, compstat_results) %>%
  arrange(snp_rank)

semi_join(Margres_cans_full, compstat_results) %>%
  arrange(snp_rank)%>%
  select(sex) %>%
  table()

semi_join(Margres_cans_full, compstat_results) %>%
  arrange(snp_rank)%>%
  select(phenotype) %>%
  table()

#### Epstein et al 2016 ####
epstein_candidates <- read.csv("Epstein2016candidates", header = FALSE, 
                               col.names = "ensembl_gene_id", stringsAsFactors = FALSE)
semi_join(epstein_candidates, sig_gwb)
semi_join(epstein_candidates, sig_can)
merge(epstein_candidates, compstat_results)

Hubert_cans <- read.csv("hubert-etal-2018-plosone-candidates.csv",
                        stringsAsFactors = FALSE, col.names = "ensembl_gene_id")
semi_join(Hubert_cans, sig_gwb)
semi_join(Hubert_cans, sig_can)
semi_join( Hubert_cans, compstat_results)

#### Wright regression host variation ###
Wright_cans <- read.csv("Wright_2017_candidategenes.csv",
                        stringsAsFactors = FALSE)
semi_join(Wright_cans, sig_gwb)
semi_join(Wright_cans, sig_can)
semi_join(Wright_cans,compstat_results)

#### Fraik GEA pre- and post- disease ####
FraikPreDisease_cans <- read.csv("Fraik_PreDisease_Candidates.csv",
                        stringsAsFactors = FALSE)
semi_join(FraikPreDisease_cans, sig_gwb)
semi_join(FraikPreDisease_cans, sig_can)
semi_join(FraikPreDisease_cans,compstat_results)

FraikPostDisease_cans <- read.csv("Fraik_PostDisease_Candidates.csv",
                                 stringsAsFactors = FALSE)
semi_join(FraikPostDisease_cans, sig_gwb)
semi_join(FraikPostDisease_cans, sig_can)
semi_join(FraikPostDisease_cans,compstat_results)


#### defining the multi-verse ####
contemporary_universe <- read.delim("all_composite_stat_results/2019-2-22/results/annotation_top100.0/composite.snps.everything.top.genes.100000bp.txt",
                                    header = FALSE, col.names = c("beg_scaffold", 
                                                                  "beg_chr_start", "beg_chr_stop", "compstat", 
                                                                  "end_scaffold", "end_chr_start", "end_chr_stop", 
                                                                  "annotation"))

head(contemporary_universe)
contemporary_universe <- str_split(contemporary_universe$annotation, ":", simplify = TRUE)[,1]

compstat_paml_universe <- intersect(allpamlresults$ensembl_gene_id, contemporary_universe)
background <- compstat_paml_universe

allpamlresults <- read.csv("all_paml_results_200707.csv")
head(allpamlresults)
table(allpamlresults$ID)

nrow(allpamlresults[allpamlresults$ID == "sig_can",])/nrow(allpamlresults[allpamlresults$ID == "nonsig_can",])
nrow(allpamlresults[allpamlresults$ID == "sig_gwb",])/nrow(allpamlresults[allpamlresults$ID == "nonsig_gwb",])

#### compare the distributions of dN/dS estimates and propotions 
allsig <- filter(allpamlresults, pval < .05)

ks_f3 <- ks.test(filter(allsig, ID=="sig_gwb")$f3,  
                 filter(allsig, ID=="sig_can")$f3)
ks_f3

ks_p2a_2b <- ks.test(filter(allsig, ID=="sig_gwb")$p2a_2b,  
                     filter(allsig, ID=="sig_can")$p2a_2b)
ks_p2a_2b

# install.packages("kSamples")
library(kSamples)
set.seed(42)
kSamples::ad.test(filter(allsig, ID=="sig_gwb")$f3,  
                  filter(allsig, ID=="sig_can")$f3, method = "exact")
# Number of samples:  2
# Sample sizes:  1963, 19
# Number of ties: 357
# AD    T.AD  asympt. P-value  sim. P-value
# version 1: 0.3466 -0.8543           0.8963        0.8050
# version 2: 0.7960 -0.2663           0.4827        0.4836

kSamples::ad.test(filter(allsig, ID=="sig_gwb")$p2a_2b,  
                  filter(allsig, ID=="sig_can")$p2a_2b, method = "exact")

# Sample sizes:  1963, 19
# Number of ties: 181
# AD    T.AD  asympt. P-value  sim. P-value
# version 1: 0.2415 -0.9916           0.9731        0.9759
# version 2: 0.2410 -0.9927           0.9730        0.9760


### build contigency table for fisher's exact test 
B <- filter(allpamlresults, pval < .05) %>%
  select(ensembl_gene_id)

library(GeneOverlap)


go.obj_comp_margres <- newGeneOverlap(listB = compositestat_candidates$ensembl_gene_id,
                         listA = B$ensembl_gene_id,
                         genome.size = length(compstat_paml_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

length(intersect(contemporary_universe, allpamlresults$ensembl_gene_id))

go.obj_historical_margres <- newGeneOverlap(listA = B$ensembl_gene_id,
                         listB = Margres_cans_full$ensembl_gene_id,
                         genome.size = length(
                           union(contemporary_universe, 
                                     allpamlresults$ensembl_gene_id)))
go.obj_historical_margres <- testGeneOverlap(go.obj_historical_margres)
print(go.obj_historical_margres)

### What if we filter these by phenotype? 
Margres_cans_case_control <- filter(Margres_cans_full, phenotype == "case_control")
go.obj <- newGeneOverlap(listA = compositestat_candidates$ensembl_gene_id,
                         listB = Margres_cans_case_control$ensembl_gene_id,
                         genome.size = length(contemporary_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

Margres_cans_age <- filter(Margres_cans_full, phenotype == "age")
go.obj <- newGeneOverlap(listA = compositestat_candidates$ensembl_gene_id,
                         listB = Margres_cans_age$ensembl_gene_id,
                         genome.size = length(contemporary_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

Margres_cans_survival <- filter(Margres_cans_full, phenotype == "survival")
go.obj <- newGeneOverlap(listA = compositestat_candidates$ensembl_gene_id,
                         listB = Margres_cans_survival$ensembl_gene_id,
                         genome.size = length(contemporary_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

Margres_cans_female <- filter(Margres_cans_full, sex == "female")
go.obj <- newGeneOverlap(listA = compositestat_candidates$ensembl_gene_id,
                         listB = Margres_cans_female$ensembl_gene_id,
                         genome.size = length(contemporary_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

Margres_cans_both <- filter(Margres_cans_full, sex == "both")
go.obj <- newGeneOverlap(listA = compositestat_candidates$ensembl_gene_id,
                         listB = Margres_cans_both$ensembl_gene_id,
                         genome.size = length(contemporary_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

Margres_cans_male <- filter(Margres_cans_full, sex == "male")
go.obj <- newGeneOverlap(listA = compositestat_candidates$ensembl_gene_id,
                         listB = Margres_cans_male$ensembl_gene_id,
                         genome.size = length(contemporary_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)


### Define an appropriate universe for this one
go.obj <- newGeneOverlap(listB = as.character(B$ensembl_gene_id),
                         listA = unique(Margres_cans_full$ensembl_gene_id),
                         genome.size = length(contemporary_universe))
go.obj <- testGeneOverlap(go.obj)
print(go.obj)
