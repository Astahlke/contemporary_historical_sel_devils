
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

######## CODEML - Branch site models #######

# a list structure would probably work here
colnames <- c("ensembl_gene_id", "likelihood", "treelength", "p0", 
              "b0", "f0", "p1", "b1", "f1", "p2", "b2", "f2", "p3", "b3", "f3")

paml_posi <- read.table(file = "paml_out/paml_posiBS_wgcan_koala_111617.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
paml_posi[1] <-  strtrim(paml_posi[,1], 18)
colnames(paml_posi) <- colnames

paml_null <- read.table(file = "paml_out/paml_nullBS_wgcan_koala_111617.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
paml_null[1] <-  strtrim(paml_null[,1], 18)
colnames(paml_null) <- colnames


#### Candidate gene lists ####
candidategene_lists <- c("compstat_candidates") 
                         # "af_compstat_candidates", 
                         # "timeseries_compstat_candidates",
                         # "af_parallel_candidates",
                         # "top_mm_candidates",
                         # "top_spat",
                          # "combined_candidates") 
                         # "oldcandidategenes")
compositestat_candidates <- read.csv("all_composite_stat_results/2019-2-22/results/annotation_top1.0/composite.snps.everything.top.genes.100000bp.txt", 
                                     fill = TRUE, sep="\t", stringsAsFactors = FALSE, header = FALSE)
# compositestat_candidates <- read.csv("compstat_results_20190114/annotation_top1.0/composite.snps.everything.top.genes.100000bp.txt", fill = TRUE, 
#                                      sep="\t", 
#                                      stringsAsFactors = FALSE, header = FALSE)
# combined_candidates <- read.csv("compstat_results_20190114/annotation_top1.0/combined.genes.100000bp.tsv", fill = TRUE, 
#                                      sep="\t", 
#                                      stringsAsFactors = FALSE)

nrow(compositestat_candidates) # 252
length(unique(compositestat_candidates$V8)) #186

compstat_candidates <- strtrim(compositestat_candidates[,'V8'], 18)
compstat_candidates <- as.data.frame(unique(compstat_candidates))
colnames(compstat_candidates) <- "ensembl_gene_id"
nrow(compstat_candidates) #186

candidates_posi <- semi_join(paml_posi, compstat_candidates)
wg_posi <- anti_join(paml_posi, compstat_candidates)
nrow(candidates_posi)

nrow(candidates_posi)/nrow(compstat_candidates) # 36% were tested


candidates_null <- semi_join(paml_null, compstat_candidates)
wg_null <- anti_join(paml_null, compstat_candidates)
nrow(wg_null)

nrow(candidates_null) - nrow(candidates_posi) #### SHOULD BE ZERO

interest_cols <- c(c("ensembl_gene_id", "likelihood", "p0", "f0", "p1", "f1", "p2", "f2", "p3", "f3"))
can_posi <- candidates_posi[,interest_cols]
can_null <- candidates_null[,interest_cols]
wg_posi <- wg_posi[,interest_cols]
wg_null <- wg_null[,interest_cols]


## f3>1 when under positive selection or constrained by model
## f3 is omega value of interest
can_posi.f3w1 <- as.data.frame(can_posi[(which(can_posi$f3>1 )),])
wg_posi.f3w1 <- as.data.frame(wg_posi[(which(wg_posi$f3>1)), ])
head(can_posi.f3w1)

### What proportion of genes have f3 sites estimated to be under positive selection? 
nrow(wg_posi.f3w1)/nrow(wg_posi) # 42%
nrow(can_posi.f3w1)/nrow(can_posi) # 38% 

wg_posi_null_f3.df <- merge(wg_posi, wg_null, by = "ensembl_gene_id") 
can_posi_null_f3.df <- merge(can_posi, can_null, by = "ensembl_gene_id")
head(wg_posi_null_f3.df)
nrow(wg_posi_null_f3.df)

## select columns for LRT
## likelihood ratio test on genes potentially under selection
LRT_cols <- c("ensembl_gene_id", "likelihood.x", "f3.x", "likelihood.y", "f3.y")
wg_posi_null_f3.df <- wg_posi_null_f3.df[,LRT_cols]
can_posi_null_f3.df <- can_posi_null_f3.df[,LRT_cols]
head(can_posi_null_f3.df)

new_LRT_cols <- c('ensembl_gene_id', 'likelihood_1', 'f3_1', 'likelihood_0', 'f3_0')

## wholegenome
colnames(wg_posi_null_f3.df) <- new_LRT_cols 
wg_lrt <- matrix(nrow = nrow(wg_posi_null_f3.df), ncol=2)
wg_lrt[,1] <-  2*(wg_posi_null_f3.df[,'likelihood_1']-wg_posi_null_f3.df[,'likelihood_0'])
wg_lrt[,2] <- pchisq(as.numeric(wg_lrt[,1]),1, lower.tail = FALSE)
colnames(wg_lrt) <- c("LRT", "pval")
head(wg_lrt)

## candidates
colnames(can_posi_null_f3.df) <- new_LRT_cols
can_lrt <- matrix(nrow = nrow(can_posi_null_f3.df), ncol=2)
can_lrt[,1] <-  2*(can_posi_null_f3.df[,'likelihood_1']-can_posi_null_f3.df[,'likelihood_0'])
can_lrt[,2] <- pchisq(as.numeric(can_lrt[,1]),1, lower.tail = FALSE)
colnames(can_lrt) <- c("LRT", "pval")
head(can_lrt)

## combine LRT with f3 omega estimates and geneID
can_LRT_f3.df <- cbind(can_posi_null_f3.df, can_lrt)
wg_LRT_f3.df <- cbind(wg_posi_null_f3.df, wg_lrt)
head(can_LRT_f3.df)

# Genes which have dN/dS rates significantly greater than 1 
## contemporary, compositie 
sig_can_LRT_f3.df <- subset(can_LRT_f3.df, pval <=.05)
nonsig_can_LRT_f3.df <- subset(can_LRT_f3.df,  pval>.05)
head(sig_can_LRT_f3.df)
nrow(sig_can_LRT_f3.df)/nrow(nonsig_can_LRT_f3.df)*100

sig_can_LRT_f3.df <- merge(can_posi[,interest_cols], sig_can_LRT_f3.df[, c('ensembl_gene_id', 'pval')], by = 'ensembl_gene_id')
nonsig_can_LRT_f3.df <- merge(can_posi[,interest_cols], nonsig_can_LRT_f3.df[, c('ensembl_gene_id', 'pval')], by = 'ensembl_gene_id')
sig_can_LRT_f3.df <- mutate(sig_can_LRT_f3.df, p2a_2b= p2+p3)
head(sig_can_LRT_f3.df)

## historical, GWB 
sig_wg_LRT_f3.df <- subset(wg_LRT_f3.df, pval <=.05)
nonsig_wg_LRT_f3.df <- subset(wg_LRT_f3.df,  pval>.05)
nrow(sig_wg_LRT_f3.df)/nrow(nonsig_wg_LRT_f3.df)*100

sig_wg_LRT_f3.df <- merge(wg_posi[,interest_cols], sig_w
                          _LRT_f3.df[, c('ensembl_gene_id', 'pval')], by = 'ensembl_gene_id')#21
nonsig_wg_LRT_f3.df <- merge(wg_posi[,interest_cols], nonsig_wg_LRT_f3.df[, c('ensembl_gene_id', 'pval')], by = 'ensembl_gene_id')#21
sig_wg_LRT_f3.df <- mutate(sig_wg_LRT_f3.df, p2a_2b= p2+p3)

sig_can_LRT_f3.df$ID <- "sig_can"
sig_wg_LRT_f3.df$ID <- "sig_wg"
nonsig_wg_LRT_f3.df$ID <- "nonsig_wg"
nonsig_can_LRT_f3.df$ID <- "nonsig_can"

colnames(nonsig_can_LRT_f3.df)
colnames(nonsig_wg_LRT_f3.df)

nonsig_can_wg_paml <- rbind(nonsig_wg_LRT_f3.df, nonsig_can_LRT_f3.df)
write.csv(nonsig_can_LRT_f3.df, "nonsig_can_LRT_f3_20190420.csv")

prop_wg_under_positive_sel <- nrow(sig_wg_LRT_f3.df)/nrow(wg_LRT_f3.df) 
prop_wg_under_positive_sel # 31.76%
prop_can_under_positive_sel <-nrow(sig_can_LRT_f3.df)/nrow(can_LRT_f3.df) 
prop_can_under_positive_sel # 27.94%

# sig_cans <- as.data.frame(sig_can_LRT_f3.df$ensembl_gene_id)
# colnames(sig_cans) <- "geneID"

# can_wg_paml <- rbind(wg_posi, can_posi)
sig_can_wg_paml <- rbind(sig_wg_LRT_f3.df, sig_can_LRT_f3.df)
nonsig_can_wg_paml <- mutate(nonsig_can_wg_paml, p2a_2b= p2+p3)

sig_nonsig_can_wg_paml <- rbind(sig_can_wg_paml, nonsig_can_wg_paml)

# write.csv(sig_nonsig_can_wg_paml, file = "all_paml_results_20190116.csv", quote=FALSE, row.names = FALSE)

allpamlresults <- read.csv("all_paml_results_20190116.csv")

# Of those genes with sites in the statistically significant f3 class 
# ggplot(filter(sig_can_wg_paml, ID=="can" | ID=="wg")) +
#   geom_density(kernel = "gaussian", adjust = .2, (aes(f3, fill = ID))) +
#   scale_fill_manual(name = c("Gene set" ), labels = c("Candidates","Whole Genome"), values = alpha(c("red", "blue"), .3)) +
#   #geom_vline(xintercept = 1, color = "black")  +
#   #labs(title = "") +
#   xlab("Class 2a Omega Estimates \n (2b Omega estimate significantly greater than 1)") +
#   #xlim(0, 1000) +
#   theme_minimal()
library(ggplot2)
plot <- ggplot(allpamlresults) +
  geom_histogram(position = "dodge", bins = 100, aes(f3, ..density.., fill = ID)) +
  # scale_fill_discrete(name = c("Gene set" ), labels = c("Candidates","Whole Genome")) +
  # geom_vline(xintercept = 1, color = "black")  +
  labs(title = "Omega estimate for site class 2") +
  # ggtitle(paste(candidategene_lists[i], "vs whole genome")) +
  # xlab("Proportion of sites in Class2, not necessarily significant Class 2 dN/dS estimates > 1") +
  theme_minimal()
plot


plot <- ggplot(sig_can_wg_paml) +
  geom_histogram(position = "dodge", bins =100, aes(p2a_2b, ..density.., fill = ID)) +
  scale_fill_discrete(name = c("Gene set" ), labels = c("Candidates","Whole Genome")) +
  # geom_vline(xintercept = 1, color = "black")  +
  # labs(title = "Omega estimate for site class 2") +
  # ggtitle(paste(candidategene_lists[i], "vs whole genome")) +
  xlab("Proportion of sites under positive selection per gene \n
       ") +
  theme_minimal()
plot


ks <- ks.test(filter(sig_can_wg_paml, ID=="sig_can")$f3,  
              filter(sig_can_wg_paml, ID=="sig_wg")$f3)
assign(paste0(candidategene_lists[i], "_ks"), ks)
# ks.test(merge(sig_can_wg_paml, regression_casecontrol)$f3,  filter(sig_can_wg_paml, ID=="wg")$f3)

install.packages("kSamples")
library(kSamples)
set.seed(2627)
kSamples::ad.test(filter(sig_can_wg_paml, ID=="sig_can")$f3,  
                  filter(sig_can_wg_paml, ID=="sig_wg")$f3, method = "exact")
kSamples::ad.test(filter(sig_can_wg_paml, ID=="sig_can")$p2a_2b,  
                  filter(sig_can_wg_paml, ID=="sig_wg")$p2a_2b, method = "simulated")

can_conting_tbl <- rbind(nrow(sig_can_LRT_f3.df), nrow(nonsig_can_LRT_f3.df))
wg_conting_tbl <- rbind(nrow(sig_wg_LRT_f3.df), nrow(nonsig_wg_LRT_f3.df))
conting_tbl <- cbind(can_conting_tbl, wg_conting_tbl)
colnames(conting_tbl) <- c("can", "wg")
rownames(conting_tbl) <- c("sig", "nonsig")
conting_tbl
assign(paste0(candidategene_lists[i], "_conting_tbl"), conting_tbl)
fish <- fisher.test(conting_tbl)
fish
assign(paste0(candidategene_lists[i], "_fish"), fish)

totals <- rbind(ngenes_candidatelist, prop_can_genes_analyzed, 
                prop_can_under_positive_sel, prop_wg_under_positive_sel, 
                fish$p.value, ks$p.value)
colnames(totals) <- paste0(candidategene_lists[i], "_totals")
rownames(totals)[5:6] <- c("fish$p.value", "ks$p.value")
assign(paste0(candidategene_lists[i], "_totals"), totals)
 # 51 have already been identified as candidates

alltotals <- c()
for (i in 1:length(candidategene_lists)){
  alltotals <- cbind(get(paste0(candidategene_lists[i], "_totals")), alltotals)
}
alltotals

write.csv(alltotals, "alltotals_candidategenelists.csv", quote=FALSE)
write.csv(compstat_candidates_sig_nonsig_can_wg_paml, "compstat_candidates_sig_nonsig_can_wg_paml.csv", quote=FALSE, row.names = FALSE)
write.csv(sig_can_wg_paml, "compstat_candidates_sig_can_wg_paml.csv", quote=FALSE, row.names = FALSE)
write.csv(sig_compstat_candidates, "sig_compstat_candidates.csv", quote=FALSE, row.names = FALSE)
nrow(sig_can_wg_paml) -26


oldcandidategenes_plot
combined_candidates_plot
timeseries_compstat_candidates_plot
compstat_candidates_plot
af_compstat_candidates_plot

shared_lg_pre <- merge(landscapegenomics_pre, sig_can_wg_paml)
shared_lg_post <- merge(landscapegenomics_post, sig_can_wg_paml)
merge(shared_lg_post, shared_lg_pre)

regression_casecontrol <- read.csv("~/Downloads/Host_wcFst_segments.csv")
regression_casecontrol$ensembl_gene_id <- str_split(regression_casecontrol$Gene_Annotation, ":", simplify = TRUE)[,1]
merge(sig_can_wg_paml, regression_casecontrol, by = "ensembl_gene_id")
sig_can_wg_paml[which(merge(sig_can_wg_paml, regression_casecontrol)$ensembl_gene_id==sig_can_wg_paml),]$ID <- "can"

nrow(merge(sig_can_wg_paml, regression_casecontrol))/ nrow(sig_can_LRT_f3.df)
## write.csv(merge(sig_can_wg_paml, regression_casecontrol, by = "geneID"), "regression_casecontrol_wg_can_100917.csv")

matches <- merge(sig_can_wg_paml, landscapegenomics_post)
for (i in 1:nrow(matches)){
  geneID <- matches[i,1]
  sig_can_wg_paml[which(sig_can_wg_paml$ensembl_gene_id==geneID),]$ID <- "can"
}

sig_candidates <- filter(sig_can_wg_paml, ID=="can")
nrow(sig_candidates)
## write.csv(sig_candidates, "sig_candidates_101817.csv")
# 
# twotailed_stats <- c()
# twotailed_pvals <- c()
# greater_stats <- c()
# greater_pvals <- c()

## resample parametric bootsrap
# for(i in 1:1000){
#   sample <- sample_n(sig_wg_LRT_f3.df, nrow(sig_can_LRT_f3.df), replace = FALSE)$p2a_2b
#   
#   wilcox_twotailed <- wilcox.test(sig_can_LRT_f3.df$p2a_2b, sample, alternative = "two.sided")
#   twotailed_stats <- rbind(twotailed_stats, wilcox_twotailed$statistic)
#   twotailed_pvals <- rbind(twotailed_pvals, wilcox_twotailed$p.value)
#   
#   wilcox_greater <- wilcox.test(sig_can_LRT_f3.df$p2a_2b, sample, alternative = "greater")
#   greater_stats <- rbind(greater_stats, wilcox_greater$statistic)
#   greater_pvals <- rbind(greater_pvals, wilcox_greater$p.value)
# }
# 
# count(twotailed_pvals<.1)$freq/1000 # reject the null of same 5% of the time, greater 13% of the time
# count(greater_pvals<.1)$freq/1000 # reject the null of same 5% of the time, greater 13% of the time
# 
# hist(twotailed_pvals)
# hist(greater_pvals)

# oppwall <- read.table("~/Code/devils/add_koala/oppwall.txt",  sep = "\t") 
# oppwallgenes <- str_split(oppwall[,1], "[.]", simplify = TRUE)[,1]
# oppwallgenes <- data.frame(oppwallgenes)
# colnames(oppwallgenes) <- "geneID"
# compare_trees <- semi_join(compare_trees, oppwallgenes)
# 
# compare_trees <- merge(old_candidates_posi, newtree_candidates_posi, by = "geneID")
# compare_trees <- mutate(compare_withkoala, "f3diff" = f3.y-f3.x, "likediff" = likelihood.y - likelihood.x, "treediff" = treelength.y-treelength.x, "p3diff"=p3.y-p3.x)
# 
# hist(x=compare_trees$p3diff, breaks = 100)
## write.csv(compare_trees, "compare_trees.csv")
