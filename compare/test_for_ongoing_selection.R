### Test for ongoing selection: contemporary & historical

library(dplyr)
library(stringr)

library(ggplot2)
library(cowplot)

#### contemporary selection ####
contemporarydir <- "~/Code/devils/all_composite_stat_results/2019-2-22/results/"

## main compstat results
compstat_results <- read.csv(paste0(contemporarydir, "annotation_top1.0/composite.snps.everything.top.genes.100000bp.txt"), 
                             row.names = NULL, sep = "\t", header = FALSE)
head(compstat_results)
compstat_results$ensembl_gene_id <- strtrim((compstat_results$V8), 18) 
length(unique(compstat_results$ensembl_gene_id)) # 186

## delta allele freq 
composite_af_results <- read.csv(paste0(contemporarydir, "annotation_top1.0/composite.snps.afchange.top.genes.100000bp.txt"), 
                                 row.names = NULL, sep = "\t", header = FALSE)
length(unique(composite_af_results$V8)) # 162

parallel_af_results <- read.csv(paste0(contemporarydir, "parallel_top5.0/top.genes.100000bp.txt"), 
                                row.names = NULL, sep = "\t", header = FALSE)
length(unique(parallel_af_results$V8)) # 174

## spatpg
spatpg <- read.csv(paste0(contemporarydir, "spatpg_significant/top.genes.100000bp.txt"), 
                   row.names = NULL, sep = "\t", header = FALSE)
length(unique(spatpg$V8)) # 68

## spatpg
mm <-  read.csv(paste0(contemporarydir, "mm_significant_3_pops/top.genes.100000bp.txt"), 
                row.names = NULL, sep = "\t", header = FALSE)
head(mm)
length(unique(mm$V8)) # 49

### how many genes are shared between these tests? 
shared_af_spatpg <- 
  strtrim(intersect(unique(parallel_af_results$V8) ,unique(spatpg$V8)), 18)
#2

shared_af_mm <- 
  strtrim(intersect(unique(parallel_af_results$V8) ,unique(mm$V8)), 18)
#13

shared_spatpg_mm <- 
  strtrim(intersect(unique(spatpg$V8) ,unique(mm$V8)), 18)
#19

#### historical selection ####
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

sig_can_wg_paml <- rbind(sig_can, sig_gwb)
sig_can_wg_paml_genes <- select(sig_can_wg_paml, ensembl_gene_id)

nonsig_can_wg_paml <- rbind(nonsig_can, nonsig_gwb)
nonsig_can_wg_paml_genes <- select(nonsig_can_wg_paml, ensembl_gene_id)

#### Compare contempororay and historical selection!! ####
nrow(sig_can_wg_paml_genes)/nrow(allpamlresults)
nrow(sig_can)/(nrow(nonsig_can)+nrow(sig_can))
nrow(sig_gwb)/(nrow(nonsig_gwb)+nrow(sig_gwb))

## magnitude of dn/ds
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

## proportion of sites
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

### statistical tests comparing distributions
# install.packages("kSamples")
library(kSamples)
set.seed(42)
allsig <- filter(allpamlresults, pval < .05)

andersondarling_f3 <- kSamples::ad.test(filter(allsig, ID=="sig_gwb")$f3,  
                  filter(allsig, ID=="sig_can")$f3, method = "exact")
andersondarling_f3

andersondarling_p2a2b <- kSamples::ad.test(filter(allsig, ID=="sig_gwb")$p2a_2b,  
                  filter(allsig, ID=="sig_can")$p2a_2b, method = "exact")
andersondarling_p2a2b

#### compare the distributions of dN/dS estimates and propotions 
ks_f3 <- ks.test(filter(allsig, ID=="sig_gwb")$f3,  
                 filter(allsig, ID=="sig_can")$f3)
ks_f3

ks_p2a_2b <- ks.test(filter(allsig, ID=="sig_gwb")$p2a_2b,  
                     filter(allsig, ID=="sig_can")$p2a_2b)
ks_p2a_2b

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

### build contigency table for fisher's exact test 
library(GeneOverlap)

go.obj_comp_historical <- newGeneOverlap(listB = compstat_results$ensembl_gene_id,
                                         listA = sig_can_wg_paml_genes$ensembl_gene_id,
                                         genome.size = length(compstat_paml_universe))
go.obj_comp_historical <- testGeneOverlap(go.obj_comp_historical)
print(go.obj_comp_historical)


