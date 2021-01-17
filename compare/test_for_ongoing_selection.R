### Test for ongoing selection: contemporary & historical
### This script takes the results from contemporary and historical tests
### and compares the results through several ways:
### 1. which genes are in both sets? straight overlap in sets
###       output: sig_contemporary_historical_paml.csv
### 2. how do the distributions of dN/dS and proportion of sites compare?
###   a. generate plots Figure3_PAMLresults.jpg
###   b. AD test
###   c. KS test
### 3. gene overlap with contingency tables
### 4. MSigDB overlap - do gene sets share a similar function?

library(tidyverse)
library(ggplot2)
library(cowplot)

### Functions
## get hgnc ids
genelist_hugo <- function(genelist) {
  
  hugo_ids = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                   filters = 'ensembl_gene_id', 
                   values = genelist, 
                   mart = mart)
  
  hugo_ids[hugo_ids == ""] <- NA
  hugo_ids <- hugo_ids[-which(is.na(hugo_ids[,'hgnc_symbol'])),]
  cat(nrow(hugo_ids)/length(genelist), "of genelist has hugo annotations")
  
  return(hugo_ids)
}

### run gsea
gsea <- function(query_name, background_name) {
  
  all_query_background <- c()
  query <- get(query_name)$hgnc_symbol
  my_background <- get(background_name)$hgnc_symbol
  m <- length(query)
  y <- c()
  
  for (j in 1:length(MSigDB_ofinterest)){
    
    db_name <- MSigDB_ofinterest[j]
    db <- MSigDB[[db_name]]
    for (i in 1:length(db)){
      cat(db_name, names(db[i]), sep = "\n") 
      background <- intersect(my_background, db[[i]])
      N <- length(background)
      x_intersect <- intersect(query, background)
      gsea <- data.frame(db_name=db_name,
                         background_name = background_name,
                         query_name = query_name,
                         geneset_name = names(db[i]),
                         N = N, # gene universe; all balls in the urn
                         m = m, # genes under selection; all white balls in urn
                         x = length(x_intersect), # genes under sel and in gene set
                         n = N - length(x_intersect), # genes in background, not under sel; black balls
                         prop_overlap = length(x_intersect)/length(background),
                         overlap = paste(x_intersect, collapse = ";"))
      # Use the dhyper built-in function for hypergeometric density
      probabilities <- mutate(gsea, pval=(dhyper(c(0:gsea$x), 
                                                 gsea$m, 
                                                 gsea$n, 
                                                 gsea$x, log = FALSE)[1]))
      y <- rbind(y, probabilities)
    }
  }
  y <- y %>% 
    mutate(p_adj = p.adjust(pval, method = "BH")) %>%
    dplyr::filter(p_adj < .05)
  
  return(y)
}

### Set up to run locally ~/Code/devils/
getwd()

#### contemporary selection ####
contemporarydir <- "all_composite_stat_results/2020-10-22/results/"

## main compstat results
compstat_results <- read.csv(paste0(contemporarydir, "annotation_top1.0/composite.snps.everything.top.genes.100000bp.txt"), 
                             row.names = NULL, sep = "\t", header = FALSE)
head(compstat_results)
compstat_results$ensembl_gene_id <- strtrim((compstat_results$V8), 18) 
length(unique(compstat_results$ensembl_gene_id)) 


#### historical selection ####
allpamlresults <- read.csv("historical/paml_LRT_21-1-13.csv")
sig_paml <- read.csv("historical/sig_paml_results_21-1-13.csv")

nrow(sig_paml)
nrow(sig_paml)/nrow(allpamlresults)

sig_can_paml <- semi_join(sig_paml, compstat_results, by ="ensembl_gene_id") %>% 
  mutate(ID = "recurrent")
nrow(sig_can_paml)

write.csv(sig_can_paml, "sig_contemporary_historical_paml.csv", row.names = F, quote = F)

## these are under historical but not contemporary candidates
sig_wg_paml <- anti_join(sig_paml, compstat_results, by ="ensembl_gene_id") %>% 
  mutate(ID = "gwb")
sig_can_wg_paml <- rbind(sig_can_paml, sig_wg_paml)

sig_can_wg_paml$ID <- as.factor(sig_can_wg_paml$ID)
levels(sig_can_wg_paml$ID) <- levels(sig_can_wg_paml$ID)[c(2,1)]

mean(sig_can_paml$f3)
mean(sig_wg_paml$f3)

min(sig_can_paml$f3)
min(sig_wg_paml$f3)

max(sig_can_paml$f3)
max(sig_wg_paml$f3)

sd(sig_can_paml$f3)
sd(sig_wg_paml$f3)

### plots! ####
## magnitude of dn/ds
omega <- ggplot(sig_can_wg_paml, aes(x= ID, y = f3, color = ID)) +
  geom_violin(show.legend = F) +
  geom_jitter(shape=16, size = 3, position=position_jitter(0.1), alpha = .5, show.legend = F) +
  scale_fill_grey() +
  scale_color_manual(values=c("#b2182b", "black")) +
  xlab("") +  
  scale_x_discrete(guide = guide_axis(n.dodge=3), labels= c("Genome-wide Background", "Contemporary Candidates")) +
  ylab("dN/dS Estimates in Site Class 2") +
  theme(text = element_text(size = 18)) +
  theme_minimal()
omega

## proportion of sites
p2ab <- ggplot(sig_can_wg_paml, aes(x= ID, y = p2a_2b, color = ID)) +
  geom_violin(show.legend = F) +
  geom_jitter(shape=16, size = 3, position=position_jitter(0.1), alpha = .5, show.legend = F) +
  scale_fill_grey() +
  scale_color_manual(values=c( "#b2182b", "black")) +
  xlab("") +
  ylab("Proportion of Sites in Class 2") +
  scale_x_discrete(guide = guide_axis(n.dodge=3), labels= c( "Genome-wide Background", "Contemporary Candidates")) +
  # labs(title = "Proportion of sites under positive selection per gene") +
  # ggtitle(paste(candidategene_lists[i], "vs whole genome")) +
  theme(text = element_text(size = 18)) +
  theme_minimal() 
p2ab

plot_grid(omega, p2ab, labels = c('A', 'B'))
ggsave("Figure3_PAMLresults.jpg", width = 7, units = "in", height = 6, dpi = 300)

#### Compare contempororay and historical selection!! ####

### statistical tests comparing distributions ####
# install.packages("kSamples")
set.seed(42)

andersondarling_f3 <- kSamples::ad.test(sig_can_paml$f3, sig_wg_paml$f3, method = "exact")
andersondarling_f3

andersondarling_p2a2b <- kSamples::ad.test(sig_can_paml$p2a_2b, sig_wg_paml$p2a_2b, method = "exact")
andersondarling_p2a2b

#### compare the distributions of dN/dS estimates and propotions 
ks_f3 <- ks.test(sig_can_paml$f3, sig_wg_paml$f3)
ks_f3

ks_p2a_2b <- ks.test(sig_can_paml$p2a_2b, sig_wg_paml$p2a_2b)
ks_p2a_2b


#### make gene lists: foreground-background ####

### build contigency table for fisher's exact test 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GeneOverlap")


#### defining the multi-verse ####
contemporary_universe <- read.delim("all_composite_stat_results/2020-10-22/results/annotation_top100.0/composite.snps.everything.top.genes.100000bp.txt",
                                    header = FALSE, col.names = c("beg_scaffold", 
                                                                  "beg_chr_start", "beg_chr_stop", "compstat", 
                                                                  "end_scaffold", "end_chr_start", "end_chr_stop", 
                                                                  "annotation"))
contemporary_universe <- str_split(contemporary_universe$annotation, ":", simplify = TRUE)[,1]

compstat_paml_universe <- intersect(allpamlresults$ensembl_gene_id, contemporary_universe)

### Gene overlap
go.obj_comp_historical <- GeneOverlap::newGeneOverlap(listB = compstat_results$ensembl_gene_id,
                                         listA = sig_can_wg_paml$ensembl_gene_id,
                                         genome.size = length(compstat_paml_universe))
go.obj_comp_historical <- GeneOverlap::testGeneOverlap(go.obj_comp_historical)
print(go.obj_comp_historical)


##### get hugo ids ######

library(biomaRt)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'sharrisii_gene_ensembl')

### Pull gene IDs 
compstat_cans_hugo <- genelist_hugo(compstat_results$ensembl_gene_id)
compstat_uni_hugo <- genelist_hugo(contemporary_universe)

paml_wg_hugo <- genelist_hugo(sig_wg_paml$ensembl_gene_id)
paml_uni_hugo <- genelist_hugo(allpamlresults$ensembl_gene_id)





#####

# install.packages("devtools")
# devtools::install_github('oganm/MSigDB')
library(MSigDB)

# names(MSigDB)
# # [1] "HALLMARK"                 
# # [3] "C2_CURATED"               
# # [5] "C4_COMPUTATIONAL"         
# # [6] "C5_GENE_ONTOLOGY"         
# # [7] "C6_ONCOGENIC_SIGNATURES"  
# # [8] "C7_IMMUNOLOGIC_SIGNATURES"

MSigDB_ofinterest <- c("HALLMARK", "C2_CURATED", "C4_COMPUTATIONAL", "C7_IMMUNOLOGIC_SIGNATURES")

genesets.df <- data.frame(foreground = c("compstat_cans_hugo", "paml_wg_hugo"), 
                          background = c("compstat_uni_hugo", "paml_uni_hugo"),
                          geneset = c("contemporary", "historical"))

for(g in 1:nrow(genesets.df)){
  assign(x = paste0(genesets.df[g, "geneset"], "_gsea"), 
         gsea(genesets.df[g, "foreground"], genesets.df[g, "background"]))
  }

historical_gsea
contemporary_gsea

all_gsea <- rbind(historical_gsea, contemporary_gsea)

all_gsea %>%
  count(query_name)

all_gsea %>%
  group_by(query_name) %>% 
  dplyr::filter(N > 10) %>% 
  count(query_name)

all_gsea <- all_gsea %>%
  group_by(query_name) %>% 
  dplyr::filter(N > 10) 

## Is 10 an appropriate threshold? Seven genes was the maximum for contemporary
contemporary_gsea %>%
  count(N)

contemporary_gsea %>%
  arrange(N)

## with a minimum of 10, there are no shared gene sets
shared <- all_gsea[which(duplicated(all_gsea$geneset_name)),]
shared
nrow(shared)

### Nothing to write!
# write.csv(semi_join(all_gsea, 
#           data.frame(geneset_name = shared$geneset_name)), 
#           file = "shared_geneset_enrichment.csv",
#           row.names = FALSE,
#           quote = FALSE)


#### permutation test ####
### do contemporary candidates share geneset enrichment less/more often than expected by chance? ###

## filter for these GWB genesets first
sig_GWB_genesets <- historical_gsea %>% 
  dplyr::select("geneset_name") 
sig_GWB_genesets <- as.character(sig_GWB_genesets$geneset_name)

subset_genesets <- function(x) {  
  sg <- x[which(names(x) %in% sig_GWB_genesets)]
  return(sg)
}

subset_sigGWB_genesets <- lapply(MSigDB, subset_genesets)

Ncandidate_hugo <- nrow(compstat_cans_hugo)

reps <- 1000
query_list <- vector("list", length = reps)

## with replacement 
for (rep in 1:reps){
  query_list[[rep]]  <- sample(compstat_uni_hugo$hgnc_symbol, size = nrow(compstat_cans_hugo))
}

y <- c()
Q <- NULL
j <- c()
i <- c()
permuted_overlap <- c()

### the function I wrote above doesn't work as well with the permutation set up. So, back to the loop.

for (Q in 1:length(query_list)){
  query_name <- paste0("permutated_withreplacement_contemporary_", Q) 
  query <- query_list[[Q]] 
  m <- length(query)
  y <- c()
  
  for (j in 1:length(MSigDB_ofinterest)){
    db_name <- MSigDB_ofinterest[j]
    db <- MSigDB[[db_name]]
    for (i in 1:length(db)){
      cat(paste0(Q, " interation; ", j, " MSigDb of ",length(MSigDB_ofinterest), "; ",   i, "th datababse of ", length(db)), sep = "\n")
      background <- intersect(compstat_uni_hugo$hgnc_symbol, db[[i]])
      N <- length(background)
      x_intersect <- intersect(query, background)
      ## probably don't need to build this data frame - just the prop overlap 
      probabilities <- data.frame(Q = Q, # iteration 
                                  N = N, # gene universe; all balls in the urn
                                  m = m, # genes under selection; all white balls in urn
                                  x = length(x_intersect), # genes under sel and in gene set
                                  n = N - length(x_intersect), # genes in background, not under sel; black balls
                                  prop_overlap = length(x_intersect)/length(background)) %>% 
        mutate(pval=(dhyper(c(0:x), m, n, x, log = FALSE)[1])) %>% 
        dplyr::select(prop_overlap, pval, Q)
      y <- rbind(probabilities, y)
    }
  }
  permuted_overlap <- y %>%
    mutate(p_adj = p.adjust(pval, method = "BH")) %>%
    dplyr::filter(p_adj < .05 & N > 10) %>% 
    count() %>%
    data.frame(., Q) %>% 
    rbind(., permuted_overlap)
}

arrange(permuted_overlap, n)
mean(permuted_overlap$n)
sd(permuted_overlap$n)

nrow(shared)
table(permuted_overlap$n > nrow(shared))
## p-value
table(permuted_overlap$n > nrow(shared))/nrow(permuted_overlap)

### calculate an appropriate binwidth based on distribution
bw <- 2 * IQR(permuted_overlap$n) / length(permuted_overlap$n)^(1/3)

ggplot(data = permuted_overlap) +
  geom_histogram(position = "dodge", aes(n), binwidth = bw) +
  geom_vline(xintercept = nrow(shared), color = "red") +
  annotate("text", x=5, y= 115,
           label=paste0("Observed = ", nrow(shared)), color = "red") +
  theme_bw()
ggsave("Supplemental_geneset_hist.jpg", width = 7, units = "in", height = 6, dpi = 300)
save(permuted_overlap, file = "permuted_overlap.RData")
