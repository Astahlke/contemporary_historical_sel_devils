library(stringr)
library(dplyr)

### run your own hypergeometric test

#### make gene lists: foreground-background ####
## are contemporary candidates enriched for some pathway?
## contemporary candidates against all targeted loci

## are contemporary candidates, also under historical positive selection, enriched?
## contemporary_historical_candidates against targted loci that were tested in PAML

## are gwb candidates (not under contemporary selection) enriched?
## Sig PAML GWB against all tested in PAML

contemporary_universe <- read.delim("all_composite_stat_results/2020-10-22/results/annotation_top100.0/composite.snps.everything.top.genes.100000bp.txt",
                       header = FALSE, col.names = c("beg_scaffold", 
                                                     "beg_chr_start", "beg_chr_stop", "compstat", 
                                                     "end_scaffold", "end_chr_start", "end_chr_stop", 
                                                     "annotation"))

head(contemporary_universe)
contemporary_universe <- str_split(contemporary_universe$annotation, ":", simplify = TRUE)[,1]

contemporary_sel <-read.delim("all_composite_stat_results/2020-10-22/results/annotation_top1.0/composite.snps.everything.top.genes.100000bp.txt",
                              header = FALSE, col.names = c("beg_scaffold", 
                                                            "beg_chr_start", "beg_chr_stop", "compstat", 
                                                            "end_scaffold", "end_chr_start", "end_chr_stop", 
                                                            "annotation"))
contemporary_sel <- str_split(contemporary_sel$annotation, ":", simplify = TRUE)[,1]

length(contemporary_sel) # 252 snps
contemporary_sel <- unique(contemporary_sel)
length(contemporary_sel) # 186 genes

paml_results <- read.csv("all_paml_results.csv")
head(paml_results)

## Are contemporary candidates under positive selection enriched for function? 
compstat_paml_universe <- semi_join(paml_results, data.frame(ensembl_gene_id = contemporary_sel))
# universe
compstat_paml_universe$ID <- "compstat_paml_universe"

# foreground
sig_historical_contemporary <- dplyr::filter(compstat_paml_universe, pval < .05)
sig_historical_contemporary$ID <- "sig_historical_contemporary"

## Are genes of the GWB under positive selection enriched for function?
historical <- anti_join(paml_results, data.frame(ensembl_gene_id = contemporary_sel))
sig_historical <- dplyr::filter(historical, pval < .05)
sig_historical$ID <- "sig_historical"
## these are a subset of the genome-wide background paml universe
paml_universe <- paml_results$ensembl_gene_id


##### get hugo ids ######

library(biomaRt)
mart <- 'ENSEMBL_MART_ENSEMBL'
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'sharrisii_gene_ensembl')

#### contemporary_universe ####
genelist <- contemporary_universe
hugo_ids = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                 filters = 'ensembl_gene_id', 
                 values = genelist, 
                 mart = mart)
head(hugo_ids)

hugo_ids[hugo_ids == ""] <- NA
hugo_ids <- hugo_ids[-which(is.na(hugo_ids[,'hgnc_symbol'])),]
hugo_ids[which(duplicated(hugo_ids$hgnc_symbol)),]
contemporary_universe_hugo <- hugo_ids
nrow(contemporary_universe_hugo)/length(contemporary_universe) #19%
write.csv(hugo_ids$mdomestica_homolog_ensembl_gene, "mdomestica_allPAMLgeneIDs.csv", col.names = NULL, row.names = FALSE)
#### contemporary_sel_hugo ####
genelist <- contemporary_sel
hugo_ids = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                 filters = 'ensembl_gene_id', 
                 values = genelist, 
                 mart = mart)

hugo_ids[hugo_ids == ""] <- NA
hugo_ids <- hugo_ids[-which(is.na(hugo_ids[,'hgnc_symbol'])),]
hugo_ids[which(duplicated(hugo_ids$hgnc_symbol)),]
contemporary_sel_hugo <- hugo_ids
nrow(contemporary_sel_hugo)/length(contemporary_sel) #60%

contemporary_hugo <- rbind(contemporary_universe_hugo, contemporary_sel_hugo)
nrow(contemporary_sel_hugo)
# write.csv(contemporary_sel_hugo, "contemporary_sel_hugo.csv", quote = FALSE, row.names = FALSE)

#### paml_universe ####
genelist <- paml_universe
hugo_ids = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                 filters = 'ensembl_gene_id', 
                 values = genelist, 
                 mart = mart)

hugo_ids[hugo_ids == ""] <- NA
hugo_ids <- hugo_ids[-which(is.na(hugo_ids[,'hgnc_symbol'])),]
hugo_ids[which(duplicated(hugo_ids$hgnc_symbol)),]
paml_universe_hugo <- hugo_ids
nrow(paml_universe_hugo)/length(paml_universe) #76%

#### sig_historical: hist_sel_wg_hugo ####
genelist <- sig_historical$ensembl_gene_id
hugo_ids = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                 filters = 'ensembl_gene_id', 
                 values = genelist, 
                 mart = mart)

hugo_ids[hugo_ids == ""] <- NA
hugo_ids <- hugo_ids[-which(is.na(hugo_ids[,'hgnc_symbol'])),]
hugo_ids[which(duplicated(hugo_ids$hgnc_symbol)),]
hist_sel_wg_hugo <- hugo_ids
head(hugo_ids)
write.csv(hist_sel_wg_hugo, "historical_sel_hugo.csv", quote = FALSE, row.names = FALSE)

#### sig_historical_contemporary: hist_sel_can_hugo ####
genelist <- sig_historical_contemporary$ensembl_gene_id
hugo_ids = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                 filters = 'ensembl_gene_id', 
                 values = genelist, 
                 mart = mart)

hugo_ids[hugo_ids == ""] <- NA
hugo_ids <- hugo_ids[-which(is.na(hugo_ids[,'hgnc_symbol'])),]
hugo_ids[which(duplicated(hugo_ids$hgnc_symbol)),]
hist_sel_can_hugo <- hugo_ids
head(hist_sel_can_hugo)
write.csv(hist_sel_can_hugo, "contemporaryhistorical_sel_hugo.csv", quote = FALSE, row.names = FALSE)

#### historical no contemporary ####
genelist <- historical$ensembl_gene_id
hugo_ids = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                 filters = 'ensembl_gene_id', 
                 values = genelist, 
                 mart = mart)

hugo_ids[hugo_ids == ""] <- NA
hugo_ids <- hugo_ids[-which(is.na(hugo_ids[,'hgnc_symbol'])),]
hugo_ids[which(duplicated(hugo_ids$hgnc_symbol)),]
historical_nocontemporary <- hugo_ids

#### compstat_paml_universe ####
genelist <- compstat_paml_universe$ensembl_gene_id
hugo_ids = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                 filters = 'ensembl_gene_id', 
                 values = genelist, 
                 mart = mart)

hugo_ids[hugo_ids == ""] <- NA
hugo_ids <- hugo_ids[-which(is.na(hugo_ids[,'hgnc_symbol'])),]
hugo_ids[which(duplicated(hugo_ids$hgnc_symbol)),]
compstat_paml_universe_hugo <- hugo_ids
#####

# nrow(hist_nosel_can_hugo)/length(hist_nosel_can) # 71%
# nrow(hist_sel_can_hugo)/length(hist_sel_can) # 46%
# nrow(paml_universe_hugo)/length(paml_universe) # 64%
# nrow(hist_sel_wg_hugo)/length(hist_sel_wg) # 55%
# nrow(contemporary_sel_hugo)/nrow(contemporary_sel) # 56%
# nrow(contemporary_universe_hugo)/length(contemporary_universe) # 16% 

# hist_nosel_hugo$ID <- "hist_nosel_hugo"
# hist_sel_wg_hugo$ID <- "hist_sel_wg_hugo"
# contemporary_universe_hugo$ID <- "contemporary_universe_hugo"
# contemporary_sel_hugo$ID <- "contemporary_sel_hugo"
# hist_sel_can_hugo$ID  <- "hist_sel_can_hugo"
# contemporary_hist_universe_hugo$ID <- "contemporary_hist_universe_hugo"
#  
# all_hugo_sets <- rbind(hist_nosel_hugo, 
#                        hist_nosel_hugo, 
#                        hist_sel_wg_hugo,
#                        hist_sel_can_hugo,
#                        contemporary_universe_hugo, 
#                        contemporary_sel_hugo, 
#                        contemporary_hist_universe_hugo)
# 

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

#### contemporary_universe ####
genelist <- contemporary_universe
hugo_ids = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                 filters = 'ensembl_gene_id', 
                 values = genelist, 
                 mart = mart)

hugo_ids[hugo_ids == ""] <- NA
hugo_ids <- hugo_ids[-which(is.na(hugo_ids[,'hgnc_symbol'])),]
hugo_ids[which(duplicated(hugo_ids$hgnc_symbol)),]
contemporary_universe_hugo <- hugo_ids


#### contemporary_sel ####
genelist <- contemporary_sel
hugo_ids = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                 filters = 'ensembl_gene_id', 
                 values = genelist, 
                 mart = mart)

hugo_ids[hugo_ids == ""] <- NA
hugo_ids <- hugo_ids[-which(is.na(hugo_ids[,'hgnc_symbol'])),]
hugo_ids[which(duplicated(hugo_ids$hgnc_symbol)),]
contemporary_sel_hugo <- hugo_ids


#### paml_universe ####
genelist <- paml_universe
hugo_ids = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                 filters = 'ensembl_gene_id', 
                 values = genelist, 
                 mart = mart)

hugo_ids[hugo_ids == ""] <- NA
hugo_ids <- hugo_ids[-which(is.na(hugo_ids[,'hgnc_symbol'])),]
hugo_ids[which(duplicated(hugo_ids$hgnc_symbol)),]
paml_universe_hugo <- hugo_ids


#### hist_sig_wg ####
genelist <- hist_sig_wg
hugo_ids = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                 filters = 'ensembl_gene_id', 
                 values = genelist, 
                 mart = mart)

hugo_ids[hugo_ids == ""] <- NA
hugo_ids <- hugo_ids[-which(is.na(hugo_ids[,'hgnc_symbol'])),]
hugo_ids[which(duplicated(hugo_ids$hgnc_symbol)),]
hist_sig_wg_hugo <- hugo_ids
head(hist_sig_wg_hugo)

#### hist_sig_can ####
genelist <- hist_sig_can
hugo_ids = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                 filters = 'ensembl_gene_id', 
                 values = genelist, 
                 mart = mart)

hugo_ids[hugo_ids == ""] <- NA
hugo_ids <- hugo_ids[-which(is.na(hugo_ids[,'hgnc_symbol'])),]
hugo_ids[which(duplicated(hugo_ids$hgnc_symbol)),]
hist_sig_can_hugo <- hugo_ids

#####

all_query_background <- c()

query_name <- "contemporary_sel_hugo"
background_name <- "contemporary_universe_hugo"

query <- get(query_name)$hgnc_symbol
my_background <- get(background_name)$hgnc_symbol
m <- length(query)
y <- c()

for (j in 1:length(MSigDB_ofinterest)){
  # j <- length(MSigDB_ofinterest)
  db_name <- MSigDB_ofinterest[j]
  # [1] "HALLMARK"                 
  # [3] "C2_CURATED"               
  # [5] "C4_COMPUTATIONAL"         
  # [8] "C7_IMMUNOLOGIC_SIGNATURES"
  db <- MSigDB[[db_name]]
  for (i in 1:length(db)){
    cat(db_name, query_name, background_name, names(db[i]), sep = "\n") 
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
write.csv(y, file = paste(query_name, background_name, "raw_geneset_enrichment_.csv", sep="_"),
          row.names = FALSE,
          quote = FALSE)
all_query_background <- rbind(y, all_query_background)


## candidates for contemporary and historical
query_name <- "hist_sel_wg_hugo" 
background_name <- "paml_universe_hugo"

query <- get(query_name)$hgnc_symbol
my_background <- get(background_name)$hgnc_symbol
m <- length(query)
y <- c()
for (j in 1:length(MSigDB_ofinterest)){
  db_name <- MSigDB_ofinterest[j]
  db <- MSigDB[[db_name]]

    for (i in 1:length(db)){
    cat(db_name, query_name, background_name, names(db[i]), sep = "\n") 
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
write.csv(y, file = paste(query_name, background_name, "raw_geneset_enrichment_.csv", sep="_"),
          row.names = FALSE,
          quote = FALSE)
all_query_background <- rbind(y, all_query_background)

## candidates for contemporary and historical
query_name <- "hist_sel_can_hugo"
background_name <- "compstat_paml_universe_hugo"

query <- get(query_name)$hgnc_symbol
my_background <- get(background_name)$hgnc_symbol
m <- length(query)
y <- c()
for (j in 1:length(MSigDB_ofinterest)){
  db_name <- MSigDB_ofinterest[j]
  db <- MSigDB[[db_name]]
  
  for (i in 1:length(db)){
    cat(db_name, query_name, background_name, names(db[i]), sep = "\n") 
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
write.csv(y, file = paste(query_name, background_name, "raw_geneset_enrichment_.csv", sep="_"),
          row.names = FALSE,
          quote = FALSE)
all_query_background <- rbind(y, all_query_background)

all_query_background %>%
  count(query_name)

all_query_background %>% 
  mutate(p_adj = p.adjust(pval, method = "BH")) %>%
  dplyr::filter(p_adj < .01 & N > 10) %>% 
  count(query_name)

write.csv(all_query_background,  file = paste("all_query_background_raw_geneset_enrichment_.csv", sep="_"),
          row.names = FALSE,
          quote = FALSE)

sig_all_query_background <- all_query_background %>% 
  mutate(p_adj = p.adjust(pval, method = "BH")) %>%
  dplyr::filter(p_adj < .01 & N > 10)

shared_hist_contemporary <- sig_all_query_background[which(duplicated(sig_all_query_background$geneset_name)),]
nrow(shared_hist_contemporary)


geneset_name <- data.frame(geneset_name=allresults[duplicated(allresults$geneset_name),]$geneset_name)
write.csv(semi_join(allresults, geneset_name), 
          file = "shared_geneset_enrichment.csv",
          row.names = FALSE,
          quote = FALSE)


#### permutation test ####
### do contemporary candidates share geneset enrichment more often than expected by chance? ###

## filter for these GWB genesets first
sig_GWB_genesets <- sig_all_query_background %>% 
  dplyr::filter(query_name == "hist_sel_wg_hugo") %>% 
  dplyr::select(geneset_name) 

sig_GWB_genesets <- as.character(sig_GWB_genesets$geneset_name)

subset_genesets <- function(x) {  
  sg <- x[which(names(x) %in% sig_GWB_genesets)]
  return(sg)
}

subset_sigGWB_genesets <- lapply(MSigDB, subset_genesets)
length(sig_GWB_genesets)
sum(unlist(lapply(subset_sigGWB_genesets, length)))
      
       
### 1. re-draw N contemporary candidates from contemporary universe
### 1.a With and without replacement
### 2. Test for gene set enrichment 
### 3. measure the number (and strength?) of shared geneset enrichment 

Ncandidate_hugo <- nrow(contemporary_sel_hugo)
sample_n( contemporary_hugo, Ncandidate_hugo)

## without replacement 

x <- sample.int(nrow(contemporary_hugo), replace = FALSE)

x_chunks <- split(x, ceiling(seq_along(x)/nrow(contemporary_sel_hugo)))
query_list <- sapply(x_chunks, function(z) contemporary_hugo$hgnc_symbol[z])
background_name <- "contemporary_background"

y <- c()

for (Q in 1: length(x_chunks)) {
  query_name <- paste0("permutated_contemporary_", Q) 
  query <- query_list[Q] 
  my_background <- contemporary_hugo$hgnc_symbol
  m <- length(query)
  for (j in 1:length(MSigDB_ofinterest)){
    db_name <- MSigDB_ofinterest[j]
    cat(db_name, query_name, "\n") 
    db <- subset_sigGWB_genesets[[db_name]]
    for (i in 1:length(db)){
    if(length(subset_sigGWB_genesets[[db_name]])<1){
        cat("skip", db_name, "\n") 
      next }
      cat(query_name, names(db[i]), i, "of", length(db), i/length(db),"\n")
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
      if( length(x_intersect) > 0){
        cat(names(db[i]), "intersect", x_intersect, "\n")
      }
      y <- rbind(y, probabilities)
    }
  }
}
head(arrange(y, desc(prop_overlap)))
without_replacement_gsea <- y

permuted_overlap <- filter(without_replacement_gsea, N > 10, pval < .001) %>%
  group_by(query_name) %>%
  count() %>%
  data.frame()

filter(without_replacement_gsea, N > 10, pval < .001) %>%
  group_by(db_name) %>%
  count() %>%
  data.frame()


### do Mark's GWAS candidates overlap with historical candidates more often than expected by chance? ###
sig_GWB_genesets <- read.csv("sig_GWB_genesets.csv", stringsAsFactors = FALSE, header = FALSE)
sig_GWB_genesets <- as.character(sig_GWB_genesets$V1)


subset_genesets <- function(x) {  
  sg <- x[which(names(x) %in% sig_GWB_genesets)]
  return(sg)
}

subset_sigGWB_genesets <- lapply(MSigDB, subset_genesets)
length(sig_GWB_genesets)
sum(unlist(lapply(subset_sigGWB_genesets, length)))


### 1. re-draw N contemporary candidates from contemporary universe
### 1.a With and without replacement
### 2. Test for gene set enrichment 
### 3. measure the number (and strength?) of shared geneset enrichment 

MSigDB_ofinterest <- c("HALLMARK", "C2_CURATED", "C4_COMPUTATIONAL", "C6_ONCOGENIC_SIGNATURES", "C7_IMMUNOLOGIC_SIGNATURES")


Ncandidate_hugo <- nrow(contemporary_sel_hugo)
sample_n( contemporary_hugo, Ncandidate_hugo)

## without replacement 

x <- sample.int(nrow(contemporary_hugo), replace = FALSE)

x_chunks <- split(x, ceiling(seq_along(x)/nrow(contemporary_sel_hugo)))
query_list <- sapply(x_chunks, function(z) contemporary_hugo$hgnc_symbol[z])
background_name <- "contemporary_background"

y <- c()

for (Q in 1: length(x_chunks)) {
  query_name <- paste0("permutated_contemporary_", Q) 
  query <- query_list[,Q] 
  my_background <- contemporary_hugo$hgnc_symbol
  m <- length(query)
  for (j in 1:length(MSigDB_ofinterest)){
    db_name <- MSigDB_ofinterest[j]
    cat(db_name, query_name, "\n") 
    # [1] "HALLMARK"                 
    # [3] "C2_CURATED"               
    # [5] "C4_COMPUTATIONAL"         
    # [7] "C6_ONCOGENIC_SIGNATURES"  
    # [8] "C7_IMMUNOLOGIC_SIGNATURES"
    db <- subset_sigGWB_genesets[[db_name]]
    for (i in 1:length(db)){
      if(length(subset_sigGWB_genesets[[db_name]])<1){
        cat("skip", db_name, "\n") 
        next }
      cat(query_name, names(db[i]), i, "of", length(db), i/length(db),"\n")
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
      if( length(x_intersect) > 0){
        cat(names(db[i]), "intersect", x_intersect, "\n")
      }
      y <- rbind(y, probabilities)
      cat("top geneset overlap", as.character(arrange(y, desc(pval))$geneset_name[1]),  "\n")
    }
  }
}
head(arrange(y, desc(prop_overlap)))
without_replacement_gsea <- y

permuted_overlap <- filter(without_replacement_gsea, N > 10, pval < .001) %>%
  group_by(query_name) %>%
  count() %>%
  data.frame()

filter(without_replacement_gsea, N > 10, pval < .001) %>%
  group_by(db_name) %>%
  count() %>%
  data.frame()




library(ggplot2)

table(permuted_overlap$n > 22)
mean(permuted_overlap$n)
sd(permuted_overlap$n)


load("permutation1000/sig_shared_results.RData")
sig_shared_results[[1]]
permuted1000_overlap <- data.frame(n = unlist(lapply(sig_shared_results, nrow)))

?IQR

bw <- 2 * IQR(permuted1000_overlap$n) / length(permuted1000_overlap$n)^(1/3)
library(ggplot2)
ggplot(data = permuted1000_overlap) +
  geom_histogram(position = "dodge", aes(n), binwidth = bw) +
  geom_vline(xintercept =22, color = "red") +
  annotate("text", x=35, y= 115,
           label="Observed = 22", color = "red") +
  # ggtitle("Number of shared gene set overlaps for 1000 permutations of contemporary set \n with replacement") +
  theme_bw()
ggsave("genesetoverlaps_permutation_hist.pdf", width = 7, units = "in", height = 6, dpi = 300)

length(which(permuted1000_overlap$n < 22))/1000


#the proportion of permutations with larger difference
(sum(permuted1000_overlap$n < abs(22)) + 1) / (length(permuted1000_overlap$n ) + 1)

reps <- 1000
query_list <- vector("list", length = reps)
## with replacement 
for (rep in 1:reps){
  query_list[[rep]]  <- sample_n(contemporary_hugo, size = nrow(contemporary_sel_hugo))$hgnc_symbol
}

Q <- NULL
for (Q in 1: length(query_list)) {
  query_name <- paste("permutated_withreplacement_contemporary_", Q) 
  query <- query_list[[Q]] 
  my_background <- contemporary_hugo$hgnc_symbol
  m <- length(query)
  # y <- c()
  for (j in 1:length(MSigDB_ofinterest)){
    db_name <- MSigDB_ofinterest[j]
    # [1] "HALLMARK"                 
    # [3] "C2_CURATED"               
    # [5] "C4_COMPUTATIONAL"         
    # [7] "C6_ONCOGENIC_SIGNATURES"  
    # [8] "C7_IMMUNOLOGIC_SIGNATURES"
    db <- MSigDB[[db_name]]
    for (i in 1:length(db)){
      cat(db_name, query_name, background_name, names(db[i]), sep = "\n") 
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
}
with_replacement_gsea <- y



query_name <- "contemporary_sel_hugo"
background_name <- "contemporary_universe_hugo"

## without replacement
out <- vector("list", nrow(contemporary_hugo) /nrow(contemporary_sel_hugo))
for (i in seq_along(means)) {
  n <- sample(100, 1)
  out[[i]] <- rnorm(n, means[[i]])
}
str(out)

## with replacement
out <- vector("list", 1000)
for (i in seq_along(means)) {
  n <- sample(100, 1)
  out[[i]] <- rnorm(n, means[[i]])
}
str(out)


query <- get(query_name)$hgnc_symbol
my_background <- get(background_name)$hgnc_symbol
m <- length(query)
y <- c()
for (j in 1:length(MSigDB_ofinterest)){
  # j <- length(MSigDB_ofinterest)
  db_name <- MSigDB_ofinterest[j]
  # [1] "HALLMARK"                 
  # [3] "C2_CURATED"               
  # [5] "C4_COMPUTATIONAL"         
  # [6] "C5_GENE_ONTOLOGY"         
  # [7] "C6_ONCOGENIC_SIGNATURES"  
  # [8] "C7_IMMUNOLOGIC_SIGNATURES"
  db <- MSigDB[[db_name]]
  for (i in 1:length(db)){
    cat(db_name, query_name, background_name, names(db[i]), sep = "\n") 
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




