
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

######## CODEML - Branch site models #######

colnames <- c("ensembl_gene_id", "likelihood", "treelength", "p0", 
              "b0", "f0", "p1", "b1", "f1", "p2", "b2", "f2", "p3", "b3", "f3")
interest_cols <- c(c("ensembl_gene_id", "likelihood")) # "p3" don't think I need this; or THIS

paml_posi <- read.table(file = "paml_out/paml_posiBS_wgcan_koala_111617.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
paml_posi[1] <-  strtrim(paml_posi[,1], 18)
colnames(paml_posi) <- colnames
paml_posi_LRT <- paml_posi[, interest_cols]


paml_null <- read.table(file = "paml_out/paml_nullBS_wgcan_koala_111617.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
paml_null[1] <-  strtrim(paml_null[,1], 18)
colnames(paml_null) <- colnames
paml_null_LRT <- paml_null[, interest_cols]

paml_LRT <- merge(paml_posi_LRT, paml_null_LRT, by = "ensembl_gene_id", suffixes = c(".posi",".null")) %>% 
  mutate(likelihood = 2*(likelihood.posi-likelihood.null))

paml_LRT$pval <- pchisq(as.numeric(paml_LRT[,"likelihood"]), 1, lower.tail = FALSE)
paml_LRT$FDR <- p.adjust(paml_LRT$pval, method = "BH")

sig_paml_genes <- subset(paml_LRT, FDR < 0.05) %>% 
  dplyr::select(ensembl_gene_id, FDR) %>% 
  merge(paml_posi, .,  by = "ensembl_gene_id") %>% 
  dplyr::filter(f3 > 1) %>% 
  mutate(p2a_2b= p2+p3)

nrow(sig_paml_genes)/nrow(paml_LRT) 

## What do I actually need to carry forward? Probably the full set of LRT and the subset of significant
write.csv(sig_paml_genes, "sig_paml_results_21-1-13.csv", quote = F)
write.csv(paml_LRT, "paml_LRT_21-1-13.csv", quote = F)