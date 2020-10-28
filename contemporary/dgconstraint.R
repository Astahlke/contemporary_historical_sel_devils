# This script calculates a constraint/repeatabulity statistic using the column
# labeled "adjp" in a set of tsv files.
#
# DETAILS
#
# Yeaman, Gerstein, Hodgins, Whitlock. Quantifying how constraints limit the diversity of viable routes to adaptation.=
# Plos Genetics. 2018. https://doi.org/10.1371/journal.pgen.1007717. 
#   For each statistic (here, af & mm), calculate c-score and p-value
#
# NOTE
#
#   The expected columns are:
#
#       scaffold    the scaffold / contig / chromosome
#       start       (0-based)
#       adjp        the p-value-like quantity to use
#
#   Other columns are ok, but will be ignored.
#

#==============================================================================#

# devtools::install_github("samyeaman/dgconstraint", build_vignettes = TRUE)
# browseVignettes(package = "dgconstraint")
library("dgconstraint")


options('scipen'=50)

scaleup <- function(x, na.rm = TRUE) 
  (-log10(x)) 

infiles <- list.files("all_composite_stat_results/2020-10-22/results", 
                      pattern =  glob2rx("*adjusted.tsv"), full.names = TRUE)

infiles
indata = vector('list', length=length(infiles))
for(i in seq_along(infiles)) {
    d = read.csv(infiles[i], sep='\t', header=TRUE, as.is=TRUE)
    indata[[i]] = d[, c('scaffold', 'start', 'adjp')]
}

# First merge the datasets
snps = unique(unlist(lapply(indata,
                            function(x) paste0(x[, 'scaffold'], '-', x[, 'start']))))
id_table  = data.frame('scaffold'=gsub('-.+', '', snps),
                       'start'=as.numeric(gsub('.+-', '', snps)))

data = vector('list', length=length(indata))
for(i in seq_along(data)) {
    tmp = merge(id_table, indata[[i]][, c('scaffold', 'start', 'adjp')],
                all=TRUE, stringsAsFactors=FALSE)
    data[[i]] = tmp[order(tmp[, 'scaffold'], tmp[, 'start']), ]
}

# Check to make sure the merging step worked
scaffolds = data[[1]][, 'scaffold']
positions = data[[1]][, 'start']
if(length(data) > 1) {
    passing = sapply(data[-1], function(x) {
                     all(x[, 'scaffold'] == scaffolds) &&
                     all(x[, 'start'] == positions)
                    })
    if(!all(passing))
        stop('Merging did not work')
}

names(data) <- list.files("all_composite_stat_results/2020-10-22/results", 
                            pattern =  glob2rx("*adjusted.tsv"), full.names = FALSE)

## somehow need to account for LD? 
df <- do.call(cbind, lapply(data, as.data.frame))
df[,1] <- paste0(df$Fentonbury.afchange.adjusted.tsv.scaffold, "_", df$Fentonbury.afchange.adjusted.tsv.start) # replace scaffold with full SNP
colnames(df)[1] <- "SNP"

df <-  df %>%
  select("SNP", contains("adjp")) %>%
  mutate_at(vars(-("SNP")), scaleup)

mm_matrix <- df %>%
  select(contains("mm")) %>%
  as.matrix()
sum(rowSums(is.na(mm_matrix))==0)/nrow(mm_matrix) # 71% have no missing data

af_matrix <- df %>%
  select(contains("af")) %>%
  as.matrix() 
sum(rowSums(is.na(af_matrix))==0)/nrow(af_matrix) # only 43% have no missing

#### no adjustment
c_score_af <- pairwise_c_chisq(af_matrix, num_permute = 10000, na.rm = T) # 3.725592
p_value_af <- allwise_p_chisq(af_matrix, num_permute = 10000, na.rm = T) # 1e-04
# Warning message: In allwise_p_chisq(af_matrix, num_permute = 10000, na.rm = T) Note: estimated p-value is lower than 1/(number of permutations)

c_score_mm <- pairwise_c_chisq(mm_matrix, num_permute = 10000, na.rm = T) # 4.974022
p_value_mm <- allwise_p_chisq(mm_matrix, num_permute = 10000, na.rm = T) # 1e-04


#### with adjustment
# adjust continuous data by adding g0 (1 − q) new entries to the dataset by 
# randomly sampling genes with replacement from the existing dataset, 
# and then applying Eq 2 (allwise_p_chisq) to this
# extended set.
# g0 = total number of genes in genome
# q = proportion that can be accurately estimated
## proportion  accurately estimated = # contemporary universe / # all genes possible
contemporary_universe <- read.delim("all_composite_stat_results/2020-10-22/results/annotation_top100.0/composite.snps.everything.top.genes.100000bp.txt",
                                    header = FALSE, col.names = c("beg_scaffold", 
                                                                  "beg_chr_start", "beg_chr_stop", "compstat", 
                                                                  "end_scaffold", "end_chr_start", "end_chr_stop", 
                                                                  "annotation"))

head(contemporary_universe)
contemporary_universe <- unique(str_split(contemporary_universe$annotation, ":", simplify = TRUE)[,1])

## all genes
library("biomaRt")
ensembl <- useMart("ensembl")
ensembl = useDataset("sharrisii_gene_ensembl",mart=ensembl)
listAttributes(ensembl)
allgenes <- getBM(attributes= 'ensembl_gene_id', 
      mart = ensembl)
g0 = nrow(allgenes) # g0 = total number of genes in genome

length(intersect(contemporary_universe, allgenes$ensembl_gene_id)) == length(contemporary_universe)
length(contemporary_universe)/nrow(allgenes) # 29 %
q = length(contemporary_universe)/nrow(allgenes) # q = proportion that can be accurately estimated 
### (or should this be number with no NA?)

n_resamps <- g0*(1-q)

af_mat_noNA <- af_matrix[complete.cases(af_matrix),]
af_mat_resamped <- as.data.frame(af_mat_noNA) %>%
  slice_sample(n = n_resamps, replace = T) %>%
  rbind(., af_mat_noNA) %>%
  as.matrix()

c_score_af <- pairwise_c_chisq(af_mat_resamped, num_permute = 10000, na.rm = F)
c_score_af  # 23.8
p_value_af <- allwise_p_chisq(af_mat_resamped, num_permute = 10000, na.rm = F) 
p_value_af # 1e-04


mm_mat_noNA <- mm_matrix[complete.cases(mm_matrix),]
mm_mat_resamped <- as.data.frame(mm_mat_noNA) %>%
  slice_sample(n = n_resamps, replace = T) %>%
  rbind(., mm_mat_noNA) %>%
  as.matrix()

c_score_mm <- pairwise_c_chisq(mm_mat_resamped, num_permute = 10000, na.rm = F)
c_score_mm  # 13.70
p_value_mm <- allwise_p_chisq(mm_mat_resamped, num_permute = 10000, na.rm = T) 
p_value_mm # ???

# convert this to binary data with some threshold
# goal:  estimate_pa(af_matrix, ndigits = 4, show.plot = T, na.rm = T)

tp <- .01

sort_mm <- sort(as.vector(mm_mat_resamped), decreasing = T)
top_mm <- sort_mm[1:floor(length(sort_mm) * tp)]
thresh <- min(top_mm)
mm_binary <- +(mm_mat_resamped >= thresh)

prop_adaptive_mm <- estimate_pa(mm_binary, ndigits = 4, show.plot = T, na.rm = T) # .34
prop_adaptive_mm # .3483

sort_af <- sort(as.vector(af_mat_resamped), decreasing = T)
top_af <- sort_af[1:floor(length(sort_af) * tp)]
thresh <- min(top_af)

af_binary <- +(af_mat_resamped >= thresh)
prop_adaptive_af <- estimate_pa(af_binary, ndigits = 4, show.plot = T, na.rm = T)
prop_adaptive_af #.164

## these results still leave much to be desired. coop and lee 2017; buffalo and coop 2020; silas tittes
## 
