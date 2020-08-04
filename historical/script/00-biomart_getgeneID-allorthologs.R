#### This generates a dataframe of random and candidate genes to later pull from in 01-getSeq
### inputs: lines 14-22; analysis_top_genes for gene IDs
### outputs: candidategenelist_061817.txt unique gene list for all analyses
###          ran_geneslist_061817.txt; 1000 random genes
###          orthologs_rangenes_061817.csv; dataframe of ensembl gene IDs for RANDOM devil, oppossum, wallaby
###          orthologs_cangenes_061817.csv; dataframe of ensembl gene IDs for CANDIDATES devil, oppossum, wallaby

library(biomaRt)
library(dplyr)

#### Create data frame of candidate genes across analyses ####

## for list of  all candidate genes 
AF_topgenes <- read.delim('AF_top.genes.100000bp.txt', header = FALSE)
AF_topgenes  <- data.frame(unique(strtrim(as.character(AF_topgenes[,'V8']), 18)))
colnames(AF_topgenes) <- 'ensembl_gene_id'

spatpg_topgenes<- read.delim('spatg_top.genes.100000bp.txt', header = FALSE)
spatpg_topgenes <- data.frame(unique(strtrim(as.character(spatpg_topgenes[,'V8']), 18)))
colnames(spatpg_topgenes) <- 'ensembl_gene_id'

MM_topgenes<- read.delim('MM_top.genes.100000bp.txt', header = FALSE)
MM_topgenes <- data.frame(unique(strtrim(as.character(MM_topgenes[,'V8']), 18)))
colnames(MM_topgenes) <- 'ensembl_gene_id'

AF_spatpg <- merge(AF_topgenes, spatpg_topgenes, by = 'ensembl_gene_id')
AF_MM <- merge(AF_topgenes, MM_topgenes, by = 'ensembl_gene_id')
spatpg_MM <- merge(spatpg_topgenes, MM_topgenes, by = 'ensembl_gene_id')
AF_spatpg_MM <- merge(AF_MM, spatpg_MM, by = 'ensembl_gene_id') #0 in all

only_spat_noAF <- anti_join(spatpg_topgenes, AF_topgenes)
only_spat_noMM <- anti_join(spatpg_topgenes, MM_topgenes)
only_spatpg <- merge(only_spat_noAF, only_spat_noMM)

only_AF_nospat <- anti_join(AF_topgenes, spatpg_topgenes)
only_AF_noMM <- anti_join(AF_topgenes, MM_topgenes)
only_AF <- merge(only_AF_nospat, only_AF_noMM)

only_MM_nospat <- anti_join(MM_topgenes, spatpg_topgenes)
only_MM_noAF <- anti_join(MM_topgenes, AF_topgenes)
only_MM <- merge(only_MM_nospat, only_MM_noAF)

alluniquecandidatesgenelist <- rbind(only_AF, only_MM, only_spatpg)


## write.table(alluniquecandidatesgenelist, "candidategenelist_061817.txt", sep = "\t", row.names = FALSE, quote = FALSE)

candidategenes <- as.data.frame(alluniquecandidatesgenelist)


##### Create data frame of random genes across analyses ####

## get all avaiable genes
allgenes <-  read.table('~/Code/devils/Inputs/ensGene.txt', header = FALSE)
nrow(unique(allgenes['V13'])) #20419 total unique genes
allg <- substr(x = (allgenes[,'V13']), 0, nchar(as.character(allgenes[1,'V13']))-2)
allg <- as.data.frame(unique(allg))
colnames(allg) <- "ensembl_gene_id"

## remove candidates already being studied here...
nrow(allg)-nrow(candidategenes) == nrow(anti_join(allg, candidategenes)) #must be true
allg <- anti_join(allg, candidategenes)
head(allg)

length(allg[,1])
nrow(allg) - nrow(candidategenes)

## randomly sample, without replacement, 
# ran_genes <- sample_n(allg, 1000, replace = FALSE)
## write.table(ran_genes, "ran_geneslist_061817.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#######

genelist <- candidategenes #allg, candidategenes

## Initialize biomart, devil dataset 

if(interactive()){
  mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = ("sharrisii_gene_ensembl")) 
}
datasets <- listDatasets(mart)[,1]

## Get orthologous gene IDs for tasmanian devil, wallaby, oppossum 
homologs <- strsplit(datasets, split = "_" )
ensembl_name <- c()
for (i in 1:length(homologs)){
  ensembl_name[i] <- paste0(homologs[[i]][1],  "_homolog_ensembl_gene")
}
sort(ensembl_name)
target_species <- ensembl_name[which( ensembl_name=="mdomestica_homolog_ensembl_gene" | ensembl_name == "neugenii_homolog_ensembl_gene") #  ensembl_name == "sharrisii_homolog_ensembl_gene" |
] # devil, oppossum, wallaby
target_species

mdomestica_orthologs <- getBM(attributes = c("ensembl_gene_id",  "mdomestica_homolog_ensembl_gene"), 
                   filters = "ensembl_gene_id", 
                   values = genelist,
                   mart = mart)
mdomestica_deduped_orthologs <- mdomestica_orthologs[-which(duplicated(mdomestica_orthologs[,1])),]

neugenii_orthologs <- getBM(attributes = c("ensembl_gene_id",   "neugenii_homolog_ensembl_gene"), 
                            filters = "ensembl_gene_id", 
                            values = genelist,
                            mart = mart)
neugenii_deduped_orthologs <- neugenii_orthologs[-which(duplicated(neugenii_orthologs[,1])),]


genelist_orthologs <- merge(mdomestica_deduped_orthologs, neugenii_deduped_orthologs)

koala_devil_orthologs <- read.csv("add_koala/KoalaDevilOrthologs_12864_2014_6686_MOESM2_ESM.csv")

koala_devil_orthologs <- koala_devil_orthologs[row.names(unique(koala_devil_orthologs['ensembl_gene_id'])),]


genelist_orthologs_addkoala <- full_join(genelist_orthologs, koala_devil_orthologs[,c("ensembl_gene_id", "koala_transcript_id")], by="ensembl_gene_id")
for (i in 1:ncol(genelist_orthologs_addkoala)){
  genelist_orthologs_addkoala[which(genelist_orthologs_addkoala[,i]==""),i] <- NA
}

genelist_orthologs_addkoala_atleastthree <- c()

for (i in 1:nrow(genelist_orthologs_addkoala)){
  if (sum(is.na(genelist_orthologs_addkoala[i,])) < 2)
    genelist_orthologs_addkoala_atleastthree <- rbind(genelist_orthologs_addkoala_atleastthree, genelist_orthologs_addkoala[i,])
}

head(genelist_orthologs_addkoala_atleastthree)
nrow(genelist_orthologs_addkoala_atleastthree) #10468 whole genome; 244 candidates

write.csv(genelist_orthologs_addkoala_atleastthree, file = "~/Code/devils/candevilsgenelist_orthologs_addkoala_atleastthree_110117.csv", row.names = FALSE) 
write.table(genelist_orthologs_addkoala_atleastthree, file = "~/Code/devils/candevilsgenelist_orthologs_addkoala_atleastthree_110117.tsv", row.names = FALSE, quote = FALSE, sep = "\t") 


devilopossumwallaby <- c("ensembl_gene_id", "mdomestica_homolog_ensembl_gene", "neugenii_homolog_ensembl_gene")

test <- matrix(data = NA, nrow = nrow(genelist_orthologs_addkoala_atleastthree), ncol = 3)

for (i in 1:nrow(genelist_orthologs_addkoala_atleastthree)){
  if (length(colnames(genelist_orthologs_addkoala_atleastthree)[which(!is.na(genelist_orthologs_addkoala_atleastthree[i,])==TRUE)])==4)
    {
    next  
  }
  test[i,] <- colnames(genelist_orthologs_addkoala_atleastthree)[which(!is.na(genelist_orthologs_addkoala_atleastthree[i,])==TRUE)]
}

test <- cbind(genelist_orthologs_addkoala_atleastthree[,1], test)


opp <- test[which(test[,3]=="mdomestica_homolog_ensembl_gene"),]
oppwall <- opp[which(opp[,4]=="neugenii_homolog_ensembl_gene"),]
oppkoala <- opp[which(opp[,4]=="koala_transcript_id"),]

wall <- test[which(test[,3]=="neugenii_homolog_ensembl_gene"),]
oppwall <- opp[which(wall[,4]=="koala_transcript_id"),]




