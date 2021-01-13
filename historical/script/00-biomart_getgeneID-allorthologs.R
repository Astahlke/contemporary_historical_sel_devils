#### This generates a dataframe of random and candidate genes to later pull from in 01-getSeq
### inputs: lines 14-22; analysis_top_genes for gene IDs
### outputs: candidategenelist_061817.txt unique gene list for all analyses
###          ran_geneslist_061817.txt; 1000 random genes
###          orthologs_rangenes_061817.csv; dataframe of ensembl gene IDs for RANDOM devil, oppossum, wallaby
###          orthologs_cangenes_061817.csv; dataframe of ensembl gene IDs for CANDIDATES devil, oppossum, wallaby


library(biomaRt)
library(dplyr)
library(stringr)

## Initialize biomart, devil dataset 

if(interactive()){
  mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = ("sharrisii_gene_ensembl")) 
}
datasets <- listDatasets(mart)[,1]
listAttributes(mart)

genelist <- getBM(
  attributes="ensembl_gene_id",
  mart = mart)

## Get orthologous gene IDs for tasmanian devil, wallaby, oppossum 
homologs <- strsplit(datasets, split = "_")
ensembl_name <- c()
for (i in 1:length(homologs)){
  ensembl_name[i] <- paste0(homologs[[i]][1],  "_homolog_ensembl_gene")
}
sort(ensembl_name)
target_species <- ensembl_name[which( ensembl_name=="mdomestica_homolog_ensembl_gene" | ensembl_name == "neugenii_homolog_ensembl_gene") #  ensembl_name == "sharrisii_homolog_ensembl_gene" |
] # devil, oppossum, wallaby
target_species

mdomestica_deduped_orthologs <- mdomestica_orthologs[-which(duplicated(mdomestica_orthologs[,1])),]


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

koala_devil_orthologs <- semi_join(koala_devil_orthologs, genelist)

genelist_orthologs_addkoala <- full_join(genelist_orthologs, koala_devil_orthologs[,c("ensembl_gene_id", "koala_transcript_id")], by="ensembl_gene_id")
for (i in 1:ncol(genelist_orthologs_addkoala)){
  genelist_orthologs_addkoala[which(genelist_orthologs_addkoala[,i]==""),i] <- NA
}

genelist_orthologs_addkoala_atleastthree <- c()

for (i in 1:nrow(genelist_orthologs_addkoala)){
  if (sum(is.na(genelist_orthologs_addkoala[i,])) < 2)
    genelist_orthologs_addkoala_atleastthree <- rbind(genelist_orthologs_addkoala_atleastthree, genelist_orthologs_addkoala[i,])
}


nrow(genelist_orthologs_addkoala_atleastthree)
write.table(genelist_orthologs_addkoala_atleastthree, file = "wholegenome_orthologtable.tsv", row.names = FALSE, quote = FALSE, sep = "\t") 

## rename for simplicity
orth <- genelist_orthologs_addkoala_atleastthree

sets <- matrix(data = NA, nrow = nrow(orth), ncol = 3)

for (i in 1:nrow(orth)){
  if (length(colnames(orth)[which(!is.na(orth[i,])==TRUE)])==4)
  {
    next  
  }
  sets[i,] <- colnames(orth)[which(!is.na(orth[i,])==TRUE)]
}

sets <- cbind(as.character(orth[,1]), sets)

allfour <- sets[which(complete.cases(sets)),]
nrow(allfour)
write.table(paste0(allfour[,1],".phylip"), file = "allfour.txt", row.names = FALSE, quote = FALSE, sep = "\t") 

opp <- sets[which(sets[,3]=="mdomestica_homolog_ensembl_gene"),]
oppwall <- opp[which(opp[,4]=="neugenii_homolog_ensembl_gene"),]
nrow(oppwall) # devil, opp, wall
write.table(paste0(oppwall[,1],".phylip"), file = "oppwall.txt", row.names = FALSE, quote = FALSE, sep = "\t") 

oppkoala <- opp[which(opp[,4]=="koala_transcript_id"),]
nrow(oppkoala) # devil, opp, koala
write.table(paste0(oppkoala[,1],".phylip"), file = "oppkoala.txt", row.names = FALSE, quote = FALSE, sep = "\t") 

wall <- sets[which(sets[,3]=="neugenii_homolog_ensembl_gene"),]
wallkoala <- wall[which(wall[,4]=="koala_transcript_id"),]
nrow(wallkoala) # devil, wallaby, koala
write.table(paste0(wallkoala[,1],".phylip"), file = "wallkoala.txt", row.names = FALSE, quote = FALSE, sep = "\t") 