
biomart  <- 'ENSEMBL_MART_ENSEMBL'
outdir <- 'ensembl_seq/'
dir.create(outdir)

can_orthologs <- "candidate_orthologtable_110117.tsv"
allg_orthologs <- "allg/wholegenome_orthologtable_110717.tsv"

orthologfile <- can_orthologs

genelist_orthologs <- read.table(orthologfile, as.is = TRUE,  na.strings = "", header=TRUE, sep = "\t",  col.names=c("sharrisii", "mdomestica", "neugenii", "koala"))

cat( 'nrow genelistorthologs:', nrow(genelist_orthologs), '\n')

#### This generates a dataframe of random and candidate genes to later pull from in 01-getSeq
### inputs:  ran_geneslist_061817.txt; 1000 random genes
###          orthologs_rangenes_061817.csv; dataframe of ensembl gene IDs for RANDOM devil, oppossum, wallaby
###          orthologs_cangenes_061817.csv; dataframe of ensembl gene IDs for CANDIDATES devil, oppossum, wallaby
### outputs: candidategenelist_orthologs_061817.txt unique gene list for all analyses


library(biomaRt)
library(stringr)


# Read the table of orthologs. I set as.is=TRUE because I don't want
cat('read',  nrow(genelist_orthologs), 'orthologs from', orthologfile, '\n')

# Set up dataset list
datasets = paste0(names(genelist_orthologs),'_gene_ensembl')

cat('eventually will pull genes from:', datasets)

skippedgenes <- c()
aalength <- c()
aalengths <- c()

# Loop over species

k <- 0

for(j in 1:3) {
#j <- 3
    # Set up biomart connection;
    dt = datasets[j]
    cat('currenlty pulling from:', dt, '\n')
    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = dt)

    k = k+1

    # Loop through genes and get sequences
    for(i in 1:nrow(genelist_orthologs)) {
         
      # Get the name of the tasmanian devil gen
        devil_gene_id = as.character(genelist_orthologs[i,1])
        
        # Open file handles for the nucleotide and aa sequences; this
        # will overwrite any existing files
        mode = 'ab'
        if(k == 1)
            mode = 'wb'
        aa_handle = file(paste0(outdir, devil_gene_id, '.aa.fasta'), mode)
        nu_handle = file(paste0(outdir, devil_gene_id, '.nuc.fasta'), mode)
        
        aa = getSequence(seqType = 'peptide',
                         #filters='ensembl_gene_id',
                         type = 'ensembl_gene_id',
                         id = genelist_orthologs[i, j],
                         mart=mart)
        aa_ID <- aa[['ensembl_gene_id']][1]
        aa_sequence <- aa[['peptide']][1]
        np = nchar(aa_sequence)
        aalength <- rbind(aalength, cbind(aa_ID, np))
        if(length(np) == 0 || is.na(np) || np == 0 || aa_sequence == 'Sequence unavailable') {
            skippedgenes <- c(skippedgenes, skippedgenes)
            cat('\t skipped', genelist_orthologs[i, j], '\n')
            close(aa_handle)
            close(nu_handle)
            next
         }
        cat('>', str_sub(aa_ID, 4, 6), str_sub(aa_ID, -5, -1), '\n', aa_sequence, '\n',         # truncate spname so that PAML doesn't get confused abouut seq ID
            sep='', file=aa_handle)
        nu = getSequence(seqType = 'coding',
                                  #filters='ensembl_gene_id',
                                  type = 'ensembl_gene_id',
                                  id = genelist_orthologs[i, j],
                                  mart=mart)
        nu_ID <- nu[['ensembl_gene_id']][1]
        nu_sequence <- nu[['coding']][1]
        cat('>', str_sub(nu_ID, 4, 6), str_sub(nu_ID, -5, -1),'\n', nu_sequence, '\n',
            sep='', file=nu_handle)
        cat('\t finished', aa_ID, 'homolog of', genelist_orthologs[i, 1], ';', (i/nrow(genelist_orthologs))*100, '% finished \n')
        close(aa_handle)
        close(nu_handle)
    }    
    cat('finished species: ', dt, '\n')
    aa_lengths <- merge(genelist_orthologs[,j], aalength)
}

# if j = 2 
skipped_meugenii <- skippedgenes
aalength_meugenii <- aalength
write.csv(skipped_meugenii, "skipped_meugenii")

# 
# for(j in 1:ncol(completecases_allg_orthologs)) {
#   # Set up biomart connection;
#   # to convert strings into factors.
#   devil_gene <- completecases_allg_orthologs[,1]
#   dt = datasets[j]
#   cat(dt, '\n')
#   mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = dt)
#   sp_name = gsub('_gene_ensembl', '', dt)
#   
#   k = k + 1
#   
#   # Loop through genes and get sequences
#   for(i in 4098){ #:length(devil_gene)) {
#     
#     # Get the name of the tasmanian devil gen
#     devil_gene_id = as.character(devil_gene[i])
#     
#     # Open file handles for the nucleotide and aa sequences; this
#     # will overwrite any existing files
#     mode = 'ab'
#     if(k == 1)
#       mode = 'wb'
#     ### for devil-wallaby comparison, added prefix
#     aa_handle = file(paste0('mars_allgenes_090717/', devil_gene_id, '.aa.fasta'), mode)
#     nu_handle = file(paste0('mars_allgenes_090717/', devil_gene_id, '.nuc.fasta'), mode)
#     
#     # Skip if missing
#     
#     if(is.na(allg_orthologs[i, j]))
#       next
#     
#     aa_sequence = getSequence(seqType = 'peptide',
#                               #filters='ensembl_gene_id',
#                               type = 'ensembl_gene_id',
#                               id = allg_orthologs[i, j],
#                               mart=mart)[['peptide']]
#     np = nchar(aa_sequence)
#     if(length(np) == 0 || np == 0 || aa_sequence == 'Sequence unavailable') {
#       close(aa_handle)
#       close(nu_handle)
#       skippedgenes <- c(skippedgenes, allg_orthologs[i, j])
#       cat('skipped', allg_orthologs[i, j], '\n')
#       next
#     }
#     # truncate spname so that PAML doesn't get confused abouut seq ID
#     cat('>', strtrim(sp_name, 8), '\n', aa_sequence[[1]], '\n',
#         sep='', file=aa_handle)
#     nu_sequence = getSequence(seqType = 'coding',
#                               #filters='ensembl_gene_id',
#                               type = 'ensembl_gene_id',
#                               id = allg_orthologs[i, j],
#                               mart=mart)
#     cat('>', strtrim(sp_name, 8), '\n', nu_sequence[[1]], '\n',
#         sep='', file=nu_handle)
#     cat('\tfinished', sp_name, 'homolog of', devil_gene_id, ';', (i/length(devil_gene))*100, '% finished \n')
#     close(aa_handle)
#     close(nu_handle)
#   }    
#   cat('finished species: ', sp_name, '\n')
# }

#which(allg_orthologs=="ENSSHAG00000004148")

write.csv(skippedgenes, "skippedgenes.csv")
