library(biomaRt)
library(stringr)

projdir <- "historical/"
ortholgIDdir <- "00-getorthologID/"

outdir <- paste0(projdir, "01-nuc_aa.fasta/ensembl_seq/")
dir.create(paste0(projdir, outdir), recursive = T, showWarnings = F)

orthologfile <- paste0(projdir, ortholgIDdir,"wholegenome_orthologtable.tsv")

genelist_orthologs <- read.table(orthologfile, as.is = TRUE,  na.strings = "", header=TRUE, sep = "\t",  col.names=c("sharrisii", "mdomestica", "neugenii", "koala"))
cat( 'nrow genelistorthologs:', nrow(genelist_orthologs), '\n')

# Read the table of orthologs. I set as.is=TRUE because I don't want
cat('read',  nrow(genelist_orthologs), 'orthologs from', orthologfile, '\n')

# Set up dataset list
biomart  <- 'ENSEMBL_MART_ENSEMBL'
datasets = paste0(names(genelist_orthologs),'_gene_ensembl')

cat('eventually will pull genes from:', datasets)


# Loop over species
skippedgenes <- c()
aalength <- c()
aalengths <- c()

k <- 0

for(j in 1:3) {
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