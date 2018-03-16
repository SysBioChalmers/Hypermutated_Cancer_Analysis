#Convert rownames of a lrt.table from ENSEMBL to gene ID

getIDgenes <- function(input){
require(biomaRt)
  #Convert ENSEMBL to gene ID
  ensembl <- useMart('ENSEMBL_MART_ENSEMBL')                      # select the Ensembl BioMart database
  ensembl <- useDataset('hsapiens_gene_ensembl', mart = ensembl)  # we want the Homo sapiens genes
  bm <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),   # these are the data we want to extract
              filters = 'ensembl_gene_id',                        # this is the type of data that we currently have
              values = rownames(input),                           # these are the genes for which symbols are desired
              mart = ensembl)
  
  gene_symbols <- bm[match(rownames(input),bm[,1]),2]
  keep <- !(gene_symbols %in% c('',NA)) & !duplicated(gene_symbols)
  # filter out removed genes from the data
  input <- input[keep,]
  rownames(input) <- gene_symbols[keep]
  return(input)
}
