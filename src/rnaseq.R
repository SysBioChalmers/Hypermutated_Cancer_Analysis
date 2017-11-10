library('TCGAbiolinks')
library("ggplot2")
library("SummarizedExperiment")
library("magrittr")
library("ggpubr")
setwd("~/Dropbox/PhD/15 - CRISPR X/CancerAnalysis")
AllPatients <- read.table(file = 'CRISPR_Cancer/data/Allpatients.tsv',sep = '\t', header = TRUE) #Patients with G529 mutations
AllPatients <- AllPatients[which(AllPatients[,2] == "Colorectal"),] #Select only COAD patients

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  #barcode = data$barcode, 
                  legacy = TRUE)

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
data <- GDCprepare(query)

#Fix factor in data
data$disease_type    <- factor(unlist(data$disease_type))                          #unlist added
data$primary_site    <- relevel(factor(unlist(data$primary_site)), ref="Stomach")  #unlist added
data$shortLetterCode <- relevel(factor(unlist(data$shortLetterCode)), ref="NT")    #unlist added
#Filter out samples of unknown primary site
df <- data[,!is.na(data$primary_site)]

RNAmat <- assay(data,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")

#equifax <- list()
#equifax[[2]] <- RNAmat[12642,]

##== GAG RETRIEVAL ==##
GAG.annot  <- getGAGgenes()
GAG.annot  <- GAG.annot[-which(duplicated(GAG.annot[,4]) == TRUE),]
GAG_13     <- paste0(GAG.annot[,3],"|",GAG.annot[,1]) #combuine column 3 and 1 in order to match the RNAseq matrix, XYLT2|64132 
IndexOfGAG <- sapply(GAG_13, function(x) which(rownames(RNAmat) == x))#find index of XYLT2
GAG_13 <- cbind(GAG_13,0) #Dirty for loop to retrieve the index of the GAG genes in the matrix
for(i in 1:length(GAG_13[,1])){
  GAG_13[i,2] <- which(GAG_13[i,1] == rownames(RNAmat))[1]
}
GAG_13 <- GAG_13[-(which(is.na(GAG_13[,2]) == TRUE)),] #Remove genes with no RNAseq data 

NewSTADmat <- as.matrix(RNAmat[as.numeric(GAG_13[,2]),]) #Extract RNA of the GAG genes

NewSTADmat_TP <- NewSTADmat[,-which(substr(data$shortLetterCode, 1, 12) == "NT")] #create a matrix without NT samples

for (i in 1:length(AllPatients[,1])){ # remove patients carrying G529 mutation from the cancer study (poorly coded)
  indexofpatient <- which(rownames(AllPatients)[i] == substr(colnames(NewSTADmat),1,12))
  if (length(indexofpatient) == 0){
    next
  }
  NewSTADmat_TP  <- NewSTADmat[,-indexofpatient] 
}

corrr <- cor.test(NewSTADmat_TP[1,], NewSTADmat_TP[2,], 
         method = "pearson")

#Only G529 Patient
NewSTADmat_G529 <- NewSTADmat[,-which(substr(data$shortLetterCode, 1, 12) == "NT")] #create a matrix without NT samples
indxG529 <- integer(0)
for (i in 1:length(AllPatients[,1])){ 
  indexofpatient <- which(rownames(AllPatients)[i] == substr(colnames(NewSTADmat),1,12))
  if (length(indexofpatient) == 0){
    next
  }
  indxG529 <- c(indxG529,indexofpatient)
  #NewSTADmat_TP  <- NewSTADmat_TP[,-indexofpatient] 
}
NewSTADmat_G529 <- NewSTADmat[,indxG529]

library(pheatmap)
col.pal <- RColorBrewer::brewer.pal(9, "Reds")
pairheatmap(log1p(NewSTADmat_G529), log1p(NewSTADmat_TP[,1:100]))
map <- pheatmap((NewSTADmat_TP[,1:20]+0.001),
                cluster_cols = T,
                color = col.pal, 
                fontsize = 6.5,
                fontsize_row=6, 
                fontsize_col = 6)#,


boxplot(NewSTADmat[1,])
corrr2 <- cor.test(NewSTADmat_G529[1:10,], NewSTADmat_TP[1:10,], 
                  method = "pearson")




boxplot(equifax,
        xlab="Cancer types")
equifax2 <- list();NewRNAmat <- ""

## RNA seq of NT vs TP vs G529 ##

equifax <- list()
equifax[[1]] <- RNAmat[12642,which(data$shortLetterCode == "NT")]
equifax[[2]] <- RNAmat[12642,which(data$shortLetterCode == "TP")]
indexofpatient <- ""
for (i in 1:length(AllPatients[,1])){ # remove patients carrying G529 mutation from the cancer study (poorly coded)
  indexofpatient <- c(indexofpatient,which(rownames(AllPatients)[i] == substr(colnames(RNAmat),1,12)))
  if (length(indexofpatient) == 0){
    next
  }
  #NewRNAmat  <- RNAmat[,indexofpatient] 
}
indexofpatient <- indexofpatient[-1]
equifax[[3]] <- RNAmat[12642,-as.numeric(indexofpatient)]
equifax[[4]] <- RNAmat[12642,as.numeric(indexofpatient)]
boxplot(equifax,
        col="grey")


