library('TCGAbiolinks')
library("ggplot2")
library("SummarizedExperiment")
library("magrittr")
library("ggpubr")

AllPatients <- read.table(file = '../data/Allpatients.tsv',sep = '\t', header = TRUE) #Patients with G529 mutations
AllPatients <- AllPatients[which(AllPatients[,2] == "Stomach"),] #Select only COAD patients

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-STAD", 
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
IndexofXylt1 <- which(rownames(RNAmat) == "XYLT2|64132") #find index of XYLT2
IndexofXylt2 <- which(rownames(RNAmat) == "XYLT1|64131") #find index of XYLT1
NewSTADmat <- as.matrix(RNAmat[IndexofXylt1:IndexofXylt2,]) #Retrieve XYLT1 and XYLT2 RNA seq data 

NewSTADmat_TP <- NewSTADmat[,-which(substr(data$shortLetterCode, 1, 12) == "NT")] #create a matrix without NT samples

for (i in 1:length(AllPatients[,1])){ # remove patients carrying G529 mutation from the cancer study (poorly coded)
  indexofpatient <- which(rownames(AllPatients)[i] == substr(colnames(NewSTADmat_TP),1,12))
  if (length(indexofpatient) == 0){
    next
  }
  NewSTADmat_TP  <- NewSTADmat_TP[,-indexofpatient] 
}

corrr <- cor.test(NewSTADmat_TP[1,], NewSTADmat_TP[2,], 
         method = "pearson")

#Only G529 Patient
NewSTADmat_G529 <- NewSTADmat[,-which(substr(data$shortLetterCode, 1, 12) == "NT")] #create a matrix without NT samples
indxG529 <- integer(0)
for (i in 1:length(AllPatients[,1])){ 
  indexofpatient <- which(rownames(AllPatients)[i] == substr(colnames(NewSTADmat_TP),1,12))
  if (length(indexofpatient) == 0){
    next
  }
  indxG529 <- c(indxG529,indexofpatient)
  #NewSTADmat_TP  <- NewSTADmat_TP[,-indexofpatient] 
}
NewSTADmat_G529 <- NewSTADmat[,indxG529]
corrr2 <- cor.test(NewSTADmat_G529[1,], NewSTADmat_G529[2,], 
                  method = "pearson")
