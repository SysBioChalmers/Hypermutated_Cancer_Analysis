#== RETRIEVE HMR2 LIST OF GENES ==#
setwd("~/Dropbox/PhD/15 - CRISPR X/CancerAnalysis")
source(file = 'CRISPR_Cancer/src/getHMRgenes.R') #Load list of genes in HMR2 GEM
HMR_ECmatrix <- getHMRgenes()
HMR_ECmatrix <- HMR_ECmatrix[,2]

#== MUTATION ==#
library('TCGAbiolinks')

#list of all the different cancer types included in the TCGA studies:
cancer_types <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH',
                  'KIRC','KIRP','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD',
                  'READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')

#== CANCER PURITY ==#
# 1. Generate purity matrix #
load("CRISPR_Cancer/data/PurityTable.rda") #Load the purity scores from tumor of several patients (Aran et al. 2015 Nat. Comm.)

threshold = 300 #!# Filter to study with more than X patients data
sample_size <- sapply(cancer_types, function(x) length(which(x == PurityTable[,2]))) #retrieve size 
PurityMatrix <- matrix(0, nrow=length(which(sample_size > threshold)), ncol = 4) #Purity Matrix which associate a purity score to cancer type
rownames(PurityMatrix) <- rownames(as.matrix(which(sample_size > threshold))) #rownames = Tumor with more than "threshold" number of patient
colnames(PurityMatrix) <- c("CPE","Q1-Q3","IQR","N") #Purity score to be calculated
PurityCPE <- list() #DataFrame of the CPE values for Plotting purposes
for(i in 1:length(PurityMatrix[,1])){ #Loop to fill the purityMatrix 
  type_cancer       <- rownames(PurityMatrix)[i] 
  newMat            <- PurityTable[which(PurityTable[,2] == type_cancer),] 
  PurityMatrix[i,1] <- median(as.numeric(gsub(",", ".", gsub("\\.", "", (newMat[,3])))))                #CPE median
  PurityMatrix[i,2] <- paste0(quantile(as.numeric(gsub(",", ".", gsub("\\.", "", (newMat[,3])))))[[2]], #IQR: Q1-Q3
                                "-",
                                quantile(as.numeric(gsub(",", ".", gsub("\\.", "", (newMat[,3])))))[[4]])
  PurityMatrix[i,3] <- IQR(as.numeric(gsub(",", ".", gsub("\\.", "", (newMat[,3])))))                   #IQR value
  PurityMatrix[i,4] <- length(newMat[,1]) #Number of patients studied 
  PurityCPE[[i]]    <- as.numeric(gsub(",", ".", gsub("\\.", "", (newMat[,3]))))
}
PurityCPE           <- PurityCPE[match(rownames(PurityMatrix[order(PurityMatrix[,1],decreasing = TRUE),]), rownames(PurityMatrix))] #Sort the list 
PurityMatrix        <- PurityMatrix[order(PurityMatrix[,1],decreasing = TRUE),] #Sort matrix with the best CPE first

# 2. Plotting the scores #
library(RColorBrewer)
boxplot(PurityCPE, col="grey58", medcol="black", whiskcol="black", staplecol="black", boxcol="black", outcol="black",  cex=0.5,
        names=rownames(PurityMatrix),
        aspect = "fill",
        main=("(CPE) Consensus measurement of Purity Estimations"))

      
# 3. Mutations Screening ==#
  #3.1 Per Cancer Types   

Purity_threshold = 0.8 #Filtering patients with a CPE < 0.8
ImpactMatrix <- list()
study_tab           <- matrix(0,nrow=length(HMR_ECmatrix), ncol=length(rownames(PurityMatrix)))     
rownames(study_tab) <- HMR_ECmatrix# Matrix with the mutation impact from each HMR gene 
colnames(study_tab) <- rownames(PurityMatrix)
for(i in 1:length(PurityMatrix[,1])){
  Spec_Matrix         <- PurityTable[which(PurityTable[,2] == rownames(PurityMatrix)[i] # Retrieve only patients from the study 
                    & as.numeric(gsub(",", ".", gsub("\\.", "",PurityTable[,3]))) >= Purity_threshold),] # And with a CPE >= threshold  
  if(length(Spec_Matrix[,1]) <= 360){                                                   # Filter study with less than 360 patients having <0.8 CPE
    next
  }
  maf                 <- GDCquery_Maf(rownames(PurityMatrix)[i], pipelines = 'mutect2') # Retrieve mutation annotation file (MAF) from Genomic Data Commons (GDC)
  
  
  for(k in 1:length(HMR_ECmatrix)){
    indexes_gene      <- which(maf$Hugo_Symbol == HMR_ECmatrix[k])
    impacts_gene      <- maf[indexes_gene,]$IMPACT
    if(length(impacts_gene) > 0){
      study_tab[k,i] <- length(which(impacts_gene == "HIGH"))/(length(Spec_Matrix[,1])) #% of the HIGH mutation in the total number of sample
                       #c(length(which(impacts_gene == "LOW")),
                       # length(which(impacts_gene == "MODIFIER")),
                       # length(which(impacts_gene == "MODERATE")),
                       # length(which(impacts_gene == "HIGH")))
    }
  }
  ImpactMatrix[[i]] <- study_tab
}  
# 4. Plot #
library(pheatmap)
col.pal <- RColorBrewer::brewer.pal(3, "Reds")
map <- pheatmap(study_tab[,1:3],
         cluster_cols = T,
         color = col.pal, 
         fontsize = 6.5,
         fontsize_row=6, 
         fontsize_col = 6)#,
pheatmap(study_tab[map$tree_row$order,][1:50,1:4],
        cluster_cols = T,
        color = col.pal, 
        fontsize = 6.5,
        #fontsize_row=6, 
        fontsize_col = 6)


  # 3.2 Per Patients #

Purity_threshold = 0.8 #Filtering patients with a CPE < 0.8
ImpactMatrix <- list()
study_tab           <- matrix(0,nrow=length(HMR_ECmatrix), ncol=length(rownames(PurityMatrix)))     
rownames(study_tab) <- HMR_ECmatrix# Matrix with the mutation impact from each HMR gene 
colnames(study_tab) <- rownames(PurityMatrix)#c("LOW","MODIFIER","MODERATE","HIGH")             # Different types of Impact established
for(i in 1:length(unique(PurityTable[,2]))){
  maf                 <- GDCquery_Maf(rownames(PurityMatrix)[i], pipelines = 'mutect2') # Retrieve mutation annotation file (MAF) from Genomic Data Commons (GDC)
  Spec_Matrix         <- PurityTable[which(PurityTable[,2] == rownames(PurityMatrix)[i] # Retrieve only patients from the study 
                                           & as.numeric(gsub(",", ".", gsub("\\.", "",PurityTable[,3]))) >= Purity_threshold),] # And with a CPE >= threshold 
  
  indexPatient <- unlist(lapply(Spec_Matrix[,1] , function(x) which(x == substr(maf$Tumor_Sample_Barcode, 1, 16)))) # Filter the maf matrix for only patients matching the Spec_Matrix
  newMaf       <- maf[indexPatient,]
  for(k in 1:length(HMR_ECmatrix)){
    indexes_gene      <- which(newMaf$Hugo_Symbol == HMR_ECmatrix[k])
    impacts_gene      <- newMaf[indexes_gene,]$IMPACT
    if(length(impacts_gene) > 0){
      study_tab[k,i] <- length(which(impacts_gene == "HIGH"))/(length(Spec_Matrix[,1]))
    }
  }
  ImpactMatrix[[i]] <- study_tab
}     

## Incidence ##
library(readxl)
#Incidence <- read_excel("~/Dropbox/PhD/15 - CRISPR X/CancerAnalysis/Annotations/Incidence.xlsx")
Incid_Tab  <- study_tab
Incidences <- c(15088,161881,14430,15139,15088,(23635+9824),137109,37419,156045,51320,(9428+18311+974+15514+14200),119739,36433,119739)
Incid_Tab  <- Incid_Tab*Incidences
Incid_Tab  <- Incid_Tab[which(rowSums(Incid_Tab) > 0),] 
Incid_Tab  <- cbind(Incid_Tab,rowSums(Incid_Tab))
Incid_Tab  <- Incid_Tab[do.call(order, as.data.frame(Incid_Tab[,length(Incid_Tab[1,])]))]
Incid_Tab  <- Incid_Tab[do.call(order, as.data.frame(Incid_Tab[,(length(Incid_Tab[1,]))])),]

levelplot(t((Incid_Tab[(length(Incid_Tab[,1])-75):(length(Incid_Tab[,1])),1:14])),
          xlab="Cancer types",
          ylab=paste0("Number of cases"),
          col.regions = heat.colors(100)[length(heat.colors(100)):1],
          aspect="fill",
          main="")
