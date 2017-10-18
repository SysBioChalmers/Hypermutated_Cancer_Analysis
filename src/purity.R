#== RETRIEVE HMR2 LIST OF GENES ==#
source(file = 'getHMRgenes.R')
HMR_ECmatrix <- getHMRgenes()
HMR_ECmatrix <- HMR_ECmatrix[,2]

#== MUTATION ==#
library('TCGAbiolinks')

#list of all the different cancer types included in the TCGA studies:
cancer_types <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH',
                  'KIRC','KIRP','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD',
                  'READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')

#== Cancer Purity ==#
library(ggplot2)
load("~/Dropbox/PhD/15 - CRISPR X/CancerAnalysis/Annotations/PurityTable.rda")

PurityMatrix <- matrix(0, nrow=length(unique(PurityTable[,2])), ncol = 4)
rownames(PurityMatrix) <- unique(unique(PurityTable[,2]))
colnames(PurityMatrix) <- c("CPE","Q1-Q3","IQR","N")
idxtoremove <- ""

for(i in 1:length(PurityMatrix[,1])){
  type_cancer       <- rownames(PurityMatrix)[i]
  if (length(which(PurityTable[,2] == type_cancer)) < 400){ #Filter to study with more than 400 patients data
    #PurityTable     <- PurityTable[-which(PurityTable[,2] == type_cancer),]
    #PurityMatrix    <- PurityMatrix[-which(rownames(PurityMatrix) == type_cancer),]
    idxtoremove      <- c(idxtoremove,i)
    next
  }
  else{
    newMat            <- PurityTable[which(PurityTable[,2] == type_cancer),]
    PurityMatrix[i,1] <- median(as.numeric(gsub(",", ".", gsub("\\.", "", (newMat[,3])))))                #CPE median
    PurityMatrix[i,2] <- paste0(quantile(as.numeric(gsub(",", ".", gsub("\\.", "", (newMat[,3])))))[[2]], #IQR: Q1-Q3
                                "-",
                                quantile(as.numeric(gsub(",", ".", gsub("\\.", "", (newMat[,3])))))[[4]])
    PurityMatrix[i,3] <- IQR(as.numeric(gsub(",", ".", gsub("\\.", "", (newMat[,3])))))                   #IQR value
    PurityMatrix[i,4] <- length(newMat[,1]) #Number of patients studied 
  }
}
idxtoremove        <-  idxtoremove[-1]
PurityMatrix       <-  PurityMatrix[-as.numeric(idxtoremove),]
PurityMatrix       <-  PurityMatrix[order(PurityMatrix[,1],decreasing = TRUE),]
PurityCPE <- list()
for(i in 1:length(PurityMatrix[,1])){
  type_cancer       <- rownames(PurityMatrix)[i]
  newMat            <- PurityTable[which(PurityTable[,2] == type_cancer),]
  PurityCPE[[i]]    <- as.numeric(gsub(",", ".", gsub("\\.", "", (newMat[,3])))) # Generate a Dataframe of the CPE values for boxplot
}
c1 <- "black" #color
c2 <- rainbow(40, alpha=0.2) #color
c3 <- rainbow(40, v=0.1) #color
boxplot(PurityCPE, col=c2, medcol=c3, whiskcol=c1, staplecol=c3, boxcol=c3, outcol=c3,  cex=0.5,
        names=rownames(PurityMatrix),
        main="(CPE) Consensus measurement of Purity Estimations ")  

#save(PurityMatrix, PurityCPE, file = "PurityFiltering.rda")


#== Mutations Screening ==#

#load(file = 'PurityFiltering.rda')
Purity_threshold = 0.8 #Filtering patients with a CPE < 0.8
ImpactMatrix <- list()
study_tab           <- matrix(0,nrow=length(HMR_ECmatrix), ncol=length(rownames(PurityMatrix)))     
rownames(study_tab) <- HMR_ECmatrix# Matrix with the mutation impact from each HMR gene 
colnames(study_tab) <- rownames(PurityMatrix)#c("LOW","MODIFIER","MODERATE","HIGH")             # Different types of Impact established
for(i in 1:length(PurityMatrix[,1])){
  maf                 <- GDCquery_Maf(rownames(PurityMatrix)[i], pipelines = 'mutect2') # Retrieve mutation annotation file (MAF) from Genomic Data Commons (GDC)
  Spec_Matrix         <- PurityTable[which(PurityTable[,2] == rownames(PurityMatrix)[i] # Retrieve only patients from the study 
                    & as.numeric(gsub(",", ".", gsub("\\.", "",PurityTable[,3]))) >= Purity_threshold),] # And with a CPE >= threshold  
  if(length(Spec_Matrix[,1]) <= 400){                                                    # Filter study with less than 400 patients having <0.8 CPE
    print("fuck")
    next
  }
#  study_tab           <- matrix(0,nrow=length(HMR_ECmatrix), ncol=4)                    # Matrix with the mutation impact from each HMR gene 
#  rownames(study_tab) <- HMR_ECmatrix
#  colnames(study_tab) <- namesOcancer#c("LOW","MODIFIER","MODERATE","HIGH")             # Different types of Impact established
  indexPatient        <- 0
  for (j in 1:length(Spec_Matrix[,1])){
    indexPatient <- c(indexPatient,which(substr(maf$Tumor_Sample_Barcode, 1, 16) == Spec_Matrix[j,1]))  # Filter the maf matrix for only patients matching the Spec_Matrix (slow way)
  }
  indexPatient <- indexPatient[-1]
  
  for(k in 1:length(HMR_ECmatrix)){
    indexes_gene      <- which(maf$Hugo_Symbol == HMR_ECmatrix[k])
    impacts_gene      <- maf[indexes_gene,]$IMPACT
    if(length(impacts_gene) > 0){
      study_tab[k,i] <- length(which(impacts_gene == "HIGH"))/(length(Spec_Matrix[,1]))
                       #c(length(which(impacts_gene == "LOW")),
                       # length(which(impacts_gene == "MODIFIER")),
                       # length(which(impacts_gene == "MODERATE")),
                       # length(which(impacts_gene == "HIGH")))
      
    }
  }
  #name_file <- paste0(rownames(PurityMatrix)[i],".csv")
  #write.csv(study_tab,file=name_file)
  ImpactMatrix[[i]] <- study_tab
}     
#save(ImpactMatrix, file = "ImpacMat_400_CPE8.rda")
ans <- study_tab*100
idx <- ans[which(rowSums(ans) > 0),] 
ans <- ans[idx,]
ans <- cbind(ans,ans[,1]+ans[,2]+ans[,3]+ans[,4])
ans <- ans[do.call(order, as.data.frame(ans[,5])),]
levelplot(t((ans[(length(ans[,1])-35):length(ans[,1])-1,1:4])),
          xlab="Cancer types",
          ylab=paste0("% High Impact"),
          col.regions = heat.colors(100)[length(heat.colors(100)):1],
          main="")



#######################################
#  All Patients with a certain CPE    #
#######################################

Purity_threshold = 0.8 #Filtering patients with a CPE < 0.8
ImpactMatrix <- list()
study_tab           <- matrix(0,nrow=length(HMR_ECmatrix), ncol=length(rownames(PurityMatrix)))     
rownames(study_tab) <- HMR_ECmatrix# Matrix with the mutation impact from each HMR gene 
colnames(study_tab) <- rownames(PurityMatrix)#c("LOW","MODIFIER","MODERATE","HIGH")             # Different types of Impact established
for(i in 1:length(unique(PurityTable[,2]))){
  maf                 <- GDCquery_Maf(rownames(PurityMatrix)[i], pipelines = 'mutect2') # Retrieve mutation annotation file (MAF) from Genomic Data Commons (GDC)
  Spec_Matrix         <- PurityTable[which(PurityTable[,2] == rownames(PurityMatrix)[i] # Retrieve only patients from the study 
                                           & as.numeric(gsub(",", ".", gsub("\\.", "",PurityTable[,3]))) >= Purity_threshold),] # And with a CPE >= threshold  
  
  #  study_tab           <- matrix(0,nrow=length(HMR_ECmatrix), ncol=4)                    # Matrix with the mutation impact from each HMR gene 
  #  rownames(study_tab) <- HMR_ECmatrix
  #  colnames(study_tab) <- namesOcancer#c("LOW","MODIFIER","MODERATE","HIGH")             # Different types of Impact established
  indexPatient        <- 0
  for (j in 1:length(Spec_Matrix[,1])){
    indexPatient <- c(indexPatient,which(substr(maf$Tumor_Sample_Barcode, 1, 16) == Spec_Matrix[j,1]))  # Filter the maf matrix for only patients matching the Spec_Matrix (slow way)
  }
  indexPatient <- indexPatient[-1]
  newMaf       <- maf[indexPatient,]
  for(k in 1:length(HMR_ECmatrix)){
    indexes_gene      <- which(newMaf$Hugo_Symbol == HMR_ECmatrix[k])
    impacts_gene      <- newMaf[indexes_gene,]$IMPACT
    if(length(impacts_gene) > 0){
      study_tab[k,i] <- length(which(impacts_gene == "HIGH"))/(length(Spec_Matrix[,1]))
      #c(length(which(impacts_gene == "LOW")),
      # length(which(impacts_gene == "MODIFIER")),
      # length(which(impacts_gene == "MODERATE")),
      # length(which(impacts_gene == "HIGH")))
    }
  }
  #name_file <- paste0(rownames(PurityMatrix)[i],".csv")
  #write.csv(study_tab,file=name_file)
  ImpactMatrix[[i]] <- study_tab
}     
#save(ImpactMatrix, file = "ImpacMat_400_CPE8.rda")

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
