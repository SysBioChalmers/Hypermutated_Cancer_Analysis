setwd("~/Dropbox/PhD/15 - CRISPR X/CancerAnalysis")
library('TCGAbiolinks')
library('readxl')
library('lattice')
library('RColorBrewer')
#frequent_mutations100 <- as.matrix(read_excel("CRISPR_Cancer/data/frequent_mutations100.xlsx")) #Top Mutations found in STAD cancer
xylt2_3cancers        <- as.matrix(read_excel("CRISPR_Cancer/data/xylt2-3cancers.xlsx"))        #xylt2 patients for 3 cancers
mvk_3cancers          <- as.matrix(read_excel("CRISPR_Cancer/data/mvk-3cancers.xlsx"))          #mvk patients for 3 cancers
cancer_types <- c('COAD','STAD','UCEC') #The 3 cancers in question

#MAKE BINARY MATRIX
ListCancer <- list()
for(i in 1:length(cancer_types)){
  maf            <- GDCquery_Maf(cancer_types[i], pipelines = 'mutect2',)
  patient_list   <- unique(substr(maf$Tumor_Sample_Barcode,1,12)) #extract the 12 first character of Tumor Sample Barcode
  patient_matrix <- matrix(0, ncol = length(patient_list), nrow = length(unique(frequent_mutations100[,3]))) #Binary matrix
  rownames(patient_matrix) <- unique(frequent_mutations100[,3]) #row = unique top frequent gene 
  colnames(patient_matrix) <- patient_list #column = all patients of the cancer study
  for(j in 1:length(patient_list)){
    gene_inv <- which(patient_list[j] == substr(maf$Tumor_Sample_Barcode,1,12)) #retrieve index of the genes mutated for patient j 
    patient_matrix[which(unique(frequent_mutations100[,3]) %in% maf$Hugo_Symbol[gene_inv]),j] <- 1 #if patient has top frequent gene, then 1
  }
  ListCancer[[i]] <- patient_matrix
}

#LEVEL PLOT#
#XYLT2#
xylt2Ind <- which(colnames(patient_matrix) %in% xylt2_3cancers[,1]) #from binary matrix which ones are G529 patients
xylt2matrix <- patient_matrix
xylt2matrix[,xylt2Ind] <- xylt2matrix[,xylt2Ind]*3 #ugly way to make those patients appear in the binary matrix without changing the info
levelplot(t(xylt2matrix),
          aspect = "fill",
          xlab = "patients",
          ylab="genes",
          main=paste0(cancer_types[i]," (Red=Xylt2 patients)"),
          col.regions = heat.colors(100)[length(heat.colors(100)):1])

levelplot(t(patient_matrix[,xylt2Ind]), #Only XYLT2 patients
          aspect = "fill",
          xlab = "patients",
          ylab="genes",
          main=paste0(cancer_types[i]," (Only Xylt2 patients)"),
          col.regions = heat.colors(100)[length(heat.colors(100)):1])
#MVK#
mvkInd <- which(colnames(patient_matrix) %in% mvk_3cancers[,1])
mvkmatrix <- patient_matrix
mvkmatrix[,mvkInd] <- mvkmatrix[,mvkInd]*3
levelplot(t(mvkmatrix),
          aspect = "fill",
          xlab = "patients",
          ylab="genes",
          main=paste0(cancer_types[i]," (Red=MVK patients)"),
          col.regions = heat.colors(100)[length(heat.colors(100)):1])

levelplot(t(patient_matrix[,mvkInd]),
          aspect = "fill",
          xlab = "patients",
          ylab="genes",
          main=paste0(cancer_types[i]," (Only MVK patients)"),
          col.regions = heat.colors(100)[length(heat.colors(100)):1])
#HEATMAP#
library(pheatmap)
pheatmap(patient_matrix)

#MafTools#
library(maftools)
maf <- read.maf('TCGA.UCEC.mutect.df91f748-56d5-41fb-927d-f1d280abe8db.DR-7.0.somatic.maf.gz') #LOAD MAF.

maf2 <- subsetMaf(maf)
mutated_gene <- getGeneSummary(maf)
genelimit      <- as.matrix(mutated_gene[1:131,1])#,(which(mutated_gene[,1] == "XYLT2"))),1])
patient_list   <- unique(substr(maf2$Matched_Norm_Sample_Barcode,1,12)) #extract the 12 first character of Tumor Sample Barcode
patient_matrix <- matrix(0, ncol = length(patient_list), nrow = length(as.matrix(genelimit))) #Binary matrix
rownames(patient_matrix) <- as.matrix(genelimit) #row = unique top frequent gene 
colnames(patient_matrix) <- patient_list #column = all patients of the cancer study
for(j in 1:length(patient_list)){
  gene_inv <- which(patient_list[j] == substr(maf2$Matched_Norm_Sample_Barcode,1,12)) #retrieve index of the genes mutated for patient j 
  patient_matrix[which(as.matrix(unique(as.matrix(genelimit))) %in% maf2$Hugo_Symbol[gene_inv]),j] <- 1 #if patient has top frequent gene, then 1
}

#XYLT2#
patient_matrix2 <- patient_matrix
whereisXylt <- which(colnames(patient_matrix2) %in% xylt2_3cancers[,1] == TRUE)
patient_matrix2[,whereisXylt] <- patient_matrix2[,whereisXylt]*2
#patient_matrix2[,which(patient_matrix2[whereisXylt,] == 2)] <- patient_matrix2[,which(patient_matrix2[whereisXylt,] == 2)]*2
col.pal <- RColorBrewer::brewer.pal(9, "Blues")
pheatmap(patient_matrix2,
         cluster_row = T,
         #cluster_cols = F,
         #annotation_col = annotation_col,
         #annotation_row = annotation_row,
         color = col.pal, 
         fontsize = 6.5,
         fontsize_row=6, 
         fontsize_col = 6)#,
         #gaps_col=50)

#MVK#
patient_matrix2 <- patient_matrix
whereisMVK <- which(colnames(patient_matrix2) %in% mvk_3cancers[,1] == TRUE)
patient_matrix2[,whereisMVK] <- patient_matrix2[,whereisMVK]*2
#patient_matrix2[,which(patient_matrix2[whereisMVK,] == 2)] <- patient_matrix2[,which(patient_matrix2[whereisMVK,] == 2)]*2
col.pal <- RColorBrewer::brewer.pal(9, "Blues")
pheatmap(patient_matrix2,
         cluster_row = T,
         #cluster_cols = F,
         # annotation_col = annotation_col,
         # annotation_row = annotation_row,
         color = col.pal, 
         fontsize = 6.5,
         fontsize_row=6, 
         fontsize_col = 6)#,
#gaps_col=50)

#COMBINED#
#patient_matrix2 <- patient_matrix
whereisDouble <- which(whereisXylt %in% whereisMVK== TRUE)
whereisXplusM <- whereisXylt[whereisDouble]
patient_matrix2[,whereisXplusM] <- patient_matrix2[,whereisXplusM]*2
#patient_matrix2[,which(patient_matrix2[whereisXylt,] == 2)] <- patient_matrix2[,which(patient_matrix2[whereisXylt,] == 2)]*2
col.pal <- RColorBrewer::brewer.pal(9, "Blues")
pheatmap(patient_matrix2,
         cluster_row = T,
         #cluster_cols = F,
         # annotation_col = annotation_col,
         # annotation_row = annotation_row,
         color = col.pal, 
         fontsize = 6.5,
         fontsize_row=6, 
         fontsize_col = 6)#,
#gaps_col=50)



plotmafSummary(maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#We will draw oncoplots for top ten mutated genes. (Removing non-mutated samples from the plot for better visualization)
oncoplot(maf,genes = as.matrix(mutated_gene[c(1:30,121),1]), removeNonMutated = TRUE)

mafmutsig <- prepareMutSig(maf)
laml.pancan = pancanComparision(mutsigResults = mafmutsig, qval = 0.1,
                               cohortName = 'COAD', inputSampleSize = 200,
                               label = 1, normSampleSize = TRUE)




#Histogram Distribution of Xylt2#
matrixHist <- matrix(0,nrow = 10,ncol=3)
colnames(matrixHist) <- cancer_types
rownames(matrixHist) <- c("0","1-5","6-10","11-15","16-20","21-25","26-30","31-35","36-40",">41")

for(j in 1:length(ListCancer)){
  #whereisXylt <- which(colnames(ListCancer[[j]]) %in% xylt2_3cancers[,1] == TRUE)
  whereisTP53 <- which(ListCancer[[j]][which(rownames(ListCancer[[j]]) == "BRAF"),] == 1)
  #patient_matrix2 <- ListCancer[[j]][,whereisXylt]
  patient_matrix2 <- ListCancer[[j]][,whereisTP53]
  for(k in 1:length(patient_matrix2[1,])){
    if(sum(patient_matrix2[,k]) == 0 ){
      matrixHist[1,j] <- matrixHist[1,j]+1
    }
    if(sum(patient_matrix2[,k]) > 0 && sum(patient_matrix2[,k]) < 6){
      matrixHist[2,j] <- matrixHist[2,j]+1
    }
    if(sum(patient_matrix2[,k]) > 5 && sum(patient_matrix2[,k]) < 11){
      matrixHist[3,j] <- matrixHist[3,j]+1
    }
    if(sum(patient_matrix2[,k]) > 10 && sum(patient_matrix2[,k]) < 16){
      matrixHist[4,j] <- matrixHist[4,j]+1
    }
    if(sum(patient_matrix2[,k]) > 15 && sum(patient_matrix2[,k]) < 21){
      matrixHist[5,j] <- matrixHist[5,j]+1
    }
    if(sum(patient_matrix2[,k]) > 20 && sum(patient_matrix2[,k]) < 26){
      matrixHist[6,j] <- matrixHist[6,j]+1
    }
    if(sum(patient_matrix2[,k]) > 25 && sum(patient_matrix2[,k]) < 31){
      matrixHist[7,j] <- matrixHist[7,j]+1
    }
    if(sum(patient_matrix2[,k]) > 30 && sum(patient_matrix2[,k]) < 36){
      matrixHist[8,j] <- matrixHist[8,j]+1
    }
    if(sum(patient_matrix2[,k]) > 35 && sum(patient_matrix2[,k]) < 41){
      matrixHist[9,j] <- matrixHist[9,j]+1
    }
    if(sum(patient_matrix2[,k]) > 41){
      matrixHist[10,j] <- matrixHist[10,j]+1
    }
  }
}

mtrxforannot <- matrix("X=0",nrow=length(CombinedCancer[1,]),ncol=2) #matrix for annotation

CombinedCancer <- cbind(ListCancer[[1]],ListCancer[[2]],ListCancer[[3]]) #combine all 3 cancer
mtrxforannot[which(colnames(CombinedCancer) %in% mvk_3cancers[,1] == TRUE),1] <- "X=1"  #which patients are the G529 patient
mtrxforannot[1:length(ListCancer[[1]][1,]),2] <- "COAD" #which patients are part of COAD cancer
mtrxforannot[length(ListCancer[[1]][1,])+1:(length(ListCancer[[1]][1,])+1)+(length(ListCancer[[2]][1,])),2] <- "STAD"
mtrxforannot[837:1366,2] <- "UCEC"
rownames(mtrxforannot) <- colnames(CombinedCancer) 
annotation_col <- data.frame(
  Patient = mtrxforannot[,1],
  Cancer = mtrxforannot[,2])

mtrxforannot2 <- matrix(0, nrow=length(CombinedCancer[,1]))
mtrxforannot2[which(rownames(CombinedCancer) == "MVK"),] = "MVK"
annotation_row <- data.frame(
  GOI = mtrxforannot2)
rownames(annotation_row) <- rownames(CombinedCancer)
col.pal <- brewer.pal(9, "Blues")
pheatmap(CombinedCancer,
         cluster_row = T,
         cluster_cols = T,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         #annotation_colors = rainbow(1:1000),
         show_colnames = F,
         color = col.pal, 
         #fontsize = 6.5,
         fontsize_row=5, 
         fontsize_col = 1,gaps_col=3)

col_groups <- substr(colnames(mat), 1, 1)


annotation_col <- data.frame()




pheatmap(test,
         #cluster_row = T,
         #cluster_cols = F,
         #annotation_row = xylt2_3cancers[,1],
         # annotation_row = annotation_row,
         color = col.pal, 
         fontsize = 6.5,
         fontsize_row=6, 
         fontsize_col = 6)#,
