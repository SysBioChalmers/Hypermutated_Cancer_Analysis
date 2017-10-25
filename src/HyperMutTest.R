library('TCGAbiolinks')
library('readxl')
library('lattice')
library('RColorBrewer')
frequent_mutations100 <- as.matrix(read_excel("data/frequent_mutations100.xlsx")) #Top Mutations found in STAD cancer
xylt2_3cancers        <- as.matrix(read_excel("data/xylt2-3cancers.xlsx"))        #xylt2 patients for 3 cancers
mvk_3cancers          <- as.matrix(read_excel("data/mvk-3cancers.xlsx"))          #mvk patients for 3 cancers
cancer_types <- c('COAD','STAD','UCEC') #The 3 cancers in question

#MAKE BINARY MATRIX
for(i in 1:length(cancer_types)){
  maf            <- GDCquery_Maf(cancer_types[i], pipelines = 'mutect2')
  patient_list   <- unique(substr(maf$Tumor_Sample_Barcode,1,12)) #extract the 12 first character of Tumor Sample Barcode
  patient_matrix <- matrix(0, ncol = length(patient_list), nrow = length(unique(frequent_mutations100[,3]))) #Binary matrix
  rownames(patient_matrix) <- unique(frequent_mutations100[,3]) #row = unique top frequent gene 
  colnames(patient_matrix) <- patient_list #column = all patients of the cancer study
  for(j in 1:length(patient_list)){
    gene_inv <- which(patient_list[j] == substr(maf$Tumor_Sample_Barcode,1,12)) #retrieve index of the genes mutated for patient j 
    patient_matrix[which(unique(frequent_mutations100[,3]) %in% maf$Hugo_Symbol[gene_inv]),j] <- 1 #if patient has top frequent gene, then 1
  }
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
        # annotation_col = annotation_col,
        # annotation_row = annotation_row,
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


