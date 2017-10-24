library('TCGAbiolinks')
library('readxl')
library('lattice')

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



