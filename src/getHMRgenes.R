getHMRgenes <- function(addnovel=T){
  load("~/Dropbox/PhD/15 - CRISPR X/CancerAnalysis/CRISPR_Cancer/data/HMRlistgeneCurated.rda")
  #ENS_list     <-  as.matrix(read.csv('HMR_ENSG.csv'))
  #ENS_ListCor  <-  as.matrix(read.csv('HMR_geneNames.csv'))
  #ENS_Matrix   <-  cbind(ENS_list,ENS_ListCor)
  
  #for(i in 1:length(indexUnknown)){
  #  foo          <- c(HMR_ECmatrix[indexUnknown[i]],
  #                    ENS_Matrix[which(ENS_Matrix[,1] == HMR_ECmatrix[indexUnknown[i]]),2])
  #  G_list2 <- rbind(G_list2,foo)
  #}
  
  HMR_list <- G_list2
  return(HMR_list)
}