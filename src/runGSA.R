#Perform a gene set analysis (GSA) starting from the lrt.table generated in multiparameters.R

#1. Input preparation
DE_res <- getIDgenes(lrt.table)
qvals  <- dplyr::select(DE_res,FDR.hmRegular)
logFC  <- dplyr::select(DE_res,logFC)
gscGO  <- loadGSC('CRISPR_Cancer/data/c2.cp.kegg.v6.0.symbols.gmt')

#2. Run GSA
gsaRes1 <- runGSA(qvals,geneSetStat="mean",gsc=gscGO,
                  nPerm=10000,gsSizeLim=c(10,800))
gsaRes2 <- runGSA(qvals,geneSetStat="median",gsc=gscGO,
                  nPerm=10000,gsSizeLim=c(10,800))
gsaRes3 <- runGSA(qvals,geneSetStat="sum",gsc=gscGO,
                  nPerm=10000,gsSizeLim=c(10,800))
#gsaRes4 <- runGSA(geneLevelStats = logFC, directions = logFC, geneSetStat = 'fgsea', #try FGSEA method
#                  gsc = gscGO, nPerm = 10000, gsSizeLim = c(10,800))
gsaRes5 <- runGSA(qvals,logFC,geneSetStat="fisher",gsc=gscGO,
                  nPerm=10000,gsSizeLim=c(10,800))
gsaRes6 <- runGSA(qvals,logFC,geneSetStat="stouffer",gsc=gscGO,
                  nPerm=10000,gsSizeLim=c(10,800))
gsaRes7 <- runGSA(qvals,logFC,geneSetStat="tailStrength",gsc=gscGO,
                  nPerm=10000,gsSizeLim=c(10,800))
gsaRes <- runGSA(geneLevelStats = qvals, directions = logFC, geneSetStat = 'stouffer', #FGSEA method is faster
                 gsc = gscGO, nPerm = 10000, gsSizeLim = c(10,200))

resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes5,gsaRes6,gsaRes7)
names(resList) <- c("mean","median","sum","fisher",
                    "stouffer","tailStrength")
ch <- consensusHeatmap(resList,cutoff=20,method="mean")
View(ch$pMat)

# (optional) export the table of results to view or analyze externally
GSAsummaryTable(gsaRes, save = TRUE, file = paste0("",'gsa_results.xls'))
writeFilesForKiwi(gsaRes, label = 'outputFileName', overwrite = FALSE)
