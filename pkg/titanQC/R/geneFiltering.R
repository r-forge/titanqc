geneFiltering <- function (eset, varDMSO = TRUE, minFC = NULL, INI = FALSE, minSampleVar = NULL) {
  
  if (varDMSO == TRUE) {
    varC <- apply(exprs(eset)[, !eset$GroupID %in% c("DMSO", "DMSO2")], 1, var)
    varD <- apply(exprs(eset)[, !eset$GroupID %in% c("DMSO")], 1, var)
    varD2 <- apply(exprs(eset)[, eset$GroupID %in% c("DMSO", "DMSO2")], 1, var)
    dropListVar <- names(varC)[varC <= varD | varC <= varD2]
  } else {
    dropListVar <- c()
  }
  
  if (!is.null(minFC)) {
    maxFoldChange <- apply(abs(exprs(eset)[, !eset$treat %in% c("DMSO", "DMSO2")]), 1, max)
    droplistFC <- names(maxFoldChange)[maxFoldChange < minFC]
  } else {
    droplistFC <- c()
  }
  
  
  if (INI == TRUE) {
    droplistFarmsLaplace <- rownames(fData(eset)[fData(eset)[, "INIcalls"] == "NI", ])
  } else {
    droplistFarmsLaplace <- c()
  }
  
  
  if (!is.null(minSampleVar)) {
    varSample <- apply(exprs(eset)[, !eset$GroupID %in% c("DMSO", "DMSO2")], 1, var)
    droplistVarSample <- names(varSample)[varSample < minSampleVar]
  } else {
    droplistVarSample <- c()
  }
  
  droplist <- unique(c(dropListVar, droplistFC, droplistFarmsLaplace, droplistVarSample))
  eset <- eset[!rownames(fData(eset)) %in% droplist, ]
  eset
}

