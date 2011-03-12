#' Plot affyPLM QC Images on a Grid Respecting the Plate Layout
#' 
#' @param x object of class 'PLMset' as produced by 'fitPLM' of the 'affyPLM' package
#' @param celFilePositions dataframe as produced by getCelFilePosition; 
#'   providing all necessary information on location on the plate of 
#'   the samples
#' @param type type of residuals to plot; currently only type 'resids'
#'   is supported
#' @param use.log logical; defaults to \code{TRUE}
#' @param standardize logical; see \code{\link[affyPLM]{PLMset-class}}
#' @param col color palette to use for the coloring of residuals
#' @param addSampleName add the sample names on top of each single
#'   well image for the plate; defaults to FALSE
#' @param sampleNameAffix affix used to comply with a legacy application
#'   where '.CEL' is added to each sample name, defaults to '' in which
#'   case sample names are expected to be the CEL file names (without the
#'   .CEL extension)
#' @return no return value; an graph is drawn to the current device
#' @seealso \code{\link[affyPLM]{PLMset-class}}
#' @author Code from the image methods of the affyPLM package by Ben Bolker, adapted
#'   to titan QC purposes by Tobias Verbeke
#' @note currently only type 'resids' from the affyPLM package is supported
#' @import affyPLM
#' @export
plmPlate <- function(x, celFilePositions, type = "resids", use.log=TRUE,
    standardize=FALSE, col=NULL, addSampleName = FALSE, sampleNameAffix = ""){
  
  
  if (type != "resids")
    stop("Currently only type 'resids' is supported")
  
  col.resids <- if (is.null(col)) pseudoPalette(low = "blue", high = "red", mid = "white") else col
  
  pm.index <- unique(unlist(indexProbes(x, "pm",row.names(coefs(x)))))
  rows <- x@nrow
  cols <- x@ncol
  pm.x.locs <- pm.index%%rows
  pm.x.locs[pm.x.locs == 0] <- rows
  pm.y.locs <- pm.index%/%rows + 1
  xycoor <- matrix(cbind(pm.x.locs,pm.y.locs),ncol=2)
  
  mm.index <-  unique(unlist(indexProbes(x, "mm",row.names(coefs(x)))))
  mm.x.locs <- mm.index%%rows
  mm.x.locs[mm.x.locs == 0] <- rows
  mm.y.locs <- mm.index%/%rows + 1
  
  xycoor2 <-matrix(cbind(mm.x.locs,mm.y.locs),ncol=2) ##xycoor## matrix(cbind(pm.x.locs,pm.y.locs+1),ncol=2)
  
  if (any(is.na(xycoor2))){
    xycoor2 <-xycoor
  }
  
  
  if (is.element(type, c("resids","pos.resids","neg.resids","sign.resids"))){
    if (any(dim(x@residuals[[1]]) ==0) & any(dim(x@residuals[[2]]) ==0)){
      stop("Sorry this PLMset does not appear to have residuals\n");
    }
    
    if (standardize & type == "resids"){
      if (x@model.description$R.model$response.variable == 0){
        resid.range <- c(-4,4)
      } else if (x@model.description$R.model$response.variable == -1){
        resid.range <- range(resid(x,standardize)[[2]])
      } else if (x@model.description$R.model$response.variable == 1){
        resid.range <- range(resid(x,standardize)[[1]])
      }
      
    } else {
      if (x@model.description$R.model$response.variable == 0){
        resid.range1 <- range(x@residuals[[1]])
        resid.range2 <- range(x@residuals[[2]])
        resid.range <- resid.range1
        resid.range[1] <- min(resid.range1 , resid.range2)
        resid.range[2] <- max(resid.range1 , resid.range2)
      } else if (x@model.description$R.model$response.variable == -1){
        resid.range <- range(x@residuals[[2]])
      } else if (x@model.description$R.model$response.variable == 1){
        resid.range <- range(x@residuals[[1]])
      }
    }
  }
  
  ### create list of residual matrices
  residsMatrixList <- vector(mode = "list", length = nrow(celFilePositions))
  
  for (iSample in seq(from = 1, to = nrow(celFilePositions))){
    
    residsmatrix <- matrix(nrow=rows, ncol=cols)
    if (standardize){
      if (x@model.description$R.model$response.variable == 0){
        residsmatrix[xycoor]<- resid(x,standardize)[[1]][,iSample]
        residsmatrix[xycoor2]<- resid(x,standardize)[[2]][,iSample]
      } else if  (x@model.description$R.model$response.variable == -1){
        residsmatrix[xycoor]<- resid(x,standardize)[[2]][,iSample]
        residsmatrix[xycoor2]<- resid(x,standardize)[[2]][,iSample]
      } else if (x@model.description$R.model$response.variable == 1){
        residsmatrix[xycoor]<- resid(x,standardize)[[1]][,iSample]
        residsmatrix[xycoor2]<- resid(x,standardize)[[1]][,iSample]
      }
    } else {
      if (x@model.description$R.model$response.variable == 0){
        residsmatrix[xycoor]<- x@residuals[[1]][,iSample]
        residsmatrix[xycoor2]<- x@residuals[[2]][,iSample]
      } else if (x@model.description$R.model$response.variable == -1){
        residsmatrix[xycoor]<- x@residuals[[2]][,iSample]
        residsmatrix[xycoor2]<- x@residuals[[2]][,iSample]
      } else if (x@model.description$R.model$response.variable == 1){
        residsmatrix[xycoor]<- x@residuals[[1]][,iSample]
        residsmatrix[xycoor2]<- x@residuals[[1]][,iSample]
      }
      
    }
    # this line
    # flips the matrix around so it is correct
    residsmatrix <- as.matrix(rev(as.data.frame(residsmatrix)))
    
    if (use.log)
      residsmatrix <- sign(residsmatrix)*log2(abs(residsmatrix)+1)
    
    residsMatrixList[[iSample]] <- residsmatrix
  }
  names(residsMatrixList) <- sampleNames(x)
  
  ### actual plotting
  
  # order position information
  celFilePositions$titanRowNumber <- as.numeric(factor(celFilePositions$titanRow, 
          levels = LETTERS[1:8]))
  
  celFilePositions <- celFilePositions[order(celFilePositions$titanColumn, 
          celFilePositions$titanRowNumber), ]
  
  # reorder list of residual matrices before plotting
  celFilePositions$sampleNameCEL <- paste(celFilePositions$sampleName, sampleNameAffix, sep = "")
  residsMatrixList <- residsMatrixList[celFilePositions$sampleNameCEL]
  
  # prepare plot layout (as not all cells on the grid will be taken)
  presenceMatrix <- matrix(0, nrow = 8, ncol = 12)
  interactionVariable <- interaction(celFilePositions$titanRowNumber, celFilePositions$titanColumn)
  for (iRow in 1:8){
    for (iCol in 1:12){
      rowPos <- which(paste(iRow, iCol, sep = ".") == interactionVariable)
      presenceMatrix[iRow, iCol] <- if (length(rowPos)) 1 else 0
    }
  }
  zeroPositions <- as.logical(1-presenceMatrix)
  plotOrder <- seq(from = 1, to = 96 - sum(zeroPositions))
  
  presenceMatrix[!zeroPositions] <- plotOrder
  
  # add space for row and column annotation
  nArrays <- length(residsMatrixList)
  presenceMatrix <- cbind((nArrays+1):(nArrays+8), presenceMatrix)
  presenceMatrix <- rbind(c(0, (nArrays+9):(nArrays+9+11)), presenceMatrix)
  
  plateLayoutLayout <- layout(presenceMatrix)
  
  
  op <- par(mar = c(0, 0, 0, 0))
  for (iMatrix in seq(length.out = nArrays)){
    if (use.log){
      image(residsMatrixList[[iMatrix]], col = col.resids, xaxt = "n",
          yaxt = "n",
          zlim = c(-max(log2(abs(resid.range)+1)), max(log2(abs(resid.range)+1))))
      if (addSampleName)
        legend("center", legend = names(residsMatrixList)[[iMatrix]], bty = "n", cex = 3)
    } else {
      image(residsMatrixList[[iMatrix]], col = col.resids, xaxt = "n",
          yaxt = "n",
          zlim = c(-max(abs(resid.range)), max(abs(resid.range))))
      if (addSampleName)
        legend("center", legend = names(residsMatrixList)[[iMatrix]], bty = "n", cex = 3)
    }  
  }
  par(op)
  
  # add row annotation
  for (i in 1:8){
    plot(1, type = "n", ann =  FALSE, axes = FALSE)
    legend("center", legend = LETTERS[i], cex = 3, bty = "n")
  }
  
  # add column annotation
  for (j in 1:12){
    plot(1, type = "n", ann =  FALSE, axes = FALSE)
    legend("center", legend = j, cex = 3, bty = "n")
  }
  
}



