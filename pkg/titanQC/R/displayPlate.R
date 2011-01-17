#' Draw a Summary Display for a Titan Plate
#' @param eset ExpressionSet object with additional plate number and position information
#'  in the pData; for the Incubator plates the column names of the pData need to be
#'  incubatorPlateNo,  incubatorColumn, incubatorRow; for Titan plates the column names of the
#'  pData need to be titanPlateNo,  titanColumn, titanRow
#' @param statistic summary statistic or information to be displayed for a well on the plate;
#' can be either a name of a pData variable or a vector of values (one for each well) respecting 
#' the row order of pData variables  
#' @param type type of plate; one of 'GeneTitan' (default) or 'Incubator'
#' @param stepSize step size used to discretize the range of summary values to be displayed; the
#' step size is reflected in the distance between ticks on the heatmap legend and each step is 
#' associated with a different color level in the heatmap 
#' @param main main title for the plot; if NULL (default) 'Plate n' is used where n is the plate
#'   number; if only one plate is present no title is displayed; if the user provides a main title
#'  (string) and there is more than one plate present, the plate numbers are added to the title
#'  between round brackets, as in 'Main Title (Plate n)'.   
#' @param legend title of the legend (character of length one); if NULL (default) no legend is 
#'   drawn 
#' @param breaks vector of break points supplied to the image function; if NULL (default) the break
#'   points are generated automatically
#' @param text name of a pData variable to be displayed; each cell of the plate grid will be annotated
#'   using the value in the corresponding pData variable
#' @param textsize cex argument to be passed to the text function which will display the values specified in the
#'  'text' argument
#' @param textcol col argument to be passed to the text function which will display the values specified in the
#'  'text' argument
#' @param color color scheme to be used for the cells of the grid;  
#' if 'heat' heat colors are used, if 'redgreen' the colors go from green (low) to black (medium) 
#' to red (high); if 'random' rainbow colors are used (typically used for discrete variables), 
#' if 'white' the background color of the cells is 'white'; alternatively the name of a pData variable
#' can be passed and will be used to set the colors of the cells; for any color scheme all missing 
#' values are set to 'black' 
#' @param filePrefix prefix used for all pdf files; one pdf file is generated for each plate
#' @return one pdf file is generated for each plate
#' @import Biobase
#' @export
displayPlate <- function(eset, 
    statistic,   # either a name of a variable in pData or a vector with names = samplenames
    type = "GeneTitan",
    stepSize = 0.01, 
    main = NULL,
    legend = NULL, 
    breaks = NULL,
    text = NULL,
    textsize = 0.7,
    textcol = "black",
    color = "heat", # redgreen, random, white, or variable from pData
    filePrefix = "display"
){
  
  # if the variable is the expression of a gene, it is attached to the pData
  if (!is.character(statistic)){
    pData(eset) <- cbind(pData(eset), stat = statistic)
    eset$stat <- as.numeric(as.character(eset$stat))
    statistic <- "stat"
  }
  
  if (type == "GeneTitan") {
    heatinfoAll <- pData(eset)[, c(statistic, "titanPlateNo", "titanColumn", "titanRow")]
    colnames(heatinfoAll) <- c(statistic, "plate", "col", "row")
  }  else {
    heatinfoAll <- pData(eset)[, c(statistic, "incubatorPlateNo", "incubatorColumn", "incubatorRow")]
    colnames(heatinfoAll) <- c(statistic, "plate", "col", "row")
  }
  
  if (!is.null(text)) {
    heatinfoAll <- cbind(heatinfoAll, text = pData(eset)[, text])
  }
  
  if (color %in% colnames(pData(eset))) {
    heatinfoAll <- cbind(heatinfoAll, color = pData(eset)[, color])
  }
  
  # !!!!!!!!!!!! add following section
  if (is.null(breaks) & !color %in% c("heat", "redgreen")) {
    heatinfoAll[, statistic] <- as.numeric(as.factor(heatinfoAll[, statistic]))
  }
    
  for (plate in unique(heatinfoAll$plate)) {
    
    heatinfo <- heatinfoAll[heatinfoAll$plate==plate,] 
 
    #dev.new(width = 12, height = 6)
    pdf(file = paste(filePrefix, "Plate", plate, ".pdf",sep=""), width = 12, height = 6) 
    if (!is.null(legend)) {
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2), 5, 2, byrow = TRUE))
      par(mar = c(1, 6, 4, 2) + 0.1)
    }
    
    #complete plate 
    heatinfo <- cbind(heatinfo, loc = paste(heatinfo$row, heatinfo$col, sep = ""))
    heatinfo$loc <- as.character(heatinfo$loc)
    completeplate <- matrix("0", ncol = 3, nrow = 96)
    rowt <- c("A", "B", "C", "D", "E", "F", "G", "H")
    colt <- c(1:12)
    for (rownr in 1:length(rowt)) {
      completeplate[((rownr - 1) * 12 + 1):(12 * rownr), 1] <- paste(rowt[rownr], colt, sep = "")
      completeplate[((rownr - 1) * 12 + 1):(12 * rownr), 2] <- rowt[rownr]
      completeplate[((rownr - 1) * 12 + 1):(12 * rownr), 3] <- colt
    }
    colnames(completeplate) <- c("loc", "Trow", "Tcol")
    heatinfo <- merge(x = heatinfo, y = completeplate, by = "loc", all.y = TRUE)
    heatinfo$plate[is.na(heatinfo$col)] <- plate
    heatinfo$row[is.na(heatinfo$row)] <- as.character(heatinfo$Trow)[is.na(heatinfo$row)]
    heatinfo$col[is.na(heatinfo$col)] <- as.character(heatinfo$Tcol)[is.na(heatinfo$col)]
    heatinfo <- heatinfo[, !colnames(heatinfo)%in%c("loc", "Trow", "Tcol")]
    
    #breakvec
    factor <- 1 / (10 ^ floor(log10(stepSize)))
    if (is.null(breaks)) {
      if (color %in% c("heat", "redgreen")) {
        mind <- floor(min(heatinfoAll[, statistic], na.rm = TRUE) * factor) / factor
        maxd <- ceiling(max(heatinfoAll[, statistic], na.rm = TRUE) * factor) / factor
        breakvecP <- unique(c(seq(mind - stepSize, maxd, by = stepSize), maxd))
      }  else {
        tmp <- sort(unique(heatinfoAll[, statistic]))
        breakvecP <- seq(min(tmp) - 0.5, max(tmp) + 0.5, by = 1)
      }
    } else {
      breakvecP <- breaks
    }
    
    breakvecP <- c(min(breakvecP) - 1, breakvecP)
    heatinfo[is.na(heatinfo[, statistic]), statistic] <- min(breakvecP)
    lbv <- length(breakvecP)
    n <- lbv - 2
    
    # color
    if (color %in% colnames(pData(eset))) {
      tmp <- unique(heatinfoAll[, c(statistic, "color")])
      tmp <- tmp[order(tmp[, statistic]), ]
      colplate <- c("black", as.character(tmp$color[!is.na(tmp$color)]))
    }
    
    if (color == "random") {
      colplate <- c("black", rainbow(n))
    }
    
    if (color == "white") {
      colplate <- c("white", rep("white", n))
    }
    
    if (color == "heat") {
      colplate <- c("black",heat.colors(n, alpha = 1)[n:1])
    }
    
    if (color == "redgreen") {
      if (n/2 == floor(n/2)){
        colplate <- c("white", rgb(red = 0, green = ((n/2):1) / (n/2), blue = 0), rgb((1:(n/2))/(n/2), green = 0, blue = 0))
      }
      else{
        colplate <- c("white", rgb(red = 0, green = (floor(n/2):1) / floor(n/2), blue = 0), rgb((0:floor(n/2)) / floor(n/2), green = 0, blue = 0))
      }
    }
     
    heatinfo <- heatinfo[order(as.numeric(heatinfo$col), heatinfo$row), ]
    hmat <- matrix(heatinfo[, statistic], nrow=8, ncol=12, byrow = FALSE)
    z <- t(hmat)
    z <- z[, 8:1]
    rownames(z) <- c(1:12)
    colnames(z) <- c("H", "G", "F", "E", "D", "C", "B", "A")
    x <- c(1:nrow(z))
    y <- c(1:ncol(z))
    image(x, y, z, axes = FALSE, xlab = "", ylab = "", col = colplate, breaks = breakvecP)
    axis(1, at = x, labels = rownames(z), las = 1, cex.axis = 1.5, tick = FALSE)
    axis(2, at = y, labels = colnames(z), las = 1, cex.axis = 1.5, tick = FALSE)
    abline(h = c(0.5,8.5), v = c(0.5, 12.5)) 
    
    if (color == "white") {
      abline(h = c(0.5:8.5), v = c(0.5:12.5))
    } else {
      abline(h = c(0.5,8.5), v = c(0.5,12.5))
     }
    
    if (!is.null(text)){
      textloc <- cbind(rep(c(8:1), 12), sort(rep(c(1:12), 8)))
      text(textloc[, 2], textloc[, 1], heatinfo$text, cex = textsize, col = textcol)
    }
    
    if (length(unique(heatinfoAll$plate))== 1) {
      title(main,cex.main = 3)
    } else {
      title(paste(main, " (Plate ", plate, ")", sep = ""), cex.main = 3)
    }
    
    if (!is.null(legend)) {
      par(mar = c(5, 6, 3, 2) + 0.1)
      image(x = c(1:(length(breakvecP) - 1) - 0.5), y = c(1), z = as.matrix(breakvecP[c(-1, -2)]), axes = FALSE, 
          xlab = legend, cex.lab = 2, ylab = "", col = colplate[-1], breaks = breakvecP[-1])
      axis(1, at = c(1:(length(breakvecP) -1 )) - 0.5, labels = breakvecP[-1], las = 1, cex.axis = 1.5)
    }
    dev.off()
  }
}

