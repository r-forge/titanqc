#' Patched version of the plot.ahull function from the alphahull package
#' This version TODO
#' @param x x
#' @param add add
#' @param do.shape do.shape 
#' @param wlines wlines
#' @param wpoints wpoints
#' @param number number
#' @param col col
#' @param xlim xlim
#' @param ylim ylim
#' @param lwd lwd
#' @param ... further arguments
#' @return plot is displayed on the current device
#' @export
plot.ahullS <- function (x, add = FALSE, do.shape = FALSE, wlines = c("none", 
        "both", "del", "vor"), wpoints = TRUE, number = FALSE, col = NULL, 
    xlim = NULL, ylim = NULL, lwd = NULL, ...) 
{
  wlines <- match.arg(wlines)
  if (is.null(class(x)) || class(x) != "ahull") {
    cat("Argument is not of class ahull.\n")
    return(invisible())
  }
  if (is.null(col)) {
    col <- c(1, 1, 1, 1, 1, 1)
  }
  else {
    col <- rep(col, length.out = 6)
  }
  if (is.null(lwd)) {
    lwd <- c(1, 1, 2)
  }
  else {
    lwd <- rep(lwd, length.out = 3)
  }
  wlines <- match.arg(wlines)
  plot.dd <- switch(wlines, none = TRUE, both = FALSE, del = FALSE, 
      vor = FALSE)
  if (do.shape) {
    plot.ashape(x$ashape.obj, add = add, wlines = wlines, 
        wpoints = wpoints, number = number, col = col[2:6], 
        xlim = xlim, ylim = ylim, lwd = lwd[1:2], ...)
  }
  else {
    if (plot.dd) {
      if (!add) {
        if (is.null(xlim)) 
          xlim <- range(x$ashape.obj$x[, 1])
        if (is.null(ylim)) 
          ylim <- range(x$ashape.obj$x[, 2])
        plot(0, 0, type = "n", xlim = xlim, ylim = ylim, 
            axes = FALSE, ...)
        axis(side = 1)
        axis(side = 2)
      }
      if (wpoints) {
        points(x$ashape.obj$x, col = col[3], ...)
      }
      if (number) {
        xoff <- 0.02 * diff(range(x$ashape.obj$x[, 1]))
        yoff <- 0.02 * diff(range(x$ashape.obj$x[, 2]))
        text(x$ashape.obj$x[, 1] + xoff, x$ashape.obj$x[, 
                2] + yoff, 1:(dim(x$ashape.obj$x)[1]), col = col[6], 
            ...)
      }
    }
    else {
      plot.delvor(x$ashape.obj$delvor.obj, add = add, wlines = wlines, 
          wpoints = wpoints, number = number, col = col[3:6], 
          lwd = lwd[1], xlim = xlim, ylim = ylim, ...)
    }
  }
  arcs <- which(x$arcs[, 3] > 0)
  if (length(arcs) > 0) {
    for (i in arcs) {
      arc(x$arcs[i, 1:2], x$arcs[i, 3], x$arcs[i, 4:5], 
          x$arcs[i, 6], col = col[1], lwd = lwd[3])
    }
  }
}


#' Generate a MA plots for all individual wells on the plate
#' @param eset ExpressionSet object as produced in the 'preprocessing' function 
#' @param filePrefix prefix that will be used to generate the output files; for plate <i>, the default
#' name will be <filePrefix>Plate<i>.png; the default file prefix is "MAPlate" 
#' @param alpha alpha value for the alphahull function (numeric of length one); when NULL no contours are
#'   plotted
#' @param gradient if TRUE a smoothscatter plot is used to depict the raw data points; if FALSE the points
#'    are plotted as is; defaults to FALSE
#' @param title name of the variable in the pData of the eset ExpresionSet that can be used to provide titles
#' on the individual well plots; if NULL no titles are plotted for the individual wells
#' @note the color codes correspond to the sampleColor column of the pData of the ExpressionSet passed to argument
#' 'eset'; if no sampleColor column is present, the data for all samples will be plotted in 'blue'
#' @return the function generates as many png files as there are plates contained in the ExpressionSet object 'eset'
#' @export
MAPlate <- function(eset, filePrefix = "MAPlate", alpha = NULL, gradient = FALSE, title = NULL) {
  
  if (is.null(pData(eset)$sampleColor))
    pData(eset)$sampleColor <- "blue"
  
  for (plate in sort(unique(eset$titanPlateNo))) {
    
    esetTemp <- eset[, eset$titanPlateNo == plate]
    
    png(file = paste(filePrefix, "Plate", plate, ".png", sep = ""), width = 12, height = 6, units = "in", res=72) 
    par(mfrow = c(8,12), mar = c(0.5, 0.5, 1, 0.5))
    x <- apply(exprs(esetTemp), 1, median)
    y <- exprs(esetTemp)
    xlim <- summary(unmatrix(x + y)/2)[c(1, 6)]
    ylim <- summary(unmatrix(y - x))[c(1, 6)]   
    esetTemp <- esetTemp[, order(esetTemp$titanRow, esetTemp$titanColumn)]  
    
    plotSmoothScat <- function(x, y, col1, col2, xlim, ylim){
      if (!is.null(alpha)){
        datamatrix <- as.data.frame(cbind(x, y))
        Q05x <- quantile(datamatrix$x, probs = 0.05)
        Q95x <- quantile(datamatrix$x, probs = 0.95)
        Q05y <- quantile(datamatrix$y, probs = 0.05)
        Q95y <- quantile(datamatrix$y, probs = 0.95)
        datamatrix <- datamatrix[ datamatrix$x < Q05x | datamatrix$x > Q95x | datamatrix$y < Q05y | datamatrix$y > Q95y, ]
        ahull.obj <- ahull(unique(datamatrix), alpha = alpha)
      }
      plot(x, y, axes = FALSE, type = "n", xlim = xlim, ylim = ylim)
      axis(1, lwd = 0, labels = FALSE)
      axis(2, lwd = 0, labels = FALSE)
      box(bty = 'l', lwd = 1.5)
      blue.grad <- colorpanel(100, "white", col1, col2)[round(seq(25,96,length=10))]
      blue.grad2 <- densCols(x, y, nbin = 128, colramp = colorRampPalette(blue.grad))
      points(x, y, pch = 20, cex = 1, col = blue.grad2)
      if (!is.null(alpha)){
        plot.ahullS(ahull.obj, add = TRUE, wpoints = FALSE, wlines = "none", lwd = 0.1)
      }
    }
    
    for (row in c("A", "B", "C", "D", "E", "F", "G", "H")) {
      for (col in c(1:12)) {
        y <- exprs(esetTemp)[,esetTemp$titanRow == row & esetTemp$titanColumn == col]
        if(sum(esetTemp$titanRow == row & esetTemp$titanColumn == col) == 0) {
          plot(0, 0, type = "n", axes = FALSE)
        } else {
          meanv <- (x + y)/2
          logr <- y - x
          col1 <- esetTemp$sampleColor[esetTemp$titanRow == row & esetTemp$titanColumn == col]
          if(gradient == TRUE) {
            basecol <- col2rgb(col1)
            basecol[basecol == max(basecol)] <- basecol[basecol == max(basecol)] - 0.30 * basecol[basecol == max(basecol)]
            col2 <- rgb((basecol[1, 1])/255, (basecol[2, 1])/255, (basecol[3, 1])/255)
          } else {
            col2 <- col1
          }
          plotSmoothScat(x = meanv, y = logr, col1 = col1, col2 = col2, xlim = xlim, ylim = ylim)
          lines(lowess(meanv, logr), col = "black", lwd=2)
          if(!is.null(title)) {
            title(pData(esetTemp)[esetTemp$titanRow == row & esetTemp$titanColumn == col, title])
          }
        }
      }
    }  
    dev.off()
  } 
}

