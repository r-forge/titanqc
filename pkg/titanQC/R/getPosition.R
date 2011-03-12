#' Get the Position on the Plate for One Single CEL File
#' @param fileName path to CEL file
#' @return data frame with columns titanPlateBarcode, titanRow, titanColumn
#' @export
getPosition <- function(fileName){
  
  celFileCon <- file(fileName, "rb")
  tmp <- readBin(celFileCon, n = 12000, what = "raw")
  # ascii only, no control characters, remove spaces
  celFileHeader <- rawToChar(tmp[!(tmp %in% as.character(as.hexmode(c(0:32, 127)))) & tmp %in% as.character(as.hexmode(0:128))])
  close(celFileCon)
  
  barcodeId <-strapply(celFileHeader, 
      "affymetrix-plate-barcode,(.{22})text",
      simplify = c)
  
  posData <- strapply(celFileHeader, 
      "affymetrix-plate-peg-wellposition(...)text",
      simplify = c)
  rowPos <- substr(posData, start = 1, stop = 1)
  colPos <- as.numeric(substr(posData, start = 2, stop = nchar(posData)))
  posDf <- data.frame(titanPlateBarcode = barcodeId, titanRow = rowPos, 
      titanColumn = colPos, row.names = fileName)
  return(posDf)
}

#' Get the Position of a Set of CEL Files on a Titan Plate
#' @param celFiles character vector of the paths of the CEL files for which the position on the
#'  plate needs to be identified
#' @return  data frame with columns titanPlateBarcode, titanPlateNo, titanRow, titanColumn
#' @note It is assumed that all CEL file names of a given plate are contiguous numbers and either
#'  all greater than or all smaller than the CEL file names (numbers) of another plate; this assumption
#'  is used to assign an integer number to a plate, i.e. a simpler identifier than the full plate barcode
#'  number
#' @export
getCelFilePosition <- function(celFiles = NULL){
  if (is.null(celFiles)){
    celFiles <- list.files(pattern = "\\.CEL$")
  }
  sampleName <- gsub("(.+)\\.CEL", "\\1", basename(celFiles))
  celFileList <- vector(mode = "list", length = length(celFiles))
  for (iFile in seq_along(celFiles)){
    celFileCon <- file(celFiles[iFile], "rb")
    tmp <- readBin(celFileCon, n = 12000, what = "raw")
    tmp <- rawToChar(tmp[!(tmp %in% as.character(as.hexmode(c(0:32, 127)))) & tmp %in% as.character(as.hexmode(0:128))])
    celFileList[[iFile]] <- tmp
    close(celFileCon)
  }
  
  barcodeIdFunction <- function(x){
    barcodeId <- strapply(x, 
        "affymetrix-plate-barcode,(.{22})text", # ID 22 characters (confirmed by Affy)
        simplify = c)
    return(barcodeId)
  }
  
  rowPosFunction <- function(x){
    posData <- strapply(x, 
        "affymetrix-plate-peg-wellposition(...)text",
        simplify = c)
    rowPos <- substr(posData, start = 1, stop = 1)
    return(rowPos)
  }
  
  colPosFunction <- function(x){
    posData <- strapply(x, 
        "affymetrix-plate-peg-wellposition(...)text",
        simplify = c)
    colPos <- as.numeric(substr(posData, start = 2, stop = nchar(posData)))  
    return(colPos)
  }
  
  rowPos <- sapply(celFileList, rowPosFunction)
  if (length(errorRowPos <- which(!(rowPos %in% LETTERS[1:8]))))
    stop(paste("No valid row position (A-H) in CEL file ", 
            paste(celFiles[errorRowPos], collapse = ", "), sep = ""))
  
  colPos <- sapply(celFileList, colPosFunction)
  if (length(errorColPos <- which(colPos < 1 | colPos > 12)))
    stop(paste("No valid column position (1-12) in CEL file ",
            paste(celFiles[errorColPos], collapse = ", "), sep = ""))
  
  barcodeId <- sapply(celFileList, barcodeIdFunction)
  
  posDf <- data.frame(sampleName = sampleName, titanPlateBarcode = barcodeId, titanRow = rowPos, 
      titanColumn = colPos, row.names = basename(celFiles), stringsAsFactors = FALSE)

  plateNo <- cbind(titanPlateBarcode = unique(posDf$titanPlateBarcode), titanPlateNo = c(1:length(unique(posDf$titanPlateBarcode)))) 
  posDf <- merge(posDf, plateNo, by = "titanPlateBarcode")
  posDf$titanPlateNo <- as.character(posDf$titanPlateNo)
		  
  return(posDf)
}
