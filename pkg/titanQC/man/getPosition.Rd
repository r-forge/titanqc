\name{getPosition}
\alias{getPosition}
\title{Get the position on the plate for one single CEL File...}
\usage{getPosition(fileName)
}
\description{Get the position on the plate for one single CEL File}
\value{data frame with columns 
\item{titanPlateBarcode}{character vector with the plate barcodes}
\item{sampleName}{character vector of sample names, concretely speaking the name of the CEL file without the .CEL extension}
\item{titanRow}{character vector indicating the row on the plate, from "A" to "H"}
\item{titanColumn}{numeric vector indicating the column on the plate, an integer between 1 and 12}
\item{titanPlateNo}{character vector with integer numbers "1" up to the number of plates}}
\arguments{\item{fileName}{path to CEL file}
}

