\name{getCelFilePosition}
\alias{getCelFilePosition}
\title{Get the Position of a Set of CEL Files on a Titan Plate...}
\usage{getCelFilePosition(celFiles)
}
\description{Get the Position of a Set of CEL Files on a Titan Plate}
\value{data frame with columns titanPlateBarcode, titanPlateNo, titanRow, titanColumn}
\note{It is assumed that all CEL file names of a given plate are contiguous numbers and either
all greater than or all smaller than the CEL file names (numbers) of another plate; this assumption
is used to assign an integer number to a plate, i.e. a simpler identifier than the full plate barcode
number}
\arguments{\item{celFiles}{character vector of the paths of the CEL files for which the position on the
plate needs to be identified}
}

