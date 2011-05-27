\name{getCelFilePosition}
\alias{getCelFilePosition}
\title{Get the position of a set of CEL files on a GeneTitan plate...}
\usage{getCelFilePosition(celFiles)
}
\description{Get the position of a set of CEL files on a GeneTitan plate}
\value{data frame with columns titanPlateBarcode, titanPlateNo, titanRow, titanColumn}
\note{It is assumed that all CEL file names of a given plate are contiguous numbers and either
all greater than or all smaller than the CEL file names (numbers) of another plate; this assumption
is used to assign an integer number to a plate, i.e. a simpler identifier than the full plate barcode
number}
\arguments{\item{celFiles}{character vector of the paths of the CEL files for which the position on the
plate needs to be identified}
}

