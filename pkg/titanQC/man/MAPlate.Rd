\name{MAPlate}
\alias{MAPlate}
\title{Generate MA plots for all individual wells on a GeneTitan plate...}
\usage{MAPlate(eset, filePrefix="MAPlate", alpha, gradient=FALSE, title)
}
\description{Generate MA plots for all individual wells on a GeneTitan plate}
\note{the color codes correspond to the sampleColor column of the pData of the ExpressionSet passed to argument
'eset'; if no sampleColor column is present, the data for all samples will be plotted in 'blue'}
\value{no return value; the function generates as many png files as there are plates contained in the ExpressionSet object 'eset'}
\references{TODO}
\arguments{\item{eset}{ExpressionSet object as produced in the 'preprocessing' function}
\item{filePrefix}{prefix that will be used to generate the output files; for plate <i>, the default
name will be <filePrefix>Plate<i>.png; the default file prefix is "MAPlate"}
\item{alpha}{alpha value for the alphahull function (numeric of length one); when NULL (default) 
no contours are plotted}
\item{gradient}{if TRUE a smoothscatter plot is used to depict the raw data points; if FALSE the points
are plotted as is; defaults to FALSE}
\item{title}{name of the variable in the pData of the eset ExpresionSet that can be used to provide titles
on the individual well plots; if NULL no titles are plotted for the individual wells}
}

