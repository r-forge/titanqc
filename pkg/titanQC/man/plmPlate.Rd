\name{plmPlate}
\alias{plmPlate}
\title{Plot affyPLM QC images on a grid respecting the plate layout...}
\usage{plmPlate(x, celFilePositions, type="resids", use.log=TRUE, standardize=FALSE,
    col, addSampleName=FALSE, sampleNameAffix="")
}
\description{Plot affyPLM QC images on a grid respecting the plate layout}
\value{no return value; a graph is drawn to the current device}
\seealso{\code{\link[affyPLM]{PLMset-class}}}
\author{Code from the image methods of the affyPLM package by Ben Bolstad, adapted
to titan QC purposes by Tobias Verbeke}
\note{currently only type 'resids' from the affyPLM package is supported}
\references{TODO}
\arguments{\item{x}{object of class 'PLMset' as produced by 'fitPLM' of the 'affyPLM' package}
\item{celFilePositions}{dataframe as produced by getCelFilePosition; this data
frame provides all necessary information on location on the plate of 
the samples}
\item{type}{type of residuals to plot; currently only type 'resids'
is supported}
\item{use.log}{logical; defaults to \code{TRUE}}
\item{standardize}{logical; see \code{\link[affyPLM]{PLMset-class}}}
\item{col}{color palette to use for the coloring of residuals; if NULL (default)
the coloring of the affyPLM package is used (pseudoPalette)}
\item{addSampleName}{add the sample names on top of each single
well image for the plate; defaults to FALSE}
\item{sampleNameAffix}{affix used to comply with a legacy application
where '.CEL' is added to each sample name, defaults to '' in which
case sample names are expected to be the CEL file names (without the
.CEL extension)}
}

