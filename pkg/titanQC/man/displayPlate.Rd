\name{displayPlate}
\alias{displayPlate}
\title{Draw a summary display for a GeneTitan plate...}
\usage{displayPlate(eset, statistic, type="GeneTitan", stepSize=0.01, main, legend, breaks,
    text, textsize=0.7, textcol="black", color="heat",
    filePrefix="displayPlate")
}
\description{Draw a summary display for a GeneTitan plate}
\note{If the profiled samples originate from an in vitro experiment, the information of the cell culture plates (Incubator plates) 
can be added to the pData of the ExpressionSet object}
\value{no return value; one pdf file is generated for each plate}
\arguments{\item{eset}{ExpressionSet object with additional plate number and position information
in the pData; for the Incubator plates (cf. note)  the column names of the pData need to be
incubatorPlateNo,  incubatorColumn, incubatorRow; for GeneTitan plates the column names of the
pData need to be titanPlateNo,  titanColumn, titanRow}
\item{statistic}{summary statistic or information to be displayed for a well on the plate;
can be either a name of a pData variable or a vector of values (one for each well) respecting 
the row order of pData variables}
\item{type}{type of plate; one of 'GeneTitan' (default) or 'Incubator'}
\item{stepSize}{step size used to discretize the range of summary values to be displayed; the
step size is reflected in the distance between ticks on the heatmap legend and each step is 
associated with a different color level in the heatmap}
\item{main}{main title for the plot; if NULL (default) 'Plate n' is used where n is the plate
number; if only one plate is present no title is displayed; if the user provides a main title
(string) and there is more than one plate present, the plate numbers are added to the title
between round brackets, as in 'Main Title (Plate n)'.}
\item{legend}{title of the legend (character of length one); if NULL (default) no legend is 
drawn}
\item{breaks}{vector of break points supplied to the image function; if NULL (default) the break
points are generated automatically}
\item{text}{name of a pData variable to be displayed; each cell of the plate grid will be annotated
using the value in the corresponding pData variable}
\item{textsize}{cex argument to be passed to the text function which will display the values specified in the
'text' argument}
\item{textcol}{col argument to be passed to the text function which will display the values specified in the
'text' argument}
\item{color}{color scheme to be used for the cells of the grid;  
if 'heat' heat colors are used, if 'redgreen' the colors go from green (low) to black (medium) 
to red (high); if 'random' rainbow colors are used (typically used for discrete variables), 
if 'white' the background color of the cells is 'white'; alternatively the name of a pData variable
can be passed and will be used to set the colors of the cells; for any color scheme all missing 
values are set to 'black'}
\item{filePrefix}{prefix used for all pdf files; one pdf file is generated for each plate; defaults to "displayPlate"}
}

