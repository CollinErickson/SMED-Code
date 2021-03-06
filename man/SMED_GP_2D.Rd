% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMED_GP_2D.R
\name{SMED_GP_2D}
\alias{SMED_GP_2D}
\title{SMED GP 2D}
\usage{
SMED_GP_2D(f, n0 = 10, n = 10, GP.package = "", opt.method = "genoud",
  contour.fit = 1, continue.option = F, max.time = NULL)
}
\arguments{
\item{f}{Function}

\item{n0}{Number of initial points}

\item{n}{Number of points to select}

\item{GP.package}{GP package}

\item{opt.method}{Optimization method. Current options are 'genoud', 'LHS'.}

\item{contour.fit}{Number of iterations after which to plot model contour plot.}

\item{continue.option}{If TRUE, you are given choice to continue sampling}

\item{max.time}{max.time for GenSA optimization for each point}
}
\value{
X, design matrix, and Y, function values at design points.
}
\description{
SMED GP 2D
}
\examples{
SMED_GP_2D(TestFunctions::banana,10,10, GP.package='laGP', contour.fit=0, max.time=.1)
}
