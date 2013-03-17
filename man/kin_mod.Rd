\name{kin_mod}
\alias{kin_mod}
\title{kinetic models}
\usage{
  kin_mod(P, inside = FALSE, plot = TRUE, plottitle,
    pnames = NULL, ...)
}
\arguments{
  \item{P}{parameters}

  \item{inside}{logical, whether the function is evaluated
  inside a fit function.}

  \item{plot}{logical, whether to make a plot}

  \item{plottitle}{the title of the plot}

  \item{pnames}{names of the parameters, nls function will
  forget them!}

  \item{...}{other parameters to be passed into plot}
}
\value{
  the calculated model values
}
\description{
  Define the model cost function for kinetic models
}
\author{
  Zhenglei Gao
}

