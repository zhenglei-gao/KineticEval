\name{kingraph}
\alias{kingraph}
\title{Function to generate graph data for the GUI to make fit plot.}
\usage{
  kingraph(fit, filename = "graphdata.txt", xlab = "Time",
    ylab = "Observed",
    xlim = c(1, 1.05) * range(fit$data$time),
    ylim = c(1, 1.05) * range((summary(fit))$data$observed, na.rm = TRUE),
    legend = TRUE, length.out = 100, ...)
}
\arguments{
  \item{fit}{An object of class 'kingui'.}

  \item{filename}{The file in which the graph data will be
  stored.}

  \item{xlim}{The plot range for x axis.}

  \item{ylim}{The plot range for y axis.}
}
\value{
  \code{NULL}
}
\description{
  A data table is generated and stored in file=filename for
  later usage.
}
\note{
  This function only writes the data table for making the
  graph.
}
\author{
  Zhenglei Gao
}
\seealso{
  \code{\link{kinplot}}
}

