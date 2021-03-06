% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network.R
\name{geneNeighborhood,network-method}
\alias{geneNeighborhood,network-method}
\alias{geneNeighborhood}
\alias{geneNeighborhood-methods}
\title{Find the neighborhood of a set of nodes.}
\usage{
\S4method{geneNeighborhood}{network}(
  net,
  targets,
  nv = 0,
  order = length(net@time_pt) - 1,
  label_v = NULL,
  ini = NULL,
  frame.color = "white",
  label.hub = FALSE,
  graph = TRUE,
  names = FALSE
)
}
\arguments{
\item{net}{a network object}

\item{targets}{a vector containing the set of nodes}

\item{nv}{the level of cutoff. Defaut to 0.}

\item{order}{of the neighborhood. Defaut to `length(net@time_pt)-1`.}

\item{label_v}{vector defining the vertex labels.}

\item{ini}{using the ``position'' function, you can
fix the position of the nodes.}

\item{frame.color}{color of the frames.}

\item{label.hub}{logical ; if TRUE only the hubs are labeled.}

\item{graph}{plot graph of the network. Defaults to `TRUE`.}

\item{names}{return names of the neighbors. Defaults to `FALSE`.}
}
\value{
The neighborhood of the targeted genes.
}
\description{
Find the neighborhood of a set of nodes.
}
\examples{

data(Selection)
data(network)
#A nv value can chosen using the cutoff function
nv=.11 
EGR1<-which(match(Selection@name,"EGR1")==1)
P<-position(network,nv=nv)

geneNeighborhood(network,targets=EGR1,nv=nv,ini=P,
label_v=network@name)

}
\references{
Jung, N., Bertrand, F., Bahram, S., Vallat, L., and
Maumy-Bertrand, M. (2014). Cascade: a R-package to study, predict and
simulate the diffusion of a signal through a temporal gene network.
\emph{Bioinformatics}, btt705.

Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M., Bertrand, F.,
Meyer, N., ... & Bahram, S. (2013). Reverse-engineering the genetic
circuitry of a cancer cell with predicted intervention in chronic
lymphocytic leukemia. \emph{Proceedings of the National Academy of
Sciences}, 110(2), 459-464.
}
\author{
Nicolas Jung, Frédéric Bertrand , Myriam Maumy-Bertrand.
}
\keyword{methods}
