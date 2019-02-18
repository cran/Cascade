\name{geneNeighborhood}
\alias{geneNeighborhood}
\alias{geneNeighborhood-methods}
\alias{geneNeighborhood,network-method}
\title{Find the neighborhood of a set of nodes.}

\description{
Find the neighborhood of a set of nodes.
}
\usage{
	geneNeighborhood(net,targets,...)
}

\arguments{
 \item{net}{a network object}
 \item{targets}{a vector containing the set of nodes}
 \item{...}{Optional arguments. See plot options. }
}

\value{
  The neighborhood of the targeted genes. 
}

\references{
Jung, N., Bertrand, F., Bahram, S., Vallat, L., and Maumy-Bertrand, M. (2014). Cascade: a R-package to study, predict and simulate the diffusion of a signal through a temporal gene network. \emph{Bioinformatics}, btt705. 
	
Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M., Bertrand, F., Meyer, N., ... & Bahram, S. (2013). Reverse-engineering the genetic circuitry of a cancer cell with predicted intervention in chronic lymphocytic leukemia. \emph{Proceedings of the National Academy of Sciences}, 110(2), 459-464.
}
\author{
Nicolas Jung, Frédéric Bertrand , Myriam Maumy-Bertrand.
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

\keyword{methods}
