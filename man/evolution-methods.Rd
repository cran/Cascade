\name{evolution}
\alias{evolution}
\alias{evolution-methods}
\alias{evolution,network-method}

\title{See the evolution of the network with change of cutoff}
\description{
See the evolution of the network with change of cutoff. 
This function may be usefull to see if the global topology is changed while increasing the cutoff. 
}
\usage{
	evolution(net,list_nv,\dots)
}

\arguments{
 \item{net}{a network object}
 \item{list_nv}{a vector of cutoff at which the network should be shown}
 \item{\dots}{Optionnal arguments:
 	\describe{ \item{gr}{a vector giving the group of each gene}
  \item{color.vertex}{a vector giving the color of each node}
  \item{fix}{logical, should the position of the node in the network be calculated once at the beginning ? Defaut to TRUE.}
  \item{taille}{vector giving the size of the plot. Default to c(2000,1000)}
 \item{\dots}{see plot function}
}
}}

\value{ A HTML page with the evolution of the network.}


\references{
Jung, N., Bertrand, F., Bahram, S., Vallat, L., and Maumy-Bertrand, M. (2014). Cascade: a R-package to study, predict and simulate the diffusion of a signal through a temporal gene network. \emph{Bioinformatics}, btt705. 
	
Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M., Bertrand, F., Meyer, N., ... & Bahram, S. (2013). Reverse-engineering the genetic circuitry of a cancer cell with predicted intervention in chronic lymphocytic leukemia. \emph{Proceedings of the National Academy of Sciences}, 110(2), 459-464.
}

\author{
Nicolas Jung, Frédéric Bertrand , Myriam Maumy-Bertrand.
}



\examples{
\donttest{
	data(network)
	sequence<-seq(0,0.2,length.out=20)
	#setwd("inst/animation")
	#evolution(network,sequence)
}
}

\keyword{methods}

