\name{analyze_network}
\alias{analyze_network}
\alias{analyze_network-methods}
\alias{analyze_network,network-method}

\title{Analysing the network}
\description{
Calculates some indicators for each node in the network.
}
\usage{
	analyze_network(Omega,nv,...)
}

\arguments{
 \item{Omega}{a network object}
 \item{nv}{the level of cutoff at which the analysis should be done}	
 \item{\dots}{label_v : (optionnal) name of the genes}
}

\value{
  A matrix containing, for each node, its betweenness,its degree, its output, its closeness.
}

\references{
Jung, N., Bertrand, F., Bahram, S., Vallat, L., and Maumy-Bertrand, M. (2014). Cascade: a R-package to study, predict and simulate the diffusion of a signal through a temporal gene network. \emph{Bioinformatics}, btt705. 
	
Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M., Bertrand, F., Meyer, N., ... & Bahram, S. (2013). Reverse-engineering the genetic circuitry of a cancer cell with predicted intervention in chronic lymphocytic leukemia. \emph{Proceedings of the National Academy of Sciences}, 110(2), 459-464.
}

\author{
Nicolas Jung, Frédéric Bertrand , Myriam Maumy-Bertrand.
}


\examples{
	data(network)
	analyze_network(network,nv=0)
}

\keyword{methods}




