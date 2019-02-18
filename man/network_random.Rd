\name{network_random}
\alias{network_random}
\title{
Generates a network.
}
\description{
Generates a network.
}
\usage{
network_random(nb, time_label, exp, init, regul, min_expr, max_expr, casc.level)
}
\arguments{
  \item{nb}{
Integer. The number of genes.
}
  \item{time_label}{
Vector. The time points measurements.
}
  \item{exp}{
The exponential parameter, as in the barabasi.game function in igraph package.
}
  \item{init}{
The attractiveness of the vertices with no adjacent edges. See barabasi.game function.
}
  \item{regul}{
A vector mapping each gene with its number of regulators.
}
  \item{min_expr}{
Minimum of strength of a non-zero link
}
  \item{max_expr}{
Maximum of strength of a non-zero link
}
  \item{casc.level}{
...}
}

\value{
A network object. 
}
\references{
Jung, N., Bertrand, F., Bahram, S., Vallat, L., and Maumy-Bertrand, M. (2014). Cascade: a R-package to study, predict and simulate the diffusion of a signal through a temporal gene network. \emph{Bioinformatics}, btt705. 
	
Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M., Bertrand, F., Meyer, N., ... & Bahram, S. (2013). Reverse-engineering the genetic circuitry of a cancer cell with predicted intervention in chronic lymphocytic leukemia. \emph{Proceedings of the National Academy of Sciences}, 110(2), 459-464.
}

\author{
Nicolas Jung, Frédéric Bertrand , Myriam Maumy-Bertrand.
}

\examples{
set.seed(1)
Net<-network_random(
	nb=100,
	time_label=rep(1:4,each=25),
	exp=1,
	init=1,
	regul=round(rexp(100,1))+1,
	min_expr=0.1,
	max_expr=2,
	casc.level=0.4
	)
plot(Net)
}

