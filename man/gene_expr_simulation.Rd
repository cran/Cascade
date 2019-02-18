\name{gene_expr_simulation}
\alias{gene_expr_simulation}
\alias{gene_expr_simulation-methods}
\alias{gene_expr_simulation,network-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates microarray data based on a given network. 
}
\description{
Simulates microarray data based on a given network. 
}
\usage{
gene_expr_simulation(network,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{A network object.}
\item{...}{\describe{
  \item{time_label}{a vector containing the time labels.}
  \item{subject}{the number of subjects}
  \item{level_peak}{the mean level of peaks. }
}}
}

\value{
A micro_array object.
}

\references{
Jung, N., Bertrand, F., Bahram, S., Vallat, L., and Maumy-Bertrand, M. (2014). Cascade: a R-package to study, predict and simulate the diffusion of a signal through a temporal gene network. \emph{Bioinformatics}, btt705. 
	
Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M., Bertrand, F., Meyer, N., ... & Bahram, S. (2013). Reverse-engineering the genetic circuitry of a cancer cell with predicted intervention in chronic lymphocytic leukemia. \emph{Proceedings of the National Academy of Sciences}, 110(2), 459-464.
}

\author{
Nicolas Jung, Frédéric Bertrand , Myriam Maumy-Bertrand.
}




\examples{
data(simul)
set.seed(1)

#We simulate gene expression according to the network Net
Msim<-gene_expr_simulation(
	network=Net,
	time_label=rep(1:4,each=25),
	subject=5,
	level_peak=200)
head(Msim)
}
