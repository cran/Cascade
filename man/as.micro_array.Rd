\name{as.micro_array}
\alias{as.micro_array}

\title{
Coerce a matrix into a micro_array object.
}
\description{
Coerce a matrix into a micro_array object.
}
\usage{
as.micro_array(M, time, subject)
}

\arguments{
  \item{M}{
A matrix. Contains the microarray measurements. Should of size N * K, with N the number of genes and K=T*P with T the number of time points, and P the number of individuals. This matrix should be created using cbind(M1,M2,...) with M1 a N*T matrix with the measurements for individual 1, M2 a N*T matrix with the measurements for individual 2.
}
  \item{time}{
A vector. The time points measurements.
}
  \item{subject}{
The number of subjects. 
}
}

\value{A micro_array object.}

\references{
Jung, N., Bertrand, F., Bahram, S., Vallat, L., and Maumy-Bertrand, M. (2014). Cascade: a R-package to study, predict and simulate the diffusion of a signal through a temporal gene network. \emph{Bioinformatics}, btt705. 
	
Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M., Bertrand, F., Meyer, N., ... & Bahram, S. (2013). Reverse-engineering the genetic circuitry of a cancer cell with predicted intervention in chronic lymphocytic leukemia. \emph{Proceedings of the National Academy of Sciences}, 110(2), 459-464.
}

\author{
Nicolas Jung, Frédéric Bertrand , Myriam Maumy-Bertrand.
}




\examples{
  if(require(CascadeData)){
	data(micro_US)
	micro_US<-as.micro_array(micro_US,time=c(60,90,210,390),subject=6)
	}
}
