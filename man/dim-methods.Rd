\name{dim}
\docType{methods}
\alias{dim}
\alias{dim-methods}
\alias{dim,micro_array-method}
\title{Dimension of the data}
\description{
Dimension of the data
}
\section{Methods}{
\describe{

\item{\code{signature(x = "micro_array")}}{
Gives the dimension of the matrix of measurements.
}
}}
\examples{
  if(require(CascadeData)){
	data(micro_US)
	micro_US<-as.micro_array(micro_US,time=c(60,90,210,390),subject=6)
	dim(micro_US)
	}
}
\keyword{methods}

