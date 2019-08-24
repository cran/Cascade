\name{inference}
\alias{inference}
\alias{inference-methods}
\alias{inference,micro_array-method}
\title{ Reverse-engineer the network}
\description{
Reverse-engineer the network.
}
\usage{
	inference(M,\dots)
}

\arguments{
 \item{M}{a micro_array object.}
  \item{\dots}{Optional arguments: 	  	\describe{ 
  	  		  	\item{tour.max=30}{maximal number of steps.}
  	\item{g=function(x) 1/x }{the new solution is choosen as (the old solution + g(x) * the new solution)/(1+g(x)) where x is the number of steps. }
  	\item{conv=10e-3}{convergence criterion.}
  	\item{cv.subjects=TRUE}{should the cross validation be done removing the subject one by one ?}
  	\item{nb.folds=NULL}{Relevant only if cv.subjects is FALSE. The number of folds in cross validation.}
  	\item{eps=10e-5}{machine zero}	
  	\item{type.inf="iterative"}{"iterative" or "noniterative" : should the algorithm be computed iteratively}	
  	}
   	}
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
\donttest{
#With simulated data
data(M)
infM <- inference(M)
str(infM)

#With selection of genes from GSE39411
data(Selection)
infSel <- inference(Selection)
str(infSel)
}
}

\keyword{methods}

