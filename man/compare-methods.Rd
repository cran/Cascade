\name{compare-methods}
\docType{methods}
\alias{compare-methods}
\alias{compare}
\alias{compare,network,network,numeric-method}
\title{
Some basic criteria of comparison between actual and inferred network.
}
\description{
Allows comparison between actual and inferred network.}
\section{Methods}{
\describe{

\item{\code{signature(Net = "network", Net_inf = "network", nv = "numeric")}}{
\describe{
  \item{Net}{
  A network object containing the actual network.
  }
  \item{Net_inf}{
  A network object containing the inferred network.
  }
  \item{nv}{
  A number that indicates at which level of cutoff the comparison should be done.
  }
}
}                        

}
}

\value{
A vector containing : sensibility, predictive positive value, and the F-score
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

#Comparing true and inferred networks
F_score=NULL

#Here are the cutoff level tested
test.seq<-seq(0,max(abs(Net_inf@network*0.9)),length.out=200)
for(u in test.seq){
	F_score<-rbind(F_score,Cascade::compare(Net,Net_inf,u))
}
matplot(test.seq,F_score,type="l",ylab="criterion value",xlab="cutoff level",lwd=2)
}
