#' @keywords internal
#' @aliases Cascade-package Cascade NULL
#' @author This package has been written by Frédéric Bertrand, Myriam
#' Maumy-Bertrand and Nicolas Jung with biological insights from Laurent
#' Vallat. Maintainer: Frédéric Bertrand <frederic.bertrand@@utt.fr>
#' 
#' @references Jung, N., Bertrand, F., Bahram, S., Vallat, L., and
#' Maumy-Bertrand, M. (2014). Cascade: a R-package to study, predict and
#' simulate the diffusion of a signal through a temporal gene network.
#' \emph{Bioinformatics}, btt705.
#' 
#' Vallat, L., Kemper, C. A., Jung, N., Maumy-Bertrand, M., Bertrand, F.,
#' Meyer, N., ... & Bahram, S. (2013). Reverse-engineering the genetic
#' circuitry of a cancer cell with predicted intervention in chronic
#' lymphocytic leukemia. \emph{Proceedings of the National Academy of
#' Sciences}, 110(2), 459-464.
"_PACKAGE"

#' @importFrom VGAM rlaplace rpareto zeta 
#' @importFrom grDevices col2rgb colorRamp dev.new grey rainbow rgb
#' @importFrom graphics abline hist legend lines matplot par rect text
#' @importFrom methods new
#' @importFrom stats aggregate cor lm loess model.matrix quantile rbinom reshape rmultinom rnorm runif var wilcox.test
#' @importFrom utils sessionInfo 
#' @importFrom cluster agnes
#' @importFrom animation ani.options saveHTML
#' @importFrom abind abind
#' @import grid
#' @import igraph 
#' @import lattice
#' @import limma
#' @import magic
#' @import nnls
#' @import splines 
#' @import stats4
#' @import survival
#' @import tnet
#' 
NULL

setGeneric("geneSelection",package="Cascade",def = function(x,y,tot.number,... ){standardGeneric("geneSelection")})
setGeneric("genePeakSelection",package="Cascade",def = function(x,peak,... ){standardGeneric("genePeakSelection")})
setGeneric("unionMicro",package="Cascade",def = function(M1,M2 ){standardGeneric("unionMicro")})
setGeneric("position",package="Cascade",def = function(net,... ){standardGeneric("position")})
setGeneric("geneNeighborhood",package="Cascade",def = function(net,targets,... ){standardGeneric("geneNeighborhood")})
setGeneric("evolution",package="Cascade",def = function(net,list_nv,... ){standardGeneric("evolution")})
setGeneric("inference",package="Cascade",def = function(M,... ){standardGeneric("inference")})
setGeneric("cutoff",package="Cascade",def = function(Omega,... ){standardGeneric("cutoff")})
setGeneric("analyze_network",package="Cascade",def = function(Omega,nv,...){standardGeneric("analyze_network")})
#setGeneric("predict",def = function(object,...){standardGeneric("predict")})
setGeneric("gene_expr_simulation",package="Cascade",def = function(network,...){standardGeneric("gene_expr_simulation")})
setGeneric("compare",package="Cascade",def = function(Net,Net_inf,nv){standardGeneric("compare")})

#' Class \code{"micro_array"}
#' 
#' The \code{"micro_array"} class
#' 
#' 
#' @name micro_array-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("micro_array", ...)}. %% ~~ describe objects here ~~
#' @keywords classes
#' @examples
#' 
#' showClass("micro_array")
#' 
#' @export
setClass(
  Class = "micro_array",
  representation(
    microarray = "matrix",
    name = "vector",
    group = c("vector", NULL),
    start_time = c("vector", NULL),
    time = c("vector", NULL),
    subject = "numeric"
  ),
  prototype = prototype(group = 0, start_time = 0),
  validity = function(object) {
    
    if(dim(object@microarray)[2] != length(object@time)*object@subject){
      
      stop("[Error: ]Number of colomns must be equal to the number of time points * the number of subject")
    }
    if(dim(object@microarray)[1] != length(object@name)&&length(object@name)!=0){
      
      stop("[Error: ] Length of the vector of names must equal to the number of genes")
    }
    
    if(dim(object@microarray)[1] != length(object@group)&&length(object@group)!=1){
      
      print(object@group)
      stop("[Error: ] Length of the vector of group must equal to the number of genes or null")
    }
    
    if(dim(object@microarray)[1] != length(object@start_time)&&length(object@start_time)!=1){
      
      
      stop("[Error: ] Length of the vector of starting time must equal to the number of genes or null")
    }
    
    
    if(object@subject<1){
      
      stop("[Error: ] There must be at least one subject")
    }	
    
  }
  
)

#' Class \code{"network"}
#' 
#' The \code{"network"} class
#' 
#' 
#' @name network-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("network", ...)}. %% ~~ describe objects here ~~
#' @keywords classes
#' @examples
#' 
#' showClass("network")
#' 
#' @export
setClass(Class = "network",
         representation(network="matrix",name="vector",F="array",convF="matrix",convO="vector",time_pt="vector")
)


#' Class \code{"micropredict"}
#' 
#' The \code{"micropredict"} class
#' 
#' 
#' @name micropredict-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("micropredict", ...)}.
#' @keywords classes
#' @examples
#' 
#' showClass("micropredict")
#' 
#' @export
setClass(Class = "micropredict",
         representation(microarray_unchanged="micro_array"
                        ,microarray_changed="micro_array"
                        ,microarray_predict="micro_array"
                        ,nv="numeric"
                        ,network="network"
                        ,targets="numeric")
)

