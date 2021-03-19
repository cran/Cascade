#' Coerce a matrix into a micro_array object.
#' 
#' Coerce a matrix into a micro_array object.
#' 
#' 
#' @param M A matrix. Contains the microarray measurements. Should of size N *
#' K, with N the number of genes and K=T*P with T the number of time points,
#' and P the number of individuals. This matrix should be created using
#' cbind(M1,M2,...) with M1 a N*T matrix with the measurements for individual
#' 1, M2 a N*T matrix with the measurements for individual 2.
#' @param time A vector. The time points measurements.
#' @param subject The number of subjects.
#' @return A micro_array object.
#' @author Nicolas Jung, Frédéric Bertrand , Myriam Maumy-Bertrand.
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
#' @examples
#' 
#'   if(require(CascadeData)){
#' 	data(micro_US)
#' 	micro_US<-as.micro_array(micro_US,time=c(60,90,210,390),subject=6)
#' 	}
#' 
#' @export
as.micro_array<-function(M,time,subject){
  
  if(is.null(row.names(M))){row.names(M)<-paste("gene",1:dim(M)[1])}
  return(new("micro_array",microarray=as.matrix(M),name=row.names(M),time=time,subject=subject,group=0,start_time=0))	
}

cv.lars1 <- function (x, y, K = 10, index, trace = FALSE, plot.it = TRUE, 
          se = TRUE, type = c("lasso", "lar", "forward.stagewise", 
                              "stepwise"), mode = c("fraction", "step"), cv.fun
          #, cv.fun.name
          , ...) 
{
#  requireNamespace("lars")
#  cat(cv.fun.name)
  type = match.arg(type)
  if (missing(mode)) {
    mode = switch(type, lasso = "fraction", lar = "step", 
                  forward.stagewise = "fraction", stepwise = "step")
  }
  else mode = match.arg(mode)
  all.folds <- cv.fun(length(y), K)
#  cat(all.folds[[1]],"\n")
  if (missing(index)) {
    index = seq(from = 0, to = 1, length = 100)
    if (mode == "step") {
      fit = lars::lars(x, y, type = type, ...)
      nsteps = nrow(fit$beta)
      maxfold = max(sapply(all.folds, length))
      nsteps = min(nsteps, length(y) - maxfold)
      index = seq(nsteps)
    }
  }
  residmat <- matrix(0, length(index), K)
  for (i in seq(K)) {
    omit <- all.folds[[i]]
    fit <- lars::lars(x[-omit, , drop = FALSE], y[-omit], trace = trace, 
                type = type, ...)
    fit <- lars::predict.lars(fit, x[omit, , drop = FALSE], mode = mode, 
                   s = index)$fit
    if (length(omit) == 1) 
      fit <- matrix(fit, nrow = 1)
    residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
    if (trace) 
      cat("\n CV Fold", i, "\n\n")
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object <- list(index = index, cv = cv, cv.error = cv.error, 
                 mode = mode)
  if (plot.it) 
    lars::plotCVLars(object, se = se)
  invisible(object)
}

lasso_reg<-function(M,Y,K,eps,cv.fun=lars::cv.folds
                    #,cv.fun.name="lars::cv.folds"
                    ){
  #  require(lars)
#  cat("lasso_reg",cv.fun.name,"\n")
  model<-try(cv.lars1(t(M),(Y),intercept=FALSE,K=K,plot.it=FALSE,eps=10^-5,cv.fun=cv.fun
                      #, cv.fun.name=cv.fun.name
                      ))
  n<-try(model$index[which(model$cv %in% min(model$cv))])
  model<-try(lars::lars(t(M),(Y),intercept=FALSE,eps=10^-5))
  repu<-try(lars::coef.lars(model,s=n,mode="fraction"))
  if(!is.vector(repu)){repu<-rep(0,dim(M)[1])}
  return(repu)
}

lasso_reg_old<-function(M,Y,K,eps){
  #  require(lars)
  model<-try(lars::cv.lars(t(M),(Y),intercept=FALSE,K=K,plot.it=FALSE,eps=10^-5))
  n<-try(model$index[which(model$cv %in% min(model$cv))])
  model<-try(lars::lars(t(M),(Y),intercept=FALSE,eps=10^-5))
  repu<-try(lars::coef.lars(model,s=n,mode="fraction"))
  if(!is.vector(repu)){repu<-rep(0,dim(M)[1])}
  return(repu)
}

F_f<-function(F,x){
  for(i in 1:dim(F)[1]){
    F[col(F)==row(F)-(i-1)]<-x[i]
  }
  #for(i in 1:dim(F)[1]){
  #if(sum(abs(F[i,]))!=0){
  #F[i,]<-F[i,]/sum(abs(F[i,]))
  #}
  #}
  #for(i in  1:dim(F)[1]){
  #if(sum(F[i,])==0){F[i,i]<-1}
  #}
  return(F)
}

sumabso<-function(x){
  d<-sum(abs(x))
  if(d==0){
    d<-1
  }
  return(d)
}

expo<-function(x,hub){	
  sum(x>=hub)/sum(x>0)
}

choice_cutoff<-function(O,nb,eps,hub,plot.g=FALSE,prop.hub=NULL){
  O<-abs(O)
  #plus petit omega sup a eps
  minO<-min(O[O>eps])
  #mediane des omegas
  qq<-quantile(O[O>eps],0.50)
  #max des omega
  maxO<-max(O)
  #cutoffs
  sequence<-seq(minO,maxO,length.out=nb)
  
  Mcut<-array(0,c(dim(O)[1],nb))
  #pour chaque valeur du cutoff on calcule 
  for(i in 1:nb){
    Mcut[,i]<-apply(O>sequence[i],1,sum)
  }
  
  expoh<-function(x){expo(x,hub)}
  hh<-apply(Mcut[,1:(dim(Mcut)[2]-1)],2,expoh)
  if(plot.g==TRUE){
    matplot(t(Mcut),type="l")
    
    dev.new()
    plot(sequence[which(hh>0)],hh[hh>0],type="l",xlab="cut off",ylab="Proportion of hubs")
  }

  auto<-min(which(hh==min(hh[hh>0])))
  if(!is.null(prop.hub)){
    if(prop.hub>min(hh[hh>0])){
      auto<-min(which(hh==min(hh[hh>prop.hub])))
      #abline(h=prop.hub)
    }
    
  }
  
  #print(Mcut[,auto])
  return(sequence[auto])
  
}

choice_cutoff_final<-function(O,nb,eps,hub,plot.g=FALSE,prop.hub){
  
  if(length(hub)>1){
    choix<-NULL
    for(i in hub){
      choix<-c(choix,choice_cutoff(O,nb,eps,i,prop.hub=prop.hub)	)
    }
    plot(hub,choix,type="l")
    lines(hub,predict(loess(choix ~ hub, span=0.75)),col="red")
    return(choix)
  }
  else{
    choix<-NULL
    for(i in prop.hub){
      choix<-c(choix,choice_cutoff(O,nb,eps,hub,prop.hub=i)	)
    }
    plot(prop.hub,choix,type="l")
    lines(prop.hub,predict(loess(choix ~ prop.hub, span=0.75)),col="red")
    return(choix)
  }
}

#simulations



#' Generates a network.
#' 
#' Generates a network.
#' 
#' 
#' @param nb Integer. The number of genes.
#' @param time_label Vector. The time points measurements.
#' @param exp The exponential parameter, as in the barabasi.game function in
#' igraph package.
#' @param init The attractiveness of the vertices with no adjacent edges. See
#' barabasi.game function.
#' @param regul A vector mapping each gene with its number of regulators.
#' @param min_expr Minimum of strength of a non-zero link
#' @param max_expr Maximum of strength of a non-zero link
#' @param casc.level ...
#' @return A network object.
#' @author Nicolas Jung, Frédéric Bertrand , Myriam Maumy-Bertrand.
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
#' @examples
#' 
#' set.seed(1)
#' Net<-network_random(
#' 	nb=100,
#' 	time_label=rep(1:4,each=25),
#' 	exp=1,
#' 	init=1,
#' 	regul=round(rexp(100,1))+1,
#' 	min_expr=0.1,
#' 	max_expr=2,
#' 	casc.level=0.4
#' 	)
#' plot(Net)
#' 
#' @export
network_random<-function(nb,time_label,exp,init,regul,min_expr,max_expr,casc.level){
  
  net<-matrix(0,nb,nb)
  net2<-net
  while(!identical(regul[which(time_label != 1)] , apply(net,2,sum)[which(time_label != 1)])) {
    for(i in 1:nb){
      if(time_label[i] !=1){
        
        if(rbinom(1,1,1-casc.level)==1){	
          reg<-which(time_label<time_label[i])}
        else{
          reg<-which(time_label==(time_label[i]-1))
        }
        if(length(reg)!=0 && sum(sum(net[,i])<regul[i])){
          pb<-apply(net[reg,],1,sum)^exp
          pb<-(pb+init)/(sum(pb+init))
          r<-rmultinom(1, 1, pb)
          net[reg[which(r==1)],i]<-1
          net2[reg[which(r==1)],i]<-runif(1,min_expr,max_expr)*(-1)^rbinom(1,1,0.5)
        }
      }
    }
  }
  length(unique(time_label))->T
  F<-array(0,c(T-1,T-1,T*(T-1)/2))
  for(i in 1:(T*(T-1)/2)){diag(F[,,i])<-1}
  N<-new("network",network=net2,name=paste("gene",1:nb),F=F,convO=0,convF=matrix(0,1,1),time_pt=1:length(unique(time_label)))
  return(N)
}
