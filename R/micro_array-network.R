#' Prediction of the gene expressions after a knock-out experience
#' \code{predict}
#' 
#' Prediction of the gene expressions after a knock-out experience
#' 
#' 
#' @aliases predict predict-methods predict,ANY-method
#' predict,micro_array-method
#' @param object a micro_array object
#' @param Omega a network object.
#' @param nv [=0] numeric; the level of the cutoff
#' @param targets [NULL] vector; which genes are knocked out?
#' @param adapt [TRUE] boolean; do not raise an error if used with vectors 
#' instead of one column matrices.
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
#' @keywords methods
#' @examples
#' 
#' data(Selection)
#' data(network)
#' #A nv value can chosen using the cutoff function
#' nv=.11
#' EGR1<-which(match(Selection@name,"EGR1")==1)
#' P<-position(network,nv=nv)
#' 
#' #We predict gene expression modulations within the network if EGR1 is experimentaly knocked-out. 
#' prediction_ko5<-predict(Selection,network,nv=nv,targets=EGR1)
#' 
#' #Then we plot the results. Here for example we see changes at time point t2:
#' plot(prediction_ko5,time=2,ini=P,label_v=Selection@name)
#' 
#' @exportMethod predict
setMethod("predict"
          ,c("micro_array")
          ,function(object
                    ,Omega
                    ,nv=0
                    ,targets=NULL
                    ,adapt=TRUE
          ){
            require(magic)
            
            micro<-object
            if(!is.null(targets)){
              micro@microarray[targets,]<-0
            }
            #groups
            groupe<-object@group
            #microarray
            M<-micro@microarray
            #number of measurements (=number of timepoints)
            T<-length(object@time)
            #genes
            gene<-1:length(groupe)
            gene2<-gene
            #first times
            supp<-seq(1,T*object@subject,T)
            #all the times
            supp2<-1:(T*object@subject)
            #all times that are not first times
            supp2<-supp2[-supp]
            #removing silenced genes
            if(!is.null(targets)){
              gene<-gene[-targets]
            }
            #weights
            O<-Omega@network
            #cutoff
            O[abs(O)<nv]<-0
            #F matrices
            F<-Omega@F
            #microarray (again)
            microP<-micro
            #predictors
            sup_pred<-rep(1:(T-1),micro@subject)+rep(seq(0,T*(micro@subject-1),T),each=T-1)
            
            
            #microarray (again)
            micro2<-micro
            #index for the F matrices
            u<-0
            for(peak in 2:(T)){
              
              #keeping predictor groups
              IND<-which(groupe[gene2]%in%1:(peak-1))
              grIND<-groupe[IND]
              #silencing targets
              if(!is.null(targets)){
                micro2@microarray[targets,]<-micro@microarray[targets,]
              }
              #predictors genes (timepoints 1..T-1)
              pred<-micro2@microarray[IND,sup_pred]
              #print(pred[IND==72,])
              #for every group
              for(k in 1:(peak-1)){
                #indices for the kth group
                ind<-which(grIND %in% k)
                #increasing the F matrix index
                u<-u+1
                #product of the F matrix and a vector
                f<-function(x){(F[,,u]%*%(x))}#/(sqrt(sum((F[,,u]%*%(x))^2)))}
                #for each subject
                for(i in 1:micro@subject){
                  pred[ind,1:(T-1)+(i-1)*(T-1)]<-t(apply(pred[ind,1:(T-1)+(i-1)*(T-1)],1,f))
                }
              }
              pred[is.na(pred)]<-0
              #Indices des genes du groupe d'interet
              IND2<-which(groupe[gene2]==(peak))
              #Pour chaque gene du groupe reponse
              for(j in IND2){
                #predj<-(pred)*O[IND,j]
                #Pour chaque gene on cherche les genes du groupe predicteur qui ont un lien non nul
                predj<-pred[O[IND,j]!=0,]
                #Si il y en a
                if(length(predj)!=0){
                  #On isole les valeurs finales dans Y (reponse)
                  Y<-micro@microarray[j,sup_pred+1]
                  if(adapt==TRUE){
                    if(!is.null(dim(predj))){
                      mm<-lm(Y~t(predj)-1)
                    }
                    else{
                      mm<-lm(Y~(predj)-1)
                    }
                    #On remplace par la valeur inferee
                    micro2@microarray[j,sup_pred+1]<-predict(mm)
                    #On update les coefficients de Omega
                    O[IND,j][O[IND,j]!=0]<-coef(mm)[]
                  }
                  #Si il n'y en a pas
                }
                else{
                  #On prend la somme ponderee des sorties 
                  #des genes exprimes aux temps precedents
                  predj<-apply((pred)*O[IND,j],2,sum)	
                  micro2@microarray[j,sup_pred+1]<-predj
                }
                
              }
            }
            #Pour pouvoir faire un plot
            micro33<-micro2
            if(!is.null(targets)){ 
              pppp<- unique(unlist(geneNeighborhood(Omega,targets,nv,graph=FALSE)))
              genes3<-gene2[-pppp]                           
            }            else{
              genes3<-gene2
            }
            if(!is.null(targets)){
              micro33@microarray[ genes3,]<-micro@microarray[genes3,]
            }		
            
            if(is.null(targets)){
              targets<- -1
            }
            subjects<-object@subject
            times<-object@time
            ntimes<-length(times)
            patients<-paste(rep("P",subjects*ntimes),rep(1:subjects,each=ntimes),sep="")
            temps<-paste(rep("T",subjects*ntimes),rep(times,subjects),sep="")
            indicateurs<-paste(patients,temps,sep="")
            expr<-rep("log(S/US)",subjects*ntimes)
            nomscol<-paste(expr,":",indicateurs)
            #nomscol<-paste(patients,timestring1)
            
            colnames(object@microarray)<-nomscol
            colnames(micro@microarray)<-nomscol
            colnames(micro33@microarray)<-nomscol
            return(new("micropredict"
                       ,microarray_unchanged=object
                       ,microarray_changed=micro
                       ,microarray_predict=micro33
                       ,nv=nv
                       ,network=Omega
                       ,targets=targets
            )
            )	
          }
)


#' Reverse-engineer the network
#' 
#' Reverse-engineer the network.
#' 
#' 
#' @aliases inference inference-methods inference,micro_array-method
#' @param M a micro_array object.
#' @param tour.max maximal number of steps. Defaults to `tour.max=30`
#' @param g the new solution is choosen as
#' (the old solution + g(x) * the new solution)/(1+g(x)) where x is the number
#' of steps. Defaults to `g=function(x) 1/x` 
#' @param conv convergence criterion. Defaults to `conv=10e-3`
#' @param cv.subjects should the cross validation be done removing the
#' subject one by one ? Defaults to `cv.subjects=TRUE`.
#' @param nb.folds Relevant only if cv.subjects is
#' FALSE. The number of folds in cross validation. Defaults to `NULL`.
#' @param eps machine zero. Defaults to `10e-5`.
#' @param type.inf "iterative" or "noniterative" : should the
#' algorithm be computed iteratively. Defaults to `"iterative"`.
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
#' @keywords methods
#' @examples
#' 
#' \donttest{
#' #With simulated data
#' data(M)
#' infM <- inference(M)
#' str(infM)
#' 
#' #With selection of genes from GSE39411
#' data(Selection)
#' infSel <- inference(Selection)
#' str(infSel)
#' }
#' 
#' @exportMethod inference
setMethod(f="inference"
          ,signature=c("micro_array")
          ,definition=function(M
                               ,tour.max=30
                               ,g=function(x){1/x}
                               ,conv=0.001
                               ,cv.subjects=TRUE
                               ,nb.folds=NULL
                               ,eps=10^-5
                               ,type.inf="iterative"
          ){
            #Package requis
            require(nnls)
            #Quelques indicateurs
            mat<-M@microarray 
            #La matrice contenant les donnees
            gr<-M@group 
            #Le vecteur des groupes
            N<-dim(mat)[1] 
            #Nombre de genes
            T<-length(unique(M@group)) 
            #Nombre de temps de mesure
            P<-M@subject 
            #Nombre de patients
            
            #La condition suivante determine le nombre de folds 
            # pour la cross validation
            if(is.null(nb.folds)){
              K<-T-1
            }
            else{
              K<-nb.folds
            }
            
            nF<-(1+T-1)/2*(T-1) 
            #nb de matrices Fab
            #Simple calcul du nombre de matrices Fab 
            # (selon le modele initial)
            
            #Initialisation des matrices Fab
            F<-array(0,c(T-1,T-1,nF)) 
            #J'ai prefere construire une matrice F 
            # Tri dimensionnelle laquelle contient toutes les matrices Fab
            # F[,,1] correspond a la matrice F12
            # F[,,2] correspond a la matrice F13
            # F[,,3] correspond a la matrice F14
            # ... (les matrices nulles ne sont pas 
            # indicees F11 par exemple)
            
            for(i in 1:nF){
              for(j in 1:(T-1)){
                F[j,j,i]<-1
              }
            } 
            #L'initialisation prend pour toutes les matrices Fab 
            #la matrice identite T-1 * T-1
            
            
            #Initialisation de la matrice Omega
            Omega<-array(0,c(N,N))
            
            #Initialisation des deux indicateurs de convergence
            convF<-rep(mean(F^2),nF)
            convO<-mean(mat^2)
            
            #Support : correspond aux colonnes de mat 
            # qui servent pour la prediction
            #Attention le modele se sert que des temps 1 a 
            # T-1 pour chaque patient pour les predicteurs
            sup_pred<-rep(1:(T-1),P)+rep(seq(0,T*(P-1),T),each=T-1)
            
            #Initialisation du nombre de tours	
            tour<-1
            
            #Si on veut faire la version non iterative
            # il faut deux passages : dans le premier on infere
            # F et dans le second on infere omega
            if(type.inf=="noniterative"){
              tour.max<-2
            }
            
            #L'algorithme commence ici
            while(tour<= tour.max && convO[length(convO)]>conv){ 
              #Condition d'arret : soit nombre de tour max atteint, 
              #soit la convergence de la matrice Omega est suffisante
              
              cat(paste("We are at step : ",tour))
              cat("\n") 
              #Pour montrer que l'algorithme est en train de calculer
              
              OmegaS<-Omega 
              #OmegaS comme sauvegarde ; necessaire pour calculer
              # la convergence
              
              u<-0 
              #u a un role essentiel, puisque qu'il determine 
              # laquelle des matrices Fab est indicee. 
              # Se referer plus haut pour en connaitre l'ordre
              
              for(peak in 2:T){ 
                #Comprendre ici que peak correspond au groupe du 
                # gene REPONSE ; en consequence, nous 
                # commencons a deux.
                
                IND<-which(gr %in% 1:(peak-1)) 
                #Ici nous cherchons les genes possiblement 
                # PREDICTEURS.
                #En effet, le gene reponse etant de groupe peak, 
                # un predicteur ne peut etre que de groupe 1 a (peak-1)
                
                grIND<-gr[IND] 				
                #Nous recuperons le groupe des individus possiblement
                # predicteurs.
                pred<-mat[IND,sup_pred]		
                #Nous creons la matrice des predicteurs 
                
                
                for(k in 1:(peak-1)) { 
                  #Cette boucle sert a transformer les 
                  # predicteurs en fonction des matrices Fab
                  #Le groupe de la variable reponse 
                  # etant peak, les groupes des predicteurs 
                  # sont 1..(peak-1)
                  
                  ind<-which(grIND %in% k) 
                  #On regarde successivement, et dans 
                  # l'ordre, tous les groupes possibles
                  
                  u<-u+1 
                  #u est initialise a 0. La premiere fois, 
                  # u vaut donc 1, peak 2, et k =1. Donc F[,,u]=F12.
                  #La deuxieme fois qu'on arrive ici 
                  # peak est passe a 3, et k vaut de nouveau 1  
                  # F[,,u]=F13
                  #La troisieme fois peak reste a 2 car la boucle 
                  # avec k n'est pas finie, et k vaut 2 donc 
                  # F[,,u]=F23
                  # ... c'est bien l'ordre dans lequel nous 
                  # avons range les F
                  
                  f<-function(x){(F[,,u]%*%(x))} 
                  #On construit une fonction generique
                  
                  for(i in 1:P){
                    pred[ind,1:(T-1)+(i-1)*(T-1)]<-t(apply(pred[ind,1:(T-1)+(i-1)*(T-1)],1,f)) 
                    #La transformation est faite ici.
                  }
                }
                #En sortant de la boucle avec k, pred 
                # contient les genes predicteurs correctement 
                # transformes. 				
                
                pred[is.na(pred)]<-0 
                #Ceci est une securite, au cas ou une matrice F 
                # deviendrait nulle
                
                #On construit ci dessous la matrice des vecteurs 
                # reponse
                Y<-mat[which(gr %in% peak),sup_pred+1]
                Omega[IND, which(gr %in% peak)]<-Omega[IND, which(gr %in% peak)]*0
                
                #Nous allons passer au Lasso
                #retenir<-lars::cv.folds 
                #Ceci permet de changer une fonction interne de lars 
                # qui s'occupe de la validation croisee. 
                
                if(cv.subjects==TRUE){
                  cv.folds1=function(n,folds){
                                      split(1:dim(pred)[2]
                                            ,rep(1:P,each=dim(pred)[2]/P))}
#                  cv.fun.name="cv.folds1"
                  } else {
                    cv.folds1=lars::cv.folds
#                    cv.fun.name="lars::cv.folds"
                  }
#                cat(cv.fun.name)
                fun_lasso<-function(x){lasso_reg(pred,x,K=K,eps,cv.fun=cv.folds1
                                                 #,cv.fun.name=cv.fun.name
                                                 )} 
                
                  # assignInNamespace("cv.folds"
                  #                   ,function(n,folds){
                  #                     split(1:dim(pred)[2]
                  #                           ,rep(1:P,each=dim(pred)[2]/P)
                  #                     )
                  #                   }
                  #                   , ns="lars")
                  Omega[IND, which(gr %in% peak)]<-apply(Y,1,fun_lasso)
#                  assignInNamespace("cv.folds",retenir, ns="lars")
                # }
                # else{
                #   Omega[IND, which(gr %in% peak)]<-apply(Y,1,fun_lasso)
                # }
                
              } 
              #fin de la boucle for avec peak ; 
              # la matrice omega est inferee
              
              co<-apply(Omega,2,sumabso)
              Omega<-t(t(Omega)/co)
              
              if(tour!=1 && type.inf=="iterative"){
                Omega<-(g(tour)*Omega+OmegaS)/(1+g(tour)) 
                #On prend seulement une partie de l'innovation
              }
              
              
              convO<-c(convO,mean(abs(Omega-OmegaS)))
              
              if( type.inf=="iterative"){
                cat(paste("The convergence of the network is (L1 norm) :", round(convO[length(convO)],5)))	
                cat("\n")	
              }
              uuu<-0
              sauvF<-F
              
              if(tour==1 && type.inf=="noniterative"){
                Omega<-Omega*0+1 
                #Tous les predicteurs sont mis a un dans le 
                # cadre de l'inference non iterative
              }
              
              #Maintenant, il s'agit d'estimer la matrice F ; 
              # on reprend la meme maniere de faire que pour Omega
              #Nous annotons les differences
              
              for(peak in 2:T){
                IND<-which(gr %in% 1:(peak-1))
                grIND<-gr[IND]
                sup_pred<-rep(1:(T-1),P)+rep(seq(0,T*(P-1),T),each=T-1)
                pred<-(mat[IND,sup_pred])
                IND2<-which(gr %in% peak)
                
                Xf<-NULL
                
                #Important : etant donne qu'ils apparaissent dans 
                # les memes equations Fij, i=1...(j-1) doivent 
                #etre estimees en meme temps
                
                for(i in 1:(peak-1)){ 
                  #Cette  boucle permet de creer la matrice des 
                  # predicteurs selon une forme pratique
                  X<-NULL
                  suma<-function(x){sum(abs(x))}
                  f<-Vectorize(
                    function(x){
                      #Multiplie la matrice des pr??dicteurs du groupe i
                      # par les omega correspondants aux genes du groupe
                      # dont les indices sont dans x et on somme les 
                      # valeurs absolues
                      apply(pred[which(grIND==i),]*Omega[IND[which(grIND==i)],x],2,suma)
                    }
                  )
                  #On applique cette fonction aux genes de peak
                  Xa<-(f(IND2))
                  Xb<-NULL
                  
                  #Boucle sur les patients
                  for(p in 1:P){
                    Q<-NULL
                    for(r in 1:(T-1)){
                      q<-Xa[1:(T-1)+(p-1)*(T-1),]
                      q<-c(rep(0,(r-1)),q[1:(length(q)-(r-1))])
                      if(r!=1){q[which((1:length(q) %% (T-1)) %in% 1:(r-1))]<-0 }
                      Q<-cbind(Q,q)
                    }
                    X<-rbind(X,Q)
                  }
                  Xf<-cbind(Xf,X)
                }
                
                
                Y<-c(t(mat[IND2,sup_pred+1]))
                pond<-rep(0,P)
                coeffi<-array(0,c(P,dim(Xf)[2]))
                
                for(pat in 1:P){
                  support<-1:length(Y)
                  enl<-(1:(length(Y)/(P))+(pat-1)*(length(Y)/(P)))
                  support<-support[-enl]
                  model<-nnls(abs(Xf[support,]),abs(Y[support]))
                  pond[pat]<-1/(mean((Xf[enl,]%*%coef(model)-Y[enl])^2) )
                  coeffi[pat,]<-coef(model)
                }
                
                model<-apply(coeffi*( pond/sum(pond)),2,sum)
                pp<-length(model)/(T-1)
                
                pk<-uuu+1
                
                for(jj in 1:pp){
                  uuu<-uuu+1
                  F[,,uuu]<-F_f(F[,,uuu],model[1:(T-1)+(jj-1)*(T-1)])
                }
              } 
              #fin de la boucle peak
              
              if(type.inf=="iterative"){
                F<-(g(tour)*F+sauvF)/(1+g(tour))
              }
              cc<-rep(0,nF)
              for(i in 1:nF){cc[i]<-mean(abs((F[,,i]/sum(F[,,i])-sauvF[,,i]/sum(sauvF[,,i]))))}
              convF<-cbind(convF,cc)
              tour<-tour+1
            }
            
            if(type.inf=="iterative"){
              plot(convO[-1],type="l")
              #dev.new()
              matplot(t(convF),type="l")
            }
            else{
              F<-sauvF   
            }
            result<-new("network"
                        ,network=Omega
                        ,name=M@name
                        ,F=F
                        ,convF=convF
                        ,convO=convO
                        ,time_pt=M@time
            )
            return(result)
          }
          
)


#' Simulates microarray data based on a given network.
#' 
#' Simulates microarray data based on a given network.
#' 
#' 
#' @aliases gene_expr_simulation gene_expr_simulation-methods
#' gene_expr_simulation,network-method
#' @param network A network object.
#' @param time_label a vector containing the time labels.
#' @param subject the number of subjects
#' @param level_peak the mean level of peaks.
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
#' data(Net)
#' set.seed(1)
#' 
#' #We simulate gene expression according to the network Net
#' Msim<-gene_expr_simulation(
#' 	network=Net,
#' 	time_label=rep(1:4,each=25),
#' 	subject=5,
#' 	level_peak=200)
#' head(Msim)
#' 
#' @exportMethod gene_expr_simulation
setMethod("gene_expr_simulation"
          ,"network"
          ,function(network
                    ,time_label=1:4
                    ,subject=5
                    ,level_peak=100
          ){
            require(VGAM)
            
            N<-network@network
            M<-matrix(0,dim(network@network)[1],length(unique(time_label))*subject)
            T<-length(unique(time_label))
            gene1<-which(time_label==1)
            supp<-seq(1,dim(M)[2],by=length(unique(time_label)))
            M[gene1,supp]<-VGAM::rlaplace(length(supp)*length(gene1),level_peak,level_peak*0.9)*(-1)^rbinom(length(supp)*length(gene1),1,0.5)
            supp<-(1:dim(M)[2])[-supp]
            M[gene1,supp]<-VGAM::rlaplace(length(supp)*length(gene1),0,level_peak*0.3) 
            
            
            for(i in 2:T){
              
              genei<-which(time_label==i)
              supp<-seq(1,dim(M)[2],by=length(unique(time_label)))
              M[genei,supp]<-VGAM::rlaplace(length(supp)*length(genei),0,level_peak*0.3)
              for(j in genei){
                for( t in 2:T){
                  M[j,supp+t-1]<-apply(N[,j]*M[,supp+(t-2)],2,sum) + 	rnorm(length(supp+t),0,50)	
                  
                }			
                
              }
              
              
            }
            
            MM<-as.micro_array(M,1:length(unique(time_label)),subject)
            MM@group<-time_label
            
            
            
            G<-predict(MM,network)@microarray_predict
            supp<-seq(1,dim(M)[2],by=length(unique(time_label)))
            G@microarray[,supp]<-M[, supp]	
            return(G)
          }
)
