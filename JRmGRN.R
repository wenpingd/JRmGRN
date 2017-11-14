#simulate ER based net, the meaning of the parameters are as following.
#p:  the number of genes.
#K:  the number of tissues / conditions to reconstruct GRN jointly, for this version, we consider K=2 only.
#sparsity:  the sparsity of elementary network.
#hubnumber:  the number of hubs.
#hubsparsity:  the sparsity of the hub columns.
#differencerate: the difference rate between two elementary networks.

ER_based_net = function(p,K,sparsity,hubnumber,hubsparsity,differencerate){
    if(sparsity<0 || sparsity>1) stop("sparsity has to take value between 0 and 1!")
    if(hubsparsity<0 || hubsparsity>1) stop("hubsparsity has to take value between 0 and 1!")
    if(hubnumber>p) stop("hubnumber cannot be larger than the number of features!")
    if(hubnumber<0) stop("hubnumber cannot be negative!")
    if(differencerate <0 ||differencerate>0.5) stop("defferencerate has to take value between 0 and 0.5")
    
    Theta=list()
    sparse.base<- rbinom(p*p,1,1-sparsity)*sample(c(-1,1),p*p,replace=TRUE)*runif(p*p,0.25,0.75)
    hubcol <- sample(1:p,hubnumber,replace=FALSE)
    hubvalue=rbinom(hubnumber*p,1,1-hubsparsity)*sample(c(-1,1),hubnumber*p,replace=TRUE)*runif(hubnumber*p,0.25,0.75)
    hubmat=matrix(0,p,p)
    hubmat[,hubcol]=hubvalue
    hubmat=(hubmat+t(hubmat))/2
    for(i in 1:K){
      newsparse=sparse.base
      differpos=sample(1:(p*p), floor(p*p*differencerate),replace=FALSE)
      newsparse[differpos]=rbinom(length(differpos),1,1-sparsity)*sample(c(-1,1),length(differpos),replace=TRUE)*runif(length(differpos),0.25,0.75)
      Theta[[i]]=matrix(data=newsparse,p,p)
      Theta[[i]][lower.tri(Theta[[i]],diag=FALSE)] <- 0
      Theta[[i]] <- Theta[[i]]+t(Theta[[i]])
      Theta[[i]]=Theta[[i]]+hubmat
      Theta[[i]]=ifelse(abs(Theta[[i]])<1e-5,0,Theta[[i]])
      diag(Theta[[i]]) <- 0
      ee <- min(eigen(Theta[[i]],only.values=T)$values)
      diag(Theta[[i]]) <- ifelse(ee < 0, -ee + 0.1, 0.1)
    }
    return(list(Theta=Theta,hubcol=sort(hubcol))) 
  }


# Main function
# S: list of sample covariance matrix
# ns: list of the sample size of each condition
#lambda1,lambda2,lambda3,lambda4: four tuning parameters

JRmGRN=function(S,ns,lambda1,lambda2,lambda3,lambda4,convergence=1e-10,maxiter=2000,
               start="cold",startValue=NULL,verbose=FALSE,penalizeD=TRUE){
  
  K=length(S)
  p=nrow(S[[1]])
  #Initialize the parameters
  # rho is the parameter for scaled Lagrangian form, we set rho=2.5 as is used in HGLASSO.
  rho=2.5
  #primal variable
  if(start=="cold"){
    Theta=list()
    Theta.tilt=list()
    Z=list()
    Z.tilt=list()
    for(k in 1:K){
      Theta[[k]]=diag(p)
      Theta.tilt[[k]]=matrix(0,p,p)
      Z[[k]]=diag(p)
      Z.tilt[[k]]=matrix(0,p,p)
    }
    V=diag(p)
    V.tilt=matrix(0,p,p)
    
    #dual variable
    W.theta=list()
    W.Z=list()
    for(k in 1:K){
      W.theta[[k]]=matrix(0,p,p)
      W.Z[[k]]=matrix(0,p,p)
    }
    W.V=matrix(0,p,p)
  }else{
    Theta=startValue$Theta
    Z=startValue$Z
    V=startValue$V
    Theta.tilt=startValue$Theta.tilt
    Z.tilt=startValue$Z.tilt
    V.tilt=startValue$V.tilt
    W.theta=startValue$W.theta
    W.Z=startValue$W.Z
    W.V=startValue$W.V
  }
  
  criteria = 1e10 	
  i <- 1  	
  # While loop for the iterations
  while(criteria > convergence && i <= maxiter){
    if(verbose){
      print(paste("interation ",i))
    }
    oldTheta=Theta
    
    #update theta(k)
    for(k in 1:K){
      Theta[[k]]=update.theta(S[[k]], Theta.tilt[[k]], W.theta[[k]],ns[[k]],rho)
    }
    
    #update Z
    Z=update.Z(Z.tilt,W.Z,lambda1,lambda2,rho,penalizeD)
    
    #update V
    V=update.V(V.tilt,W.V,lambda3,lambda4,rho)
    
    #update Theta.tilt, Z.tilt, V
    res=update.tilt(Theta,W.theta,Z,W.Z,V,W.V,rho)
    Theta.tilt=res[[1]]
    Z.tilt=res[[2]]
    V.tilt=res[[3]]
    
    #udate W.theta
    for(k in 1:K){
      W.theta[[k]]=W.theta[[k]]+Theta[[k]]-Theta.tilt[[k]]
    }
    
    #update W.Z
    for(k in 1:K){
      W.Z[[k]]=W.Z[[k]]+Z[[k]]-Z.tilt[[k]]
    }
    
    #update W.V
    W.V=W.V+V-V.tilt
    
    criteria <- max(sum((Theta[[1]]-oldTheta[[1]])^2)/sum((oldTheta[[1]])^2),sum((Theta[[2]]-oldTheta[[2]])^2)/sum((oldTheta[[2]])^2))
    criteria
    i=i+1
  }
  V=ifelse(abs(V)<1e-5, 0 , V)
  for(k in 1:K){
    Z[[k]]=ifelse(abs(Z[[k]])<1e-5, 0, Z[[k]])
  }
  if(i>maxiter){
    warning(paste("the real exit converge criteria is ",criteria, " instead of ",convergence ))
  }
  obj=JGLwH.obj(Theta,S,Z,V,ns,lambda1,lambda2,lambda3,lambda4)
  a=abs(V)
  b=apply(a,2,sum)
  hubs=sort(which(b!=0))
  return(list(Theta=Theta,Z=Z,V=V,hubs=hubs,ns=ns,Theta.tilt=Theta.tilt, Z.tilt=Z.tilt, 
              V.tilt=V.tilt, W.theta=W.theta, W.Z=W.Z, W.V=W.V,obj=obj,criteria=criteria))
}
update.theta=function(S,Theta.tilt,W.theta,n,rho){
  
  A=Theta.tilt-W.theta-n/rho*S
  a=eigen(A)
  d=a$values
  D=diag((d+sqrt(d^2+4*n/rho))/2)
  U=a$vectors
  Theta <- U%*%D%*%t(U)
  return (Theta)
}

update.Z=function(Z.tilt,W.Z,lambda1,lambda2,rho,penalizeD){
  K=length(Z.tilt)
  p=nrow(Z.tilt[[1]])
  lam1=penalty.as.matrix(lambda1,p,penalize.diagonal=penalizeD)
  if(K==2){
    A1=Z.tilt[[1]]-W.Z[[1]]
    A2=Z.tilt[[2]]-W.Z[[2]]
    Z1=ifelse(A1>A2+2*lambda2/rho, A1-lambda2/rho, ifelse(A2>A1+2*lambda2/rho,A1+lambda2/rho, (A1+A2)/2))
    Z2=ifelse(A1>A2+2*lambda2/rho, A2+lambda2/rho, ifelse(A2>A1+2*lambda2/rho,A2-lambda2/rho, (A1+A2)/2))
    Z1=sign(Z1)*pmax(abs(Z1)-lam1/rho,0)
    Z2=sign(Z2)*pmax(abs(Z2)-lam1/rho,0)
    return(list(Z1,Z2))
  }
}

update.V=function(V.tilt,W.V,lambda3,lambda4,rho){
  p=nrow(V.tilt)
  V=matrix(0,p,p)
  for(j in 1:p){
    Cj=V.tilt[,j]-W.V[,j]
    Scj=sign(Cj)*pmax(abs(Cj)-lambda3/rho,0)
    a=max(1-lambda4/(rho*sqrt(sum(Scj^2))),0)
    V[,j]=a*Scj
  }
  return(V)
}

update.tilt=function(Theta, W.theta, Z, W.Z, V, W.V, rho){
  K=length(Theta)
  p=nrow(Theta[[1]])
  Gam.s=matrix(0,p,p)
  for(k in 1:K){
    Gam.s=Gam.s+(Theta[[k]]+W.theta[[k]])-(Z[[k]]+W.Z[[k]])-(V+W.V)-t(V+W.V)
  }
  Gam.s=Gam.s*rho/(4*K+2)
  Gam=list()
  Theta.tilt=list()
  Z.tilt=list()
  for(k in 1:K){
    Gam[[k]]=rho/2*((Theta[[k]]+W.theta[[k]])-(Z[[k]]+W.Z[[k]])-(V+W.V)-t(V+W.V))-2*Gam.s
    Theta.tilt[[k]]=(Theta[[k]]+W.theta[[k]])-Gam[[k]]/rho
    Z.tilt[[k]]=(Z[[k]]+W.Z[[k]])+Gam[[k]]/rho
  }
  V.tilt=(V+W.V)+2*Gam.s/rho
  return(list(Theta.tilt,Z.tilt,V.tilt))
}

JGLwH.obj <- function(Theta,S,Z,V,ns,lambda1,lambda2,lambda3,lambda4){
  p = dim(S[[1]])[1]
  K = length(S)
  lam1 = penalty.as.matrix(lambda1,p,penalize.diagonal=FALSE)
  lam2 = penalty.as.matrix(lambda2,p,penalize.diagonal=TRUE)	
  lam3 = penalty.as.matrix(lambda3,p,penalize.diagonal=TRUE)	
  crit = 0
  for(k in 1:K)
  {
    # add log det that was entered as an argument, or else calculate it
    crit = crit-ns[[k]]*log(det(Theta[[k]]))+ns[[k]]*sum(S[[k]]*Theta[[k]])+sum(lam1*abs(Z[[k]])) 
    for(kp in k:length(Theta))
    {
      crit = crit + sum(lam2*abs(Z[[k]]-Z[[kp]]))
    }
  }
  crit=crit+sum(lam3*abs(V))
  tempV <- sum(sqrt(apply(V^2,2,sum)))
  crit=crit+lambda4*tempV
  return(crit)
  
}

penalty.as.matrix <-
  function(lambda,p,penalize.diagonal)
  {
    # for matrix penalties:  check dim and symmetry:
    if(is.matrix(lambda))
    {
      if(sum(lambda!= t(lambda))>0) {stop("error: penalty matrix is not symmetric")}
      if(sum(abs(dim(lambda)-p))!=0 ) {stop("error: penalty matrix has wrong dimension")}
    }
    # for scalar penalties: convert to matrix form:
    if(length(lambda)==1) {lambda=matrix(lambda,p,p)}
    # apply the penalize.diagonal argument:
    if(!penalize.diagonal) {diag(lambda)=0}
    return(lambda)
  }  
AIC.JGLwH = function(S,ns, jglwh.res,c=0.5){
  p = dim(S[[1]])[1]
  K = length(S)
  n=mean(unlist(ns))
  aic= -ns[[1]]*log(det(jglwh.res$Theta[[1]]))+ns[[1]]*sum(S[[1]]*jglwh.res$Theta[[1]])
  -ns[[2]]*log(det(jglwh.res$Theta[[2]]))+ns[[2]]*sum(S[[2]]*jglwh.res$Theta[[2]])
  
  a1=jglwh.res$Z[[1]]
  a2=jglwh.res$Z[[2]]  
  aic=aic+2*(sum(((a1!=0)+(a2!=0))!=0)-p)/2
  aic=aic+2*(length(jglwh.res$hubs))
  b=jglwh.res$V+t(jglwh.res$V)
  diag(b)=0
  aic=aic+2*(c*(sum(b!=0)/2-length(jglwh.res$hubs)))
  return (aic)
}
BIC.JGLwH = function(S,ns, jglwh.res,c=0.5){
  p = dim(S[[1]])[1]
  K = length(S)
  n=mean(unlist(ns))
  bic= -ns[[1]]*log(det(jglwh.res$Theta[[1]]))+ns[[1]]*sum(S[[1]]*jglwh.res$Theta[[1]])
  -ns[[2]]*log(det(jglwh.res$Theta[[2]]))+ns[[2]]*sum(S[[2]]*jglwh.res$Theta[[2]])
  
  a1=jglwh.res$Z[[1]]
  a2=jglwh.res$Z[[2]]  
  bic=bic+log(n)*(sum(((a1!=0)+(a2!=0))!=0)-p)/2
  bic=bic+log(n)*(length(jglwh.res$hubs))
  b=jglwh.res$V+t(jglwh.res$V)
  diag(b)=0
  bic=bic+log(n)*(c*(sum(b!=0)/2-length(jglwh.res$hubs)))
  return (bic)
}

create_grid_mul=function(lam.min,lam.max,n){
  u=(lam.max/lam.min)^(1/(n-1))
  lam.grid=lam.min*(u)^(0:(n-1))
  lam.grid
}

check_sim=function(net,ans,S,ns){
  num_true_hubs=length(net$hubcol)
  num_est_hubs=length(ans$hubs)
  TPR_hubs=length(intersect(ans$hubs,net$hubcol))/length(net$hubcol)
  FPR_hubs=length(setdiff(ans$hubs,net$hubcol))/(nrow(net$Theta[[1]])- length(net$hubcol))
  
  
  
  K=length(net$Theta)
  p=nrow(net$Theta[[1]])
  true_edges=0
  false_edges=0
  est_true_edges=0
  est_false_edges=0
  for(i in 1:K){
    true_edges=true_edges+(sum(net$Theta[[i]]!=0)-p)/2
    false_edges=false_edges+p*(p-1)/2-(sum(net$Theta[[i]]!=0)-p)/2
    
    est_edge=(sum((ans$Z[[i]]+ans$V+t(ans$V))!=0)-p)/2
    est_t_edge=(sum(((ans$Z[[i]]+ans$V+t(ans$V))!=0)*(net$Theta[[i]]!=0))-p)/2
    est_true_edges=est_true_edges+est_t_edge
    est_false_edges=est_false_edges+(est_edge-est_t_edge)
  }
  TPR_edges=est_true_edges/true_edges
  FPR_edges=est_false_edges/false_edges
  
  c(num_true_hubs,num_est_hubs,TPR_hubs,FPR_hubs,TPR_edges,FPR_edges)
}

