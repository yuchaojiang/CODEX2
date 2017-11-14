
# EM algorithm
# Expectation function
Estep = function(Ytemp,h,beta0,g,mu,pi,Ntemp,fGC){
    Zhat=rep(NA,length(Ytemp))
    for(i in 1:length(Ytemp)){
        lambda1=Ntemp[i]*fGC[i]*exp(beta0+g%*%h[i,]+mu)
        lambda0=Ntemp[i]*fGC[i]*exp(beta0+g%*%h[i,])
        prob1=dpois(Ytemp[i],lambda=lambda1)*pi
        prob0=dpois(Ytemp[i],lambda=lambda0)*(1-pi)
        if((prob1+prob0)==0){
            if(abs(Ytemp[i]-lambda1)>=abs(Ytemp[i]-lambda0)){
                Zhat[i]=0
            } else{
                Zhat[i]=1
            }
        } else{
            Zhat[i]=prob1/(prob1+prob0)
        }
    }
    return(Zhat)
}

# Maximization function
Mstep = function (Ytemp,h,Zhat,Ntemp,fGC,mu,...){
    alpha=mu*Zhat
    offset.temp=alpha+log(Ntemp)+log(fGC)
    glm.fit=glm(Ytemp~h,offset=offset.temp,family=poisson)
    
    glm.coefficients=glm.fit$coefficients
    beta0=glm.coefficients[1]
    g=glm.coefficients[2:length(glm.coefficients)]
    
    pi=sum(Zhat)/length(Ytemp)
    return(list(beta0=beta0,g=g,pi=pi,mu=mu))
}

