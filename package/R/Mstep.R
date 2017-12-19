Mstep = function (Ytemp, h, Zhat, Ntemp, fGC, mu){
  alpha=mu*Zhat
  offset.temp=alpha+log(Ntemp)+log(fGC)
  glm.fit=glm(Ytemp~h,offset=offset.temp,family=poisson)
  
  glm.coefficients=glm.fit$coefficients
  beta0=glm.coefficients[1]
  g=glm.coefficients[2:length(glm.coefficients)]
  
  pi=sum(Zhat)/length(Ytemp)
  return(list(beta0=beta0, g=g, pi=pi, mu=mu))
}