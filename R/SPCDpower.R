SPCDPower=function(n=NULL, power=NULL, p,w=0.5, placeboProp=.66, drop = 0, alpha = 0.025,
          effect_size = rep(NULL, 2)){
  if(!is.null(effect_size[1])) p[1,]=find_probs(p[2,],effect_size)
  if(is.null(n)+is.null(power)!=1)
    print("n or power must be specified but not both")
  #descrete
  npred=!is.null(n)
  if(!npred) {den=qnorm(power)+qnorm(1-alpha);n=1}
  C2=p[1,1];C3=p[2,1];C4=p[1,2];C5=p[2,2]#from spreadsheet
  C6=1-drop
  C7=w;C8=placeboProp;
  if(is.finite(n)) C9=n else C9=1
  ncen= (C4 + C5*(-1 + C7) + (C2 - C3)*C7 - C4*C7)/
    sqrt((((-1 + C2)*C2*C7^2)/(-1 + C8) +
            (2*(((-1 + C4)*C4*(-1 + C7)^2)/((-1 + C3)*C6) -
                  ((-1 + C3)*C3*C7^2)/4))/C8 +
            (2*(((-1 + C5)*C5*(-1 + C7)^2)/((-1 + C3)*C6) -
                  ((-1 + C3)*C3*C7^2)/4))/C8)/C9)
  if(npred) pow0=pnorm(ncen-qnorm(1-alpha))
  else n0=as.integer(ceiling(1/(ncen/den)^2))
  a=placeboProp/2
  k=qnorm(1-p[2,1])
  mu=k-qnorm(1-p)
  contr=c(w,-w,(1-w),-(1-w))
  var=1/c(n*(1-2*a),2*n*a,rep((1-drop)*(1-p[2,1])*n*a,2))
  if (is.null(effect_size[1]))
    ncen1=sum(mu*contr)/sqrt(sum(contr*contr*var))
  else ncen1=(w*effect_size[1]+(1-w)*effect_size[2])/sqrt(sum(contr*contr*var))
  #ncen1a=sum(mu*contr)/sqrt(w^2/(2*a*n*(1-2*a))+2*(1-p[2,1])*(1-w)^2/(a*n*(1-p[2,1])^2))
  firstNonCent=(mu[1,1]-mu[2,1])/sqrt(sum(var[1:2]))
  secondNonCent=(mu[1,2]-mu[2,2])/sqrt(sum(var[3:4]))
  pow1a=function(nn,correct)
    1-pnorm(qnorm(1-alpha*correct)-firstNonCent*sqrt(nn))*
    pnorm(qnorm(1-alpha*correct)-secondNonCent*sqrt(nn))
  if(npred) pow1=pnorm(ncen1-qnorm(1-alpha))
  else n1=as.integer(ceiling(1/(ncen1/den)^2))
  ncen2=(p[1,1]-p[2,1])/sqrt(2*p[1,1]*(1-p[1,1])/n+2*p[1,2]*(1-p[1,2])/n)
  if(npred) pow2=pnorm(ncen2-qnorm(1-alpha))
  else n2=as.integer(ceiling(1/(ncen2/den)^2))
  if (is.null(effect_size[1]))
    ncen3=(mu[1,1]-mu[2,1])/sqrt(4/n)
  else ncen3=effect_size[1]/sqrt(4/n)
  if(npred) pow3=pnorm(ncen3-qnorm(1-alpha))
  else n3=as.integer(ceiling(1/(ncen3/den)^2))
  if(npred) {pow5=pow1a(1,1);pow6=pow1a(1,.5)}
  else {
    n5=ceiling(uniroot(function(nn) pow1a(nn,1)-power,interval=c(1,1000))$root)
    n6=ceiling(uniroot(function(nn) pow1a(nn,.5)-power,interval=c(1,1000))$root)
  }
  if(npred) return(c(n=n,SPCD_Response=pow0,
                     SPCD_Continuous=pow1,
                     Conventional_Response=pow2,
                     Conventional_Continous=pow3,
                     uncorrected=pow5,
                     corrected=pow6
  ))
  else return(c(power=power,SPCD_Response=n0,
                SPCD_Continuous=n1,
                Conventional_Response=n2,
                Conventional_Continuous=n3,
                uncorrected=n5,
                corrected=n6)
  )
}

find_probs=function(p,effect_sizes){
  k=qnorm(1-p[1])
  return(c(p11=1-pnorm(k-effect_sizes[1]),p12=1-pnorm(qnorm(1-p[2])-effect_sizes[2])))
}
