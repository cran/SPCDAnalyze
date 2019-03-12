createTestData = function(times = c(0:11),
                          transition = 5,
                          rx1 = 0.5,
                          rx2 = 1,
                          slope = 0.5,
                          error = 0.3,
                          n = c(50, 50, 50)) {
  nn = sum(n)
  df1 = data.frame(
    intercepts = rnorm(nn, 0, 1),
    slopes1 = rnorm(nn,slope, 1),
    slopes2 = rnorm(nn,slope, 1),
    rx = c(rep(1, n[1]), rep(2, n[2]), rep(3, n[3])),
    ID = 1:nn,
    transition = rep(transition, nn)
  )
  data1 = plyr::ddply(df1, "ID", function(xdf1) {
    pred1 = times <= transition
    if(xdf1$rx==3) mult=rx1 else mult=0
    y1 = xdf1$intercepts + (xdf1$slopes1 + mult)*times[pred1]
    if(xdf1$rx!=1) mult=rx2 else mult=0
    y2 = y1[length(y1)] + (xdf1$slopes2 + mult)*(times[!pred1] - transition)
    vs=data.frame(
      time = times,
      y = c(y1, y2) + rnorm(length(times), 0, error)
    )
    return(vs)
  },.inform=TRUE)
  combData = merge(df1[, c(5, 4, 6)],data1)
  return(combData)
}
#ADDS data as to whether patient is a placebo non-responder to the perpatient file

placeboNonResponder = function(data, ID = 'ID', k = 0) {
  vs = plyr::ddply(data, c('ID'), function(ts) {
    if ((ts[ts$time == ts$transition, 'y'] - ts[1, 'y']) > k)
      outp = TRUE
    else
      outp = FALSE
    return(c(nonResponder = outp))
  })
  return(merge(data, vs))
}


runSlopes = function(combData,times='time') {
  #set start to 0
  combData$newTime = plyr::ddply(combData,"ID", function(data)
    data.frame(V1=data[, times] - min(data[, times])))$V1
  #variable 'trt' supplied by calling program, 0=Placebo 1-Active
  combData$rxx=ifelse(combData$trt=='Placebo',0,1)
  mod1 = lme4::lmer(y ~ newTime + newTime:rxx + (1 + newTime |
                                             ID), data = combData)
  theta1 = summary(mod1)$coefficients
  #returns treatment effect, standard error and t-score
  return(theta1[dim(theta1)[1], ])
}

runGls=function(combData,times='time'){
  #Recoding times by their interger values.
  mx=max(combData[,times])+.1
  combData$newTimes=cut(combData[,times],breaks=unique(c(combData[,times]-.1,mx)),labels=FALSE)
   #sets time 0 to placebo
    combData$trt1=ifelse((combData$newTimes==1)|(combData$trt=='Placebo'),0,1)
    mod1=summary(nlme::gls(y~-1+as.factor(newTimes)+trt1,
                           correlation=nlme::corSymm(form=~newTimes|ID),
                           weights=nlme::varIdent(form=~1|newTimes),
         na.action=na.omit,data=combData,control = list(singular.ok = TRUE)))
    m=length(mod1$coefficients)
#contrasts=summary(multcomp::glht(mod1,linfct=matrix(c(0,0,0,-1,0,0,0,1),1,8),alternative='two.sided',
#                                 test=adjusted(type="none")))
outp=c(mod1$coefficients[m],
        sqrt(mod1$varBeta[m,m]),mod1$coefficients[m]/sqrt(mod1$varBeta[m,m]))
names(outp)<-c("Estimate","Std. Error","t value")
return(outp)
}

SPCDcontinuous = function(combData,recordID,
                          times,
                          group,
                          transition,
                          nonResponder,
                          outcome,
                          runmod = runSlopes,w=0.5)
{
  #Analyze Phase 1
  #isoloate first phase
  #Code to deal with issue where column numbers are used
  namesCombData=names(combData)
  args=as.list(match.call())
  nargs=names(args)
for (i in 3:8) {if (is.numeric(args[[i]])) eval(parse(text=paste(nargs[[i]],"='",namesCombData[args[[i]]],"'",sep='')))}
  #fix problem that y is assumed to be the outcome name.
  combData$y=combData[,outcome]
  combData$ID=combData[,recordID]
  combData=combData[,c(recordID,times,group,transition,nonResponder,outcome)]
  names(combData)[c(1,2,5,6)]<-c("ID",'time','nonResponder','y')
  data = combData[combData[, 'time'] <= combData[, transition], ]
  data$trt = ifelse(data[, group] == 3, 'Active', 'Placebo')
  phase1 = runmod(data)
  #isoloate first phase
  nu1=length(unique(data[,"ID"]))-2
  data = combData[(combData[, 'time'] >= combData[, transition]) &
                    (combData[, 'nonResponder'] == TRUE) &
                    (combData[, group] != 3), ]#set start to zero
  data$trt = ifelse(data[,group] == 2, 'Active', 'Placebo')
  phase2 = runmod(data)
  nu2=length(unique(data[,"ID"]))-2
  outp = c(phase1,
           phase2,
           w*phase1[1] +(1-w)*phase2[1],
           sqrt((w*phase1[2]) ^ 2 + ((1-w)*phase2[2]) ^ 2))
  outp = c(outp, outp[7] / outp[8])
  nup=outp[8]^4/sum( (c(w,1-w)*outp[c(2,5)])^4/c(nu1,nu2))
  p = c(
    pf(outp[3] ^ 2, 1, nu1, lower.tail = FALSE),
    pf(outp[6] ^ 2, 1, nu1, lower.tail = FALSE),
    pf(outp[9] ^ 2, 1, nup, lower.tail = FALSE)
  )
  outp = c(outp, p)
  names(outp) <- c(
    'Phase1 Estimate',
    'Phase1 SE',
    't-value',
    'Phase2 Estimate',
    'Phase2 SE',
    't-value',
    'Pooled Estimate',
    'Pooled SE',
    't-value',
    'Phase1 p',
    'Phase2 p',
    'Pooled p'
  )
  return(outp)
}
createBinaryData=function(p1,q1,p2,q2,s,a,n){
  n1=n2=round(a*n)
  n3=n-n1-n2
  np1=rbinom(1,n1,q1)
  np2=rbinom(1,n2,q1)
  np3=rbinom(1,n3,p1)
  n12=rbinom(1,(n1-np1),s)
  n22=rbinom(1,(n2-np2),s)
  m12=rbinom(1,n12,q2)
  m22=rbinom(1,n22,p2)
  return(matrix(
    c(np1,n1-np1,m12,n12-m12,np2,n2-np2,m22,n22-m22,np3,n3-np3,NA,NA),
    3,4,byrow=TRUE))
}


SPCDbinary=function(results,w=.5){
  f=function(p,n) p*(1-p)/n
  rand1=rowSums(results[,1:2])
  rand2=rowSums(results[,3:4])
  p1=results[3,1]/rand1[3]
  q1=sum(results[1:2,1])/sum(rand1[-3])
  p2=results[2,3]/rand2[2]
  q2=results[1,3]/rand2[1]
  h1=p1-q1
  h2=p2-q2
  h=w*h1+(1-w)*h2
  seh1=sqrt(f(p1,rand1[3])+f(q1,sum(rand1[-3])))
  seh2=sqrt(f(p2,rand2[2])+f(q2,rand2[1]))
  seh=sqrt( f(p1,rand1[3])*w^2 +
            f(q1,sum(rand1[1:2]))*w^2/2 +
            f(q2,rand2[1])*(1-w)^2+
            f(p2,rand2[2])*(1-w)^2)
  outp=c(h1,seh1,h1/seh1,h2,seh2,h2/seh2,h,seh,h/seh)
  p = c(
    pchisq(outp[3] ^ 2, 1, lower.tail = FALSE),
    pchisq(outp[6] ^ 2, 1, lower.tail = FALSE),
    pchisq(outp[9] ^ 2, 1, lower.tail = FALSE)
  )
  outp=c(outp,p)
  names(outp)<-c(
    'Phase1 Estimate',
    'Phase1 SE',
    't-value',
    'Phase2 Estimate',
    'Phase2 SE',
    't-value',
    'Pooled Estimate',
    'Pooled SE',
    't-value',
    'Phase1 p',
    'Phase2 p',
    'Pooled p'
  )
    return(outp)
  }












