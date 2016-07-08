

stratadisc=function(Treatment,Outcome,Matrix){
  check=as.matrix(unique(Matrix[,Outcome]))
  if(nrow(check) > 2 & sum(check[,1]) != 1)stop("Warning: Output variable must be discrete for this function. For continuous outputs, please use stratacont().")
  a=ncol(Matrix)-1
  for(i in 1:a){
    h=as.matrix(unique(Matrix[,i]))
    if(nrow(h) > 5)message("We advise using input variables with at most 5 categories.")
    if(nrow(h) > 5)dd=sprintf("Column %s contains more than 5 categories" , i)
    if(nrow(h) > 5)message(dd)
  }
  g=stratifydisc(Treatment,Outcome,Matrix)
  l=weightdisc(g)
  l3=summarydisc(g)
  dispdisc(g)
  k=l3
  message("Stratification complete.")
  print(l)
  return(k)

}


stratacont=function(Treatment,Outcome,Matrix){
  check=as.matrix(unique(Matrix[,Outcome]))
  if(nrow(check) <= 2 )stop("Warning: Output variable must be continuous for this function. For discrete outputs, please use stratadisc().")
  a=ncol(Matrix)-1
  for(i in 1:a){
    h=as.matrix(unique(Matrix[,i]))
    if(nrow(h) > 5)message("We advise using input variables with at most 5 categories.")
    if(nrow(h) > 5)dd=sprintf("Column %s contains more than 5 categories" , i)
    if(nrow(h) > 5)message(dd)
  }
  g=stratifycont(Treatment,Outcome,Matrix)
  l=weightcont(g)
  l3=summarycont(g)
  dispcont(g)
  k=l3
  message("Stratification complete.")
  print(l)
  return(k)
}



stratifydisc=function(Treatment,Outcome,Matrix){
  Matrix=as.data.frame(Matrix)
  l=as.data.frame(cbind(Matrix[,Treatment],Matrix[,Outcome]))
  r=c(colnames(Matrix))
  p=c(r[Treatment],r[Outcome])
  colnames(l)=p
  if(Treatment < Outcome){
    Matrix=Matrix[,-Outcome]
    Matrix=Matrix[,-Treatment]
  }
  if(Outcome < Treatment){
    Matrix=Matrix[,-Treatment]
    Matrix=Matrix[,-Outcome]
  }
  Matrix=cbind(Matrix,l)
  for (i in 1:ncol(Matrix)){
  if(colnames(Matrix)[1] == "V1")colnames(Matrix)[1]="Var1"
  }
  z=ncol(Matrix)
  t=ddply(Matrix,(1:z),nrow)
  w=z-2
  km=ncol(t)-3
  mm=ddply(t,(1:km),nrow)
  a=ncol(mm)
  b=a-1
  t=as.matrix(t)
  p=matrix(ncol=b,nrow=nrow(mm))
  pi=matrix(ncol=8,nrow=nrow(mm))
  for(i in 1:ncol(pi)){
    pi[,i]=0
  }
  colnames(p)=colnames(Matrix[,(1:b)])
  mm=as.matrix(mm)
  count=1
  for(SI in 1:nrow(p)){
    p[SI,]=mm[SI,(1:b)]
    qwer=as.numeric(mm[,(b+1)])
    f=qwer[SI]
    temp=matrix(nrow=f,ncol=ncol(t))
    for (k in 1:nrow(t)){
      if(identical(t[k,(1:b)],mm[SI,(1:b)]) == TRUE){
        temp[count,]=t[k,]
        count=count+1

      }
    }
    count=1
    gg=ncol(temp)-2
    gg2=ncol(temp)-1
    gg3=ncol(temp)
    a1=as.numeric(temp[,gg])
    a2=as.numeric(temp[,gg2])
    a3=as.numeric(temp[,gg3])

    kop=cbind(a1,a2,a3)

    for(l in 1:nrow(kop)){
      if(kop[l,1]==1 & kop[l,2]==1) pi[SI,7]=kop[l,3]
      if(kop[l,1]==1 & kop[l,2]==0) pi[SI,8]=kop[l,3]
      if(kop[l,1]==0 & kop[l,2]==1) pi[SI,3]=kop[l,3]
      if(kop[l,1]==0 & kop[l,2]==0) pi[SI,4]=kop[l,3]
    }
    pi[SI,1]=pi[SI,3]+pi[SI,4]
    pi[SI,5]=pi[SI,7]+pi[SI,8]
  }
  for ( i in 1:nrow(pi)){
    if ( pi[i,1] > 0 & pi[i,5] > 0)pi[i,6] = 1
    else pi[i,6] = 0
  }
  for ( i in 1:nrow(pi)){
    if ( pi[i,6] == 1) pi[i,2] = pi[i,5] / pi[i,1]
    else pi[i,2] = 0
  }
  colnames(pi) = c("T=0", "Wk0" , "Y=1" , "Y=0", "T=1" , "Wk1" , "Y=1", "Y=0")
  for ( i in 1:nrow(pi)){
    if(pi[i,2] == 0)pi[i,1]=NA
  }

  p=cbind(p,pi)
  p=na.omit(p)
  return(p)
}





weightdisc = function(mat){
  cola=ncol(mat)
  colb=ncol(mat)-7
  mat=as.matrix(mat[,(colb:cola)])
  sss=sum(mat[,5])
  l=matrix(nrow=nrow(mat),ncol=ncol(mat))
  for(i in 1:ncol(l)){
    l[,i]=as.numeric(mat[,i])
  }
  mat=l
  mat[,3]=mat[,3]*mat[,2]
  mat[,4]=mat[,4]*mat[,2]
  n= 0
  a = 0
  b=0
  c=0
  d=0
  run=0
  run2=0
  E = 0
  V = 0
  O = 0
  v1=0
  v2=0
  v3=0
  v4=0
  for( i in 1:nrow(mat)){
    n = mat[i,3] + mat[i,4]  + mat[i,7] + mat[i,8]
    a = mat[i,7]
    b = mat[i,8]
    c = mat[i,3]
    d = mat[i,4]
    run = run + (a*d) /n
    run2 = run2 + (b*c)/n
    v1 = v1 + (((a*d)/n) * ((a+d)/n))
    v2 = v2 + (((b*c)/n) * ((b+c)/n))
    v3 = v3 + (((a*d)/n) * ((b+c)/n)) + (((b*c)/n) * ((a+d)/n))
    O=O+a
    E = E + ((mat[i,7] + mat[i,8])*(mat[i,7] + mat[i,3]))/ (mat[i,3] + mat[i,4]  + mat[i,7] + mat[i,8])
    V = V + ((mat[i,7] + mat[i,8])*(mat[i,3] + mat[i,8])*(mat[i,7]+mat[i,3])*(mat[i,8]+mat[i,4])) / (mat[i,3] + mat[i,4]  + mat[i,7] + mat[i,8])^2*(mat[i,3] + mat[i,4]  + mat[i,7] + mat[i,8]-1)
  }
  v4 = v1/(2*run^2) + v2/(2*run2^2) + v3/(2*run*run2)
  CHI= (abs(O - E) - 0.5)^2 / V
  l=pchisq(CHI,df=1)
  result = run / run2
  k=matrix(nrow=5,ncol=1)
  rownames(k)=c("Odds Ratio Of Impact Of Treatment on Outcome","Mantel Haenszel P-Value","No. Of Cases Matched" ,"95% C.I Upper Bound" , "95% C.I Lower Bound")
  k[1,1]=round(result,3)
  k[2,1]=round(l,4)
  k[3,1]=sss
  k[4,1]=exp(log(k[1,1]) + pnorm(0.975,0,1)*sqrt(v4))
  k[5,1]=exp(log(k[1,1]) - pnorm(0.975,0,1)*sqrt(v4))
  return(k)
}








oddsdisc=function(mat){

cola=ncol(mat)
colb=ncol(mat)-7
t=as.matrix(mat[,(colb:cola)])
t2=as.matrix(mat[,-(colb:cola)])
l=matrix(nrow=nrow(t),ncol=ncol(t))
for(i in 1:ncol(l)){
  l[,i]=as.numeric(t[,i])
}
t=l

en=matrix(nrow=1,ncol=1)
for( kk in 1:ncol(t2)){
  d=na.omit(as.numeric(unique(t2[,kk])))
  if ( length(d) == 0){
    s=unique(t2[,kk])
    s=as.matrix(s,byrow=T)
    en=rbind(en,s)
  }
  else{
    s=colnames(t2)[kk]
    a=as.matrix(s)
    en=rbind(en,s)
  }
}
en=na.omit(en)
results=matrix(nrow=nrow(en),ncol=2)
results[,2]=rep(1)
colnames(results)=c("Before Weighting" , "After Weighting")
count=1
for ( kk in 1:ncol(t2)){
  d=na.omit(as.numeric(unique(t2[,kk])))
  if ( length(d) == 0 ){
    s=unique(t2[,kk])

    for ( i in 1:length(s)){
      dd=as.matrix(t2[,kk])
      for (o in 1:nrow(dd)){
        if (dd[o,1] != s[i])dd[o,1]=0
        if (dd[o,1] == s[i])dd[o,1]=1
      }
      dd=as.matrix(as.numeric(dd[,1]))
      a=t[,5]*dd[,1]
      b=t[,1]*dd[,1]
      a1=sum(a)/sum(t[,5])
      b1=sum(b)/sum(t[,1])
      f=a1/b1
      results[count,1]=f
      count=count+1
    }
  }
  if (length(d) != 0){
    dd=as.matrix(t2[,kk])
    dd=as.matrix(as.numeric(dd[,1]))
    a=t[,5]*dd[,1]
    b=t[,1]*dd[,1]
    a1=sum(a)/sum(t[,5])
    b1=sum(b)/sum(t[,1])
    f=b1/a1
    results[count,1]=f
    count=count+1
  }

}
rownames(results)=en[,1]
results=na.omit(results)
return(results)
}





dispdisc=function(mat){
  mat=oddsdisc(mat)
time=c(1:nrow(mat))
betagal.abs=c(mat[,1])
cell.density=c(mat[,2])
h=max(mat[,1])
h=h+1
plot(time, betagal.abs, pch=17, axes=F, ylim=c(0,h), xlab="", ylab="",col="black",main="Balance Of Covariates")
axis(2, ylim=c(0,1),col="black")
mtext("Odds",side=2,line=2.5)
box()

par(new=T)

plot(time, cell.density, pch=15,  xlab="", ylab="", ylim=c(0,h), axes=F, type="b", col="red")

bb=nrow(mat)
axis(1,1:bb,rownames(mat),las=2)

legend("topright", legend=c("Before Weighting","After Weighting"), pch=c(17,15))
}


stratifycont=function(Treatment,Outcome,Matrix) {
  Matrix=as.data.frame(Matrix)
  l=as.data.frame(cbind(Matrix[,Treatment],Matrix[,Outcome]))
  r=c(colnames(Matrix))
  p=c(r[Treatment],r[Outcome])
  colnames(l)=p
  if(Treatment < Outcome){
    Matrix=Matrix[,-Outcome]
    Matrix=Matrix[,-Treatment]
  }
  if(Outcome < Treatment){
    Matrix=Matrix[,-Treatment]
    Matrix=Matrix[,-Outcome]
  }
  Matrix=cbind(Matrix,l)
  d=as.data.frame(Matrix)
  for (i in 1:ncol(d)){
  if(colnames(d)[1] == "V1")colnames(d)[1]="Var1"
  }
d=un1(d)
 a=makeit(d)
 gg=un3(a,d)
 if (nrow(gg) == 0)stop("Stratification could not be completed, strata too sparsely populated")
f=colnames(d)
w=ncol(d)-2
f=f[(1:w)]
w=ncol(gg)
q=c(1:w)
w=length(f)
q[(1:w)]=f
colnames(gg)=q
gg=clean(gg)
return(gg)
}




weightcont=function(Matrix){
  resd=matrix(nrow=5)
  a=ncol(Matrix)-5
  b=ncol(Matrix)
  r=Matrix[,(a:b)]
  f=matrix(nrow=nrow(r))
  r=cbind(r,f)
  t=matrix(nrow=nrow(r),ncol=ncol(r))
  for(i in 1:ncol(t)){
    t[,i]=as.numeric(r[,i])
  }
  r=t
  sum1=sum(r[,1])
  r[,2]=(r[,1]/sum1)*r[,2]
  r[,4]=r[,1]/r[,4]
  sum2=sum(r[,4])
  r[,4]=r[,4]/sum2
  r[,5]=r[,5]*r[,4]
  res=sum(r[,2])-sum(r[,5])
  a=r[,2]-r[,5]
  f=t.test(a,mu=0)
  resd[1,1]=res
  resd[2,1]=as.numeric(f[3])
  resd[3,1]=round(sum(r[,1]),1)
  resd[4,1]=f$conf.int[2]*resd[3,1]
  resd[5,1]=f$conf.int[1]*resd[3,1]
  rownames(resd)=c("Average Of Cases - Average Of Controls" , "t-test p-value", "No. Of Cases Matched", "95% C.I Upper Bound", "95% C.I Lower Bound")
  return(resd)

}


oddscont=function(mat){

  cola=ncol(mat)
  colb=ncol(mat)-5
  t=as.matrix(mat[,(colb:cola)])
  t2=as.matrix(mat[,-(colb:cola)])
  l=matrix(nrow=nrow(t),ncol=ncol(t))
  for(i in 1:ncol(l)){
    l[,i]=as.numeric(t[,i])
  }
  t=l

  en=matrix(nrow=1,ncol=1)
  for( kk in 1:ncol(t2)){
    d=na.omit(as.numeric(unique(t2[,kk])))
    if ( length(d) == 0){
      s=unique(t2[,kk])
      s=as.matrix(s,byrow=T)
      en=rbind(en,s)
    }
    else{
      s=colnames(t2)[kk]
      a=as.matrix(s)
      en=rbind(en,s)
    }
  }
  en=na.omit(en)
  results=matrix(nrow=nrow(en),ncol=2)
  results[,2]=rep(1)
  colnames(results)=c("Before Weighting" , "After Weighting")
  count=1
  for ( kk in 1:ncol(t2)){
    d=na.omit(as.numeric(unique(t2[,kk])))
    if ( length(d) == 0 ){
      s=unique(t2[,kk])

      for ( i in 1:length(s)){
        dd=as.matrix(t2[,kk])
        for (o in 1:nrow(dd)){
          if (dd[o,1] != s[i])dd[o,1]=0
          if (dd[o,1] == s[i])dd[o,1]=1
        }
        dd=as.matrix(as.numeric(dd[,1]))
        a=t[,4]*dd[,1]
        b=t[,1]*dd[,1]
        a1=sum(a)/sum(t[,4])
        b1=sum(b)/sum(t[,1])
        f=a1/b1
        results[count,1]=f
        count=count+1
      }
    }
    if (length(d) != 0){
      dd=as.matrix(t2[,kk])
      dd=as.matrix(as.numeric(dd[,1]))
      a=t[,4]*dd[,1]
      b=t[,1]*dd[,1]
      a1=sum(a)/sum(t[,4])
      b1=sum(b)/sum(t[,1])
      f=a1/b1
      results[count,1]=f
      count=count+1
    }

  }
  rownames(results)=en[,1]
  results=na.omit(results)
  return(results)
}

dispcont=function(mat){
  mat=oddscont(mat)
  time=c(1:nrow(mat))
  betagal.abs=c(mat[,1])
  cell.density=c(mat[,2])
  h=max(mat[,1])
  h=h+1
  plot(time, betagal.abs, pch=17, axes=F, ylim=c(0,h), xlab="", ylab="",col="black",main="Balance Of Covariates")
  axis(2, ylim=c(0,1),col="black")
  mtext("Odds",side=2,line=2.5)
  box()

  par(new=T)

  plot(time, cell.density, pch=15,  xlab="", ylab="", ylim=c(0,h), axes=F, type="b", col="red")

  bb=nrow(mat)
  axis(1,1:bb,rownames(mat),las=2)
  legend("topright", legend=c("Before Weighting","After Weighting"), pch=c(17,15))
}





gerom=function(d,rer) {
  z=ncol(d)-2
  f=merge(d, rer[,(1:z)], all.y=TRUE)
  l=as.matrix(f)
  k=as.matrix(rer)
  gg=matrix(nrow=1,ncol=(z+6))


  for (w in 1:nrow(k)){
    lala=l
    bb=matrix(nrow=nrow(lala),data=c(1:nrow(lala)))
    for(ww in 1:z){
      for ( i in 1:nrow(lala)){
        if(lala[i,ww] != k[w,ww])bb[i,1]=NA
        if(lala[i,ww] != k[w,ww])lala[i,1]=NA

      }
      lala=na.omit(lala)
      bb=na.omit(bb)
    }


    for(i in 1:nrow(bb)){
      num=bb[i,1]
      l[num,1]=NA
    }
    l=na.omit(l)

    p=lala
    p2=p
    for(i in 1:nrow(p)){
      if(p[i,(z+1)] == 0)p[i,1]=NA
      else(p2[i,1] = NA)
    }
    p=na.omit(p)
    p2=na.omit(p2)
    ll=as.numeric(p[,(z+2)])
    ll2=as.numeric(p2[,(z+2)])
    a=matrix(nrow=1,ncol=(z+6),data=rep(0))
    a[,(1:z)]=k[w,(1:z)]
    a[1,(z+1)]=nrow(p)

  a[1,(z+2)]=mean(ll)
    a[1,(z+3)]=0
    a[1,(z+4)]=nrow(p2)

  a[1,(z+5)]=mean(ll2)
    a[1,(z+6)]=0
    gg=rbind(gg,a)

  }
  gg=na.omit(gg)
  return(gg)
}

makeit=function(d){
  z=ncol(d)
  z=z-2
  h=z+1
  t=ddply(d,(1:z),nrow)
  if(nrow(t) < 100){
    a=matrix(nrow=1,ncol=2)
    a[1,1]=1
    a[1,2]=nrow(t)
    return(a)
  }
  else {
  m=trunc(nrow(t)/50) + 1
  a=matrix(nrow=m,ncol=2)
  a[1,1]=1
  a[1,2]=50
  for(i in 2:nrow(a)){
    a[i,1]=a[i-1,1]+50
    a[i,2]=a[i-1,2]+50
  }
  c=nrow(a)
  a[c,2]=(nrow(t))
  return(a)
  }
}

un1=function(d){
  d=as.data.frame(d)
  z=ncol(d)
  z=z-2
  t=ddply(d,(1:z),nrow)
  h=ncol(t)
  for(i in 1:nrow(t)){
    if(t[i,h] == 1)t[i,h]=NA
  }
  t=na.omit(t)
  f=merge(d, t[,(1:z)], all.y=TRUE)
  d=f
  return(d)
}

un3=function(a,d){
  z=ncol(d)
  z=z-2
  h=z+1
  tt=ddply(d,(1:z),nrow)
  gg=matrix(nrow=1,ncol=(z+6))
  for(i in 1:nrow(a)){
    rer=tt[(a[i,1]:a[i,2]),]
    ly=gerom(d,rer)
    gg=rbind(gg,ly)

  }
  gg=na.omit(gg)
  return(gg)
}







clean=function(g){
  r=ncol(g)-5
  r1=ncol(g)
  f=as.matrix(g[,(r:r1)])
  for(i in 1:nrow(f)){
    if(f[i,2] == "NaN")g[i,1]=NA
    if(f[i,5]=="NaN")g[i,5]=NA
  }
  g=na.omit(g)
  colnames(g)[(r:r1)]=c("Cases","Avg.","Sd(NA)","Controls","Avg.","Sd(NA)")
  return(g)
}



summarydisc=function(mat){
  colb=ncol(mat)-8
  t=as.matrix(mat[,(1:colb)])
  c=0
  for(i in 1:ncol(t)){
    u=as.matrix(unique(t[,i]))
    if(nrow(u) != 2)c=c+1
  }
  if ( c > 0){
    ll=matrix(nrow=c)
    c=1
    for(i in 1:ncol(t)){
      u=as.matrix(unique(t[,i]))
      if(nrow(u) != 2)ll[c,1]=i
      if(nrow(u) != 2)c=c+1
    }
    gg=matrix(nrow=nrow(t),ncol=nrow(ll))
    for(pp in 1:nrow(ll)){
      a=ll[pp,1]
      for(i in 1:nrow(t)){
        y=sprintf("%s %s", colnames(t)[a] , t[i,a])
        gg[i,pp]=y
      }
      t[1,a]=NA
    }
    t=t(t)
    t=na.omit(t)
    t=t(t)
    results=matrix(nrow=nrow(t))
    for( i in 1:nrow(results)){
      y="_"
      for(kk in 1:ncol(gg)){
        y=sprintf("%s %s" , y,gg[i,kk])
      }
      for(kk in 1:ncol(t)){
        if (t[i,kk] != 0) y=sprintf("%s, %s" , y,colnames(t)[kk])
      }
      results[i,1]=y
    }
  }
  if ( c == 0) {
    results=matrix(nrow=nrow(t))
    for( i in 1:nrow(results)){
      y="_"
      for(kk in 1:ncol(t)){
        if (t[i,kk] != 0) y=sprintf("%s, %s" , y,colnames(t)[kk])
      }
      results[i,1]=y
    }
  }
  if ( c > 0 & c == ncol(t)){
    ll=matrix(nrow=c)
    c=1
    for(i in 1:ncol(t)){
      u=as.matrix(unique(t[,i]))
      if(nrow(u) != 2)ll[c,1]=i
      if(nrow(u) != 2)c=c+1
    }
    gg=matrix(nrow=nrow(t),ncol=nrow(ll))
    for(pp in 1:nrow(ll)){
      a=ll[pp,1]
      for(i in 1:nrow(t)){
        y=sprintf("%s %s", colnames(t)[a] , t[i,a])
        gg[i,pp]=y
      }
      results=matrix(nrow=nrow(t))
      for( i in 1:nrow(results)){
        y="_"
        for(kk in 1:ncol(gg)){
          y=sprintf("%s %s" , y,gg[i,kk])
        }
        results[i,1]=y
      }
    }
  }
  colnames(results)=c("Strata Characteristics")
  cola=ncol(mat)
  colb=ncol(mat)-7
  mat=as.matrix(mat[,(colb:cola)])
  colnames(mat)=c(" No. Of Controls", "Control Weights" , "Controls Y=1", "Controls Y=0", "No. Of Cases", "Case Weights", "Cases Y=1","Cases Y=0")
  mat=cbind(results,mat)
  f=as.matrix(mat[,(2:5)])
  f=f[,-4]
  f2=as.matrix(mat[,6:9])
  f2=f2[,-4]
  mat=cbind(results,f2,f)
  return(mat)
}


summarycont=function(mat){
  colb=ncol(mat)-6
  t=as.matrix(mat[,(1:colb)])
  c=0
  for(i in 1:ncol(t)){
    u=as.matrix(unique(t[,i]))
    if(nrow(u) != 2)c=c+1
  }
  if ( c > 0 & c != ncol(t)){
    ll=matrix(nrow=c)
    c=1
    for(i in 1:ncol(t)){
      u=as.matrix(unique(t[,i]))
      if(nrow(u) != 2)ll[c,1]=i
      if(nrow(u) != 2)c=c+1
    }
    gg=matrix(nrow=nrow(t),ncol=nrow(ll))
    for(pp in 1:nrow(ll)){
      a=ll[pp,1]
      for(i in 1:nrow(t)){
        y=sprintf("%s %s", colnames(t)[a] , t[i,a])
        gg[i,pp]=y
      }
      t[1,a]=NA
    }
    t=t(t)
    t=na.omit(t)
    t=t(t)
    results=matrix(nrow=nrow(t))
    for( i in 1:nrow(results)){
      y="_"
      for(kk in 1:ncol(gg)){
        y=sprintf("%s %s" , y,gg[i,kk])
      }
      for(kk in 1:ncol(t)){
        if (t[i,kk] != 0) y=sprintf("%s, %s" , y,colnames(t)[kk])
      }
      results[i,1]=y
    }
  }
  if ( c == 0) {
    results=matrix(nrow=nrow(t))
    for( i in 1:nrow(results)){
      y="_"
      for(kk in 1:ncol(t)){
        if (t[i,kk] != 0) y=sprintf("%s, %s" , y,colnames(t)[kk])
      }
      results[i,1]=y
    }
  }
  if ( c > 0 & c == ncol(t)){
    ll=matrix(nrow=c)
    c=1
    for(i in 1:ncol(t)){
      u=as.matrix(unique(t[,i]))
      if(nrow(u) != 2)ll[c,1]=i
      if(nrow(u) != 2)c=c+1
    }
    gg=matrix(nrow=nrow(t),ncol=nrow(ll))
    for(pp in 1:nrow(ll)){
      a=ll[pp,1]
      for(i in 1:nrow(t)){
        y=sprintf("%s %s", colnames(t)[a] , t[i,a])
        gg[i,pp]=y
      }
      results=matrix(nrow=nrow(t))
      for( i in 1:nrow(results)){
        y="_"
        for(kk in 1:ncol(gg)){
          y=sprintf("%s %s" , y,gg[i,kk])
        }
        results[i,1]=y
      }
    }
  }
  colnames(results)=c("Strata Characteristics")
  cola=ncol(mat)
  colb=ncol(mat)-5
  mat=as.matrix(mat[,(colb:cola)])
  mat[1,3]=NA
  mat[1,6]=NA
  mat=t(mat)
  mat=na.omit(mat)
  mat=t(mat)
  colnames(mat)=c(" No. Of Cases", "Weighted Average Response" , "No. Of Controls", "Weighted Average Response")
  mat=cbind(results,mat)
  return(mat)
}



