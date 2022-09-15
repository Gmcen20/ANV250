vec<-function(x,p){
  l=nchar(p)
  v=0
  n=0
  for(i in 1:l){
    if(substr(p,i,i)==x){
      n=n+1
    }
    v[i]=n
  }
  return(v)
}

num<-function(v){
  l=length(v)
  return(v[l])
}

theta<-function(v){
  l=length(v)
  s=sum(v)
  if(num(v)==0){return (0)}
  else{return(s/l)}
}

mdis<-function(v){
  s=sum(v)
  if(num(v)==0){return (0)}
  else{return(s/num(v))}
}

d2<-function(v){
  l=length(v)
  m=theta(v)
  s=0
  for(i in 1:l){
    s=s+(v[i]-m)^2
  }
  if(num(v)==0){return (0)}
  else{return(s/(num(v)^2))}
}

cova<-function(v1,v2){
  l=length(v1)
  m1=theta(v1)
  m2=theta(v2)
  s=0
  for(i in 1:l){
    s=s+(v1[i]-m1)*(v2[i]-m2)
  }
  if(num(v1)*num(v2)==0){return(0)}
  else{return(s/(num(v1)*num(v2)))}
}

aa=c('A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V')
setwd("/Users/guanmengcen/Desktop/virus/giant_virus/fasta_new")
#??????sequence
mydata <- readLines("Coronaviridae.fasta")

id <- which(grepl(">",mydata))
n <- length(id)
pro <- list()
tmp <- ""
for(i in 1:(n-1)){
  tmp <- ""
  for(j in (id[i]+1):(id[i+1]-1)){
    tmp <- paste(tmp,mydata[j],sep = "")
  }
  pro[[i]] <- tmp
}
tmp <- ""
for(j in (id[n]+1):length(mydata)){
  tmp <- paste(tmp,mydata[j],sep = "")
}
pro[[n]] <- tmp



#natural vector
t=n
re=matrix(0,t,250)
for(k in 1:t){
  p=pro[[k]]
  v=matrix(0,20,nchar(p))
  for(i in 1:20){
    v[i,]=vec(aa[i],p)
  }
  for(i in 1:20){
    re[k,i]=num(v[i,])
    re[k,i+20]=mdis(v[i,])
    re[k,i+40]=d2(v[i,])
  }
  u=61
  for(i in 1:19){
    for(j in (i+1):20){
      re[k,u]=cova(v[i,],v[j,])
      u=u+1
    }
  }
}
write.csv(re,file="Coronaviridae.csv")