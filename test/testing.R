library(semiSupNMF)
library(reshape)
library(ggplot2)
set.seed(1)
m=15;n=20
centro=cbind(rpois(m,20)+(rpois(m,2)*100),rpois(m,20)+(rpois(m,2)*100))
dimnames(centro)=list(letters[1:m],LETTERS[1:2])
x=runif(n);props=cbind(x,1-x)
X=centro%*%t(props)
ismaxof=apply(centro,1,which.max)
maskofmax=t(sapply(ismaxof,function(i)(1:2) ==i))
markers=split(names(ismaxof),ismaxof)


mDeconv(X,markers)


#true stuff



library(semiSupNMF)
library(reshape)
library(ggplot2)

load("test/testdata.Rdata")



subssupprop=mDeconv(subsimulexp,selUPgeneL,addWhiteComp=F)

mean(abs(subsimulp-subssupprop$Hprop[,colnames(subsimulp)]))
mean(sqrt((subsimulp-subssupprop$Hprop[,colnames(subsimulp)])^2))
sort(diag(cor(subsimulp,subssupprop$Hprop[,colnames(subsimulp)])))

df=data.frame(melt(subsimulp),pred=melt(subssupprop$Hprop[,colnames(subsimulp)])[,3])
ggplot(df,aes(x=pred,y=value))+geom_point()+facet_wrap(~X2)+ geom_abline(slop=1)




selcomp=names(selUPgeneL)[1:10]
add=NULL;
add=unlist(selUPgeneL)
subssuppropMC=mDeconv(subsimulexp,selUPgeneL[selcomp],addG = ,addWhiteComp=T)

mean(abs(subsimulp[,selcomp]-subssuppropMC$Hprop[,selcomp]))
mean(sqrt((subsimulp[,selcomp]-subssuppropMC$Hprop[,selcomp])^2))
sort(diag(cor(subsimulp[,selcomp],subssuppropMC$Hprop[,selcomp])))


df=data.frame(melt(subsimulp[,selcomp]),pred=melt(subssuppropMC$Hprop[,selcomp])[,3])
ggplot(df,aes(x=pred,y=value))+geom_point()+facet_wrap(~X2)+ geom_abline(slop=1)





XEXP=subsimulexp;gL=selUPgeneL[selcomp];addG=unlist(selUPgeneL);addWhiteComp=T;
maxIter=2000L;nthreads=1;printInfo=F
