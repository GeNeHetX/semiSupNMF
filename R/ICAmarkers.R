
.internaldoica = function(datanorm,k = 4,silent=F,maxiter = 1000000,eps = 10^-6,addmsg=""){
  X=as.matrix(datanorm)
  resICA=NULL

  try({

    if(!silent)print(paste("Do ICA",addmsg))
    resICA <- JADE::JADE(X, n.comp = k, maxiter = maxiter, eps = eps)
    rownames(resICA$A) <- colnames(X)
    rownames(resICA$S) <- rownames(X)


    colnames(resICA$S)=paste("ICA",1:ncol(resICA$S),sep="")
    colnames(resICA$A)=paste("ICA",1:ncol(resICA$A),sep="")

  })

  resICA
}



#' doICA: Independent Compoent analysis for marker gene identification
#'
#' @param X A matrix of gene expression
#' @param k the number of components to extract
#'
#' @return a list:
#'           - initICA : the internals of the ICA
#'           - dirS: the gene weights S oriented towards the direction most likely corresponding to cell type
#'           - dirA: the sample weights A oriented towards the direction most likely corresponding to cell type
#'           - orient: a dataframe with data used for component orientation
#'
#' @examples
doICA=function(X,k){
  coreica= semiSupNMF:::.internaldoica(normalize(X,"gc"),k)
  orientdf=do.call(rbind,apply(coreica$S,2,function(x){
    ex=range(cumsum(x[order(-abs(x))]))
    data.frame(peak=ex[which.max(abs(ex))],revpeak=ex[which.min(abs(ex))],
               low=min(ex),high=max(ex),extreme=sign(x[which.max(abs(x))]))
  }))

  dirS=t(t(coreica$S)*orientdf$extreme)
  dirA=t(t(coreica$A)*orientdf$extreme)
  list(initICA=coreica,dirS=dirS,dirA=dirA,orient=orientdf)
}



#' getMarkerFromICA
#'
#' @param anica  An ICA object obtained using doICA
#' @param thresh the threshold to apply to S to select component specific markers
#'
#' @return a list of gene sets
getMarkerFromICA=function(anica,thresh=1){
  Xs=t(apply((anica$dirS),1,sort,decreasing=T))
  ig=apply(anica$dirS[which( (Xs[,1]-Xs[,2]) > thresh ),],1,which.max)
  ign=colnames(anica$dirS)[ig]
  gMarkL=split(names(ig),ign)

  min(rowMaxs(anica$dirS[names(ig),]))

  nmiss=sum(sapply(gMarkL[colnames(anica$dirS)],length)==0)


  print(paste("Of a total of",ncol(anica$dirS),"components,",nmiss,"are missing. (",
              paste(colnames(anica$dirS)[which(sapply(gMarkL[colnames(anica$dirS)],length)==0)],collapse=","),")"))


  gMarkL
}


corscoreGS=function(anexp,gsmark,cormeth="pearson"){
  y=sapply(gsmark,function(gs){
    comgs=intersect(gs,rownames(anexp))
    xc=cor(t(anexp[comgs,]),method=cormeth)
    return(narm(xc[upper.tri(xc)]))
  })
  y[order(sapply(y,mean),decreasing=T)]
}




checkMarkGene=function(anexp,gsmark){
  newgmark=lapply(gsmark,intersect,y=rownames(anexp))
  dupg=unlist(newgmark)[which(duplicated(unlist(newgmark)))]
  newgmark=lapply(newgmark,setdiff,y=dupg)
  newgmark[which(sapply(newgmark,length)>=3)]
}
