#' Title
#'
#' @param XEXP A gene expression matrix. Should be non-negative values (log2 normalized counts, counts)
#' @param gL A list of gene ID/names, all found in the rownames of XEXP. Each entry in the
#' @param addG A set of gene ID/names,
#' @param addWhiteComp
#'
#' @return
#' @export
#'
#' @examples
mDeconv <- function(XEXP,gL,addG=NULL,addWhiteComp=T,
                    maxIter=2000L,nthreads=1,printInfo=F,testInternalUnlog=FALSE){

  #Run everything in log2 scale
  doLog=TRUE

  gL=lapply(gL,intersect,y=rownames(XEXP))
  if(min(sapply(gL,length))<2){
    stop("At least 2 marker genes  (arg: gL) should be present in each signature and present in the expression data")
  }
  if(max(table(unlist(gL)))>1){
    stop("Each gene in MArker gene list (arg: gL) should only be asscoiated to one phenotype")
  }

  if(length(gL)>100){
    stop("More than 100 gene lists seems unreasonable.")
  }

  if(length(unique(names(gL))) < length(gL)){
    names(gL)= paste0("geneL",combn(letters,2,paste,collapse=""))[1:length(gL)]
  }


  allselg=unique(unlist(c(gL,addG)))
  MARK=do.call(cbind,lapply(gL,function(x){allselg %in% x}))
  rownames(MARK)=allselg
  colnames(MARK)=names(gL)


  if(addWhiteComp){
    MARK=cbind(MARK,"additional"=rep(F,nrow(MARK)))
  }

  selg=rownames(MARK)

  X=as.matrix(XEXP[selg,])


  if(testInternalUnlog){
   X=semiSupNMF:::.qcpm(X,testInternalUnlog,doLog)
  }
  initwmax=(MARK*matrixStats::rowMaxs(X)) #+ ((abs(MARK-1)*matrixStats::rowMins(X)))

  nmark=colSums(MARK)
  nmark[which(nmark==0)]=1
  initHraw=(t(MARK)%*%X)/  nmark
  if(addWhiteComp){
    initHprop=t(t(initHraw)/colSums(initHraw))
    initHprop[nrow(initHprop),]=0.01
  }else{
    initHprop=t(t(initHraw)/colSums(initHraw))
  }



  if(printInfo){
    print("init W")
    print(dim(initwmax))
    print("init H")
    print(dim(initHprop))

    print("markers W")
    print(dim(MARK))
    print("markers H")
    print(dim(matrix(F,ncol(MARK),ncol(X))))
    print("K")
    print(ncol(MARK))

  }


  mynmf=semiSupNMF:::semiSupNMF(A=X, k = ncol(MARK),
                                init = list(W=initwmax,H=initHprop), mask = list(W=MARK,H=matrix(F,ncol(MARK),ncol(X))),
                                check.k = TRUE,max.iter = maxIter, rel.tol = 1e-4, n.threads = nthreads, trace = 0.1,
                                verbose = 0, show.warning = TRUE, inner.rel.tol = 1e-9)



  W=mynmf$W
  colnames(W)=colnames(MARK)
  Ht=t(mynmf$H)
  colnames(Ht)=colnames(MARK)
  # normfac=(sapply(names(mark),\(icn){(sum((Wp[mark[[icn]],icn]) ))}))#/colSums(W[,names(mark)])
  #
  # Hn=t(t(Hp[,1:length(mark)])*normfac)
  # Hpn = MASS::ginv(diag(rowSums(Hpn)))%*% Hpn
  # D = diag(rowSums(mynmf$H))
  # NewH = MASS::ginv(D)%*% mynmf$H
  # NewW= W %*%D
  # NewHtprop=t(NewH)*(1/colSums(NewH))
  # colnames(NewHtprop)=colnames(Ht)
  #
  # Wnorm=semiSupNMF:::.qcpm(W[,which(colSums(W)>0)], unlog=F,outlog=isLogged2)
  # normWISP=NULL
  # try({normWISP=semiSupNMF:::.qWISProutine(Wnorm,X)})

  # list(H=Ht,W=W,Hp=NewHtprop,Wp=NewW,Wnorm=Wnorm,normWISP=normWISP,NMF=mynmf,expUsed=X)


  normwh=normByWSum(mynmf,mean(colSums(X)))

  list(H=Ht,W=W,normH=normwh$H,normW=normwh$W,NMF=mynmf,expUsed=X,markers=MARK,CVH=sd(rowSums(Ht))/mean(rowSums(Ht)))
}


#' normByWSum normalize NMF H by sum of values in W. Highly experimental
#'
#' @param nmf
#' @param sumto
#'
#' @return normalized H and W
#' @export
#'
#' @examples
normByWSum=function(nmf,sumto=10000){
  H=nmf$H
  W=nmf$W
  D = (colSums(W))/sumto
  iD=1/D
  iD[which(is.infinite(iD))]=0
  list(H=t(D* H),W=t(t(W) * iD))
}

#' Normalize H (supposedly from NMF) to get sum to 1
#'
#' @param H H from NMF
#' @param normFactor a vector of normalization fators (H columns are multiplied by normFactor before sumto1)
#'
#' @return H with lines sum to 1
#' @export
#'
#' @examples
normSumTo1=function(H,normFactor=NULL){
  if(!is.null(normFactor)){
    if(length(normFactor)==ncol(H)){
      H=t(t(H)*normFactor)
    }
  }
  H*(1/rowSums(H))
}





.qcpm=function(x,unlog=F,outlog=F){
  if(unlog){
    x=((2^(x))-1)
  }
  y=edgeR::cpm(x, normalized.lib.sizes = TRUE,log = F)
  if(outlog){y=log2(y+1)}
  y
}


.qWISProutine=function(W,expdat){

  W=W[,which(colSums(W)>0)]

  if(any(rownames(expdat) != rownames(W))){
    stop(paste0("In the WISP routine, expression data and gene weight matrix should",
                " have exactly the same gene/row names in the same order"))
  }


  wisped=NULL
  subcomps=setdiff(colnames(W),"additional")
  try({wisped=semiSupNMF:::.qWISP(expdat,W[,subcomps],scaling="none",sum_LessThanOne=T)},silent=T)
  cpt=0
  while(is.null(wisped) & length(subcomps)>=2 & cpt <4) {
    cpt=cpt+1
    try({
      subcomps=setdiff(subcomps,subcomps[which.min(colMeans(W[,subcomps]))])
      if(length(  subcomps)>=2){
        wisped=semiSupNMF:::.qWISP(expdat,W[,subcomps],scaling="none",sum_LessThanOne=T)
      }
    },silent=T)
  }

  if(is.null(wisped)){  return(NULL) }

  return(
    list(
      prop  =wisped[,subcomps],
      scores=wisped[,setdiff(colnames(wisped),subcomps)]
    )
  )
}

.qWISP=function (data, centro, scaling = c("none", "scale", "center")[1],
                 sum_LessThanOne = TRUE)
{
  cutoff_ttest_weights=0.05
  g = intersect(rownames(data), rownames(centro))
  data = data[g, ]
  centro = as.matrix(centro[g, ])
  if (scaling == "scale") {
    data = scale(data, scale = TRUE)
    centro = scale(centro, scale = TRUE)
  }
  if (scaling == "center") {
    data = scale(data, scale = FALSE)
    centro = scale(centro, scale = FALSE)
  }
  nJ = ncol(centro)
  A = centro
  nSubj = ncol(data)
  mixCoef = matrix(0, nSubj, nJ)
  rownames(mixCoef) = colnames(data)
  colnames(mixCoef) = colnames(centro)
  Amat = cbind(rep(-1, nJ), diag(nJ))
  b0vec = c(-1, rep(0, nJ))
  if (sum_LessThanOne == TRUE) {
    meq = 0
  }
  else {
    meq = 1
  }
  output = data.frame(t(apply(data, 2, function(y) {
    obs = which(!is.na(y))
    Dmat = t(centro[obs, ]) %*% centro[obs, ]
    diag(Dmat) <- diag(Dmat) + 1e-08
    mixCoef = quadprog::solve.QP(Dmat, t(centro[obs, ]) %*%
                                   y[obs], Amat, b0vec, meq = meq)$sol
    B = as.matrix(y)
    coeff = round(mixCoef, 4)
    names(coeff) = colnames(centro)
    vBeta = matrix(coeff)
    dSigmaSq = sum((B - A %*% vBeta)^2)/(nrow(A) - ncol(A))
    dTotalSq = sum((B)^2)/(nrow(A))
    dModelSq = sum((A %*% vBeta)^2)/(ncol(A))
    mVarCovar = try(dSigmaSq * chol2inv(chol((t(A) %*% A)+10^-8 )))
    Adjusted.R.squared = round((dTotalSq - dSigmaSq)/dTotalSq,
                               3)
    Ftest = dModelSq/dSigmaSq
    p.Ftest = stats::pf(q = Ftest, df1 = ncol(A), df2 = nrow(A) -
                          ncol(A), lower.tail = FALSE)
    ng = nrow(centro)
    if (!is.character(mVarCovar)) {
      vStdErr = sqrt(diag(mVarCovar))
      CI.inf = sapply(1:nJ, function(j) {
        round(coeff[j] - (stats::qt(1 - cutoff_ttest_weights/2,
                                    ng - nJ) * vStdErr[j]), 2)
      })
      CI.sup = sapply(1:nJ, function(j) {
        round(coeff[j] + (stats::qt(1 - cutoff_ttest_weights/2,
                                    ng - nJ) * vStdErr[j]), 2)
      }
      )
      tvalue = sapply(1:nJ, function(j) {
        coeff[j]/vStdErr[j]
      })
      p.Ttest = sapply(1:nJ, function(j) {
        stats::pt(abs(coeff[j]/vStdErr[j]), df = ng -
                    nJ, lower = FALSE) * 2
      }
      )
    }
    else {
      mVarCovar = NA
      vStdErr = CI.inf = CI.sup = p.Ttest = rep(NA, nJ)
    }
    names(CI.inf) = paste("CI.inf.", colnames(centro), sep = "")
    names(CI.sup) = paste("CI.sup.", colnames(centro), sep = "")
    names(p.Ttest) = paste("Pvalue.", colnames(centro), sep = "")
    coeff.filtered = coeff
    coeff.filtered[p.Ttest > cutoff_ttest_weights] = 0
    names(coeff.filtered) = paste(names(coeff), ".filtered",
                                  sep = "")
    dist.Obs.Model = round(sqrt(sum(((centro %*% coeff) -
                                       y)^2)), 2)
    c(coeff, dist.Obs.Model = dist.Obs.Model, Ftest.pvalue = p.Ftest,
      Adjusted.R.squared = Adjusted.R.squared)
  }
  )))
  output$topWeightedClass = colnames(centro)[apply(output[,1:nJ], 1, which.max)]
  output$deltaTopWeights = apply(output[, 1:nJ], 1, function(z) abs(diff(sort(z,
                                                                              decreasing = TRUE)[1:2]))
  )

  return(output)
}




