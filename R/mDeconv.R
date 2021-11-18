#' Title
#'
#' @param XEXP
#' @param gL
#' @param addG
#' @param addWhiteComp
#'
#' @return
#' @export
#'
#' @examples
mDeconv <- function(XEXP,gL,addG=NULL,addWhiteComp=F,
                    maxIter=2000L,nthreads=1,printInfo=F){


  gL=lapply(gL,intersect,y=rownames(XEXP))
  if(min(sapply(gL,length))<2){
    stop("At least 2 marker genes should be present in each signature and present in the expression data")
  }
  if(max(table(unlist(gL)))>1){
    stop("Each gene in MArker gene list should only be asscoiated to one phenotype")
  }


  allselg=unique(unlist(c(gL,addG)))
  MARK=do.call(cbind,lapply(gL,function(x){allselg %in% x}))
  rownames(MARK)=allselg
  colnames(MARK)=names(gL)


  if(addWhiteComp){
    MARK=cbind(MARK,"additional"=rep(F,nrow(MARK)))
  }

  selg=rownames(MARK)

  X=XEXP[selg,]
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
  D = diag(rowSums(mynmf$H))
  NewH = MASS::ginv(D)%*% mynmf$H
  NewW= W %*%D
  NewHtprop=t(NewH)*(1/colSums(NewH))
  colnames(NewHtprop)=colnames(Ht)
  NewHtprop
  list(H=Ht,W=W,Hprop=NewHtprop,Wprop=NewW,NMF=mynmf)

}






