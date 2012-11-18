cmBE <- function(N, pbar,px,n,c,type=c("LTPD","AOQL"))
{ type = match.arg(type);
e2 <- function(N, cm, pbar, n, k,n2_,c_) 100*(n*cm+(N-n)*(1-OC(pbar,n,k)))/(N-(N-n2_)*LAtHy(pbar,n2_,c_,N));
LAtHy <-function(p,n_,c_,N) phyper(q=c_, m=N*p, n=N*(1-p), k=n_);
cm_limit =  uniroot(function(cm) {
    if (type == "LTPD"){
planMerenim=planLTPD(N=N,pt=px,pbar=pbar,cm=cm)}
 if (type == "AOQL"){
planMerenim=planAOQL(N=N,pbar=pbar,pL=px,cm=cm)}
e2(N, cm,   pbar,n(planMerenim),k(planMerenim),n,c) - 100}, c(1, 9),tol=.Machine$double.eps)$root;
return(cm_limit)}
