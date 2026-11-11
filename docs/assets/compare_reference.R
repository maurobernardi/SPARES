# Reproducible accuracy and timing comparisons, on matched objectives/scales.
# Run: Rscript compare_reference.R [output-directory]
library(spares)
stopifnot(requireNamespace("glmnet",quietly=TRUE),requireNamespace("gglasso",quietly=TRUE))
set.seed(20261111)
n<-100;p<-12;groups<-rep(1:6,each=2)
X<-matrix(rnorm(n*p),n,p)
# Within-group orthonormalization; no intercept/centering in either implementation.
for(g in unique(groups))X[,groups==g]<-qr.Q(qr(X[,groups==g]))*sqrt(n)
y<-as.numeric(X%*%c(1,-.8,.6,0,rep(0,p-4))+rnorm(n,sd=.5))
lambda<-c(.25,.12,.06,.03)
ctl<-list(maxit=30000,abstol=1e-10,reltol=1e-10,print.out=FALSE)
pilot<-as.numeric(qr.coef(qr(X),y));wi<-exp(-abs(pilot))
wg<-vapply(unique(groups),function(g)exp(-sqrt(sum(pilot[groups==g]^2))),0)
ours<-list(
 lasso=function()lmreg(X,y,lambda=lambda,control=ctl),
 glasso=function()lmreg(X,y,penalty='glasso',groups=groups,lambda=lambda,control=ctl),
 adalasso=function()lmreg(X,y,penalty='adalasso',profile='exponential',lambda=lambda,control=ctl),
 adaglasso=function()lmreg(X,y,penalty='adaglasso',groups=groups,profile='exponential',lambda=lambda,control=ctl))
refs<-list(
 lasso=function()glmnet::glmnet(X,y,lambda=lambda,intercept=FALSE,standardize=FALSE,thresh=1e-12),
 glasso=function()gglasso::gglasso(X,y,group=groups,lambda=lambda,intercept=FALSE,loss='ls',eps=1e-12),
 # glmnet rescales penalty.factor to sum p: compensate its lambda accordingly.
 adalasso=function()glmnet::glmnet(X,y,lambda=lambda*mean(wi),penalty.factor=wi,intercept=FALSE,standardize=FALSE,thresh=1e-12),
 adaglasso=function()gglasso::gglasso(X,y,group=groups,pf=wg,lambda=lambda,intercept=FALSE,loss='ls',eps=1e-12))
rows<-list()
for(m in names(ours)) {
 a<-suppressWarnings(ours[[m]]());b<-refs[[m]]();bb<-as.matrix(b$beta)
 stopifnot(ncol(bb)==length(lambda),all(a$convergence==1))
 ta<-median(replicate(5,system.time(suppressWarnings(ours[[m]]()))[['elapsed']]))
 tb<-median(replicate(5,system.time(refs[[m]]())[['elapsed']]))
 rows[[m]]<-data.frame(method=m,reference=if(grepl('glasso',m))'gglasso'else'glmnet',
  max_coefficient_error=max(abs(a$beta-bb)),max_prediction_error=max(abs(X%*%(a$beta-bb))),
  spares_seconds=ta,reference_seconds=tb)
}
results<-do.call(rbind,rows);print(results,row.names=FALSE)
stopifnot(all(results$max_prediction_error<1e-4))
args<-commandArgs(trailingOnly=TRUE);out<-if(length(args))args[1]else'.'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
utils::write.csv(results,file.path(out,'reference-comparison.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'benchmark-session.txt'))
# Timings include the public wrapper and retained metadata, exclude DoF,
# use five warm runs, and are illustrative rather than a speed ranking.
