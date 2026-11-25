library(disdat)
library(glmnet)
dir.create("results/ridge",recursive=TRUE,showWarnings=FALSE)
regions<-c("CAN","NZ","SWI")
auc<-function(y,p){n<-sum(y==1);m<-sum(y==0);if(!n||!m)return(NA_real_);(sum(rank(p)[y==1])-n*(n+1)/2)/(n*m)}
set.seed(20261125)
out<-list()
for(r in regions){
 pa<-disPa(r);env<-disEnv(r);stopifnot(identical(pa[,1:4],env[,1:4]))
 pred<-intersect(disPredictors(r),names(env));pred<-pred[vapply(env[pred],is.numeric,logical(1))]
 pred<-setdiff(pred,"vegsys")
 x0<-as.matrix(env[pred]);ok<-apply(x0,1,function(z)all(is.finite(z)))&is.finite(pa$x)&is.finite(pa$y)
 cuts<-unique(as.numeric(quantile(pa$x[ok],seq(0,1,length.out=6))))
 fold<-cut(pa$x,cuts,include.lowest=TRUE,labels=FALSE)
 # Fixed candidate cohort: all species with at least 100 presences and absences.
 for(s in names(pa)[-(1:4)]){
 y<-pa[[s]];valid<-ok&!is.na(y)&y%in%c(0,1)
 if(sum(y[valid]==1)<100||sum(y[valid]==0)<100)next
 for(k in c(1,3,5)){
 pool<-which(valid&fold!=k);te<-which(valid&fold==k)
 if(length(pool)<2000||length(unique(y[te]))<2)next
 for(rep in 1:3){
 order<-sample(pool)
 for(n in c(100,500,2000)){
 tr<-head(order,n)
 if(min(table(factor(y[tr],levels=0:1)))<5)next
 mu<-colMeans(x0[tr,,drop=FALSE]);sdv<-apply(x0[tr,,drop=FALSE],2,sd);keep<-is.finite(sdv)&sdv>0
 if(!any(keep))next
 x<-scale(x0[,keep,drop=FALSE],center=mu[keep],scale=sdv[keep])
 for(structure in c("linear","quadratic")){
 z<-if(structure=="quadratic")cbind(x,x^2)else x
 for(penalty in c("none","ridge")){
 w<-character()
 p<-tryCatch(withCallingHandlers({
 if(penalty=="none"){
 fit<-glm.fit(cbind(1,z[tr,,drop=FALSE]),y[tr],family=binomial(),control=glm.control(maxit=50))
 if(!fit$converged)stop("nonconvergence")
 cf<-fit$coefficients;cf[is.na(cf)]<-0
 plogis(cbind(1,z[te,,drop=FALSE])%*%cf)
 }else{
 # Fixed penalty sensitivity, not tuned using held-out data.
 fit<-glmnet(z[tr,,drop=FALSE],y[tr],family="binomial",alpha=0,lambda=.01,standardize=TRUE)
 as.numeric(predict(fit,z[te,,drop=FALSE],s=.01,type="response"))
 }
 },warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart("muffleWarning")}),error=function(e){w<<-c(w,conditionMessage(e));rep(NA_real_,length(te))})
 good<-all(is.finite(p))
 out[[length(out)+1L]]<-data.frame(region=r,species=s,fold=k,rep=rep,n=n,train_presence=sum(y[tr]==1),structure=structure,penalty=penalty,auc=if(good)auc(y[te],p)else NA_real_,success=good,warnings=paste(unique(w),collapse=";"))
 }
 }
 }
 }
 }
 }
}
a<-do.call(rbind,out);write.csv(a,"results/ridge/metrics.csv",row.names=FALSE)
a$model<-paste(a$structure,a$penalty,sep="_")
b<-reshape(a[c("region","species","fold","rep","n","model","auc")],idvar=c("region","species","fold","rep","n"),timevar="model",direction="wide")
b$delta_none<-b$auc.quadratic_none-b$auc.linear_none
b$delta_ridge<-b$auc.quadratic_ridge-b$auc.linear_ridge
write.csv(b,"results/ridge/paired.csv",row.names=FALSE)
# Strictly matched cohort across all sample sizes and four models.
b<-b[complete.cases(b[,c("delta_none","delta_ridge")]),]
key<-interaction(b$region,b$species,b$fold,b$rep,drop=TRUE)
b<-b[ave(b$n,key,FUN=function(x)length(unique(x)))==3,]
stopifnot(nrow(b)>0)
c<-aggregate(cbind(delta_none,delta_ridge)~region+species+n,b,mean)
print(aggregate(cbind(delta_none,delta_ridge)~region+n,c,mean))
print(aggregate(cbind(delta_none,delta_ridge)~n,c,mean))
write.csv(c,"results/ridge/species_matched.csv",row.names=FALSE)
