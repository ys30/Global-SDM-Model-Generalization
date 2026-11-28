# Descriptive coefficient uncertainty: resample whole species within regions.
dir.create("results/mechanism_uncertainty",recursive=TRUE,showWarnings=FALSE)
b<-read.csv("inputs/matched_features.csv")
vars<-c("log_n","logit_prev","log_geo","log_env","outside_fraction")
models<-c("quadratic","ridge_linear","ridge_quadratic","xgb_shallow","gam_low")
keys<-c("region","group","species","fold")
stopifnot(nrow(b)==876,!anyDuplicated(b[keys]))
b$id<-interaction(b$region,b$group,b$species,drop=TRUE)
stopifnot(length(unique(b$id))==198)
X<-model.matrix(~factor(region)+log_n+logit_prev+log_geo+log_env+outside_fraction,b)
Y<-as.matrix(b[paste0("auc.",models)])-b$auc.linear
w<-b$w
stopifnot(all(is.finite(X)),all(is.finite(Y)),all(w>0))
fit<-function(w){qr.solve(crossprod(X,w*X),crossprod(X,w*Y))}
ref<-fit(w)
ids<-split(as.character(b$id),b$region)
ids<-lapply(ids,unique)
set.seed(20261128)
boot<-array(NA_real_,c(length(vars),length(models),4000))
for(k in 1:4000){
 cnt<-table(unlist(lapply(ids,function(z)sample(z,length(z),replace=TRUE))))
 mult<-as.numeric(cnt[as.character(b$id)]);mult[is.na(mult)]<-0
 cf<-tryCatch(fit(w*mult),error=function(e)NULL)
 if(!is.null(cf))boot[,,k]<-cf[match(vars,colnames(X)),,drop=FALSE]
}
out<-list()
for(j in seq_along(models))for(i in seq_along(vars)){
 z<-boot[i,j,];z<-z[is.finite(z)];stopifnot(length(z)>=3900)
 q<-quantile(z,c(.025,.975))
 out[[length(out)+1]]<-data.frame(model=models[j],term=vars[i],coefficient=ref[match(vars,colnames(X)),j],lower=q[1],upper=q[2],valid_boot=length(z))
}
o<-do.call(rbind,out)
loo<-read.csv("inputs/coefficients.csv")
l<-loo[loo$excluded!="none",]
ll<-aggregate(coefficient~model+term,l,function(x)c(min=min(x),max=max(x)))
ll<-data.frame(model=ll$model,term=ll$term,loo_min=ll$coefficient[,1],loo_max=ll$coefficient[,2])
o<-merge(o,ll,by=c("model","term"))
print(o)
write.csv(o,"results/mechanism_uncertainty/intervals.csv",row.names=FALSE)
# These are conditional species-resampling intervals, not site-independent inference.
