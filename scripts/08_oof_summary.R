auc <- function(y,p) { n<-as.double(sum(y==1)); m<-as.double(sum(y==0)); if(!n||!m)return(NA_real_); (sum(rank(p)[y==1])-n*(n+1)/2)/(n*m) }
# Regression test: both rank numerator and pair denominator exceed int32.
stopifnot(abs(auc(c(rep(0,50000),rep(1,50000)),1:100000)-1)<1e-12)
if(!exists("o",inherits=FALSE)) o<-read.csv(gzfile("results/oof/predictions.csv.gz"),colClasses=c(region="character",group="character",species="character",scheme="character",model="character"))
stopifnot(all(o$y %in% 0:1),all(is.finite(o$p)),all(o$p>=0 & o$p<=1),all(o$prevalence>=0 & o$prevalence<=1))
# Match every prediction across all three models before computing metrics.
id<-c("region","group","species","scheme","fold","row_id","y","prevalence")
stopifnot(!anyDuplicated(o[c(id,"model")]))
w<-reshape(o[c(id,"model","p")],idvar=id,timevar="model",direction="wide")
w<-w[complete.cases(w),]
cols<-c("p.linear","p.quadratic","p.xgb_shallow","prevalence")
# Pooled decomposition: exact same-species positive-negative pair counts.
calc<-function(d){
 n<-as.double(sum(d$y==1));m<-as.double(sum(d$y==0))
 z<-split(d,interaction(d$group,d$species,drop=TRUE))
 pairs<-vapply(z,function(q)as.double(sum(q$y==1))*as.double(sum(q$y==0)),numeric(1))
 do.call(rbind,lapply(cols,function(p){
 aa<-vapply(z,function(q)auc(q$y,q[[p]]),numeric(1))
 within<-sum(pairs*aa,na.rm=TRUE); total<-auc(d$y,d[[p]])*n*m
 data.frame(model=p,rows=nrow(d),species=length(z),macro=mean(aa,na.rm=TRUE),pooled=total/(n*m),within_pair_auc=within/sum(pairs),cross_pair_auc=if(n*m>sum(pairs))(total-within)/(n*m-sum(pairs))else NA,cross_pair_share=1-sum(pairs)/(n*m))
 }))
}
res<-do.call(rbind,lapply(split(w,interaction(w$region,w$scheme,drop=TRUE)),function(d)cbind(region=d$region[1],scheme=d$scheme[1],calc(d))))
stopifnot(all(is.finite(res$pooled)), all(abs(res$pooled-((1-res$cross_pair_share)*res$within_pair_auc+res$cross_pair_share*res$cross_pair_auc))<1e-10), all(res$cross_pair_auc>=-1e-10 & res$cross_pair_auc<=1+1e-10))
write.csv(res,"results/oof/pooling.csv",row.names=FALSE)
print(res,row.names=FALSE)
# Compare schemes only for identical records that received all three model scores in both.
key<-c("region","group","species","row_id","y")
v<-merge(w[w$scheme=="spatial",],w[w$scheme=="random",],by=key,suffixes=c(".spatial",".random"))
stopifnot(!anyDuplicated(v[key]))
sp<-do.call(rbind,lapply(split(v,interaction(v$region,v$group,v$species,drop=TRUE)),function(d){
 do.call(rbind,lapply(cols[1:3],function(p)data.frame(region=d$region[1],group=d$group[1],species=d$species[1],model=p,n=nrow(d),spatial=auc(d$y,d[[paste0(p,".spatial")]]),random=auc(d$y,d[[paste0(p,".random")]]))))
}))
write.csv(sp,"results/oof/matched_scheme.csv",row.names=FALSE)
print(aggregate(cbind(spatial,random)~region+model,sp,mean))
