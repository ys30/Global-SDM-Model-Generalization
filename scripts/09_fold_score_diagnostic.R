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

# Separate within-fold discrimination from between-fold score shifts.
d<-do.call(rbind,lapply(split(w,interaction(w$region,w$scheme,w$group,w$species,w$fold,drop=TRUE)),function(q){
 if(length(unique(q$y))<2)return(NULL)
 do.call(rbind,lapply(cols,function(p)data.frame(region=q$region[1],scheme=q$scheme[1],group=q$group[1],species=q$species[1],fold=q$fold[1],model=p,auc=auc(q$y,q[[p]]),n=nrow(q))))
}))
stopifnot(all(abs(d$auc[d$model=="prevalence"]-.5)<1e-12))
sp<-aggregate(auc~region+scheme+group+species+model,d,mean)
print(aggregate(auc~region+scheme+model,sp,mean))
write.csv(d,"results/oof/within_fold_auc.csv",row.names=FALSE)
# Rank transform within each species and fold: diagnostic only, not deployment scores.
for(p in cols[1:3])w[[paste0(p,".rank")]]<-ave(w[[p]],interaction(w$region,w$scheme,w$group,w$species,w$fold,drop=TRUE),FUN=function(x)(rank(x,ties.method="average")-.5)/length(x))
a<-do.call(rbind,lapply(split(w,interaction(w$region,w$scheme,drop=TRUE)),function(q)do.call(rbind,lapply(cols[1:3],function(p)data.frame(region=q$region[1],scheme=q$scheme[1],model=p,raw=auc(q$y,q[[p]]),rank_normalized=auc(q$y,q[[paste0(p,".rank")]]))))))
print(a,row.names=FALSE)
write.csv(a,"results/oof/rank_diagnostic.csv",row.names=FALSE)
