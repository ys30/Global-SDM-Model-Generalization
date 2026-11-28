# Paired novelty selector inference and held-region delta prediction diagnostics.
dir.create("results/novelty_summary",recursive=TRUE,showWarnings=FALSE)
d<-read.csv("inputs/fold_predictions.csv")
p<-read.csv("inputs/predicted_deltas.csv")
keys<-c("region","group","species","fold")
versions<-c("M1_size","M2_prevalence","M3_class_counts","M4_geo","M5_env","M6_combined")
stopifnot(setequal(unique(d$version),versions),!anyDuplicated(d[c(keys,"version")]),!anyDuplicated(p[c(keys,"version","model")]))
b<-reshape(d[c(keys,"version","auc_selected","auc_ridge_l")],idvar=keys,timevar="version",direction="wide")
stopifnot(nrow(b)==876,all(complete.cases(b)),nrow(d)==876*6)
metrics<-character()
for(v in versions){
 stopifnot(all(b[[paste0("auc_ridge_l.",v)]]==b$auc_ridge_l.M1_size))
 n<-paste0(v,"_minus_ridge");b[[n]]<-b[[paste0("auc_selected.",v)]]-b$auc_ridge_l.M1_size;metrics<-c(metrics,n)
}
for(v in versions[4:6])for(ref in versions[1:2]){
 n<-paste0(v,"_minus_",ref);b[[n]]<-b[[paste0("auc_selected.",v)]]-b[[paste0("auc_selected.",ref)]];metrics<-c(metrics,n)
}
sp<-aggregate(b[metrics],b[keys[1:3]],mean);stopifnot(nrow(sp)==198)
rr<-aggregate(sp[metrics],sp["region"],mean)
set.seed(20261128)
ci<-do.call(rbind,lapply(metrics,function(m){
 z<-split(sp[[m]],sp$region);bs<-sapply(z,function(x)replicate(4000,mean(sample(x,length(x),replace=TRUE))))
 v<-vapply(z,mean,0.0);q<-quantile(rowMeans(bs),c(.025,.975))
 loo<-sapply(names(v),function(r)mean(v[names(v)!=r]))
 data.frame(contrast=m,n=nrow(sp),mean=mean(v),lower=q[1],upper=q[2],loo_min=min(loo),loo_max=max(loo),region_wins=sum(v>0),species_wins=sum(sp[[m]]>0))
}))
print(ci)
# Evaluate model-minus-linear prediction errors, excluding tautological linear reference.
p<-p[p$model!="linear",];stopifnot(nrow(p)==876*6*5,all(is.finite(p$pred_delta)),all(is.finite(p$obs_delta)))
p$se<-(p$pred_delta-p$obs_delta)^2;p$ae<-abs(p$pred_delta-p$obs_delta)
p$bias<-p$pred_delta-p$obs_delta
pk<-c("version","region","group","species","model")
ps<-aggregate(p[c("se","ae","bias")],p[pk],mean)
pr<-aggregate(ps[c("se","ae","bias")],ps[c("version","region","model")],mean)
pg<-aggregate(pr[c("se","ae","bias")],pr["version"],mean);pg$rmse<-sqrt(pg$se)
print(pg)
write.csv(ci,"results/novelty_summary/contrasts.csv",row.names=FALSE)
write.csv(rr,"results/novelty_summary/region.csv",row.names=FALSE)
write.csv(pr,"results/novelty_summary/prediction_region_model.csv",row.names=FALSE)
write.csv(pg,"results/novelty_summary/prediction_global.csv",row.names=FALSE)
