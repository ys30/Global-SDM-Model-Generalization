# Compare selector variants on exactly the same evaluation cells.
dir.create("results/selection_variants",recursive=TRUE,showWarnings=FALSE)
d<-read.csv("inputs/fold_predictions.csv")
keys<-c("region","group","species","fold")
versions<-c("M1_size","M2_prevalence","M3_class_counts")
stopifnot(setequal(unique(d$version),versions),!anyDuplicated(d[c(keys,"version")]))
b<-reshape(d[c(keys,"version","auc_selected","auc_ridge_l")],idvar=keys,timevar="version",direction="wide")
stopifnot(all(complete.cases(b)),nrow(b)*3==nrow(d))
stopifnot(all(b$auc_ridge_l.M1_size==b$auc_ridge_l.M2_prevalence),all(b$auc_ridge_l.M1_size==b$auc_ridge_l.M3_class_counts))
b$M2_minus_M1<-b$auc_selected.M2_prevalence-b$auc_selected.M1_size
b$M3_minus_M1<-b$auc_selected.M3_class_counts-b$auc_selected.M1_size
for(v in versions)b[[paste0(v,"_minus_ridge")]]<-b[[paste0("auc_selected.",v)]]-b$auc_ridge_l.M1_size
metrics<-c("M2_minus_M1","M3_minus_M1",paste0(versions,"_minus_ridge"))
sp<-aggregate(b[metrics],b[keys[1:3]],mean)
rr<-aggregate(sp[metrics],sp["region"],mean)
set.seed(20261127)
ci<-do.call(rbind,lapply(metrics,function(m){
 z<-split(sp[[m]],sp$region);bs<-sapply(z,function(x)replicate(4000,mean(sample(x,length(x),replace=TRUE))))
 v<-vapply(z,mean,0.0);q<-quantile(rowMeans(bs),c(.025,.975))
 loo<-sapply(names(v),function(r)mean(v[names(v)!=r]))
 data.frame(contrast=m,n=nrow(sp),mean=mean(v),lower=q[1],upper=q[2],loo_min=min(loo),loo_max=max(loo),region_wins=sum(v>0),species_wins=sum(sp[[m]]>0))
}))
print(rr);print(ci)
write.csv(sp,"results/selection_variants/species.csv",row.names=FALSE)
write.csv(rr,"results/selection_variants/region.csv",row.names=FALSE)
write.csv(ci,"results/selection_variants/intervals.csv",row.names=FALSE)
