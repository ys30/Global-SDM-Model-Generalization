# Descriptive paired uncertainty: species resampling within region, not independent-site inference.
dir.create("results/uncertainty",recursive=TRUE,showWarnings=FALSE)
set.seed(20261126)
regions<-c("AWT","CAN","NSW","NZ","SA","SWI")
readpa<-function(r){d<-read.csv(file.path("inputs/pa",paste0("corrected-shard-",r),"species_macro.csv"));stopifnot(all(d$region==r));d}
pa<-do.call(rbind,lapply(regions,readpa))
ne<-read.csv("inputs/nested/species.csv")
bu<-read.csv("inputs/buffer/species.csv")
stopifnot(nrow(pa)==198,nrow(ne)==183,!anyDuplicated(pa[c("region","group","species")]),!anyDuplicated(ne[c("region","group","species")]))
contrasts<-list()
add<-function(d,label,v){stopifnot(length(v)==nrow(d));data.frame(region=d$region,group=d$group,species=d$species,contrast=label,delta=v)}
for(m in c("quadratic","xgb_shallow","ridge_linear","ridge_quadratic","gam_low"))contrasts[[length(contrasts)+1]]<-add(pa,paste0("PA_",m,"_minus_linear"),pa[[paste0("auc.",m)]]-pa$auc.linear)
for(m in c("ridge_l","ridge_q","xgb"))contrasts[[length(contrasts)+1]]<-add(ne,paste0("Nested_",m),ne[[paste0("delta_",m)]])
for(m in unique(bu$model))for(v in c("delta_size","delta_spatial","delta_total")){d<-bu[bu$model==m,];contrasts[[length(contrasts)+1]]<-add(d,paste("Buffer",m,v,sep="_"),d[[v]])}
a<-do.call(rbind,contrasts);B<-4000L
summary<-list();loo<-list()
for(n in unique(a$contrast)){
 d<-a[a$contrast==n,];stopifnot(all(is.finite(d$delta)),!anyDuplicated(d[c("region","group","species")]))
 z<-split(d$delta,d$region);stopifnot(length(z)==6)
 boots<-sapply(z,function(v)replicate(B,mean(sample(v,length(v),replace=TRUE))))
 est<-vapply(z,mean,0.0)
 for(r in names(z)){
 v<-z[[r]];ci<-quantile(boots[,r],c(.025,.975))
 summary[[length(summary)+1]]<-data.frame(contrast=n,scope=r,n=length(v),mean=mean(v),median=median(v),win=mean(v>0),lower=ci[1],upper=ci[2])
 }
 for(scope in c("region_weighted","species_weighted")){
 w<-if(scope=="region_weighted")rep(1/6,6)else lengths(z)/sum(lengths(z))
 ci<-quantile(as.vector(boots%*%w),c(.025,.975))
 summary[[length(summary)+1]]<-data.frame(contrast=n,scope=scope,n=nrow(d),mean=sum(est*w),median=NA_real_,win=NA_real_,lower=ci[1],upper=ci[2])
 }
 for(r in names(z))loo[[length(loo)+1]]<-data.frame(contrast=n,held_out=r,mean_remaining=mean(est[names(est)!=r]))
}
s<-do.call(rbind,summary);l<-do.call(rbind,loo)
write.csv(a,"results/uncertainty/paired_species.csv",row.names=FALSE)
write.csv(s,"results/uncertainty/descriptive_intervals.csv",row.names=FALSE)
write.csv(l,"results/uncertainty/leave_region_out_sensitivity.csv",row.names=FALSE)
print(s[s$scope=="region_weighted",],row.names=FALSE)
cat("Intervals conditional on observed regions and species resampling; shared sites and folds preclude confirmatory independence claims.\n")
