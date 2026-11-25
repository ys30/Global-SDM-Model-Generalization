b<-read.csv("results/pilot/paired_auc.csv")
b<-b[is.finite(b$delta),]
a<-aggregate(delta~region+group+species,b,mean)
set.seed(20261125)
summary<-do.call(rbind,lapply(split(a,a$region),function(d){
 x<-d$delta; boot<-replicate(2000,mean(sample(x,length(x),replace=TRUE)))
 data.frame(region=d$region[1],species=length(x),mean_delta=mean(x),median_delta=median(x),win_rate=mean(x>0),ci_low=unname(quantile(boot,.025)),ci_high=unname(quantile(boot,.975)))
}))
write.csv(summary,"results/pilot/species_macro_summary.csv",row.names=FALSE)
print(summary,row.names=FALSE)
m<-read.csv("results/pilot/glm_fold_metrics.csv")
k<-m[m$model=="linear",c("region","group","species","fold","n_train","n_test","test_presence")]
d<-merge(b,k,by=c("region","group","species","fold"))
d$sample_bin<-cut(d$n_train,c(0,100,500,2000,10000,Inf),include.lowest=TRUE)
print(aggregate(delta~sample_bin,d,mean))
write.csv(d,"results/pilot/paired_sample_diagnostics.csv",row.names=FALSE)
