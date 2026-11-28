# Aggregate corrected nested CV after all NZ shards complete.
dir.create("results/nested_revised",recursive=TRUE,showWarnings=FALSE)
regions<-c("AWT","CAN","NSW","NZ","SA","SWI")
paths<-c(file.path("inputs/main",paste0("nested-path-",c("AWT","NSW","SA","SWI"))),"inputs/can",file.path("inputs/nz",paste0("nested-nz-shard-",1:8)))
read_all<-function(f)do.call(rbind,lapply(paths,function(p){stopifnot(file.exists(file.path(p,f)));read.csv(file.path(p,f),stringsAsFactors=FALSE)}))
assigned<-do.call(rbind,lapply(1:8,function(i){a<-read.csv(file.path("inputs/nz",paste0("nested-nz-shard-",i),"assigned_taxa.csv"));stopifnot(all(a$shard==i));a}))
stopifnot(nrow(assigned)==52,!anyDuplicated(assigned$species))
audit<-read_all("eligibility.csv");outer<-read_all("outer_metrics.csv");inner<-read_all("inner_scores.csv")
keys<-c("region","group","species","fold")
stopifnot(setequal(unique(audit$region),regions),!anyDuplicated(audit[keys]),!anyDuplicated(outer[c(keys,"model")]),!anyDuplicated(inner[c(keys,"candidate")]))
stopifnot(setequal(assigned$species,unique(audit$species[audit$region=="NZ"])),nrow(audit[audit$region=="NZ",])==52*5)
b<-reshape(outer[c(keys,"model","auc")],idvar=keys,timevar="model",direction="wide")
models<-c("ridge_l_fixed","ridge_l_tuned","ridge_q_fixed","ridge_q_tuned","xgb_fixed","xgb_tuned")
ac<-paste0("auc.",models);stopifnot(all(ac%in%names(b)))
b<-b[complete.cases(b[ac]),]
sp<-aggregate(b[ac],b[keys[1:3]],mean)
for(m in c("ridge_l","ridge_q","xgb"))sp[[paste0("delta_",m)]]<-sp[[paste0("auc.",m,"_tuned")]]-sp[[paste0("auc.",m,"_fixed")]]
cols<-grep("^delta_",names(sp),value=TRUE)
rr<-aggregate(sp[c(ac,cols)],sp["region"],mean)
print(rr);print(table(sp$region));print(with(audit,table(region,status)));print(colMeans(rr[cols]));print(colMeans(sp[cols]))
set.seed(20261127)
ci<-do.call(rbind,lapply(cols,function(m){
 z<-split(sp[[m]],sp$region);bs<-sapply(z,function(x)replicate(4000,mean(sample(x,length(x),replace=TRUE))))
 v<-vapply(z,mean,0.0);q<-quantile(rowMeans(bs),c(.025,.975))
 data.frame(contrast=m,n=nrow(sp),mean=mean(v),lower=q[1],upper=q[2])
}))
print(ci)
# Compare the original and revised fixed/tuned scores on common cells. Check sample sizes.
old<-read.csv("inputs/old/outer_metrics.csv")
j<-merge(outer,old,by=c(keys,"model"),suffixes=c("_new","_old"))
stopifnot(!anyDuplicated(j[c(keys,"model")]),all(j$n_train_new==j$n_train_old),all(j$n_test_new==j$n_test_old))
j$delta_solver<-j$auc_new-j$auc_old
print(aggregate(delta_solver~region+model,j,mean))
for(n in c("assigned","audit","outer","inner","sp","rr","ci","j"))write.csv(get(n),file.path("results/nested_revised",paste0(n,".csv")),row.names=FALSE)
