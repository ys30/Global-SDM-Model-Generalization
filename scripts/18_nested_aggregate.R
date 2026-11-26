# Aggregate completed regional nested tuning artifacts.
dir.create("results/nested_combined",recursive=TRUE,showWarnings=FALSE)
regions<-c("AWT","CAN","NSW","NZ","SA","SWI")
read_region<-function(r,f){
 p<-file.path("inputs",paste0("nested-full-",r),f)
 stopifnot(file.exists(p))
 a<-read.csv(p,stringsAsFactors=FALSE);stopifnot(nrow(a)>0,all(a$region==r));a
}
sp<-do.call(rbind,lapply(regions,read_region,f="species_macro.csv"))
audit<-do.call(rbind,lapply(regions,read_region,f="eligibility.csv"))
outer<-do.call(rbind,lapply(regions,read_region,f="outer_metrics.csv"))
stopifnot(!anyDuplicated(sp[c("region","group","species")]),!anyDuplicated(outer[c("region","group","species","fold","model")]))
for(m in c("ridge_l","ridge_q","xgb"))sp[[paste0("delta_",m)]]<-sp[[paste0("auc.",m,"_tuned")]]-sp[[paste0("auc.",m,"_fixed")]]
cols<-grep("^delta_",names(sp),value=TRUE)
print(aggregate(sp[cols],sp["region"],mean))
print(table(sp$region));print(with(audit,table(region,status)))
print(aggregate(sp[cols],list(scope="species_weighted"),mean))
rr<-aggregate(sp[cols],sp["region"],mean);print(colMeans(rr[cols]))
write.csv(sp,"results/nested_combined/species.csv",row.names=FALSE)
write.csv(rr,"results/nested_combined/region.csv",row.names=FALSE)
write.csv(audit,"results/nested_combined/eligibility.csv",row.names=FALSE)
write.csv(outer,"results/nested_combined/outer_metrics.csv",row.names=FALSE)
