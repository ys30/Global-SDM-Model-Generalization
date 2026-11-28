# Descriptive within-region associations, not causal effects.
dir.create("results/novelty_mechanism",recursive=TRUE,showWarnings=FALSE)
regions<-c("AWT","CAN","NSW","NZ","SA","SWI")
keys<-c("region","group","species","fold")
a<-read.csv("inputs/pa/fold_metrics.csv")
f<-do.call(rbind,lapply(regions,function(r)read.csv(file.path("inputs/novelty",paste0("novelty-",r),"features.csv"))))
models<-c("linear","quadratic","ridge_linear","ridge_quadratic","xgb_shallow","gam_low")
ac<-paste0("auc.",models)
stopifnot(!anyDuplicated(a[c(keys,"model")]),!anyDuplicated(f[keys]))
b<-reshape(a[c(keys,"model","auc","n_train","n_test")],idvar=keys,timevar="model",direction="wide")
b<-b[complete.cases(b[ac]),]
b<-merge(b,f,by=keys,all.x=TRUE)
stopifnot(nrow(b)==876,all(b$n_train==b$n_train.linear),all(b$n_test==b$n_test.linear))
b$log_n<-log(b$n_train)
b$logit_prev<-qlogis(b$train_presence/b$n_train)
b$log_geo<-log1p(b$geo_median)
b$log_env<-log1p(b$env_median)
vars<-c("log_n","logit_prev","log_geo","log_env","outside_fraction")
stopifnot(all(is.finite(as.matrix(b[vars]))))
# Within-region weighted centering removes between-region mean differences.
# Equal region and species weights, equal fold weights within species.
b$id<-interaction(b$region,b$group,b$species,drop=TRUE)
b$w<-1/as.numeric(table(b$id)[b$id])
for(r in regions)b$w[b$region==r]<-b$w[b$region==r]/sum(b$w[b$region==r])
b$w<-b$w/sum(b$w)
for(v in vars){
 for(r in regions){i<-b$region==r;b[i,v]<-b[i,v]-weighted.mean(b[i,v],b$w[i])}
 b[[v]]<-b[[v]]/sqrt(sum(b$w*b[[v]]^2))
}
out<-list();diag<-list()
for(m in models[-1]){
 b$delta<-b[[paste0("auc.",m)]]-b$auc.linear
 for(drop in c("none",regions)){
 d<-if(drop=="none")b else b[b$region!=drop,]
 fit<-lm(delta~factor(region)+log_n+logit_prev+log_geo+log_env+outside_fraction,data=d,weights=w)
 cf<-coef(fit)[vars]
 out[[length(out)+1]]<-data.frame(model=m,excluded=drop,term=vars,coefficient=unname(cf))
 if(drop=="none"){
  base<-lm(delta~factor(region)+log_n+logit_prev,data=d,weights=w)
  diag[[length(diag)+1]]<-data.frame(model=m,r2_base=summary(base)$r.squared,r2_full=summary(fit)$r.squared,partial_r2=1-sum(d$w*resid(fit)^2)/sum(d$w*resid(base)^2),condition=kappa(model.matrix(fit)))
 }
 }
}
o<-do.call(rbind,out);q<-do.call(rbind,diag)
print(o[o$excluded=="none",]);print(q)
write.csv(o,"results/novelty_mechanism/coefficients.csv",row.names=FALSE)
write.csv(q,"results/novelty_mechanism/fit_diagnostics.csv",row.names=FALSE)
write.csv(b,"results/novelty_mechanism/matched_features.csv",row.names=FALSE)
