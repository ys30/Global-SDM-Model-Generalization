# Leave-one-region-out model selection using training-only covariates.
# Fold-level features are computed without held-out labels; evaluation uses matched model outcomes.
dir.create("results/selection_novelty",recursive=TRUE,showWarnings=FALSE)
regions<-c("AWT","CAN","NSW","NZ","SA","SWI")
models<-c("linear","quadratic","ridge_linear","ridge_quadratic","xgb_shallow","gam_low")
a<-read.csv("inputs/pa/fold_metrics.csv");stopifnot(setequal(unique(a$region),regions))
keys<-c("region","group","species","fold")
stopifnot(!anyDuplicated(a[c(keys,"model")]))
# n_train is outcome-independent apart from missing response availability. Do not use test prevalence or AUC as features.
b<-reshape(a[c(keys,"model","n_train","auc")],idvar=keys,timevar="model",direction="wide")
ac<-paste0("auc.",models);nc<-paste0("n_train.",models)
b<-b[complete.cases(b[c(ac,nc)]),]
stopifnot(all(apply(b[nc],1,function(x)length(unique(x)))==1))
b$log_n<-log(b[[nc[1]]])
features<-do.call(rbind,lapply(regions,function(r){
 d<-read.csv(file.path("inputs/novelty",paste0("novelty-",r),"features.csv"))
 stopifnot(all(d$region==r),!anyDuplicated(d[keys]));d
}))
b<-merge(b,features,by=keys,all.x=TRUE)
stopifnot(all(is.finite(b$train_presence)),all(b$n_train==b[[nc[1]]]),all(b$train_presence>0),all(b$train_presence<b$n_train))
b$logit_prev<-qlogis(b$train_presence/b$n_train)
b$log_pos<-log(b$train_presence)
b$log_neg<-log(b$n_train-b$train_presence)
b$log_geo<-log1p(b$geo_median)
b$log_env<-log1p(b$env_median)
stopifnot(nrow(b)==876,all(is.finite(b$log_geo)),all(is.finite(b$log_env)),all(is.finite(b$outside_fraction)))

# Fit simple ridge regression of paired model-minus-linear AUC on log training size.
# Balance regions and taxa; training-only centering; fixed penalty, no test-driven tuning.
out<-list();coefout<-list();predout<-list()
for(version in c("M1_size","M2_prevalence","M3_class_counts","M4_geo","M5_env","M6_combined"))for(r in regions){
 tr<-b[b$region!=r,];te<-b[b$region==r,]
 key<-interaction(tr$region,tr$group,tr$species,drop=TRUE)
 w<-1/as.numeric(table(key)[key])
 for(rr in unique(tr$region))w[tr$region==rr]<-w[tr$region==rr]/sum(w[tr$region==rr])
 w<-w/sum(w)
 vars<-switch(version,M1_size="log_n",M2_prevalence=c("log_n","logit_prev"),M3_class_counts=c("log_pos","log_neg"),M4_geo=c("log_n","logit_prev","log_geo"),M5_env=c("log_n","logit_prev","log_env","outside_fraction"),M6_combined=c("log_n","logit_prev","log_geo","log_env","outside_fraction"))
 t<-as.matrix(tr[vars]);v<-as.matrix(te[vars])
 mu<-colSums(w*t);sdv<-sqrt(colSums(w*sweep(t,2,mu)^2));stopifnot(all(sdv>0))
 x<-sweep(sweep(t,2,mu),2,sdv,"/");z<-sweep(sweep(v,2,mu),2,sdv,"/")
 for(k in seq_along(vars))z[,k]<-pmin(max(x[,k]),pmax(min(x[,k]),z[,k]))
 x<-cbind(1,x);z<-cbind(1,z)
 d<-as.matrix(tr[ac])-tr$auc.linear
 beta<-solve(crossprod(x,w*x)+diag(c(0,rep(.1,length(vars)))),crossprod(x,w*d))
 pred<-z%*%beta
 best_global<-which.max(colSums(w*as.matrix(tr[ac])))
 ix<-max.col(pred,ties.method="first")
 obs<-as.matrix(te[ac]);oracle<-apply(obs,1,max)
 for(j in seq_along(models))predout[[length(predout)+1]]<-data.frame(version=version,te[keys],model=models[j],pred_delta=pred[,j],obs_delta=obs[,j]-te$auc.linear)
 result<-data.frame(version=version,te[keys],auc_ridge_l=te$auc.ridge_linear,n_train=te[[nc[1]]],selected=models[ix],global_model=models[best_global],auc_selected=obs[cbind(seq_len(nrow(te)),ix)],auc_global=obs[,best_global],auc_linear=te$auc.linear,auc_oracle=oracle)
 out[[length(out)+1]]<-result
 coefout[[length(coefout)+1]]<-data.frame(version=version,held_out=r,model=rep(models,each=nrow(beta)),term=rep(c("intercept",vars),length(models)),coefficient=as.vector(beta))
}
d<-do.call(rbind,out)
d$gain_global<-d$auc_selected-d$auc_global;d$gain_linear<-d$auc_selected-d$auc_linear;d$regret<-d$auc_oracle-d$auc_selected
d$gain_ridge_l<-d$auc_selected-d$auc_ridge_l
metrics<-c("gain_ridge_l","auc_selected","auc_global","auc_linear","auc_oracle","gain_global","gain_linear","regret")
sp<-aggregate(d[metrics],d[c("version","region","group","species")],mean)
rr<-aggregate(sp[metrics],sp[c("version","region")],mean)
print(rr);print(aggregate(rr[metrics],rr["version"],mean));print(with(d,table(version,region,selected)))
write.csv(d,"results/selection_novelty/fold_predictions.csv",row.names=FALSE)
write.csv(sp,"results/selection_novelty/species.csv",row.names=FALSE)
write.csv(rr,"results/selection_novelty/region.csv",row.names=FALSE)
write.csv(do.call(rbind,coefout),"results/selection_novelty/coefficients.csv",row.names=FALSE)
write.csv(do.call(rbind,predout),"results/selection_novelty/predicted_deltas.csv",row.names=FALSE)
cat("Conditional matched benchmark; oracle uses held-out outcomes only as an unattainable reference. No independent external dataset claim.\n")
