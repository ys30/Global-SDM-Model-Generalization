# Leave-one-region-out model selection using training-only covariates.
# Fold-level features are computed without held-out labels; evaluation uses matched model outcomes.
dir.create("results/selection_revised",recursive=TRUE,showWarnings=FALSE)
regions<-c("AWT","CAN","NSW","NZ","SA","SWI")
models<-c("linear","quadratic","ridge_linear","ridge_quadratic","xgb_shallow","gam_low")
a<-read.csv("inputs/fold_metrics.csv");stopifnot(setequal(unique(a$region),regions))
keys<-c("region","group","species","fold")
stopifnot(!anyDuplicated(a[c(keys,"model")]))
# n_train is outcome-independent apart from missing response availability. Do not use test prevalence or AUC as features.
b<-reshape(a[c(keys,"model","n_train","auc")],idvar=keys,timevar="model",direction="wide")
ac<-paste0("auc.",models);nc<-paste0("n_train.",models)
b<-b[complete.cases(b[c(ac,nc)]),]
stopifnot(all(apply(b[nc],1,function(x)length(unique(x)))==1))
b$log_n<-log(b[[nc[1]]])
# Fit simple ridge regression of paired model-minus-linear AUC on log training size.
# Balance regions and taxa; training-only centering; fixed penalty, no test-driven tuning.
out<-list();coefout<-list()
for(r in regions){
 tr<-b[b$region!=r,];te<-b[b$region==r,]
 key<-interaction(tr$region,tr$group,tr$species,drop=TRUE)
 w<-1/as.numeric(table(key)[key])
 for(rr in unique(tr$region))w[tr$region==rr]<-w[tr$region==rr]/sum(w[tr$region==rr])
 w<-w/sum(w)
 mu<-sum(w*tr$log_n);sdv<-sqrt(sum(w*(tr$log_n-mu)^2));stopifnot(sdv>0)
 x<-cbind(1,(tr$log_n-mu)/sdv)
 z<-cbind(1,(te$log_n-mu)/sdv)
 # Clamp to training feature range, to avoid unconstrained meta-model extrapolation.
 z[,2]<-pmin(max(x[,2]),pmax(min(x[,2]),z[,2]))
 d<-as.matrix(tr[ac])-tr$auc.linear
 beta<-solve(crossprod(x,w*x)+diag(c(0,.1)),crossprod(x,w*d))
 pred<-z%*%beta
 best_global<-which.max(colSums(w*as.matrix(tr[ac])))
 ix<-max.col(pred,ties.method="first")
 obs<-as.matrix(te[ac]);oracle<-apply(obs,1,max)
 result<-data.frame(te[keys],n_train=te[[nc[1]]],selected=models[ix],global_model=models[best_global],auc_selected=obs[cbind(seq_len(nrow(te)),ix)],auc_global=obs[,best_global],auc_linear=te$auc.linear,auc_oracle=oracle)
 out[[length(out)+1]]<-result
 coefout[[length(coefout)+1]]<-data.frame(held_out=r,model=models,intercept=beta[1,],slope=beta[2,],center=mu,scale=sdv)
}
d<-do.call(rbind,out)
d$gain_global<-d$auc_selected-d$auc_global;d$gain_linear<-d$auc_selected-d$auc_linear;d$regret<-d$auc_oracle-d$auc_selected
metrics<-c("auc_selected","auc_global","auc_linear","auc_oracle","gain_global","gain_linear","regret")
sp<-aggregate(d[metrics],d[c("region","group","species")],mean)
rr<-aggregate(sp[metrics],sp["region"],mean)
print(rr);print(colMeans(rr[metrics]));print(with(d,table(region,selected)))
write.csv(d,"results/selection_revised/fold_predictions.csv",row.names=FALSE)
write.csv(sp,"results/selection_revised/species.csv",row.names=FALSE)
write.csv(rr,"results/selection_revised/region.csv",row.names=FALSE)
write.csv(do.call(rbind,coefout),"results/selection_revised/coefficients.csv",row.names=FALSE)
cat("Conditional matched benchmark; oracle uses held-out outcomes only as an unattainable reference. No independent external dataset claim.\n")
