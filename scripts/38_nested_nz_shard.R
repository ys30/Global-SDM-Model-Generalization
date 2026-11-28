library(disdat)
library(xgboost)
library(sf)
library(glmnet)
source('scripts/ridge_path.R')

set.seed(20261125)
dir.create("results/nested_nz_shard",recursive=TRUE,showWarnings=FALSE)
regions <- list(AWT=c("bird","plant"),CAN=NA_character_,NSW=c("ba","db","nb","ot","ou","rt","ru","sr"),NZ=NA_character_,SA=NA_character_,SWI=NA_character_)
auc <- function(y,p) { n<-as.double(sum(y==1)); m<-as.double(sum(y==0)); if(!n||!m)return(NA_real_); (sum(rank(p)[y==1])-n*(n+1)/2)/(n*m) }
cats<-list(AWT=character(),CAN="ontveg",NSW="vegsys",NZ=c("age","toxicats"),SA=character(),SWI="calc")
ap<-function(y,p){
 if(!sum(y==1))return(NA_real_)
 i<-order(p,decreasing=TRUE); y<-y[i];p<-p[i]
 end<-c(which(diff(p)!=0),length(p)); tp<-cumsum(y==1)[end]
 sum(diff(c(0,tp))/sum(y==1)*tp/end)
}

# Fixed pilot grid; outer test labels never enter tuning.
grid<-expand.grid(depth=c(2L,4L),rounds=c(100L,250L),eta=.05)
lambdas<-c(.001,.01,.1)
prep<-function(env,continuous,categorical,tr,te){
 supported<-rep(TRUE,length(te))
 for(p in categorical)supported<-supported & env[[p]][te]%in%unique(env[[p]][tr])
 te<-te[supported]
 x<-as.data.frame(env[continuous]);names(x)<-paste0("v",seq_along(continuous))
 mu<-vapply(x[tr,,drop=FALSE],mean,0.0);sdv<-vapply(x[tr,,drop=FALSE],sd,0.0)
 keep<-is.finite(sdv)&sdv>0
 x<-as.data.frame(scale(x[,keep,drop=FALSE],center=mu[keep],scale=sdv[keep]))
 num<-names(x)
 for(j in seq_along(categorical)){
 lev<-sort(unique(env[[categorical[j]]][tr]))
 if(length(lev)>1)for(h in seq_along(lev)[-1])x[[paste0("cat",j,"_",h)]]<-as.numeric(env[[categorical[j]]]==lev[h])
 }
 if(!ncol(x))return(NULL)
 list(x=as.matrix(x[tr,,drop=FALSE]),z=as.matrix(x[te,,drop=FALSE]),q=as.matrix(cbind(x,x[num]^2)[tr,,drop=FALSE]),qz=as.matrix(cbind(x,x[num]^2)[te,,drop=FALSE]),te=te)
}
safe<-function(f,n)tryCatch(suppressWarnings({p<-as.numeric(f());if(length(p)!=n||any(!is.finite(p)))stop("invalid predictions");p}),error=function(e)rep(NA_real_,n))
fitridge<-function(d,y,lambda,q=FALSE){
 x<-if(q)d$q else d$x;z<-if(q)d$qz else d$z
 f<-ridge_path(x,y,lambda)
 as.numeric(predict(f,z,s=lambda,type="response"))
}
fitxgb<-function(d,y,j){
 f<-xgb.train(params=list(objective="binary:logistic",max_depth=grid$depth[j],eta=grid$eta[j],nthread=1,seed=20261125),data=xgb.DMatrix(d$x,label=y),nrounds=grid$rounds[j],verbose=0)
 predict(f,xgb.DMatrix(d$z))
}
out<-list();tuning<-list();audit<-list()
target<-Sys.getenv("TARGET_REGION")
stopifnot(target %in% names(regions))
for(r in target)for(g in regions[[r]]){
 group<-if(is.na(g))NULL else g;label<-if(is.na(g))"all"else g
 pa<-disPa(r,group);env<-disEnv(r,group);stopifnot(identical(pa[,1:4],env[,1:4]))
 pred<-disPredictors(r);categorical<-intersect(cats[[r]],pred);continuous<-setdiff(pred,categorical)
 ok<-complete.cases(env[pred])&is.finite(pa$x)&is.finite(pa$y)
 for(p in continuous)ok<-ok&is.finite(env[[p]])
 pts<-st_as_sf(pa[ok,],coords=c("x","y"),crs=disCRS(r))
 ll<-st_coordinates(st_transform(pts,4326));lon<-atan2(mean(sin(ll[,1]*pi/180)),mean(cos(ll[,1]*pi/180)))*180/pi;lat<-mean(ll[,2])
 xy<-st_coordinates(st_transform(pts,sprintf("+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m",lat,lon)))
 u<-unique(xy);if(nrow(u)<5)next
 centers<-kmeans(u,centers=5,nstart=10)$centers
 fold<-rep(NA_integer_,nrow(pa));fold[ok]<-max.col(-sapply(1:5,function(j)rowSums(sweep(xy,2,centers[j,])^2)),ties.method="first")
 # Deterministic pilot selection independent of performance, up to 3 taxa per survey group.
 taxa<-sort(names(pa)[-(1:4)])
 taxa<-taxa[vapply(taxa,function(s){y<-pa[[s]][ok];sum(y==1,na.rm=TRUE)>=10&&sum(y==0,na.rm=TRUE)>=10},logical(1))]
 # Partition sorted eligible taxa only after the original NZ spatial folds are constructed.
 shard<-as.integer(Sys.getenv("SHARD"));stopifnot(r=="NZ",shard%in%1:8)
 taxa<-taxa[(seq_along(taxa)-1L)%%8L+1L==shard]
 write.csv(data.frame(species=taxa,shard=shard),"results/nested_nz_shard/assigned_taxa.csv",row.names=FALSE)
 for(s in taxa)for(k in 1:5){
 y<-pa[[s]];valid<-ok&!is.na(y)&y%in%0:1
 tr<-which(valid&fold!=k);te<-which(valid&fold==k)
 if(sum(y[tr]==1)<5||sum(y[tr]==0)<5||length(unique(y[te]))<2){audit[[length(audit)+1]]<-data.frame(region=r,group=label,species=s,fold=k,inner_n=0L,status="outer_ineligible");next}
 # Inner held-out blocks are the four remaining outer spatial blocks.
 scores<-list();inner_n<-0L
 for(j in setdiff(1:5,k)){
 it<-tr[fold[tr]!=j];iv<-tr[fold[tr]==j]
 stopifnot(!length(intersect(c(it,iv),te)))
 if(sum(y[it]==1)<5||sum(y[it]==0)<5||length(unique(y[iv]))<2)next
 d<-prep(env,continuous,categorical,it,iv)
 if(is.null(d)||length(unique(y[d$te]))<2)next
 inner_n<-inner_n+1L
 z<-c(vapply(lambdas,function(l)auc(y[d$te],safe(function()fitridge(d,y[it],l),length(d$te))),0.0),
 vapply(lambdas,function(l)auc(y[d$te],safe(function()fitridge(d,y[it],l,TRUE),length(d$te))),0.0),
 vapply(seq_len(nrow(grid)),function(h)auc(y[d$te],safe(function()fitxgb(d,y[it],h),length(d$te))),0.0))
 scores[[length(scores)+1]]<-z
 }
 if(inner_n<2){audit[[length(audit)+1]]<-data.frame(region=r,group=label,species=s,fold=k,inner_n=inner_n,status="insufficient_inner_folds");next}
 m<-do.call(rbind,scores);common<-complete.cases(m)
 if(sum(common)<2){audit[[length(audit)+1]]<-data.frame(region=r,group=label,species=s,fold=k,inner_n=inner_n,status="candidate_failure");next}
 avg<-colMeans(m[common,,drop=FALSE]);a<-which.max(avg[1:3]);b<-which.max(avg[4:6]);c<-which.max(avg[7:10])
 tuning[[length(tuning)+1]]<-data.frame(region=r,group=label,species=s,fold=k,inner_n=sum(common),candidate=1:10,inner_auc=avg)
 d<-prep(env,continuous,categorical,tr,te)
 if(is.null(d)||length(unique(y[d$te]))<2){audit[[length(audit)+1]]<-data.frame(region=r,group=label,species=s,fold=k,inner_n=sum(common),status="outer_unsupported");next}
 fs<-list(ridge_l_fixed=function()fitridge(d,y[tr],.01),ridge_l_tuned=function()fitridge(d,y[tr],lambdas[a]),ridge_q_fixed=function()fitridge(d,y[tr],.01,TRUE),ridge_q_tuned=function()fitridge(d,y[tr],lambdas[b],TRUE),xgb_fixed=function()fitxgb(d,y[tr],1),xgb_tuned=function()fitxgb(d,y[tr],c))
 for(model in names(fs)){
 p<-safe(fs[[model]],length(d$te))
 out[[length(out)+1]]<-data.frame(region=r,group=label,species=s,fold=k,model=model,n_train=length(tr),n_test=length(d$te),auc=auc(y[d$te],p),lambda_l=lambdas[a],lambda_q=lambdas[b],xgb_depth=grid$depth[c],xgb_rounds=grid$rounds[c])
 }
 audit[[length(audit)+1]]<-data.frame(region=r,group=label,species=s,fold=k,inner_n=sum(common),status="evaluated")
 }
 cat("Completed",r,label,"\n")
 if(length(out))write.csv(do.call(rbind,out),"results/nested_nz_shard/outer_metrics.csv",row.names=FALSE)
 if(length(tuning))write.csv(do.call(rbind,tuning),"results/nested_nz_shard/inner_scores.csv",row.names=FALSE)
 if(length(audit))write.csv(do.call(rbind,audit),"results/nested_nz_shard/eligibility.csv",row.names=FALSE)
}
a<-do.call(rbind,out);stopifnot(nrow(a)>0)
write.csv(a,"results/nested_nz_shard/outer_metrics.csv",row.names=FALSE)
write.csv(do.call(rbind,tuning),"results/nested_nz_shard/inner_scores.csv",row.names=FALSE)
write.csv(do.call(rbind,audit),"results/nested_nz_shard/eligibility.csv",row.names=FALSE)
b<-reshape(a[c("region","group","species","fold","model","auc")],idvar=c("region","group","species","fold"),timevar="model",direction="wide")
b<-b[complete.cases(b),];cols<-grep("^auc\\.",names(b),value=TRUE)
stopifnot(nrow(b)>0)
sp<-aggregate(b[cols],b[c("region","group","species")],mean)
print(aggregate(sp[cols],sp["region"],mean));print(table(sp$region))
write.csv(sp,"results/nested_nz_shard/species_macro.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"results/nested_nz_shard/session_info.txt")
