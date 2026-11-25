library(disdat)
library(xgboost)
library(sf)
set.seed(20261125)
dir.create("results/oof",recursive=TRUE,showWarnings=FALSE)
regions <- list(AWT=c("bird","plant"),CAN=NA_character_,NSW=c("ba","db","nb","ot","ou","rt","ru","sr"),NZ=NA_character_,SA=NA_character_,SWI=NA_character_)
auc <- function(y,p) { n<-sum(y==1); m<-sum(y==0); if(!n||!m)return(NA_real_); (sum(rank(p)[y==1])-n*(n+1)/2)/(n*m) }
out<-list(); oof<-list()
for(r in names(regions)) for(g in regions[[r]]){
 group<-if(is.na(g))NULL else g; label<-if(is.na(g))"all" else g
 pa<-disPa(r,group); env<-disEnv(r,group)
 stopifnot(identical(pa[,1:4],env[,1:4]))
 pred<-intersect(disPredictors(r),names(env)); stopifnot(length(pred)>0)
 # Pilot deliberately uses only numeric predictors; categorical metadata review pending.
 pred<-pred[vapply(env[pred],is.numeric,logical(1))]
 pred<-setdiff(pred,c("vegsys"))
 stopifnot(length(pred)>0)
 ok<-complete.cases(env[pred]) & is.finite(pa$x) & is.finite(pa$y)
 for(p in pred)ok<-ok & is.finite(env[[p]])
 # Project each survey to local azimuthal equidistant coordinates in metres.
 pts<-st_as_sf(pa[ok,],coords=c("x","y"),crs=disCRS(r))
 ll<-st_coordinates(st_transform(pts,4326))
 lon<-atan2(mean(sin(ll[,1]*pi/180)),mean(cos(ll[,1]*pi/180)))*180/pi
 lat<-mean(ll[,2])
 xy<-st_coordinates(st_transform(pts,sprintf("+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m",lat,lon)))
 # Five contiguous Voronoi regions; identical locations always share a fold.
 u<-unique(xy); if(nrow(u)<5)next
 centers<-kmeans(u,centers=5,nstart=10)$centers
 distmat<-sapply(1:5,function(j)rowSums(sweep(xy,2,centers[j,])^2))
 fold<-rep(NA_integer_,nrow(pa));fold[ok]<-max.col(-distmat,ties.method="first")
 spatial_fold<-fold
 # Random assignment by identical raw coordinate pairs, never individual repeat visits.
 ids<-paste(pa$x[ok],pa$y[ok],sep="|")
 uid<-unique(ids)
 random_fold<-rep(NA_integer_,nrow(pa))
 random_fold[ok]<-sample(rep(1:5,length.out=length(uid)))[match(ids,uid)]
 for(scheme in c("spatial","random")){
 fold<-if(scheme=="spatial")spatial_fold else random_fold
 for(s in names(pa)[-(1:4)]){
 y<-pa[[s]]; valid<-ok & !is.na(y) & y %in% c(0,1)
 if(sum(y[valid]==1)<10||sum(y[valid]==0)<10)next
 for(k in 1:5){
 tr<-which(valid & fold!=k); te<-which(valid & fold==k)
 if(sum(y[tr]==1)<5||sum(y[tr]==0)<5||length(unique(y[te]))<2)next
 x<-as.data.frame(env[pred]); names(x)<-paste0("v",seq_along(pred))
 mu<-vapply(x[tr,,drop=FALSE],mean,0.0); sdv<-vapply(x[tr,,drop=FALSE],sd,0.0)
 keep<-is.finite(sdv)&sdv>0
 if(!any(keep))next
 x<-as.data.frame(scale(x[,keep,drop=FALSE],center=mu[keep],scale=sdv[keep]))
 for(model in c("linear","quadratic","xgb_shallow")){
 d<-x
 if(model=="quadratic")for(p in names(x))d[[paste0(p,"_sq")]]<-x[[p]]^2
 d$y<-y
 warnings<-character()
 fit<-tryCatch(withCallingHandlers(if(model=="xgb_shallow") xgb.train(params=list(objective="binary:logistic",max_depth=2,eta=.05,nthread=1,subsample=1,colsample_bytree=1),data=xgb.DMatrix(as.matrix(x[tr,,drop=FALSE]),label=y[tr]),nrounds=100,verbose=0) else glm(y~.,data=d[tr,,drop=FALSE],family=binomial(),control=glm.control(maxit=50)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
 p<-if(inherits(fit,"error"))rep(NA_real_,length(te)) else tryCatch(suppressWarnings(if(model=="xgb_shallow")predict(fit,xgb.DMatrix(as.matrix(x[te,,drop=FALSE]))) else predict(fit,newdata=d[te,,drop=FALSE],type="response")),error=function(e)rep(NA_real_,length(te)))
 success<-!inherits(fit,"error") && (model=="xgb_shallow" || isTRUE(fit$converged)) && all(is.finite(p))
 if(success) oof[[length(oof)+1L]]<-data.frame(region=r,group=label,species=s,scheme=scheme,fold=k,row_id=te,y=y[te],model=model,p=as.numeric(p),prevalence=mean(y[tr]))
 out[[length(out)+1L]]<-data.frame(scheme=scheme,region=r,group=label,species=s,fold=k,model=model,n_train=length(tr),n_test=length(te),test_presence=sum(y[te]==1),auc=if(success)auc(y[te],p)else NA_real_,brier=if(success)mean((y[te]-p)^2)else NA_real_,success=success,warnings=paste(unique(warnings),collapse=";"))
 }
 }
 }
}
}
a<-do.call(rbind,out);write.csv(a,"results/oof/glm_fold_metrics.csv",row.names=FALSE)

o<-do.call(rbind,oof)
write.csv(o,gzfile("results/oof/predictions.csv.gz"),row.names=FALSE)
# Match every prediction across all three models before computing metrics.
id<-c("region","group","species","scheme","fold","row_id","y","prevalence")
stopifnot(!anyDuplicated(o[c(id,"model")]))
w<-reshape(o[c(id,"model","p")],idvar=id,timevar="model",direction="wide")
w<-w[complete.cases(w),]
cols<-c("p.linear","p.quadratic","p.xgb_shallow","prevalence")
# Pooled decomposition: exact same-species positive-negative pair counts.
calc<-function(d){
 n<-sum(d$y==1);m<-sum(d$y==0)
 z<-split(d,interaction(d$group,d$species,drop=TRUE))
 pairs<-vapply(z,function(q)sum(q$y==1)*sum(q$y==0),numeric(1))
 do.call(rbind,lapply(cols,function(p){
 aa<-vapply(z,function(q)auc(q$y,q[[p]]),numeric(1))
 within<-sum(pairs*aa,na.rm=TRUE); total<-auc(d$y,d[[p]])*n*m
 data.frame(model=p,rows=nrow(d),species=length(z),macro=mean(aa,na.rm=TRUE),pooled=total/(n*m),within_pair_auc=within/sum(pairs),cross_pair_auc=if(n*m>sum(pairs))(total-within)/(n*m-sum(pairs))else NA,cross_pair_share=1-sum(pairs)/(n*m))
 }))
}
res<-do.call(rbind,lapply(split(w,interaction(w$region,w$scheme,drop=TRUE)),function(d)cbind(region=d$region[1],scheme=d$scheme[1],calc(d))))
write.csv(res,"results/oof/pooling.csv",row.names=FALSE)
print(res,row.names=FALSE)
# Compare schemes only for identical records that received all three model scores in both.
key<-c("region","group","species","row_id","y")
v<-merge(w[w$scheme=="spatial",],w[w$scheme=="random",],by=key,suffixes=c(".spatial",".random"))
stopifnot(!anyDuplicated(v[key]))
sp<-do.call(rbind,lapply(split(v,interaction(v$region,v$group,v$species,drop=TRUE)),function(d){
 do.call(rbind,lapply(cols[1:3],function(p)data.frame(region=d$region[1],group=d$group[1],species=d$species[1],model=p,n=nrow(d),spatial=auc(d$y,d[[paste0(p,".spatial")]]),random=auc(d$y,d[[paste0(p,".random")]]))))
}))
write.csv(sp,"results/oof/matched_scheme.csv",row.names=FALSE)
print(aggregate(cbind(spatial,random)~region+model,sp,mean))
