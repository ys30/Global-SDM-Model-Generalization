library(disdat)
library(xgboost)
library(sf)
set.seed(20261125)
dir.create("results/spatial_ml",recursive=TRUE,showWarnings=FALSE)
regions <- list(AWT=c("bird","plant"),CAN=NA_character_,NSW=c("ba","db","nb","ot","ou","rt","ru","sr"),NZ=NA_character_,SA=NA_character_,SWI=NA_character_)
auc <- function(y,p) { n<-sum(y==1); m<-sum(y==0); if(!n||!m)return(NA_real_); (sum(rank(p)[y==1])-n*(n+1)/2)/(n*m) }
out<-list()
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
 out[[length(out)+1L]]<-data.frame(region=r,group=label,species=s,fold=k,model=model,n_train=length(tr),n_test=length(te),test_presence=sum(y[te]==1),auc=if(success)auc(y[te],p)else NA_real_,brier=if(success)mean((y[te]-p)^2)else NA_real_,success=success,warnings=paste(unique(warnings),collapse=";"))
 }
 }
 }
}
a<-do.call(rbind,out);write.csv(a,"results/spatial_ml/glm_fold_metrics.csv",row.names=FALSE)
b<-reshape(a[,c("region","group","species","fold","model","auc")],idvar=c("region","group","species","fold"),timevar="model",direction="wide")
b$delta<-b$auc.quadratic-b$auc.linear
b$delta_xgb<-b$auc.xgb_shallow-b$auc.linear
c<-b[complete.cases(b[,c("auc.linear","auc.quadratic","auc.xgb_shallow")]),]
sp<-aggregate(cbind(auc.linear,auc.quadratic,auc.xgb_shallow)~region+group+species,c,mean)
print(aggregate(cbind(auc.linear,auc.quadratic,auc.xgb_shallow)~region,sp,mean))
print(table(sp$region))
write.csv(sp,"results/spatial_ml/species_macro.csv",row.names=FALSE)
write.csv(b,"results/spatial_ml/paired_auc.csv",row.names=FALSE)
print(with(a,table(region,model,success)))
print(aggregate(delta~region,b,mean))
