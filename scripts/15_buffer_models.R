library(disdat)
library(xgboost)
library(sf)
library(glmnet)
library(RANN)
set.seed(20261125)
dir.create("results/buffer_models",recursive=TRUE,showWarnings=FALSE)
regions <- list(AWT=c("bird","plant"),CAN=NA_character_,NSW=c("ba","db","nb","ot","ou","rt","ru","sr"),NZ=NA_character_,SA=NA_character_,SWI=NA_character_)
auc <- function(y,p) { n<-as.double(sum(y==1)); m<-as.double(sum(y==0)); if(!n||!m)return(NA_real_); (sum(rank(p)[y==1])-n*(n+1)/2)/(n*m) }
cats<-list(AWT=character(),CAN="ontveg",NSW="vegsys",NZ=c("age","toxicats"),SA=character(),SWI="calc")
ap<-function(y,p){
 if(!sum(y==1))return(NA_real_)
 i<-order(p,decreasing=TRUE); y<-y[i];p<-p[i]
 end<-c(which(diff(p)!=0),length(p)); tp<-cumsum(y==1)[end]
 sum(diff(c(0,tp))/sum(y==1)*tp/end)
}
out<-list(); coverage<-list()
for(r in names(regions)) for(g in regions[[r]]){
 group<-if(is.na(g))NULL else g; label<-if(is.na(g))"all" else g
 pa<-disPa(r,group); env<-disEnv(r,group)
 stopifnot(identical(pa[,1:4],env[,1:4]))
 pred<-intersect(disPredictors(r),names(env)); stopifnot(length(pred)>0)
 categorical<-intersect(cats[[r]],pred)
 continuous<-setdiff(pred,categorical)
 stopifnot(all(vapply(env[continuous],is.numeric,logical(1))))
 ok<-complete.cases(env[pred]) & is.finite(pa$x) & is.finite(pa$y)
 for(p in continuous)ok<-ok & is.finite(env[[p]])
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
 coords<-matrix(NA_real_,nrow(pa),2);coords[ok,]<-xy
 for(s in names(pa)[-(1:4)]){
 y<-pa[[s]]; valid<-ok & !is.na(y) & y %in% c(0,1)
 if(sum(y[valid]==1)<10||sum(y[valid]==0)<10)next
 for(k in 1:5){
 tr<-which(valid & fold!=k); te<-which(valid & fold==k)
 if(sum(y[tr]==1)<5||sum(y[tr]==0)<5||length(unique(y[te]))<2)next
 base_tr<-tr; base_te<-te
 # Fix test records across buffer levels using the largest buffer's category support.
 dist<-nn2(coords[base_te,,drop=FALSE],coords[base_tr,,drop=FALSE],k=1)$nn.dists[,1]/1000
 maxtr<-base_tr[dist>100]
 if(sum(y[maxtr]==1)<5||sum(y[maxtr]==0)<5)next
 supported_all<-rep(TRUE,length(base_te))
 for(p in categorical)supported_all<-supported_all & env[[p]][base_te] %in% unique(env[[p]][maxtr])
 base_te<-base_te[supported_all]
 if(length(unique(y[base_te]))<2)next
 for(buffer in c(0,25,50,100)){
 tr<-base_tr[dist>buffer];te<-base_te
 stopifnot(min(dist[dist>buffer])>buffer)
 # Fit all transformations on training records only. Unsupported test levels
 # are excluded identically for all models and explicitly counted.
 old_te<-length(te)
 supported<-rep(TRUE,length(te))
 for(p in categorical)supported<-supported & env[[p]][te] %in% unique(env[[p]][tr])
 te<-te[supported]
 coverage[[length(coverage)+1L]]<-data.frame(region=r,group=label,species=s,fold=k,buffer_km=buffer,test_before=old_te,test_after=length(te),unseen_category=old_te-length(te))
 if(!length(te)||length(unique(y[te]))<2)next
 x<-as.data.frame(env[continuous]); names(x)<-paste0("v",seq_along(continuous))
 mu<-vapply(x[tr,,drop=FALSE],mean,0.0); sdv<-vapply(x[tr,,drop=FALSE],sd,0.0)
 keep<-is.finite(sdv)&sdv>0
 x<-as.data.frame(scale(x[,keep,drop=FALSE],center=mu[keep],scale=sdv[keep]))
 numeric_names<-names(x)
 for(j in seq_along(categorical)){
 p<-categorical[j]; lev<-sort(unique(env[[p]][tr]))
 if(length(lev)>1)for(h in seq_along(lev)[-1])x[[paste0("cat",j,"_",h)]]<-as.numeric(env[[p]]==lev[h])
 }
 if(!ncol(x))next
 for(model in c("linear","quadratic","xgb_shallow","ridge_linear","ridge_quadratic")){
 d<-x
 if(model %in% c("quadratic","ridge_quadratic"))for(p in numeric_names)d[[paste0(p,"_sq")]]<-x[[p]]^2
 d$y<-y
 warnings<-character()
 fit<-tryCatch(withCallingHandlers(if(startsWith(model,"ridge_")) glmnet(as.matrix(d[tr,setdiff(names(d),"y"),drop=FALSE]),y[tr],family="binomial",alpha=0,lambda=.01,standardize=TRUE) else if(model=="xgb_shallow") xgb.train(params=list(objective="binary:logistic",max_depth=2,eta=.05,nthread=1,subsample=1,colsample_bytree=1),data=xgb.DMatrix(as.matrix(x[tr,,drop=FALSE]),label=y[tr]),nrounds=100,verbose=0) else glm(y~.,data=d[tr,,drop=FALSE],family=binomial(),control=glm.control(maxit=50)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
 p<-if(inherits(fit,"error"))rep(NA_real_,length(te)) else tryCatch(suppressWarnings(if(startsWith(model,"ridge_")) as.numeric(predict(fit,as.matrix(d[te,setdiff(names(d),"y"),drop=FALSE]),s=.01,type="response")) else if(model=="xgb_shallow")predict(fit,xgb.DMatrix(as.matrix(x[te,,drop=FALSE]))) else predict(fit,newdata=d[te,,drop=FALSE],type="response")),error=function(e)rep(NA_real_,length(te)))
 success<-!inherits(fit,"error") && (model=="xgb_shallow" || startsWith(model,"ridge_") || isTRUE(fit$converged)) && all(is.finite(p))
 out[[length(out)+1L]]<-data.frame(region=r,group=label,species=s,fold=k,buffer_km=buffer,model=model,n_train=length(tr),n_test=length(te),train_presence=sum(y[tr]==1),test_presence=sum(y[te]==1),median_geo_km=median(nn2(coords[tr,,drop=FALSE],coords[te,,drop=FALSE],k=1)$nn.dists[,1])/1000,auc=if(success)auc(y[te],p)else NA_real_,ap=if(success)ap(y[te],p)else NA_real_,brier=if(success)mean((y[te]-p)^2)else NA_real_,success=success,warnings=paste(unique(warnings),collapse=";"))
 }
 }
 }
}
}
write.csv(do.call(rbind,coverage),"results/buffer_models/category_coverage.csv",row.names=FALSE)
print(aggregate(cbind(test_before,test_after,unseen_category)~region,do.call(rbind,coverage),sum))
a<-do.call(rbind,out);write.csv(a,"results/buffer_models/glm_fold_metrics.csv",row.names=FALSE)
models<-c("linear","quadratic","xgb_shallow","ridge_linear","ridge_quadratic")
b<-reshape(a[,c("region","group","species","fold","buffer_km","model","auc")],idvar=c("region","group","species","fold","buffer_km"),timevar="model",direction="wide")
cols<-paste0("auc.",models)
c<-b[complete.cases(b[,cols]),]
key<-interaction(c$region,c$group,c$species,c$fold,drop=TRUE)
c<-c[ave(c$buffer_km,key,FUN=function(x)length(unique(x)))==4,]
stopifnot(nrow(c)>0)
sp<-aggregate(c[,cols],c[c("region","group","species","buffer_km")],mean)
print(aggregate(sp[,cols],sp[c("region","buffer_km")],mean))
print(with(sp,table(region,buffer_km)))
print(with(a,table(region,model,success)))
write.csv(sp,"results/buffer_models/species_macro.csv",row.names=FALSE)
write.csv(c,"results/buffer_models/strict_matched_auc.csv",row.names=FALSE)
write.csv(b,"results/buffer_models/paired_auc.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"results/buffer_models/session_info.txt")
