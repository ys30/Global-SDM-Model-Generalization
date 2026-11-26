library(disdat)
library(sf)
library(glmnet)
set.seed(20261125)
dir.create("results/pa_ridge_replay",recursive=TRUE,showWarnings=FALSE)
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
target<-Sys.getenv("TARGET_REGION");stopifnot(target %in% names(regions))
for(r in target) for(g in regions[[r]]){
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
 set.seed(20261125L+match(r,names(regions))*100L+match(g,regions[[r]]))
 centers<-kmeans(u,centers=5,nstart=10)$centers
 distmat<-sapply(1:5,function(j)rowSums(sweep(xy,2,centers[j,])^2))
 fold<-rep(NA_integer_,nrow(pa));fold[ok]<-max.col(-distmat,ties.method="first")
 for(s in names(pa)[-(1:4)]){
 y<-pa[[s]]; valid<-ok & !is.na(y) & y %in% c(0,1)
 if(sum(y[valid]==1)<10||sum(y[valid]==0)<10)next
 for(k in 1:5){
 tr<-which(valid & fold!=k); te<-which(valid & fold==k)
 if(sum(y[tr]==1)<5||sum(y[tr]==0)<5||length(unique(y[te]))<2)next
 # Fit all transformations on training records only. Unsupported test levels
 # are excluded identically for all models and explicitly counted.
 old_te<-length(te)
 supported<-rep(TRUE,length(te))
 for(p in categorical)supported<-supported & env[[p]][te] %in% unique(env[[p]][tr])
 te<-te[supported]
 coverage[[length(coverage)+1L]]<-data.frame(region=r,group=label,species=s,fold=k,test_before=old_te,test_after=length(te),unseen_category=old_te-length(te))
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

 for(spec in c("linear","quadratic")){
 d<-x
 if(spec=="quadratic")for(p in numeric_names)d[[paste0(p,"_sq")]]<-x[[p]]^2
 xx<-as.matrix(d[tr,,drop=FALSE]);zz<-as.matrix(d[te,,drop=FALSE])
 mu<-colMeans(xx);sdv<-apply(xx,2,sd);keep<-is.finite(sdv)&sdv>1e-10
 xs<-scale(xx[,keep,drop=FALSE],center=mu[keep],scale=sdv[keep]);zs<-scale(zz[,keep,drop=FALSE],center=mu[keep],scale=sdv[keep])
 yy<-y[tr];lam<-.01
 loss<-function(b){eta<-drop(b[1]+xs%*%b[-1]);mean(pmax(eta,0)+log1p(exp(-abs(eta)))-yy*eta)+lam*sum(b[-1]^2)/2}
 grad<-function(b){eta<-drop(b[1]+xs%*%b[-1]);e<-plogis(eta)-yy;c(mean(e),drop(crossprod(xs,e))/length(yy)+lam*b[-1])}
 opt<-optim(c(qlogis(mean(yy)),rep(0,ncol(xs))),loss,grad,method="BFGS",control=list(maxit=2000,reltol=1e-12))
 path<-sort(unique(c(exp(seq(log(1),log(1e-5),length.out=100)),lam)),decreasing=TRUE)
 for(method in c("single_original","path_default","path_explicit","optim_explicit")){
 warn<-character();f<-NULL
 p<-tryCatch(withCallingHandlers({
 if(method=="optim_explicit")drop(opt$par[1]+zs%*%opt$par[-1])else{
 explicit<-method=="path_explicit";a<-if(explicit)xs else xx;b<-if(explicit)zs else zz
 f<-if(method=="single_original")glmnet(a,yy,family="binomial",alpha=0,lambda=lam,standardize=TRUE)else glmnet(a,yy,family="binomial",alpha=0,lambda=path,standardize=!explicit,thresh=1e-12,maxit=1000000)
 as.numeric(predict(f,b,s=lam,type="link"))
 }
 },warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e){warn<<-c(warn,conditionMessage(e));rep(NA_real_,length(te))})
 good<-length(p)==length(te)&&all(is.finite(p))
 out[[length(out)+1]]<-data.frame(region=r,group=label,species=s,fold=k,spec=spec,method=method,n_train=length(tr),n_test=length(te),train_presence=sum(yy),auc=if(good)auc(y[te],plogis(p))else NA_real_,auc_link=if(good)auc(y[te],p)else NA_real_,score_sd=if(good)sd(p)else NA_real_,optim_code=opt$convergence,optim_max_grad=max(abs(grad(opt$par))),glmnet_jerr=if(is.null(f))NA_integer_ else f$jerr,warnings=paste(unique(warn),collapse=";"))
 }
 }
 }
 }
 if(length(out))write.csv(do.call(rbind,out),"results/pa_ridge_replay/metrics.csv",row.names=FALSE)
}
a<-do.call(rbind,out);stopifnot(nrow(a)>0)
b<-reshape(a[c("region","group","species","fold","spec","method","auc")],idvar=c("region","group","species","fold","spec"),timevar="method",direction="wide")
cols<-grep("^auc\\.",names(b),value=TRUE);b<-b[complete.cases(b[cols]),]
b$delta_path<-b$auc.path_default-b$auc.single_original
b$delta_solver<-b$auc.path_explicit-b$auc.optim_explicit
sp<-aggregate(b[c(cols,"delta_path","delta_solver")],b[c("region","group","species","spec")],mean)
print(aggregate(sp[c(cols,"delta_path","delta_solver")],sp[c("region","spec")],mean))
print(aggregate(abs(delta_path)~region+spec,b,max))
print(aggregate(abs(delta_solver)~region+spec,b,max))
print(aggregate(optim_max_grad~region+spec,a,max))
old<-read.csv("inputs/glm_fold_metrics.csv");old<-old[old$model%in%c("ridge_linear","ridge_quadratic"),];old$spec<-sub("ridge_","",old$model)
j<-merge(b,old[c("region","group","species","fold","spec","auc")],by=c("region","group","species","fold","spec"))
j$replay_error<-j$auc.single_original-j$auc
print(aggregate(abs(replay_error)~region+spec,j,max))
write.csv(j,"results/pa_ridge_replay/paired_replay.csv",row.names=FALSE)
write.csv(sp,"results/pa_ridge_replay/species.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"results/pa_ridge_replay/session_info.txt")
