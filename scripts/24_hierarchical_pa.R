library(disdat)
library(sf)
library(lme4)
library(glmnet)
set.seed(20261125)
dir.create("results/hierarchical",recursive=TRUE,showWarnings=FALSE)
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

 taxa<-names(pa)[-(1:4)]
 for(k in 1:5){
 tr<-which(ok&fold!=k);te<-which(ok&fold==k)
 if(!length(tr)||!length(te))next
 supported<-rep(TRUE,length(te))
 for(p in categorical)supported<-supported&env[[p]][te]%in%unique(env[[p]][tr])
 te<-te[supported];if(!length(te))next
 x<-as.data.frame(env[continuous]);names(x)<-paste0("v",seq_along(continuous))
 mu<-vapply(x[tr,,drop=FALSE],mean,0.0);sdv<-vapply(x[tr,,drop=FALSE],sd,0.0);keep<-is.finite(sdv)&sdv>0
 x<-as.data.frame(scale(x[,keep,drop=FALSE],center=mu[keep],scale=sdv[keep]));num<-names(x)
 for(j in seq_along(categorical)){
 lev<-sort(unique(env[[categorical[j]]][tr]))
 if(length(lev)>1)for(h in seq_along(lev)[-1])x[[paste0("cat",j,"_",h)]]<-as.numeric(env[[categorical[j]]]==lev[h])
 }
 # Include training-eligible taxa only; test outcomes never determine training inclusion.
 eligible<-taxa[vapply(taxa,function(s){y<-pa[[s]][tr];sum(y==1,na.rm=TRUE)>=5&&sum(y==0,na.rm=TRUE)>=5},logical(1))]
 if(length(eligible)<2||!ncol(x))next
 long<-function(ids){
 do.call(rbind,lapply(eligible,function(s){
 y<-pa[[s]][ids];v<-!is.na(y)&y%in%0:1
 data.frame(x[ids[v],,drop=FALSE],y=y[v],taxon=factor(s,levels=eligible),row_id=ids[v])
 }))
 }
 train<-long(tr);test<-long(te)
 for(spec in c("linear","quadratic")){
 a<-train;b<-test
 if(spec=="quadratic")for(p in num){a[[paste0(p,"_sq")]]<-a[[p]]^2;b[[paste0(p,"_sq")]]<-b[[p]]^2}
 vars<-setdiff(names(a),c("y","taxon","row_id"))
 # Common fixed slopes, partially pooled taxon intercepts; not random slopes.
 f<-as.formula(paste("y~",paste(vars,collapse="+"),"+(1|taxon)"))
 warn<-character()
 fit<-tryCatch(withCallingHandlers(glmer(f,data=a,family=binomial(),nAGQ=1,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
 good<-!inherits(fit,"error")
 singular<-if(good)isSingular(fit)else NA
 messages<-if(good)paste(fit@optinfo$conv$lme4$messages,collapse=";")else conditionMessage(fit)
 opt<-if(good)fit@optinfo$conv$opt else NA
 pred<-if(good)tryCatch(predict(fit,newdata=b,type="response"),error=function(e)rep(NA_real_,nrow(b)))else rep(NA_real_,nrow(b))
 coverage[[length(coverage)+1]]<-data.frame(region=r,group=label,fold=k,spec=spec,n_taxa=length(eligible),n_train=nrow(a),n_test=nrow(b),singular=singular,opt_code=if(length(opt))opt else NA,convergence=messages,warnings=paste(unique(warn),collapse=";"))
 for(s in eligible){
 ia<-which(a$taxon==s);ib<-which(b$taxon==s)
 if(length(unique(b$y[ib]))<2)next
 xx<-as.matrix(a[ia,vars,drop=FALSE]);zz<-as.matrix(b[ib,vars,drop=FALSE])
 p0<-tryCatch(suppressWarnings({m<-glm.fit(cbind(1,xx),a$y[ia],family=binomial(),control=glm.control(maxit=50));if(!m$converged)stop("nonconvergence");cf<-m$coefficients;cf[is.na(cf)]<-0;plogis(cbind(1,zz)%*%cf)}),error=function(e)rep(NA_real_,length(ib)))
 p1<-tryCatch(as.numeric(predict(glmnet(xx,a$y[ia],family="binomial",alpha=0,lambda=.01),zz,s=.01,type="response")),error=function(e)rep(NA_real_,length(ib)))
 for(model in c("independent_glm","independent_ridge","hierarchical_intercept")){
 p<-switch(model,independent_glm=p0,independent_ridge=p1,hierarchical_intercept=pred[ib])
 success<-all(is.finite(p))&& (model!="hierarchical_intercept"|| (good&&isTRUE(opt==0)&&!nzchar(messages)))
 out[[length(out)+1]]<-data.frame(region=r,group=label,species=s,fold=k,spec=spec,model=model,n_train=length(ia),n_test=length(ib),auc=if(success)auc(b$y[ib],p)else NA_real_,ap=if(success)ap(b$y[ib],p)else NA_real_,brier=if(success)mean((b$y[ib]-p)^2)else NA_real_,success=success)
 }
 }
 }
 cat("Completed",r,label,k,"\n")
 }
 if(length(out))write.csv(do.call(rbind,out),"results/hierarchical/metrics.csv",row.names=FALSE)
 if(length(coverage))write.csv(do.call(rbind,coverage),"results/hierarchical/fit_audit.csv",row.names=FALSE)
}
a<-do.call(rbind,out);stopifnot(nrow(a)>0)
b<-reshape(a[c("region","group","species","fold","spec","model","auc")],idvar=c("region","group","species","fold","spec"),timevar="model",direction="wide")
cols<-grep("^auc\\.",names(b),value=TRUE);b<-b[complete.cases(b[cols]),]
stopifnot(nrow(b)>0)
sp<-aggregate(b[cols],b[c("region","group","species","spec")],mean)
print(aggregate(sp[cols],sp[c("region","spec")],mean));print(with(sp,table(region,spec)))
print(with(a,table(spec,model,success)))
write.csv(sp,"results/hierarchical/species_macro.csv",row.names=FALSE)
write.csv(b,"results/hierarchical/matched.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"results/hierarchical/session_info.txt")
