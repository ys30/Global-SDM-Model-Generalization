library(disdat)
library(maxnet)
library(glmnet)
source('scripts/ridge_path.R')
library(xgboost)
dir.create("results/maxnet_sensitivity",recursive=TRUE,showWarnings=FALSE)
set.seed(20261125)
regions<-list(AWT=c("bird","plant"),CAN=NA_character_,NSW=c("ba","db","nb","ot","ou","rt","ru","sr"),NZ=NA_character_,SA=NA_character_,SWI=NA_character_)
cats<-list(AWT=character(),CAN="ontveg",NSW="vegsys",NZ=c("age","toxicats"),SA=character(),SWI="calc")
auc<-function(y,p){n<-as.double(sum(y==1));m<-as.double(sum(y==0));if(!n||!m)return(NA_real_);(sum(rank(p)[y==1])-n*(n+1)/2)/(n*m)}
ap<-function(y,p){i<-order(p,decreasing=TRUE);y<-y[i];p<-p[i];end<-c(which(diff(p)!=0),length(p));tp<-cumsum(y==1)[end];sum(diff(c(0,tp))/sum(y==1)*tp/end)}
grid<-expand.grid(feature=c("l","lq","lqh","lqhpt"),regmult=c(.5,1,2,4),stringsAsFactors=FALSE)
grid$model<-paste0("maxnet_",grid$feature,"_r",grid$regmult)
models<-grid$model
out<-list();coverage<-list()
target<-Sys.getenv("TARGET_REGION");stopifnot(target %in% names(regions))
for(r in target){
 po<-disPo(r);bg<-disBg(r);pred<-disPredictors(r)
 sid<-intersect(c("species","spid","species_id","sp","taxon"),names(po));stopifnot(length(sid)==1)
 catnames<-intersect(cats[[r]],pred);cont<-setdiff(pred,catnames)
 xykey<-function(d)paste(d$x,d$y,sep="|")
 for(g in regions[[r]]){
 group<-if(is.na(g))NULL else g;label<-if(is.na(g))"all"else g
 pa<-disPa(r,group);env<-disEnv(r,group);stopifnot(identical(pa[,1:4],env[,1:4]))
 # Remove background locations coinciding with any test location, and with focal presences.
 bg0<-bg[!xykey(bg)%in%xykey(pa),,drop=FALSE]
 for(sp in names(pa)[-(1:4)]){
 pres<-po[!is.na(po[[sid]])&as.character(po[[sid]])==sp,,drop=FALSE]
 pres<-pres[!xykey(pres)%in%xykey(pa),,drop=FALSE]
 back<-bg0[!xykey(bg0)%in%xykey(pres),,drop=FALSE]
 train<-rbind(pres[,pred,drop=FALSE],back[,pred,drop=FALSE]);y<-c(rep(1,nrow(pres)),rep(0,nrow(back)))
 test<-env[,pred,drop=FALSE];yt<-pa[[sp]]
 ok<-complete.cases(train);te<-complete.cases(test)&!is.na(yt)&yt%in%0:1
 for(p in cont){ok<-ok&is.finite(train[[p]]);te<-te&is.finite(test[[p]])}
 train<-train[ok,,drop=FALSE];y<-y[ok]
 if(sum(y==1)<10||sum(y==0)<10)next
 before<-sum(te)
 for(p in catnames){
 lev<-sort(unique(as.character(train[[p]])));te<-te&as.character(test[[p]])%in%lev
 train[[p]]<-factor(as.character(train[[p]]),levels=lev);test[[p]]<-factor(as.character(test[[p]]),levels=lev)
 }
 coverage[[length(coverage)+1]]<-data.frame(region=r,group=label,species=sp,po_n=sum(y==1),bg_n=sum(y==0),test_before=before,test_after=sum(te),test_presence=sum(yt[te]==1),test_absence=sum(yt[te]==0))
 if(length(unique(yt[te]))<2)next
 test<-test[te,,drop=FALSE];yt<-yt[te]
 keep<-vapply(train,function(z)length(unique(z))>1,logical(1));train<-train[,keep,drop=FALSE];test<-test[,keep,drop=FALSE]
 numeric_names<-intersect(cont,names(train))
 for(p in numeric_names){mu<-mean(train[[p]]);sdv<-sd(train[[p]]);train[[p]]<-(train[[p]]-mu)/sdv;test[[p]]<-(test[[p]]-mu)/sdv}
 # Full-rank treatment coding for GLM and ridge; MaxEnt receives actual factors.
 xx<-model.matrix(~.,train)[,-1,drop=FALSE];zz<-model.matrix(~.,test,contrasts.arg=attr(model.matrix(~.,train),"contrasts"))[,-1,drop=FALSE]
 stopifnot(identical(colnames(xx),colnames(zz)))
 for(model in models){
 warnings<-character()
 p<-tryCatch(withCallingHandlers({
 if(startsWith(model,"maxnet_")){
 i<-match(model,grid$model);cls<-grid$feature[i]
 f<-maxnet.formula(y,train,classes=cls)
 fit<-maxnet(y,train,f=f,regmult=grid$regmult[i],addsamplestobackground=FALSE)
 predict(fit,test,type="link",clamp=TRUE)
 }else{
 x<-xx;z<-zz
 if(model%in%c("glm_q","ridge_q")){x<-cbind(x,as.matrix(train[numeric_names])^2);z<-cbind(z,as.matrix(test[numeric_names])^2)}
 if(startsWith(model,"glm_")){
 fit<-glm.fit(cbind(1,x),y,family=binomial(),control=glm.control(maxit=50));if(!fit$converged)stop("nonconvergence")
 cf<-fit$coefficients;cf[is.na(cf)]<-0;as.numeric(cbind(1,z)%*%cf)
 }else if(startsWith(model,"ridge_")){
 fit<-ridge_path(x,y,.01);as.numeric(predict(fit,z,s=.01,type="link"))
 }else{
 fit<-xgb.train(params=list(objective="binary:logistic",max_depth=2,eta=.05,nthread=1),data=xgb.DMatrix(x,label=y),nrounds=100,verbose=0)
 predict(fit,xgb.DMatrix(z))
 }
 }
 },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e){warnings<<-c(warnings,conditionMessage(e));rep(NA_real_,length(yt))})
 good<-length(p)==length(yt)&&all(is.finite(p))
 out[[length(out)+1]]<-data.frame(region=r,group=label,species=sp,model=model,n_presence=sum(y==1),n_background=sum(y==0),n_test=length(yt),feature=grid$feature[match(model,grid$model)],regmult=grid$regmult[match(model,grid$model)],score_sd=if(good)sd(as.numeric(p))else NA_real_,auc=if(good)auc(yt,p)else NA_real_,ap=if(good)ap(yt,p)else NA_real_,success=good,warnings=paste(unique(warnings),collapse=";"))
 }
 }
 }
 cat("Completed",r,"\n")
}
a<-do.call(rbind,out);write.csv(a,"results/maxnet_sensitivity/metrics.csv",row.names=FALSE)
write.csv(do.call(rbind,coverage),"results/maxnet_sensitivity/coverage.csv",row.names=FALSE)
b<-reshape(a[c("region","group","species","model","auc")],idvar=c("region","group","species"),timevar="model",direction="wide")
cols<-paste0("auc.",models);stopifnot(all(cols%in%names(b)))
c<-b[complete.cases(b[,cols]),];cat("All-16 matched taxa:",nrow(c),"\n")
print(table(c$region));if(nrow(c))print(aggregate(c[,cols],c["region"],mean));print(with(a,table(region,model,success)))
write.csv(c,"results/maxnet_sensitivity/matched_species.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"results/maxnet_sensitivity/session_info.txt")
