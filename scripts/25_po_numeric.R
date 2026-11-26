library(disdat)
library(glmnet)
dir.create("results/po_numeric",recursive=TRUE,showWarnings=FALSE)
set.seed(20261125)
regions<-list(AWT=c("bird","plant"),CAN=NA_character_,NSW=c("ba","db","nb","ot","ou","rt","ru","sr"),NZ=NA_character_,SA=NA_character_,SWI=NA_character_)
cats<-list(AWT=character(),CAN="ontveg",NSW="vegsys",NZ=c("age","toxicats"),SA=character(),SWI="calc")
auc<-function(y,p){n<-as.double(sum(y==1));m<-as.double(sum(y==0));if(!n||!m)return(NA_real_);(sum(rank(p)[y==1])-n*(n+1)/2)/(n*m)}
ap<-function(y,p){i<-order(p,decreasing=TRUE);y<-y[i];p<-p[i];end<-c(which(diff(p)!=0),length(p));tp<-cumsum(y==1)[end];sum(diff(c(0,tp))/sum(y==1)*tp/end)}
out<-list();coverage<-list();features<-list()
for(r in c("CAN","AWT")){
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

 for(spec in c("linear","quadratic")){
 x<-xx;z<-zz
 if(spec=="quadratic"){q<-as.matrix(train[numeric_names])^2;qz<-as.matrix(test[numeric_names])^2;colnames(q)<-colnames(qz)<-paste0(numeric_names,"_sq");x<-cbind(x,q);z<-cbind(z,qz)}
 stopifnot(identical(colnames(x),colnames(z)))
 mu<-colMeans(x);sdv<-apply(x,2,sd)
 for(j in seq_len(ncol(x)))features[[length(features)+1]]<-data.frame(region=r,group=label,species=sp,spec=spec,feature=colnames(x)[j],train_sd=sdv[j],train_min=min(x[,j]),train_max=max(x[,j]),test_min=min(z[,j]),test_max=max(z[,j]),outside=mean(z[,j]<min(x[,j])|z[,j]>max(x[,j])),max_abs_z=if(sdv[j]>0)max(abs((z[,j]-mu[j])/sdv[j]))else NA_real_)
 keep<-is.finite(sdv)&sdv>1e-10
 xs<-scale(x[,keep,drop=FALSE],center=mu[keep],scale=sdv[keep]);zs<-scale(z[,keep,drop=FALSE],center=mu[keep],scale=sdv[keep])
 # Hold lambda and objective fixed; independent numerical solver in standardized space.
 lam<-.01
 loss<-function(b){eta<-drop(b[1]+xs%*%b[-1]);mean(pmax(eta,0)+log1p(exp(-abs(eta)))-y*eta)+lam*sum(b[-1]^2)/2}
 grad<-function(b){eta<-drop(b[1]+xs%*%b[-1]);e<-plogis(eta)-y;c(mean(e),drop(crossprod(xs,e))/length(y)+lam*b[-1])}
 opt<-optim(c(qlogis(mean(y)),rep(0,ncol(xs))),loss,grad,method="BFGS",control=list(maxit=2000,reltol=1e-12))
 path<-sort(unique(c(exp(seq(log(1),log(1e-5),length.out=100)),lam)),decreasing=TRUE)
 for(method in c("single_default","path_default","path_explicit","optim_explicit")){
 f<-NULL;warn<-character()
 p<-tryCatch(withCallingHandlers({
 if(method=="optim_explicit")drop(opt$par[1]+zs%*%opt$par[-1])else{
 explicit<-method=="path_explicit";a<-if(explicit)xs else x;b<-if(explicit)zs else z
 f<-glmnet(a,y,family="binomial",alpha=0,lambda=if(method=="single_default")lam else path,standardize=!explicit,thresh=1e-12,maxit=1000000)
 as.numeric(predict(f,b,s=lam,type="link"))
 }
 },warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e){warn<<-c(warn,conditionMessage(e));rep(NA_real_,length(yt))})
 good<-length(p)==length(yt)&&all(is.finite(p))
 out[[length(out)+1]]<-data.frame(region=r,group=label,species=sp,spec=spec,method=method,auc=if(good)auc(yt,p)else NA_real_,score_sd=if(good)sd(p)else NA_real_,score_unique=if(good)length(unique(p))else NA_integer_,optim_code=opt$convergence,optim_max_grad=max(abs(grad(opt$par))),glmnet_jerr=if(is.null(f))NA_integer_ else f$jerr,warnings=paste(unique(warn),collapse=";"))
 }
 }
 }
 }
}
a<-do.call(rbind,out);f<-do.call(rbind,features)
write.csv(a,"results/po_numeric/metrics.csv",row.names=FALSE)
write.csv(f,"results/po_numeric/features.csv",row.names=FALSE)
print(aggregate(cbind(auc,score_sd)~region+spec+method,a,mean))
print(aggregate(max_abs_z~region+feature,f,max))
print(head(f[order(f$max_abs_z,decreasing=TRUE),],20))
print(aggregate(optim_max_grad~region+spec,a,max))
writeLines(capture.output(sessionInfo()),"results/po_numeric/session_info.txt")
