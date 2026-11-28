library(disdat)
library(sf)
library(RANN)
set.seed(20261125)
dir.create("results/novelty",recursive=TRUE,showWarnings=FALSE)
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
 coords<-matrix(NA_real_,nrow(pa),2);coords[ok,]<-xy
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


 # Distances use unlabeled test covariates and training sites only.
 geo<-nn2(coords[tr,,drop=FALSE],coords[te,,drop=FALSE],k=1)$nn.dists[,1]/1000
 if(!length(numeric_names))next
 z<-as.matrix(x[,numeric_names,drop=FALSE])
 envd<-nn2(z[tr,,drop=FALSE],z[te,,drop=FALSE],k=1)$nn.dists[,1]/sqrt(ncol(z))
 lo<-apply(z[tr,,drop=FALSE],2,min);hi<-apply(z[tr,,drop=FALSE],2,max)
 novel<-rowSums(sweep(z[te,,drop=FALSE],2,lo,"<")|sweep(z[te,,drop=FALSE],2,hi,">"))>0
 out[[length(out)+1]]<-data.frame(region=r,group=label,species=s,fold=k,n_train=length(tr),n_test=length(te),train_presence=sum(y[tr]==1),geo_median=median(geo),geo_p90=unname(quantile(geo,.9)),env_median=median(envd),env_p90=unname(quantile(envd,.9)),outside_fraction=mean(novel))
 }
 }
}
a<-do.call(rbind,out);stopifnot(nrow(a)>0,!anyDuplicated(a[c("region","group","species","fold")]))
write.csv(a,"results/novelty/features.csv",row.names=FALSE)
print(aggregate(a[c("geo_median","env_median","outside_fraction")],a["region"],median));print(table(a$region))
