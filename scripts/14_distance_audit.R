library(disdat)
library(sf)
library(RANN)
dir.create("results/distance",recursive=TRUE,showWarnings=FALSE)
set.seed(20261125)
regions<-list(AWT=c("bird","plant"),CAN=NA_character_,NSW=c("ba","db","nb","ot","ou","rt","ru","sr"),NZ=NA_character_,SA=NA_character_,SWI=NA_character_)
cats<-list(AWT=character(),CAN="ontveg",NSW="vegsys",NZ=c("age","toxicats"),SA=character(),SWI="calc")
out<-list(); buf<-list(); foldout<-list()
for(r in names(regions))for(g in regions[[r]]){
 group<-if(is.na(g))NULL else g;label<-if(is.na(g))"all"else g
 pa<-disPa(r,group);env<-disEnv(r,group);stopifnot(identical(pa[,1:4],env[,1:4]))
 pred<-disPredictors(r);cont<-setdiff(pred,cats[[r]])
 ok<-complete.cases(env[pred])&is.finite(pa$x)&is.finite(pa$y)
 for(p in cont)ok<-ok&is.finite(env[[p]])
 pts<-st_as_sf(pa[ok,],coords=c("x","y"),crs=disCRS(r))
 ll<-st_coordinates(st_transform(pts,4326))
 lon<-atan2(mean(sin(ll[,1]*pi/180)),mean(cos(ll[,1]*pi/180)))*180/pi
 lat<-mean(ll[,2])
 xy<-st_coordinates(st_transform(pts,sprintf("+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m",lat,lon)))
 u<-unique(xy);if(nrow(u)<5)next
 centers<-kmeans(u,centers=5,nstart=10)$centers
 f<-max.col(-sapply(1:5,function(j)rowSums(sweep(xy,2,centers[j,])^2)),ties.method="first")
 ids<-match(paste(xy[,1],xy[,2],sep="|"),unique(paste(xy[,1],xy[,2],sep="|")))
 rf<-sample(rep(1:5,length.out=max(ids)))[ids]
 foldout[[length(foldout)+1]]<-data.frame(region=r,group=label,row_id=which(ok),spatial=f,random=rf)
 for(scheme in c("spatial","random")){
 ff<-if(scheme=="spatial")f else rf
 for(k in 1:5){
 tr<-which(ff!=k);te<-which(ff==k);if(!length(tr)||!length(te))next
 dist<-nn2(xy[tr,,drop=FALSE],xy[te,,drop=FALSE],k=1)$nn.dists[,1]/1000
 quant<-quantile(dist,c(0,.1,.25,.5,.75,.9,1))
 # Environmental nearest-neighbour distance, scaled only on training records.
 x<-as.matrix(env[ok,cont,drop=FALSE]);mu<-colMeans(x[tr,,drop=FALSE]);sdv<-apply(x[tr,,drop=FALSE],2,sd)
 keep<-is.finite(sdv)&sdv>0
 x<-scale(x[,keep,drop=FALSE],center=mu[keep],scale=sdv[keep])
 ed<-nn2(x[tr,,drop=FALSE],x[te,,drop=FALSE],k=1)$nn.dists[,1]/sqrt(ncol(x))
 out[[length(out)+1]]<-data.frame(region=r,group=label,scheme=scheme,fold=k,n_train=length(tr),n_test=length(te),min_km=quant[1],p10_km=quant[2],p25_km=quant[3],median_km=quant[4],p75_km=quant[5],p90_km=quant[6],max_km=quant[7],median_env_distance=median(ed))
 if(scheme!="spatial")next
 td<-nn2(xy[te,,drop=FALSE],xy[tr,,drop=FALSE],k=1)$nn.dists[,1]/1000
 for(b in c(0,5,10,25,50,100)){
 tt<-tr[td>b]
 for(sp in names(pa)[-(1:4)]){
 y<-pa[[sp]][ok]
 if(sum(y==1,na.rm=TRUE)<10||sum(y==0,na.rm=TRUE)<10)next
 np<-sum(y[tt]==1,na.rm=TRUE);na<-sum(y[tt]==0,na.rm=TRUE)
 buf[[length(buf)+1]]<-data.frame(region=r,group=label,species=sp,fold=k,buffer_km=b,n_train=length(tt),train_presence=np,train_absence=na,test_presence=sum(y[te]==1,na.rm=TRUE),test_absence=sum(y[te]==0,na.rm=TRUE),eligible=np>=5&na>=5&sum(y[te]==1,na.rm=TRUE)>0&sum(y[te]==0,na.rm=TRUE)>0)
 }
 }
 }
 }
}
a<-do.call(rbind,out);b<-do.call(rbind,buf)
write.csv(a,"results/distance/fold_distances.csv",row.names=FALSE)
write.csv(b,"results/distance/buffer_feasibility.csv",row.names=FALSE)
write.csv(do.call(rbind,foldout),"results/distance/fold_assignments.csv",row.names=FALSE)
print(aggregate(cbind(median_km,median_env_distance)~region+scheme,a,median))
print(aggregate(eligible~region+buffer_km,b,mean))
