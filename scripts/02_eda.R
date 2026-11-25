library(disdat)
dir.create("results/eda",recursive=TRUE,showWarnings=FALSE)
regions <- list(AWT=c("bird","plant"),CAN=NA_character_,NSW=c("ba","db","nb","ot","ou","rt","ru","sr"),NZ=NA_character_,SA=NA_character_,SWI=NA_character_)
groups <- list(); predictors <- list(); folds <- list()
key <- function(d) do.call(paste,c(lapply(d[,1:4,drop=FALSE],as.character),sep="\034"))
for(r in names(regions)) for(g in regions[[r]]) {
 group <- if(is.na(g)) NULL else g
 label <- if(is.na(g)) "all" else g
 pa <- disPa(r,group); env <- disEnv(r,group)
 stopifnot(identical(names(pa)[1:4],names(env)[1:4]))
 kp <- key(pa); ke <- key(env)
 unique_keys <- !anyDuplicated(kp) && !anyDuplicated(ke)
 matched <- unique_keys && length(kp)==length(ke) && setequal(kp,ke)
 xy <- intersect(c("x","y"),names(pa))
 pred <- names(env)[-(1:4)]
 groups[[length(groups)+1]] <- data.frame(region=r,group=label,pa_rows=nrow(pa),env_rows=nrow(env),species=ncol(pa)-4,predictors=length(pred),unique_keys=unique_keys,matched=matched,same_order=identical(kp,ke),crs=paste(disCRS(r),collapse=" "),duplicate_coordinates=if(length(xy)==2)sum(duplicated(pa[,xy])) else NA_integer_)
 for(p in pred){
 z <- env[[p]]
 predictors[[length(predictors)+1]] <- data.frame(region=r,group=label,predictor=p,type=paste(class(z),collapse="/"),missing=sum(is.na(z)),unique_values=length(unique(z[!is.na(z)])),nonfinite=if(is.numeric(z))sum(!is.finite(z)&!is.na(z)) else NA_integer_)
 }
 if(!matched || length(xy)!=2) next
 env <- env[match(kp,ke),,drop=FALSE]
 ok <- complete.cases(env[,pred,drop=FALSE]) & complete.cases(pa[,xy,drop=FALSE])
 for(p in pred) if(is.numeric(env[[p]])) ok <- ok & is.finite(env[[p]])
 # Diagnostic partition only: five contiguous x-coordinate bands, ties never split.
 xx <- pa$x
 cuts <- unique(as.numeric(quantile(xx[ok],probs=seq(0,1,length.out=6),na.rm=TRUE)))
 if(length(cuts)<2) next
 f <- cut(xx,breaks=cuts,include.lowest=TRUE,labels=FALSE)
 for(s in names(pa)[-(1:4)]) for(k in seq_len(length(cuts)-1)){
 y <- pa[[s]]; eligible <- ok & !is.na(y) & y %in% c(0,1) & !is.na(f)
 tr <- eligible & f!=k; te <- eligible & f==k
 folds[[length(folds)+1]] <- data.frame(region=r,group=label,species=s,fold=k,train_presence=sum(y[tr]==1),train_absence=sum(y[tr]==0),test_presence=sum(y[te]==1),test_absence=sum(y[te]==0),complete_rows=sum(eligible),eligible_auc=sum(y[tr]==1)>=5 & sum(y[tr]==0)>=5 & sum(y[te]==1)>0 & sum(y[te]==0)>0)
 }
}
a <- do.call(rbind,groups); b <- do.call(rbind,predictors); c <- do.call(rbind,folds)
write.csv(a,"results/eda/group_audit.csv",row.names=FALSE)
write.csv(b,"results/eda/predictor_audit.csv",row.names=FALSE)
write.csv(c,"results/eda/fold_feasibility.csv",row.names=FALSE)
print(a,row.names=FALSE)
print(aggregate(eligible_auc~region,c,mean))
stopifnot(all(a$matched))
