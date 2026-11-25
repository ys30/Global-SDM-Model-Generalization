library(disdat)
dir.create("results/po_audit",recursive=TRUE,showWarnings=FALSE)
regions<-list(AWT=c("bird","plant"),CAN=NA_character_,NSW=c("ba","db","nb","ot","ou","rt","ru","sr"),NZ=NA_character_,SA=NA_character_,SWI=NA_character_)
cats<-list(AWT=character(),CAN="ontveg",NSW="vegsys",NZ=c("age","toxicats"),SA=character(),SWI="calc")
inventory<-list(); metadata<-list(); schemas<-list()
for(r in names(regions)){
 po<-disPo(r);bg<-disBg(r);pred<-disPredictors(r)
 print(list(region=r,po_names=names(po),bg_names=names(bg),predictors=pred))
 stopifnot(all(pred %in% names(po)),all(pred %in% names(bg)))
 # Resolve species identifier explicitly; fail on ambiguous schema.
 candidates<-intersect(c("species","spid","species_id","sp","taxon"),names(po))
 if(length(candidates)!=1)stop("Unresolved species identifier: ",r)
 sid<-candidates[1]
 xykey<-function(d)paste(d$x,d$y,sep="|")
 for(g in regions[[r]]){
 group<-if(is.na(g))NULL else g; label<-if(is.na(g))"all"else g
 pa<-disPa(r,group);env<-disEnv(r,group)
 stopifnot(identical(pa[,1:4],env[,1:4]),all(pred %in% names(env)))
 for(p in pred){
 vals<-function(d)paste(head(sort(unique(as.character(d[[p]][!is.na(d[[p]])]))),20),collapse="|")
 metadata[[length(metadata)+1L]]<-data.frame(region=r,group=label,predictor=p,documented_categorical=p %in% cats[[r]],po_class=paste(class(po[[p]]),collapse="/"),pa_class=paste(class(env[[p]]),collapse="/"),po_unique=length(unique(na.omit(po[[p]]))),pa_unique=length(unique(na.omit(env[[p]]))),po_missing=sum(is.na(po[[p]])),bg_missing=sum(is.na(bg[[p]])),pa_missing=sum(is.na(env[[p]])),po_values=if(p %in% cats[[r]])vals(po)else "",pa_values=if(p %in% cats[[r]])vals(env)else "")
 }
 for(sp in names(pa)[-(1:4)]){
 q<-po[!is.na(po[[sid]]) & as.character(po[[sid]])==sp,,drop=FALSE]
 y<-pa[[sp]]
 inventory[[length(inventory)+1L]]<-data.frame(region=r,group=label,species=sp,po_records=nrow(q),po_unique_coords=length(unique(xykey(q))),bg_records=nrow(bg),pa_presence=sum(y==1,na.rm=TRUE),pa_absence=sum(y==0,na.rm=TRUE),po_pa_overlap=length(intersect(xykey(q),xykey(pa))),bg_pa_overlap=length(intersect(xykey(bg),xykey(pa))),eligible=nrow(q)>=10 & sum(y==1,na.rm=TRUE)>0 & sum(y==0,na.rm=TRUE)>0)
 }
 }
}
a<-do.call(rbind,inventory);b<-do.call(rbind,metadata)
write.csv(a,"results/po_audit/species_inventory.csv",row.names=FALSE)
write.csv(b,"results/po_audit/predictor_metadata.csv",row.names=FALSE)
print(aggregate(cbind(po_records,po_unique_coords,po_pa_overlap,eligible)~region,a,sum))
print(b[b$documented_categorical,c("region","group","predictor","po_class","pa_class","po_values","pa_values")],row.names=FALSE)
stopifnot(nrow(a)==226,all(a$po_records>0))
