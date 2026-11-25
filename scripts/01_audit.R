dir.create("results/audit", recursive=TRUE, showWarnings=FALSE)
library(disdat)
regions <- list(AWT=c("bird","plant"), CAN=NA_character_, NSW=c("ba","db","nb","ot","ou","rt","ru","sr"), NZ=NA_character_, SA=NA_character_, SWI=NA_character_)
inventory <- list(); schemas <- list()
schema <- function(x,r,g,t) data.frame(region=r,group=g,table=t,rows=nrow(x),column=names(x),class=vapply(x,function(z) paste(class(z),collapse="/"),""),missing=vapply(x,function(z) sum(is.na(z)),integer(1)))
for(r in names(regions)){
 for(t in c("po","bg")){
 x <- if(t=="po") disPo(r) else disBg(r)
 schemas[[length(schemas)+1L]] <- schema(x,r,"all",t)
 }
 for(g in regions[[r]]){
  group <- if(is.na(g)) NULL else g
  label <- if(is.na(g)) "all" else g
  pa <- disPa(r,group); env <- disEnv(r,group)
  schemas[[length(schemas)+1L]] <- schema(pa,r,label,"pa")
  schemas[[length(schemas)+1L]] <- schema(env,r,label,"env")
  stopifnot(ncol(pa)>4)
  for(s in names(pa)[-(1:4)]){
   y <- pa[[s]]
   valid <- !is.na(y) & y %in% c(0,1)
   n <- sum(valid); p <- sum(y[valid]==1); a <- sum(y[valid]==0)
   inventory[[length(inventory)+1L]] <- data.frame(region=r,group=label,species=s,rows=length(y),presence=p,absence=a,missing=sum(is.na(y)),invalid=sum(!is.na(y)&!valid),prevalence=if(n>0)p/n else NA,eligible_initial=p>=10 & a>=10)
  }
 }
}
i <- do.call(rbind,inventory); s <- do.call(rbind,schemas)
write.csv(i,"results/audit/species_inventory.csv",row.names=FALSE)
write.csv(s,"results/audit/table_schema.csv",row.names=FALSE)
write.csv(aggregate(cbind(presence,absence,missing,invalid,eligible_initial)~region,i,sum),"results/audit/region_counts.csv",row.names=FALSE)
write.csv(as.data.frame(table(i$region)),"results/audit/species_by_region.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"results/audit/session_info.txt")
print(table(i$region)); print(nrow(i))
if(any(i$invalid>0)) warning("Nonbinary labels detected; review inventory before modelling")
if(nrow(i)!=226) warning("Candidate species count differs from published 226; reconcile before modelling")
