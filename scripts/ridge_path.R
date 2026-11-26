# Stable ridge path with explicit target and convergence checks.
ridge_path<-function(x,y,lambda){
 grid<-sort(unique(c(exp(seq(log(max(1,lambda)),log(min(1e-5,lambda)),length.out=100)),lambda)),decreasing=TRUE)
 f<-glmnet(x,y,family="binomial",alpha=0,lambda=grid,standardize=TRUE,thresh=1e-12,maxit=1000000)
 if(f$jerr!=0||!any(abs(f$lambda-lambda)<1e-12))stop("ridge path incomplete")
 f
}
