trend_measure<-function(y_t, alpha=0.05){
  require(dlm)
  buildFun <- function(x) {
    m<- dlmModPoly(2, dV=exp(x[1]), dW= exp(c(x[2], x[3]))); 
    return(m)
  }
  fit<- dlmMLE(y_t, parm=c(0,0,0), build= buildFun)
  mod<- buildFun(fit$par)
  mod$m0[1] <- y_t[1]
  y_smooth<-dlmSmooth(y_t, mod)
  sde<-sqrt(sapply(dlmSvd2var(y_smooth$U.S, y_smooth$D.S), FUN=function(x){x[2,2]}))
  t<- qnorm(1-alpha/2)
  trend_upper <- y_smooth$s[,2] + t* sde
  trend_lower <- y_smooth$s[,2] - t* sde
  sig_diff_zero<-!((0> trend_lower) & (0< trend_upper))
  positive_size<- length(which(sig_diff_zero & (y_smooth$s[,2]>0)))
  negative_size <- length(which(sig_diff_zero & (y_smooth$s[,2]<0)))
  return(list("trend_or_not" = sig_diff_zero, "positive score" = min(1,positive_size/length(y_t)), 
              "neg score" = min(1,negative_size/length(y_t))))
}

# MCMC
trend_measure_baysian<-function(y_t, alpha=0.05, nsample = 1000, discountFactor = 0.05){
  require(dlm)
  buildFun <- function(x) {
    m<- dlmModPoly(2, dV=exp(x[1]), dW= exp(c(x[2], x[3]))); 
    return(m)
  }
  fit<- dlmMLE(y_t, parm=c(0,0,0), build= buildFun)
  mod<- buildFun(fit$par)
  nobs <- NROW(y_t)
  p <- nrow(mod$W)
  mod$m0[1] <- y_t[1]
  y_filter <- dlmFilter(y_t, mod)
  y_smooth<-dlmSmooth(y_t, mod)
  y.center<-y_t- tcrossprod(y_smooth$s[-1, , drop=FALSE], mod$FF)
  shape_par_y <-2
  shape_par_theta <- 2
  SSy <- crossprod(y.center)
  discountFactor<- 0.01
  rate.y <- 0.5*SSy *discountFactor
  shape.y <-  0.5* nobs * discountFactor
  # shape.y <- shape_par_y
  # df_y <- shape_par_y/(0.5*nobs)
  # rate.y <- 0.5 * SSy* df_y
  
  theta.center <- y_smooth$s[-1, , drop = FALSE] - tcrossprod(y_smooth$s[-(nobs + 1), , drop = FALSE], mod$GG)
  SStheta <- drop(sapply(1:p, function(i) crossprod(theta.center[, i])))
  shape.theta <- 0.5* rep(nobs, p) * discountFactor
  rate.theta <- 0.5* SStheta * discountFactor
  # shape.theta <- rep(shape_par_theta, p)
  # df_theta <- shape_par_theta/(0.5*nobs)
  # rate.theta <- 0.5 * SStheta * df_theta
  BSamples<-dlmGibbsDIG(y_t, mod = mod, shape.y = shape.y, rate.y= rate.y, shape.theta = shape.theta, rate.theta = rate.theta, thin = 2, n.sample = nsample, save.states = T, progressBar = F)
  betas<-sapply((dim(BSamples$theta)[3]/2):dim(BSamples$theta)[3], function(x){BSamples$theta[,2,x]})
  l_p <- alpha/2
  u_p <- 1-l_p
  quants<-sapply(1:nrow(betas), function(x){quantile(betas[x,], c(l_p, u_p))})
  trend_upper <- quants[2,]
  trend_lower <- quants[1,]
  sig_diff_zero<-!((0> trend_lower) & (0< trend_upper))
  positive_size<- length(which(sig_diff_zero & (trend_lower>0)))
  negative_size <- length(which(sig_diff_zero & (trend_upper<0)))
  return(list("trend_or_not" = sig_diff_zero, "positive score" = min(1,positive_size/length(y_t)), 
              "neg score" = min(1,negative_size/length(y_t))))
}

trend_measure_lm <- function(y_t){
  mdl<-lm(yt~t, data = data.frame(yt=y_t, t=1:length(y_t)))
  return(score<-summary(mdl)$r.squared)
}



# r2ww <- function(x){
#   SSe <- sum(x$w*(x$resid)^2); #the residual sum of squares is weighted
#   observed <- x$resid+x$fitted;   
#   SSt <- sum(x$w*(observed-weighted.mean(observed,x$w))^2)#the total sum of squares is weighted  
#   value <- 1-SSe/SSt;
#   return(value);
# }
# 
# 
# r2w <- function(x){
#    SSe <- sum((x$w*x$resid)^2); #the residual sum of squares is weighted
#    observed <- x$resid+x$fitted;
#    SSt <- sum((x$w*observed-mean(x$w*observed))^2); #the total sum of squares is weighted      
#    value <- 1-SSe/SSt;
#    return(value);
# }
