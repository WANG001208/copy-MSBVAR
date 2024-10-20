### msvar.R -- Estimates the MLE for MSVAR and MS univariate
### models.

### hregime在得到fph后进行加权MLE
### 得到fph是通过blkopt，即EM algorithm算出来的。所谓的block optimization，就是先固定除mu以外的其他参数，优化mu，然后再优化别人，类似于EM algorithm，但不完全一样.

# 20120109 : Initial version by Ryan Davis


msvar <- function(Y, p, h, niterblkopt=10)
{

# Switching indicator: 'IAH' totally switching
# fixed for now
    indms <- 'IAH'

    n <- nrow(Y)
    m <- ncol(Y)

    # check value of h
    if(h<2) stop("h should be an integer >=2")

# Now do a baseline, non-regime model using szbvar() since this
# gives us all of the inputs we need for later.
# n.b. the below specification is equivalent to MLE
#      (but if mu5 or mu6 does not equal zero, then problems)
    init.model <- szbvar(ts(Y), p,
                         lambda0=1, lambda1=1, lambda3=1, lambda4=1,
                         lambda5=1, mu5=0, mu6=0, prior=2)

    # print(init.model$ar.coefs)
    # print(init.model$intercept)

# set initial parameters for blockwise optimization
# initial Q
    Qhat.start <- (1-(h*0.1/(h-1)))*diag(h) + matrix(0.1/(h-1), h, h)

# array for storage
    thetahat.start1 <- array(NA, c(m, 1+m*p+m, h))
    thetahat.start2 <- array(NA, c(m, 1+m*p+m, h))
    thetahat.start3 <- array(NA, c(m, 1+m*p+m, h))

# set intercept and AR coef initial values all to zero
# 这里和YaweiZhao的原代码不同，我在这里进行了修改
# /Here is different than the original code by Yawei Zhao, I made some adjustment here
# 用0作为初始值并不合适，而且在面对高维度的MSM问题时，问题是非凸的，需要多个initial value 以争取达到最好的值(先试一下用szbvar的结果如何)
# /It is inappropriate to use 0 as initial value, in addition, when we are facing high-dimensional MSM problems, the problem is non-convex, we need sample several different intial values to get the best we could get.(Firstly, let us try the results of using the value of szbvar as initial value.)
    for (j in 1:h){
        for (i in 1:m){
            start=1+(i-1)*m*p
            end=i*m*p
            thetahat.start1[i,1:(m*p),j] <- init.model$ar.coefs[start:end]
            # print(init.model$ar.coefs[start:end])
            thetahat.start1[i,1+(m*p),j] <- init.model$intercept[i]+(j-1)/100 # Adjust by a small amount over regimes for convinence in optimization.
            # print(i)
            # print(init.model$intercept[i])
            }
        }
    
    thetahat.start2[,1:(1+m*p),] <- 0
    thetahat.start2[1,1,] <- 0.4
    thetahat.start2[2,2,] <- 0.4

    print(thetahat.start2)

    thetahat.start3[,1:(1+m*p),] <- 0
    thetahat.start3[1,1,] <- -0.4
    thetahat.start3[2,2,] <- -0.4

    print(thetahat.start3)

# set sigma initial values
# first, get residuals from initial model
# dummy obs are appended, so adjust for those
    res.im <- init.model$residuals[(m+2):n,]
    sig2.start <- (1/n)*crossprod(res.im, res.im)

# the sig2 starting values need to be different though,
# so adjust these by a small amount over regimes
    for (i in 1:h) { thetahat.start1[,(1+m*p+1):(1+m*p+m),i] <-
                         sig2.start}
    for (i in 1:h) { thetahat.start2[,(1+m*p+1):(1+m*p+m),i] <-
                         sig2.start}
    for (i in 1:h) { thetahat.start3[,(1+m*p+1):(1+m*p+m),i] <-
                     sig2.start}

    blkopt_est1 <- blkopt(Y=Y, p=p, thetahat.start=thetahat.start1,
                         Qhat.start=Qhat.start, niter=niterblkopt,
                         indms)

    blkopt_est2 <- blkopt(Y=Y, p=p, thetahat.start=thetahat.start2,
                         Qhat.start=Qhat.start, niter=niterblkopt,
                         indms)

    blkopt_est3 <- blkopt(Y=Y, p=p, thetahat.start=thetahat.start3,
                     Qhat.start=Qhat.start, niter=niterblkopt,
                     indms)

# now, setup hreg, adjusting for dummies
    if (blkopt_est2$llfval[niterblkopt+1] <= blkopt_est1$llfval[niterblkopt+1] & blkopt_est3$llfval[niterblkopt+1] <= blkopt_est1$llfval[niterblkopt+1]){
        
        hreg <- hregime.reg2.mle(h, m, p, TT=(n-p), fp=blkopt_est1$fpH, init.model)
        output <- list(init.model=init.model,
               hreg=hreg,
               Q=blkopt_est1$Qhat,
               fp=blkopt_est1$fpH,
               m=m, p=p, h=h,
               llfval=blkopt_est1$llfval,
               DirectBFGSLastSuccess=blkopt_est1$DirectBFGSLastSuccess)
        } else if(blkopt_est1$llfval[niterblkopt+1] <= blkopt_est2$llfval[niterblkopt+1] & blkopt_est3$llfval[niterblkopt+1] <= blkopt_est2$llfval[niterblkopt+1]){
        
        hreg <- hregime.reg2.mle(h, m, p, TT=(n-p), fp=blkopt_est2$fpH, init.model)
        output <- list(init.model=init.model,
               hreg=hreg,
               Q=blkopt_est2$Qhat,
               fp=blkopt_est2$fpH,
               m=m, p=p, h=h,
               llfval=blkopt_est2$llfval,
               DirectBFGSLastSuccess=blkopt_est2$DirectBFGSLastSuccess)
        } else{
        hreg <- hregime.reg2.mle(h, m, p, TT=(n-p), fp=blkopt_est3$fpH, init.model)
        output <- list(init.model=init.model,
               hreg=hreg,
               Q=blkopt_est3$Qhat,
               fp=blkopt_est3$fpH,
               m=m, p=p, h=h,
               llfval=blkopt_est3$llfval,
               DirectBFGSLastSuccess=blkopt_est3$DirectBFGSLastSuccess)
        }

    # param_opt里面的参数顺序，取决于llf.msvar的需求
    # 因为原始的llf.msar比较奇怪，这里把llf.msar也改了一下，加上optstr="all"，并且
    param_opt <- array(NA, c(1+m*p+m, m, h))
    param_opt[1:m*p,,] <- output$hreg$Bk[1:m*p,,]
    param_opt[m*p+1,,] <- output$hreg$Bk[1+m*p,,]
    # 这里需要注意，在R的index里面，:比+有更高的优先计算级，所以当我们的index设计加法运算时，必须加上括号
    print(param_opt[(2+m*p):(m+m*p+1),,])
    print(output$Q)
    param_opt[(2+m*p):(m+m*p+1),,] <- output$hreg$Sigmak

    output_theta <- array(NA, c((1+m*p+m),m,h))
    print(output_theta[1:(m*p+1),,])
    print(output$hreg$Bk)
    output_theta[1:(m*p+1),,] <- output$hreg$Bk
    output_theta[(1+m*p+1):(1+m*p+m),,] <- output$hreg$Sigmak

    print("let us get X and Y")
    Y <- init.model$Y[(m+1+1):nrow(init.model$Y),]
    X <- init.model$X[(m+1+1):nrow(init.model$X),]
    print("got X and Y")
    
    print(p)
    print(param_opt)
    print(p)
    print(output_theta)
    print(output$Q)

    output_theta = aperm(output_theta, c(2, 1, 3))
    # optim_result <- fdHess(pars=param_opt, fun=llf_msar, Y=Y, X=X, p=p, theta=output_theta,Q=output$Q, optstr='all', ms.switch=indms)$Hessian
    optim_result <- optim(par=param_opt, fn=llf_msar, gr=NULL, Y=Y, X=X, p=p, theta=output_theta,Q=output$Q, optstr='all', ms.switch=indms,hessian=TRUE)$Hessian
    
    print("got optim_result")
    std <- sqrt(abs(diag(solve(optim_result))))
    print("got std")
    print(std)
    output$hessian <- std
    class(output) <- "MSVAR"

return(output)

} # end mlemsvar() function


# Ryan adjusted several items from the original hregime.reg2 function
hregime.reg2.mle <- function(h, m, p, TT, fp, init.model)
{

    # Storage
    tmp <- vector(mode="list", length=h)
    Bk <- array(0, c(m*p+1, m, h))
    Sigmak <- array(0, c(m,m,h))
    df <- apply(fp, 2, sum)
    e <- array(0, c(TT, m, h))
    Y <- init.model$Y[(m+1+1):nrow(init.model$Y),]
    X <- init.model$X[(m+1+1):nrow(init.model$X),]

    # Loops to compute
    # 1) sums of squares for X and Y
    # 2) B(k) matrices
    # 3) Residuals
    # 4) Sigma(k) matrices

    for(i in 1:h)
    {
        # Note how the dummy obs. are appended to the moment matrices
        Sxy <- crossprod(X, diag(fp[,i]))%*%Y + crossprod(init.model$X[1:(m+1),], init.model$Y[1:(m+1),])
        Sxx <- crossprod(X, diag(fp[,i]))%*%X + crossprod(init.model$X[1:(m+1),])

        # Compute the regression coefficients
        hstar <- Sxx# + init.model$H0
        #Bk[,,i] <- solve(hstar,
       #                  (Sxy + init.model$H0[,1:m]))
        Bk[,,i] <- solve(hstar,Sxy,tol=1e-100)

        # Compute residuals and Sigma (based on Krolzig)

        # Get the full residuals -- need these for filtering
        e[,,i] <- Y - X%*%Bk[,,i]

        #Sigmak[,,i] <- (init.model$S0 + crossprod(e[,,i],diag(fp[,i]))%*%e[,,i])/df[i]
        Sigmak[,,i] <- (crossprod(e[,,i],diag(fp[,i]))%*%e[,,i])/df[i]

        # Save the moments
        tmp[[i]] <- list(Sxy=Sxy, Sxx=Sxx) #, ytmp=ytmp, xtmp=xtmp)
    }

    return(list(Bk=Bk, Sigmak=Sigmak, df=df, e=e, moment=tmp))
}

llf_msar <- function(param_opt, Y, X, p, theta, Q, optstr, ms.switch) {

  print("Start llf_msar")
  m <- ncol(Y)
  n <- nrow(Y) + p
  h <- nrow(Q)

  # initially assign values from theta
  beta0 <- array(theta[,1,], c(m,1,h))
  betap <- NULL
  if (p > 0) betap <- array(theta[,2:(1+m*p),], c(m,m*p,h))
  sig2  <- array(theta[,(1+m*p+1):ncol(theta),], c(m,m,h))
  Qhat  <- Q

  # now choose the parameter over which we are optimizing
  if (optstr=='beta0') {
    beta0 <- array(param.opt, c(m,1,h))
  } else if (optstr=='betap') {
    if (p > 0) betap <- array(param.opt, c(m,m*p,h))
  } else if (optstr=='sig2') {
      sig2  <- array(NA, c(m,m,h))
      # number of distinct: m*(m+1)/2
      nd <- (m*(m+1)/2)
      for (i in 1:h) {
        low <- param.opt[(1+(i-1)*nd):(nd+(i-1)*nd)]
        sig2[,,i] <- xpnd(low, nrow=m)  # user-defined function below
      }
  } else if (optstr=='Qhat') {
    # only passing in first h-1 columns, so add column
    Qhat <- matrix(param.opt, nrow=h, ncol=h-1)
    Qhat <- cbind(Qhat, 1-rowSums(Qhat))
  } else if (optstr=='all'){
    # passing all the estimation in param.opt
    Qhat <- matrix(param.opt[(2+m*p):(m+m*p+1),,], nrow=h, ncol=h-1)
    Qhat <- cbing(Qhat, 1-rowSums(Qhat))

    beta0 <- array(param.opt[(m*p+1),,],c(m,1,h))
    betap <- array(param.opt[1:m*p,,],c(m,m*p,h))
    sig2 <- array(param.opt[(1+m*p+1):(1+m*p+m),,],c(m,m,h))
  }

  # numerical checks on Q matrix
  # prevents elements in Q from going negative or greater than 1
  # if that occurs during optimization, then just set to previous Q
  if ( (min(Qhat) <= 0.0001) || (max(Qhat) >= 0.9999 )) Qhat <- Q

  # numerical checks on sig2
  # prevents elements in sig2 from going negative
  # if that occurs during optimization, then just set to previous sig2
  if ( (min(sig2) <= 0.0001) ) sig2 <- array(theta[,(1+m*p+1):ncol(theta),], c(m,m,h))
  # constrain max(off-diagonal) to be less than min(diagonal)
  blnUseOld <- FALSE
  if (m>1) {
    blnSetPast = 0
    for (i in 1:h) {
      if (max(sig2[,,i][lower.tri(sig2[,,i], diag=FALSE)]) > min(diag(sig2[,,i]))) blnSetPast = 1
    }
    if (blnSetPast == 1) blnUseOld <- TRUE
  }
  if (blnUseOld==TRUE) sig2 <- array(theta[,(1+m*p+1):ncol(theta),], c(m,m,h))

  # by default, everything switches, so now adjust by
  # assigning values that do not switch to first state
  # intercept only
  if (ms.switch=='I') {
    if (p > 0) betap <- array(betap[,,1], c(m,m*p,h))
    sig2  <- array(sig2[,,1], c(m,m,h))
  } else if (ms.switch=='H') { # heteroskedastic
    beta0 <- array(beta0[,,1], c(m,1,h))
    if (p > 0) betap <- array(betap[,,1], c(m,m*p,h))
  } else if (ms.switch=='A') {
    beta0 <- array(beta0[,,1], c(m,1,h))
    sig2  <- array(sig2[,,1], c(m,m,h))
  } else if (ms.switch=='IA') { # homoskedastic
    sig2  <- array(sig2[,,1], c(m,m,h))
  } else if (ms.switch=='IH') {
    if (p > 0) betap <- array(betap[,,1], c(m,m*p,h))
  } else if (ms.switch=='AH') {
    beta0 <- array(beta0[,,1], c(m,1,h))
  }

  print("Value assigned sucessfully")

  #############################################
  # Filtering section
  # Note: there is both Fortran and native R
  #       code in the package (they parallel).
  # Use Fortran for speed, but R for
  # pedagogical purposes.
  #############################################

  # First, obtain residuals
  # e[,,i] is an (n-p) x m matrix containing conditional means
  # for regime i, i=1,2,...,h (e's third dimension is regime)
  e <- array(NA, c(n-p, m, h))
  betahat <- array(NA, c(m,1+m*p,h))
  betahat[,1,] <- beta0
  if (p>0) betahat[,2:(1+m*p),] <- betap
  # adjust for univariate vs. multivariate case
  if (m>1) {
    for (i in 1:h) { e[,,i] <- Y - X %*% t(betahat[,,i]) }
  } else {
    for (i in 1:h) { e[,,i] <- Y - X %*% betahat[,,i] }
  }
  print("got residuals")

  # Using residuals, get filtered regime probabilities
  # HamFilt <- filter.Hamresid(e, sig2.it, Qhat.it)  # R code
  HamFilt <- .Fortran("HamiltonFilter",
         bigt=as.integer(n),
         m = as.integer(m), p = as.integer(p), h = as.integer(h),
         e = e,
         sig2 = sig2,
         Qhat = Qhat,
         f = double(1),
         filtprSt = matrix(0,as.integer(n-p),as.integer(h))
         )

  print("got HamFilt")
  f <- HamFilt$f

  return(-f) # optim() minimizes negative

}
