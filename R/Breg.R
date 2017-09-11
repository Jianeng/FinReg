library(glmnet)
library(spcov)

# Sigma_inverse_half ------------------
Sigma_inverse_half <- function(sigma)
{
  # get inverse half of sigma by SVD

  s = svd(sigma)
  sigmaD <- 1/(sqrt(s$d))
  sigmaU <- s$u

  sigmaU%*%diag(sigmaD)%*%t(sigmaU)
}

# Get_coef ------------------
coef.Breg <- function(model,lambda,X,Y,type="var")
{
  # recover coefficient and residual covariance from glmnet model
  if(type == "var")
  {
    k = dim(X)[2]
    IK = diag(rep(1,k))
    LI = kronecker(X,IK)
    C1 = coef(model)[-1,10]
    Z = Y - LI%*%C1
    Z1 = t(matrix(Z,k))
    resid_cov = var(Z1)
    beta = matrix(C1,k,k)
  }

  if(type == "sur")
  {
    m = dim(Y)[2]
    beta = model$beta[,10]
    res = Y - predict(model,newx=X,s=lambda)
    resid_cov = var(matrix(res,ncol=m))
  }

  return(list("beta"=beta,"recov"=resid_cov))
}

# var.estimation ------------------
Breg <- function(x,y=NULL,lambda,gamma,
                     alpha=1,type = "var",
                     tol=10^-7,Time=10,r=10^-14)
{
  # estimate beta and sigma with fixed tune parameters
  # type: var, sur

  # data dimension
  T = dim(x)[1] # number of observations
  k = dim(x)[2] # dimension of observations

  # remove the mean
  x.mean = colMeans(x)
  x = sweep(x,2,x.mean)

  # X & Y setting
  if(type == "var")
  {
    X = x[1:(T-1),]
    Y = t(x)[(k+1):(k*T)]
  }

  if(type == "sur")
  {
    m = dim(y)[2]
    y.mean = colMeans(y)
    y = sweep(y,2,y.mean)
    X = bdiag(rep(list(x),m))
    Y = c(y)
  }

  # lambda for glmnet
  l = seq(exp(20),lambda,length.out = 10)

  # initial estimates
  sigma_temp = diag(k)
  if(type == "var")
    beta_temp = matrix(1,k,k)
  if(type == "sur")
    beta_temp = rep(1,k*m)


  # VAR iteration ------------------------------
  for (t in 1:Time)
  {
    cat("VAR iterate step:", t,"\n")

    # update estimates
    beta = beta_temp
    sigma = sigma_temp

    # set up equation
    sigma_half = Sigma_inverse_half(sigma)
    if(type == "var")
    {
      I = diag(rep(1,times=T-1))
      newX = kronecker(X,sigma_half)
      newY = kronecker(It,sigma_half) %*% Y
    }

    if(type == "sur")
    {
      omega_half = kronecker(sigma_half,diag(T))
      newX = omega_half %*% X
      newY = omega_half %*% Y
    }

    # estimate beta
    model = glmnet(newX,newY,lambda=l,alpha=alpha,intercept=FALSE)

    # get coefficients
    res_1 = coef.Breg(model,lambda,X,Y,type=type)
    beta_temp = res_1$beta
    sigma_temp = res_1$recov

    # estimate sigma (penalize correlation matrix)
    if(gamma!=0)
    {
      corr = cov2cor(sigma_temp)
      half_var = diag(diag(sigma_temp)^(.5))

      if(rcond(corr)<r)
      {
        print("residual correlation matrix is singular!\n")
        exit()
      }

      penalty = gamma - diag(gamma,k)
      cor_temp = spcov(diag(diag(corr)),corr,lambda = penalty,
                   step.size = step.size)
      sigma_temp = half_var %*% temp$Sigma %*% half_var
    }

    change.beta = norm(beta-beta_temp,type="M")
    change.sigma = norm(sigma-sigma_temp,type="M")

    if(change.beta < tol & change.sigma < tol)
      break
  }

  return(list("beta"=beta_temp,"sigma"=sigma_temp))
}
