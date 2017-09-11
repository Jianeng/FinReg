# cvCom -----
cvCom <- function(x,lambda,gamma,alpha,T,width,Time)
{
  L = T-width
  prediction = rep(NA,L)

  for(t in 1:L)
  {
    newx = x[t:(t+width-1),]
    beta = Breg(x = newx,
                lambda = lambda,gamma = gamma,
                alpha = alpha,Time = Time)$beta
    prediction[t] = mean(newx[,1]) + beta[1,]%*%(newx[width,]-colMeans(newx))
  }

  mean((x[(1+width):T,1] - prediction)^2)
}

# icCom -----
icCom <- function(x,type,lambda,gamma,alpha,T,Time)
{
  beta = Breg(x = x,
              lambda = lambda,gamma = gamma,
              alpha=alpha,Time=Time)$beta
  x.mean = colMeans(x)
  prediction = (x[1:(T-1),] - x.mean)%*%beta[1,] + x.mean[1]

  P = sum(beta[1,]!=0)
  Likelihood = (T-1)*log(mean((x[2:T,1] - prediction)^2))

  ifelse(type == "aic",
         Likelihood + 2*P,
         Likelihood + log(T-1)*P)
}

# cv.map -----
cv_Breg <- function(x,alpha=1,type="cv",width=NULL,lambda,gamma,tol=10^-3,Time=30)
{
  # find tune parameters by cross-validation
  # width: width of the moving window
  # type: cv, aic, bic

  T = dim(x)[1]
  k = dim(x)[2]
  map = matrix(NA,length(lambda),length(gamma))

  for(i in 1:length(lambda))
    for(j in 1:length(gamma))
    {
      cat("alpha=",alpha,"  lambda_range:",i,"  gamma_range",j,"\n")

      if(type=="cv")
        map[i,j] = cvCom(x,lambda[i],gamma[j],alpha,T,width,Time=Time)
      else
        map[i,j] = icCom(x,type=type,lambda[i],gamma[j],alpha,T,Time = Time)
    }


  # find the "best" tune parameters by MSE
  min_index = which(map == min(map), arr.ind = TRUE)

  return(list("lambda"=lambda[min_index[1]],
              "gamma"=gamma[min_index[2]],
              "map"=map))
}

