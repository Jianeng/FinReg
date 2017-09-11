path <- function(x,y=NULL,lambda,gamma,r.type = "var",p.type = "coef",alpha=1,Time=5)
{
  # r.type: var, group, sur
  # p.type: coef, pred, corr,impul

  N.l = length(lambda)
  N.g = length(gamma)
  T = dim(x)[1]
  k = dim(x)[2]
  m = dim(y)[2]
  x.mean = colMeans(x)
  if(!is.null(y))
    y.mean = colMeans(y)

  # result format setting
  if(p.type == "coef")
  {
    if(r.type == "var")
      path_data = matrix(NA,N.l,k)

    if(r.type == "sur" | r.type == "group")
      path_data = matrix(NA,N.l,k*m)
  }

  if(p.type == "pred")
  {
    if(r.type == "var")
      path_data = rep(NA,N.l)

    if(r.type == "sur" | r.type == "group")
      path_data = matrix(NA,N.l,m)
  }

  if(p.type == "corr")
  {
    if(r.type == "var")
      path_data = matrix(NA,N.g,k)

    if(r.type == "sur" | r.type == "group")
      path_data = array(NA,dim=c(m,m,N.g))
  }

  if(p.type == "impul")
  {
    path_data = matrix(NA,N.g,k)
  }

  # coefficient path
  if(p.type == "coef")
  {
    if(r.type == "var")
    {
      for (i in 1:N.l)
      {
        fit = Breg(x = x,
                       lambda = lambda[i],gamma = gamma,
                       alpha = alpha,Time = Time)
        path_data[i,] = fit$beta[1,]
      }
    }

    if(r.type == "sur")
    {
      for (i in 1:N.l)
      {
        fit = Breg(x = x,y = y,
                       lambda = lambda[i],gamma = gamma,
                       alpha=alpha,Time = Time)
        path_data[i,] = fit$beta
      }
    }

    if(r.type == "group")
    {
      model = glmnet::glmnet(x,y,lambda = lambda,family = "mgaussian",alpha=alpha)
      for(i in 1:m)
        path_data[,1:k+(i-1)*k] = as.matrix(t(model$beta[[i]]))
    }
  }

  # prediction path
  if(p.type == "pred")
  {
    if(r.type == "var")
    {
      for (i in 1:N.l)
      {
        fit = Breg(x = x,
                       lambda = lambda[i],gamma = gamma,
                       alpha=alpha,Time = Time)
        path_data[i] = x.mean[1] + (fit$beta%*%(x[T,]-x.mean))[1]
      }
    }

    if(r.type == "sur")
    {
      for (i in 1:N.l)
      {
        fit = Breg(x = x,y = y,
                       lambda = lambda[i],gamma = gamma,
                       alpha=alpha,Time = Time)
        path_data[i,]= (x[T,]-x.mean) %*% matrix(fit$beta,nrow=k) + y.mean
      }
    }

    if(r.type == "group")
    {
      model = glmnet::glmnet(x,y,lambda = lambda,family = "mgaussian",alpha=alpha)
      prediction = predict(model,newx = t(x[T,]),s=lambda)
      dim(prediction) = c(m,N.l)
      path_data = prediction
    }

  }

  # correlation path
  if(p.type == "corr")
  {
    if(r.type == "var")
    {
      for (i in 1:N.g)
      {
        fit = Breg(x = x,
                       lambda = lambda,gamma = gamma[i],
                       alpha=alpha,Time = Time)
        path_data[i,] = cov2cor(fit$sigma)[1,]
      }
    }

    if(r.type == "sur")
    {
      for (i in 1:N.g)
      {
        fit = Breg(x = x,y=y,
                       lambda = lambda,gamma = gamma[i],
                       alpha=alpha,Time = Time)
        path_data[,,i] = cov2cor(fit$sigma)
      }
    }

    if(r.type == "group")
    {
      model = glmnet::glmnet(x,y,lambda = 0,family = "mgaussian")
      res = y - predict(model,newx = x, s = 0)[,,1]
      resid_cov = var(matrix(res,ncol=m))
      for(i in 1:N.g)
      {
        fit = spcov::spcov(diag(diag(resid_cov)),resid_cov,step.size=10,
                      thr.inner=10^-7,lambda = gamma[i]-diag(gamma[i],m))
        path_data[,,i] = cov2cor(fit$Sigma)
      }
    }

  }

  # 1-step impulse path
  if(p.type == "impul")
  {
    if(r.type == "var")
    {
      for (i in 1:N.g)
      {
        fit = Breg(x = x,
                       lambda = lambda,gamma = gamma[i],
                       alpha=alpha,Time = Time)
        path_data[i,]=(fit$beta%*%t(chol(fit$sigma)))[1,]
      }
    }

    if(r.type == "sur")
    {
      for (i in 1:N.g)
      {
        fit = Breg(x = x,y=y,
                       lambda = lambda,gamma = gamma[i],
                       alpha=alpha,Time = Time)
        path_data[i,]=(fit$beta%*%t(chol(fit$sigma)))[1,]
      }
    }

    if(r.type == "group")
    {
      model = glmnet::glmnet(x,y,lambda = 0,family = "mgaussian")
      res = y - predict(model,newx = x, s = 0)[,,1]
      resid_cov = var(matrix(res,ncol=m))
      for(i in 1:N.g)
      {
        fit = spcov::spcov(diag(diag(resid_cov)),resid_cov,step.size=10,
                    thr.inner=10^-7,lambda = gamma[i]-diag(gamma[i],m))
        path_data[i,]=(fit$beta%*%t(chol(fit$sigma)))[1,]
      }
    }

  }

  return(path_data)
}
