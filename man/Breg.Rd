\name{Breg}
\alias{Breg}
\title{
Regularize Bayesian predictive regression
}
\description{
Given regularization parameters, this function uses an iterative method to estimate the VAR(1) model, 

\deqn{X_t = \beta X_{t-1} + ep_t
}
where all eigenvalues of \eqn{\beta} are smaller than one in modulus. The noise vector is assumed to be multivariate normally distributed with mean 0 and covariance \eqn{\Sigma}.

or SUR model

\deqn{y_i = X^T \beta_{i} + \epsilon_{i}
}
where \eqn{y_i} is N x 1, \eqn{\beta_i} is K x 1 and \eqn{\epsilon_i} is N x 1; \eqn{X^T} is the design matrix. i = 1,2, ..., m.
}
\usage{
Breg(x,y,lambda,gamma,alpha,type)
}

\arguments{
  \item{x}{For \code{type = "var"}, an N x K matrix containing observations of the VAR(1) model; 
  
  For \code{type = "sur"}, an N x K design matrix.}
  \item{y}{For \code{type = "var"}, NULL; 
  
  For \code{type = "sur"}, an N x M response matrix. }
  \item{lambda}{Regularization parameter for \eqn{\beta} estimation. Must be a nonnegative scalar.}
  \item{gamma}{Regularization parameter for \eqn{\Sigma} estiamtion. Must be a nonnegative scalar. The penalty is defined as \eqn{||P*\Sigma||_1} where \eqn{P} is a K × K matrix of all 1 and with 0 on the diagonal to ensure the positive-definiteness.}
  \item{alpha}{Elastic-net mixing parameter as in \code{glmnet}, with \eqn{0 \le \alpha \le 1}. The coefficient penalty is defined as
  \deqn{\lambda*{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1}.}
  \code{alpha=1} is the lasso penalty, and \code{alpha=0} is the ridge penalty.}
  \item{type}{Model type: \code{"var"}(default), \code{"sur"}.}
  \item{tol}{Convergence threshold for the iterative estimation. The convergence criteria is defined as 
  \deqn{||\beta^(i)-\beta^(i-1)||_m \le tol & ||\Sigma^(i)-\Sigma^(i-1)||_m \le tol}  where 'm' denotes the maximum modulus of all the elements.}
  \item{Time}{Maximum number of iterations.}
}

\value{
\item{beta}{The estimate of coefficients. 

For \code{type = "var"}, an K x K matrix; 

For \code{type = "sur"}, a vector of length K*M.}

\item{sigma}{The estimate of noise variance matrix.

For \code{type = "var"}, an K x K matrix.

For \code{type = "sur"}, an M x M matrix.}
}
\references{
Guanhao Feng and Nicholas G. Polson (2016), Regularizing Bayesian Predictive Regressions. \url{https://arxiv.org/abs/1606.01701}
}
\seealso{
\code{\link{glmnet},\link{spcov}}
}
\examples{
## simulating data
library(mvtnorm)
beta <- matrix(c(0.9,-0.1,-0.1,0.8),2,2)
sigma <- matrix(c(1,-0.5,-0.5,1),2,2)

x <- t(rep(0,2))
for (i in 1:1000)
  x<-rbind(x,t(beta\%*\%x[i,])+rmvnorm(1,mean=rep(0,2),sigma))

x <- x[-1,]

## estimating the model
Breg(x,lambd=0,gamma=0,alpha=1,tol=10^-3,Time=30)
}

