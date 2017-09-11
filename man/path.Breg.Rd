\name{path}
\alias{path}

\title{
Regularizaiton path of Bayesian predictive regression 
}
\description{
This functions computes the Bayesian regularization path for coefficient, prediction, correlation and 1-step ahead impulse, given a sequence of lambda or gamma.
}
\usage{
path(x,y=NULL,lambda,gamma,r.type,p.type,alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{For \code{r.type = "var"}, an N x K matrix.
  
  For \code{r.type = "sur"} or \code{r.type = "group"}, an N x K design matrix.}
  
  \item{y}{For \code{r.type = "var"}, NULL; 
  
  For \code{r.type = "sur"} or \code{r.type = "group"}, an N x M response matrix. }
  
  \item{r.type}{Model type: \code{"var"}(default), \code{"sur"}, \code{"group"}.}
  \item{p.type}{Path type: \code{"coef"}(default), \code{"pred"}, \code{"corr"}, \code{"impul"}.}
  \item{...}
}
\details{
Only one of \code{lambda} and \code{gamma} can be an increasing sequence and another one is a non-negative scalar (typical 0).

For \code{p.type = "coef","pred"}, \code{lambda} is a sequence.

For \code{p.type = "corr","impul"}, \code{gamma} is a sequence.
}
\value{
\item{path}{A vector/matrix/array of regularization result. The first dimension size equals to the length of parameter sequence.}
}

\seealso{
\code{\link{Breg}}
}

\examples{
## simulating data 
library(mvtnorm)
beta <- matrix(c(0.9,-0.1,-0.1,0.8),2,2)
sigma <- matrix(c(1,-0.5,-0.5,1),2,2)

x <- t(rep(0,2))
for (i in 1:1000)
  x <- rbind(x,t(beta\%*\%x[i,])+rmvnorm(1,mean=rep(0,2),sigma))
x <- x[-1,]

## set the regularization parameter lambda
lambda = seq(0,0.1,by=0.01)

## compute regularization path
path(x=x,gamma=0,lambda=lambda,p.type="pred",r.type="var")
}