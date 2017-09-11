\name{cv.Breg}
\alias{cv.Breg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select regularization parameters by cross-validation, AIC and BIC
}

\usage{
cv.Breg(x,alpha=1,type,width,lambda,gamma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{\code{"cv"}(default), \code{"aic"}, \code{"bic"}.}
  \item{width}{Length of the rolling window in computing 1-step ahead prediction errors when \code{type = "cv"}.}
  \item{...}
}

\details{
This function uses VAR(1) model.
}
\value{
\item{lambda}{Selected lambda.}
\item{gamma}{Selected gamma.}
\item{map}{MSE/AIC/BIC grid calculated to choose the lambda and gamma from.}
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

## set regularization parameters
gamma <- seq(0,0.1,by=0.05)
lambda <- seq(0,0.1,by=0.05)

## cross validation
cv.Breg(x,alpha=1,type="cv",width=500,lambda,gamma)
}

