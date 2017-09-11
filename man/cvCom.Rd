\name{cvCom}
\alias{cvCom}
\title{
Compute MSE of predictions in a moving window
}
\description{
This function is used in \code{cv.Breg} when \code{type = "cv"}. Estimation an VAR(1) model in a moving window and make a 1-step ahead prediction. The MSE of prediciton is computed.
}
\usage{
cvCom(x,lambda,gamma,alpha,T,width,k,Time)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{VAR(1) observations in a selected period.}

  \item{width}{Bandwidth of the moving window.}

  \item{...}
}

\seealso{
\code{\link{Breg}}
}
