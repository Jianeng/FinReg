\name{cv_Breg}
\alias{cv_Breg}
\alias{cv.Breg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select regularization parameters by cross-validation, AIC and BIC
}
\description{
This function uses VAR(1) model.
}
\usage{
cv_Breg(x,alpha=1,type="cv",width=NULL,lambda,gamma,tol=10^-3,Time=30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{See \code{Breg}.}
  \item{type}{\code{"cv"}(default), \code{"aic"},\code{"bic"}.}
  \item{width}{Length of the rolling window in computing 1-step ahead prediction errors when \code{type = "cv"}.}
  \item{lambda}{See \code{Breg}.}
  \item{gamma}{See \code{Breg}.}
  \item{alpha}{See \code{Breg}.}
  \item{tol}{See \code{Breg}.}
  \item{Time}{See \code{Breg}.}
}


\value{
\item{lambda}{Selected lambda.}
\item{gamma}{Selected gamma.}
\item{map}{MSE/AIC/BIC grid calculated to choose the lambda and gamma from.}
}
\seealso{
\code{\link{Breg}}
}


