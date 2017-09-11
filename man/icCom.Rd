\name{icCom}
\alias{icCom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute AIC and BIC of VAR(1) model
}
\description{
This function calculates the AIC and BIC of estimated VAR(1) model.
}
\usage{
icCom(x,type,lambda,gamma,alpha,T,Time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{x matrix as in \code{Breg}.}
\item{type}{"aic", "bic".}
  \item{lambda}{See \code{Breg}.}
  \item{gamma}{See \code{Breg}.}
  \item{alpha}{See \code{Breg}.}
  \item{T}{Number of observations.}
  \item{Time}{See \code{Breg}.}
}

\seealso{
\code{\link{Breg}}
}
