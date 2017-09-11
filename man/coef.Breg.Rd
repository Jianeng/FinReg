\name{coef.Breg}
\alias{coef.Breg}
\title{
Extract VAR(1) and SUR coefficient and residual covariance from glmnet model
}
\description{
This function extracts the coefficient matrix \eqn{\beta} and residual covariance matrix \eqn{\Sigma} of VAR(1) and SUR model from an \code{glmnet} object.
}
\usage{
coef.Breg(model,lambda,X,Y,type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{An \code{glmnet} model.}

  \item{lambda}{Regularization parameter used in the estimation of \eqn{\beta}.}

  \item{X,Y}{\code{Breg} input.}

  \item{type}{Model type: \code{"var"}(default), \code{"sur"}.}
}


\value{
\item{beta}{Coefficient matrix/vector.}
\item{recov}{Residual covariance matrix.}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{coef.glmnet}, \code{Breg}
}

