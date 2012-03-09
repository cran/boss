\name{updateR2}
\Rdversion{1.1}
\alias{updateR2} 

\title{
%%  ~~function to do ... ~~
Internal function to update choleski decomposition.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is a modified version of the function `updateR()' function from the 'lars' package. This version can handle rank >= 1 updates at a time.
}
\usage{
updateR2(xnew, R = NULL, xold, eps = .Machine$double.eps)
}
\arguments{
	\item{xnew}{ variable(s) to be added to choleski decomposition}
	\item{R}{ choleski decomposition of crossprod(xold)}
	\item{xold}{ old design matrix}
	\item{eps}{ limit of what to consider a zero}
}
\references{
  Trevor Hastie and Brad Efron (2012). lars: Least Angle Regression, Lasso and Forward Stagewise. R package version 1.1.
  http://CRAN.R-project.org/package=lars
}
\author{
Arend Voorman
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{boss.fit}}
}
%\keyword{ ~htest }
%\keyword{ ~regression }
