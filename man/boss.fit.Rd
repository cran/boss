\name{boss.fit}
\Rdversion{1.1}
\alias{boss.fit} 
\alias{boss.fit.fgee} 
\alias{boss.fit.fgee.GxE} 
\alias{boss.fit.flm}
\alias{boss.fit.flm.GxE}
\alias{boss.fit.flr}
\alias{boss.fit.flr.GxE}
\alias{boss.fit.fgee.lr}
\alias{boss.fit.fgee.lr.GxE}
\alias{boss.fit.lmm}
\alias{boss.fit.lmm.GxE}
\alias{boss.fit.smcg}
\alias{boss.fit.smvg}
\alias{updateR}
\alias{updateR2}
\title{
%%  ~~function to do ... ~~
fit a model with genotype
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Uses an object from 'boss.set' and calls an appropriate fitter function to add genotype to the model
}
\usage{
boss.fit(g, init, thresh = 1e-7, robust=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
A genotype variable, assumed to be in the same order as the outcome variable used in `boss.set()'.  For the 'swap' method, this can be a matrix, 
}
  \item{init}{
%%     ~~Describe \code{init} here~~
An object from `boss.set()'
}
 \item{thresh}{
 p-value threshold beyond which further iteratively re-weighted least squares iterations should be performed.  Ignored for OLS and linear GEE, where estimates are fully converged, also ignored for "swap" method.
}
\item{robust}{
 logical.  Whether or not robust standard errors should be used.  Defaults to TRUE.  Ignored for mixed models, where they do not exist, and "swap" method.
} 
 \item{\dots}{
%%     ~~Describe \code{\dots} here~~
optional arguments to be passed to the fitter functions if further iterations are desired.
}
}
\value{
  \item{beta.main}{estimated coefficient for the main effect of genotype}
  \item{v.main}{the estimated variance of beta.main}
  \item{beta.inter}{when interaction term is supplied, the estimated coefficient for the interaction with genotype}
  \item{v.inter}{when interaction term is supplied, the estimated variance of beta.inter}
  \item{cov.inter}{when interaction term is supplied, the estimated covariance between beta.inter and beta.main}
}

\author{
Arend Voorman
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{boss.set}}
\code{\link{glm}}
\code{\link{geese}}
\code{\link{lmer}}
}
\examples{
#generate data:
n <- 500
nsnps <- 1000
snpMat <- replicate(nsnps,rbinom(n,2,.2))
colnames(snpMat) <- 1e6+(1:nsnps)
rownames(snpMat) <- 1:n

#generate 3 observations per person with random intercept and covariates:
id <- as.numeric(gl(n,3))
y <- rnorm(3*n) + rnorm(n,sd=.3)[id]
data <- data.frame(y=y, x1=rnorm(3*n,mean=y), x2 =rnorm(3*n), site=sample(6,n,replace=TRUE)[id], id=id )

# set up a GEE with AR-1 correlation structure, using "chol" method:
init1 <- boss.set(y~x1+x2+factor(site), id=data$id, type = "gee", method = "chol", corstr= "ar1", data=data)

# set up a linear mixed model (lmm) with random intercept using "swap" method:
init2 <- boss.set(y~x1+x2+factor(site) + (1|id), type = "lmm", method = "swap", data=data)

#set up an interaction with logistic regression, which must use the "chol" method:
y <- rbinom(n,1,0.5) 
data.bin <- data.frame(y=y, x1=rnorm(n,mean=y), x2 =rnorm(n), e = rnorm(n), site=sample(6,n,replace=TRUE)[id], id=1:n)

init3 <- boss.set(y~x1+x2+factor(site)+e, E.name = "e", id=id, family=binomial(), type="glm", data = data.bin)

##fit the model (at one locus):
g <- snpMat[1,]

###GEE
#one step approximation:
boss.fit(g[id], init1, robust=TRUE)

#full iteration
boss.fit(g[id], init1, robust=TRUE, thresh=1)


##LMM
#one step approximation:
boss.fit(g[id], init2)

#full iteration (uses lmer, since "swap" method does not permit further iterations currently)
mod <- summary(lmer(update(init2$formula,.~.+g[id]),data=data))
list(beta.main = mod@fixef[9], v.main=mod@vcov[9,9], chi2 = mod@fixef[9]^2/mod@vcov[9,9])

##Logistic Reg.
boss.fit(g[id], init3)

#full iteration
boss.fit(g[id], init3, thresh = 1)

#for swap method, matrix arguments can be used to simultaneously fit many models:
results <- boss.fit(t(snpMat[,id]),init2)
}
%\keyword{ ~htest }
%\keyword{ ~regression }
